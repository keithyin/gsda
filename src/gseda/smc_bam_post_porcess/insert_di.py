#!/usr/bin/env python3
"""Insert bases recorded in the smc `di` tag back into the consensus sequence.

The smc output BAM carries a `di` tag (a string) that records bases which have
subread support but are missing from the smc consensus sequence (`S`).

This script parses `di`, inserts each recorded base at its position in both
the sequence and its quality vector (using the provided `phreq` as the qual),
and recomputes the whole-read quality as stored in the `rq` tag.

The input and output BAMs are NOT alignment BAMs (no CIGAR / reference
coordinates).

`di` tag format (authoritative, per insert_di.md):
    "pos,base,frac,depth,phreq;pos,base,frac,depth,phreq;..."
    - pos:    int   insertion position in the smc seq
    - base:   ACGT  base to insert
    - frac:   float support fraction
    - depth:  int   support depth
    - phreq:  float -10*log(1-p), a raw phred error score (no further scaling)

`pos` may repeat, meaning multiple supported bases land at the same position;
they are all inserted in their `di` order.

With `--poly-only` a unit is inserted only when its base equals `seq[pos]`
(i.e. it extends an existing homopolymer run); non-matching units are kept
out, including trailing ones at `pos == len(seq)` which have no base to
compare against.

With `--q-cut Q` (an integer quality) a unit is inserted only when its own
`phreq` reaches `Q` (`phreq >= Q`, compared against the raw, unrounded phred
score).  The two filters compose: a unit has to pass both.  `Q` is also baked
into the default output name, `<input>.with_di.q<Q>.bam`.

`rq` (read quality, probability form) is recomputed as the mean per-base
correct rate:  rq = mean(1 - 10**(-q/10))  — matching
`gseda/file_format_cvt/fastq2bam.py:calc_rq`.
"""

import argparse
import array
import math
import pathlib
import sys
from collections import OrderedDict
from typing import List, Tuple

import pysam
from tqdm import tqdm

_BASES = "ACGT"


def parse_di(di_string: str) -> List[Tuple[int, str, float, int, float]]:
    """Parse a `di` tag into a list of (pos, base, frac, depth, phreq).

    Malformed units (wrong field count, bad base, non-numeric) are skipped,
    matching the tolerant style of `reinsert_model.read_di_tag`.
    """
    if not di_string:
        return []
    entries = []
    for segment in di_string.split(";"):
        segment = segment.strip()
        if not segment:
            continue
        parts = segment.split(",")
        if len(parts) != 5:
            continue
        try:
            pos = int(parts[0])
            base = parts[1].strip().upper()
            frac = float(parts[2])
            depth = int(parts[3])
            phreq = float(parts[4])
        except (ValueError, IndexError):
            continue
        if base not in _BASES:
            continue
        entries.append((pos, base, frac, depth, phreq))
    return entries


def _clamp(x: int, lo: int = 0, hi: int = 50) -> int:
    return max(lo, min(hi, x))


def calc_rq(quals) -> float:
    """Compute the whole-read correctness rate (probability form).

    rq = mean(1 - 10**(-q/10)) over bases, rounded to 6 decimals — identical
    to `gseda/file_format_cvt/fastq2bam.py:calc_rq`.
    """
    if quals is None or len(quals) == 0:
        return 0.0
    acc = sum(1.0 - math.pow(10.0, -q / 10.0) for q in quals) / len(quals)
    return round(acc, 6)


def filter_poly_only(
    seq: str,
    entries: List[Tuple[int, str, float, int, float]],
) -> List[Tuple[int, str, float, int, float]]:
    """Keep only entries whose base equals the base already at that position.

    A `di` unit is a homopolymer extension only when its base matches
    `seq[pos]`; anything else is a non-poly insertion and is dropped.  Units
    with `pos >= len(seq)` have no base to compare against and are dropped too.
    """
    n = len(seq)
    return [
        entry for entry in entries
        if 0 <= entry[0] < n and entry[1] == seq[entry[0]]
    ]


def filter_q_cut(
    entries: List[Tuple[int, str, float, int, float]],
    q_cut: int,
) -> List[Tuple[int, str, float, int, float]]:
    """Keep only entries whose `phreq` reaches the integer quality `q_cut`.

    `phreq` is already a raw phred score, so it is compared as-is (against the
    unrounded value) and the cut is inclusive: `phreq == q_cut` is kept.  A
    `q_cut` of 0 keeps everything, since `-10*log(1-p)` is never negative.
    """
    return [entry for entry in entries if entry[4] >= q_cut]


def insert_at_positions(
    seq: str,
    quals,
    entries: List[Tuple[int, str, float, int, float]],
) -> Tuple[str, "array.array"]:
    """Insert `di` bases at their positions, keeping `di` ordering.

    For each position, the recorded bases are inserted *at* that index
    (i.e. before `seq[pos]`); multiple entries at the same position are
    appended in their `di` order.  A `pos == len(seq)` group lands at the end.

    Returns (new_seq, new_qual) where new_qual is an array('B').
    """
    by_pos: "OrderedDict[int, List[Tuple[str, float]]]" = OrderedDict()
    for (pos, base, _frac, _depth, phreq) in entries:
        by_pos.setdefault(pos, []).append((base, phreq))

    out_bases: List[str] = []
    out_qual: List[int] = []

    n = len(seq)
    original_quals = list(quals) if quals is not None else [0] * n

    for i in range(n):
        for (base, phreq) in by_pos.get(i, ()):
            out_bases.append(base)
            out_qual.append(_clamp(round(float(phreq))))
        out_bases.append(seq[i])
        out_qual.append(original_quals[i])

    # Trailing inserts at pos == len(seq).
    for (base, phreq) in by_pos.get(n, ()):
        out_bases.append(base)
        out_qual.append(_clamp(round(float(phreq))))

    return "".join(out_bases), array.array("B", out_qual)


def process_bam(
    input_bam: str,
    output_bam: str,
    threads: int = 40,
    show_progress: bool = True,
    poly_only: bool = False,
    q_cut: int = 0,
) -> Tuple[int, int]:
    """Read `input_bam`, insert `di` bases, recompute `rq`, write `output_bam`.

    With `poly_only`, only `di` bases matching the consensus base at their
    position are inserted (see `filter_poly_only`).  `q_cut` additionally
    requires a unit's own `phreq` to reach that value (see `filter_q_cut`).

    Returns (n_modified, n_total).
    """
    n_modified = 0
    n_total = 0

    with pysam.AlignmentFile(
        input_bam, mode="rb", check_sq=False, threads=threads
    ) as in_f, pysam.AlignmentFile(
        output_bam, mode="wb", check_sq=False, threads=threads, header=in_f.header
    ) as out_f:

        iterator = in_f.fetch(until_eof=True)
        if show_progress:
            iterator = tqdm(iterator, desc=f"insert_di {input_bam}")

        for record in iterator:
            n_total += 1
            if record.has_tag("di"):
                di = record.get_tag("di")
                entries = parse_di(di)
                if entries:
                    seq = record.query_sequence or ""
                    if q_cut > 0:
                        entries = filter_q_cut(entries, q_cut)
                    if poly_only:
                        entries = filter_poly_only(seq, entries)
                    if entries:
                        quals = record.query_qualities
                        new_seq, new_qual = insert_at_positions(seq, quals, entries)
                        record.query_sequence = new_seq
                        record.query_qualities = new_qual
                        record.set_tag("rq", calc_rq(new_qual), value_type="d")
                        n_modified += 1
            out_f.write(record)

    return n_modified, n_total


def _default_output_path(input_bam: str, q_cut: int = 0) -> str:
    """`<stem>.with_di.q<N>.bam` — the cut is always part of the name.

    Tagging unconditionally means a run at `q_cut=0` (no filter) and one at
    `q_cut=20` never overwrite each other, and every output states the
    threshold that produced it.
    """
    p = pathlib.PurePosixPath(input_bam)
    stem = p.name.rsplit(".", maxsplit=1)[0]
    return str(p.with_name(f"{stem}.with_di.q{q_cut}.bam"))


def main_cli():
    parser = argparse.ArgumentParser(
        description="Insert `di`-tag bases back into the smc consensus "
                    "sequence, recompute per-read quality, update `rq`, "
                    "and write a new (unmapped) BAM."
    )
    parser.add_argument("input_bam", help="Input smc BAM (with `di` tags)")
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output BAM path (default: <input>.with_di.q<Q>.bam, where Q is "
             "the --q-cut value)",
    )
    parser.add_argument(
        "--threads", type=int, default=40,
        help="Number of threads for BAM I/O (default: 40)",
    )
    parser.add_argument(
        "--no-progress", action="store_true", help="Disable progress bar",
    )
    parser.add_argument(
        "--poly-only", action="store_true",
        help="Only insert `di` bases that equal the consensus base at their "
             "position (homopolymer-run extensions); drop the other units",
    )
    parser.add_argument(
        "--q-cut", type=int, default=0, metavar="Q",
        help="Minimum integer quality Q a `di` unit's phred score must reach "
             "to be inserted (inclusive; default: 0, i.e. no quality filter). "
             "Also tags the default output name as *.with_di.q<Q>.bam",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    out = args.output or _default_output_path(args.input_bam, args.q_cut)

    n_mod, n_tot = process_bam(
        input_bam=args.input_bam,
        output_bam=out,
        threads=args.threads,
        show_progress=not args.no_progress,
        poly_only=args.poly_only,
        q_cut=args.q_cut,
    )
    print(f"Processed {n_tot} reads; modified {n_mod}; output: {out}")


if __name__ == "__main__":
    main_cli()
