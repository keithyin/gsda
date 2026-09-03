#!/usr/bin/env python3
"""Print reads (channels) whose `di` tag contains units with phreq in a range.

The smc output BAM carries a `di` tag: a ";"-separated list of units
`pos,base,frac,depth,phreq` recording bases which have subread support but
are missing from the smc consensus sequence.  This script scans the BAM and,
for every `di` unit whose `phreq` (rounded to int) falls in the inclusive
range given by `--phreq LO:HI` (default 30:45), prints the read name (channel)
and the full `di` info for that unit.

The input BAM is NOT an alignment BAM (unmapped smc consensus), so it is
opened with `check_sq=False` and iterated via `fetch(until_eof=True)`.
"""

from __future__ import annotations

import argparse
import sys

import pysam
from tqdm import tqdm

try:
    from .insert_di import parse_di
except ImportError:  # pragma: no cover - only if run as a script
    from insert_di import parse_di


def _parse_phreq_range(range_str: str) -> tuple[int, int]:
    """Parse a `lo:hi` (inclusive) phred range, e.g. `30:33`.

    A bare value (`45`) is treated as the degenerate range `45:45`.
    """
    parts = range_str.split(":")
    if len(parts) == 1:
        lo = hi = int(parts[0])
    elif len(parts) == 2:
        lo, hi = int(parts[0]), int(parts[1])
    else:
        raise argparse.ArgumentTypeError(
            f"invalid range {range_str!r}: expected `LO:HI`"
        )
    if lo > hi:
        raise argparse.ArgumentTypeError(
            f"invalid range {range_str!r}: lower bound exceeds upper bound"
        )
    return lo, hi


def dump_di_phreq(
    input_bam: str,
    phreq_range: tuple[int, int] = (30, 45),
    threads: int = 4,
    show_progress: bool = True,
) -> tuple[int, int]:
    """Print reads whose `di` tag has a unit with phreq in `phreq_range`.

    The range is inclusive on both ends; `phreq` is rounded to int first.
    Returns (n_matched_reads, n_matched_units).
    """
    phreq_lo, phreq_hi = phreq_range
    n_matched_reads = 0
    n_matched_units = 0

    channels = []

    with pysam.AlignmentFile(
        input_bam, mode="rb", check_sq=False, threads=threads
    ) as in_f:
        iterator = in_f.fetch(until_eof=True)
        # if show_progress:
        #     iterator = tqdm(iterator, desc=input_bam)

        for record in iterator:
            if not record.has_tag("di"):
                continue
            entries = parse_di(record.get_tag("di"))
            rq = record.get_tag("rq")
            hits = [
                (pos, base, frac, depth, phreq)
                for (pos, base, frac, depth, phreq) in entries
                if phreq_lo <= round(float(phreq)) <= phreq_hi
            ]
            if not hits:
                continue

            n_matched_reads += 1
            for (pos, base, frac, depth, phreq) in hits:
                n_matched_units += 1
                if pos > 10 and pos < 20:
                    print(
                        f"{record.query_name}\t{pos}\t{base}\t"
                        f"{frac:.4f}\t{depth}\t{phreq}\t{rq}"
                    )
                    if rq < 0.999:
                        channels.append(record.get_tag("ch"))
    print("channels: {}".format(",".join(list(map(str, channels)))))
    return n_matched_reads, n_matched_units


def main_cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Print channels (read names) and `di` info for `di` units "
                    "whose phreq falls in the target range (default 30:45)."
    )
    parser.add_argument("input_bam", help="Input smc BAM (with `di` tags)")
    parser.add_argument(
        "--phreq", type=_parse_phreq_range, default=(30, 45), metavar="LO:HI",
        help="Inclusive phreq range to match, as `LO:HI` (e.g. `30:33`). "
             "A bare value (e.g. `45`) matches that single value "
             "(default: 30:45).",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for BAM I/O (default: 4).",
    )
    parser.add_argument(
        "--no-progress", action="store_true", help="Disable progress bar",
    )

    if len(sys.argv) == 1 and argv is None:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args(argv)
    n_reads, n_units = dump_di_phreq(
        args.input_bam,
        phreq_range=args.phreq,
        threads=args.threads,
        show_progress=not args.no_progress,
    )
    print(
        f"matched: {n_units} `di` units (phreq in [{args.phreq[0]}, {args.phreq[1]}]) "
        f"across {n_reads} reads",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
