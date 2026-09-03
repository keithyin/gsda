#!/usr/bin/env python3
"""Check whether `di` units sit in homopolymer context of the smc consensus.

Verdict rule (as specified): a `di` unit (pos, base) is *correct* when the
inserted base matches the consensus at the insertion point, or at the position
right before it:

    base == seq[pos]   or   base == seq[pos-1]

Both cases mean the base joins a run of the same base that already exists in
the consensus — `seq[pos]` extends that run rightwards, `seq[pos-1]` leftwards.
Anything else is reported as *wrong*.

The rule is looser than "is in a homopolymer", though: a single adjacent copy
(`run == 1`) satisfies it without being a homopolymer in the usual sense.  So
each unit is additionally scored by the length of the run of `base` that is
contiguous with the insertion point, and the report breaks the correct units
down by it.  `run >= 2` is the strict reading of the question.

`di` tag format (authoritative, per insert_di.md):
    "pos,base,frac,depth,phreq;pos,base,frac,depth,phreq;..."
    `pos` may repeat; it is an index into the smc seq, and `pos == len(seq)`
    means a trailing insertion (only `seq[pos-1]` can then be compared).

`prev` is checked further when it is the sole reason a unit passes, since that
is the one case the rule can excuse for the wrong reason.  A contest at `pos`
is expected to be behind it (see `prev_supported_by_multi`), so every such unit
is classified as *explained* — an earlier unit at the same pos carries a
different base — or *anomalous* — none does.  A unit matching both `seq[pos]`
and `seq[pos-1]` sits in a run the consensus already shows, so it is counted in
`both` and needs no such excuse.

The input BAM is NOT an alignment BAM (unmapped smc consensus), so it is opened
with `check_sq=False` and iterated via `fetch(until_eof=True)`.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
from collections import Counter

import pysam
from tqdm import tqdm

# --- reuse the tolerant `di` parser from the sibling module --------------------
try:
    from .insert_di import parse_di
except ImportError:  # pragma: no cover - only if run as a standalone script
    from insert_di import parse_di


# ---------------------------------------------------------------------------
# The check
# ---------------------------------------------------------------------------

def classify(seq: str, pos: int, base: str) -> tuple[bool, bool, int]:
    """Score one `di` unit against the consensus `seq`.

    Returns (at, prev, run):
      at   — `base == seq[pos]`      (insertion point itself matches)
      prev — `base == seq[pos-1]`    (position before it matches)
      run  — number of `base` copies contiguous with the insertion point,
             i.e. left-run + right-run; 0 exactly when neither `at` nor `prev`
             holds, so `run > 0` is the verdict and `run >= 2` the strict
             homopolymer reading of it.

    Positions outside the sequence are never correct: `pos < 0`, and
    `pos > len(seq)`, leave nothing to compare against (a trailing unit at
    `pos == len(seq)` still has `seq[pos-1]`, so it is allowed).
    """
    n = len(seq)
    if pos < 0 or pos > n:
        return False, False, 0

    at = pos < n and seq[pos] == base
    prev = pos > 0 and seq[pos - 1] == base

    left = pos - 1
    while left >= 0 and seq[left] == base:
        left -= 1
    right = pos
    while right < n and seq[right] == base:
        right += 1
    return at, prev, right - left - 1


def prev_supported_by_multi(base: str, earlier_bases) -> bool:
    """Whether `base == seq[pos-1]` is excused by a co-inserted unit at the pos.

    A repeated `pos` means the site is contested: several subreads support
    bases the consensus dropped.  Only one of them can be read as extending
    rightwards from `seq[pos]`, so the rest can only satisfy the rule
    leftwards, against `seq[pos-1]` — `prev` is then the expected outcome of
    the collision rather than evidence of a homopolymer.  When no earlier unit
    shares the `pos`, nothing contested it, and a `prev` that does not also
    satisfy `at` is anomalous.

    `earlier_bases` is the bases of the units listed before this one at the
    same `pos`; a differing base is required, a repeat of the same base says
    nothing about the slot.
    """
    return any(b != base for b in earlier_bases)


def context(seq: str, pos: int, flank: int = 8) -> str:
    """Consensus bases around the insertion point, with `|` marking it."""
    n = len(seq)
    p = max(0, min(pos, n))
    return seq[max(0, p - flank):pos] + "|" + seq[pos:p + flank]


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------

def check_bam(
    input_bam: str,
    min_np: int | None = None,
    min_rq: float | None = None,
    q_cut: float = 20.0,
    n_examples: int = 5,
    dump_bad: str | None = None,
    threads: int = 4,
    show_progress: bool = True,
) -> tuple[Counter, list[dict]]:
    """Scan the smc consensus BAM and tally the verdict for every `di` unit.

    A record is kept only when its `np` tag is >= min_np and its `rq` tag is
    >= min_rq; an absent tag fails that filter (min_* = None disables it).

    Returns (stats, examples) where `examples` holds up to n_examples wrong
    units; if `dump_bad` is set, every wrong unit is appended there as TSV.
    """
    st: Counter = Counter()
    examples: list[dict] = []

    bad_out = None
    if dump_bad:
        bad_out = open(dump_bad, "w")
        bad_out.write("read\tpos\tbase\tphreq\tfrac\tdepth\trq\tnp\tcontext\n")

    channels = []
    
    try:
        with pysam.AlignmentFile(
            input_bam, mode="rb", check_sq=False, threads=threads
        ) as in_f:
            iterator = in_f.fetch(until_eof=True)
            # if show_progress:
            #     iterator = tqdm(iterator,
            #                     desc=f"check {pathlib.PurePosixPath(input_bam).name}")
            

            for record in iterator:
                st["n_read"] += 1
                if min_np is not None and (
                        not record.has_tag("np") or int(record.get_tag("np")) < min_np):
                    continue
                if min_rq is not None and (
                        not record.has_tag("rq") or float(record.get_tag("rq")) < min_rq):
                    continue
                if not record.has_tag("di"):
                    continue

                units = parse_di(record.get_tag("di"))
                if not units:
                    continue

                st["n_read_di"] += 1
                seq = record.query_sequence or ""
                r_rq = float(record.get_tag("rq")) if record.has_tag(
                    "rq") else float("nan")
                r_np = int(record.get_tag("np")
                           ) if record.has_tag("np") else -1
                n_wrong = 0
                # bases of the units already passed, per `pos` — the contested
                # sites that can excuse a `prev` match
                bases_at: dict[int, list[str]] = {}

                for (pos, base, frac, depth, phreq) in units:
                    at, prev, run = classify(seq, pos, base)

                    # Only a `prev` that fails `at` needs excusing: a base
                    # matching both sides sits in a run the consensus already
                    # shows.  The other way it can pass is a contest at `pos`.
                    if prev and not at:
                        prev_expl = prev_supported_by_multi(
                            base, bases_at.get(pos, ()))
                        st["prev_alone_expl" if prev_expl else
                           "prev_alone_anom"] += 1
                    bases_at.setdefault(pos, []).append(base)

                    """
                        "5556-6-65-14-----4565666-666-5
                        "CTGG-A-CG-TT-----GTACTCC-CAG-C
                        "CTGG-A-CG-TT-----GTACTCC-CAG-C
                        "CTGG-ACCG-TTGTA--GTACT-C-CAG-C
                        "CTGG-A-CG-TT-----GTACTCC-CAG-C
                        "CTGG-A-CG-TT-----GTACTCC-CAG-C
                        "CTGG-A-CG-TTGTAAAGTACTCC-CAGCC
                        "CTGG-A-CA-TT-----GTACTCC-CAG-C
                        "CTGG-A-CGTTT-----GTACTCCACA--C
                        "CTGG-A-CG-TT-----G-ACTCC-CAG-C
                        "CTGGCA-CG-TT-----GTACTCC-CAG-C
                        
                        在当前 SMC 逻辑前，为什么会存在 prev = true, 上述给了一个解释。是当一个位点可能存在多个插入的时候
                        这个是符合预期的
                    
                    """
                    if prev and pos < 20:
                        ch = record.get_tag("ch")
                        if ch not in channels:
                            channels.append(ch)

                        if at:
                            verdict = "benign (`at` too — real run)"
                        else:
                            verdict = "explained" if prev_expl else "ANOMALOUS"
                        print(
                            f"{record.query_name}. {(pos, base, frac, depth, phreq)}"
                            f" {verdict}")

                    st["n_unit"] += 1
                    coh = "hi" if phreq > q_cut else "lo"
                    st[f"{coh}_n"] += 1
                    st[f"run_{min(run, 12)}"] += 1

                    if run > 0:
                        st["n_ok"] += 1
                        st[f"{coh}_ok"] += 1
                        if run >= 2:
                            st["n_ok_run2"] += 1
                            st[f"{coh}_run2"] += 1
                        st["n_at"] += int(at)
                        st["n_prev"] += int(prev)
                        st["n_both"] += int(at and prev)
                    else:
                        n_wrong += 1
                        st["n_wrong"] += 1
                        st[f"{coh}_wrong"] += 1
                        st[f"wrong_base_{base}"] += 1
                        if bad_out is not None:
                            bad_out.write(
                                "{}\t{}\t{}\t{:.2f}\t{:.4f}\t{}\t{:.4f}\t{}\t{}\n".format(
                                    record.query_name, pos, base, phreq, frac,
                                    depth, r_rq, r_np, context(seq, pos)))
                        if len(examples) < n_examples:
                            examples.append({
                                "name": record.query_name, "pos": pos, "base": base,
                                "phreq": phreq, "frac": frac, "depth": depth,
                                "rq": r_rq, "np": r_np, "rlen": len(seq),
                                "ctx": context(seq, pos),
                            })

                if n_wrong:
                    st["n_read_wrong"] += 1
                else:
                    st["n_read_ok"] += 1
    finally:
        if bad_out is not None:
            bad_out.close()
    print("{}".format(",".join(map(str, channels))))
    return st, examples


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def report(st: Counter, examples: list[dict], q_cut: float,
           min_np: int | None, min_rq: float | None) -> None:
    """Print the verdict table, the run-length breakdown and wrong examples."""
    W = 74
    print("=" * W)
    print("di HOMOPLYMER-CONTEXT CHECK")
    print("  correct  <=>  base == seq[pos]  or  base == seq[pos-1]")
    filt = "none"
    if min_np is not None or min_rq is not None:
        parts = []
        if min_np is not None:
            parts.append(f"np>={min_np}")
        if min_rq is not None:
            parts.append(f"rq>={min_rq}")
        filt = ", ".join(parts)
    print(f"  read filter: {filt}")

    def pct(x: int, y: int) -> str:
        return f"{100.0 * x / y:6.2f}%" if y else "     -"

    print("-" * W)
    print("  records scanned       {:>14,}".format(st["n_read"]))
    print("  records with di       {:>14,}".format(st["n_read_di"]))
    print("    all units correct   {:>14,}  {}".format(
        st["n_read_ok"], pct(st["n_read_ok"], st["n_read_di"])))
    print("    >=1 unit wrong      {:>14,}  {}".format(
        st["n_read_wrong"], pct(st["n_read_wrong"], st["n_read_di"])))
    print()
    print("  di units              {:>14,}".format(st["n_unit"]))
    print("    correct             {:>14,}  {}".format(
        st["n_ok"], pct(st["n_ok"], st["n_unit"])))
    print("      seq[pos]   only   {:>14,}".format(st["n_at"] - st["n_both"]))
    print(
        "      seq[pos-1] only   {:>14,}".format(st["n_prev"] - st["n_both"]))
    print("      both               {:>13,}".format(st["n_both"]))
    print("    WRONG               {:>14,}  {}".format(
        st["n_wrong"], pct(st["n_wrong"], st["n_unit"])))
    if st["n_wrong"]:
        comp = " ".join(
            "{}={:,}".format(b, st[f"wrong_base_{b}"]) for b in "ACGT"
            if st[f"wrong_base_{b}"])
        print("      wrong units by base: " + comp)

    print("-" * W)
    print("  run length of `base` contiguous with the insertion point")
    print("  (run==1 satisfies the rule but is a lone copy, not a homopolymer)")
    tot_run = st["n_ok"] + st["n_wrong"]
    for k in range(13):
        lbl = f"  run {k}" if k < 12 else "  run >=12"
        n_run = st[f"run_{k}"]
        if not n_run:
            continue
        note = ""
        if k == 0:
            note = "<- wrong (no neighbouring copy)"
        elif k == 1:
            note = "<- correct, single copy"
        print("{:<12}{:>10,}  {:>9}  {}".format(
            lbl, n_run, pct(n_run, tot_run), note))
    print("  correct AND in a run >=2: {:,}  {}".format(
        st["n_ok_run2"], pct(st["n_ok_run2"], st["n_unit"])))

    print("-" * W)
    print("  `prev` units (base == seq[pos-1]) — excused by an earlier `di` unit")
    print("  at the same pos carrying a different base?  (see prev_supported_by_multi)")
    n_pe = st["prev_alone_expl"]
    n_pa = st["prev_alone_anom"]
    print("  {:<20}{:>11}{:>11}{:>10}".format(
        "", "explained", "anomalous", "anom%"))
    print("  {:<18}{:>11,}{:>11,}{:>10}".format(
        "prev alone", n_pe, n_pa, pct(n_pa, n_pe + n_pa).strip()))
    print("  `at` matches too (real run, needs no excuse): {:,}".format(st["n_both"]))

    print("-" * W)
    print("  {:<12}{:>10}{:>12}{:>12}{:>12}".format(
        "cohort", "units", "correct%", "run>=2%", "wrong%"))
    for lbl, key in ((f">Q{q_cut:g}", "hi"), (f"<=Q{q_cut:g}", "lo")):
        n_coh = st[f"{key}_n"]
        print("  {:<12}{:>10,}{:>12}{:>12}{:>12}".format(
            lbl, n_coh,
            pct(st[f"{key}_ok"], n_coh).strip(),
            pct(st[f"{key}_run2"], n_coh).strip(),
            pct(st[f"{key}_wrong"], n_coh).strip()))
    print("=" * W)

    if not examples:
        return
    print()
    print("WRONG `di` UNITS — first {}:".format(len(examples)))
    for ex in examples:
        print("  {}  pos={} base={} Q{:.1f} frac={:.3f} depth={} np={} rq={:.4f} len={}".format(
            ex["name"][:48], ex["pos"], ex["base"], ex["phreq"], ex["frac"],
            ex["depth"], ex["np"], ex["rq"], ex["rlen"]))
        print("      consensus: ..." + ex["ctx"] + "...")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify that `di` units from an smc consensus BAM land in "
                    "homopolymer context: a unit is correct when its base "
                    "matches seq[pos] or seq[pos-1], wrong otherwise."
    )
    parser.add_argument(
        "input_bam", help="Input smc consensus BAM (unmapped).")
    parser.add_argument(
        "--min-np", type=int, default=3,
        help="Keep records with np >= --min-np (default: 3; 0 keeps all).",
    )
    parser.add_argument(
        "--min-rq", type=float, default=None,
        help="Keep records with rq >= --min-rq (e.g. 0.9). Default: no rq filter.",
    )
    parser.add_argument(
        "--q-cut", type=float, default=20.0,
        help="phreq threshold of the high/low cohort split (default: 20).",
    )
    parser.add_argument(
        "--examples", type=int, default=5,
        help="How many wrong units to print to stdout (default: 5).",
    )
    parser.add_argument(
        "--dump-bad", type=str, default=None,
        help="Also write every wrong unit to this TSV path (default: off).",
    )
    parser.add_argument(
        "--threads", type=int, default=40,
        help="Threads for BAM I/O (default: 40).",
    )
    parser.add_argument(
        "--no-progress", action="store_true",
        help="Disable the tqdm progress bar.",
    )

    args_to_parse = sys.argv[1:] if argv is None else argv
    if len(args_to_parse) == 0:
        parser.print_help(sys.stderr)
        return 1

    args = parser.parse_args(args_to_parse)
    min_np = args.min_np if args.min_np > 0 else None

    st, examples = check_bam(
        args.input_bam,
        min_np=min_np,
        min_rq=args.min_rq,
        q_cut=args.q_cut,
        n_examples=args.examples,
        dump_bad=args.dump_bad,
        threads=args.threads,
        show_progress=not args.no_progress,
    )

    if st["n_unit"] == 0:
        print("No `di` units found — nothing to check.", file=sys.stderr)
        return 1

    report(st, examples, q_cut=args.q_cut, min_np=min_np, min_rq=args.min_rq)
    if args.dump_bad:
        print("\nwrong units ->", args.dump_bad)
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
