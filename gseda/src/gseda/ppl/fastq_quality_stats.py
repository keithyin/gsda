#!/usr/bin/env python3
"""Compute Q20/Q30 read percentages for every FASTQ file in a directory.

Read Q is calculated via accuracy conversion:
  for each read:  avg_acc = mean(1 - 10^(-qi/10))  across all bases
                  read_Q = -10 * log10(1 - avg_acc)
  Q20% = fraction of reads with read_Q >= 20
  Q30% = fraction of reads with read_Q >= 30
"""
from __future__ import annotations

import argparse
import glob
import math
import os
import sys

import pysam
from tqdm import tqdm
import polars as pl


def _phred_to_error_prob(q: int | float) -> float:
    """Phred Q -> base error probability: Pe = 10^(-Q/10)."""
    return 10.0 ** (-q / 10.0)


def _error_prob_to_phred(pe: float) -> float:
    """Error probability -> Phred Q: Q = -10 * log10(Pe)."""
    if pe <= 0:
        return float("inf")
    if pe >= 1:
        return 0.0
    return -10.0 * math.log10(pe)


def _phred_to_accuracy(q: int | float) -> float:
    """Phred Q -> base-calling accuracy: Pa = 1 - 10^(-Q/10)."""
    return 1.0 - _phred_to_error_prob(q)


def _accuracy_to_phred(pa: float) -> float:
    """Accuracy -> Phred Q."""
    pe = 1.0 - pa
    return _error_prob_to_phred(pe)


def _quality_to_phreds(qual) -> list[int]:
    """Convert pysam quality (string or array) to list of Phred scores (Phred+33)."""
    if isinstance(qual, (list, tuple)):
        return [int(q) for q in qual]
    # pysam FastxRecord.quality is a Phred+33 encoded string
    return [ord(c) - 33 for c in qual]


def read_q_stats(filepath: str, rq_thr: float = 0.0) -> tuple[int, int, int]:
    """Return (total_reads, q20_reads, q30_reads) for a FASTQ file.

    Read-level Q is computed by converting each base's Phred score to
    accuracy, averaging across the read, then converting back to Phred.
    """
    total = q20 = q30 = 0
    with pysam.FastxFile(filepath) as fh:
        for rec in tqdm(fh, desc=filepath, unit="read"):
            if rec.quality is None or len(rec.quality) == 0:
                continue
            phreds = _quality_to_phreds(rec.quality)
            accs = [_phred_to_accuracy(q) for q in phreds]
            avg_acc = sum(accs) / len(accs)
            read_q = _accuracy_to_phred(avg_acc)
            if read_q < rq_thr:
                continue
            total += 1
            if read_q >= 20:
                q20 += 1
            if read_q >= 30:
                q30 += 1
    return total, q20, q30


def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Q20/Q30 read percentages per FASTQ file in a directory.",
    )
    parser.add_argument("fastq_dir", help="Directory containing FASTQ files")
    parser.add_argument("--gz", action="store_true", help="Include gzipped FASTQ files")
    parser.add_argument("--rq-thr", type=float, default=0.0,
                        help="read Q threshold, reads below this value are filtered out (default: 0.0, no filter)")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.fastq_dir):
        print(f"Error: {args.fastq_dir} is not a directory", file=sys.stderr)
        return 1

    exts = ["*.fastq", "*.fq"]
    if args.gz:
        exts += ["*.fastq.gz", "*.fq.gz"]

    rows = []
    for ext in exts:
        for fpath in sorted(glob.glob(os.path.join(args.fastq_dir, ext))):
            if os.path.isfile(fpath):
                total, q20, q30 = read_q_stats(fpath, rq_thr=args.rq_thr)
                if total == 0:
                    continue
                q20_pct = q20 / total * 100
                q30_pct = q30 / total * 100
                rows.append({
                    "file": fpath,
                    "total_reads": total,
                    "Q20_pct": round(q20_pct, 4),
                    "Q30_pct": round(q30_pct, 4),
                })

    if rows:
        df = pl.DataFrame(rows)
        print(df)
    if args.rq_thr > 0:
        print(f"\n[filter: excluding reads with read Q < {args.rq_thr}]", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
