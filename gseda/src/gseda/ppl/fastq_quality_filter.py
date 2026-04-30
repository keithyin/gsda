#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Filter FASTQ reads by quality and plot read length distribution.

Accepts one or more FASTQ files (supports glob patterns), filters reads
by minimum average quality score, and produces a length distribution plot.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
import glob as glob_module

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def _phred_to_error_prob(q: int) -> float:
    """Phred Q -> base error probability: Pe = 10^(-Q/10)."""
    return 10.0 ** (-q / 10.0)


def _error_prob_to_phred(pe: float) -> float:
    """Error probability -> Phred Q: Q = -10 * log10(Pe)."""
    if pe <= 0:
        return float("inf")
    if pe >= 1:
        return 0.0
    return -10.0 * math.log10(pe)


def _phred_to_accuracy(q: int) -> float:
    """Phred Q -> base-calling accuracy: Pa = 1 - 10^(-Q/10)."""
    return 1.0 - _phred_to_error_prob(q)


def _accuracy_to_phred(pa: float) -> float:
    """Accuracy -> Phred Q: convert accuracy back to Phred score."""
    pe = 1.0 - pa
    return _error_prob_to_phred(pe)


def _qual_str_to_phreds(qual_str: str) -> list[int]:
    """FASTQ quality string -> list of Phred scores (Phred+33)."""
    return [ord(c) - 33 for c in qual_str]


def parse_fastq(paths: list[str]) -> pl.DataFrame:
    """Read FASTQ files (supports glob patterns) and return a Polars DataFrame."""
    names, seqs, quals = [], [], []

    for path in paths:
        expanded = glob_module.glob(path)
        if not expanded:
            logging.warning("no files matched: %s", path)
            continue

        for fpath in expanded:
            with open(fpath, "r") as fh:
                while True:
                    header = fh.readline().strip()
                    if not header:
                        break
                    seq = fh.readline().strip()
                    fh.readline()  # + line
                    qual = fh.readline().strip()
                    if len(header) < 1 or not seq or not qual:
                        continue
                    names.append(header)
                    seqs.append(seq)
                    quals.append(qual)

    return pl.DataFrame({
        "name": names,
        "seq": seqs,
        "qual": quals,
    })


def filter_by_min_quality(df: pl.DataFrame, min_q: int) -> pl.DataFrame:
    """Keep reads where every base has quality >= min_q."""
    if df.is_empty():
        return df

    all_qual = df["qual"].to_list()
    mask = [min(_qual_str_to_phreds(q)) >= min_q for q in all_qual]
    return df.filter(pl.Series(mask))


def filter_by_avg_quality(df: pl.DataFrame, min_avg: float) -> pl.DataFrame:
    """Keep reads whose average Phred quality >= min_avg.

    Method: convert each base's Phred score to accuracy, average the
    accuracies across the read, then convert back to Phred.
    """
    if df.is_empty():
        return df

    all_qual = df["qual"].to_list()
    mask = []
    for q in all_qual:
        phreds = _qual_str_to_phreds(q)
        accuracies = [_phred_to_accuracy(qi) for qi in phreds]
        avg_acc = sum(accuracies) / len(accuracies)
        avg_q = _accuracy_to_phred(avg_acc)
        mask.append(avg_q >= min_avg)
    return df.filter(pl.Series(mask))


def plot_length_dist(df: pl.DataFrame, output: str):
    """Draw and save the read length distribution histogram."""
    lengths = df.select(pl.col("seq").str.len_bytes().alias("length"))

    if not lengths.is_empty():
        arr = lengths["length"].to_numpy()
        min_l = int(arr.min())
        max_l = int(arr.max())
        mean_len = arr.mean()
        median_len = np.median(arr)
        total = len(arr)
    else:
        min_l = max_l = 0
        mean_len = median_len = 0
        total = 0
        
    min_l = 0
    max_l = 4000

    plt.figure(figsize=(12, 6))
    if total > 0:
        sns.histplot(data=pl.DataFrame({"length": arr}), x="length",
                     bins=max(1, max_l - min_l + 1), stat="count",
                     kde=False, color="#4c72b0", edgecolor="white")

    plt.xlabel("Read Length (bp)")
    plt.ylabel("Count")
    plt.title("Filtered Read Length Distribution")

    stats_text = (
        f"Total: {total} reads\n"
        f"Mean: {mean_len:.1f} bp\n"
        f"Median: {median_len:.1f} bp"
    )
    plt.text(0.02, 0.95, stats_text, transform=plt.gca().transAxes,
             fontsize=11, verticalalignment="top",
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    plt.close()
    logging.info("plot saved to %s", output)


def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Filter FASTQ reads by quality and plot length distribution.",
    )
    parser.add_argument("fastq", nargs="+",
                        help="FASTQ file(s), supports glob patterns")
    parser.add_argument("-o", "--output", default="fastq_length_dist.png",
                        help="output plot filename (default: fastq_length_dist.png)")
    parser.add_argument("--min-avg-q", type=float, default=None,
                        help="minimum average Phred quality per read")
    parser.add_argument("--min-q", type=int, default=None,
                        help="minimum Phred quality for every base in a read")
    parser.add_argument("--length-range", type=str, default=None,
                        help="filter by read length range, e.g. '100:500'")
    args = parser.parse_args(argv)

    if args.min_avg_q is None and args.min_q is None:
        logging.warning("no quality filter specified; all reads will be included")

    df = parse_fastq(list(args.fastq))
    logging.info("loaded %d reads", len(df))

    if args.min_q is not None:
        df = filter_by_min_quality(df, args.min_q)
        logging.info("after min-q=%d filter: %d reads", args.min_q, len(df))

    if args.min_avg_q is not None:
        df = filter_by_avg_quality(df, args.min_avg_q)
        logging.info("after avg-q>=%.1f filter: %d reads", args.min_avg_q, len(df))

    if args.length_range:
        lo, hi = args.length_range.split(":")
        lo, hi = int(lo), int(hi)
        mask = pl.Series([lo <= len(s) <= hi for s in df["seq"].to_list()])
        df = df.filter(mask)
        logging.info("after length %d:%d filter: %d reads", lo, hi, len(df))

    plot_length_dist(df, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
