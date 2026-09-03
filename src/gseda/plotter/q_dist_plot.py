#!/usr/bin/env python3
"""Plot Q-value distribution at raw and normalized positions from BAM or FASTQ files.

Two subplots (both heatmaps):
  - Left:  Q-score at raw (0-indexed) positions within each read
  - Right: Q-score at normalized positions (0.0 to 1.0 relative to read length)
"""
from __future__ import annotations

import argparse
import math
import os
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from tqdm import tqdm


def detect_format(path: str) -> str:
    """Auto-detect file format from extension."""
    p = pathlib.Path(path)
    ext = p.suffix.lower()
    if ext == ".bam":
        return "bam"
    if ext in (".fq", ".fastq"):
        return "fastq"
    # Try without extension for gzipped files
    stem = p.stem
    if stem.endswith((".gz", ".xz", ".bz2")):
        stem = stem.rsplit(".", 1)[0]
    if stem.lower().endswith(".bam"):
        return "bam"
    if stem.lower().endswith((".fq", ".fastq")):
        return "fastq"
    return "auto"


def extract_from_bam(
    bam_path: str,
    threads: int = 1,
    min_rl: int = 1,
    max_rl: int | None = None,
) -> tuple[list[tuple[int, int]], list[tuple[float, int]]]:
    """Extract (position, Q) pairs from BAM, return raw and normalized lists."""
    raw_pairs: list[tuple[int, int]] = []
    norm_pairs: list[tuple[float, int]] = []

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=threads) as bam:
        for read in tqdm(bam.fetch(until_eof=True), desc=f"reading {bam_path}"):
            if read.is_unmapped:
                continue
            quals = read.query_qualities
            L = len(quals)
            if L < min_rl:
                continue
            if max_rl is not None and L > max_rl:
                continue
            for pos, q in enumerate(quals):
                raw_pairs.append((pos, q))
                norm_pairs.append((pos / (L - 1) if L > 1 else 0.0, q))

    return raw_pairs, norm_pairs


def extract_from_fastq(
    fastq_path: str,
    min_rl: int = 1,
    max_rl: int | None = None,
) -> tuple[list[tuple[int, int]], list[tuple[float, int]]]:
    """Extract (position, Q) pairs from FASTQ, return raw and normalized lists."""
    raw_pairs: list[tuple[int, int]] = []
    norm_pairs: list[tuple[float, int]] = []

    with pysam.FastxFile(fastq_path) as fh:
        for entry in tqdm(fh, desc=f"reading {fastq_path}", unit="read"):
            if entry.quality is None:
                continue
            quals = [ord(c) - 33 for c in entry.quality]
            L = len(quals)
            if L < min_rl:
                continue
            if max_rl is not None and L > max_rl:
                continue
            for pos, q in enumerate(quals):
                raw_pairs.append((pos, q))
                norm_pairs.append((pos / (L - 1) if L > 1 else 0.0, q))

    return raw_pairs, norm_pairs


def _plot_heatmap(
    ax: plt.Axes,
    x_values: list[float | int],
    y_values: list[int],
    x_min: float,
    x_max: float,
    x_bins: int,
    y_min: int = 0,
    y_max: int = 40,
    title: str = "",
    x_label: str = "",
    y_label: str = "Q-score",
) -> None:
    """Draw a 2D heatmap (x vs Q-score, color = count) with log color scale."""
    x_edges = np.linspace(x_min, x_max, x_bins + 1)
    y_edges = np.arange(y_min, y_max + 2)  # +1 so max bin edge > y_max
    X, Y, Z = np.histogram2d(x_values, y_values, bins=[x_edges, y_edges])

    # Apply log1p for better visual dynamic range
    Z = np.log1p(Z)
    # Replace zeros with -1 so they show as background color
    Z = np.ma.masked_where(Z == 0, Z)

    dx = (x_max - x_min) / x_bins
    dy = y_edges[1] - y_edges[0]

    im = ax.pcolormesh(
        X, Y, Z,
        shading="auto",
        cmap="viridis",
        vmin=0,
        vmax=np.ma.max(Z) * 0.95 if np.ma.max(Z) > 0 else 1,
    )
    fig = ax.get_figure()
    fig.colorbar(im, ax=ax, label="log1p(count)")

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(axis="both", labelsize=11)


def plot_q_distribution(
    raw_pairs: list[tuple[int, int]],
    norm_pairs: list[tuple[float, int]],
    output_path: str,
    dpi: int = 300,
) -> str:
    """Create figure with two heatmaps and save to output_path."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 7))

    max_pos = max(p for p, _ in raw_pairs) if raw_pairs else 0
    # Adaptive x bins: cap at 200 for readability
    x_bins_raw = min(max_pos + 1, 200) if max_pos > 0 else 1
    x_bins_norm = 200

    _plot_heatmap(
        ax1,
        x_values=[p for p, _ in raw_pairs],
        y_values=[q for _, q in raw_pairs],
        x_min=-0.5,
        x_max=max_pos + 0.5,
        x_bins=x_bins_raw,
        title="Raw Position",
        x_label="Position in Read",
    )

    _plot_heatmap(
        ax2,
        x_values=[p for p, _ in norm_pairs],
        y_values=[q for _, q in norm_pairs],
        x_min=0.0,
        x_max=1.0,
        x_bins=x_bins_norm,
        title="Normalized Position",
        x_label="Normalized Position (0 = start, 1 = end)",
    )

    fig.suptitle(f"Q-value Distribution  |  Total bases: {len(raw_pairs):,}", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return output_path


def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plot Q-value distribution at raw and normalized positions.",
    )
    parser.add_argument("input_file", help="Path to BAM or FASTQ file")
    parser.add_argument("-o", "--output", default="q_dist.png", help="Output PNG path (default: q_dist.png)")
    parser.add_argument("--min-read-len", type=int, default=1, help="Minimum read length (default: 1)")
    parser.add_argument("--max-read-len", type=int, default=None, help="Maximum read length (default: no limit)")
    parser.add_argument("--threads", type=int, default=1, help="Threads for BAM (default: 1)")
    parser.add_argument("--format", choices=["auto", "bam", "fastq"], default="auto", help="Force format (default: auto)")
    parser.add_argument("--dpi", type=int, default=300, help="Output DPI (default: 300)")
    args = parser.parse_args(argv)

    fmt = args.format if args.format != "auto" else detect_format(args.input_file)

    if fmt == "bam":
        raw, norm = extract_from_bam(args.input_file, threads=args.threads,
                                     min_rl=args.min_read_len, max_rl=args.max_read_len)
    elif fmt == "fastq":
        raw, norm = extract_from_fastq(args.input_file, min_rl=args.min_read_len,
                                       max_rl=args.max_read_len)
    else:
        print(f"Error: could not auto-detect format for {args.input_file}", file=sys.stderr)
        return 1

    if not raw:
        print("Warning: no bases found matching criteria.", file=sys.stderr)
        return 1

    out = plot_q_distribution(raw, norm, args.output, dpi=args.dpi)
    print(f"Saved to {out} ({len(raw):,} bases)")
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
