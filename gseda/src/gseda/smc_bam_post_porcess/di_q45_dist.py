#!/usr/bin/env python3
"""Distance distribution of Q45 `di` units from the ends of the smc consensus.

For every `di` unit (pos,base,frac,depth,phreq) in the smc consensus BAM whose
`phreq` rounds to the target (default 45), compute how far the insertion
position `pos` is from the *nearest* end of the read:

    dist = min(pos, len(seq) - pos)

(i.e. distance from the start if the unit is nearer the 5' end, distance from
the tail otherwise).  The resulting distances are plotted as a histogram.

The input BAM is NOT an alignment BAM (unmapped smc consensus), so it is
opened with `check_sq=False` and iterated via `fetch(until_eof=True)`.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pysam  # noqa: E402
from tqdm import tqdm  # noqa: E402

# --- project palette (validated, see references/palette.md) -------------------
COLOR_DIST = "#2a78d6"   # distance to nearest end (categorical slot 1)
SURFACE = "#fcfcfb"
GRID = "#e1e0d9"
INK = "#0b0b0b"

# --- reuse the tolerant `di` parser from the sibling module --------------------
try:
    from .insert_di import parse_di
except ImportError:  # pragma: no cover - only if sibling is missing
    def parse_di(di_string):
        """Fallback 5-field tolerant parser (mirrors insert_di.parse_di)."""
        if not di_string:
            return []
        out = []
        for seg in str(di_string).split(";"):
            seg = seg.strip()
            if not seg:
                continue
            parts = seg.split(",")
            if len(parts) != 5:
                continue
            try:
                phreq = float(parts[4])
            except (ValueError, IndexError):
                continue
            out.append((int(parts[0]), parts[1].upper(),
                        float(parts[2]), int(parts[3]), phreq))
        return out


# ---------------------------------------------------------------------------
# Collection
# ---------------------------------------------------------------------------

def collect_end_distances(
    input_bam: str,
    target_phreq: int = 45,
    threads: int = 4,
    show_progress: bool = True,
) -> tuple[list[int], int, int]:
    """Scan the smc consensus BAM; return distances of matching `di` units.

    Returns (distances, n_read, n_matched_units) where distances are
    min(pos, len(seq) - pos) for every `di` unit whose phreq rounds to
    target_phreq.
    """
    distances: list[int] = []
    n_read = 0
    n_matched = 0

    with pysam.AlignmentFile(input_bam, mode="rb", check_sq=False, threads=threads) as in_f:
        iterator = in_f.fetch(until_eof=True)
        if show_progress:
            iterator = tqdm(iterator, desc=f"collect {pathlib.PurePosixPath(input_bam).name}")

        for record in iterator:
            n_read += 1
            if not record.has_tag("di"):
                continue
            seq_len = len(record.query_sequence or "")
            if seq_len == 0:
                continue

            for (pos, _base, _frac, _depth, phreq) in parse_di(record.get_tag("di")):
                if round(float(phreq)) != target_phreq:
                    continue
                distances.append(min(pos, seq_len - pos))
                n_matched += 1

    return distances, n_read, n_matched


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_histogram(
    distances: list[int],
    out_png: str,
    target_phreq: int = 45,
    log_y: bool = False,
    dpi: int = 150,
) -> str:
    """Plot the distance samples as a histogram (one bin per nt) and save PNG."""
    hi = max(distances) if distances else 0
    edges = np.arange(0, hi + 2, 1)

    fig = plt.figure(figsize=(9, 5.5), facecolor=SURFACE)
    ax = fig.add_subplot(1, 1, 1, facecolor=SURFACE)

    ax.hist(
        distances,
        bins=edges,
        color=COLOR_DIST,
        alpha=0.8,
        label="distance to nearest end",
        edgecolor=SURFACE,
    )

    if log_y:
        ax.set_yscale("log")
    ax.set_xscale("linear")

    ax.set_xlabel("Distance from nearest read end (nt)", color=INK, fontsize=11)
    ax.set_ylabel("di units" + (" (log scale)" if log_y else ""), color=INK, fontsize=11)
    ax.set_title(f"Q{target_phreq} di units: distance from start vs tail of smc consensus",
                 color=INK, fontsize=12, fontweight="bold")

    ax.tick_params(axis="both", colors=INK, labelsize=9, which="both")
    ax.grid(True, which="major", axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)

    ax.legend(loc="upper right", frameon=False, fontsize=10,
              facecolor="none", edgecolor="none")
    for text in ax.get_legend().get_texts():
        text.set_color(INK)

    plt.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return out_png


def _default_output(input_bam: str) -> str:
    p = pathlib.PurePosixPath(input_bam)
    stem = p.name
    if stem.endswith(".bam"):
        stem = stem[:-4]
    return str(p.with_name(f"{stem}.di_q45_dist.png"))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plot the distribution of Q45 `di` units' distance from "
                    "the nearest end (start or tail) of the smc consensus."
    )
    parser.add_argument("input_bam", help="Input smc consensus BAM (unmapped).")
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output PNG path (default: <input>.di_q45_dist.png).",
    )
    parser.add_argument(
        "--phreq", type=int, default=45,
        help="Target phreq to match (default: 45).",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for BAM I/O (default: 4).",
    )
    parser.add_argument(
        "--linear", action="store_true",
        help="Use a linear (counts) y-axis. Default: log scale.",
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
    out = args.output or _default_output(args.input_bam)

    distances, n_read, n_matched = collect_end_distances(
        args.input_bam,
        target_phreq=args.phreq,
        threads=args.threads,
        show_progress=not args.no_progress,
    )

    if not distances:
        print(f"No `di` units with phreq={args.phreq} found — nothing to plot.",
              file=sys.stderr)
        return 1

    saved = plot_histogram(distances, out, target_phreq=args.phreq, log_y=not args.linear)
    d = np.asarray(distances)
    print(f"records: {n_read:,}, matched di units (phreq={args.phreq}): {n_matched:,}")
    print(f"distance: min={d.min()}, median={int(np.median(d))}, "
          f"mean={d.mean():.1f}, max={d.max()}")
    print(f"within 50 nt of an end: {int((d <= 50).sum()):,} "
          f"({(d <= 50).mean() * 100:.1f}%)")
    print(f"saved: {saved}")
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
