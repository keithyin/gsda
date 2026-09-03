#!/usr/bin/env python3
"""Compare the quality of bases *missing* from the smc consensus (`di` tag)
against the quality of the bases that *are* present in the smc consensus seq.

The smc consensus BAM carries a `di` tag: bases which have subread support but
are absent from the smc consensus sequence.  Each `di` unit is
`pos,base,frac,depth,phreq` (phreq = -10*log(1-p), a raw phred error score);
the `phreq` field is the quality of that missing base.  This script plots, as
a histogram with a legend, the distribution of:

  - `di` qual    -> the `phreq` of every `di` unit
  - `smc seq qual` -> the per-base quality of the consensus sequence

Reads are optionally filtered on the whole-read `rq` and `np` tags before
sampling.

The input (and any related) BAM is NOT an alignment BAM — it is the unmapped
smc consensus, so we open it with `check_sq=False` and use `fetch(until_eof=True)`.

Colors: two categorical series, first two slots of the validated project
palette (blue / orange), colorblind-safe for two hues.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pysam  # noqa: E402
from tqdm import tqdm  # noqa: E402

# --- project palette (validated, see references/palette.md) -------------------
COLOR_DI = "#2a78d6"     # di phreq     (categorical slot 1)
COLOR_SEQ = "#eb6834"    # smc seq qual (categorical slot 2)
SURFACE = "#fcfcfb"
GRID = "#e1e0d9"
INK = "#0b0b0b"

# --- reuse the tolerant `di` parser from the sibling module --------------------
try:
    from .insert_di import parse_di  # 5-field: pos,base,frac,depth,phreq
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

def collect_qual_samples(
    input_bam: str,
    min_rq: float | None = None,
    min_np: int | None = None,
    threads: int = 4,
    show_progress: bool = True,
):
    """Scan the smc consensus BAM; return per-base quality samples.

    Returns (di_values, seq_values, n_read, n_filtered) where
    di_values   = [phreq of every `di` unit],
    seq_values  = [every consensus-seq base quality],
    n_read      = total records seen,
    n_filtered  = records dropped by the rq/np filters.

    A record is kept only if its `rq` (float) is >= min_rq and its `np` (int)
    is >= min_np; a tag that is absent is treated as unfiltered-for that value.
    """
    di_values: list[int] = []
    seq_values: list[int] = []
    n_read = 0
    n_filtered = 0

    with pysam.AlignmentFile(input_bam, mode="rb", check_sq=False, threads=threads) as in_f:
        iterator = in_f.fetch(until_eof=True)
        if show_progress:
            iterator = tqdm(
                iterator, desc=f"collect {pathlib.PurePosixPath(input_bam).name}")

        for record in iterator:
            n_read += 1

            # --- read-level filter on rq / np -----------------------------------
            if min_rq is not None:
                if not record.has_tag("rq") or float(record.get_tag("rq")) < min_rq:
                    n_filtered += 1
                    continue
            if min_np is not None:
                if not record.has_tag("np") or int(record.get_tag("np")) < min_np:
                    n_filtered += 1
                    continue

            # --- smc seq qual (per base) ----------------------------------------
            seq_qual = record.query_qualities
            if seq_qual is not None:
                seq_values.extend(int(q) for q in seq_qual)

            # --- di qual (phreq of each `di` unit) -------------------------------
            if record.has_tag("di"):
                for _pos, _base, _frac, _depth, phreq in parse_di(record.get_tag("di")):
                    int_q = int(round(float(phreq)))
                    if int_q >= 0:
                        di_values.extend([int_q])

    return di_values, seq_values, n_read, n_filtered


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _bin_edges(di_values, seq_values):
    """Integer phred edges covering both series (one bin per phred value)."""
    all_vals = list(di_values) + list(seq_values)
    hi = max(all_vals) if all_vals else 0
    # edges: 0..hi+1 -> bins [0,1),[1,2),...,[hi,hi+1)
    return np.arange(0, hi + 2, 1)


def plot_histogram(
    di_values: list[int],
    seq_values: list[int],
    out_png: str,
    log_y: bool = False,
    dpi: int = 150,
) -> str:
    """Plot the two quality samples as one overlaid histogram and save PNG."""
    edges = _bin_edges(di_values, seq_values)

    fig = plt.figure(figsize=(9, 5.5), facecolor=SURFACE)
    ax = fig.add_subplot(1, 1, 1, facecolor=SURFACE)

    ax.hist(
        [di_values, seq_values],
        bins=edges,
        histtype="bar",
        alpha=0.65,
        color=[COLOR_DI, COLOR_SEQ],
        label=["di phreq", "smc seq qual"],
        edgecolor=SURFACE,
    )

    if log_y:
        ax.set_yscale("log")
    ax.set_xscale("linear")

    ax.set_xlabel("Base quality (phred)", color=INK, fontsize=11)
    ax.set_ylabel("Reads" + (" (log scale)" if log_y else ""),
                  color=INK, fontsize=11)
    ax.set_title("di phreq vs smc consensus seq-qual distribution", color=INK,
                 fontsize=12, fontweight="bold")

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
    return str(p.with_name(f"{stem}.di_vs_seq.png"))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plot a two-series histogram of `di` (missing-base) phreq "
                    "vs smc consensus seq-qual from an smc consensus BAM."
    )
    parser.add_argument(
        "input_bam", help="Input smc consensus BAM (unmapped).")
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output PNG path (default: <input>.di_vs_seq.png).",
    )
    parser.add_argument(
        "--min-rq", type=float, default=None,
        help="Keep records with rq >= --min-rq (e.g. 0.9). Default: no rq filter.",
    )
    parser.add_argument(
        "--min-np", type=int, default=5,
        help="Keep records with np >= --min-np (e.g. 5). Default: no np filter.",
    )
    parser.add_argument(
        "--threads", type=int, default=40,
        help="Threads for BAM I/O (default: 40).",
    )
    parser.add_argument(
        "--linear", action="store_true",
        help="Use a linear (counts) y-axis. Default: log scale, which keeps both "
             "series legible given their very different total sample counts.",
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

    di_vals, seq_vals, n_read, n_filtered = collect_qual_samples(
        args.input_bam,
        min_rq=args.min_rq,
        min_np=args.min_np,
        threads=args.threads,
        show_progress=not args.no_progress,
    )

    if not di_vals and not seq_vals:
        print("No quality samples found — nothing to plot.", file=sys.stderr)
        return 1

    filters = []
    if args.min_rq is not None:
        filters.append(f"rq >= {args.min_rq}")
    if args.min_np is not None:
        filters.append(f"np >= {args.min_np}")
    fl = f" [filter: {', '.join(filters)}]" if filters else ""

    saved = plot_histogram(di_vals, seq_vals, out, log_y=not args.linear)
    print(f"records: {n_read}, filtered out: {n_filtered}{fl}")
    print(f"di samples: {len(di_vals):,}, seq samples: {len(seq_vals):,}")
    print(f"saved: {saved}")

    gsda_tmp_data_dir = pathlib.Path("/root/projects/gsda/tmp-data-dir")
    if gsda_tmp_data_dir.exists():
        print(f"copy {saved} to {gsda_tmp_data_dir / saved}")
        shutil.copy(saved, gsda_tmp_data_dir / saved)
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
