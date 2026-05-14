import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import logging
from typing import List

import polars as pl
import matplotlib

matplotlib.use("Agg")


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)

DEFAULT_METRIC_DIR = "/data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/channel_subreads"


def scan_aggr_csvs(metric_dir: str) -> pl.DataFrame:
    """
    Scan all *-hp-aggr.csv files under metric_dir and concat them.
    """
    dfs = []
    count = 0

    for dirpath, _, filenames in os.walk(metric_dir):
        for fname in filenames:

            if not fname.endswith("-hp-aggr.csv"):
                continue

            fpath = os.path.join(dirpath, fname)

            try:
                df = pl.read_csv(fpath, separator="\t")
                dfs.append(df)
                count += 1

            except Exception as e:
                logging.warning("Failed reading %s: %s", fpath, e)

    if not dfs:
        logging.warning("No aggr CSV files found under %s", metric_dir)
        return pl.DataFrame()

    logging.info("Scanned %d aggr CSV files", count)

    return pl.concat(dfs, how="vertical_relaxed")


def sample_values(vals, max_points=500):
    """
    Random sample for scatter overlay.
    """
    vals = np.asarray(vals)

    if len(vals) <= max_points:
        return vals

    idx = np.random.choice(
        len(vals),
        size=max_points,
        replace=False,
    )

    return vals[idx]


def plot_ratio_distribution_violin(
    metric_dir: str,
    motifs: List[str],
    outpath: str = "ratio_distribution_violin.png",
    tag: str = "pure+mixed",
    abs_max: int = 2,
    true_cnt_min: int = 1,
    scatter_sample_size: int = 500,
):
    """
    Violin plot + vertical scatter overlay.

    Scatter 不使用 jitter。
    """

    df = scan_aggr_csvs(metric_dir)

    if df.is_empty():
        logging.warning("Empty dataframe")
        return outpath

    # ==========================================================
    # Filter
    # ==========================================================
    df = df.filter(
        pl.col("motif").is_in(motifs)
        & (pl.col("tag") == tag)
        & ((pl.col("called") - pl.col("true_cnt")).abs() <= abs_max)
        & (pl.col("true_cnt") >= true_cnt_min)
    )

    total_rows = len(df)

    logging.info("Rows after filtering: %d", total_rows)

    if total_rows == 0:
        logging.warning("No rows after filtering")
        return outpath

    # ==========================================================
    # Group
    # ==========================================================
    # ratio_within_motif_tag
    grouped = (
        df.group_by(["motif", "true_cnt", "called"])
        .agg(
            [
                pl.col("ratio").alias("ratios"),
            ]
        )
        .sort(["motif", "true_cnt", "called"])
    )

    motif_list = (
        grouped.select("motif")
        .unique()
        .get_column("motif")
        .to_list()
    )

    motif_list = sorted(motif_list)
    n = len(motif_list)

    cols = min(3, n)
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(cols * 7, rows * 5),
        squeeze=False,
    )

    axes = axes.flatten()

    # ==========================================================
    # Per motif
    # ==========================================================
    for ax_idx, motif in enumerate(motif_list):

        ax = axes[ax_idx]

        sub = (
            grouped.filter(pl.col("motif") == motif)
            .sort(["true_cnt", "called"])
        )

        if sub.is_empty():
            ax.set_visible(False)
            continue

        x_positions = list(range(len(sub)))

        x_labels = [
            f"{t}/{c}"
            for t, c in zip(
                sub["true_cnt"].to_list(),
                sub["called"].to_list(),
            )
        ]

        all_ratios = sub["ratios"].to_list()

        all_ratios_np = [
            np.asarray(r.to_list() if hasattr(r, "to_list") else r)
            for r in all_ratios
        ]

        valid = [
            (x, lbl, arr)
            for x, lbl, arr in zip(
                x_positions,
                x_labels,
                all_ratios_np,
            )
            if len(arr) > 0
        ]

        if not valid:
            ax.set_visible(False)
            continue

        x_positions, x_labels, all_ratios_np = zip(*valid)

        # ======================================================
        # Violin plot
        # ======================================================
        vp = ax.violinplot(
            all_ratios_np,
            positions=x_positions,
            widths=0.8,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        for body in vp["bodies"]:
            body.set_alpha(0.35)

        if "cmedians" in vp:
            vp["cmedians"].set_linewidth(1.5)

        # ======================================================
        # Scatter overlay (NO jitter)
        # ======================================================
        median_list = []

        for pos, vals in zip(
            x_positions,
            all_ratios_np,
        ):

            vals = np.asarray(vals)

            if len(vals) == 0:
                median_list.append(np.nan)
                continue

            med = np.median(vals)

            median_list.append(med)

            vals_sample = sample_values(
                vals,
                max_points=scatter_sample_size,
            )

            x_vals = np.full(
                len(vals_sample),
                pos,
            )

            ax.scatter(
                x_vals,
                vals_sample,
                s=6,
                alpha=0.25,
                linewidths=0,
            )

        # ======================================================
        # Annotate
        # ======================================================
        for pos, vals, med in zip(
            x_positions,
            all_ratios_np,
            median_list,
        ):

            n_pts = len(vals)

            ax.text(
                pos,
                med + 0.015,
                f"m={med:.3f}\nn={n_pts}",
                ha="center",
                va="bottom",
                fontsize=6,
            )

        # ======================================================
        # Axis style
        # ======================================================
        ax.set_xticks(x_positions)

        ax.set_xticklabels(
            x_labels,
            fontsize=7,
            rotation=45,
            ha="right",
        )

        ax.set_title(
            motif,
            fontsize=11,
        )

        ax.set_xlabel("true_cnt / called")

        ax.set_ylabel("ratio_within_motif_tag")

        ax.grid(axis="y", alpha=0.3)

    # Hide unused axes
    for ax_idx in range(n, len(axes)):
        axes[ax_idx].set_visible(False)

    fig.suptitle(
        (
            f"ratio_within_motif_tag distribution "
            f"(rows={total_rows}, motifs={motifs})"
        ),
        fontsize=12,
    )

    fig.tight_layout()

    fig.savefig(
        outpath,
        dpi=180,
        bbox_inches="tight",
    )

    plt.close(fig)

    logging.info("Saved plot to %s", outpath)

    return outpath


def main_cli():

    parser = argparse.ArgumentParser(
        prog="metric_ratio_distribution_violin"
    )

    parser.add_argument(
        "metric_dir",
        nargs="?",
        default=None,
        help="channel_subreads directory",
    )

    parser.add_argument(
        "--motifs",
        nargs="+",
        default=None,
        help="motifs to analyze",
    )

    parser.add_argument(
        "--tag",
        default="pure+mixed",
        help="tag filter",
    )

    parser.add_argument(
        "--abs-max",
        type=int,
        default=2,
        help="abs(true_cnt-called) <= abs_max",
    )

    parser.add_argument(
        "--true-cnt-min",
        type=int,
        default=1,
        help="minimum true_cnt",
    )

    parser.add_argument(
        "--scatter-sample-size",
        type=int,
        default=10000000000000,
        help="max scatter points per group",
    )

    parser.add_argument(
        "--outpath",
        default="ratio_distribution_violin.png",
        help="output image path",
    )

    args = parser.parse_args()

    metric_dir = args.metric_dir or DEFAULT_METRIC_DIR

    motifs = args.motifs or ["(C)5", "(C)4"]

    plot_ratio_distribution_violin(
        metric_dir=metric_dir,
        motifs=motifs,
        outpath=args.outpath,
        tag=args.tag,
        abs_max=args.abs_max,
        true_cnt_min=args.true_cnt_min,
        scatter_sample_size=args.scatter_sample_size,
    )


if __name__ == "__main__":
    main_cli()
