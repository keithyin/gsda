#!/usr/bin/env python3

from pathlib import Path
import argparse
import math
import re
import os

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib
from tqdm import tqdm
from multiprocessing import Pool
import random
matplotlib.use("Agg")


MOTIF_RE = re.compile(r"^\(([ACGT])\)(\d+)$")


def motif_length(motif: str) -> int:
    m = MOTIF_RE.match(motif)
    return int(m.group(2)) if m else 999999


def load_csv(filepath):
    filepath = Path(filepath)
    df = (
        pl.read_csv(filepath, separator="\t")
        .with_columns(pl.lit(filepath.stem).alias("sample"))
        .filter(
            (pl.col("called") >= pl.col("true_cnt") - 1)
            & (pl.col("called") <= pl.col("true_cnt") + 1)
        )
    )
    return df


def load_csv_wrapper(params):
    filepath, tag = params
    return load_csv(filepath, tag)


def load_all_csv(input_dir: Path) -> pl.DataFrame:
    dfs = []

    all_csv_files = []

    for metric_dir in tqdm(input_dir.rglob("*filtered-metric"), desc="walking ... "):
        if not metric_dir.is_dir():
            continue
        for csv_file in metric_dir.glob("*.gsmm2-hp-aggr.csv"):
            all_csv_files.append(csv_file)

    print(f"NumCsvFiles:{len(all_csv_files)}")
    with Pool(processes=10) as pool:
        dfs = []
        for df in tqdm(pool.imap_unordered(load_csv, all_csv_files, chunksize=1), desc=f"loading ...  "):
            dfs.append(df)

    if not dfs:
        raise RuntimeError(
            f"No *.gsmm2-hp-aggr.csv found under {input_dir}"
        )

    return pl.concat(dfs, how="vertical_relaxed")


def plot_delta_distribution(
    df: pl.DataFrame,
    out_dir: Path,
    dpi: int,
    tag: str,
    ncols: int,
    quantiles=[5, 10, 25, 50, 75, 90, 95],
):
    """
    delta =
        ratio_within_motif_tag(called=true_cnt-1)
        -
        ratio_within_motif_tag(called=true_cnt)

    Histogram + Quantile Lines
    """

    delta_df = (
        df
        .filter(
            (pl.col("called") == pl.col("true_cnt"))
            |
            (pl.col("called") == pl.col("true_cnt") - 1)
        )
        .with_columns(
            pl.when(
                pl.col("called") == pl.col("true_cnt") - 1
            )
            .then(pl.lit("minus1"))
            .otherwise(pl.lit("true"))
            .alias("called_type")
        )
        .pivot(
            values="ratio_within_motif_tag",
            index=[
                "sample",
                "motif",
                "true_base",
            ],
            on="called_type",
        )
        .with_columns(
            (
                pl.col("minus1")
                - pl.col("true")
            ).alias("delta")
        )
        .drop_nulls("delta")
    )

    if len(delta_df) == 0:
        print(f"[WARN] tag={tag}: no delta data")
        return

    motif_infos = (
        delta_df
        .select(
            [
                "true_base",
                "motif",
            ]
        )
        .unique()
        .to_dicts()
    )

    motif_infos = sorted(
        motif_infos,
        key=lambda x: (
            x["true_base"],
            motif_length(x["motif"]),
        ),
    )

    nplots = len(motif_infos)

    ncols = min(
        ncols,
        nplots,
    )

    nrows = (
        nplots + ncols - 1
    ) // ncols

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(
            ncols * 5,
            nrows * 4,
        ),
        squeeze=False,
        sharex=True,
    )

    axes = axes.flatten()

    bins = np.linspace(
        -1,
        1,
        21,
    )

    #
    # quantile颜色
    #
    color_cycle = [
        "gray",
        "royalblue",
        "red",
        "royalblue",
        "gray",
    ]

    for ax, info in zip(
        axes,
        motif_infos,
    ):

        true_base = info["true_base"]
        motif = info["motif"]

        values = (
            delta_df
            .filter(
                (pl.col("true_base") == true_base)
                &
                (pl.col("motif") == motif)
            )
            .get_column("delta")
            .to_numpy()
        )

        if len(values) == 0:
            ax.axis("off")
            continue

        #
        # histogram
        #
        ax.hist(
            values,
            bins=bins,
            density=True,
            alpha=0.7,
            color="tab:blue",
            edgecolor="black",
            linewidth=0.5,
        )

        #
        # quantiles
        #
        qvals = np.percentile(
            values,
            quantiles,
        )

        legend_handles = []

        for idx, (q, qv) in enumerate(
            zip(quantiles, qvals)
        ):

            if idx < len(color_cycle):
                color = color_cycle[idx]
            else:
                color = f"C{idx}"

            line = ax.axvline(
                qv,
                color=color,
                linewidth=2,
                alpha=0.9,
                linestyle="--",
            )

            legend_handles.append(
                (
                    line,
                    f"P{q}={qv:.3f}",
                )
            )

        #
        # delta=0参考线
        #
        ax.axvline(
            0,
            color="black",
            linestyle=":",
            linewidth=1.5,
        )

        ax.set_xlim(
            -1,
            1,
        )

        ax.grid(
            axis="y",
            alpha=0.3,
        )

        median = np.median(values)

        ax.set_title(
            (
                f"{true_base}:{motif}\n"
                f"n={len(values)}\n"
                f"median={median:.3f}"
            ),
            fontsize=9,
        )

        #
        # legend
        #
        ax.legend(
            [h[0] for h in legend_handles],
            [h[1] for h in legend_handles],
            fontsize=7,
            loc="upper right",
            framealpha=0.8,
        )

    #
    # 多余subplot关闭
    #
    for ax in axes[nplots:]:
        ax.axis("off")

    fig.supxlabel(
        "Δ = ratio(true_cnt-1) - ratio(true_cnt)"
    )

    fig.supylabel(
        "Density"
    )

    fig.suptitle(
        (
            f"{tag}\n"
            "Histogram of Delta Distribution"
        ),
        fontsize=18,
    )

    fig.tight_layout()

    outfile = (
        out_dir
        / (
            "delta-histogram-"
            f"{tag.replace('+', 'And')}.jpg"
        )
    )

    fig.savefig(
        outfile,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"[INFO] Saved {outfile}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot homopolymer ratio distributions"
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Root directory",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory",
    )

    parser.add_argument(
        "--ncols",
        type=int,
        default=4,
        help="Number of subplot columns",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image DPI",
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("[INFO] Loading data...")
    tags = ["pure", "mixed", "pure+mixed"]
    all_df = load_all_csv(input_dir)
    for tag in tags:
        df = all_df.filter(pl.col("tag") == pl.lit(tag))
        print(f"[INFO] tag={tag} Rows: {len(df):,}")
        true_bases = (
            df["true_base"]
            .unique()
            .sort()
            .to_list()
        )

        print(f"[INFO] Bases: {true_bases}")
        #  ncols=args.ncols
        # plot_delta_distribution(df, output_dir,
        #                         dpi=args.dpi, tag=tag)

        # plot_delta_distribution_v2(df, output_dir,
        # dpi=args.dpi, tag=tag, ncols=args.ncols)

        plot_delta_distribution(df, output_dir,
                                   dpi=args.dpi, tag=tag, ncols=args.ncols)

    print("[INFO] Done")


if __name__ == "__main__":
    main()
