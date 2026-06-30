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


def plot_true_base(
    df_base: pl.DataFrame,
    true_base: str,
    out_dir: Path,
    ncols: int,
    dpi: int,
):

    print(f"[INFO] Plotting {true_base}. rows:{df_base.shape}")

    motifs = sorted(
        df_base["motif"].unique().to_list(),
        key=motif_length,
    )

    print(f"motifs:{motifs}")
    if not motifs:
        return

    ncols = min(ncols, len(motifs))
    nrows = (len(motifs) + ncols - 1) // ncols

    figsize = (ncols * 7, nrows * 5)

    print(f"nrows={nrows}, ncols={ncols}, figsize:{figsize}")

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize,
        squeeze=False,
    )

    print("subplots")

    axes = axes.flatten()

    rng = np.random.default_rng(1234)

    for ax, motif in zip(axes, motifs):
        print("line: 101")
        sub = (
            df_base
            .filter(pl.col("motif") == motif)
            .sort("called")
        )
        print("line: 107")

        if len(sub) == 0:
            continue

        true_cnt = int(sub["true_cnt"][0])

        calleds = [true_cnt - 1, true_cnt, true_cnt + 1]

        data = []
        positions = []
        labels = []

        for pos, called in enumerate(calleds, start=1):
            print("line: 122")
            ratios = (
                sub
                .filter(pl.col("called") == called)
                .get_column("ratio_within_motif_tag")
                .to_list()
            )

            print(f"motif:{motif}, called:{called}, cnt:{len(ratios)}")

            # random.shuffle(ratios)
            # ratios = ratios[:5000]

            if not ratios:
                continue

            data.append(ratios)
            positions.append(pos)
            labels.append(str(called))

        if not data:
            continue

        vp = ax.violinplot(
            data,
            positions=positions,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        for body, pos in zip(vp["bodies"], positions):
            if pos == 1:
                body.set_alpha(0.5)
                body.set_color("tab:red")
            elif pos == 2:
                body.set_alpha(0.5)
                body.set_color("tab:blue")
            else:
                body.set_alpha(0.5)
                body.set_color("tab:green")

        for pos, values in zip(positions, data):
            jitter = rng.normal(0, 0.04, len(values))

            if pos == 1:
                color = "tab:red"
            elif pos == 2:
                color = "tab:blue"
            else:
                color = "tab:green"

            ax.scatter(
                np.full(len(values), pos) + jitter,
                values,
                s=8,
                alpha=0.6,
                color=color,
            )

        ax.set_title(motif)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)

        ax.set_xlabel("called")
        ax.set_ylabel("ratio_within_motif_tag")

        ax.set_ylim(0, 1)
        ax.grid(axis="y", alpha=0.3)

    for ax in axes[len(motifs):]:
        ax.axis("off")

    fig.suptitle(
        f"true_base = {true_base}",
        fontsize=18,
    )

    fig.tight_layout()

    outfile = out_dir / f"{true_base}.jpg"

    fig.savefig(
        outfile,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"[INFO] Saved {outfile}")


def plot_all_bases(
    df: pl.DataFrame,
    out_dir: Path,
    ncols: int,
    dpi: int,
    tag: str,
    cumulative=False,
):
    motif_infos = (
        df.select(["true_base", "motif"])
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

    ncols = min(ncols, nplots)
    nrows = (nplots + ncols - 1) // ncols

    figsize = (ncols * 7, nrows * 5)

    print(
        f"[INFO] nplots={nplots}, "
        f"nrows={nrows}, "
        f"ncols={ncols}"
    )

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)

    axes = axes.flatten()

    rng = np.random.default_rng(1234)

    for ax, info in zip(axes, motif_infos):

        true_base = info["true_base"]
        motif = info["motif"]

        sub = (
            df.filter(
                (pl.col("true_base") == true_base)
                & (pl.col("motif") == motif)
            )
            .sort("called")
        )

        if len(sub) == 0:
            continue

        true_cnt = int(sub["true_cnt"][0])

        calleds = [
            true_cnt - 1,
            true_cnt,
            true_cnt + 1,
        ]

        data = []
        positions = []
        labels = []

        ratio_field = "ratio" if cumulative else "ratio_within_motif_tag"

        for pos, called in enumerate(calleds, start=1):

            ratios = (
                sub.filter(pl.col("called") == called)
                .get_column(ratio_field)
                .to_list()
            )

            if not ratios:
                continue

            data.append(ratios)
            positions.append(pos)
            labels.append(str(called))

        if not data:
            continue

        vp = ax.violinplot(
            data,
            positions=positions,
            showmeans=True,
            showextrema=False,
            quantiles=[[0.1, 0.25, 0.5, 0.75, 0.9]] * len(data)
        )
        
        vp["cquantiles"].set_color("black")
        vp["cquantiles"].set_linewidth(2)
        vp["cquantiles"].set_linestyle("--")

        for body, pos in zip(vp["bodies"], positions):
            if pos == 1:
                body.set_color("tab:red")
            elif pos == 2:
                body.set_color("tab:blue")
            else:
                body.set_color("tab:green")

            body.set_alpha(0.5)

        for pos, values in zip(positions, data):

            jitter = rng.normal(
                0,
                0.04,
                len(values),
            )

            color = (
                "tab:red"
                if pos == 1
                else "tab:blue"
                if pos == 2
                else "tab:green"
            )

            ax.scatter(
                np.full(len(values), pos) + jitter,
                values,
                s=5,
                alpha=0.3,
                color=color,
            )

        ax.set_title(
            f"{true_base}:{motif}",
            fontsize=10,
        )

        ax.set_xticks(positions)
        ax.set_xticklabels(labels)

        ax.set_ylim(0, 1)

    for ax in axes[nplots:]:
        ax.axis("off")

    fig.suptitle(
        "Homopolymer Ratio Distribution",
        fontsize=20,
    )

    fig.tight_layout()

    tag = tag.replace("+", "And")
    cum_str = "cum" if cumulative else ""
    outfile = out_dir / f"all-bases-{tag}-{cum_str}.jpg"

    fig.savefig(
        outfile,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"[INFO] Saved {outfile}")


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
        for cumlative in [True, False]:
            plot_all_bases(df, output_dir, ncols=args.ncols,
                           dpi=args.dpi, tag=tag, cumulative=cumlative)

    print("[INFO] Done")


if __name__ == "__main__":
    main()
