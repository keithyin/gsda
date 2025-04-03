import polars as pl
import math
import os
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def set_polars_env():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "100"
    os.environ["POLARS_FMT_STR_LEN"] = "50"


def phreq2q(v):
    return 1 - math.pow(10, v / -10)


def q2phreqExpr(name, new_name=None):
    if new_name is None:
        new_name = name
    return (
        -10
        * (
            1
            - pl.when(pl.col(name) > (1 - 1e-10))
            .then(1 - 1e-10)
            .otherwise(pl.col(name))
        ).log10()
    ).alias(new_name)


def AccuracyExpr(phreq_threshold):
    threshold = phreq2q(phreq_threshold)
    # threshold = phredq_threshold
    return (
        pl.col("rq").ge(threshold).and_(pl.col("iy").ge(threshold)).sum()
        / pl.col("rq").ge(threshold).sum()
    ).alias(f"≥Q{phreq_threshold}Accuracy")


def RecallExpr(phreq_threshold):
    threshold = phreq2q(phreq_threshold)
    # threshold = phredq_threshold
    return (
        pl.col("rq").ge(threshold).and_(pl.col("iy").ge(threshold)).sum()
        / pl.col("iy").ge(threshold).sum()
    ).alias(f"≥Q{phreq_threshold}Recall")


def plot_rq_iy_scatter(df: pl.DataFrame, o_prefix=None):
    figure = plt.figure(figsize=(10, 10))
    axs = figure.add_subplot(1, 1, 1)
    plt.sca(axs)
    plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    df = df.to_pandas()
    sns.scatterplot(df, x="phreq_rq", y="phreq_iy", ax=axs)
    sns.kdeplot(
        df,
        x="phreq_rq",
        y="phreq_iy",
        ax=axs,
        cmap="RdYlGn",
        fill=True,
        alpha=0.5,
        label="Density",
        cbar=True,
    )

    axs.set_xticks(list(range(0, 60, 2)))
    axs.set_yticks(list(range(0, 60, 2)))
    axs.set_xlabel("PredictedChannelQ", fontdict={"size": 16})
    axs.set_ylabel("EmpericalChannelQ", fontdict={"size": 16})
    perfect_line = pl.DataFrame(
        {
            "x": list(range(0, 60)),
            "y": list(range(0, 60)),
        }
    )

    sns.lineplot(
        perfect_line.to_pandas(), x="x", y="y", ax=axs, color="blue", linestyle="--"
    )
    fname = "rq_iy_scatter.png"
    if o_prefix is not None:
        fname = f"{o_prefix}-{fname}"
    figure.savefig(fname=fname)
    

def plot_predQ20_but_lower_Q20(df: pl.DataFrame, o_prefix=None):
    figure = plt.figure(figsize=(10, 10))
    axs = figure.add_subplot(1, 1, 1)
    plt.sca(axs)
    plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    df = df.filter(pl.col("phreq_rq") >= pl.lit(20))\
        .filter(pl.col("phreq_iy") < pl.lit(20))\
        .select([pl.col("phreq_iy")])\
        .to_pandas()
    
    axs.set_xlim((10, 20))
    axs.set_xticks(np.linspace(10, 20, num=21))
    
    sns.histplot(df, x="phreq_iy", ax=axs, bins=100)

    axs.set_xlabel("TrueQ", fontdict={"size": 16})
    axs.set_ylabel("Count", fontdict={"size": 16})
    plt.title("True Q-score distribution for predicted Q >= 20")
    fname = "pred_lt_q20_but_lower_q20.png"
    if o_prefix is not None:
        fname = f"{o_prefix}-{fname}"
    figure.savefig(fname=fname)
    
    pass


def stat(fname: str, o_prefix):
    df = pl.read_csv(fname, separator="\t")
    df = df.with_columns(
        [
            pl.when(pl.col("rq") > 0.99999)
            .then(0.99999)
            .otherwise(pl.col("rq"))
            .alias("rq"),
            pl.when(pl.col("iy") > 0.99999)
            .then(0.99999)
            .otherwise(pl.col("iy"))
            .alias("iy"),
        ]
    ).with_columns([q2phreqExpr("rq", "phreq_rq"), q2phreqExpr("iy", "phreq_iy")])

    plot_rq_iy_scatter(df=df, o_prefix=o_prefix)
    plot_predQ20_but_lower_Q20(df=df, o_prefix=o_prefix)

    stat_res = df.select(
        (pl.col("rq") - pl.col("iy")).mean().alias("ME(pred-real)(rq-iy)"),
        (pl.col("rq") - pl.col("iy")).median().alias("MedE"),
        (pl.col("rq") - pl.col("iy")).abs().mean().alias("MAE"),
        ((pl.col("rq") - pl.col("iy")).abs() / pl.col("iy")).mean().alias("MAPE"),
        AccuracyExpr(13),
        RecallExpr(13),
        AccuracyExpr(20),
        RecallExpr(20),
        AccuracyExpr(25),
        RecallExpr(25),
        AccuracyExpr(30),
        RecallExpr(30),
    )
    print(stat_res)

    # .with_columns([pl.col("phreq_rq").mul(1/3).cast(pl.Int32).mul(3)])\
    stats2 = (
        df.with_columns([pl.col("phreq_rq").cast(pl.Int32)])
        .group_by("phreq_rq")
        .agg(
            [
                (pl.col("rq") - pl.col("iy")).mean().alias("ME(pred-real)(rq-iy)"),
                (pl.col("rq") - pl.col("iy")).median().alias("MedE"),
                (pl.col("rq") - pl.col("iy")).abs().mean().alias("MAE"),
                ((pl.col("rq") - pl.col("iy")).abs() / pl.col("iy"))
                .mean()
                .alias("MAPE"),
            ]
        )
        .sort(["phreq_rq"], descending=[False])
    )
    print(stats2)

    stats3 = (
        df.with_columns(
            [
                q2phreqExpr("rq", "phreq_rq").cast(pl.Int32),
            ]
        )
        .group_by("phreq_rq")
        .agg(
            [
                pl.col("iy").mean().alias("avg"),
                pl.col("iy").median().alias("median"),
                pl.col("iy").std().alias("std"),
                pl.col("phreq_iy").min().alias("min"),
                pl.quantile("phreq_iy", quantile=0.05).alias("percent_5"),
                pl.quantile("phreq_iy", quantile=0.20).alias("percent_20"),
                pl.quantile("phreq_iy", quantile=0.50).alias("percent_50"),
                pl.quantile("phreq_iy", quantile=0.80).alias("percent_80"),
                pl.quantile("phreq_iy", quantile=0.95).alias("percent_95"),
                pl.col("phreq_iy").max().alias("max"),
                pl.len().alias("cnt")
            ]
        )
        .with_columns(
            [
                q2phreqExpr("avg"),
                q2phreqExpr("median"),
            ]
        )
        .select(
            pl.col("phreq_rq"),
            pl.col("avg"),
            pl.col("median"),
            pl.concat_str(
                pl.col("min").round(2),
                pl.col("percent_5").round(2),
                pl.col("percent_20").round(2),
                pl.col("percent_50").round(2),
                pl.col("percent_80").round(2),
                pl.col("percent_95").round(2),
                pl.col("max").round(2),
                separator=", ",
            ).alias("Percent_0_5_20_50_80_95_100"),
            pl.col("cnt")
        )
        .sort(["phreq_rq"], descending=[False])
    )
    print(stats3)


def main(args):
    stat(args.fact_table, args.o_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("rq_iy_diff")
    parser.add_argument("fact_table", metavar="fact_aligned_bam_bam_basic")
    parser.add_argument("--o-prefix", metavar="o-prefix", default=None, dest="o_prefix")

    args_ = parser.parse_args()

    set_polars_env()
    main(args_)
