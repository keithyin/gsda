import polars as pl
import math
import os
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import pathlib
import sys
cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))

print(f"cur_dir:{cur_dir}")
print(f"prev_dir:{prev_dir}")
print(f"prev_prev_dir:{prev_prev_dir}")
print(f"sys.path={sys.path}")
from fact_table_ana import pred_baseq_and_emp_q  # noqa: E402


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

    smc_path = pathlib.Path(args.smc_bam)
    smc_root = smc_path.parent
    smc_stem = smc_path.stem

    o_prefix = f"{smc_root}/{smc_stem}.aligned"

    gsmm2_cmd = f"gsmm2 align -q {args.smc_bam} -t {args.ref} --noMar -p {o_prefix}"
    if args.np_range is not None:
        gsmm2_cmd += f" --np-range {args.np_range}"
    if args.rq_range is not None:
        gsmm2_cmd += f" --rq-range {args.rq_range}"
        

    print(f"running {gsmm2_cmd} ")
    subprocess.check_call(gsmm2_cmd, shell=True)

    gsetl_o_dir = f"{o_prefix}-gsetl"
    gsda_cmd = f"gsetl --outdir {gsetl_o_dir} aligned-bam --bam {o_prefix}.bam --ref-file {args.ref} --factRecordStat 0 --factRefLocusInfo 0 --factErrorQueryLocusInfo 0 --factBaseQStat 1 --factPolyInfo 0"
    print(f"running {gsda_cmd}")
    subprocess.check_call(gsda_cmd, shell=True)

    fact_table_path = f"{gsetl_o_dir}/fact_aligned_bam_bam_basic.csv"
    stat_o_prefix = f"{gsetl_o_dir}/{smc_stem}"

    print(f"CHANNEL_Q_ANALYSIS")
    stat(fact_table_path, stat_o_prefix)

    print("BASE_Q_ANALYSIS")

    baseq_ana_args = {
        "fact_table": f"{gsetl_o_dir}/fact_baseq_stat.csv",
        "o_prefix": None
    }
    pred_baseq_and_emp_q.main(argparse.Namespace(**baseq_ana_args))


if __name__ == "__main__":
    parser = argparse.ArgumentParser("rq_iy_diff")
    parser.add_argument("--smc-bam", metavar="smc_bam", dest="smc_bam")
    parser.add_argument("--ref", metavar="ref.fasta")
    parser.add_argument("--np-range", default=None, type=str, dest="np_range", help="1:3,5,7:9 means [[1, 3], [5, 5], [7, 9]]. only valid for bam input that contains np field")
    parser.add_argument("--rq-range", default=None, type=str, dest="rq_range", help="0.9:1.1 means 0.9<=rq<=1.1. only valid for bam input that contains rq field")

    args_ = parser.parse_args()

    set_polars_env()
    main(args_)
