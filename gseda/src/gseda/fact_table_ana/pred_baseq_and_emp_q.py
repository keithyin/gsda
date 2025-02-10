import pysam
import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)

import utils
import polars as pl
import argparse
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import argparse


import polars_init


def main(args):
    polars_init.polars_env_init()

    # plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    df = pl.read_csv(args.fact_table, separator="\t")
    df_shift = df.with_columns([pl.col("baseq")])
    df = df.with_columns(
        [
            (pl.col("eq") / (pl.col("eq") + pl.col("diff") + pl.col("ins"))).alias(
                "emp_rq"
            )
        ]
    ).with_columns([utils.q2phreq_expr("emp_rq", "emp_phreq")])
    figure = plt.figure(figsize=(10, 10))
    axs = figure.add_subplot(1, 1, 1)
    plt.sca(axs)
    plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    sns.scatterplot(df.to_pandas(), x="baseq", y="emp_phreq", ax=axs)
    axs.set_xticks(list(range(0, 60, 2)))
    axs.set_yticks(list(range(0, 60, 2)))
    axs.set_xlabel("PredictedBaseQ", fontdict={"size": 16})
    axs.set_ylabel("EmpericalBaseQ", fontdict={"size": 16})
    perfect_line = pl.DataFrame(
        {
            "x": list(range(0, 60)),
            "y": list(range(0, 60)),
        }
    )

    sns.lineplot(
        perfect_line.to_pandas(), x="x", y="y", ax=axs, color="blue", linestyle="--"
    )

    print(df.head(60))
    fname = "baseq2empq.png"
    if args.o_prefix is not None:
        fname = f"{args.o_prefix}-{fname}"
    figure.savefig(fname=fname)
    print(f"check image {fname}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="",
        usage="""
    gsmm2 align -q $query_file -t $ref_file -p outputbam_prefix --noMar
    gsetl --outdir $outdir aligned-bam --bam $aligned_bam --ref-file $ref_file
    python preq-baseq-and-emp-q.py $outdir/fact_baseq_stat.csv
""",
    )
    parser.add_argument("fact_table", metavar="fact_baseq_stat")
    parser.add_argument("--o-prefix", metavar="o-prefix", default=None, dest="o_prefix")

    main(parser.parse_args())
