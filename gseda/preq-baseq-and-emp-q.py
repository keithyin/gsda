import pysam
import utils
import polars as pl
import argparse
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from numba import jit
import polars_init


def main(args):
    # plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    df = pl.read_csv(args.data, separator="\t")
    df = df.with_columns(
        [
            (pl.col("eq") / (pl.col("eq") + pl.col("diff") + pl.col("ins"))).alias(
                "emp_rq"
            )
        ]
    ).with_columns([utils.q2phreq_expr("emp_rq", "emp_phreq")])
    figure = plt.figure(figsize=(20, 10))
    axs = figure.add_subplot(1, 1, 1)
    plt.sca(axs)
    plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    sns.scatterplot(df.to_pandas(), x="baseq", y="emp_phreq", ax=axs)
    perfect_line = pl.DataFrame(
        {
            "x": list(range(0, 60)),
            "y": list(range(0, 60)),
        }
    )

    sns.lineplot(
        perfect_line.to_pandas(), x="x", y="y", ax=axs, color="blue", linestyle="--"
    )

    print(df.head(10))
    figure.savefig(fname="baseq2empq.png")


if __name__ == "__main__":
    polars_init.polars_env_init()
    params = {
        "data": "/data/ccs_data/ccs_eval2024q3/Ludaopei/subread_bak/output-all/analysis/fact_baseq_stat.csv",
    }

    main(argparse.Namespace(**params))
