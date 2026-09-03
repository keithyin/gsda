import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import polars as pl
import os
import sys
import numpy as np

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)
import polars_init  # noqa: E402
import utils  # noqa: E402


def _cdf(cnt: np.ndarray) -> np.ndarray:
    total = cnt.sum()
    if total == 0:
        return np.zeros_like(cnt, dtype=float)
    return np.cumsum(cnt) / total


def plot_ins_vs_non_ins_baseq_dist(df: pl.DataFrame, o_prefix: str):
    """Distribution of predicted baseQ for insertion bases vs non-insertion bases,
    drawn as a shared histogram plus one CDF curve per category."""
    ins_dist = (
        df.group_by("baseq").agg(pl.col("ins").sum().alias("cnt")).sort("baseq"))
    nonins_dist = (
        df.group_by("baseq").agg(
            (pl.col("eq") + pl.col("diff")).sum().alias("cnt")).sort("baseq"))

    ins_bq = ins_dist["baseq"].to_list()
    nonins_bq = nonins_dist["baseq"].to_list()
    max_baseq = max(ins_bq + nonins_bq + [0])
    baseq_axis = np.arange(0, max_baseq + 1, dtype=float)

    def to_axis(d, col):
        m = dict(zip(d["baseq"].to_list(), d[col].to_list()))
        return np.array([m.get(b, 0) for b in baseq_axis], dtype=float)

    ins_cnt = to_axis(ins_dist, "cnt")
    nonins_cnt = to_axis(nonins_dist, "cnt")

    # Normalize each group to its own fraction, so group-size differences
    # don't distort the visual comparison of baseQ distributions.
    ins_total = ins_cnt.sum()
    nonins_total = nonins_cnt.sum()
    ins_prob = ins_cnt / ins_total if ins_total > 0 else ins_cnt
    nonins_prob = nonins_cnt / nonins_total if nonins_total > 0 else nonins_cnt

    figure = plt.figure(figsize=(10, 10))
    axs = figure.add_subplot(1, 1, 1)
    plt.sca(axs)
    plt.grid(True, linestyle=":", linewidth=0.5, color="gray")

    width = 0.8
    half = width / 2
    axs.bar(
        baseq_axis - half / 2, ins_prob, half,
        label="Insertion base", color="#d62728")
    axs.bar(
        baseq_axis + half / 2, nonins_prob, half,
        label="Non-insertion base", color="#1f77b4")

    axs.set_xlim((0, max_baseq + 1))
    axs.set_xticks(np.arange(0, max_baseq + 1, 2))
    axs.set_xlabel("PredictedBaseQ", fontdict={"size": 16})
    axs.set_ylabel("Proportion", fontdict={"size": 16})

    axs2 = axs.twinx()
    axs2.plot(
        baseq_axis, _cdf(ins_cnt), color="#d62728",
        linestyle="--", label="Insertion CDF")
    axs2.plot(
        baseq_axis, _cdf(nonins_cnt), color="#1f77b4",
        linestyle="--", label="Non-insertion CDF")
    axs2.set_ylim((0, 1.02))
    axs2.set_ylabel("CDF", fontdict={"size": 16})

    h1, l1 = axs.get_legend_handles_labels()
    h2, l2 = axs2.get_legend_handles_labels()
    axs.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=12)

    plt.title("BaseQ distribution: insertion vs non-insertion", fontdict={"size": 16})

    fpath = f"{o_prefix}.ins-vs-nonins-baseq-dist.png"
    figure.savefig(fname=fpath)
    print(f"check image {fpath}")
    return fpath


def main(args):
    polars_init.polars_env_init()

    # plt.grid(True, linestyle=":", linewidth=0.5, color="gray")
    fact_table_path = pathlib.Path(args.fact_table)
    df = pl.read_csv(args.fact_table, separator="\t")

    if args.union_fwd_rev:
        df = df.with_columns(
            pl.col("refname").str.replace(r"___fwd$|___rev$", ""))

        df = df.group_by(["refname", "baseq"]).agg(
            [pl.col("eq").sum(), pl.col("diff").sum(), pl.col(
                "ins").sum(), pl.col("del").sum(), pl.col("depth").sum()]
        ).sort(by=["refname", "baseq"], descending=[False, False])

    df = df.with_columns(
        [
            (pl.col("eq") / (pl.col("eq") + pl.col("diff") + pl.col("ins") + pl.col("del"))).alias(
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

    print(df)

    summary = df.select([
        (pl.col("depth").filter(pl.col("baseq") >= 20).sum() /
         pl.col("depth").sum()).alias("baseq20Ratio"),
        (pl.col("depth").filter(pl.col("baseq") >= 30).sum() /
         pl.col("depth").sum()).alias("baseq30Ratio"),
        (pl.col("depth").filter(pl.col("baseq") >= 35).sum() /
         pl.col("depth").sum()).alias("baseq35Ratio"),
        (pl.col("depth").filter(pl.col("baseq") >= 40).sum() /
         pl.col("depth").sum()).alias("baseq40Ratio"),
    ])
    print("Base Q Summary:\n", summary.transpose(
        include_header=True, header_name="name", column_names=["value"]
    ))

    metric = (df.filter(pl.col("depth") > 10000)
              .select([(pl.col("baseq") - pl.col("emp_phreq")).pow(2.0).alias("SquareError"),
                       (pl.col("baseq") - pl.col("emp_phreq")
                        ).abs().alias("AbsError"),
                       ((pl.col("emp_phreq") - pl.col("baseq")).abs() /
                        pl.col("emp_phreq")).alias("ape"),
                       (2 * (pl.col("baseq")-pl.col("emp_phreq")).abs() /
                        (pl.col("baseq")+pl.col("emp_phreq"))).alias("sape")
                       ])
              .select([
                  pl.col("SquareError").mean().alias("MSE"),
                  pl.col("SquareError").mean().sqrt().alias("RMSE"),
                  pl.col("AbsError").mean().alias("MAE"),
                  pl.col("AbsError").median().alias("MedAE"),
                  pl.col("ape").mean().alias("MAPE"),
                  pl.col("sape").mean().alias("sMAPE"),

              ])
              )

    metric = metric.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )
    print("Base Q Metric: \n", metric)

    baseq2emp_baseq_fpath = f"{args.o_prefix}.baseq2empq.png"
    figure.savefig(fname=baseq2emp_baseq_fpath)
    print(f"check image {baseq2emp_baseq_fpath}")

    baseq_cnt = (df.filter(pl.col("depth") > 10000)
                 .select([pl.col("baseq"), pl.col("depth")]))
    baseq_cnt = baseq_cnt.to_pandas()

    # TODO baseq distribution
    figure = plt.figure(figsize=(10, 10))
    axs = figure.add_subplot(1, 1, 1)
    plt.sca(axs)
    plt.grid(True, linestyle=":", linewidth=0.5, color="gray")
    sns.barplot(baseq_cnt, x="baseq", y="depth",
                ax=axs, order=list(range(0, 50)),)

    axs.set_xticks(list(range(0, 50, 2)))
    axs.set_xlabel("PredictedBaseQ", fontdict={"size": 16})
    axs.set_ylabel("Count", fontdict={"size": 16})

    baseq_dist_fpath = f"{args.o_prefix}.baseq-dist.png"
    figure.savefig(fname=baseq_dist_fpath)
    print(f"check image {baseq_dist_fpath}")

    ins_vs_nonins_fpath = plot_ins_vs_non_ins_baseq_dist(
        df=df, o_prefix=args.o_prefix)

    return [baseq2emp_baseq_fpath, baseq_dist_fpath, ins_vs_nonins_fpath]


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
    parser.add_argument("--o-path", metavar="o-path",
                        default=None, dest="o_path")
    parser.add_argument("--union-fwd-rev", action="store_true", default=False)

    main(parser.parse_args())
