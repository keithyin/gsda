import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import polars as pl
import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)
import polars_init  # noqa: E402
import utils  # noqa: E402


def main(args):
    polars_init.polars_env_init()

    # plt.grid(True, linestyle=":", linewidth=0.5, color="gray")
    fact_table_path = pathlib.Path(args.fact_table)
    out_dir = fact_table_path.parent

    df = pl.read_csv(args.fact_table, separator="\t")
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

    fname = "baseq2empq.png"
    fpath = out_dir.joinpath(fname)
    print(df.head(60))

    summary = df.select([
        (pl.col("depth").filter(pl.col("baseq") >= 20).sum() /
         pl.col("depth").sum()).alias("baseq20Ratio"),
        (pl.col("depth").filter(pl.col("baseq") >= 30).sum() /
         pl.col("depth").sum()).alias("baseq30Ratio"),
        (pl.col("depth").filter(pl.col("baseq") >= 35).sum() /
         pl.col("depth").sum()).alias("baseq35Ratio"),
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
    

    if args.o_path is not None:
        fpath = args.o_path

    figure.savefig(fname=fpath)
    print(f"check image {fpath}")


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

    main(parser.parse_args())
