import argparse
import os

import polars as pl


def set_polars_env():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "100"
    os.environ["POLARS_FMT_STR_LEN"] = "50"


def stat(fact_table: str):
    df = pl.read_csv(fact_table, separator="\t")
    quantiles = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00]
    res = df.select(
        [pl.col("iy").quantile(q).alias(f"p{int(q * 100)}") for q in quantiles]
    )
    print(f"rows: {df.height}  iy min: {df['iy'].min():.6f}  iy max: {df['iy'].max():.6f}")
    print(res)

    res_by_np = (
        df.group_by("np")
        .agg(
            [pl.col("iy").quantile(q).alias(f"p{int(q * 100)}") for q in quantiles]
            + [pl.len().alias("cnt")]
        )
        .sort("np")
    )
    print("\nby np:")
    print(res_by_np)


def main(args):
    stat(args.fact_table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("fact_aligned_bam_basic_ana")
    parser.add_argument(
        "fact_table",
        metavar="fact_aligned_bam_bam_basic.csv",
        nargs="?",
        default="/data1/ccs_data/20260717-rna-report/20260731_240601Y0004_Run0006_low_concentration_c/smc-rna-1k-3h-gsetl/fact_aligned_bam_bam_basic.csv",
    )
    args_ = parser.parse_args()

    set_polars_env()
    main(args_)
