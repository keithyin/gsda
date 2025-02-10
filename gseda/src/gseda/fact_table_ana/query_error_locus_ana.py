import polars as pl
import argparse

import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)

import polars_init


def error_channel_analysis(
    fact_aligned_bam_bam_basic: str, fact_error_query_locus_info: str
):

    basic = pl.read_csv(fact_aligned_bam_bam_basic, separator="\t")
    query_error_locus = pl.read_csv(fact_error_query_locus_info, separator="\t")

    basic = basic.filter((pl.col("np") == 7).and_(pl.col("iy") < 0.999)).select(
        [
            pl.col("qname"),
            pl.col("np"),
            pl.col("rq"),
            pl.col("iy"),
            pl.col("qlen"),
            pl.col("fwd"),
        ]
    )

    df = (
        basic.join(query_error_locus, on="qname", how="inner")
        .with_columns(
            [
                pl.when(pl.col("fwd"))
                .then(pl.col("qstart"))
                .otherwise(pl.col("qlen") - pl.col("qend"))
                .alias("qstart"),
                pl.when(pl.col("fwd"))
                .then(pl.col("qend"))
                .otherwise(pl.col("qlen") - pl.col("qstart"))
                .alias("qend"),
            ]
        )
        .sort(by=["qname", "rstart"], descending=[True, False])
    )
    print(df.head(100))


def main(args):
    error_channel_analysis(args.bam_basic, args.error_query_locus)
    pass


if __name__ == "__main__":
    polars_init.polars_env_init()

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("bam_basic")
    parser.add_argument("error_query_locus")
    main(parser.parse_args())
