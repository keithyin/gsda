import polars as pl
import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)

import polars_init
import argparse


def variant_calling(df: pl.DataFrame):
    v = (
        df.with_columns([(pl.col("eq") + pl.col("diff")).alias("eq_and_diff")])
        .with_columns(
            [
                (pl.col("eq") / pl.col("eq_and_diff")).alias("eq_in_eqdiff"),
                (pl.col("diff") / pl.col("eq_and_diff")).alias("diff_in_eqdiff"),
            ]
        )
        .filter((pl.col("eq_in_eqdiff") > 0.2).and_(pl.col("diff_in_eqdiff") > 0.2))
        .filter(pl.col("diffDetail").str.split(",").list.len() < 2)
        .shape[0]
    )

    print("variant calling ratio: {} / {} = {}".format(v, df.shape[0], v / df.shape[0]))


def del_calling(df: pl.DataFrame):
    v = (
        df.with_columns([(pl.col("del") / pl.col("depth")).alias("del_ratio")])
        .filter((pl.col("del_ratio") > 0.4))
        .shape[0]
    )

    print("del ratio: {} / {} = {}".format(v, df.shape[0], v / df.shape[0]))


def ana(filepath: str):
    df = pl.read_csv(filepath, separator="\t")
    print(df.head(2))

    df = (
        df.with_columns(
            [
                (
                    pl.col("eq")
                    / (pl.col("eq") + pl.col("diff") + pl.col("ins") + pl.col("del"))
                ).alias("eq_rate"),
                (pl.col("eq") / pl.col("depth")).alias("eq_rate2"),
            ]
        )
        # .filter(pl.col("curIsHomo").eq(0).and_(pl.col("nextIsHomo").eq(0)))
        .sort(by=["eq_rate2"], descending=[False])
    )

    print(df.head(100))

    print(df.select((pl.col("eq_rate2") < 0.5).sum() / pl.len()))

    variant_calling(df=df)
    del_calling(df=df)

    # print(df.head(200))


if __name__ == "__main__":
    polars_init.polars_env_init()

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("fp", metavar="fact_aligned_bam_ref_locus_info.csv")
    args = parser.parse_args()
    ana(filepath=args.fp)
    pass
