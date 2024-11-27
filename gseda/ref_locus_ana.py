import polars as pl
import polars_init
import argparse


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
        .filter(pl.col("curIsHomo").eq(0).and_(pl.col("nextIsHomo").eq(0)))
        .sort(by=["eq_rate2"], descending=[False])
    )

    print(df.head(200))

    print(df.select((pl.col("eq_rate2") < 0.5).sum() / pl.len()))

    # print(df.head(200))


if __name__ == "__main__":
    polars_init.polars_env_init()

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("fp", metavar="fact_aligned_bam_ref_locus_info.csv")
    args = parser.parse_args()
    ana(filepath=args.fp)
    pass
