import polars as pl
import polars_init


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

    print(df.head(100))


if __name__ == "__main__":
    polars_init.polars_env_init()

    fp = "/data/ccs_data/case-study/kangweishiji-q20-low-ratio/bc1-np3-smc2ref-analysis/fact_aligned_bam_ref_locus_info.csv"
    ana(filepath=fp)
    pass
