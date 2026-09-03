import polars as pl
import os


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "300"
    os.environ["POLARS_FMT_STR_LEN"] = "10000"


def main():
    fname = "/data/gsmm2-debug-data/low-identity-ecoli/20250801_250302Y0004_Run0007_called.valid-metric/20250801_250302Y0004_Run0007_called.valid.gsmm2_aligned_metric_fact.csv"

    data = pl.read_csv(fname, separator="\t")
    print(len(data))
    data = data.with_columns([(pl.col("match") + pl.col("misMatch") + pl.col(
        "ins") + pl.col("homoIns") + pl.col("del") + pl.col("homoDel")).alias("alignedSpan")]).with_columns(
            [
                (pl.col("match") / pl.col("alignedSpan")).alias("identity"),
                (pl.col("ins") / pl.col("alignedSpan")).alias("NHinsRate"),
             ]
            )
    data = data.filter(pl.col("identity").gt(0.88).and_(pl.col("NHinsRate").gt(0.05)))
    data = data.sort(by="qname", descending=False)
    

    print(data.head(10).select([pl.col("qname"), pl.col("oriAlignInfo"), pl.col("identity"), pl.col("NHinsRate")]))
    align_info_txts = data.head(10).select([pl.col("oriAlignInfo")]).to_pandas()["oriAlignInfo"].tolist()
    align_info_txts = map(lambda txt: txt.replace(";", ";\n"), align_info_txts)
    align_info_txt = "\n\n".join(align_info_txts)
    print(align_info_txt)
    print(len(data))
    

    pass


if __name__ == "__main__":
    polars_env_init()
    main()
    pass
