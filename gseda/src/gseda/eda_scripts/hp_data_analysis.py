import polars as pl
import os

def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "3000"
    os.environ["POLARS_FMT_STR_LEN"] = "100"

def main():
    df = pl.read_csv(
        source="/data/ccs_data/speed-test/1M-1h-ecoli/smicing.smc_all_reads-metric/smicing.smc_all_reads.gsmm2-hp-fact.csv", separator="\t")
    df = df.filter(pl.col("tag") == "mixed").with_columns(
        [pl.col("qname").str.extract(r"/(\d+)", group_index=1).alias("channel")])
    
    wanted_channels = df.select([pl.col("qname"), pl.col("channel")]).unique().sort(
        by=[pl.col("channel")], descending=[False]).head(20)
    
    joined_result = df.join(wanted_channels, on="qname", how="inner")
    print(joined_result)
    
    wanted_channels = ",".join(wanted_channels.to_pandas()["channel"].tolist())

    print(wanted_channels)

    pass


if __name__ == "__main__":
    polars_env_init()
    main()
