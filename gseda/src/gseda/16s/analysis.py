import polars as pl
import os


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "300"
    os.environ["POLARS_FMT_STR_LEN"] = "100"


def main():
    # fname = "/data/ccs_data/qunfeng-2025/20250528_250301Y0001_Run0002_called_adapter.polish2.smc_all_reads-metric/20250528_250301Y0001_Run0002_called_adapter.polish2.smc_all_reads.gsmm2_aligned_metric_fact.csv"
    fname = "/data/ccs_data/qunfeng-2025/20250528_250301Y0001_Run0002_called_adapter.polish3.smc_all_reads-metric/20250528_250301Y0001_Run0002_called_adapter.polish3.smc_all_reads.gsmm2_aligned_metric_fact.csv"
    
    data = pl.read_csv(fname, separator="\t")
    data = data.filter(pl.col("rname").str.starts_with(
        "Akkermansia_muciniphila_16s"))
    
    data = data.filter(pl.col("rname").str.starts_with(
        "Akkermansia_muciniphila_16s_3"))
    
    print(len(data))
    print(data.group_by("rname").agg([pl.len()]))
    selected_data = data.select([
        pl.col("qname"), 
        pl.col("identity"),
        pl.col("ins-C"),
        pl.col("homoIns-C"),
        pl.col("del-C"),
        pl.col("homoDel-C"),
        pl.col("ins-G"),
        pl.col("homoIns-G"),
        pl.col("del-G"),
        pl.col("homoDel-G"),
        ])
    print(selected_data)
    print(selected_data.select([pl.col("identity").median()]))
    
    channels = selected_data.select([pl.col("qname").str.split("/").list.last()])
    channels.write_csv(file="/data/ccs_data/qunfeng-2025/expected_channels.csv", include_header=False)

    pass


if __name__ == "__main__":
    polars_env_init()
    main()
    pass
