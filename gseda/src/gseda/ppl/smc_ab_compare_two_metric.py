
import polars as pl
import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
print(cur_dir)
sys.path.append(cur_dir)
import env_prepare  # noqa: E402


def main():

    env_prepare.polars_env_init()

    exp_a_metric = "/data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter_metric/20260311_240601Y0012_Run0003.newmodel.stage1.baseline.smc_all_reads.aligned.metric.csv"
    exp_b_metric = "/data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter_metric/20260311_240601Y0012_Run0003.newmodel.stage1.norna-nofwdrev.noleftalign.smc_all_reads.aligned.metric.csv"

    df_a = pl.read_csv(exp_a_metric, separator="\t")
    df_b = pl.read_csv(exp_b_metric, separator="\t")

    df_a = df_a.filter(pl.col("subreadPasses") > 7)
    df_b = df_b.filter(pl.col("subreadPasses") > 7)

    df = df_a.join(df_b, on="channel_id", how="inner").filter(
        pl.col("concordance") > pl.col("concordance_right")).filter(
            pl.col("refStart") == pl.col("refStart_right")).filter(
                pl.col("refEnd") == pl.col("refEnd_right"))

    print(df_a.head(10))
    print(df_b.head(10))

    print(df.head(10))


if __name__ == "__main__":
    main()
