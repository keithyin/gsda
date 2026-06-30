import polars as pl

exp_a_metric = "/data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter_metric/20260311_240601Y0012_Run0003.newmodel.stage1.baseline.smc_all_reads.aligned.metric.csv"
exp_b_metric = "/data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter_metric/20260311_240601Y0012_Run0003.newmodel.stage1.norna-nofwdrev.noleftalign.smc_all_reads.aligned.metric.csv"


df_a = pl.read_csv(exp_a_metric, separator="\t")
df_b = pl.read_csv(exp_b_metric, separator="\t")

print(df_a.head(10))
print(df_b.head(10))
