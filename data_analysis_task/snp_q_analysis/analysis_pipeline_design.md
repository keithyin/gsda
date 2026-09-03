本项目目标用来分析 SNP 位点的 Q 值分布是怎么样的

1. @/data1/ccs_data/20260805-xunming/run1,  @/data1/ccs_data/20260805-xunming/run2 中是两次运行得到的 fastq 文件，相同的 barcode 数据进行合并
2. @/data1/ccs_data/20260805-xunming 中的 fasta 是 barcode样本 的 参考序列。 将 合并的 barcode 样本，进行 Q30 过滤后 比对到 参考序列上，然后输出 SNP 位点。
3. 统计SNP位点的 Q 值分布。

上述分析的每一步 都需要一个分析脚本承接，将分析脚本写到 @/root/projects/gsda/data_analysis_task/snp_q_analysis 中，实现之后对当前数据进行分析