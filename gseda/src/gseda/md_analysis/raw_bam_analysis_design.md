
目标：对于不同 run 的 molecular dynamic. 将其 一些统计的柱状图画到一个 histgram 上。用来观察不同 run 之间的 molecular dynamic 的差异

原始数据产出：
    调用 `md_analysis -i /data1/ccs_data/20260717-rna-report/20260529_250802Y0002_Run0001_called.bam  raw-bam` 这个会生成两个 csv 文件，具体什么样的csv文件，以及其命名规则是怎么的，请阅读代码 @/root/projects/md_analysis

不同 bam 会生成不同的 csv 文件，将这些csv文件，相同类别的 画在同一张图上 。用来进行比对 run 间差异。
