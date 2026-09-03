关于 16S consensus 后疑似 SNP、INDEL 位点的分析

1. 先从 扩增子的整个流程中将存在疑似提取出来。然后使用以下脚本看一下情况

```
python /root/projects/gsda/gseda/src/gseda/16s/ref_range_query_base_q_dist.py \
    --ref Group_0/barcodes_reads_cons_gen_amplicon/Consensus/Sequences/Group_0_Adaptor-barcode277-0.consensus.fasta \
    --query 20260831_240601Y0014_Run0005_called_barcode_v4/demuxed/Adaptor-barcode277-0.fastq \
    --start 1091 --end 1092
```

```
会有一个这样的输出，就能看到比较关注的区域的 不同情况的 Q值分布。还会有一个图
Per-op baseQ summary:
  eq         count=     958  mean=30.7  range=4-45
  mismatch   count=       2  mean=31.0  range=26-36
  insertion  count=      11  mean=21.1  range=6-45
  deletion   count=      40  mean=27.9  range=14-45


还会输出这个信息，如果想看更精细的 subreads 粒度的信息，可以看一下 asrtc 的结果
  deletion (40):
    20260831_240601Y0014_Run0005/483326/32-3542
    20260831_240601Y0014_Run0005/262877/32-3586
    20260831_240601Y0014_Run0005/109577/32-3546
    483326,262877,109577,206527,99365,379533,361048,334925,399390,456325,352224
```


2. 执行以下命令来查看 asrtc 的结果

```bash
asrtc \
    --ref-fa Group_0/barcodes_reads_cons_gen_amplicon/Consensus/Sequences/Group_0_Adaptor-barcode277-0.consensus.fasta  \
    -q 20260831_240601Y0014_Run0005_called_demuxed_v4_new_model.bam \
    -t 20260831_240601Y0014_Run0005_called_barcode_v4/demuxed/Adaptor-barcode277-0.fastq  \
    -p barcode277-newmodel \
    --channel-whitelist 483326,262877,109577,206527,99365,379533,361048,334925361589 \
    --ref-range 1080:1099 \
    --tn-delim "/" \
    --ch-idx 1
```

3. 开始人眼看一下感兴趣的区域，是什么个情况

## 两个 basecalling 模型 A/B 对比（对 consensus 的 identity）

当拿到两套 barcode 拆分后的 CCS reads（例如 baseline 模型 vs 优化后的模型），想客观比较哪个模型
basecalling 更准时，用 `ab_consensus_identity.py`。它把每个 barcode 的 reads 比对到该 barcode 对应的
consensus 序列，复用 `sequencing_report_v2` 的口径（gsmm2-aligned-metric），输出池化指标 + 逐 barcode 比较。

```bash
conda run -n py38 python /root/projects/gsda/gseda/src/gseda/16s/ab_consensus_identity.py \
    --run-a  <baseline 的 demuxed fastq 目录> --label-a baseline \
    --run-b  <新模型 的 demuxed fastq 目录>   --label-b new_model \
    --consensus <Consensus/Sequences 目录> \
    --prefix Group_0_ --suffix .consensus.fasta \
    --rq-range 0.99:1.1 \
    --outdir /tmp/ab_eval
```

- 配对规则：`Adaptor-barcodeXXX-Y.fastq` ↔ `{prefix}Adaptor-barcodeXXX-Y{suffix}`，
  只保留 A、B 都有且 consensus 非空的 barcode（严格逐 barcode 一一对应）。
- `--rq-range 0.99:1.1`：只统计自评高置信（rq≥0.99）的 reads，与各模型"自报高置信"集合对比。
- 输出：`outdir/per_barcode_identity.csv`（逐 barcode A/B/diff）、`outdir/summary.txt`（池化指标 + A/B 胜出数）。

## SUSPICIOUS_SITES 新旧模型对比（降噪 vs 敏感性）

当拿到同一源数据的两组 consensus suspicious_sites 输出（旧模型 vs 新模型），量化新模型相对旧模型
在**降噪**与**敏感性**上的差异，用 `suspicious_sites_compare.py`。方法定义见同目录
`suspicious_sites_model_comparison.md`（D1 位点总量 / D2 命中样本数 / D3 逐位点分类 / D4 归因画像）。

```bash
conda run -n py38 python /root/projects/gsda/gseda/src/gseda/16s/suspicious_sites_compare.py \
    --new <新模型 Suspicious_sites 目录> \
    --old <旧模型 Suspicious_sites 目录> \
    --format txt \
    --outdir /tmp/susp_cmp
```

- 两个目录各含一批 `<sample>.consensus.suspicious_sites.{txt,csv}`（`.txt` tab 分隔、`.csv` 逗号
  分隔，同一数据的两种格式）；`--format` 指定用哪一种，**只取其一不叠加**。字段按表头名解析
  （真实文件 20 列，含 `Name`/`Mut_context`）。
- 以 `(样本, Position)` 为键分三类：
  - **CLEANED**：OLD 有、NEW 无 → 降噪
  - **NEW_SUS**：NEW 有、OLD 无 → 新爆
  - **SHIFTED**：同一样本内 indel 坐标整体重定位（邻位成簇、类型一致、深度量级一致），
    单独标注、**不计入净增减**。
- 每个差异位点给归因标签：低比例 SNP / 超低深度 indel / 重复同聚区 indel / 其它。
- 输出：`outdir/suspicious_sites_comparison.txt`（D1–D4 汇总 + 逐样本逐位点明细），同时打印到 stdout。
