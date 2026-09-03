# ref_range_query_base_q_dist.py

参考区域 query baseQ 四类分布分析工具。

在 `align_mismatch_dist_and_q_dist.py`（只统计 mismatch 位点 baseQ）的基础上，圈中 reference 上的一块区域 `[start, end)`，统计该区域内四类 CIGAR 操作对应的 query baseQ 分布。

## 处理流程

1. **读取参考与 query**
   - `--ref`：单序列 FASTA，读取后校验只有 1 条序列；
   - `--query`：按 `.bam` 后缀自动选择读取方式
     - FASTQ → `read_fastq`，rq 为 `None`，后续由 `readsq_from_baseq` 计算；
     - BAM → `read_bam`（pysam），rq 从 BAM tag `"rq"`（错误概率）换算：`-10 * log(1 - rq)`（自然对数，与 FASTQ 路径的 `log10` 口径不同）。

2. **区域校验**
   - 必须满足 `0 <= start < end <= len(ref)`，否则报错退出。

3. **rq-thr 过滤**
   - `--rq-thr`（默认 `20.0`）：低于该值的 read 直接跳过，计入 `filtered_by_rq`。
   - FASTQ 路径用 `readsq_from_baseq` 从 baseq 算 readsq（对每个碱基的错误概率 `10^(-q/10)` 求均值再转回 PHRED，**不能**直接对 baseq 求算术均值）。

4. **mappy 比对**
   - `mappy.Aligner(seq=ref_seq, extra_flags=67108864, k=11, w=1, best_n=10, n_threads=1)`。
   - 负链 hit 由 `revcomp(query)` / `qual[::-1]` 统一为正链坐标，`q_st = len(query) - q_en`。

5. **三重阈值过滤**（命中 `is_primary` 的 hit）
   - `qcov = (q_en - q_st) / len(qseq) >= --min-qcov`（默认 `0.5`）
   - `rcov = (r_en - r_st) / len(ref)  >= --min-rcov`（默认 `0.5`）
   - `identity = CIGAR 中 `=` 列数 / 总列数 > --min-identity`（默认 `0.95`）
   - 未通过 → 计入 `filtered_by_cov_id`。

6. **baseQ 提取**（`extract_baseq_by_op`）

   走 CIGAR，按操作类型分别收集 baseq：

   | CIGAR op | 含义 | 取值规则 |
   |---|---|---|
   | `=`  | eq | 该列 query base 的 baseq（逐列推进 `q_pos`、`r_pos`） |
   | `X`  | mismatch | 该列 query base 的 baseq（逐列推进） |
   | `I`  | insertion | 插入的**每个** query base 各取一个 baseq；ref 不动，用 `r_pos` 判窗口 |
   | `D`  | deletion | 取 del 两端 query base 的最小值 `min(qual[q_pos-1], qual[q_pos])`，重复**裁剪后**的 del 长度 `k` 次；del 位于 read 边界（`q_pos==0` 或 `q_pos==len(qual)`）时整条排除 |

   **边界收缩**（对所有 op 统一适用，不限 mismatch）：
   - 用 `_homo_run_from_start(ref, r_st)` 和 `_homo_run_from_end(ref, r_en)` 量出两端同聚物（poly）run 长度；
   - 内侧不足 2 个但外侧紧邻碱基同类（边界从中间切开一个 poly）时按 2 处理（宁可多剔边界碱基）；
   - 防止两端 poly run 收缩重叠（read 极短或整段 poly），超过 `aln_span` 的部分先减 `n_start` 再减 `n_end`；
   - 有效分析窗口 = 用户区域 `[start, end)` ∩ core 区 `[core_st, core_en)`，只有落在窗口内的 baseq 才收集。

7. **绘图**（`plot_baseq_distributions`）
   - grouped bar chart：x 轴 = baseQ bin（宽 5），每个 bin 4 根柱子（eq / mismatch / insertion / deletion）；
   - 每根柱高度 = **该 op 类型内部**该 baseQ bin 的占比（各类各自归一，比例之和 = 1），保证分布形状可比；
   - 空 bin / 无数据时输出 `No baseQ values to plot.`；
   - 输出文件：`{output-prefix}_baseq_distributions.png`。

8. **stdout 概览**
   - `Aligned: N / M (qcov>=…, rcov>=…, ident>…)`
   - `filtered_by_cov_id: …`
   - `filtered_by_rq: …`
   - 每类 op：`count`、`mean`、`range`（min–max）。

## 用法

```bash
python ref_range_query_base_q_dist.py \
    --ref ref.fa --query reads.fq \
    --start 0 --end 1000 \
    --output-prefix out \
    --min-identity 0.99
```

## 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `--ref` | str | 必填 | 参考 FASTA（单序列） |
| `--query` | str | 必填 | query FASTQ 或 BAM（按 `.bam` 后缀自动选择） |
| `--start` | int | 必填 | 参考区域起点（0-based, inclusive） |
| `--end` | int | 必填 | 参考区域终点（0-based, exclusive） |
| `--output-prefix` | str | `refrange_baseq` | 输出文件前缀 |
| `--min-identity` | float | `0.95` | 最低 identity（CIGAR `=` 占比） |
| `--min-qcov` | float | `0.5` | 最低 query 覆盖度 |
| `--min-rcov` | float | `0.5` | 最低参考覆盖度 |
| `--rq-thr` | float | `20.0` | 最低 read quality（PHRED），低于该值的 read 过滤 |

## 依赖

- `mappy`、`matplotlib`、`tqdm`
- `pysam`（仅 BAM 路径）
- 同仓库 `gseda.align_ana.mappy_ext`（提供 `revcomp`、`parse_cigar`）

## 与 `align_mismatch_dist_and_q_dist.py` 的差异

| 项 | 旧脚本 | 本脚本 |
|---|---|---|
| 分析范围 | 全参考 | 指定区域 `[start, end)` |
| 收集 op | 仅 mismatch | `=` / `X` / `I` / `D` 四类 |
| 统计口径 | 每个错配位点 1 个 baseq | insertion 每 base 各取；deletion 取两端 min 重复 del 长度 |
| 绘图 | 三联 subplot（位置分布 / 分位点 / 标准差） | grouped bar chart（每类内部归一） |
