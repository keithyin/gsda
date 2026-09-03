# gsetl

`gsetl`（genomic sequencing ETL）是一个用 **Rust** 编写的测序数据 **ETL（Extract–Transform–Load）** 工具集。它的核心任务是读取 **BAM** 文件（含参考序列的比对结果，以及未比对的原始 basecall），并生成一批结构化的 **事实表（fact table，CSV）**，用于下游的测序质量 / 比对质量分析。

它分为两大类输入：

- **已比对 BAM（aligned）**：read 已经比对到参考序列上，工具逐 record 解析 CIGAR，统计匹配 / 错配 / 插入 / 删除、比对一致性（concordance）、等长同源重复（homopolymer）区域行为、碱基质量分布等。
- **未比对 BAM（non-aligned）**：原始 basecall（带 per-base 时序信息，如 dwell time / arrival time / capture rate），工具按 channel（read）聚合出各类基础统计，以及首尾 N 个碱基的稳定性统计。

> 本仓库通过 **git subtree** 镜像在主仓库 `gsda` 的 `gsetl/` 目录下维护（见文末 [与主仓库的关系](#与主仓库的关系)）。

---

## 目录

- [特性](#特性)
- [构建](#构建)
- [CLI 总览](#cli-总览)
- [子命令：`aligned-bam`](#子命令aligned-bam)
  - [输出文件](#aligned-bam-输出文件)
  - [关键概念](#aligned-bam-关键概念)
- [子命令：`non-aligned-bam`](#子命令non-aligned-bam)
- [子命令：`non-aligned-bam-seq-n-stats`](#子命令non-aligned-bam-seq-n-stats)
- [代码结构](#代码结构)
- [测试](#测试)
- [依赖](#依赖)
- [与主仓库的关系](#与主仓库的关系)

---

## 特性

- 单二进制、多子命令，`clap` derive 驱动的命令行。
- 已比对 BAM 的多张事实表**并行生成**（`thread::scope` + 多个进度条）。
- 支持通过 **BED（confidence regions）** 限定只统计可信区域内的位置，并通过 **VCF（known variants）** 排除已知位点，避免把真实变异误判为“错误”。
- 参考序列自动扩展出 `name` / `name___fwd` / `name___rev`（反向互补）三套，以支持**链特异性（fwd/rev）**的比对解析。
- 未比对 BAM 支持按 CPU 核数多线程解码、流式聚合，内存友好。

## 构建

```bash
cargo build --release
# 二进制位于 target/release/gsetl
./target/release/gsetl --help
```

## CLI 总览

全局参数（所有子命令都接受）：

| 参数 | 说明 |
| --- | --- |
| `--outdir <DIR>` | 输出目录。不存在会自动创建；已存在时不加 `-f` 会沿用。 |
| `-f, --force` | 若 `--outdir` 已存在，先整个删除再重建（清空旧输出）。 |

无论运行哪个子命令，输出目录里都会先生成一个 **`meta.txt`**，记录当前二进制版本（`CARGO_PKG_VERSION`）与完整命令行，便于溯源。

```
gsetl --outdir <DIR> [-f] <subcommand> [子命令参数]
```

子命令：

```
aligned-bam                       # 已比对 BAM → 多张事实表（核心功能）
non-aligned-bam                   # 未比对 BAM → 逐 channel 基础统计表
non-aligned-bam-seq-n-stats       # 未比对 BAM → 首尾 N 碱基稳定性统计
```

---

## 子命令：`aligned-bam`

对一个已比对 BAM + 参考序列，**并行**产出多张事实表。

```bash
gsetl --outdir out/ \
  aligned-bam \
  --bam        my.bam \
  --ref-file   my.fasta \
  --hcregion   confident_regions.bed \
  --vcf        known_variants.vcf \
  --rnames     A B C \
  --useSeco --useSupp
```

### 参数

| 参数 | 必填 | 说明 |
| --- | --- | --- |
| `--bam <PATH>` | 是 | 输入 BAM（需已索引，工具按参考区间 `fetch`）。 |
| `--ref-file <PATH>` | 是 | 参考序列，支持 `fa/fasta/fna`、`fq/fastq`、`bam` 三种格式。 |
| `--hcregion <BED>` | 否 | **可信区域** BED。给定后只有落在这些区域内的位置才会被统计。 |
| `--vcf <VCF>` | 否 | **已知变异位点** VCF。命中这些位点的位置会被排除（不计为错配/错误）。 |
| `--rnames <NAME>...` | 否 | 只处理指定参考名（可多个），默认全部。 |
| `--useSeco` | 否 | 是否纳入 secondary alignment（默认排除）。 |
| `--useSupp` | 否 | 是否纳入 supplementary alignment（默认排除）。 |
| `--factRecordStat <0/1>` | 否 | 是否生成 `fact_record_stat`，默认 `1`。 |
| `--factRefLocusInfo <0/1>` | 否 | 是否生成 `fact_ref_locus_info`，默认 `1`。 |
| `--factBamBasic <0/1>` | 否 | 是否生成 `fact_bam_basic`，默认 `1`。 |
| `--factErrorQueryLocusInfo <0/1>` | 否 | 是否生成 `fact_error_query_locus_info`，默认 `1`。 |
| `--factBaseQStat <0/1>` | 否 | 是否生成 `fact_baseq_stat`，默认 `1`。 |
| `--factPolyInfo <0/1>` | 否 | 是否生成 `fact_poly_info`，默认 `1`。 |

各 `fact*` 开关可单独关掉某张表，按需裁剪输出。

### 过滤（audit）

每条 record 在统计前先经过统一的 `audit` 过滤：

- 未比对（unmapped）→ 跳过；
- secondary → 除非加了 `--useSeco`，否则跳过；
- supplementary → 除非加了 `--useSupp`，否则跳过。

### `aligned-bam` 输出文件

| 文件 | 粒度 | 内容 |
| --- | --- | --- |
| `meta.txt` | run | 版本 + 完整命令行。 |
| `fact_aligned_bam_bam_basic.csv` | 每 record | 基础比对信息：`qname refname channel np rq iy ec rstart rend qstart qend qlen fwd ori_start ori_end`（含参考起止、query 起止、query 长度、正/反向、以及扩展字段 `ori_start/ori_end`）。 |
| `fact_aligned_bam_record_stat.csv` | 每 record | 比对质量汇总：`qname queryCoverage concordance concordanceQv matchBp mismatchBp nonHpInsertionBp nonHpDeletionBp hpInsertionBp hpDeletionBp ignoreBp`。 |
| `fact_aligned_bam_ref_locus_info.csv` | 每参考位点 | `refname pos eq diff ins del depth curBase nextBase curIsHomo nextIsHomo aroundBases diffDetail insDetail`，其中 `diffDetail`/`insDetail` 以 `A:3,T:1` 形式记录差异碱基/插入序列的具体计数，`aroundBases` 以 `前文[当前]后文` 给出 10bp 上下文。 |
| `fact_error_query_locus_info.csv` | 每错误片段 | `qname qstart qend rstart rend qseq rseq`，把一段连续的错配/缺口区间抽取出来，并给出 query/ref 两侧带 `-` 缺口的对齐片段（错误段用 `[...]` 高亮）。 |
| `fact_baseq_stat.csv` | 每参考 × baseq | `refname baseq eq diff ins del depth`，按碱基质量值分桶统计匹配/错配/插入/删除。 |
| `fact_poly_info.csv` | 每 record × 每同源重复区 | `rName qName rStart rEnd qStart qEnd rBase rRepeats qRepeats qSeq qClean`，记录该 read 在每个参考同源重复区上读到的序列与重复数（用于分析重复区扩张/收缩）。 |

### `aligned-bam` 关键概念

- **concordance（一致性）** = `match / (match + mismatch + ins + del)`；`concordanceQv` 是把错误率换算成 PHRED 式质量值（`-10·log10(1 - concordance)`，下限截断）。
- **hp / non-hp indel（同源重复型 / 非同源重复型）**：一个插入/删除若涉及同一碱基的连续重复，记为 `hp`（homopolymer），否则 `nonHp`。这对判断重复区特有的 indel 错误很有用。
- **可信区域（`--hcregion`）与已知变异（`--vcf`）**：`hc_regions` 是“白名单区域”（区域外不统计），`hc_variants` 是“黑名单位点”（命中则不计为差异）。两者叠加用于把“真实生物学差异”和“测序/比对错误”区分开。
- **链特异性扩展**：每个参考名会被展开成 `name`、`name___fwd`、`name___rev`（反向互补），以支持按链分别解析比对。

---

## 子命令：`non-aligned-bam`

对未比对 BAM 逐 channel 解析 per-base 时序/质量信息并输出统计表。

```bash
gsetl --outdir out/ \
  non-aligned-bam \
  --bam raw_called.bam \
  -o out/raw_called.non_aligned_fact.csv \
  -t 16
```

### 参数

| 参数 | 必填 | 说明 |
| --- | --- | --- |
| `--bam <PATH>` | 是 | 输入 BAM。 |
| `-o, --o <PATH>` | 否 | 输出文件路径，默认 `<outdir>/<bam 去扩展名>.non_aligned_fact.csv`。 |
| `-t, --t <N>` | 否 | BAM 解码线程数，默认物理核数。 |

### 输出

`{...}.non_aligned_fact.csv`，每个 channel 输出 A/C/G/T 四行：

```
qname  base  base_cnt  dw_sum  ar_sum  cr_mean  cq  oe
```

- `dw_sum` / `ar_sum`：该 channel 上各碱基 dwell time / arrival time 之和；
- `cr_mean`：capture rate 均值；
- `cq` / `oe`：channel 级指标（缺失时记为 `-1`）。

> 说明：per-base 的 `dw`/`ar`/`cr` 来自 BAM tag；当某 record 缺少相应 tag 时会以 `-1` 占位。

---

## 子命令：`non-aligned-bam-seq-n-stats`

评估未比对 read **首尾 N 个碱基**的时序稳定性（用于判断 read 两端信号质量是否退化）。

```bash
gsetl --outdir out/ \
  non-aligned-bam-seq-n-stats \
  --bam raw_called.bam \
  -n 20 \
  --length-percentile-thr 5
```

### 参数

| 参数 | 必填 | 说明 |
| --- | --- | --- |
| `--bam <PATH>` | 是 | 输入 BAM。 |
| `-n, --n <N>` | 是 | 取首尾各多少个碱基。 |
| `--length-thr <N>` | 二选一 | 长度阈值，仅统计长度 ≥ 该值的 read。 |
| `--length-percentile-thr <P>` | 二选一 | 长度百分位阈值（0–100），据此计算长度阈值。 |

> `--length-thr` 与 `--length-percentile-thr` 至少要给一个，否则报错退出。

### 输出

- `{...}.seq_n_stats.csv`：`qname` + 首/尾 N 的 dwell / arrival 的 `median` 与 `mean` 共 8 列；
- `{...}.first-n.fasta` / `{...}.last-n.fasta`：每个 read 的首 N / 尾 N 碱基序列。

实现上用 `crossbeam` 有界 channel 做“读 BAM → 多线程统计 → 汇聚写出”的流水线，读端线程数与统计线程数都按核数自适应。

---

## 代码结构

```
gsetl/
├── Cargo.toml
└── src/
    ├── main.rs                       # 入口：解析 CLI、建输出目录、写 meta.txt、分发子命令
    ├── cli.rs                        # clap 定义：Cli 与各子命令参数
    ├── utils.rs                      # 通用工具：快速选择算法求中位数
    ├── poly_n.rs                     # 同源重复区（homopolymer）识别与位置关系
    ├── aligned_bam_etl/              # 已比对 BAM 的各张事实表
    │   ├── mod.rs                    #   FastaData（含 fwd/rev 扩展）、audit、BED/VCF 加载
    │   ├── fact_bam_basic.rs         #   每 record 基础比对信息
    │   ├── fact_record_stat.rs       #   每 record 比对质量汇总（多线程 CIGAR 解析）
    │   ├── fact_ref_locus_info.rs    #   每参考位点 eq/diff/ins/del/depth + 上下文
    │   ├── fact_error_query_locus_info.rs  # 错误片段抽取 + 双侧对齐展示
    │   ├── fact_baseq_stat.rs        #   按 baseq 分桶的 eq/diff/ins/del
    │   └── fact_poly_info.rs         #   同源重复区读入序列统计
    └── non_aligned_bam_etl/          # 未比对 BAM
        ├── mod.rs                    #   模块声明
        ├── fact_bam_basic.rs         #   逐 channel base 统计（dw/ar/cr/cq/oe）
        └── seq_n_stats.rs           #   首尾 N 碱基稳定性统计（含聚合实验代码）
```

## 测试

```bash
cargo test
```

部分单测（如 `fact_poly_info` 的 CIGAR 解析用例）是纯逻辑测试；`non_aligned_bam_etl` 下个别测试依赖本机绝对路径的真实数据文件，缺失时以 `eprintln` + 提前返回处理，不会阻断其余测试。

## 依赖

- `clap`（derive）：命令行解析
- `rust-htslib`：BAM 读写 / 索引 `fetch`
- `gskits`：内部工具库（fastx 读取、BED/VCF 解析、BAM record 扩展、DNA 工具、进度条封装）
- `crossbeam`：线程间有界 channel（生产者–消费者流水线）
- `polars`（lazy）：数据聚合（用于聚合实验代码路径）
- `ndarray`、`num_cpus`、`rand`：数值 / 并发 / 随机（含确定种子采样）
- `indicatif`：进度条

## 与主仓库的关系

本仓库（`keithyin/gsetl`）通过 **git subtree** 挂载在主仓库 `gsda` 的 `gsetl/` 目录下。日常可在 `gsda` 工作区内直接编辑 `gsetl/` 下的代码，再通过 subtree 与本仓库同步：

```bash
# 在 gsda 中把 gsetl/ 的提交推送到本仓库
git subtree push --prefix=gsetl gsetl main

# 从本仓库拉回最新 gsetl 内容
git subtree pull --prefix=gsetl git@github.com:keithyin/gsetl.git main --squash
```

---

# change log

## 0.4.0
* remove useless column
* bug fix
