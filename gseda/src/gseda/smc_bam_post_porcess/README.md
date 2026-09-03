# smc_bam_post_porcess

对 smc 共识（consensus）BAM 中 `di` 字段做后处理与质量分析的脚本集合。

## 背景：`di` 字段

smc 输出的 BAM 包含一个 `di` 字段。它是 string，但实际内容是一个 list，用 `;` 分隔。
每个单元由 `pos,base,frac,depth,phreq` 构成：

```
pos,base,frac,depth,phreq;pos,base,frac,depth,phreq;...
```

- `pos`    — int，插入位置（在 smc seq 中的下标）
- `base`   — ACGT，要插入的 base
- `frac`   — float，支撑比例（support fraction）
- `depth`  — int，支撑深度
- `phreq`  — float，`= -10*log(1-p)`，原始 phred 错误分（无其它缩放）

**语义**：`di` 的每个元素表示 subreads 在该位点有支撑，但 smc 输出的 seq 上并没有该 base
（典型场景是 SMN 低估了均聚物长度）。`pos` 可能重复，表示多个受支撑 base 落在同一位置。

> **重要**：本目录所有脚本读写的 BAM **都不是比对结果 BAM**，而是未映射（unmapped）的
> smc 共识序列。因此一律用 `pysam.AlignmentFile(..., check_sq=False)` 打开，并用
> `fetch(until_eof=True)` 遍历（没有 CIGAR / 参考坐标）。

### `di` 质量分析的结论（见 `di_hi_phred_eda.png`）

>Q20 的高分 `di` base 其实分**两类**：

| 组 | n | 特征 |
|---|---|---|
| 自然尾部 Q20–39 | ~87,751 | lq-med 32、depth-med 8、frac-med 0.27，属正常的长尾 |
| **Q45 饱和尖峰** | ~21,448 | lq-med 28（反而更低）、depth 5–9（94%）、frac-med 0.28 |

所有 >Q20（n≈109,199，占全部 di 的 2.19%）的共同点：
- 均聚物上下文占比高（88.6% vs 低分组的 83.6%）
- 中位共识 run 长度 = 2 拷贝（run≥2 占 62.6%）→ 正是"漏掉一个拷贝"
- support count（中位 2 subreads）、frac（中位 0.28，最低 ~0.20，从未 =1.0）
- local qual（中位 31，反而低于 Q≤20 组的 37）

即：高分 `di` base 绝大多数发生在**均聚物**上，subreads 稳定地支持比共识多一个拷贝，
而 Q45 是一个饱和上限（hard cap，在 Q40 附近断崖）。

## 脚本清单

| 脚本 | 作用 | 输入 → 输出 |
|---|---|---|
| `insert_di.py` | 将 `di` 中的 base 插回 smc seq（qual 用 `phreq`），重算整读质量并更新 `rq` | `*.bam` → `*.with_di.q<Q>.bam` |
| `di_vs_seq_qual.py` | 画 `di` 的 `phreq` 与 smc 共识 seq-qual 的质量分布对比直方图（带 legend，可 rq/np 过滤） | `*.bam` → `*.di_vs_seq.png` |
| `di_q45_dist.py` | 统计命中目标 phreq（默认 45）的 `di` 单元距离读两端（start/tail）的距离分布 | `*.bam` → `*.di_q45_dist.png` |
| `dump_di_phreq.py` | 打印 `di` 中存在 phreq 落在指定区间（默认 30:45）单元的读名（channel）及其 `di` 信息 | `*.bam` → stdout / TSV |
| `extract_di_examples.py` | 抽取并可视化"均聚物上下文里高分 `di` base"的 6 个具体读示例 | `*.bam` → stdout + `/tmp/di_examples.tsv` |
| `eda_di_qual.py` | 直接基于 BAM（或 `di_units.npz`）打印 >Q20 与 ≤Q20 两组指标对比，并画 2×3 图 | `*.bam`（或 `*.npz`）→ stdout + PNG |
| `di_homopolymer_check.py` | 校验 `di` 是否只落在 homopolymer 上下文：`base == seq[pos]` 或 `base == seq[pos-1]` 判为正确，否则错误 | `*.bam` → stdout（可选 TSV） |

- `di_hi_phred_eda.png` — 上面结论的成品图（`di_hi_phred_eda.png`），保留在目录中供参考。
- `insert_di.md` / `di_vs_seq_qual.md` — 对应脚本的原始设计说明。

## 用法

各脚本均为独立可执行脚本（`python <script>.py`）。除 `extract_di_examples.py`
（输入 BAM 路径硬编码在文件顶部 `IN`，运行前按需修改）外，输入 BAM 都作为首个
位置参数传入。绘图/产物默认写在输入文件同目录或 `/tmp/`，均可用 `-o` 覆盖。

### insert_di.py

把 `di` 的 base 插回共识 seq，qual 取 `phreq`，重算整读 `rq` 后写出新 BAM。

```bash
python insert_di.py <input.bam> [--output out.bam] [--threads 40] [--no-progress]
                    [--poly-only] [--q-cut 20]
```

- 输出默认 `<input>.with_di.q<Q>.bam`（`<Q>` 即 `--q-cut` 取值，不加 `--q-cut` 时为 `q0`），
  这样不同阈值的多次运行不会互相覆盖，产物文件名自带其过滤条件。
- 多个相同 `pos` 的单元按 `di` 出现顺序依次插入；`pos == len(seq)` 的单元落在末尾。
- `--poly-only`：只插入与插入位点共识 base 相同（`base == seq[pos]`，即延长均聚物 run）的
  `di` 单元，其余单元丢弃；`pos >= len(seq)` 的末尾单元没有可比对的 base，也一并丢弃。
  实测（ecoli 样本前 20 万读）保留 83.5% 的 `di` 单元。
- `--q-cut Q`：`Q` 为整数质量值，只插入自身 `phreq >= Q` 的 `di` 单元（`phreq` 已是原始 phred
  分值，直接与未取整的原值比较、含等号；默认 `0` 即不过滤）。与 `--poly-only` 可叠加，两个条件都要满足。
  注意 `di` 的 `phreq` 分布严重偏低（ecoli 样本前 5 万读、45.9 万个单元：中位数 5、最大 45，
  `>=20` 仅 2.39%、`>=30` 0.67%、`>=40` 与 `>=45` 均为 0.39%），阈值取高会把插入几乎全删掉。
- 无 `di` 或 `di` 解析为空的读原样写回。
- `rq` 重算公式 `rq = mean(1 - 10**(-q/10))`，与 `gseda/file_format_cvt/fastq2bam.py:calc_rq` 一致。

### di_vs_seq_qual.py

对比"缺失 base（`di`）质量"与"共识 seq 质量"的分布。

```bash
python di_vs_seq_qual.py <input.bam> [-o out.png] [--min-rq 0.9] [--min-np 3]
                          [--threads 4] [--linear] [--no-progress]
```

- 输出默认 `<input>.di_vs_seq.png`。
- y 轴默认 log scale（两组样本量差异大）；`--linear` 切换线性。
- `--min-rq` / `--min-np` 按读级 `rq` / `np` tag 过滤（缺 tag 视为不满足过滤条件）。

### di_q45_dist.py

命中目标 phreq 的 `di` 单元，计算其插入位置到最近读端的距离 `min(pos, len(seq)-pos)` 并画分布。

```bash
python di_q45_dist.py <input.bam> [-o out.png] [--phreq 45]
                      [--threads 4] [--linear] [--no-progress]
```

- 输出默认 `<input>.di_q45_dist.png`；统计 min / median / mean / max 及"距端 ≤50nt 占比"。

### dump_di_phreq.py

打印 `di` 单元 phreq 落在目标区间内的读名与信息。

```bash
python dump_di_phreq.py <input.bam> [--phreq 30:45] [--threads 4] [--no-progress]
```

- `--phreq LO:HI` 为**闭区间**（如 `30:45` 表示 phreq ∈ [30, 45]），默认 `30:45`；也可只写
  单个值（`--phreq 45` 等价于 `45:45`）。比较前 phreq 先四舍五入到整数。
  `LO > HI` 或格式非法（如 `1:2:3`）会直接报错退出。
- 仅当单元 `pos` 落在 `(10, 20)` 区间时才 print（`pos > 10 and pos < 20`）——这是调试期
  留下的过滤，且**只作用于输出、不作用于计数**，因此 stderr 汇总里的单元数会大于实际打印行数。
  汇总计数走 stderr。

### extract_di_examples.py

抽取 6 个"高分 `di` base 位于均聚物"的具体读，打印可视示例并导出 TSV。

```bash
python extract_di_examples.py
```

- 过滤：`np >= 3`、`rq >= 0.9`、`di` 单元 `phred > 20`、且该位点存在均聚物上下文。
- 输出：stdout 可视化 + `/tmp/di_examples.tsv`。输入 BAM 路径硬编码在文件顶部 `IN`。

### eda_di_qual.py

对 `di` 单元做 >Q20 与 ≤Q20 两组对比（打印指标 + 2×3 图）。**直接基于 BAM 构建
数据**，也可从预先生成的 `di_units.npz` 载入。

```bash
# 直接从 BAM 跑（推荐）
python eda_di_qual.py <input.bam> [--min-np 3] [--min-rq 0.9]
                      [-o out.png] [--threads 4] [--dpi 140] [--no-progress]

# 或从已有的 .npz 载入
python eda_di_qual.py --npz /tmp/di_units.npz
```

- 提供 `input_bam` 或 `--npz` 二选一；两者都给时以 `--npz` 优先。
- 数据 = 每个 `di` 单元一行：单元 5 字段（`phreq,pos,base,frac,depth`）+
  读/上下文特征（`ndi`、`rq`、`lq`、`homo`）。`lq` 取插入点两侧 4 个现有共识
  base（`seq[pos-2..pos+1]`，端点截断）的 qual 均值。
- `--min-np`（默认 3）/ `--min-rq` 按读级 tag 过滤（缺 tag 视为不满足）。
- 输出 PNG 默认为输入文件同目录下的 `<输入名去扩展名>.di_q_dist.png`
  （如 `run.with_di.bam` → `run.with_di.di_q_dist.png`；用 `--npz` 时前缀取 npz 名；
  用 `-o` 可完全覆盖该路径）。
- 控制台额外打印 `di` 的 Q 值（phred）分布直方图，按 `floor(Q)` 分箱，
  条形取对数刻度（详见 `print_phred_distribution`）。
- `.npz` 为可选中间缓存，不再是必需输入。

### di_homopolymer_check.py

校验 `di` 是否只落在 homopolymer 上下文。**判定规则**：单元 `(pos, base)` 满足

```
base == seq[pos]   或   base == seq[pos-1]      → 正确；否则 → 错误
```

```bash
python di_homopolymer_check.py <input.bam> [--min-np 3] [--min-rq 0.9] [--q-cut 20]
                               [--examples 5] [--dump-bad bad.tsv] [--threads 4] [--no-progress]
```

- 两种匹配的含义相同：插入的 base 接上了一段已存在的同 base run
  （`seq[pos]` 向右延长，`seq[pos-1]` 向左延长）。
- **该规则比"位于均聚物"更宽松**：`run == 1`（只有一侧一个拷贝）也算正确，但它不是通常意义的
  homopolymer。所以脚本同时统计"与插入点连续的同 base 拷贝数 `run`"（`classify()` 从插入点向两侧
  走到 base 不同为止），`run >= 2` 才是严格口径；报告里两个口径都给。
- `pos < 0` 或 `pos > len(seq)` 无邻位可比 → 判错；`pos == len(seq)`（末尾单元）只有 `seq[pos-1]` 可比。
- `--min-np`（默认 3，传 `0` 关闭）/ `--min-rq` 为读级过滤，缺 tag 视为不满足；
  `--dump-bad` 才落盘全部错误单元（默认只打印前 `--examples` 条）。
- 输出仅 stdout 报告；无图。

**实测**（`...adapter.withdel.3.smc_all_reads.bam`，np≥3，541,552 读 / 6,910,205 个单元）：
正确 86.78%、错误 13.22%、严格口径（正确且 `run ≥ 2`）48.72%、含错误单元的读占 46.29%。
`withdel.2` 给出几乎相同的全局数字（86.78% / 13.22%）。

即 **`di` 并非只出现在 homopolymer 区域**：约 1/8 的单元两侧都没有同 base 拷贝，
另有 38% 只是"单拷贝邻位匹配"。但**按 phred 分组的结论不稳定**，不要据单次运行下判断：

| BAM | >Q20 单元数 | >Q20 正确率 | ≤Q20 正确率 |
|---|---|---|---|
| `withdel.2` | 263,035 | 95.36% | 86.44% |
| `withdel.3` | 42,331 | 70.74% | 86.88% |

高分单元的数量与构成在两版之间变化很大（.3 里 >Q20 少了 84%），方向因此相反。
把 `--q-cut` 提到 40（`withdel.3`）更能看出问题：这一档共 6,540 个单元，只有 **47.71%** 正确。
`--q-cut` 取 39 / 41 / 43 / 44 / 44.5 得到的单元数完全相同（6,540），说明这批就是 **Q45 饱和尖峰本身**
——即最高分的那批 `di` 里有一半以上，插入点两侧都没有同 base 拷贝，既不符合规则也谈不上均聚物。
（注：不带 `withdel` 的 `...adapter.smc_all_reads.bam` 根本没有 `di` tag——`di` 由
deletion-aware 的 smc 流程产生。）

## 依赖

`pysam`、`numpy`、`matplotlib`、`tqdm`。绘图脚本统一用项目调色板（`#2a78d6` 蓝 /
`#eb6834` 橙，双色色盲安全）。

## 设计说明

各脚本刻意保持独立、不过度抽象；`di_vs_seq_qual.py` / `di_q45_dist.py` / `dump_di_phreq.py`
通过 `from .insert_di import parse_di` 复用同一套容错的 5 字段 `di` 解析器（字段数不对、base
非法、非数字的单元直接跳过），sibling 缺失时有等价 fallback。
`di_homopolymer_check.py` 同样复用该解析器（相对导入失败时退回同名模块导入），
但把"是否延长已有 run"的判断做在自己的 `classify()` 里，与 `insert_di.filter_poly_only()`
（只看 `seq[pos]`）和 `eda_di_qual.py` 的粗粒度 `homo` 标记互不依赖。
