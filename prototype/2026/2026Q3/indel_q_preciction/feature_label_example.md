# CCS 三类错误：特征 & Label 可视化示例

> 配套 [design.md](./design.md) 的图示：用一个连续的 MSA 例子说明三类模型（SNP / deletion / insertion）的**特征抽取**与 **label 定义**。

## 0. 前提（来自 design.md）

- **ref = truth**（Sanger 同分子参考）；**SMC** = consensus 碱基；**subreads** = 原始读段。
- 三模型按**列类型**分诊，互斥、独立二分类：

| | 候选列（该模型要判的位点） | label = 1（错误）的含义 |
|---|---|---|
| **A (SNP)** | M 列（SMC、ref 都有碱基） | SMC 碱基 ≠ ref 碱基 |
| **B (del)** | SMC 为 `-`，且 subreads 有非 gap 覆盖 | ref 有碱基（SMC 漏 call 了真碱基） |
| **C (ins)** | M+I 列（SMC 有碱基） | 是 I 列（ref 无碱基，SMC 多 call） |

## 1. 一个完整的 MSA 窗口

```
         col:     1   2   3   4   5   6   7   8   9
   ref(truth):   A   C   C   G   T   -   A   G   -
   SMC:          A   C   A   -   T   T   A   G   -
   subread1:     A   C   A   G   T   T   A   G   G
   subread2:     A   C   C   G   T   T   A   G   G
   subread3:     A   C   C   -   T   -   A   G   -
   subread4:     A   -   A   G   T   T   A   G   G
   subread5:     A   C   A   G   T   -   A   G   -
   ─────────────────────────────────────────────────────
   列类型         M   M   M   D   M   I   M   M   D      M=SMC/ref均有碱基  I=SMC有/ref无  D=SMC为-
   Model A 候选  ✓   ✓   ✓   ·   ✓   ·   ✓   ✓   ·      候选 = M 列
   Model B 候选  ·   ·   ·   ✓   ·   ·   ·   ·   ✓      候选 = SMC为- 且 subread有非gap覆盖
   Model C 候选  ✓   ✓   ✓   ·   ✓   ✓   ✓   ✓   ·      候选 = SMC 有碱基 (M+I)
```

这 9 列里，每类错误的正/负例各至少一个：**col3** 是 SNP 错误（正例），**col4** 是 deletion 错误（正例）、**col9** 是其假阳性对照（负例），**col6** 是 insertion 错误（正例）。

## 2. 每个模型的 label 怎么定

**Model A（SNP）** — 只看 M 列，`SMC碱基 ≠ ref碱基` → 1：

| col | SMC | ref | 不一致? | label |
|---|---|---|---|---|
| 1 | A | A | ✗ | 0 |
| 2 | C | C | ✗ | 0 |
| **3** | **A** | **C** | **✓** | **1 ← 唯一的 SNP 错误** |
| 5 | T | T | ✗ | 0 |
| 7 | A | A | ✗ | 0 |
| 8 | G | G | ✗ | 0 |

**Model B（deletion）** — 候选 = SMC 为 `-` 且有 subread 覆盖；`ref 有碱基`（真漏 call）→ 1：

| col | SMC | subread 非gap覆盖 | ref 有碱基? | label | 含义 |
|---|---|---|---|---|---|
| **4** | - | 4/5 (G×4) | G → 有 | **1** | SMC 漏 call 了真碱基 G |
| 9 | - | 3/5 (G×3) | - → 无 | 0 | subreads 的假插入，SMC 正确跳过 |

**Model C（insertion）** — 候选 = SMC 有碱基；`是 I 列`（ref 无碱基）→ 1：

| col | 列类型 | SMC | ref | label |
|---|---|---|---|---|
| 1 | M | A | A | 0 |
| 2 | M | C | C | 0 |
| 3 | M | A | C | 0 |
| 5 | M | T | T | 0 |
| **6** | **I** | **T** | **-** | **1 ← SMC 多 call 了一个真序列没有的碱基** |
| 7 | M | A | A | 0 |
| 8 | M | G | G | 0 |

> 注意 col3：它是 Model A 的正例（SNP 错误），但在 Model C 里因为列类型是 M → label 0。**正/负例 = 列类型 + 一致性，先分诊再定标，这正是三个二分类"良定义"的原因。**

## 3. 特征怎么从窗口抽出来（以 col3 为 Model A 目标，窗口 k=2）

对目标位点 t=col3，取窗口 `[t-2, t+2]` = col1–5：

**① 每列 subread 碱基计数** → (2k+1)×5 = **25 维**

| col | A | C | G | T | - |
|---|---|---|---|---|---|
| 1 | 5 | 0 | 0 | 0 | 0 |
| 2 | 0 | 4 | 0 | 0 | 1 |
| **3 (t)** | **3** | **2** | 0 | 0 | 0 |
| 4 | 0 | 0 | 4 | 0 | 1 |
| 5 | 0 | 0 | 0 | 5 | 0 |

**② SMC 碱基 one-hot（t 处）**：A=1, C=0, G=0, T=0

**③ 上下文标量（假设值，示意）**：

| 特征 | 值 | 说明 |
|---|---|---|
| 均聚长度（SMC 侧） | 1 | t 处为孤立 A，两侧是 C / `-` |
| 到最近 gap 距离 | 1 | t+1（col4）就是 gap |
| 覆盖深度 | 5 | 5 条 subreads |
| pass 数 N | 12 | 示例假设 |

**④ Model B 独有特征**（以 col4 为目标）：`支撑度 = 非gap subread数 / depth = 4/5 = 0.8`（顺带：col9 的支撑度是 3/5 = 0.6 —— 正负例在支撑度上可分离，正是验证路径 §5.1 想看的分布）。

Model C 的窗口计数（目标 col6，窗口 col4–8）同理：

| col | A | C | G | T | - |
|---|---|---|---|---|---|
| 4 | 0 | 0 | 4 | 0 | 1 |
| 5 | 0 | 0 | 0 | 5 | 0 |
| **6 (t)** | 0 | 0 | 0 | **3** | **2** |
| 7 | 5 | 0 | 0 | 0 | 0 |
| 8 | 0 | 0 | 5 | 0 | 0 |

（t 处 SMC=T，one-hot；T 在 subreads 占 3/5 多数 → 模型会把它判成"像真的碱基"，而 truth 说没有 → label=1。）

## 4. 小结

| 模型 | 候选列 | 特征 | label=1 的定义 | 本示例 |
|---|---|---|---|---|
| A | M 列 | 窗口 25 维计数 + SMC one-hot + 3 个标量 | SMC ≠ ref | col3 |
| B | SMC 为-且有覆盖 | 同上 **+ 支撑度** | ref 有碱基（漏 call） | col4（正）/ col9（负） |
| C | M+I 列 | 同 A | 是 I 列（ref 无碱基） | col6 |

一个完整训练样本 = **「窗口特征向量」+「该模型在该列的 0/1 label」**；同一列在不同模型里可以是不同样本（如 col3 是 A 的正例、C 的负例）。
