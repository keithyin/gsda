# CCS 三类错误概率预估 — 方案设计（讨论稿）

> 状态：**方案合理性讨论已收敛，未进入实施**。本文是讨论结果的汇总，
> 供后续验证 / 迭代 / 实施作为基线。关联：[`demand.md`](./demand.md)、[`research.md`](./research.md)、[`ont_approach.md`](./ont_approach.md)。

---

## 0. 一句话方案

把 CCS（SMC consensus）相对 reference 的每个位点错误，按**列类型**拆成三个独立二分类模型：
SNP（M 列）、deletion（subreads 有支撑但 SMC 未 call 碱基的位点）、insertion（M+I 列），
用 **medaka 风格的 MSA 窗口特征 + GBDT** 各自预测校准概率，输出 `Q_snp / Q_ins / Q_del` 三个通道。

---

## 1. 三个模型的定义

统一前提：**参考 = 该分子的真实序列（Sanger 同分子参考）**，因此"consensus 相对 ref 的差异"即错误。

| 模型 | 位置集合 | label（正/负） | 输出 |
|---|---|---|---|
| **A (SNP)** | M 列（SMC 有碱基 且 ref 有碱基） | SMC 碱基 ≠ ref 碱基 → 1；一致 → 0 | `P(snp)` → `Q_snp` |
| **B (deletion)** | subreads 有支撑、但 SMC 未 call 碱基的位点（SMC 列为 `-` 且有 subread 非 gap 覆盖） | 该位点能与 ref 构成 eq（ref 有碱基、SMC 漏 call）→ 1；否则 0 | `P(del)` → `Q_del` |
| **C (insertion)** | M 列 + I 列（SMC 有碱基，ref 该位无碱基/有碱基） | 是 I 列（ref 无碱基）→ 1；是 M 列 → 0 | `P(ins)` → `Q_ins` |

说明：
- 三类错误在比对空间住在**不同列类型**，互斥，三个二分类是良定义的；
- 三个模型**独立训练**，输出**三个独立 Q 通道**，不做归一化（`P(correct)+P(snp)+P(ins)+P(del)≠1` 是可接受的权衡）；
- 不输出修正后的序列（只出概率），**不做序列重写**。

---

## 2. 特征设计（MSA 窗口 + GBDT）

对目标位点 t，取窗口 [t−k, t+k]，特征来自另一程序产出的 MSA 行（SMC、reference、subreads…）：

- **每列 subread 碱基计数**：A/C/G/T/- 各计几票 → (2k+1)×5 维
- **SMC 碱基** one-hot（目标位点，或窗口内）
- **上下文标量**：均聚长度（SMC 侧）、到最近 gap 距离、覆盖深度、pass 数 N
- **Model B 额外**：该位点 subreads 支撑度（非 gap 的 subread 数 / depth）

模型形态：**GBDT**（XGBoost/LightGBM），输出经 **Platt/isotonic 校准** 后转 Phred：
`Q = -10·log10(p)`（注意 p 要 clip，避免 log(0)）。

---

## 3. 已知软肋与权衡（讨论中明确）

1. **Model C 信息量有限**（原理层面）：
   - Sanger=truth 下，I 列恒为插入错误、M 列恒非插入 → 正负例即列类型；
   - 模型学到的是"像不像 I 列"，对"真实插入被比对成 mismatch"的错位情况**没有训练监督**，抓不到；
   - 它的价值 = 校准的连续 P(ins)（而非硬 0/1）。要抓错位插入需另设 label（如比对不确定样本）。
2. **三个独立概率不归一**：接受。每个通道含义是"该位点是该类错误的可信度"，下游勿当概率分布用。
3. **列类型依赖比对**：模型按列类型分诊，但 SMC↔ref 比对在均聚区不可靠 → 错位会分诊错模型。窗口特征可部分缓解，不根治；结构差距 vs medaka 式 seq2seq。
4. **ref=truth 的适用范围**：真实 SNP / 真实 STR 长度多态会被误标为错误，训练数据须排除。
5. **工程注意**：类别不均衡（M 海量、I/D 稀少）→ class weight / 欠采样；评估用 reliability curve + per-channel P/R，勿只看 AUC。

---

## 4. 数据管线待定点（尚未决策）

- **Model B 候选来源**（已定）：subreads 有支撑、但 SMC 未 call 碱基的位点
  （SMC 列为 `-` 且有 subread 非 gap 覆盖），由 MSA 直接提取，不依赖 SMC bam 的 `di` tag。
- 注意：候选的「subreads 有支撑」只看 subreads 侧；ref 该位是否有碱基只用于定 label，两者分开。

---

## 5. 验证路径（后续）

用现有数据做**诊断**，先验证"三类错误的特征可分离性"：
1. 取若干记录，在候选位点上按 label 分组统计 subread 证据分布（真 deletion（ref 有碱基）vs 假阳性（ref 无碱基）的 subread 支撑度、均聚上下文差异）；
2. 用少量手工特征（majority ratio、均聚长度）画分布，看 A/B/C 是否可分；
3. 若可分，再决定是否进入 GBDT 实施。

---

## 6. 结论

方案**成立、可落地**。三个软肋中：
- Model C 信息量有限、独立概率不归一 → 接受的设计权衡；
- 列类型依赖比对 → 结构性限制，窗口特征缓解但不根治。

作为"每位置三类错误的可信度"第一版，站得住；更进一步再考虑 joint 模型或 seq2seq。
