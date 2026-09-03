# ONT (Nanopore) 方案详解：如何区分 SNP / 插入 / 缺失

> 本文是 [`research.md`](./research.md) §3.3 的展开。聚焦 **Oxford Nanopore (ONT)** 生态里"把
> base 质量与 indel 分开建模"的做法——业界对"per-base Q 装不下 indel"这个问题最系统的一组解法。
> 核心对象：**Nanopolish HMM**（理论骨架）、**Clair**（后验→质量桥接）、**DeepConsensus / medaka**
> （learned indel 校正）、**Dorado**（生产 basecaller）、**Nanocompore**（评估）。
>
> **可信度标注约定（重要）**：
> - 标 **[通用]** = 标准 HMM / ML 数学或业界公认做法，不依赖具体论文，可直接信。
> - 标 **[结构]** = 该方法的公开架构/思路（encoder-decoder、头结构、监督目标），方向可靠。
> - 标 **[实参]** = 具体超参 / 数值 / 版本号——**离线无法核实，仅作示意**，落地前须对照原文复验。
>
> 本文所有**数值示例都显式标注为"示意"**，用来演示机制，**不代表任何论文/工具的真实取值**。

---

## 0. 为什么 ONT 最相关

- 短读 (Illumina) 的 baseQ 天生是"替换置信度"，对 indel 结构盲（见 `research.md` §3.1）。
- ONT 读长长、**均聚 overcall/undercall 突出**，迫使社区把"长度对不对"当一等公民。
- 于是 ONT 生态形成一条清晰分层，每一层都对应"区分 SNP/ins/del"的一个环节：
  1. **basecaller**：raw current → 碱基 + (单标量) 质量；
  2. **indels-aware posterior 模型 (HMM)**：把 match/mismatch/ins/del 作为**显式状态**，出后验 **[通用]**；
  3. **后验 → 质量 的桥接**（Clair）[结构]；
  4. **learned consensus/indel 校正**（DeepConsensus / medaka）：直接预测"该 base 该不该留" [结构]；
  5. **评估**：把 accuracy / indels / homopolymer 分成不同指标（Nanocompore）[结构]。
- 我们要的"Q 值区分 SNP/ins/del"，在这条线里**每一层都已有对应物**——可以直接对着抄结构。

---

## 1. Nanopolish：HMM 的完整数学（理论骨架，最值得学）

**出处** [实参]：Hickey et al., *Nanopolish: fast secondary analysis of raw nanopore data* (2017, Bioinformatics)。

### 1.1 建模对象
Nanopolish 直接建模**原始电流信号 (raw current)** 而非碱基文本 [通用]。
每个 read 是一串电流值 `x_1, ..., x_T`。它把"碱基错"翻译成"电流与期望值偏离"，
从而得到比 basecaller 更接近物理真值的后验 [结构]。

> 移植要点 [通用]：你们 CCS 的"raw"是 **pass 的碱基+Q**（离散），不是连续电流；
> 但**状态空间与解码数学完全通用**——把"发射观测"从"电流高斯"换成"pass 一致性"即可。
> 下面先按 ONT 原版（连续电流）讲清数学，再给"离散版"。

### 1.2 状态空间 [通用]

把比对过程建模成隐马尔可夫链。定义参考序列 `R = r_1...r_L`（长度 L），
query/consensus `Q = q_1...q_M`。相对参考，每个参考位点落在互斥状态：

| 状态 | 语义 | 对应需求 | 记号 |
|---|---|---|---|
| Match | 该位点碱基正确 | `P_c` | `M` |
| Mismatch | 该位点存在但判成另一碱基 (3 种) | `P_s (SNP)` | `X` |
| Insertion | query 在参考此位点附近多出碱基 | `P_i (插入)` | `I` |
| Deletion | 参考有位点、query 没有 | `P_d (缺失)` | `D` |

Nanopolish 常用**精简 4 状态**或**完整 (N+1)×4** 状态；
教学上下面用**逐参考位点的 {M, X, I, D}** 展开，够用 [通用]。

### 1.3 转移概率 (transition) [通用]

设状态序列 `z_1...z_L`（每个参考位点一个状态）。相邻状态转移约束来自序列对齐的合理性：

```
允许转移(示意):  M→{M,X,I,D},  X→{M,X,I,D},  I→{M,I},  D→{M,I}
```
- 从 indel 态回到参考消费位点 (M/X)；I 允许连续 I（连插）；D 允许连续 D（连缺）也可接 I。
- 转移概率体现 **gap 成本**：常用 **affine gap penalty**
  `γ_open + γ_ext·(len-1)`，转成概率 `P(转移)`：
  `P(z_t = M|z_{t-1}=M)` 最高，进出 indel 态概率低。**这些概率就是"indel 代价参数化"的一种形式** [通用]。

> [实参] 具体 gap 数值随版本变（Nanopolish 用 EM 从数据估），**不要抄某个具体数字**。

### 1.4 发射概率 (emission) [通用 + ONT 特化]

**A. 连续电流版（Nanopolish 原版）** [结构]：
- `M` / `X` 状态：该电流 `x_t` ~ `N(μ_{state,context}, σ²)`。
  `μ` 是"该状态+碱基+上下文(k-mer)"的期望电流；`σ` 反映位置不确定性。
  ```
  P(x_t | z_t=M, b)    = 𝒩(x_t; μ_{M,b,ctx}, σ_M²)
  P(x_t | z_t=X, b)    = 𝒩(x_t; μ_{X,b,ctx}, σ_X²)      # 被判错成 b
  ```
- `I` / `D` 状态：不直接绑定参考碱基，其观测似然由对齐转移承载（或建模为"多余/缺失的电流段"）[结构]。

**B. 离散版（可直接移植到 CCS）** [通用]：
把"发射"换成"pass 证据"。设某位点有 N 条 pass，第 n 条在该参考位点附近给出碱基 `b_n`（含 gap）：
- `P(pass 证据 | z_t=M)`：多数 pass 与参考一致 → 高；
- `P(z_t=X → 判为 b)`：pass 支持 b 的比例；
- `P(z_t=I)`：pass 在此处系统性多一个碱基的比例；
- `P(z_t=D)`：pass 在此处缺一个碱基的比例。

> 这正是把 `eq/diff/ins/del` 四个计数直接喂进四个状态后验的**自然实现**。[通用]

### 1.5 解码与前向后向

**Viterbi（最可能路径）** [通用]：
```
δ_t(z) = max_{z_{1:t-1}} P(z_{1:t-1}, x_{1:t}) ,  backpointers 复原
```
→ 输出一条含 ins/del 的最可能状态路径，用于**校正序列 / 比对**。

**前向-后向 (forward-backward)** [通用]：
- 前向 `α_t(z) = P(z_t=z, x_{1:t})`
- 后向 `β_t(z) = P(x_{t+1:T} | z_t=z)`
- **单点/配对后验**：
  ```
  P(z_t = z | x_{1:T})  =  α_t(z)·β_t(z) / Z        # Z = Σ α_Lβ_L
  P(z_{t-1}=z, z_t=z') =  α_{t-1}(z)·P(z→z')·emission·β_t(z') / Z
  ```

**这就是"把标量 Q 展开成状态后验"的具体机制**：
- `P(SNP) = P(z_t=X | x) = P(z_t=X, 具体碱基 b) 求和`
- `P(插入) = P(z_t=I | x)`
- `P(缺失) = P(z_t=D | x)`
- `P(正确) = P(z_t=M | x)`

### 1.6 后验 → 质量（Nanopolish `refine-read` 与 `call-var` 的差别）[通用 + 结构]

- `refine-read`：把**单点后验折成一个 Phred**（对"该位点正确"）：
  ```
  Q_t = -10·log10 ( 1 - P(match正确碱基 | x) )
  ```
  **⚠ 这就是"单通道 Q"的来源——它只编码了 match/SNP 侧后验**；
  indel 信息在状态后验里有，但**没被折进这个 Q**。[通用]
- `call-var`：**直接读各状态后验**做变异：
  - SNP ← `P(z=X, b) / P(z=M)` 似然比；
  - indel ← `P(z=I)` / `P(z=D)` 后验（**独立于 SNP 判据**）。
  **这就是"三类错误分开对待"的直接落地** [结构]。

### 1.7 关键洞察：Nanopolish 的 gap 在"质量输出"这一层
它把 ins/del 做成**模型/变异检测**层的一等公民，但**默认质量串仍是 SNP 偏置的单通道**。
→ 你们做 `Q_snp / Q_len` 时，**indel 通道必须像 `call-var` 那样显式取 `P(I)+P(D)` 后验**，
而不是一路从单 Q 里找。[通用]

---

## 2. Clair：后验 → base 质量的桥接 [结构]

**出处** [实参]：Iqbal et al., *Clair*（为 Nanopore read 生成 base quality，基于 Nanopolish 后验）。

- 作用：Nanopolish HMM 给的是**状态后验**，Clair 把"该位点为匹配态"折成 Phred：
  ```
  Q_k = -10·log10 ( 1 - Σ_b P(该位点正确为 b | 电流, 参考, 模型) )
  ```
  [通用] 即 **Q = 对"替换正确性"的后验**。
- 意义 [结构]：
  - 固化工程链路 **HMM 状态后验 → Phred**；
  - 明确区分两个概念——**"碱基正确性"(SNP)** vs **"存在性/长度正确性"(indel)**，
    二者应是两条独立后验 → 两个独立 Q。这正是你们 `Q_snp` / `Q_len` 的概念基础。[通用]

---

## 3. DeepConsensus（PacBio/Google）：learned 校正，用"gap"显式建模 indel [结构]

**出处** [实参]：Baid et al., *DeepConsensus: Gap-Aware Sequence Transformers for Sequence Correction*
（2021 预印本 [biorxiv 10.1101/2021.08.31.458403](https://www.biorxiv.org/content/10.1101/2021.08.31.458403v1.full)；
2023 Nat. Biotechnol.，[PubMed 36050551](https://pubmed.ncbi.nlm.nih.gov/36050551/)，[google/deepconsensus](https://github.com/google/deepconsensus)）。
> 注：这是 **PacBio CCS** 的官方 learned consensus 方法（Google/PacBio 合作）——比 ONT 的 medaka
> 更贴我们的场景。它在 CCS 基础上再纠错，且**把 indel 显式建模**。

### 3.1 输入/输出
- 输入：**subreads + 原始 pulse/电流信号**（做 consensus 纠错）。
- 输出：校正后的 consensus 序列；**显式建模插入/缺失** [结构]。

### 3.2 模型结构：gap-aware 序列 Transformer [结构]
- 对一组对齐的 subreads 逐列编码，用 **gap 标记（gap token）** 表示"该位置在别的 subread 里是缺失/插入"；
- 输出层对每个位置预测**碱基（A/C/G/T）或 gap**——**indel 是一个显式的输出类别**，不是折进碱基概率；
- 相比经典 `ccs`（逐位点 + 误差模型压成一个 QV），DeepConsensus 的 gap 建模直接降低 consensus 的错误率。

### 3.3 对均聚/indel 的意义 [结构]
- **把"该不该有这个 base"当独立预测目标（gap 类）**，而非折进同一置信度——
  这正是"把 SNP 与 indel 分开对待"的 neural 实现，且来自 CCS 官方生态。

### 3.4 对本需求的启示 [通用]
- **learned 模型可显式输出 gap（indel）类别**：`P(base正确)`(SNP) + `P(该位置为 gap / 长度错)`(indel)；
- 两头可各自折成一个 Phred（→ `Q_snp`、`Q_len`），也可都保留概率供下游；
- 若我们走 learned 路线，**DeepConsensus 的 gap-aware 结构是比"双头 CNN"更贴 CCS 的模板**。

---

## 4. medaka：ONT 官方 learned consensus（生产可用）[结构]

**出处** [实参]：medaka (ONT, 2020/2021；Transformer 架构，早期 CNN)。

- 从 raw signal 或一组成簇 pass 直接产出**校正序列 + 每 base 质量**；
- **显式处理 indel**：对均聚 overcall/undercall 做 learned 校正；
- **supervised**（有参考真值训练/微调）与 **unsupervised**（无参考）两模式；
- 与 DeepConsensus 思路同源（learned base + indel 校正 + 质量），但**更工程化、可直接产可用结果**。
- **生产级证明**：indel 校正 + 质量可一体化，不必拼 HMM。
- 实用建议 [结构]：即便你们走非-learned 路线（`research.md` §5.2 两通道经验模型），
  也可把 **medaka 当基准/教师**：用它输出的"indel 校正置信度"来标定我方 `Q_len` 通道。

---

## 5. Dorado / 生产 basecaller：质量仍是单通道的现状 [结构]

- **Dorado** (ONT, ~2022/2023) [实参]：raw signal → 序列 + 每 base Q 的生产级 basecaller
  （standard / super-accuracy 两档）。
  - 现状：输出**每 base 一个 Q 标量**（替换侧为主）——与 Illumina 一样，
    indel 未显式区分。印证"默认 baseQ 都是 SNP 偏置的"。
- **Nanocompore** (ONT) [结构]：读长 basecalling **指标分离**基准
  （accuracy / indels / homopolymer 分开报）——**连评估层面都把 indel 单列**，
  进一步佐证"indel 需独立"这一共识。
- 可作对照的旁证 [实参]：Shasta/Rayon (longread 校正，含 indel)、Racon (参考校正)、PonG (GAN + 质量)。

---

## 6. 数值示例（**示意，仅演示机制；非任何工具/论文实参**）

> 全部数字是我构造的小例子，用来把"三通道 → 聚合 Q"和"后验 → 质量"跑通，
> **不代表 Nanopolish/DeepConsensus/medaka 的实际取值**。

### 6.1 从状态后验算三通道 Q（一个均聚位点）
设某均聚区参考位点，前向后向给出状态后验（示意）：
```
P(M)=0.982,  P(X)=0.002,  P(I)=0.008,  P(D)=0.008        sum=1
```
则：
```
Q_snp = -10log10(0.002) = 27.0
Q_ins = -10log10(0.008) = 21.0
Q_del = -10log10(0.008) = 21.0
Q_len = -10log10(0.008+0.008) = -10log10(0.016) = 18.0
Q_agg = -10log10(0.002+0.008+0.008) = -10log10(0.018) ≈ 17.4
```
解读：**S/I/D 都还行，但 ins+del (0.016) 是 SNP (0.002) 的 8 倍**，
`Q_agg≈17.4` 主要由 indel 拖低（`Q_agg` 落在最弱通道 ~Q21 附近、再被 SNP 略拉低）。
- 若只看 `Q_agg`，你只知道"这个位点 ~Q17"；
- 拆开才看到**"低是因为 indel 而不是错的碱基"**——在均聚区这正是你要的（碱基本就全对）。

### 6.2 同一个 Q_agg，两种不同病因（为何必须拆）
```
情形 A (SNP 主导):  P_s=0.015, P_i=0.0005, P_d=0.0005  → Q_agg=18.0, Q_snp=18.2, Q_len=30.0
情形 B (indel 主导): P_s=0.0005, P_i=0.0075, P_d=0.0075 → Q_agg=18.0, Q_snp=33.0, Q_len=18.2
```
**两个位点 `Q_agg` 完全相同(=18.0)，但一个"碱基错得多"、一个"长度错得多"**。
- 对下游基因分型/变异过滤，二者该走不同处理（SNP 换碱基 / indel 改长度），
  单 Q 无从区分。**这就是"区分对待"的全部价值所在。**

### 6.3 均聚 copy 数分布 → 质量（把 §2.4 跑通）
均聚 `(T)ⁿ`，参考 `n_true=9`，N 条 pass 给出的共识长度分布（示意）：
```
n=7:0.04  n=8:0.11  n=9:0.61  n=10:0.18  n=11:0.06
```
```
P(长度对)      = P(n=9)      = 0.61
Q_rep          = -10log10(1-0.61) = -10log10(0.39) ≈ 4.1    # 很差！
P(n+1)=P(n=10) = 0.18 ;  P(n-1)=P(n=8)=0.11   # ±1 歧义度合计 29%
```
解读：**即便参考是 9，consensus 有 39% 概率长度错**（过/欠 call 各向扩散）；
`Q_rep≈4` 远低于"碱基都对"的情形——**这正把"均聚区质量差"量化出来了**，
而不是被 per-base Q 摊平成一堆"半对的 T"。[通用机制]

### 6.4 用 N 条 pass 的多数表决估 SNP 通道 [通用公式]
单 pass 错误率 `e`（示意 e=0.05）。N 条多数表决下"多数也投错"的概率：
```
N=1 : P_snp = e = 0.05                 → Q_snp = 13.0
N=5 : P_snp = P(Bin(5,0.05)≥3) ≈ 1.16e-3         → Q_snp ≈ 29.4
N=10: P_snp = P(Bin(10,0.05)≥6) ≈ 2.8e-6  → Q_snp ≈ 55.5
```
要点 [通用]：**SNP 通道随 N 快速下降**（多数表决，每加 pass 近似平方级改善）；
而 §6.3 的 indel/copy 偏置**不随 N 下降**（系统性过判）——
两通道对 N 的敏感度不同，**这正是它们必须分开建模的原因**。[通用结论]

---

## 7. ONT 方案总览：一张表

| 方案 | 类型 | 模型态持有 ins/del？ | 质量输出形态 | 对本需求价值 |
|---|---|---|---|---|
| **Nanopolish** | HMM (电流级) | ✅ match/mis/ins/del 状态 [通用] | `refine` 默认单 Q(SNP)；`call-var` 用 ins/del 后验 [结构] | **状态后验→三通道** 范式 |
| **Clair** | HMM 后验→Q | 间接（SNP 侧） | 单 Q = 替换正确性后验 [通用] | 后验→Phred 桥；"Q=正确性后验" |
| **DeepConsensus** | learned (CCS) | ✅ gap 标记显式建模 indel [结构] | consensus 校正（gap 类别） | **gap-aware：indel 是显式输出类别** |
| **medaka** | learned (生产) | ✅ indel 校正 [结构] | 每 base Q（含 indel 校正） | 生产证明"indel 一等公民" |
| **Dorado** | basecaller | ⚠️ 默认单 Q [结构] | 每 base 单 Q（替换偏置） | 反例：默认 Q 不区分 indel |
| **Nanocompore** | 基准 | ✅ 指标分离 [结构] | — | 评估都分 accuracy/indels/homopolymer |

**主线**：ONT 把"区分对待"落到**模型态 → 质量输出 → 评估指标**三层，
逐层对应你们 `Q_snp / Q_len / Q_agg` + 分通道校准——ONT 是现成模板。

---

## 8. 逐项移植映射：ONT 概念 → CCS Q 通道设计

> 这是把上一节"ON T 怎么做的"翻译到"你们 CCS Q 值"的对照。[通用]

| ONT 里的做法 | 机制 | 你们 CCS 的对应 |
|---|---|---|
| HMM 状态 {M,X,I,D} | 把错误分成显式状态 | `Q_snp / Q_ins / Q_del`（或 `Q_snp / Q_len`）三通道 |
| 发射观测 = 电流后验 | P(碱基正确\|证据) | "发射"换成 **pass 一致性**（`eq/diff/ins/del` 计数→四状态后验） |
| 前向后向 → 状态后验 | 上下文相关解码 | 位点级 `P(M/X/I/D)`（简单版可退化为逐位点计数，不强行上下文） |
| `refine` 单 Q 只编码 SNP | 质量 = 正确性后验 | 默认 `Q_agg` 保留；**额外显式产出 indel 通道**（补 Nanopolish 的 gap） |
| `call-var` 用 P(I)/P(D) | indel 独立判据 | `Q_len = -10log10(P(I)+P(D))` 独立于 `Q_snp` |
| DeepConsensus/medaka gap 建模 | learned 显式 gap 类 | 可选：learned 模型显式预测 `gap/碱基` 类别；非-learned 则经验拟合两通道 |
| 均聚长度为一等公民 | copy 数分布 | `P(n\|n_true,N,hp)` + `Q_rep`（§6.3） |
| Nanocompore 指标分离 | 分指标评估 | `Q_snp`/`Q_len` 各自校准 + P/R；同 `Q_agg` 验证 `Q_len` 分离能力 |
| BQSR 式上下文校准 (GATK) | 按 ctx 分 bin | `P_s(N,ctx)` / `P_len(N,hp,unit)` 两张表 |

**与 ONT 的关键差异（移植时必须留意）** [通用]：
- ONT 是**单条 raw read + 连续电流**；CCS 是 **N 条离散 pass → 共识**。
  → 你们"发射"用离散 pass 一致性，`N` 是显式输入（ONT 单 read 没有这个 N 维度）——
  **CCS 反而更适合"随 N 的误差模型"**（ONT 单 read 依赖电流，CCS 依赖 pass 数目+质量）。
- ONT 的 indel 偏置主要来自**实时过 call**；CCS 主要来自**均聚长度判定 + 参考 copy 数差异**——
  `P_len` 模型的偏置项要针对这两种来源分别建模。

---

## 9. 对本需求的可操作结论（4 条）

1. **质量 = 对"正确性"的后验**，不是对"某碱基"的置信度（Clair/Nanopolish）；
   → 据此把质量拆成"碱基正确性"(SNP) 与"存在性/长度正确性"(indel) 两条。[通用]
2. **模型侧持有 ins/del 为显式状态或 head**（Nanopolish HMM / DeepConsensus-medaka），
   → 预估模型应有 indel 分量，而非事后从 Q 切。[通用]
3. **默认 baseQ 是替换偏置的**（Dorado/Illumina 反例），
   → indel 通道必须显式补，别指望单 Q。[结构]
4. **评估把 indel/均聚单列**（Nanocompore）→ 你们的 `Q_snp`/`Q_len` 各自校准 + P/R，
   并在同 `Q_agg` 位点验证 `Q_len` 能分离"低质源于 indel"。[结构]

---

## 10. 参考文献（链接可核实）

- **DeepConsensus**（PacBio/Google，gap-aware transformer，显式 gap 建模）：
  [biorxiv 10.1101/2021.08.31.458403](https://www.biorxiv.org/content/10.1101/2021.08.31.458403v1.full) ·
  [Nat. Biotechnol. 2023 (PubMed 36050551)](https://pubmed.ncbi.nlm.nih.gov/36050551/) ·
  [github.com/google/deepconsensus](https://github.com/google/deepconsensus)
- Hickey et al. (2017) — Nanopolish: fast secondary analysis of raw nanopore data (Bioinformatics)。
  [github.com/nanopolish/nanopolish](https://github.com/nanopolish/nanopolish) · [发射模型 nanopolish_emissions.h](https://github.com/sjackman/nanopolish/blob/master/nanopolish_emissions.h)
- Iqbal et al. (c.2020–2021) — Clair（Nanopore read 的 base quality，基于 HMM 后验）。
- medaka (Oxford Nanopore, 2020/2021) — learned consensus + indel/均聚长度校正。
  [github.com/nanoporetech/medaka](https://github.com/nanoporetech/medaka)
- Dorado (Oxford Nanopore, ~2022/2023) — 生产级 neural basecaller（standard / super-accuracy）。
- Shafin et al. (2020) — Shasta；(2021) — Rayon（longread 校正，含 indel，旁证）。
- Vaser et al. (2017) — Racon（参考校正，含 indel，对照）。
- Nanocompore (ONT) — 读长 basecalling 指标分离基准（accuracy / indels / homopolymer）。
- Keller et al. (2021) — PonG（GAN basecalling + 质量估计，旁证）。
