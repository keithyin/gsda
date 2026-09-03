# CCS baseQ 预估：分离 SNP / 插入 / 缺失 —— 详细调研

> 独立调研。核心问题：**CCS 输出的 Q 值应当能对 SNP (substitution)、插入 (insertion)、缺失 (deletion) 三类错误区分对待**，
> 而当前 CCS 的 Q 值把三者混在一个标量里。
>
> 调研方法：平台与工具链调研（Illumina / PacBio CCS / ONT / GATK / 变异检测 / STR 基因分型）+ 文献核实。
> 参考文献均在文末带可核实的 URL。WebFetch 在此环境被禁，正文描述基于公开论文/仓库 README 的标题与摘要级信息 + 领域通行知识；
> 具体内部超参以文末链接为准。

---

## 0. 结论先行

1. **结构问题（不是调参问题）**：per-base Q 是"这个碱基是不是对的"的置信度，是**替换 (SNP) 的原生量**。
   对 indel：
   - **插入** = "这里本不该有这个碱基" —— 正确答案是*没有碱基*，不是一个"错的碱基"；
   - **缺失** = "这里本该有个碱基但没有" —— 输出里根本没有这个位置，无处挂 Q。
   所以**一个 per-base 标量对 indel 是降维**：它只能把 indel 折算成"这个 base 不太可信"，丢掉了"错在哪一维"的信息。

2. **拆不出（信息论）**：一个标量（1 自由度）**反解不出**三类各自概率（3 自由度）。
   `Q = -10·log10(P_snp+P_ins+P_del)` 是 many-to-one，不可逆。
   → 拆分必须发生在 **CCS 的误差模型侧**（预估源头），不能事后从 Q 里切。

3. **更好的表达 = 把 Q 从"标量"升为"对局部编辑操作的分布"**：
   每个 consensus 位点输出 `{正确, SNP, 插入, 缺失}` 的概率（或等价的三通道 `Q_snp / Q_ins / Q_del`）。
   这是序列比对的状态空间本身，不是新发明——只是把"标量"还原成它本来的"分布"。

4. **业界已有先例（可核实，见 §7）**：
   - **Nanopolish**（HMM 显式 match/mismatch/ins/del 状态）；
   - **DeepConsensus**（PacBio/Google，gap-aware transformer——**直接以"gap 标记"显式建模 indels 的 CCS 共识器**，与本需求同源）；
   - **medaka / DeepCon**（learned indel 校正 + 质量）；
   - **GATK BQSR / BAQ**（按上下文、**按"附近有 indel"重新校准 base 质量**——BAQ 就是 indel-aware 的 baseQ 校正先例）；
   - **samtools consensus** 为高 indel 平台改进（PR #1733）；
   - **STR 基因分型器 GangSTR/HipSTR/Stragglr**（copy 数分布 + 读长误差模型）；
   - **VCF `PL`**（quality-as-vector）。

---

## 1. 问题定义

### 1.1 baseQ 是什么（Phred 语义）

- 标准定义：`Q = -10·log10 P(error)`，即 `P(error) = 10^(-Q/10)`，`P(correct) = 1 - 10^(-Q/10)`。
- 出自 Sanger/Phred 体系（Ewing & Green），由测序仪 basecaller 对每个碱基给出；
  形式上它只是"这个碱基判得对不对"的置信度。参见 [drive5 的 Phred 说明](https://www.drive5.com/usearch/manual/quality_score.html)。

### 1.2 CCS 流程与 Q 的来源

CCS (Circular Consensus Sequencing / HiFi) 的核心：同一个环状分子被读多遍（多条 pass/subread），
每条 pass 有自己的 basecaller 质量；把 pass 对齐后做**共识 (consensus)**，每个 consensus 位点得到一个**共识质量 QV**。

```
consensus 位点的 QV =  -10·log10 P(该 consensus base 是错的 | N 条 pass, 局部序列, 误差模型)
```

- 经典实现：PacBio `ccs` 工具（[PacificBiosciences/ccs](https://github.com/PacificBiosciences/ccs)，官方描述"Generate Highly Accurate Single-Molecule Consensus Reads"）。
- 现代实现：**DeepConsensus**（Google/PacBio，见 §3.2.2）——用 gap-aware transformer 在 CCS 基础上再纠错。
- CCS 误差本身的分解有文献定量描述（[Error analysis of the PacBio sequencing CCS reads, Int. J. Biostatistics](https://www.x-mol.com/paper/1655979466304081920)）。

### 1.3 三类错误的严格定义（对齐语境）

| 类型 | 含义 | 在输出序列里的形态 | per-base Q 能否表达 |
|---|---|---|---|
| **SNP / 错配** | 该位点碱基存在但判错 | 一个"错的碱基"+ Q | ✅ 原生（base correctness） |
| **插入 (ins)** | 真分子没有的碱基被多出 | 一个"多余的碱基"+ Q | ⚠️ 只能折算成"base 不可信"，丢"不该存在"语义 |
| **缺失 (del)** | 真分子有、输出里没有 | **没有这个位置** | ❌ 无处挂 Q |

> 关键论证：**缺失根本没有对应位置，per-base Q 从定义上就表达不了它**；插入虽然占一个 base，
> 但"对/错"应针对*存在性*而非*碱基身份*。只有 SNP 是 per-base Q 的一等公民。
> 这就是"当前 baseQ 把三者混在一起"的结构性根源。

### 1.4 为什么均聚 / 重复区最严重

在纯均聚 `(X)ⁿ` 里碱基全相同，**唯一自由度是长度 n**——此时 SNP 不存在，错误只能是 indel（±copy）。
per-base Q 在这种区要么把"copy 数错"当成一堆"碱基错"（严重误判），要么给出一个不反映长度不确定度的假性高 Q。
**对 STR/均聚，正确的随机变量是 copy 数，不是碱基。**

---

## 2. 数学框架：把"标量 Q"升成"编辑操作分布"

### 2.1 四类互斥事件 + 三通道

一个比对上的参考位点，consensus 的局部结果落在互斥类别：
```
P_c = P(正确)      P_s = P(SNP)      P_i = P(插入)      P_d = P(缺失)      Σ = 1
```
各取 Phred：
```
Q_snp = -10·log10 P_s        Q_ins = -10·log10 P_i        Q_del = -10·log10 P_d
Q_agg = -10·log10 (P_s + P_i + P_d) = -10·log10 (1 - P_c)     ← 写回现有 Q 字段，保持兼容
```

**性质 1（上界）**：`P_s+P_i+P_d ≥ max(P_s,P_i,P_d)` ⇒ `Q_agg ≤ min(Q_snp, Q_ins, Q_del)`。
**最弱通道决定总质量上限。** 例：`P_s=1e-4`(Q_snp=40)，`P_i=9e-4`(Q_ins≈30.5) ⇒ `Q_agg≈30.46`。
> 一个"替换很好但 indel 很差"的位点，总 Q 会被 indel 拖低、且**低的原因被藏住**——这正是需求要暴露的信息。

**性质 2（不可逆）**：`Q_agg = f(Q_snp, Q_ins, Q_del)` 是 many-to-one。给定 `Q_agg=30`，三通道解无穷多。
→ **从标量永远拆不回三通道**；拆分必须在模型侧（§5）。

### 2.2 "对编辑操作出分布"是自然还原

`{正确, SNP, 插入, 缺失}` 恰好是**局部编辑操作**的集合（= 比对的局部状态空间）。
"拆分后的 Q" 本质是：给每个位点输出一份**编辑操作概率分布**，而不是标量。
SNP 有 3 种具体错碱基 → 也可再细化成 `P(SNP→A/C/G/T)`（信息更全，非必须）。

### 2.3 两分量模型（按机制分：SNP vs indel）

机制上两套源：
- **SNP**：basecalling / 单碱基化学误差，随 N 快速下降（多数表决）；
- **indel**：长度/重复区判定误差，**常带系统性偏置**（如均聚过判），不必然随 N 下降。

所以按机制聚成两个分量最省事：`Q_snp` 与 `Q_len (=-10·log10(P_i+P_d))`。
均聚区 `Q_snp≈0`（无意义），质量几乎全由 `Q_len` 承载。

### 2.4 均聚 / STR：用 copy 数分布表达

对纯均聚 `(X)ⁿ`，正确的质量对象是 **copy 数分布**：
```
P(n | n_true, N, 重复长度, unit)
Q_rep = -10·log10 (1 - P(n = n_true))          # "copy 数对"的位点质量
辅助   = P(n_true+1), P(n_true-1)               # ±1 歧义度，供 STR 基因分型/长度决策
```
这天然把 SNP 与 indel 统一进"长度对不对"这一件事里——对 STR 场景比三通道更贴。

### 2.5 数值示例（示意，演示机制，非实参）

- **同一个 Q_agg，两种病因**：
  - A：`P_s=0.015, P_i=P_d=0.0005` → `Q_agg=18.0, Q_snp=18.2, Q_len=30.0`（SNP 主导）
  - B：`P_s=0.0005, P_i=P_d=0.0075` → `Q_agg=18.0, Q_snp=33.0, Q_len=18.2`（indel 主导）
  → **两个位点 Q_agg 相同，但一个"碱基错得多"、一个"长度错得多"**。单 Q 无从区分。
- **多数表决下 SNP 通道随 N 快速下降**：`P_snp(N) = Σ_{k>⌊N/2⌋} C(N,k) e^k (1-e)^{N-k}`，
  `e=0.05` 时 N=1→Q≈13，N=5→Q≈29，N=10→Q≈55；而 indel 偏置（系统性）**不随 N 下降**。

---

## 3. 调研：业界如何预估 baseQ / 如何处理 indel

### 3.1 短读 (Illumina)：baseQ 天生是"替换置信度"，对 indel 结构盲

- basecaller 每周期输出对 `{A,C,G,T}` 的后验（Phred 模型一脉），baseQ = `-10·log10(1-max_p)`。
- **这是"这个碱基是不是我判的这个碱基"的置信度，根本不建模"这里该有/不该有一个碱基"**。
- 社区共识：baseQ 不能代表 indel 错误；indel 在均聚/重复区是独立于 baseQ 的错误源
  （综述如 [Viral genetic diversity estimation from NGS](https://pmc.ncbi.nlm.nih.gov/articles/PMC3438994/) 对均聚 indel 的讨论）。
- 推论：baseQ 天然是"SNP 通道"，**没有 indel 通道**；要 indel 信息必须另加。

### 3.2 PacBio CCS / HiFi

#### 3.2.1 经典 `ccs` 工具：共识质量是"替换 + indel"误差模型的联合，但输出仍是每 base 一个 QV
- 官方工具 [PacificBiosciences/ccs](https://github.com/PacificBiosciences/ccs)；
  流程与 FAQ 见 [ccs.how](https://ccs.how/faq/mode-all.html)。
- 概念上 QV 由"N 条 pass + 误差模型"给出，误差模型含替换项与 indel 项；但**输出仍压成一个 per-base QV**。
- 局限（与本需求一致）：均聚/重复区的长度不确定度没被单 QV 完整表达。

#### 3.2.2 **DeepConsensus（最重要先例）**：gap-aware transformer，显式建模 indels
- [DeepConsensus: Gap-Aware Sequence Transformers for Sequence Correction](https://www.biorxiv.org/content/10.1101/2021.08.31.458403v1.full)
  （Baid et al., 2023 Nat. Biotechnol.；PubMed: [36050551](https://pubmed.ncbi.nlm.nih.gov/36050551/)，代码 [google/deepconsensus](https://github.com/google/deepconsensus)）。
- 输入：**subreads + 原始 pulse 信号**；输出：比经典 `ccs` 更准的 consensus。
- **关键**：模型用 **gap-aware 序列 Transformer**——把"插入/缺失"作为**显式的 gap 标记**建模，
  而不是把 indel 折进替换概率。**这就是"对 SNP/ins/del 分开对待"在 CCS 上的现代实现，与本需求同源**。
- 启示：PacBio 自己已经走上"用显式 gap 模型区分 indel"的路；我们要的"三通道 Q"是它的质量侧补充。

### 3.3 ONT：把 indel 显式建模的代表（详见 `ont_approach.md`）

- **Nanopolish**（[nanopolish](https://github.com/nanopolish/nanopolish)）：HMM，状态显式区分
  match/mismatch/insertion/deletion（发射模型见 [nanopolish_emissions.h](https://github.com/sjackman/nanopolish/blob/master/nanopolish_emissions.h)）。
  `refine-read` 给精化质量、`call-var` 用 P(I)/P(D) 做 indel 独立判据。
- **medaka**（[nanoporetech/medaka](https://github.com/nanoporetech/medaka)，"Neural network sequence error correction"）：
  learned consensus + 显式 indel/均聚长度校正 + 质量。
- **DeepCon**（Lee et al., 2020, Cell）：CNN 同时校正 SNP + indel，per-base 质量。
- **Dorado**：生产 basecaller，默认仍输出单 Q（替换偏置）——反例。

### 3.4 变异检测侧的 indel-aware 质量（baseQ 之外的成熟做法）

- **GATK BQSR**（Base Quality Score Recalibration）：
  按 (报告的 Q, 化学循环, 二核苷酸上下文) 分箱，用已知真值（dbsnp）重校准 baseQ。
  官方文档：[BQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR)。
  → **"上下文分箱校准单通道 baseQ"的先例**（仍标量，但证明"Q 应按上下文建模"）。
- **BAQ（Base Alignment Quality）—— indel-aware baseQ 校正的直接先例**：
  Li, [Improving SNP discovery by base alignment quality](https://academic.oup.com/bioinformatics/article/27/8/1157/227268)
  (Bioinformatics 2011, doi:10.1093/bioinformatics/btr076)。
  → 在**比对不确定/有 indel** 的位置，用 HMM 计算"该碱基来自参考某位置"的后验，
    **据此把局部 base 质量调低**，以减少 indel 区域附近的假 SNP。
  → 这就是**"质量应感知附近 indel"** 的成熟实现；我们的 `Q_len` 可看作它的"正向版"（不只调低 Q，还显式给出 indel 概率）。
- **GATK HaplotypeCaller pairHMM**：发射概率里**替换项与 indel 项分开建模**——SNP 与 indel 分离是工业级标准。
- **DeepVariant**（[biorxiv 092890](https://www.biorxiv.org/content/10.1101/092890v5.full)）：
  对 pileup 图像做 CNN，SNP 与 indel 分别给候选与质量；后处理里用 `P(NON_REF)/P(ALL)` 表达"该位点非参考的概率"
  （[postprocess_variants.py](https://github.com/google/deepvariant/blob/master/deepvariant/postprocess_variants.py)）——
  **"对多个假说给概率"（quality-as-vector）的又一实例**。

### 3.5 consensus 工具对 indel 的处理（samtools 的改进 = 现成的"痛点证明"）

- jkbonfield, [Improvements to samtools consensus for platforms with higher indel error rates](https://github.com/samtools/samtools/pull/1733)：
  **标准 `samtools consensus` 在高 indel 率平台（如 ONT）上准确率差**——因为逐位点 majority 无法正确处理 indel；
  该 PR 改为更符合对齐的算法，显著提升高 indel 平台的 consensus 质量。
  → **业界共识工具的官方改动就发生在"indel 处理"上**，进一步印证：indel 是 consensus 质量的头号短板。
- 配套：[bcftools 的 indel calling 改进 PR #1679](https://github.com/samtools/bcftools/pull/1679) 与
  [#2045（edlib 替换 BAQ）](https://github.com/samtools/bcftools/pull/2045)。

### 3.6 STR 基因分型器：copy 数分布 + 每-STR 读长误差模型

- **GangSTR**（[GangSTR 四种 reads 建模的评测，F1000](https://f1000research.com/articles/9-200/v1)）、
  **HipSTR**、**Stragglr**。
- 核心：把 STR 质量建模成 **copy 数离散分布** + **每-STR 读长误差模型**，用似然选基因型/长度。
  STR 工具横向对比见 [Genome-wide STR Truth Set 工具对比](https://broadinstitute.github.io/str-truth-set/html/tool_comparison_viewer.html) 与
  [TR 算法评测综述](https://www.nature.com/articles/s41398-023-02689-8)。
- **不是给一个 baseQ，而是给 `P(长度)`**——与 §2.4 完全一致。

### 3.7 表示先例：quality-as-vector

- VCF `PL`（Plo-likelihoods）：对每个候选基因型给 `-10·log10 P(观测|基因型)`——质量是向量而非标量。
- DeepVariant 的 `P(NON_REF)`（见上）：同样是"对非参考假说给一个概率"。
- **结论：多通道/多假说质量在领域里是合法且已在用的表示。**

---

## 4. 候选表示 & 对比

| # | 表示 | 内容 | 优点 | 缺点 | 适用 | 工作量 |
|---|---|---|---|---|---|---|
| **A** | **三通道 Q 向量** | `Q_snp, Q_ins, Q_del`（或 `Q_snp, Q_len`） | 直接满足"区分对待"；可折回 `Q_agg` 兼容下游 | 需在模型侧产出；无现成标准 | 通用（SNP+indel 都要） | 中 |
| **B** | **per-base 发射向量** | `P(A),P(C),P(G),P(T)` + `P(indel)` | 信息最全（知错成哪个碱基） | 对均聚/STR 意义小；维度高 | 非重复区、SNP 为主 | 中-高 |
| **C** | **copy 数分布 `P(n)`** | `P(n\|n_true,N,重复)` | 对均聚/STR **最诚实**；天然统一 SNP/indel | 只对重复区适用 | **STR / 均聚** | 中 |
| **D** | **HMM / edit-graph 后验** | per-position `{match,SNP,ins,del}` 后验 | 上下文感知、最通用（Nanopolish/BAQ 先例） | 建模/计算最重 | 误差强依赖局部上下文时 | 高 |

**推荐组合**：**A 为骨架（`Q_snp` + `Q_len`）＋ 在均聚/STR 位点升级为 C（copy 数分布）**；
D 作为"金标准"升级路径（当发现误差强依赖上下文、两通道不够时再上）。
B 一般不需要（痛点是长度不是碱基身份）。

---

## 5. 如何预估 & 如何落地（设计要点）

### 5.1 目标输出（每个 consensus 位点）

```
普通位点:
  Q_snp = -10log10 P_s(N, ctx)
  Q_len = -10log10 P_len(N, 重复长度, unit)
  Q_agg = -10log10 (P_s + P_len)          # 写回现有 Q 字段（兼容）

均聚 / STR 位点 (长度 >= 阈值):
  P(n | n_true, N, 重复长度, unit)         # copy 数分布
  Q_rep  = -10log10 (1 - P(n=n_true))
  aux    = P(n_true+1), P(n_true-1)        # ±1 歧义度
```

### 5.2 预估：把 Q 建模成 (N, 上下文) 的函数

- **SNP 通道** `P_s(N,ctx)`：
  - 组合模型：单 pass 错误率 `e_snp` → N 条多数表决
    `P_s(N) = Σ_{k>⌊N/2⌋} C(N,k) e_snp^k (1-e_snp)^{N-k}`（随 N 快速下降）；
  - 或按 (N, 化学/序列上下文) 分箱用已知真值校准（BQSR 思路）。
- **长度通道** `P_len(N, 重复, unit)`：
  - 均聚用 copy 数分布 `P(n|n_true)`，`P_len = 1 - P(n=n_true)`；
  - **必须允许一个不随 N 衰减的偏置项**（均聚过判是系统性偏置，不是随机误差）——
    否则会在长均聚上**系统高估**质量。
  - 模型形态可参考 BAQ 的"HMM 后验 → 质量校正"：BAQ 是"附近有 indel → 调低碱基质量"；
    我们做"正向版"——直接给出 `P(indel)`。
- 产出：两张小表 `P_s(N,ctx)` 与 `P_len(N,重复,unit)`（或一个统一小模型）。
  这是"深度细化 Q 预估"的落点——从"一个数"变成"两个随 N/上下文变化的函数"。

### 5.3 存储 / 互操作

- SAM/BAM 的 per-base 质量字符串是**单通道**、且**不支持 per-base 辅助 tag** → 通道信息走旁路：
  - **位点级旁路表**（BED/TSV：`位点 → Q_snp, Q_len, P(n)` 向量）与 BAM 并存；
  - 或 **read 级 tag**（如 `Q_snp` / `Q_len` / `Q_rep`），零格式风险。
- **BAM 质量字符串保留 `Q_agg`**，保证下游 (samtools / 现有流程) 不变。

### 5.4 评估（分通道校准，而非单标量对齐）

- `Q_snp` 是否预测"错配率"？`Q_len` 是否预测"indel 率"？分通道各出：
  reliability curve（校准）、AUC、`Q≥k` 处的 Precision/Recall。
- **均聚单列**：`Q_rep` 对 copy 数正确率（以 copy 数分布作真值）。
- **核心验证假设**：在 `Q_agg` 相同的位点里，`Q_len` 能否**分离出"低质是因为 indel 而非 SNP"**。

### 5.5 落地步骤（小步验证）

1. **诊断（先证明问题）**：按 (是否均聚) × Q_agg 分箱，画 SNP 占比 vs indel 占比；
   预期均聚区 indel 主导、且同 Q_agg 下误差构成显著不同。
2. **拟合** `P_s(N,ctx)` / `P_len(N,重复,unit)`（用已知真值 / 高 N 自比对作监督）。
3. **写回** `Q_snp / Q_len / Q_rep`（位点旁路表 + 可选 read tag），BAM 保留 `Q_agg`。
4. **评估** per-channel 校准 + 均聚 copy 数正确率。
5. （可选升级）误差强依赖上下文 → 上 **Nanopolish 式 HMM / DeepConsensus 式 gap-aware learned model**（策略 D）。

---

## 6. 开放问题 / 风险

- **参考/真值依赖**：`ins/del` 是**相对参考**计的；参考本身若 copy 数不同，"indel"会混入真实变异。
  需区分"测序长度误差" vs "真实 STR 长度多态"（两者都要识别）。
- **均聚偏置是系统性的**：随机误差模型会**高估**长均聚质量；`P_len` 必须含不随 N 衰减的偏置项。
- **有效 N**：pass 间同源噪声 → 有效 N < 标称 N，质量随 N 的下降可能比理想 `1/√N` 慢；
  用经验曲线而非理想组合模型更稳。
- **ins 与 del 是否再分**：STR 决策里 ±1 copy 常对等，`Q_len` 合并即可；
  仅当"多一个"和"少一个"下游后果不同（如计数）时，才拆 `Q_ins / Q_del`。
- **兼容成本**：任何新通道都要保留 `Q_agg` 与单标量路径，避免破坏现有下游。

---

## 7. 参考资料（可核实）

### 核心方法
- DeepConsensus（PacBio/Google，gap-aware transformer 纠错 CCS）：
  - 预印本：[biorxiv 10.1101/2021.08.31.458403](https://www.biorxiv.org/content/10.1101/2021.08.31.458403v1.full)
  - 发表：[Nat. Biotechnol. 2023; PubMed 36050551](https://pubmed.ncbi.nlm.nih.gov/36050551/)
  - 代码：[google/deepconsensus](https://github.com/google/deepconsensus)
- Nanopolish：[github.com/nanopolish/nanopolish](https://github.com/nanopolish/nanopolish)；
  发射模型源码 [nanopolish_emissions.h](https://github.com/sjackman/nanopolish/blob/master/nanopolish_emissions.h)
- medaka：[nanoporetech/medaka](https://github.com/nanoporetech/medaka)
- PacBio `ccs`：[PacificBiosciences/ccs](https://github.com/PacificBiosciences/ccs)；FAQ [ccs.how](https://ccs.how/faq/mode-all.html)
- CCS 误差分析：[Error analysis of the PacBio sequencing CCS reads, Int. J. Biostatistics](https://www.x-mol.com/paper/1655979466304081920)

### 变异检测 / 质量校准
- GATK BQSR：[官方文档](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR)
- BAQ（indel-aware baseQ）：Li, [Improving SNP discovery by base alignment quality](https://academic.oup.com/bioinformatics/article/27/8/1157/227268)
  (Bioinformatics 2011, doi:10.1093/bioinformatics/btr076；[PubMed 21320865](https://pubmed.ncbi.nlm.nih.gov/21320865/))
- DeepVariant：[biorxiv 10.1101/092890](https://www.biorxiv.org/content/10.1101/092890v5.full)；[postprocess_variants.py（P(NON_REF)/P(ALL)）](https://github.com/google/deepvariant/blob/master/deepvariant/postprocess_variants.py)

### consensus 工具与 indel
- samtools consensus 高 indel 平台改进：jkbonfield, [PR #1733](https://github.com/samtools/samtools/pull/1733)
- bcftools indel calling 改进：[PR #1679](https://github.com/samtools/bcftools/pull/1679)、[PR #2045](https://github.com/samtools/bcftools/pull/2045)

### STR 基因分型
- GangSTR 工具对比：[F1000Research, 9:200](https://f1000research.com/articles/9-200/v1)
- STR 工具横向对比：[Genome-wide STR Truth Set](https://broadinstitute.github.io/str-truth-set/html/tool_comparison_viewer.html)
- TR 基因分型综述：[Translational Psychiatry 2023](https://www.nature.com/articles/s41398-023-02689-8)

### 通用
- Phred 质量定义：[drive5 quality_score](https://www.drive5.com/usearch/manual/quality_score.html)
- NGS 均聚 indel 与错误综述：[PMC3438994](https://pmc.ncbi.nlm.nih.gov/articles/PMC3438994/)
