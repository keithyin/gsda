# 需求文档：SUSPICIOUS_SITES 新旧模型对比 —— 分析方法

- 日期：2026-09-03
- 目的：给定同一源数据的两组 consensus 输出（旧模型 vs 新模型），定义一套可复用、
  可回归的对比方法，量化新模型相对旧模型在**降噪**与**敏感性（sensitivity）**上的差异。
- 产出：一份结构化的对比报告（位点数、样本数、逐位点分类），而非一次性结论。

## 1. 输入

| 组 | 路径 |
|---|---|
| NEW（新模型） | `/data1/ccs_data/20260901-yiqiao/20260831_240601Y0014_Run0005_called_barcode_v4_new_model/Consensus/Consensus/Suspious_sites` |
| OLD（旧模型） | `/data1/ccs_data/20260901-yiqiao/20260831_240601Y0014_Run0005_called_barcode_v4/Consensus/Consensus/Suspious_sites` |

约定：
- 同一份源数据，仅算法不同 → 差异可归因于算法。
- 每样本产出 `<name>.consensus.suspicious_sites.{txt,csv}`；`.txt` 为 tab 分隔、
  `.csv` 为逗号分隔，二者是**同一份数据的两种格式**（仅数字精度差异）。**计数只取其一**，不叠加。
- 样本名格式：`Adaptor-barcode<N>-<M>.consensus.suspicious_sites`，`<M>`（0/1）为同一 barcode
  的两个 cluster/侧，视为**独立样本**分别统计。

### 字段（18 列）
`Clust, Position, Consensus, Depth, A_num, C_num, G_num, T_num, Del_num, Ins_num,
Top1_base, Freq_of_top1_base, Top2_base, Freq_of_top2_base, Quality, Homo_num,
Sus_types, insert_seq`

## 2. 分析维度（要回答的问题）

对比方法需覆盖以下维度，每一项都给出 NEW / OLD 两个数值及差值：

### D1. 位点总量
所有样本位点数之和（数据行数，去掉 header）。→ 衡量整体爆点多少。

### D2. 命中样本数
位点数 > 0 的样本文件数；并列出**仅一组命中**的样本（NEW-only / OLD-only）。→ 衡量
影响面（多少个样本被"爆出"）。

### D3. 逐位点分类（核心）
以 `(样本, Position)` 为键，将位点差异分为三类，逐类列出明细：
- **CLEANED**：OLD 有、NEW 无 —— 变干净（新模型少报）。
- **NEW_SUS**：NEW 有、OLD 无 —— 新爆（新模型多报）。
- **SHIFTED**：同一物理区域的 indel 集合坐标整体偏移（位置编号变化、深度/长度近似、
  类型一致）—— 坐标重定位，**不视为净增减**，单独标注以免误判。

  判 SHIFTED 的经验规则：同一样本内，CLEANED 与 NEW_SUS 位点数接近、`Sus_types` 相同、
  位置呈连续区间且深度量级一致 → 判定为重定位，从净增减中剔除。

### D4. 位点属性画像（对每个差异位点记录，用于归因）
对每个 CLEANED / NEW_SUS 位点，记录：`Consensus / Depth / Top1(频率) / Top2(频率) /
Homo_num / Sus_types / insert_seq 前缀`。并据此对每个差异位点给出**归因标签**：
- 低比例 SNP：高 Depth、Top1≈100%、Top2 为单碱基低频（如 ~5%）→ 关注是否漏检真实低 VAF。
- 超低深度 indel：Depth 极低、Top2 为 `-`/单事件 → 多为噪声，降噪合理。
- 重复/同聚区 indel：Homo_num 高、位置邻近成串 → 坐标易漂移，对应 SHIFTED。

## 3. 输出格式

报告按样本组织，每个样本给出：
```
### <sample>   OLD=<n> NEW=<m>  cleaned=<a> newly_susp=<b> shifted=<c>
  [CLEANED]    pos=... Cons=... Depth=... Top1=..(..) Top2=..(..) Types=..  ins=..   #归因标签
  [NEW_SUS]    ...
  [SHIFTED]    OLD区间 a–b (D=..) → NEW区间 c–d (D=..)
```
末尾给出 D1/D2 汇总与 CLEANED/NEW_SUS 的归因标签分布。

## 4. 验收标准
- [ ] 方法对任意两组 `Suspicious_sites` 目录可直接运行（参数化 NEW / OLD 路径）。
- [ ] 输出含 D1、D2、D3、D4 全部内容，逐位点可追溯到原始字段。
- [ ] CLEANED / NEW_SUS / SHIFTED 三类互斥且完整覆盖所有差异位点（无遗漏、无重叠）。
- [ ] 差异归因（低比例 SNP / 超低深度 indel / 重复区漂移）有明确判据，不靠人工目测。
- [ ] 可固化为脚本（Python 优先），支持后续 run 一键回归对比。

## 5. 附：统计口径
- 位点计数 = 数据行数（去 header），空行忽略。
- 样本计数 = 位点数 > 0 的样本文件数。
- 每组只用一种分隔格式（txt 或 csv）。
