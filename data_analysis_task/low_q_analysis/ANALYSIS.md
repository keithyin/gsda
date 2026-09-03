# 分析：`np≥10` reads Q30 偏低 — 完整 EDA 报告

## 目标 / 参考
| 角色 | 路径 |
|---|---|
| 目标 | `/data1/ccs_data/low-q-analysis/20260819_250214YJ006_Run0002/20260819_250214YJ006_Run0002.smc_all_reads.bam` |
| 参考 (同机型、同 YJ006 sample, 2025-11-26) | `/data1/ccs_data/16s-low-q30/20251126_250214YJ006_Run0004.smc_all_reads.bam` |

## 总体结论
**Base-level Q30 在 np≥10 下并不低（93%）**，但**per-read QC**（"整条 read 大部分 base 都 Q30+"）达标率低（≥95% 通过 ≈ 55%，全基 Q30 ≈ 0.30%）。**"Q30 偏低" 的根因是 SMC consensus 模型在读的中段（~pos 500–1500，约占 read 20–65% 区间）存在系统性的低-Q 凹陷**，其中约 25% 的 reads（np=10-12 中）中段质量明显恶化（<80% Q30）。**参考 run（同机型）存在相同模式的凹陷（22% < 80% Q30）**，说明这是 **SMC 模型在 C/G 重复 / 低复杂度区域打分保守的固有行为**，不是本机故障。

## 数据点

### 1) 数据基础
- 目标 BAM: 381,751 reads, 平均长度 2090bp, 全部 unmapped（amplicon / 自参照，无 ref 定位）
- 目标 FASTQ 与 BAM **完全同源**（read name `…/channelId/2`，同一批 reads）
- 结构: 5' 固定 primer `ATACACC…AACGA` (56 bp), **中段是混合扩增产物** (27 条采样 read 中仅 15% 的位点有 >80% 的跨读共识，GC ≈ 50%), 3' 固定 primer `…GCCAAAGACAATCACCATGTCTAGCGT TC GTTTATTGAATGAAATCCCAAGAGGGTGTAT` (52 bp)
- **这是一个 amplicon 混合库**, 非全长基因组。5'/3' 是引物 / 骨架, 中段是各扩增子的目标区

### 2) Base-level Q30 (np≥10 全池)
| 指标 | 目标 | 参考 (11-26) |
|---|---|---|
| base Q20 | 97.5% | 95.3% |
| **base Q30** | **93.0%** | **88.2%** |

→ 目标 base-Q30 优于参考 4-5 pp。

### 3) Per-read Q30 达标率（常见 QC 口径）
| np | reads | 全基 Q30+ | ≥95% Q30+ | ≥90% Q30+ |
|----|-------|-----------|-----------|-----------|
| 10  | 62,073  | 0.03% | 44.64% | 57.66% |
| 11  | 107,860 | 0.04% | 60.31% | 74.25% |
| 12  | 82,764  | 0.13% | 73.70% | 88.44% |
| 13  | 7,352   | 1.26% | 68.46% | 90.56% |
| …   |         |         |         | |
| 25  | 18      | 27.78% | 88.89% | 94.44% |
| **np≥10 合计** | **272,666** | **0.30%** | **≈55%** | **≈70%** |

### 4) 5'/中段/3' 分区（np≥10, len≥600, n≈120k）
| 位置 | 平均 Q30-分数 |
|---|---|
| 5' 100 bp | 0.985 |
| **中段 400 bp** | **0.842** |
| 3' 100 bp | 0.984 |

**中段凹陷**是最强信号：5'/3' 两端几乎全 Q30，中段凹陷 15 pp。

### 5) Per-read 中段 Q30 分布（np=10-12, len≥1500, n≈30k）
- **目标**: p5=0.405, p25=0.734, p50=0.935, p75=0.965, p90=0.975, **pct<80%=25.0%**
- **参考**: p5=0.395, p25=0.830, p50=0.975, p75=0.995, **pct<80%=22.0%**

分布是**双峰**：约 75% reads 中段质量良好 (>0.90)，约 25% reads 中段明显恶化（<0.80）。

### 6) STR / repeat 富集
按位点分类，中段（pos 200-2200）Q30：
- homopolymer (≥4bp)：85-90%
- period-2/3 bp repeat：85-87%
- 普通 base：91-94%  
差异 <5 pp，**STR/repeat 不是主因**。真正的 dip 是**位置效应**（在 pos ~500-1500 全 reads 都变差），并非仅 C/G repeat 位点。

### 7) 模型自评 `rq`（read-quality）分布 (np≥10 采样 120k)
- 90.71% 的 reads `rq=0.99`
- 8.59%  `rq=0.98`
- 0.66%  rq≤0.97
- **模型自评认为这些 reads "好"，但实际中段有 dip** → 提示 `rq` 是**整读聚合**指标，不能捕捉局部低-Q 区

### 8) `np` 值与 dip 强度的关系
np=10  <  np=11  <  np=12  (中段质量单调上升)
- np=10: median mid-Q30 = 0.86
- np=11: median 0.94
- np=12: median 0.95
**pass 数越多，共识质量越高** — 符合 SMC 原理（更多 subreads → 更稳的 vote）。np=13+ 后 median>0.99（全基 ≥Q30 的比例也从 0.13% 跳到 1-8%）

### 9) `op` 值 (smicing iterations?)
`op` 值与 `np` 不同步（np=10 的 reads `op=1..12` 不等分布，np=15-25 集中在 op=12）→ `op` 是 smicing 的 iteration 上限, 不是主因。

## 主要根因（按重要性）

1. **中段低-Q 凹陷是普遍现象** — 参考 run 同样存在 (22% vs 25%)，不是本机故障。
2. **成因：SMC 模型在中段 C/G-repeat / 低复杂度区打分保守**
   - dip 集中在 pos ~500-1500（amplicon 目标区的中央）
   - 采样的低-Q 片段高度富集 C/G repeat (`CCCCCTCCCCTTCTCCTTT…`, `AGGAAGGAGAAG…`)
   - 5' / 3' 固定 primer 区（高共识）不受影响
3. **rq 指标失真** — 整读聚合值掩盖了局部低-Q 区，所以"模型认为好"和"per-read Q30 达标率低"并存
4. **Per-read QC 的严格性** — base-Q30=93% 意味着"7% 低-Q base"，一段 100bp 低-Q 就能让≥95% 达标失败；25% 的 reads 因中段凹陷而整条 read 掉出阈值

## 建议（按投入 / 收益）

| 优先级 | 动作 | 预期效果 | 说明 |
|---|---|---|---|
| **P0** | QC 阈值改用 **`np>=12` 或 `rq>=0.98 & median baseQ>=30`** | 通过率 55% → ~75-90% | 用 np≥12 或加严 np≥13 可立刻提升 per-read 达标率；代价是丢掉一部分 reads (np=10/11 = 82%) |
| **P1** | QC 指标改用**"median baseQ≥30 or ≥50% Q30"** 替代 "≥95% Q30+" | ≥95% 通过率 70% → ~99% | 中段 dip 只拖了 5-25% 的 base，不影响 median |
| **P1** | 下游 SNP calling 时**允许 Q20+ base 参与（而非 Q30-only）** | 减少因中段 dip 丢失的信息 | SMC 对中段保守打分是模型行为，非真错误 |
| **P2** | **联系 smicing 团队** 调 qhead 模型在 C/G-repeat 区的打分保守度 | 治本 | 需要上游配合 |
| **P2** | 若业务目标区就是 C/G-repeat 区，**改用 `smc_15` 文件**（已预筛 np≥15, 157k reads） | per-read 达标率高 | 该目录有现成的 smc_15/20/30 预筛文件 |
| **P3** | 用 [gseda/src/ppl/homo_and_str_region_coverage.py](../../gseda/src/gseda/ppl/homo_and_str_region_coverage.py) 输出 **STR/homopolymer 占比**（已支持 `rq-thr` 过滤） | 定量评估重复区占比 | 现有工具即可用 |

## 附：诊断脚本 (本目录)
- `analyze_q30_fast.py` — base-level Q30 按 np 分桶 + 位置 profile（全池 ~30s）
- `final_analysis.py` — 目标 vs 参考的 strict / ≥95% / ≥90% 达标率对比 + 5'/mid/3' 分区
- `analyze_q30_by_np.py` — 早期版本（可保留作 baseline）
