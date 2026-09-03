#!/usr/bin/env python3
r"""
Re-insertion 阈值学习模型 — 基于数据学习 poly-N 缺失碱基的插入阈值。

背景:
    smc_all_reads.bam 里每条 read 带 `di` tag，格式:
        "pos,base,support_ratio,depth;pos,base,support_ratio,depth;..."
    每个条目标记 consensus 序列上一个可疑的缺失位点。当该位点落在 poly-N 区
    且 support_ratio 超过阈值时，`smc_bam_post_process` 会在该位置插入 base。

    当前阈值是手工设定的 (base, repeat) → support_ratio 查表。本工具用真实
    数据**学习**每个 (base, repeat) 格子的最优阈值 (F1-max)，输出与
    `smc_bam_post_process` 完全兼容的阈值 JSON。

Ground truth (label):
    把 consensus 序列比对到 reference，以 ref 为准:
        - ref 该位点有碱基而 consensus 缺失 (ref run > consensus run) → emit(1)
        - ref 也没有 → gap(0)
    多碱基插入: 记录 insert_cnt = ref_run - cons_run (>1 完整处理)。

用法:
    reinsert_model input.smc_all_reads.bam reference.fa [--aligned-bam x.bam]
        [--max-records N] [--val-fraction 0.2] [--output-dir .] [--no-plot]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pysam

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

_BASES = "ACGT"
_SEP = "=" * 70

# minimum samples per (base, hp) cell for meaningful threshold
_MIN_GROUP = 10
# minimum positive (emit) samples in cell to avoid tiny-sample F1=1.0
_MIN_POS = 5

# --- 分格分策略参数 ---
# 长 poly-N: 重复数 >= 该值时视为强信号 (poly-N 越长 SMC 越易漏 call, 低阈值多插)
_LONG_POLY_MIN_REPEAT = 5
# 长 poly-N: 且 emit 占比 >= 该值才触发低阈值 (排除 rep>=5 但 ref 无长 poly 的纯 gap 格子)
_LONG_POLY_EMIT_FRAC = 0.50
# 长 poly-N 的固定低阈值 (support_ratio > thr 即插)
_LONG_POLY_FIXED_THR = 0.20
# 可分格子: train F1 需达到该值 (过低说明 support_ratio 不可分, 学的是噪声)
_MIN_SEPARABLE_F1 = 0.40
# 可分格子: 相对 always-emit baseline 的 F1 提升需达到该值
_MIN_SEPARABLE_LIFT = 0.10
# 可分格子: 正样本数下限, 防止 tiny-sample F1 过拟合 (如 6 个正样本撑起的 F1=0.8)
_MIN_SEPARABLE_POS = 50


# ---------------------------------------------------------------------------
# PolyN 语义：与 /usr/bin/smc_bam_post_process 完全一致
# ---------------------------------------------------------------------------

def find_poly_n_repeat(seq: str, pos: int, base: str) -> Optional[int]:
    """Find consecutive `base` count started at position `pos`.

    与 smc_bam_post_process.find_poly_n_repeat 完全一致，保证训练/推理同语义。
    """
    if pos >= len(seq):
        return None
    if seq[pos] != base:
        return 0
    right = pos
    while right < len(seq) and seq[right] == base:
        right += 1
    return right - pos


class PolyNThresholdMap:
    """解析阈值 JSON，与 smc_bam_post_process.PolyNThresholdMap 兼容。

    仅用于验证生成的 JSON 可被现有工具解析。
    """

    def __init__(self, flat: Dict[Tuple[str, int], float]):
        self.flat = flat

    @classmethod
    def from_dict(cls, raw: Dict[str, float]) -> "PolyNThresholdMap":
        flat: Dict[Tuple[str, int], float] = {}
        for key, val in raw.items():
            try:
                rep = int(key.split(")")[1])
                base = key[1]
            except (ValueError, IndexError):
                continue
            flat[(base, rep)] = val
        return cls(flat=flat)

    @classmethod
    def from_json(cls, path: str) -> "PolyNThresholdMap":
        with open(path) as fh:
            return cls.from_dict(json.load(fh))

    def lookup(self, base: str, repeat: int) -> Optional[float]:
        return self.flat.get((base, repeat))

    def get_threshold(self, repeat_count: int, base: str) -> Optional[float]:
        """largest <= n 策略，与 smc_bam_post_process 一致。"""
        reps = sorted((r for (b, r), _ in self.flat.items() if b == base))
        if not reps:
            return None
        best_val = 0.6
        for rep in reps:
            if rep <= repeat_count:
                best_val = self.flat[(base, rep)]
            else:
                break
        return best_val


# ---------------------------------------------------------------------------
# 数据加载
# ---------------------------------------------------------------------------

def load_reference(path: str) -> Tuple[str, str]:
    """加载参考序列，返回 (ref_name, ref_seq)。"""
    name = ""
    seq = ""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq += line
    return name, seq


def read_di_tag(record: pysam.AlignedSegment) -> List[Tuple[int, str, float, int]]:
    """解析 di tag: [(pos, base, support_ratio, depth), ...]。"""
    if not record.has_tag("di"):
        return []
    di_string = record.get_tag("di")
    entries = []
    for segment in di_string.split(";"):
        segment = segment.strip()
        if not segment:
            continue
        parts = segment.split(",")
        if len(parts) != 4:
            continue
        try:
            pos = int(parts[0])
            base = parts[1].upper()
            support = float(parts[2])
            depth = int(parts[3])
            if base in _BASES:
                entries.append((pos, base, support, depth))
        except (ValueError, IndexError):
            continue
    return entries


# ---------------------------------------------------------------------------
# CIGAR → ref 坐标映射
# ---------------------------------------------------------------------------

def build_query_to_ref(cigartuples, ref_start: int) -> Dict[int, int]:
    """构建 query-pos → ref-coord 的映射。

    处理 M/=/X (双端消耗)、I (query 端)、D (ref 端)、S (softclip)。
    Insertion (I) 位置的 query base 映射到当前 ref 坐标 (ref 没有对应碱基)。
    返回 {qpos: refcoord}；只包含比对到 ref 的 query 位置。
    """
    mapping: Dict[int, int] = {}
    q = 0
    r = ref_start
    for op, length in cigartuples:
        if op in (0, 7, 8):  # M / = / X
            for _ in range(length):
                mapping[q] = r
                q += 1
                r += 1
        elif op == 1:  # I: query only
            q += length
        elif op == 2:  # D: ref only
            r += length
        elif op == 4:  # S: softclip, query only
            q += length
        # op 5 (H) / 6 (P) 不影响坐标
    return mapping


def ref_run_length(ref_seq: str, coord: int, base: str) -> int:
    """ref 在 coord 处向右连续 base 的个数。"""
    if coord < 0 or coord >= len(ref_seq):
        return 0
    if ref_seq[coord] != base:
        return 0
    cnt = 0
    while coord < len(ref_seq) and ref_seq[coord] == base:
        cnt += 1
        coord += 1
    return cnt


def run_total_length(seq: str, pos: int, base: str) -> int:
    """pos 处 base 的 run 总长度 (向前 + 向后)。

    Label 提取用 run-total 而非单点 run，是因为 homopolymer 比对会自由选择
    placement：ref 与 consensus 的 run 可能错位对齐，导致单点比较
    (ref[refcoord] == base) 判不出"consensus 比 ref 少碱基"。
    用整条 run 总长度对比 (ref run > consensus run → emit) 天然鲁棒。
    """
    if pos < 0 or pos >= len(seq) or seq[pos] != base:
        return 0
    left = 0
    i = pos
    while i >= 0 and seq[i] == base:
        i -= 1
        left += 1
    right = 0
    i = pos + 1
    while i < len(seq) and seq[i] == base:
        i += 1
        right += 1
    return left + right


@dataclass
class Sample:
    """一条 di 位点的特征 + label。"""

    base: str
    repeat: int          # consensus 中从 pos 向右连续 base 个数
    support_ratio: float # di 中的证据强度
    depth: int           # di 中的 depth（仅诊断，不进模型）
    label: int           # 1=emit, 0=gap
    insert_cnt: int      # ref_run - cons_run (多碱基插入支持)
    pos_frac: float      # pos / len(seq)，诊断用
    qname: str           # 来源 read，诊断用
    refcoord: int        # 映射到的 ref 坐标，诊断用


def extract_samples(
    aligned_bam: str,
    source_bam: str,
    ref_seq: str,
    max_records: Optional[int] = None,
    verbose: bool = True,
) -> Tuple[List[Sample], dict]:
    """从 aligned bam + source bam 提取 di 位点的特征和 label。

    di tag 在 source bam 中，CIGAR/ref 坐标在 aligned bam 中，按 query_name 关联。

    Returns:
        (samples, stats)
        stats: {n_aligned, n_with_di, n_skipped_no_mapping, n_skipped_not_poly,
                n_skipped_edge, n_emit, n_gap}
    """
    # 先读 source bam 的 di tag（内存中保留 query_name → di 列表）
    # 用 dict 存 query_name → [(pos, base, sup, dep)]
    source_di: Dict[str, List[Tuple[int, str, float, int]]] = {}
    with pysam.AlignmentFile(source_bam, "rb", check_sq=False) as af:
        for i, record in enumerate(af.fetch(until_eof=True)):
            if max_records is not None and i >= max_records:
                break
            entries = read_di_tag(record)
            if entries:
                source_di[record.query_name] = entries
    if verbose:
        print(f"  读取 source bam: {len(source_di)} 条带 di 的 read")

    stats = {
        "n_aligned": 0,
        "n_with_di": 0,
        "n_skipped_no_mapping": 0,
        "n_skipped_not_poly": 0,
        "n_skipped_edge": 0,
        "n_emit": 0,
        "n_gap": 0,
    }

    samples: List[Sample] = []
    with pysam.AlignmentFile(aligned_bam, "rb", check_sq=False) as af:
        for record in af.fetch(until_eof=True):
            qname = record.query_name
            if qname not in source_di:
                continue
            stats["n_aligned"] += 1
            entries = source_di[qname]
            stats["n_with_di"] += 1

            if record.flag & 4:  # unmapped
                continue

            cons_seq = record.query_sequence
            if cons_seq is None:
                continue
            qlen = len(cons_seq)
            if qlen == 0:
                continue

            # 处理 secondary/supplementary: 用 primary 的比对
            mapping = build_query_to_ref(record.cigartuples, record.reference_start)

            for (pos, base, support, depth) in entries:
                if pos >= qlen:
                    continue

                repeat = find_poly_n_repeat(cons_seq, pos, base)
                if repeat is None or repeat == 0:
                    stats["n_skipped_not_poly"] += 1
                    continue

                if pos not in mapping:
                    stats["n_skipped_no_mapping"] += 1
                    continue

                refcoord = mapping[pos]
                # Label: 比较整条 run 总长度 (鲁棒于 homopolymer placement)
                rrun = run_total_length(ref_seq, refcoord, base)
                crun = run_total_length(cons_seq, pos, base)
                if crun == 0:
                    stats["n_skipped_not_poly"] += 1
                    continue

                if rrun > crun:
                    label = 1
                    insert_cnt = rrun - crun
                    stats["n_emit"] += 1
                else:
                    label = 0
                    insert_cnt = 0
                    stats["n_gap"] += 1

                samples.append(
                    Sample(
                        base=base,
                        repeat=repeat,
                        support_ratio=support,
                        depth=depth,
                        label=label,
                        insert_cnt=insert_cnt,
                        pos_frac=pos / qlen if qlen else 0.0,
                        qname=qname,
                        refcoord=refcoord,
                    )
                )

    return samples, stats


# ---------------------------------------------------------------------------
# 阈值搜索 (F1-max, 按 (base, repeat) 格子)
# ---------------------------------------------------------------------------

def _precision_recall(tp, fp, fn):
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    return prec, rec


def search_cell_threshold(samples: List[Sample]) -> Optional[dict]:
    """在单个 (base, repeat) 格子内搜索最优 support_ratio 阈值 (F1-max)。

    预测规则: support_ratio > thr → emit(1)。与 smc_bam_post_process 的
    'support_ratio <= threshold → skip' 逻辑一致。
    """
    n_pos = sum(1 for s in samples if s.label == 1)
    n_neg = len(samples) - n_pos
    if n_pos == 0 or n_neg == 0 or n_pos < _MIN_POS:
        return None

    # 只搜索支持度取值 (含 0.0 和 1.0)，与现有工具的边界一致
    ratios = sorted(set(s.support_ratio for s in samples))

    best_f1 = -1.0
    best_thr = None
    best_tp = best_fp = best_fn = 0

    for thr in ratios:
        tp = sum(1 for s in samples if s.support_ratio > thr and s.label == 1)
        fp = sum(1 for s in samples if s.support_ratio > thr and s.label == 0)
        fn = n_pos - tp
        prec, rec = _precision_recall(tp, fp, fn)
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
        if f1 > best_f1:
            best_f1 = f1
            best_thr = thr
            best_tp, best_fp, best_fn = tp, fp, fn

    if best_thr is None:
        return None

    return {
        "ratio_threshold": best_thr,
        "f1": round(best_f1, 4),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "tp": best_tp,
        "fp": best_fp,
        "fn": best_fn,
    }


def group_by_cell(samples: List[Sample]) -> Dict[Tuple[str, int], List[Sample]]:
    """按 (base, repeat) 分组。"""
    cells: Dict[Tuple[str, int], List[Sample]] = {}
    for s in samples:
        cells.setdefault((s.base, s.repeat), []).append(s)
    return cells


def _always_emit_f1(cell: List[Sample]) -> float:
    """always-emit 策略的 F1 (全部预测 emit)。"""
    n = len(cell)
    n_pos = sum(1 for s in cell if s.label == 1)
    if n == 0:
        return 0.0
    prec = n_pos / n
    return 2 * prec / (prec + 1.0) if prec > 0 else 0.0


def find_thresholds(samples: List[Sample]) -> Tuple[Dict[Tuple[str, int], dict], dict]:
    """分格分策略搜索阈值。

    每个 (base, repeat) 格子按数据可分性分策略。核心原则：**只有 support_ratio
    真正可分的格子才产出学习阈值**；其余保守处理 (推理端 get_threshold 兜底)。

    策略:
    1. **长 poly-N 强信号** (emit 占比 >= 50% 且 support_ratio 确实可分):
       长 poly-T 长尾这类格子 emit 几乎 100%, support_ratio 上 F1-max 达到或
       接近 always-emit 的最优 (lift 大或全 emit), 用固定低阈值多插。
    2. **可分格子**: F1-max 相对 always-emit baseline 显著提升 (lift 大)
       → F1-max 学阈值。
    3. **不可分格子** (emit 与 gap 的 support_ratio 重叠, F1-max 不优于 baseline):
       不产出阈值。这类格子无论怎么切 support_ratio 都插不准, 学习阈值只是
       在噪声上拟合; 推理端由 get_threshold 的 "largest <= n" 兜底保守处理。

    Returns:
        (lookup, meta): lookup[(base,repeat)] -> {ratio_threshold, f1, ...}
            meta[(base,repeat)] -> {"strategy": "long_poly_fixed"
                                    |"separable"|"inseparable", "reason": str}
    """
    lookup: Dict[Tuple[str, int], dict] = {}
    meta: Dict[Tuple[str, int], dict] = {}
    cells = group_by_cell(samples)

    for key in sorted(cells):
        base, repeat = key
        cell = cells[key]
        if len(cell) < _MIN_GROUP:
            meta[key] = {"strategy": "inseparable",
                         "reason": f"样本不足 (<{_MIN_GROUP})"}
            continue

        n_pos = sum(1 for s in cell if s.label == 1)
        n_neg = len(cell) - n_pos
        emit_frac = n_pos / len(cell)
        base_f1 = _always_emit_f1(cell)

        # --- 长 poly-N 强信号: 固定低阈值多插 ---
        # 条件: 重复数 >= _LONG_POLY_MIN_REPEAT 且 emit 占比足够高。
        #   poly-N 越长 SMC 越易漏 call, 长尾格子 emit 占比接近 100%,
        #   "有任何支持就插"是可靠策略。用 repeat 而非 emit 占比做主判据,
        #   避免 run-total label 的假 emit (A2/C3/C4 rep 2-4) 触发低阈值;
        #   同时要求 emit 占比不低, 排除 rep>=5 但 ref 无长 poly (纯 gap,
        #   如 C5/G6 emit=0) 的格子。
        if repeat >= _LONG_POLY_MIN_REPEAT and emit_frac >= _LONG_POLY_EMIT_FRAC:
            lookup[key] = {
                "ratio_threshold": _LONG_POLY_FIXED_THR,
                "f1": round(base_f1, 4),
                "n_positive": n_pos,
                "n_negative": n_neg,
                "tp": n_pos,
                "fp": n_neg,
                "fn": 0,
            }
            meta[key] = {
                "strategy": "long_poly_fixed",
                "reason": (f"长 poly-N (repeat={repeat}>=5, emit={emit_frac:.1%}), "
                           f"固定低阈值 {_LONG_POLY_FIXED_THR}"),
            }
            continue

        # --- 常规格子: 可分性判断 ---
        if n_pos < _MIN_POS:
            meta[key] = {"strategy": "inseparable",
                         "reason": f"emit 样本不足 (<{_MIN_POS})"}
            continue
        if n_neg == 0:
            meta[key] = {"strategy": "inseparable", "reason": "无 gap 样本"}
            continue

        result = search_cell_threshold(cell)
        if result is None:
            meta[key] = {"strategy": "inseparable", "reason": "F1-max 无有效阈值"}
            continue

        lift = result["f1"] - base_f1
        if (n_pos >= _MIN_SEPARABLE_POS
                and result["f1"] >= _MIN_SEPARABLE_F1
                and lift >= _MIN_SEPARABLE_LIFT):
            lookup[key] = result
            meta[key] = {
                "strategy": "separable",
                "reason": (f"train F1={result['f1']:.3f} "
                           f"(baseline always-emit F1={base_f1:.3f}, lift={lift:+.3f})"),
            }
        else:
            meta[key] = {
                "strategy": "inseparable",
                "reason": (f"train F1={result['f1']:.3f} 低于可分阈值 "
                           f"(emit={emit_frac:.1%}, n_pos={n_pos}, "
                           f"baseline F1={base_f1:.3f}, lift={lift:+.3f})"),
            }

    return lookup, meta


# ---------------------------------------------------------------------------
# 导出阈值 JSON (兼容格式 + gap-fill)
# ---------------------------------------------------------------------------

def build_threshold_json(
    lookup: Dict[Tuple[str, int], dict], samples: List[Sample]
) -> Dict[str, float]:
    """把 (base, repeat) 格子阈值转成 {"(A)5": 0.25} 格式。

    只输出策略实际学出的格子。未配置的格子由 smc_bam_post_process 的
    get_threshold "largest <= n" 语义兜底: 对某个 base 未配置的 repeat 取
    已配置的最大 <= n 的阈值。因此 G 只配置 G1/G5 时, G2-4 自动用 G1 的阈值,
    G6+ 自动用 G5 的阈值 —— 无需显式 gap-fill。
    """
    result: Dict[str, float] = {}
    for (base, repeat), entry in sorted(lookup.items()):
        result[f"({base}){repeat}"] = entry["ratio_threshold"]
    return result


# ---------------------------------------------------------------------------
# 评估
# ---------------------------------------------------------------------------

def evaluate(preds: List[int], labels: List[int]) -> dict:
    """混淆矩阵与衍生指标 (正类 = emit=1)。"""
    n = len(labels)
    tp = sum(1 for p, y in zip(preds, labels) if p == 1 and y == 1)
    fp = sum(1 for p, y in zip(preds, labels) if p == 1 and y == 0)
    fn = sum(1 for p, y in zip(preds, labels) if p == 0 and y == 1)
    tn = n - tp - fp - fn
    prec, rec = _precision_recall(tp, fp, fn)
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    acc = (tp + tn) / n if n > 0 else 0.0
    return {
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "precision": prec, "recall": rec, "f1": f1, "accuracy": acc,
    }


def apply_lookup(samples: List[Sample], lookup: Dict[Tuple[str, int], dict]):
    """对每个样本按 (base, repeat) 查表预测。

    Returns:
        (preds, covered): preds[i]=1/0 或 None (无阈值格子); covered[i]=bool
    """
    preds: List[Optional[int]] = []
    covered: List[bool] = []
    for s in samples:
        entry = lookup.get((s.base, s.repeat))
        if entry is None:
            preds.append(None)
            covered.append(False)
        else:
            preds.append(1 if s.support_ratio > entry["ratio_threshold"] else 0)
            covered.append(True)
    return preds, covered


def format_eval_table(name: str, ev: dict) -> str:
    return (
        f"{name:<20} {ev['tp']:>4d} {ev['fp']:>4d} {ev['fn']:>4d} {ev['tn']:>4d} "
        f"{ev['precision']:>9.4f} {ev['recall']:>7.4f} {ev['f1']:>6.4f} "
        f"{ev['accuracy']:>8.4f}"
    )


# ---------------------------------------------------------------------------
# 可视化
# ---------------------------------------------------------------------------

def plot_scatter(
    samples: List[Sample],
    lookup: Dict[Tuple[str, int], dict],
    output_path: str,
):
    """support_ratio × repeat 散点图，emit/gap 着色 + 阈值线。"""
    if not HAS_MPL:
        return
    _BASE_COLORS = {"G": "#E63946", "A": "#219EBC",
                    "C": "#2EC4B6", "T": "#F77F00"}

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_yscale("log", base=2)
    max_repeat = max((s.repeat for s in samples), default=1)
    y_max_log = min(max(5, max_repeat), 100)
    ax.set_ylim(bottom=0.5, top=y_max_log)
    ax.set_yticks([1, 2, 4, 8, 16, 32, 64])

    total = len(samples)
    s_size = max(3, min(8, 200 / max(total, 1) ** 0.4))
    layer_alpha = min(0.5, 0.3 / max(total ** 0.3, 1))

    for base in sorted(_BASE_COLORS):
        indices = [i for i, s in enumerate(samples) if s.base == base]
        if not indices:
            continue
        color = _BASE_COLORS[base]
        for label, alpha_mul in [(1, 1.2), (0, 0.8)]:
            sub = [i for i in indices if samples[i].label == label]
            if not sub:
                continue
            x = [samples[i].support_ratio for i in sub]
            y = [samples[i].repeat for i in sub]
            suffix = "emit" if label else "gap"
            ax.scatter(
                x, y,
                c=color + ("FF" if label else "80"),
                s=s_size, alpha=layer_alpha * alpha_mul,
                label=f"{base}-{suffix} n={len(sub)}",
            )

    if lookup:
        best_key = max(lookup.keys(), key=lambda k: lookup[k]["f1"])
        b, h = best_key
        r_thr = lookup[best_key]["ratio_threshold"]
        ax.axvline(x=r_thr, color="#457B9D", linewidth=1.2, linestyle="--",
                   label=f"({b}, hp={h}): ratio > {r_thr:.3f}")
        metric = (f"F1={lookup[best_key]['f1']:.4f} | 查表阈值: "
                  f"base={b}, hp={h}, ratio > {r_thr:.3f}")
    else:
        metric = ""
        ax.text(0.5, 0.95, "No threshold data", transform=ax.transAxes,
                ha="center", va="top", fontsize=12, color="gray")

    ax.set_title(f"Re-insertion Decision Scatter Plot\n{metric}", fontsize=13)
    ax.set_xlabel("support_ratio (evidence strength)", fontsize=11)
    ax.set_ylabel("poly-N repeat length (log scale)", fontsize=11)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    ax.grid(True, alpha=0.3, which="both")
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 比对
# ---------------------------------------------------------------------------

def align_to_ref(source_bam: str, ref_fasta: str, threads: int, out_prefix: str) -> str:
    """用 gsmm2 把 source_bam 比对到 ref，返回 aligned bam 路径。"""
    aligned_bam = out_prefix + ".bam"
    cmd = [
        "gsmm2", "--threads", str(threads),
        "align",
        "-q", source_bam,
        "-t", ref_fasta,
        "-p", out_prefix,
    ]
    print("  运行比对:", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"gsmm2 align 失败 (exit {proc.returncode})")
    return aligned_bam


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main_cli(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Re-insertion 阈值学习 — 基于 ref 比对学习 poly-N 插入阈值",
    )
    parser.add_argument("input_bam", help="带 di tag 的 SMC consensus bam")
    parser.add_argument("reference", help="参考序列 fasta")
    parser.add_argument("--aligned-bam", default=None,
                        help="已比对到 ref 的 bam (跳过比对步骤)")
    parser.add_argument("--align-threads", type=int, default=32,
                        help="gsmm2 比对线程数 (default: 32)")
    parser.add_argument("--max-records", type=int, default=None,
                        help="最大处理记录数 (default: 全部)")
    parser.add_argument("--val-fraction", type=float, default=0.2,
                        help="验证集比例 (default: 0.2, 在 train 上搜阈值、val 上评估)")
    parser.add_argument("--output-dir", default=".",
                        help="输出目录 (PNG + JSON) (default: 当前目录)")
    parser.add_argument("--no-plot", action="store_true",
                        help="不生成散点图")
    parser.add_argument("--keep-alignment", action="store_true",
                        help="保留中间比对 bam (默认清理)")
    args = parser.parse_args(argv)

    if not (0 <= args.val_fraction < 1):
        print("[ERROR] --val-fraction 必须在 [0, 1) 区间", file=sys.stderr)
        return 1

    os.makedirs(args.output_dir, exist_ok=True)
    ref_name, ref_seq = load_reference(args.reference)
    print(f"参考序列: {ref_name} (len={len(ref_seq)})")

    # Step 1: 比对 (或复用)
    tmpdir = None
    if args.aligned_bam is not None:
        aligned_bam = args.aligned_bam
        print(f"使用已比对 bam: {aligned_bam}")
    else:
        tmpdir = tempfile.mkdtemp(prefix="reinsert_model_")
        out_prefix = os.path.join(tmpdir, "aligned")
        print("步骤 1: 比对到 reference ...")
        aligned_bam = align_to_ref(
            args.input_bam, args.reference, args.align_threads, out_prefix)
    if not os.path.exists(aligned_bam):
        print(f"[ERROR] 未找到比对结果: {aligned_bam}", file=sys.stderr)
        return 1

    # Step 2: 特征 + label 提取
    print("步骤 2: 提取特征与 label ...")
    samples, stats = extract_samples(
        aligned_bam, args.input_bam, ref_seq,
        max_records=args.max_records,
    )
    n = len(samples)
    print(f"  提取 {n} 个 di 位点特征")
    for k, v in stats.items():
        if k.startswith("n_") and k not in ("n_aligned", "n_with_di"):
            print(f"    {k}: {v}")

    if n == 0:
        print("[WARN] 无有效样本，退出。", file=sys.stderr)
        return 0

    n_emit = sum(1 for s in samples if s.label == 1)
    print(f"  label 分布: emit={n_emit} ({100.0*n_emit/n:.2f}%), "
          f"gap={n - n_emit} ({100.0*(n-n_emit)/n:.2f}%)")

    # Step 3: train/val 划分
    if args.val_fraction > 0:
        n_val = max(1, int(n * args.val_fraction))
        # 按记录分块划分，避免同 read 的位点跨 train/val
        # 先按 qname 去重排序，分块
        qnames = sorted({s.qname for s in samples})
        n_q = len(qnames)
        # 按 qname 数量切分
        split_idx = int(n_q * (1 - args.val_fraction))
        train_qnames = set(qnames[:split_idx])
        train_samples = [s for s in samples if s.qname in train_qnames]
        val_samples = [s for s in samples if s.qname not in train_qnames]
        print(f"  train/val 划分: train={len(train_samples)}, val={len(val_samples)} "
              f"({len(train_qnames)}/{n_q - len(train_qnames)} reads)")
    else:
        train_samples = samples
        val_samples = []

    # Step 4: 阈值搜索 (train)
    print("步骤 3: 搜索阈值 (train) ...")
    lookup, threshold_meta = find_thresholds(train_samples)
    if not lookup:
        print("[WARN] 未找到有效阈值 (样本太稀/标签太偏)。", file=sys.stderr)
        if HAS_MPL and not args.no_plot:
            plot_path = os.path.join(args.output_dir, "reinsert_scatter.png")
            plot_scatter(samples, {}, plot_path)
            print(f"散点图已保存至: {plot_path}")
        return 0

    # Step 5: 导出阈值 JSON
    threshold_json = build_threshold_json(lookup, train_samples)
    json_path = os.path.join(args.output_dir, "reinsert_thresholds.json")
    with open(json_path, "w") as fh:
        json.dump(threshold_json, fh, indent=4)
    print(f"\n阈值 JSON 已保存至: {json_path}")
    print(f"  共 {len(threshold_json)} 个格子")

    # 阈值决策依据 (meta) 落盘
    meta_path = os.path.join(args.output_dir, "reinsert_thresholds.meta.json")
    meta_out = {
        f"({base}){repeat}": {
            "strategy": m["strategy"],
            "reason": m["reason"],
            "ratio_threshold": lookup.get((base, repeat), {}).get("ratio_threshold"),
        }
        for (base, repeat), m in sorted(threshold_meta.items())
    }
    with open(meta_path, "w") as fh:
        json.dump(meta_out, fh, indent=4, ensure_ascii=False)
    print(f"阈值决策依据已保存至: {meta_path}")

    # 打印查表阈值表
    print(f"\n{_SEP}")
    print("--- 按 (碱基, repeat) 查表阈值 (分格分策略) ---")
    print(f"{'Base':<6} {'repeat':>7s} {'ratio_thr':>10s} {'F1':>6s} "
          f"{'+pos':>6s} {'-neg':>6s} {'TP':>4s} {'FP':>4s} {'FN':>4s} 策略")
    for (base, repeat) in sorted(lookup.keys()):
        t = lookup[(base, repeat)]
        strat = threshold_meta.get((base, repeat), {}).get("strategy", "?")
        print(f"{base:<6} {repeat:>7d} {t['ratio_threshold']:>10.4f} "
              f"{t['f1']:>6.4f} {t['n_positive']:>6d} {t['n_negative']:>6d} "
              f"{t['tp']:>4d} {t['fp']:>4d} {t['fn']:>4d}  {strat}")

    # 打印被判定为不可分的格子 (推理端会走保守路径)
    n_insep = sum(1 for m in threshold_meta.values()
                  if m["strategy"] == "inseparable")
    if n_insep:
        print(f"\n{_SEP}")
        print(f"--- 不可分格子 (不产出阈值, 推理端 get_threshold 保守兜底) --- "
              f"共 {n_insep} 个")
        for key in sorted(threshold_meta):
            if threshold_meta[key]["strategy"] == "inseparable":
                print(f"  {key[0]}, repeat={key[1]}: {threshold_meta[key]['reason']}")

    # Step 6: 评估 (val)
    if val_samples:
        print(f"\n{_SEP}")
        print(f"--- 评估 (val, n={len(val_samples)}) ---")
        val_preds, val_covered = apply_lookup(val_samples, lookup)
        val_covered_idx = [i for i, c in enumerate(val_covered) if c]
        n_cov = len(val_covered_idx)
        total = len(val_samples)
        if n_cov == 0:
            print(f"coverage = 0/{total} (0.0%)，无有阈值格子，跳过评估。")
        else:
            cov_labels = [val_samples[i].label for i in val_covered_idx]
            cov_preds = [val_preds[i] for i in val_covered_idx]
            model_eval = evaluate(cov_preds, cov_labels)
            base_eval = evaluate([0] * n_cov, cov_labels)  # baseline: 全部不插

            print(f"{'方法':<20} {'TP':>4s} {'FP':>4s} {'FN':>4s} {'TN':>4s} "
                  f"{'Precision':>9s} {'Recall':>7s} {'F1':>6s} {'Accuracy':>8s}")
            print(format_eval_table("baseline (SMC全gap)", base_eval))
            print(format_eval_table("reinsert 模型", model_eval))
            print(f"coverage = {n_cov}/{total} ({100.0*n_cov/total:.1f}%) "
                  f"[有阈值格子占全部样本的比例]")

            # 分格评估
            print(f"\n{_SEP}")
            print("--- 分格评估 (val): 模型 vs baseline ---")
            print(f"{'Base':<6} {'repeat':>7s} {'N':>7s} {'方法':<20} "
                  f"{'TP':>4s} {'FP':>4s} {'FN':>4s} {'TN':>4s} "
                  f"{'Precision':>9s} {'Recall':>7s} {'F1':>6s} {'Accuracy':>8s}")
            cell_groups: Dict[Tuple[str, int], List[int]] = {}
            for i in val_covered_idx:
                s = val_samples[i]
                key = (s.base, s.repeat)
                cell_groups.setdefault(key, []).append(i)
            for key in sorted(cell_groups):
                idx = cell_groups[key]
                c_labels = [val_samples[i].label for i in idx]
                c_preds = [val_preds[i] for i in idx]
                m = evaluate(c_preds, c_labels)
                b = evaluate([0] * len(idx), c_labels)
                print(f"{key[0]:<6} {key[1]:>7d} {len(idx):>7d} "
                      + format_eval_table("baseline (SMC全gap)", b))
                print(f"{'':<6} {'':>7} {'':>7} "
                      + format_eval_table("reinsert 模型", m))

    # Step 7: 可视化
    if HAS_MPL and not args.no_plot:
        plot_path = os.path.join(args.output_dir, "reinsert_scatter.png")
        plot_scatter(samples, lookup, plot_path)
        print(f"\n散点图已保存至: {plot_path}")

    # 清理
    if tmpdir is not None and not args.keep_alignment:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main_cli())
