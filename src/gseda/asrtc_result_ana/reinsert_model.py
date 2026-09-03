#!/usr/bin/env python3
"""Re-insertion 决策模型 — SMC gap 位点预测：emit base or keep gap。

纯阈值规则模型，不依赖机器学习框架。

特征（每个 SMC gap 位点）:
    majority_base_ratio: subreads 中票数最多的非-gap base / 覆盖该位点的 subreads 数量（depth）
    smc_hp_length:     SMC consensus 在该位点左右连续相同碱基的个数
                        (左右碱基相同则取和, 不同则以右边为准)

Label（基于 reference，仅用于生成训练数据/评估）:
    ref[pos] 是 ACGT → emit (1),  SMC 漏 call
    ref[pos] 是 '-'  → gap  (0),  SMC 正确 del

用法:
    python reinsert_model.py input.asrtc.txt [--max-records N] [--output-dir ./figures]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VALID_BASES = frozenset("ACGT")
_SEP = "=" * 70


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_records(filepath: str, max_records: int | None) -> list[dict]:
    """读取 JSONL 文件，返回 msa_seqs 列表。"""
    records: list[dict] = []
    with open(filepath, "r") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            if max_records is not None and i >= max_records:
                break
            records.append(json.loads(line))
    return records


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

def _compute_smc_hp(smc: str, pos: int) -> tuple[int, str]:
    """计算 SMC consensus 在 pos 位的 homopolymer 长度和碱基类型。

    从 pos 向左右扫描（跳过 '-'），找紧邻两侧各自连续相同碱基的个数。
    若左右碱基相同 → 取和；若不同 → 以右边为准；若一边有值另一边无 → 用有值的那边。

    Returns:
        (hp_length, base_char)  e.g. (6, 'G') or (3, 'G') or (0, '')
    """
    # --- Scan right: find first non-gap character ---
    r_char = ""
    r_idx = pos + 1
    while r_idx < len(smc) and smc[r_idx] == "-":
        r_idx += 1
    if r_idx < len(smc) and smc[r_idx] in _VALID_BASES:
        r_char = smc[r_idx]

    # --- Scan left: find first non-gap character ---
    l_char = ""
    l_idx = pos - 1
    while l_idx >= 0 and smc[l_idx] == "-":
        l_idx -= 1
    if l_idx >= 0 and smc[l_idx] in _VALID_BASES:
        l_char = smc[l_idx]

    # --- Compute runs for each character ---
    def _run_right(start: int, target: str) -> int:
        length = 0
        i = start
        while i < len(smc):
            if smc[i] == "-":
                i += 1
                continue
            if smc[i] != target:
                break
            length += 1
            i += 1
        return length

    def _run_left(start: int, target: str) -> int:
        length = 0
        i = start
        while i >= 0:
            if smc[i] == "-":
                i -= 1
                continue
            if smc[i] != target:
                break
            length += 1
            i -= 1
        return length

    if l_char and r_char:
        if l_char == r_char:
            # Same base on both sides — combine
            left_run = _run_left(pos, l_char)
            right_run = _run_right(r_idx, r_char)
            hp_len = left_run + right_run
            return hp_len, l_char
        else:
            # Different bases — use right side
            if r_char:
                right_run = _run_right(r_idx, r_char)
                return right_run, r_char
            elif l_char:
                left_run = _run_left(pos, l_char)
                return left_run, l_char
    elif r_char:
        # Only right side has a base
        right_run = _run_right(r_idx, r_char)
        return right_run, r_char
    elif l_char:
        # Only left side has a base
        left_run = _run_left(pos, l_char)
        return left_run, l_char

    return 0, ""


def _extract_features_one(rec: dict) -> dict:
    """处理单条记录，返回该记录内所有 SMC gap 位点的特征 dict。

    模块级函数：ProcessPoolExecutor 的 worker 需要可 pickle。
    """
    ratios = []
    hp_lens = []
    labels = []
    hp_bases = []

    msa = rec["msa_seqs"]
    smc = msa[1]
    ref = msa[2]
    subreads = msa[3:]

    # Filter out all-gap subreads (useless noise)
    subreads = [s for s in subreads if not all(c == "-" for c in s)]

    n_total_sub = len(subreads)
    if n_total_sub == 0 or len(smc) != len(ref):
        return {
            "ratios": ratios,
            "hp_lens": hp_lens,
            "labels": labels,
            "hp_bases": hp_bases,
        }

    def _subread_covers(subread: str, pos: int) -> bool:
        """判断 subread 是否覆盖给定位点。

        subread 必须在该位点的**左右两侧都有非-gap碱基**（不含pos本身）
         才能认为它跨越并覆盖了该 gap 位点。如果任一侧全是 '-'，则不覆盖。
        """
        if subread[pos] != '-':
            return True

        left_non_gap = any(
            s != "-" for s in subread[:pos]) if pos > 0 else False
        right_non_gap = any(
            s != "-" for s in subread[pos + 1:]) if pos < len(subread) - 1 else False

        return left_non_gap and right_non_gap

    for pos in range(len(smc)):
        if smc[pos] != "-":
            continue

        # Feature 1: depth = number of subreads covering this position
        covering = [s for s in subreads if _subread_covers(s, pos)]
        n_depth = len(covering)

        # Feature 2: majority base ratio from covering subreads
        non_gap = [s[pos] for s in covering if s[pos] in _VALID_BASES]
        n_evidence = len(non_gap)
        if n_evidence == 0:
            continue

        # Feature 3: SMC homopolymer length + base at this gap pos
        hp_len, hp_base_at_pos = _compute_smc_hp(smc, pos)

        # Filter out sites with no homopolymer context on either side
        # (hp_base == "" → 孤立 gap，不进入训练/评估样本)
        if not hp_base_at_pos:
            continue

        # Label from reference at gap position
        ref_base = ref[pos]
        # 0=gap (correct del), 1=emit (missed call)
        label = 0 if ref_base == "-" else 1

        ratio = n_evidence / n_depth if n_depth > 0 else 0.0
        if ratio < 0.1:
            continue
        
        ratios.append(n_evidence / n_depth if n_depth > 0 else 0.0)
        hp_lens.append(hp_len)
        labels.append(label)
        hp_bases.append(hp_base_at_pos)

    return {
        "ratios": ratios,
        "hp_lens": hp_lens,
        "labels": labels,
        "hp_bases": hp_bases,
    }


def _merge_features(results: list[dict]) -> dict:
    """把多条记录的特征 dict 按记录顺序拼接成全局特征 dict。"""
    merged = {"ratios": [], "hp_lens": [], "labels": [], "hp_bases": []}
    for r in results:
        for k in merged:
            merged[k].extend(r[k])
    return merged


def _extract_features(records: list[dict], workers: int | None = None) -> dict:
    """遍历所有记录的 msa_seqs，提取 SMC gap 位点的特征和 label。

    workers == 1 → 串行；workers > 1 → 多进程并行（按记录分块）。
    并行/串行产出的特征顺序完全一致。

    Returns:
        dict with keys:
            - ratios: list[float]     majority base ratio at each SMC gap pos (evidence / depth)
            - hp_lens: list[int]      smc_hp_length at each SMC gap pos
            - labels: list[int]       1=emit (ref has base), 0=gap (ref is '-')
            - hp_bases: list[str]        SMC homopolymer base at each gap pos
    """
    n = len(records)

    def _serial() -> dict:
        results = [_extract_features_one(rec)
                   for rec in tqdm(records, desc="reading ...")]
        return _merge_features(results)

    if workers == 1 or n <= 1:
        return _serial()

    chunksize = max(1, min(50, (n // workers) + 1))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(
            tqdm(executor.map(_extract_features_one, records, chunksize=chunksize),
                 total=n, desc="reading ..."))
    return _merge_features(results)


# ---------------------------------------------------------------------------
# Threshold search (F1-max)
# ---------------------------------------------------------------------------

def _precision_recall(tp, fp, fn):
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    return prec, rec


# minimum samples per (base, hp) cell for meaningful threshold
_MIN_GROUP = 10
# minimum positive (emit) samples in cell to avoid tiny-sample F1=1.0
_MIN_POS = 5


_WS_RATIOS: list = []
_WS_LABELS: list = []


def _init_threshold_worker(ratios: list, labels: list) -> None:
    """ProcessPoolExecutor initializer：把只读特征数组共享给所有 worker。

    避免每个任务都 pickle 一次大数组；任务只需传 (base, h_thr, sub_indices)。
    """
    global _WS_RATIOS, _WS_LABELS
    _WS_RATIOS = ratios
    _WS_LABELS = labels


def _search_cell_threshold(task: tuple) -> tuple | None:
    """在单个 (base, hp_len) 格子内搜索最优 ratio_thr（F1-max）。

    模块级 worker：ProcessPoolExecutor 需要可 pickle。
    task = (base, h_thr, sub_indices)；ratios/labels 取自进程初始化的 _WS_RATIOS/_WS_LABELS。
    Returns:
        (base, h_thr, {"ratio_threshold", "f1", "n_positive", "n_negative"})
        或 None（该格子无有效阈值）。
    """
    base, h_thr, sub_indices = task
    ratios = _WS_RATIOS
    labels = _WS_LABELS

    cell_n = len(sub_indices)
    n_pos = sum(1 for i in sub_indices if labels[i] == 1)
    n_neg = cell_n - n_pos
    if n_pos == 0 or n_neg == 0 or n_pos < _MIN_POS:
        return None

    # Only search over ratio_thr (hp_thr is fixed as current loop variable)
    sub_ratios = sorted(set(ratios[i] for i in sub_indices))

    best_local_f1 = -1.0
    best_local_r_thr = None

    for r_thr in sub_ratios:
        tp = sum(
            1 for i in sub_indices if ratios[i] > r_thr and labels[i] == 1)
        fp = sum(
            1 for i in sub_indices if ratios[i] > r_thr and labels[i] == 0)
        fn = n_pos - tp

        prec, rec = _precision_recall(tp, fp, fn)
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0

        if f1 > best_local_f1:
            best_local_f1 = f1
            best_local_r_thr = r_thr

    if best_local_r_thr is None:
        return None

    return (base, h_thr, {
        "ratio_threshold": round(best_local_r_thr, 4),
        "f1": round(best_local_f1, 4),
        "n_positive": n_pos,
        "n_negative": n_neg,
    })


def _find_thresholds(features: dict, workers: int | None = None) -> dict | None:
    """对每个 **碱基 × hp_thr** 组合，搜索最优 ratio_thr（F1-max）。

    workers == 1 → 串行；workers > 1 → 多进程并行（按 (base, hp_len) 格子分块）。
    并行/串行产出的 lookup_table / summary_by_base 完全一致（聚合按任务顺序进行）。

    Returns:
        dict with keys:
            - lookup_table: {(base, hp_thr): {"ratio_threshold": float, "f1": float,
                                               "n_positive": int, "n_negative": int}}
                每个 (碱基, hp_len) 的最优 ratio_thr + 该格子正负样本数
            - summary_by_base: {base: {"n_unique_hp": int, "best_f1": float, "hp_values": [int, ...]}}
    """
    ratios = features["ratios"]
    hp_lens = features["hp_lens"]
    labels = features["labels"]
    hp_bases = features["hp_bases"]
    n = len(labels)

    if n == 0:
        return None

    n_emit = sum(labels)
    n_gap = n - n_emit

    # Edge case: all samples same class — no meaningful threshold
    if n_emit == 0 or n_gap == 0:
        classes = "gap" if n_emit == 0 else "emit"
        return {
            "lookup_table": {},
            "summary_by_base": {},
            "_warning": f"所有标签均为 {classes}，无法搜索阈值",
        }

    # --- Group by hp base → list of global indices ---
    groups: dict[str, list[int]] = {}  # base -> [global indices]
    for i in range(n):
        b = hp_bases[i]
        if b not in groups:
            groups[b] = []
        groups[b].append(i)

    # --- 任务构建（串行、廉价）：base 级过滤 + 按 hp_len 分格子 ---
    tasks: list[tuple] = []        # (base, h_thr, sub_indices)
    summary_bases: list[str] = []  # 通过 base 级过滤、需写 summary 的碱基
    for base in sorted(groups.keys()):
        indices_list = groups[base]  # global indices belonging to this base
        n_grp = len(indices_list)
        if n_grp < _MIN_GROUP:
            continue

        grp_emit = sum(1 for i in indices_list if labels[i] == 1)
        grp_gap = n_grp - grp_emit
        if grp_emit == 0 or grp_gap == 0:
            continue

        summary_bases.append(base)

        # Sub-group by hp_len within this base — O(n_grp) lookup dict instead of list scan
        cells: dict[int, list[int]] = {}  # hp_len -> [global indices]
        for i in indices_list:
            h = int(hp_lens[i])
            if h not in cells:
                cells[h] = []
            cells[h].append(i)

        for h_thr in sorted(cells.keys()):
            sub_indices = cells[h_thr]
            if len(sub_indices) < _MIN_GROUP:
                continue
            tasks.append((base, h_thr, sub_indices))

    # --- 并行/串行搜索每个格子的最优 ratio_thr ---
    if workers == 1 or len(tasks) <= 1:
        _init_threshold_worker(ratios, labels)  # 串行路径同样需要共享数组
        results = [_search_cell_threshold(t)
                   for t in tqdm(tasks, desc="searching thresholds ...")]
    else:
        chunksize = max(1, min(50, (len(tasks) // workers) + 1))
        with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_init_threshold_worker,
                initargs=(ratios, labels)) as executor:
            results = list(
                tqdm(executor.map(_search_cell_threshold, tasks, chunksize=chunksize),
                     total=len(tasks), desc="searching thresholds ..."))

    # --- 聚合（按任务顺序，与串行结果逐项一致）---
    lookup_table: dict[tuple[str, int], dict] = {}
    best_f1_by_base: dict[str, float] = {}
    hp_values_by_base: dict[str, list[int]] = {}
    for r in results:
        if r is None:
            continue
        base, h_thr, entry = r
        lookup_table[(base, h_thr)] = entry
        if entry["f1"] > best_f1_by_base.get(base, -1.0):
            best_f1_by_base[base] = entry["f1"]
        hp_values_by_base.setdefault(base, []).append(h_thr)

    summary_by_base: dict[str, dict] = {}
    for base in summary_bases:
        summary_by_base[base] = {
            "n_unique_hp": len(hp_values_by_base.get(base, [])),
            "best_f1": round(best_f1_by_base.get(base, -1.0), 4),
            "hp_values": sorted(set(hp_values_by_base.get(base, []))),
        }

    return {
        "lookup_table": lookup_table,
        "summary_by_base": summary_by_base,
    }


# ---------------------------------------------------------------------------
# Evaluation (baseline vs ratio lookup model)
# ---------------------------------------------------------------------------

def _apply_lookup(features: dict, lookup: dict) -> tuple[list, list]:
    """对每个样本按 (hp_base, hp_len) 查表预测。

    Returns:
        (preds, covered): 与 labels 等长的两个 list。
            preds[i]  = 1 if ratio > thr else 0（有阈值格子）
                        None（该格子被跳过，无阈值）
            covered[i] = True 表示该样本有可用阈值、参与评估
    """
    ratios = features["ratios"]
    hp_lens = features["hp_lens"]
    hp_bases = features["hp_bases"]
    n = len(features["labels"])

    preds: list = []
    covered: list = []
    for i in range(n):
        entry = lookup.get((hp_bases[i], int(hp_lens[i])))
        if entry is None:
            preds.append(None)
            covered.append(False)
        else:
            preds.append(1 if ratios[i] > entry["ratio_threshold"] else 0)
            covered.append(True)
    return preds, covered


def _evaluate(preds: list, labels: list) -> dict:
    """给定等长的预测与真实标签，计算混淆矩阵与衍生指标（正类 = emit=1）。"""
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


def _evaluate_by_cell(features: dict, preds: list, covered: list) -> list:
    """按 (base, hp_len) 格子分别评估模型 vs baseline。

    covered 样本按 (hp_base, hp_len) 分组（= lookup 的键）。
    Returns:
        list[tuple[(str, int), int, dict, dict]]  (key, n_cell, model_eval, base_eval)
    """
    hp_bases = features["hp_bases"]
    hp_lens = features["hp_lens"]
    labels = features["labels"]

    cell_groups: dict[tuple[str, int], list[int]] = {}
    for i in range(len(labels)):
        if not covered[i]:
            continue
        key = (hp_bases[i], int(hp_lens[i]))
        cell_groups.setdefault(key, []).append(i)

    rows = []
    for key in sorted(cell_groups):
        idx = cell_groups[key]
        cell_labels = [labels[i] for i in idx]
        cell_preds = [preds[i] for i in idx]
        model_eval = _evaluate(cell_preds, cell_labels)
        base_eval = _evaluate([0] * len(idx), cell_labels)
        rows.append((key, len(idx), model_eval, base_eval))
    return rows


# ---------------------------------------------------------------------------
# Visualization (scatter plot)
# ---------------------------------------------------------------------------

def _plot_scatter(features: dict, thresholds: dict | None, output_path: str):
    """绘制重插入决策散点图（按 majority base 着色 + emit/gap 叠加）。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 未安装，跳过可视化输出。", file=sys.stderr)
        return

    ratios = features["ratios"]
    hp_lens = features["hp_lens"]
    labels = features["labels"]
    maj_bases = features["hp_bases"]
    lookup: dict[tuple[str, int], dict] = (
        thresholds.get("lookup_table") if thresholds else {}) or {}

    # Use log scale on Y to handle extreme hp_len values; floor at 1
    max_hp = max(hp_lens) if hp_lens else 1
    y_max_log = min(max(5, max_hp), 100)

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_yscale("log", base=2)
    ax.set_ylim(bottom=0.5, top=y_max_log)
    ax.set_yticks([1, 2, 4, 8, 16, 32, 64])
    y_labels = [str(p) for p in ax.get_yticks()]

    total = len(labels)
    s_size = max(3, min(8, 200 / max(total, 1) ** 0.4))

    # Color palette by base, overlaid with emit/gap semantics
    _BASE_COLORS = {"G": "#E63946", "A": "#219EBC",
                    "C": "#2EC4B6", "T": "#F77F00"}
    # Global alpha for all sub-base layers (emit points slightly brighter)
    layer_alpha = min(0.5, 0.3 / max(total ** 0.3, 1))

    base_layers: dict[str, list[int]] = {}
    for i in range(len(labels)):
        b = maj_bases[i]
        if b not in base_layers:
            base_layers[b] = []
        base_layers[b].append(i)

    # Plot per-base layers (only those with enough points to matter)
    _MIN_LAYER = 5
    for base in sorted(base_layers.keys()):
        indices = base_layers[base]
        if len(indices) < _MIN_LAYER:
            continue
        x = [ratios[i] for i in indices]
        y = [hp_lens[i] for i in indices]
        # Split into emit/gap sub-layers within this base
        for label, color in [(1, _BASE_COLORS.get(base, "#8D99AE") + "FF"), (0, _BASE_COLORS.get(base, "#8D99AE") + "80")]:
            sub_x = [x[j]
                     for j in range(len(indices)) if labels[indices[j]] == label]
            sub_y = [y[j]
                     for j in range(len(indices)) if labels[indices[j]] == label]
            cls_name = "emit" if label else "gap"
            ax.scatter(sub_x, sub_y, c=color, s=s_size, alpha=layer_alpha * (1.2 if label else 0.8),
                       label=f"{base}-{cls_name} n={len(sub_x)}")

    # Threshold lines from best lookup entry
    if lookup:
        best_key = max(lookup.keys(), key=lambda k: lookup[k]["f1"])
        b, h_thr = best_key
        r_thr = lookup[best_key]["ratio_threshold"]
        ax.axvline(x=r_thr, color="#457B9D", linewidth=1.2, linestyle="--",
                   label=f"({b}, hp={h_thr}): ratio > {r_thr:.3f}")

        m_f1 = lookup[best_key]["f1"]
        metric_str = (f"F1={m_f1:.4f} | 查表阈值: base={b}, hp={h_thr}, "
                      f"ratio > {r_thr:.3f}")
    else:
        ax.text(0.5, 0.95, "No threshold data (uniform labels)", transform=ax.transAxes,
                ha="center", va="top", fontsize=12, color="gray")
        metric_str = ""

    ax.set_title(
        f"Re-insertion Decision Scatter Plot\n{metric_str}", fontsize=13)
    ax.set_xlabel("majority_base_ratio (evidence strength)", fontsize=11)
    ax.set_ylabel("smc_hp_length (log scale)", fontsize=11)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    ax.grid(True, alpha=0.3, which="both")
    ax.yaxis.set_ticklabels(y_labels)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Re-insertion 决策模型 — SMC gap 位点预测 emit base or keep gap",
    )
    parser.add_argument("input", help="ASRTC JSONL 文件路径")
    parser.add_argument("--max-records", type=int, default=None,
                        help="最大处理记录数 (default: 全部)")
    parser.add_argument("--workers", type=int, default=None,
                        help="并行提取特征的工作进程数 (default: min(CPU核数, 8); 1=串行)")
    parser.add_argument("--find-workers", type=int, default=None,
                        help="并行搜索阈值的工作进程数 (default: 沿用 --workers; 1=串行)")
    parser.add_argument("--output-dir", default=".",
                        help="输出目录 (PNG + CSV) (default: 当前目录)")
    args = parser.parse_args(argv)

    records = _load_records(args.input, args.max_records)
    if not records:
        print("[ERROR] 未加载到任何记录。", file=sys.stderr)
        return 1

    # Phase 1: Extract features (multicore via --workers)
    workers = args.workers if args.workers else min(os.cpu_count() or 1, 8)
    print(f"正在提取特征 (workers={workers})...")
    features = _extract_features(records, workers=workers)
    n = len(features["labels"])
    print(f"  提取 {n} 个 SMC gap 位点特征")

    if n == 0:
        print("[WARN] 无有效 SMC gap 位点，退出。", file=sys.stderr)
        return 0

    # Phase 2: Find best thresholds (multicore via --find-workers)
    print("正在搜索最佳阈值...")
    find_workers = args.find_workers if args.find_workers else workers
    thresholds = _find_thresholds(features, workers=find_workers)

    warning = None
    if thresholds and "_warning" in thresholds:
        warning = thresholds.pop("_warning")

    if warning:
        print(f"\n[WARN] {warning}")
        plot_path = f"{args.output_dir}/reinsert_scatter.png"
        _plot_scatter(features, None, plot_path)
        print(f"散点图已保存至: {plot_path}")
        return 0

    if thresholds and thresholds.get("lookup_table"):
        lookup = thresholds["lookup_table"]

        print(f"\n{_SEP}")
        print("--- 按 (碱基, hp_len) 查表阈值 (F1-max ratio_thr) ---")
        print(
            f"{'Base':<6} {'hp=':>4s} {'ratio_thr':>10s} {'F1':>6s} {'+pos':>6s} {'-neg':>6s}")
        for (base, h_thr) in sorted(lookup.keys()):
            t = lookup[(base, h_thr)]
            print(f"{base:<6} {h_thr:>4d} {t['ratio_threshold']:>10.4f} {t['f1']:>6.4f} "
                  f"{t['n_positive']:>6d} {t['n_negative']:>6d}")

        # Best lookup entry
        best_key = max(
            lookup.keys(), key=lambda k: lookup[k]["f1"]) if lookup else None
        if best_key:
            b, h = best_key
            bf = lookup[best_key]["f1"]
            rt = lookup[best_key]["ratio_threshold"]
            print(f"\n最佳查表项: ({b}, hp={h}) → ratio > {rt:.4f}  (F1={bf:.4f})")

        # --- Evaluation: baseline (current SMC = all-gap) vs ratio lookup model ---
        preds, covered = _apply_lookup(features, lookup)
        covered_idx = [i for i, c in enumerate(covered) if c]
        n_cov = len(covered_idx)
        total = len(features["labels"])

        print(f"\n{_SEP}")
        print(f"--- 评估指标 (仅 covered 样本) ---")
        if n_cov == 0:
            print(f"coverage = 0/{total} (0.0%)，无有阈值格子，跳过评估。")
        else:
            model_preds = [preds[i] for i in covered_idx]
            cov_labels = [features["labels"][i] for i in covered_idx]
            model_eval = _evaluate(model_preds, cov_labels)
            base_eval = _evaluate([0] * n_cov, cov_labels)

            print(f"{'方法':<20} {'TP':>4s} {'FP':>4s} {'FN':>4s} {'TN':>4s} "
                  f"{'Precision':>9s} {'Recall':>7s} {'F1':>6s} {'Accuracy':>8s}")
            print(f"{'baseline (SMC全gap)':<20} {base_eval['tp']:>4d} {base_eval['fp']:>4d} "
                  f"{base_eval['fn']:>4d} {base_eval['tn']:>4d} "
                  f"{base_eval['precision']:>9.4f} {base_eval['recall']:>7.4f} "
                  f"{base_eval['f1']:>6.4f} {base_eval['accuracy']:>8.4f}")
            print(f"{'ratio 查表模型':<20} {model_eval['tp']:>4d} {model_eval['fp']:>4d} "
                  f"{model_eval['fn']:>4d} {model_eval['tn']:>4d} "
                  f"{model_eval['precision']:>9.4f} {model_eval['recall']:>7.4f} "
                  f"{model_eval['f1']:>6.4f} {model_eval['accuracy']:>8.4f}")

            cov_pct = 100.0 * n_cov / total if total > 0 else 0.0
            print(f"coverage = {n_cov}/{total} ({cov_pct:.1f}%)  "
                  f"[有阈值格子占全部样本的比例]")

            # --- 分格评估 (base × hp_len): 模型 vs baseline ---
            print(f"\n{_SEP}")
            print("--- 分格评估 (base × hp_len): 模型 vs baseline ---")
            print(f"{'Base':<6} {'hp=':>4s} {'N':>7s} {'方法':<20} "
                  f"{'TP':>4s} {'FP':>4s} {'FN':>4s} {'TN':>4s} "
                  f"{'Precision':>9s} {'Recall':>7s} {'F1':>6s} {'Accuracy':>8s}")
            for (base, h_thr), n_cell, m, b in _evaluate_by_cell(
                    features, preds, covered):
                for name, ev in (("baseline (SMC全gap)", b),
                                 ("ratio 查表模型", m)):
                    print(f"{base:<6} {h_thr:>4d} {n_cell:>7d} {name:<20} "
                          f"{ev['tp']:>4d} {ev['fp']:>4d} {ev['fn']:>4d} {ev['tn']:>4d} "
                          f"{ev['precision']:>9.4f} {ev['recall']:>7.4f} "
                          f"{ev['f1']:>6.4f} {ev['accuracy']:>8.4f}")
    else:
        print("[WARN] 无法找到有效阈值。")

    # Phase 3: Visualize (threshold lines from best single lookup entry)
    if HAS_MPL:
        plot_path = f"{args.output_dir}/reinsert_scatter.png"
        _plot_scatter(features, thresholds, plot_path)
        print(f"\n散点图已保存至: {plot_path}")
    else:
        print("\n[WARN] matplotlib 未安装，跳过可视化。")

    return 0


if __name__ == "__main__":
    raise SystemExit(main_cli())
