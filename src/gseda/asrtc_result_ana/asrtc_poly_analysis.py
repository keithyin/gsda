#!/usr/bin/env python3
"""ASRTC MSA gap ratio + re-insertion 分析 — 纯 MSA（SMC + subreads），无需 reference。

两个分析维度：

1. SMC-call base 位点的 gap_ratio
    SMC 已 call ACGT，subreads 中 gap 比例是多少？
    → 用于判断 SMC call 的可信度（gap_ratio 越高越可疑）

2. SMC-gap 位点的 non_gap_ratio（re-insertion 证据强度）
    SMC 是 "-"（del），subreads 中有多少非-gap base？
    → 用于 re-insertion：non_gap_ratio >= T 时，插入 base 的可信度
    同时统计 subreads majority base 的一致性。

全程不依赖 reference，适用于推理阶段。

数据格式（msa_seqs 列表索引）:
    [0]  quality string      (digits 0-9 + #，跳过)
    [1]  SMC consensus       DNA bases ACGT-
    [2]  reference           DNA bases ACGT-  (present in JSONL, unused)
    [3:] subreads            DNA bases ACGT-

用法:
    python asrtc_poly_analysis.py input.asrtc.txt
    python asrtc_poly_analysis.py input.asrtc.txt --max-records 10000
    python asrtc_poly_analysis.py input.asrtc.txt --output json -o gap_ratio.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_BASES = "ACGT"       # valid DNA bases
_VALID_SMC = frozenset(_BASES + "-")  # allowed chars in SMC consensus


# ---------------------------------------------------------------------------
# Data loading & helpers
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


def _collect_gap_stats(records: list[dict]) -> dict:
    """遍历所有记录的 msa_seqs，收集 gap_ratio + re-insertion 统计。

    不依赖 reference — 适用于推理阶段。

    Returns:
        dict with keys:
            - all_ratios: list[float]  SMC-call base 位点的 gap_ratio
            - by_base: {base: [gap_ratios]}  SMC-call base 按碱基分类
            - consensus_agree: list[tuple[gap_ratio, bool]]
              SMC-call 位点上 subreads majority 与 SMC 是否一致
            - n_smc_calls: int SMC 已 call base 的位点数
            - reinsert: dict  SMC-gap 位点的 re-insertion 分析
                - n_positions: int
                - non_gap_ratios: list[float]  subreads 中证据强度 (non_gaps / n_sub)
                - by_majority_base: {base: [evidence_strengths]}
                - agree_by_evidence_bin: list[tuple[bin_label, total, agreed]]
            - n_records: int number of records processed
    """
    all_ratios = []
    by_base: dict[str, list[float]] = {b: [] for b in _BASES}
    consensus_agree: list[tuple[float, bool]] = []

    n_smc_calls = 0
    n_records = len(records)

    # re-insertion accumulators
    re_n_pos = 0
    non_gap_ratios: list[float] = []
    re_by_majority: dict[str, list[float]] = {b: [] for b in _BASES}
    re_agree_bins: dict[str, list[int]] = {}

    for rec in records:
        msa = rec["msa_seqs"]
        smc = msa[1]
        _ref = msa[2]  # present in JSONL but unused for analysis
        subreads = msa[3:]
        n_sub = len(subreads)
        if n_sub == 0:
            continue

        for pos in range(len(smc)):
            base = smc[pos]

            if base not in _VALID_SMC:
                # Skip quality chars, other non-DNA symbols in consensus
                continue

            if base != "-":
                # === SMC-call base 位点 (existing logic, kept as-is) ===
                n_smc_calls += 1
                n_gaps = sum(1 for s in subreads if s[pos] == "-")
                gap_ratio = n_gaps / n_sub

                all_ratios.append(gap_ratio)
                by_base[base].append(gap_ratio)

                non_gap = [s[pos] for s in subreads if s[pos] != "-"]
                if non_gap:
                    majority_base, majority_count = Counter(non_gap).most_common(1)[0]
                    agreed = (majority_base == base and majority_count > len(non_gap) / 2)
                    consensus_agree.append((gap_ratio, agreed))

            else:
                # === SMC-gap 位点 (new re-insertion analysis) ===
                re_n_pos += 1
                non_gap = [s[pos] for s in subreads if s[pos] != "-"]
                n_evidence = len(non_gap)
                if n_evidence == 0:
                    # all subreads also gap → no evidence to insert
                    continue

                evidence_str = n_evidence / n_sub
                non_gap_ratios.append(evidence_str)

                # Majority base among non-gap subreads at this position
                maj_counts = Counter(non_gap).most_common(1)[0]
                maj_base = maj_counts[0]
                if maj_base in _BASES:
                    re_by_majority[maj_base].append(evidence_str)

                # Agreement with subread majority at this evidence bin
                bucket = min(int(evidence_str * 10), 9)
                lo, hi = f"{bucket * 10}%", f"{(bucket + 1) * 10}%"
                key = f"[{lo}, {hi})"
                if key not in re_agree_bins:
                    re_agree_bins[key] = [0, 0]
                re_agree_bins[key][0] += 1
                # "agreed" for reinsertion means subread majority > half of non-gap
                # (i.e., there IS a clear majority)
                if maj_counts[1] > n_evidence / 2:
                    re_agree_bins[key][1] += 1

    return {
        "all_ratios": all_ratios,
        "by_base": by_base,
        "consensus_agree": consensus_agree,
        "n_smc_calls": n_smc_calls,
        "reinsert": {
            "n_positions": re_n_pos,
            "non_gap_ratios": non_gap_ratios,
            "by_majority_base": re_by_majority,
            "agree_by_evidence_bin": dict(sorted(re_agree_bins.items())),
        },
        "n_records": n_records,
    }


# ---------------------------------------------------------------------------
# Analysis & formatting
# ---------------------------------------------------------------------------

def _percentiles(ratios: list[float], percs: list[int] | None = None) -> dict[float, float]:
    """返回 gap_ratio 的指定分位数。"""
    if not ratios:
        return {}
    sorted_r = sorted(ratios)
    n = len(sorted_r)
    if percs is None:
        percs = [10, 25, 50, 75, 90, 95, 99]
    result: dict[float, float] = {}
    for p in percs:
        idx = min(int(p / 100 * n), n - 1)
        result[p] = sorted_r[idx]
    return result


def _histogram(ratios: list[float], nbins: int = 20) -> list[tuple[str, int, float]]:
    """按 gap_ratio 区间画直方图，返回 [(bin_label, count, pct), ...]。"""
    if not ratios:
        return []
    bins = [i / nbins for i in range(nbins + 1)]
    hist: list[int] = [0] * nbins
    for r in ratios:
        bucket = min(int(r * nbins), nbins - 1)
        hist[bucket] += 1
    n = len(ratios)
    result: list[tuple[str, int, float]] = []
    for i in range(nbins):
        lo, hi = bins[i], bins[i + 1]
        label = f"[{lo:.2f}, {hi:.2f})" if i < nbins - 1 else f"[{lo:.2f}, 1.00]"
        result.append((label, hist[i], hist[i] / n * 100))
    return result


def _format_text(stats: dict) -> str:
    """生成可读文本报告。"""
    all_r = stats["all_ratios"]
    lines: list[str] = []
    sep = "=" * 70

    lines.append(sep)
    lines.append("ASRTC MSA gap_ratio + re-insertion 分析")
    lines.append(sep)
    lines.append(f"记录数:          {stats['n_records']:,}")
    lines.append(f"SMC-call位点:    {stats['n_smc_calls']:,}")
    ri = stats["reinsert"]
    n_revid = len(ri["non_gap_ratios"])
    lines.append(f"SMC-gap位点(有证据):  {n_revid:,}/{ri['n_positions']:,}")
    if all_r:
        gapped_n = sum(1 for r in all_r if r > 0)
        lines.append(f"有gap的位点:   {gapped_n:,} ({gapped_n / len(all_r) * 100:.1f}%)")
    lines.append("")

    # --- SMC-call base gap_ratio distribution ---
    if not all_r:
        lines.append("[WARN] 无 SMC 已 call base 的位点。")
    else:
        percs = _percentiles(all_r, [1, 5, 10, 25, 50, 75, 90, 95, 99])
        mean_g = sum(all_r) / len(all_r)

        lines.append("--- SMC-call base gap ratio 分布 ---")
        lines.append(f"  Mean: {mean_g:.4f}")
        lines.append(f"  P1/P5/P10/P25/P50/P75/P90/P95/P99:")
        vals_str = "/".join(f"P{p}={percs[p]:.3f}" for p in [1, 5, 10, 25, 50, 75, 90, 95, 99])
        lines.append(f"    {vals_str}")

        # Histogram
        hist = _histogram(all_r)
        max_count_bin = max((c for _, c, _ in hist), default=1) if hist else 1
        lines.append("")
        lines.append("  直方图 (gap_ratio → SMC-call位点占比):")
        for label, count, pct in hist:
            bar_len = int(count / max(max_count_bin, 1) * 40) if all_r and count > 0 else 0
            bar = "#" * bar_len if count > 0 else ""
            lines.append(f"  {label:>13s}: {count:>8,} ({pct:>5.1f}%) |{bar}")

        # Per-base distribution
        lines.append("")
        lines.append("--- 按 SMC call base 分类 ---")
        for base in _BASES:
            vals = stats["by_base"][base]
            if not vals:
                continue
            bpercs = _percentiles(vals, [25, 50, 75, 90, 95])
            bmean = sum(vals) / len(vals) if vals else 0
            lines.append(f"\n  SMC={base} (n={len(vals):,}):")
            for p in [25, 50, 75, 90, 95]:
                lines.append(f"    P{p}: {bpercs[p]:.3f}")
            lines.append(f"    Mean: {bmean:.4f}")

        # Subread consensus agreement by gap_ratio bin
        lines.append("")
        lines.append("--- subreads majority 与 SMC 一致率 (按 gap ratio bin) ---")
        agree_bins: dict[str, list[int]] = {}
        for gap_r, agreed in stats["consensus_agree"]:
            bucket = min(int(gap_r * 10), 9)
            lo, hi = f"{bucket * 10:2d}%", f"{(bucket + 1) * 10:2d}%"
            key = f"[{lo}, {hi})"
            if key not in agree_bins:
                agree_bins[key] = [0, 0]
            agree_bins[key][0] += 1
            if agreed:
                agree_bins[key][1] += 1
        for key in sorted(agree_bins.keys()):
            total, matched = agree_bins[key]
            rate = matched / total * 100 if total else 0
            lines.append(f"  {key:>13s}: {matched:>6,}/{total:>6,} ({rate:5.1f}%)")

        # Recommendation thresholds for SMC-call
        lines.append("")
        lines.append(sep)
        t90 = percs.get(90, 0)
        t95 = percs.get(95, 0)
        lines.append(f"  gap_ratio <= {t90:.3f} (P90): "
                      f"{sum(1 for r in all_r if r <= t90):>6,}/{len(all_r):>6,} 位点可用")
        lines.append(f"  gap_ratio <= {t95:.3f} (P95): "
                      f"{sum(1 for r in all_r if r <= t95):>6,}/{len(all_r):>6,} 位点可用")

    # --- Re-insertion analysis ---
    lines.append("")
    lines.append(sep)
    lines.append("--- re-insertion 阈值参考 (SMC gap 位点) ---")
    ri_non_gap = ri["non_gap_ratios"]
    if not ri_non_gap:
        lines.append("  无 SMC-gap 位点有 subreads 证据，无法分析。")
    else:
        percs_ri = _percentiles(ri_non_gap, [1, 5, 10, 25, 50, 75, 90, 95, 99])
        mean_ri = sum(ri_non_gap) / len(ri_non_gap)
        lines.append(f"\n  non_gap_ratio (证据强度) 分布:")
        lines.append(f"    Mean: {mean_ri:.4f}")
        lines.append(f"    P1/P5/P10/P25/P50/P75/P90/P95/P99:")
        vals_str_ri = "/".join(f"P{p}={percs_ri[p]:.3f}" for p in [1, 5, 10, 25, 50, 75, 90, 95, 99])
        lines.append(f"      {vals_str_ri}")

        # Histogram of evidence strength
        hist_ri = _histogram(ri_non_gap)
        max_count_ri = max((c for _, c, _ in hist_ri), default=1) if hist_ri else 1
        lines.append("")
        lines.append("  直方图 (证据强度 → 位点占比):")
        for label, count, pct in hist_ri:
            bar_len = int(count / max(max_count_ri, 1) * 40) if ri_non_gap and count > 0 else 0
            bar = "#" * bar_len if count > 0 else ""
            lines.append(f"  {label:>13s}: {count:>8,} ({pct:>5.1f}%) |{bar}")

        # Per majority base
        lines.append("")
        lines.append("  按 subreads majority base 分类:")
        for base in _BASES:
            vals = ri["by_majority_base"][base]
            if not vals:
                continue
            bpercs = _percentiles(vals, [25, 50, 75, 90, 95])
            bmean = sum(vals) / len(vals) if vals else 0
            lines.append(f"\n    majority={base} (n={len(vals):,}):")
            for p in [25, 50, 75, 90, 95]:
                lines.append(f"      P{p}: {bpercs[p]:.3f}")
            lines.append(f"      Mean: {bmean:.4f}")

        # Agreement by evidence bin
        lines.append("")
        lines.append("  subreads majority >一半一致率 (按证据区间):")
        for key in sorted(ri["agree_by_evidence_bin"].keys()):
            total, matched = ri["agree_by_evidence_bin"][key]
            rate = matched / total * 100 if total else 0
            lines.append(f"  {key:>13s}: {matched:>6,}/{total:>6,} ({rate:5.1f}%)")

        # Recommended thresholds
        lines.append("")
        t25 = percs_ri.get(25, 0)
        t50 = percs_ri.get(50, 0)
        t75 = percs_ri.get(75, 0)
        lines.append(f"  non_gap_ratio >= {t25:.3f} (P25): "
                      f"{sum(1 for r in ri_non_gap if r >= t25):>6,}/{len(ri_non_gap):>6,} 位点可用")
        lines.append(f"  non_gap_ratio >= {t50:.3f} (P50): "
                      f"{sum(1 for r in ri_non_gap if r >= t50):>6,}/{len(ri_non_gap):>6,} 位点可用")
        lines.append(f"  non_gap_ratio >= {t75:.3f} (P75): "
                      f"{sum(1 for r in ri_non_gap if r >= t75):>6,}/{len(ri_non_gap):>6,} 位点可用")

    lines.append(sep)

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def _format_json(stats: dict) -> str:
    """生成 JSON 格式报告。"""
    all_r = stats["all_ratios"]
    ri = stats["reinsert"]
    non_gap_ratios = ri["non_gap_ratios"]

    report: dict = {
        "summary": {
            "n_records": stats["n_records"],
            "n_smc_calls": stats["n_smc_calls"],
            "n_smc_gaps_with_evidence": len(non_gap_ratios),
        },
        "smc_call_gap_ratio": {},
        "smc_call_by_base": {},
        "consensus_agreement_by_gap_bin": [],
        "reinsert_analysis": {},
    }

    if all_r:
        # gap ratio percentiles
        pcts = _percentiles(all_r, [1, 5, 10, 25, 50, 75, 90, 95, 99])
        report["smc_call_gap_ratio"] = {
            "mean": round(sum(all_r) / len(all_r), 4),
            "percentiles": {str(k): round(v, 4) for k, v in pcts.items()},
            "n_gapped": sum(1 for r in all_r if r > 0),
        }

        # Per-base gap ratio
        for base in _BASES:
            vals = stats["by_base"][base]
            if vals:
                report["smc_call_by_base"][base] = {
                    "count": len(vals),
                    "mean": round(sum(vals) / len(vals), 4),
                    "percentiles": {str(k): round(v, 4) for k, v in _percentiles(vals).items()},
                }

        # Consensus agreement by gap bin
        agree_bins: dict[str, list[int]] = {}
        for gap_r, agreed in stats["consensus_agree"]:
            bucket = min(int(gap_r * 10), 9)
            lo, hi = f"{bucket * 10}%", f"{(bucket + 1) * 10}%"
            if (lo, hi) not in agree_bins:
                agree_bins[(lo, hi)] = [0, 0]
            agree_bins[(lo, hi)][0] += 1
            if agreed:
                agree_bins[(lo, hi)][1] += 1
        for (lo, hi), (total, matched) in sorted(agree_bins.items()):
            report["consensus_agreement_by_gap_bin"].append({
                "gap_ratio_range": f"{lo} - {hi}",
                "total_positions": total,
                "smc_matches_majority": matched,
                "agreement_rate": round(matched / total * 100, 2) if total else 0,
            })

    # Re-insertion analysis (always included, may be empty)
    if non_gap_ratios:
        percs_ri = _percentiles(non_gap_ratios, [1, 5, 10, 25, 50, 75, 90, 95, 99])
        report["reinsert_analysis"] = {
            "n_positions_with_evidence": len(non_gap_ratios),
            "evidence_strength_distribution": {
                "mean": round(sum(non_gap_ratios) / len(non_gap_ratios), 4),
                "percentiles": {str(k): round(v, 4) for k, v in percs_ri.items()},
            },
            "by_majority_base": {},
            "evidence_agreement_by_bin": [],
        }
        for base in _BASES:
            vals = ri["by_majority_base"][base]
            if vals:
                report["reinsert_analysis"]["by_majority_base"][base] = {
                    "count": len(vals),
                    "mean": round(sum(vals) / len(vals), 4),
                    "percentiles": {str(k): round(v, 4) for k, v in _percentiles(vals).items()},
                }

        for key in sorted(ri["agree_by_evidence_bin"].keys()):
            total, matched = ri["agree_by_evidence_bin"][key]
            report["reinsert_analysis"]["evidence_agreement_by_bin"].append({
                "evidence_range": key,
                "total_positions": total,
                "majority_consensus": matched,
                "agreement_rate": round(matched / total * 100, 2) if total else 0,
            })

    return json.dumps(report, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="ASRTC MSA gap_ratio + re-insertion 分析 — 纯 MSA，无需 reference",
    )
    parser.add_argument("input", help="ASRTC JSONL 文件路径")
    parser.add_argument("--max-records", type=int, default=None,
                        help="最大处理记录数 (default: 全部)")
    parser.add_argument("-o", "--output", choices=["text", "json"], default="text",
                        help="输出格式 (default: text)")
    args = parser.parse_args(argv)

    records = _load_records(args.input, args.max_records)
    if not records:
        print("[ERROR] 未加载到任何记录。", file=sys.stderr)
        return 1

    stats = _collect_gap_stats(records)

    if args.output == "json":
        print(_format_json(stats))
    else:
        print(_format_text(stats))

    return 0


if __name__ == "__main__":
    raise SystemExit(main_cli())
