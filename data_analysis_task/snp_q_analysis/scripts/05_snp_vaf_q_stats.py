#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""分析变体位点（SNP/INS/DEL），仅保留 VAF（变异等位基因频率）> 阈值的稳定位点，
并输出这些位点的 Q 值基本统计信息。

VAF = 该位点变体事件数 / n_reads（full-length 比对 read 数）。
同一 read 在某位点至多贡献 1 个事件，故事件数 = 携带该变体的 read 数。

输入为 03_variant_calling.py 输出的变体 TSV（列：barcode read_name ref_pos(1-based)
ref_base alt_base base_q strand read_len variant_type）。

输出：
    {out_prefix}_{TYPE}_sites.tsv    逐位点表：ref_pos ref_base alt_base n_events vaf q_min q_p5 q_p10 q_p25 q_p50 q_p75 q_max q_mean
    {out_prefix}_{TYPE}_overall.tsv  保留位点全部事件的总基本统计（N/min/P1..P99/max/mean/std）
    {out_prefix}_{TYPE}_q_hist.png   保留位点事件 Q 值直方图

用法:
    python3 05_snp_vaf_q_stats.py results/variant/Barcode01.variants.tsv \
        --type DEL --n-reads 94516 --out-prefix results/vaf_q_stats/Barcode01
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)

SITE_Q_COLUMNS = ["ref_pos", "ref_base", "alt_base", "n_events", "vaf",
                  "q_min", "q_p5", "q_p10", "q_p25", "q_p50", "q_p75", "q_max", "q_mean"]
OVERALL_COLUMNS = ["group", "N", "min"] + [f"P{p}" for p in (1, 5, 10, 25, 50, 75, 90, 95, 99)] \
    + ["max", "mean", "std"]


def load_snp_by_pos(paths: list[str], vtype: str):
    """读取变体 TSV，按 ref_pos 聚合等位基因与全部 base_q，仅保留指定类型。

    Returns:
        dict: ref_pos -> {"alleles": Counter[(ref_base, alt_base)], "n": int, "qs": list[int]}
    """
    by_pos: dict[int, dict] = {}
    for path in paths:
        with open(path, "r") as f:
            header = f.readline()
            if header and not header.startswith("barcode"):
                f.seek(0)
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[8] != vtype:
                    continue
                rp = int(parts[2])
                rec = by_pos.setdefault(rp, {"alleles": Counter(), "n": 0, "qs": []})
                rec["alleles"][(parts[3], parts[4])] += 1
                rec["n"] += 1
                rec["qs"].append(int(parts[5]))
    return by_pos


def dominant_allele(rec: dict) -> tuple[str, str]:
    """返回该位点出现次数最多的 (ref_base, alt_base)。"""
    if rec["alleles"]:
        return rec["alleles"].most_common(1)[0][0]
    return ("-", "-")


def site_q_stats(qs: list[int]) -> list[float]:
    """返回 [q_min, q_p5, q_p10, q_p25, q_p50, q_p75, q_max, q_mean]。"""
    arr = np.asarray(qs, dtype=float)
    qs_val = np.percentile(arr, [5, 10, 25, 50, 75])
    return [float(arr.min()), *[float(q) for q in qs_val],
            float(arr.max()), float(arr.mean())]


def overall_stats(qs: list[int]) -> list[float]:
    """返回 [N, min, P1..P99, max, mean, std]，列序同 OVERALL_COLUMNS[1:]。"""
    arr = np.asarray(qs, dtype=float)
    if arr.size == 0:
        return [0.0] * (1 + 2 + 9 + 2)
    pcts = np.percentile(arr, [1, 5, 10, 25, 50, 75, 90, 95, 99])
    return [float(arr.size), float(arr.min()),
            *[float(q) for q in pcts],
            float(arr.max()), float(arr.mean()), float(arr.std(ddof=0))]


def plot_histogram(qs: list[int], title: str, out_path: str) -> None:
    """绘制保留位点事件 Q 值直方图。"""
    arr = np.asarray(qs, dtype=int)
    lo, hi = int(arr.min()), int(arr.max())
    plt.figure(figsize=(12, 6))
    plt.hist(arr, bins=list(range(lo, hi + 2)), color="#4c72b0",
             edgecolor="white", align="left")
    plt.xlabel("Variant base quality (PHRED)")
    plt.ylabel("Variant event count")
    plt.title(title)
    plt.xticks(range(lo, hi + 1, max(1, (hi - lo) // 25)))
    stats_text = f"N = {arr.size:,}\nmean = {arr.mean():.1f}\nmedian = {np.median(arr):.0f}"
    plt.text(0.02, 0.95, stats_text, transform=plt.gca().transAxes,
             va="top", ha="left", fontsize=11,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    log.info("直方图已保存: %s", out_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="变体位点 Q 值统计（仅保留 VAF > 阈值的位点）。",
    )
    parser.add_argument("snp_tsv", nargs="+", help="变体 TSV 路径，支持通配符")
    parser.add_argument("--type", choices=["SNP", "INS", "DEL"], default="SNP",
                        help="变体类型 (默认 SNP)")
    parser.add_argument("--n-reads", type=int, required=True,
                        help="full-length 比对 read 数（VAF 分母）")
    parser.add_argument("--vaf-threshold", type=float, default=0.05,
                        help="VAF 阈值，保留 VAF > 阈值的位点 (默认 0.05)")
    parser.add_argument("--out-prefix", default="results/vaf_q_stats/snp_vaf_q_stats",
                        help="输出前缀")
    args = parser.parse_args()
    vtype = args.type

    paths = []
    for pat in args.snp_tsv:
        paths.extend(glob.glob(pat))
    paths = [p for p in paths if os.path.exists(p)]
    if not paths:
        log.error("未找到变体 TSV 文件")
        sys.exit(1)

    by_pos = load_snp_by_pos(paths, vtype)
    total_sites = len(by_pos)
    log.info("%s 有事件的位点数: %d", vtype, total_sites)

    kept = []
    for rp in sorted(by_pos):
        n = by_pos[rp]["n"]
        vaf = n / args.n_reads
        if vaf > args.vaf_threshold:
            kept.append((rp, by_pos[rp], vaf))
    log.info("%s VAF > %.2f 的位点数: %d", vtype, args.vaf_threshold, len(kept))

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)

    # 逐位点表
    sites_path = f"{args.out_prefix}_{vtype}_sites.tsv"
    with open(sites_path, "w") as f:
        f.write("\t".join(SITE_Q_COLUMNS) + "\n")
        for rp, rec, vaf in kept:
            q_min, q_p5, q_p10, q_p25, q_p50, q_p75, q_max, q_mean = site_q_stats(rec["qs"])
            ref_base, alt_base = dominant_allele(rec)
            f.write("\t".join([
                str(rp), ref_base, alt_base, str(rec["n"]), f"{vaf:.4f}",
                f"{q_min:.0f}", f"{q_p5:.0f}", f"{q_p10:.0f}", f"{q_p25:.0f}",
                f"{q_p50:.0f}", f"{q_p75:.0f}", f"{q_max:.0f}", f"{q_mean:.2f}",
            ]) + "\n")
    log.info("逐位点表已保存: %s", sites_path)

    # 总体统计
    all_qs = [q for _, rec, _ in kept for q in rec["qs"]]
    overall_path = f"{args.out_prefix}_{vtype}_overall.tsv"
    with open(overall_path, "w") as f:
        f.write("\t".join(OVERALL_COLUMNS) + "\n")
        f.write("\t".join(["combined"] + [str(v) for v in overall_stats(all_qs)]) + "\n")
    log.info("总体统计已保存: %s", overall_path)

    # 直方图
    if all_qs:
        plot_histogram(all_qs, f"{vtype} baseQ distribution (VAF > {args.vaf_threshold})",
                       f"{args.out_prefix}_{vtype}_q_hist.png")
    else:
        log.warning("无保留位点，跳过直方图")

    # 控制台打印
    print(f"\n=== {vtype} VAF > {args.vaf_threshold} 的位点数: {len(kept)} "
          f"(总位点 {total_sites}) ===")
    print(f"{'rp':>4} {'ref':>3} {'alt':>3} {'n_events':>9} {'VAF%':>6} "
          f"{'q_min':>5} {'q5':>4} {'q10':>4} {'q25':>4} {'q50':>4} "
          f"{'q75':>4} {'q_max':>5} {'q_mean':>7}")
    for rp, rec, vaf in kept:
        q_min, q_p5, q_p10, q_p25, q_p50, q_p75, q_max, q_mean = site_q_stats(rec["qs"])
        ref_base, alt_base = dominant_allele(rec)
        alt = alt_base if len(alt_base) <= 3 else alt_base[:3] + ".."
        print(f"{rp:4d} {ref_base:>3} {alt:>3} {rec['n']:9d} {100*vaf:6.1f} "
              f"{q_min:5.0f} {q_p5:4.0f} {q_p10:4.0f} {q_p25:4.0f} {q_p50:4.0f} "
              f"{q_p75:4.0f} {q_max:5.0f} {q_mean:7.2f}")

    ov = overall_stats(all_qs)
    print(f"\n=== {vtype} 总体基本统计 (N={int(ov[0]):,}) ===")
    print(f"{'N':>5} = {int(ov[0]):,}")
    for col, v in zip(OVERALL_COLUMNS[2:], ov[1:]):
        print(f"{col:>5} = {v:.2f}")


if __name__ == "__main__":
    main()
