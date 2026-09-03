#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3：统计变体位点（SNP/INS/DEL）的 Q 值分布。

读取一个或多个 03_variant_calling.py 输出的变体 TSV，按 --type 过滤后汇总
该类型位点的碱基质量 (base_q) 分布，输出：
    - {out_prefix}_q_dist.png      碱基质量直方图（合并 + 各 barcode）
    - {out_prefix}_q_summary.tsv   分布统计（N/min/P1/P5/.../P99/max/mean/std）
    - {out_prefix}_q_per_pos.png   沿参考位点的 baseQ 分位数分布

用法:
    python3 04_snp_q_dist.py 'results/variant/*.variants.tsv' \
        --type INS --out-prefix results/qdist/Barcode01_INS_qdist
"""

from __future__ import annotations

import argparse
import glob
import logging
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)

PERCENTILES = [1, 5, 10, 25, 50, 75, 90, 95, 99]


def read_snp_tsv(path: str, vtype: str) -> list[tuple[str, int, int]]:
    """读取变体 TSV，返回 (barcode, ref_pos_1based, base_q)，仅保留 variant_type==vtype 的行。"""
    rows = []
    with open(path, "r") as f:
        header = f.readline()
        if header and not header.startswith("barcode"):
            # 无表头，第一行也是数据
            f.seek(0)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[8] != vtype:
                continue
            rows.append((parts[0], int(parts[2]), int(parts[5])))
    return rows


def stats_rows(base_q: list[int]) -> list[float]:
    """计算 base_q 的 (N,min,p1,...,p99,max,mean,std)。"""
    arr = np.asarray(base_q, dtype=float)
    if arr.size == 0:
        return [0.0] * (1 + 2 + len(PERCENTILES) + 2)
    qs = np.percentile(arr, PERCENTILES)
    return [float(arr.size), float(arr.min()),
            *[float(q) for q in qs],
            float(arr.max()), float(arr.mean()), float(arr.std(ddof=0))]


STAT_COLUMNS = (["group", "N", "min"]
                + [f"P{p}" for p in PERCENTILES]
                + ["max", "mean", "std"])


def plot_histogram(base_q: list[int], title: str, out_path: str) -> None:
    """绘制 baseQ 直方图。"""
    arr = np.asarray(base_q, dtype=int)
    lo = int(arr.min())
    hi = int(arr.max())
    bins = list(range(lo, hi + 2))
    plt.figure(figsize=(12, 6))
    plt.hist(arr, bins=bins, color="#4c72b0", edgecolor="white", align="left")
    plt.xlabel("Variant base quality (PHRED)")
    plt.ylabel("Variant event count")
    plt.title(title)
    plt.xticks(range(lo, hi + 1, max(1, (hi - lo) // 25)))
    mean_q = float(arr.mean())
    median_q = float(np.median(arr))
    stats_text = (f"N = {arr.size:,}\n"
                  f"mean = {mean_q:.1f}\n"
                  f"median = {median_q:.0f}")
    plt.text(0.02, 0.95, stats_text, transform=plt.gca().transAxes,
             va="top", ha="left", fontsize=11,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    log.info("直方图已保存: %s", out_path)


def plot_per_position(base_q_by_pos: dict[int, list[int]], out_path: str) -> None:
    """沿参考位点绘制 baseQ 分位数线。"""
    positions = sorted(base_q_by_pos)
    if not positions:
        log.warning("无数据，跳过 per-position 图")
        return
    centers = list(positions)
    quantile_fracs = [0.05, 0.25, 0.50, 0.75, 0.95]
    q_labels = ["P5", "P25", "P50", "P75", "P95"]
    colors = ["purple", "darkorange", "green", "blue", "red"]

    plt.figure(figsize=(14, 5))
    for frac, label, color in zip(quantile_fracs, q_labels, colors):
        y = []
        for p in positions:
            vals = sorted(base_q_by_pos[p])
            if len(vals) >= 2:
                idx = int(len(vals) * frac)
                y.append(vals[min(idx, len(vals) - 1)])
            else:
                y.append(float("nan"))
        plt.plot(centers, y, linestyle="--", linewidth=1.2, alpha=0.8,
                 color=color, label=label)
    mean_vals = [sum(v) / len(v) if v else float("nan")
                 for v in (base_q_by_pos[p] for p in positions)]
    plt.plot(centers, mean_vals, color="darkred", linestyle="-",
             linewidth=1.5, alpha=0.9, label="mean")
    plt.xlabel("Reference position (1-based)")
    plt.ylabel("Variant base quality (PHRED)")
    plt.title("Variant base quality along reference")
    plt.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    log.info("per-position 图已保存: %s", out_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="统计变体位点 Q 值分布。")
    parser.add_argument("snp_tsv", nargs="+", help="变体 TSV 路径，支持通配符")
    parser.add_argument("--type", choices=["SNP", "INS", "DEL"], default="SNP",
                        help="变体类型 (默认 SNP)")
    parser.add_argument("--out-prefix", default="q_dist", help="输出前缀")
    parser.add_argument("--title", default=None, help="图标题 (默认自动)")
    args = parser.parse_args()
    vtype = args.type
    title = args.title or f"{vtype} baseQ distribution"

    tsv_paths = []
    for pat in args.snp_tsv:
        tsv_paths.extend(glob.glob(pat))
    tsv_paths = [p for p in tsv_paths if os.path.exists(p)]
    if not tsv_paths:
        log.error("未找到变体 TSV 文件")
        sys.exit(1)

    all_rows: list[tuple[str, int, int]] = []
    for p in sorted(tsv_paths):
        rows = read_snp_tsv(p, vtype)
        log.info("%s: %d %s 事件", p, len(rows), vtype)
        all_rows.extend(rows)
    log.info("总 %s 事件: %d", vtype, len(all_rows))

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)

    # 按 barcode 分组
    by_barcode: dict[str, list[int]] = {}
    by_pos: dict[int, list[int]] = {}
    for barcode, rpos, bq in all_rows:
        by_barcode.setdefault(barcode, []).append(bq)
        by_pos.setdefault(rpos, []).append(bq)

    all_bq = [r[2] for r in all_rows]

    # 统计表
    stat_path = f"{args.out_prefix}_q_summary.tsv"
    with open(stat_path, "w") as f:
        f.write("\t".join(STAT_COLUMNS) + "\n")
        f.write("\t".join(["combined"] + [str(v) for v in stats_rows(all_bq)]) + "\n")
        for barcode in sorted(by_barcode):
            f.write("\t".join([barcode] + [str(v) for v in stats_rows(by_barcode[barcode])]) + "\n")
    log.info("统计表已保存: %s", stat_path)

    # 直方图（合并 + 各 barcode）
    plot_histogram(all_bq, title, f"{args.out_prefix}_q_dist.png")
    for barcode, bq in sorted(by_barcode.items()):
        plot_histogram(bq, f"{barcode} {vtype} baseQ distribution",
                       f"{args.out_prefix}_{barcode}_q_dist.png")

    # per-position 分布
    plot_per_position(by_pos, f"{args.out_prefix}_q_per_pos.png")

    # 打印合并汇总
    s = stats_rows(all_bq)
    print(f"\n=== {vtype} baseQ distribution summary (combined) ===")
    print(f"{'N':>6} = {int(s[0]):,}")
    for n, v in zip(STAT_COLUMNS[2:], s[1:]):
        print(f"{n:>6} = {v:.2f}")


if __name__ == "__main__":
    main()
