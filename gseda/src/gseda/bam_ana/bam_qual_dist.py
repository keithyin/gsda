"""绘制 BAM 文件中 qual（per-base Q）的分布图。

支持一次输入多个 bam，所有分布画在同一张图上（柱状图并排），用 legend 区分。
纵轴为 prob（归一化到总和为 1），不是 count；下半图是对应的累积概率分布（CDF）。
"""

import argparse
import os
import pathlib
from multiprocessing import cpu_count
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pysam
import seaborn as sns
from tqdm import tqdm
import pathlib

MAX_Q = 60  # PHRED 质量值上限
READ_BIN_WIDTH = 0.5  # read 级别分布的分箱宽度


def collect_qual(
    bam_file: str, level: str, max_q: int = MAX_Q
) -> Tuple[np.ndarray, np.ndarray]:
    """单次遍历 BAM，统计 qual 的直方图。

    Args:
        bam_file: bam 路径
        level: "base" 统计每个碱基的 qual；"read" 统计每条 read 的平均 qual
        max_q: 质量值上限，超出部分并入最后一个 bin

    Returns:
        (edges, counts)  edges 长度为 len(counts) + 1
    """
    assert level in ("base", "read")
    bin_width = 1.0 if level == "base" else READ_BIN_WIDTH
    n_bins = int(np.ceil((max_q + 1) / bin_width))
    counts = np.zeros(n_bins, dtype=np.int64)

    n_no_qual = 0
    with pysam.AlignmentFile(
        filename=bam_file, mode="rb", check_sq=False, threads=cpu_count() // 2
    ) as reader:
        desc = f"reading {os.path.basename(bam_file)} ({level} qual)"
        for record in tqdm(reader.fetch(until_eof=True), desc=desc):
            qual = record.query_qualities
            if qual is None:
                n_no_qual += 1
                continue

            if level == "base":
                q = np.frombuffer(qual, dtype=np.int8).astype(np.int64)
                # 负值（无质量 '-'）与超上限的一律夹到 [0, max_q]
                np.clip(q, 0, max_q, out=q)
                idx = q
            else:
                q = np.frombuffer(qual, dtype=np.int8)
                if len(q) == 0:
                    continue
                mean_q = np.clip(q.astype(np.float64).mean(), 0, max_q)
                idx = np.array([int(mean_q / bin_width)], dtype=np.int64)

            counts += np.bincount(idx, minlength=n_bins)

    if n_no_qual:
        print(f"⚠️  {os.path.basename(bam_file)}: {n_no_qual} 条 record 无 qual，已跳过")

    edges = np.arange(n_bins + 1, dtype=np.float64) * bin_width
    return edges, counts


def weighted_median_from_hist(edges: np.ndarray, counts: np.ndarray) -> float:
    """由直方图求加权中位数（bin 内线性插值）"""
    total = counts.sum()
    cum = np.cumsum(counts)
    half = total / 2.0
    i = min(int(np.searchsorted(cum, half, side="left")), len(counts) - 1)
    prev = cum[i - 1] if i > 0 else 0
    lo, hi = float(edges[i]), float(edges[i + 1])
    if counts[i] == 0:
        return (lo + hi) / 2
    return float(lo + (half - prev) / counts[i] * (hi - lo))


def plot_qual_dist(
    qual_hists: Dict[str, Tuple[np.ndarray, np.ndarray]],
    level: str,
    output_file: str,
    log_y: bool = False,
    figsize: Tuple[float, float] = (12, 9),
    dpi: int = 300,
):
    """柱状图画 prob 分布，下方共享 x 轴画累积概率（CDF）"""
    series = []
    for name, (edges, counts) in qual_hists.items():
        total = counts.sum()
        if total == 0:
            print(f"⚠️  {name}: 无有效 qual 数据，跳过")
            continue
        prob = counts / total  # 画 prob，不画 count
        centers = (edges[:-1] + edges[1:]) / 2
        mean_q = float((prob * centers).sum())
        median_q = weighted_median_from_hist(edges, counts)
        series.append((name, edges, prob, np.cumsum(prob), mean_q, median_q))
        print(
            f"📊 {name}: n={int(total):,} | mean Q={mean_q:.2f} | median Q={median_q:.2f}")

    if not series:
        print("❌ 所有 bam 都没有 qual 数据，无法画图")
        return

    sns.set_theme(style="whitegrid", font_scale=1.1)
    fig, (ax_bar, ax_cdf) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.06},
    )
    colors = sns.color_palette("deep", len(series))
    bin_width = float(series[0][1][1] - series[0][1][0])
    bar_width = bin_width / len(series) * 0.85
    xmax = max(float(s[1][np.flatnonzero(s[2] > 0)[-1] + 1]) for s in series)

    for idx, (name, edges, prob, cdf, mean_q, median_q) in enumerate(series):
        # 多文件在 bin 内并排（居中），避免柱子互相遮挡
        offset = (bin_width - len(series) * bar_width) / 2 + idx * bar_width
        ax_bar.bar(
            edges[:-1] + offset,
            prob,
            width=bar_width,
            align="edge",
            color=colors[idx],
            alpha=0.85,
            linewidth=0,
            label=f"{name}\n    mean Q={mean_q:.1f}, median Q={median_q:.1f}",
        )
        ax_bar.axvline(mean_q, color=colors[idx], linestyle="--", linewidth=1.4)
        ax_bar.axvline(median_q, color=colors[idx],
                       linestyle="-.", linewidth=1.4)

        # 用 steps 而不是 stairs：stairs 会在两端补一条回到 baseline 的竖线
        ax_cdf.plot(np.append(edges[:-1], edges[-1]), np.append(cdf, cdf[-1]),
                    drawstyle="steps-post", color=colors[idx], linewidth=1.8)
        ax_cdf.plot(median_q, 0.5, "o", color=colors[idx], markersize=6)

    ymax = max(s[2].max() for s in series)
    if log_y:
        ax_bar.set_yscale("log")
        ax_bar.set_ylim(top=ymax * 1.1, bottom=ymax / 1e6)
        ax_bar.set_ylabel("Probability (log)")
    else:
        # 顶部留白，避免最高的柱子顶到边框压住 legend
        ax_bar.set_ylim(bottom=0, top=ymax * 1.1)
        ax_bar.set_ylabel("Probability")
    ax_bar.set_title(
        f"Quality Distribution ({level} level, probability) — mean: dashed, median: dash-dot"
    )

    ax_cdf.axhline(0.5, color="grey", linestyle=":", linewidth=1)
    ax_cdf.set_ylabel("Cumulative Probability")
    ax_cdf.set_ylim(0, 1.03)
    ax_cdf.set_xlabel(
        "Base Quality (PHRED)" if level == "base" else "Mean Read Quality (PHRED)")
    # ● 标注中位数位置（● = 该曲线 CDF 与 0.5 的交点）
    ax_cdf.plot([], [], "o", color="grey", label="median (CDF = 0.5)")
    ax_cdf.legend(frameon=True, loc="lower right")

    ax_bar.legend(frameon=True, loc="upper right", fontsize=9)
    ax_bar.set_xlim(0, xmax)
    ax_bar.xaxis.set_major_locator(plt.MaxNLocator(20))

    _save_figure(fig, output_file, dpi)
    print(f"✅ 分布图已保存: {os.path.abspath(output_file)}")


def _save_figure(fig: plt.Figure, output_file: str, dpi: int):
    output_dir = os.path.dirname(os.path.abspath(output_file))
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def default_output_path(bams: List[str], level: str) -> str:
    p = pathlib.Path(bams[0])
    if len(bams) > 1:
        return str(p.parent / f"{p.stem}.multi_{len(bams)}_qual_dist.{level}.png")
    return str(p.parent / f"{p.stem}.qual_dist.{level}.png")


def main(args):
    import shutil
    qual_hists = {}
    for bam in args.bams:
        name = os.path.basename(bam)
        if name in qual_hists:  # 同名不同目录时加序号区分
            name = f"{args.bams.index(bam)}-{name}"
        qual_hists[name] = collect_qual(bam, level=args.level, max_q=args.max_q)

    output_file = args.output or default_output_path(args.bams, args.level)
    plot_qual_dist(
        qual_hists,
        level=args.level,
        output_file=output_file,
        log_y=args.log_y,
        dpi=args.dpi,
    )
    
    gsda_tmp_dir = pathlib.Path("/root/projects/gsda/tmp-data-dir")
    if gsda_tmp_dir.exists():
        shutil.copy(output_file, gsda_tmp_dir / pathlib.Path(output_file).name)
    
    


def main_cli():
    parser = argparse.ArgumentParser(
        prog="bam-qual-dist",
        description="Plot qual (PHRED) distribution of BAM files, y axis is probability",
    )
    parser.add_argument("bams", nargs="+", type=str, help="One or more BAM file paths")
    parser.add_argument(
        "--level",
        type=str,
        default="base",
        choices=["base", "read"],
        help="base: 逐碱基 qual；read: 每条 read 的平均 qual。Default: base",
    )
    parser.add_argument(
        "--max-q",
        type=int,
        default=MAX_Q,
        help=f"统计的质量值上限，更高的 qual 归入该 bin。Default: {MAX_Q}",
        dest="max_q",
    )
    parser.add_argument(
        "--log-y", action="store_true", help="纵轴取 log 尺度", dest="log_y"
    )
    parser.add_argument("--output", type=str, default=None, help="输出图片路径")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    main(args=args)


if __name__ == "__main__":
    main_cli()
