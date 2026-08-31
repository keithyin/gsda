"""绘制 BAM 文件中 qual（per-base Q）的分布图。

支持一次输入多个 bam，所有分布画在同一张图上，用 legend 区分。
纵轴为 prob（归一化到总和为 1），不是 count。
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

MAX_Q = 93  # PHRED 质量值上限
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


def plot_qual_dist(
    qual_hists: Dict[str, Tuple[np.ndarray, np.ndarray]],
    level: str,
    output_file: str,
    log_y: bool = False,
    figsize: Tuple[float, float] = (12, 7),
    dpi: int = 300,
):
    """将多个 bam 的 prob 分布叠加画在一张图上"""
    sns.set_theme(style="whitegrid", font_scale=1.1)
    fig, ax = plt.subplots(figsize=figsize)
    colors = sns.color_palette("deep", len(qual_hists))

    ymax = 0.0
    xmax = 0.0
    for idx, (name, (edges, counts)) in enumerate(qual_hists.items()):
        total = counts.sum()
        if total == 0:
            print(f"⚠️  {name}: 无有效 qual 数据，跳过")
            continue
        prob = counts / total  # 画 prob，不画 count
        ymax = max(ymax, prob.max())
        nonzero = np.flatnonzero(counts)
        xmax = max(xmax, float(edges[nonzero[-1] + 1]))
        mean_q = float(
            (prob * (edges[:-1] + edges[1:]) / 2).sum())
        ax.stairs(
            prob,
            edges,
            label=f"{name}  (mean Q={mean_q:.1f})",
            color=colors[idx],
            linewidth=1.8,
            fill=False,
        )
        print(f"📊 {name}: bins_sum={int(counts.sum()):,} | mean Q={mean_q:.2f}")

    if log_y:
        ax.set_yscale("log")
        ax.set_ylim(top=ymax, bottom=ymax / 1e6)
    ax.set_xlabel("Base Quality (PHRED)" if level == "base" else "Mean Read Quality (PHRED)")
    ax.set_ylabel("Probability")
    ax.set_title(f"Quality Distribution ({level} level, probability)")
    ax.legend(frameon=True, loc="upper right")
    ax.set_xlim(0, xmax)
    ax.xaxis.set_major_locator(plt.MaxNLocator(20))

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
        help="统计的质量值上限，更高的 qual 归入该 bin。Default: 93",
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
