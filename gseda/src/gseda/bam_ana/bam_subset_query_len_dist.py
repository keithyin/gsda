import pysam
from tqdm import tqdm
import os
import argparse
import numpy as np
from multiprocessing import cpu_count
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
from typing import Dict, List, Tuple


def extract_multi_rq_length_counts(
    bam_file: str,
    channel_tag: str,
    rq_thresholds: List[float]
) -> Dict[float, Dict[int, int]]:
    """
    单次遍历 BAM，同时统计多个 rq 阈值的长度分布

    返回: {threshold: {length: count}}
    """
    assert channel_tag in ("ch", "zm")
    # 按阈值降序排序（优化：高阈值先检查，但无法短路因阈值非连续）
    thresholds_sorted = sorted(set(rq_thresholds), reverse=True)
    accumulators = {f"{t:.5f}": {} for t in thresholds_sorted}

    with pysam.AlignmentFile(
        filename=bam_file, mode="rb", check_sq=False, threads=cpu_count() // 2
    ) as reader:
        desc = f"Processing {os.path.basename(bam_file)} for {len(thresholds_sorted)} rq thresholds"
        for record in tqdm(reader.fetch(until_eof=True), desc=desc):
            # 计算 rq 值
            if record.has_tag("rq"):
                rq = float(record.get_tag("rq"))
            elif record.has_tag("cq"):
                rq = 1 - 10 ** (float(record.get_tag("cq")) / -10)
            else:
                rq = 1.0  # 无质量标签视为完美质量

            # 获取序列长度
            length = len(record.query_sequence)

            # 为所有满足条件的阈值累加计数
            for threshold in thresholds_sorted:
                if rq >= threshold:
                    threshold = f"{threshold:.5f}"
                    accumulators[threshold][length] = accumulators[threshold].get(
                        length, 0) + 1

    return accumulators


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """安全计算加权中位数，避免内存爆炸"""
    sorted_idx = np.argsort(values)
    sorted_vals = values[sorted_idx]
    sorted_weights = weights[sorted_idx]
    cumsum = np.cumsum(sorted_weights)
    cutoff = cumsum[-1] / 2.0
    return float(sorted_vals[np.searchsorted(cumsum, cutoff)])


def plot_single_threshold(
    length_counts: Dict[int, int],
    bam_file: str,
    min_rq: float,
    output_file: str,
    bins: int = 1000,
    figsize: tuple = (10, 8),
    dpi: int = 300
):
    """单阈值分布：双子图（线性+对数坐标）"""
    if not length_counts:
        raise ValueError("未提取到任何有效序列，请检查 BAM 文件或过滤参数")

    # 转换为数组
    lengths = np.array(list(length_counts.keys()))
    counts = np.array(list(length_counts.values()))

    # 计算统计量
    total_reads = counts.sum()
    mean_len = np.average(lengths, weights=counts)
    median_len = weighted_median(lengths, counts)

    # 创建分箱
    min_len, max_len = lengths.min(), lengths.max()
    bins_linear = np.linspace(min_len, max_len, bins + 1)
    bins_log = np.logspace(np.log10(max(1, min_len)),
                           np.log10(max_len), bins + 1)

    # 绘图
    sns.set_theme(style="whitegrid", font_scale=1.2)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharey=False)

    # 线性坐标
    sns.histplot(
        x=lengths, weights=counts, bins=bins_linear,
        ax=ax1, color="steelblue", edgecolor="white", linewidth=0.5
    )
    ax1.set_title(f"Sequence Length Distribution (Linear Scale)\n"
                  f"BAM: {os.path.basename(bam_file)} | rq ≥ {min_rq}\n"
                  f"Total reads: {total_reads:,} | Mean: {mean_len:.1f} bp | Median: {median_len:.1f} bp",
                  fontsize=13, fontweight="bold")
    ax1.set_xlabel("Sequence Length (bp)", fontsize=12)
    ax1.set_ylabel("Read Count", fontsize=12)
    ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # 对数坐标
    sns.histplot(
        x=lengths, weights=counts, bins=bins_log,
        ax=ax2, color="coral", edgecolor="white", linewidth=0.5
    )
    ax2.set_xscale('log')
    ax2.set_title("Sequence Length Distribution (Log Scale)",
                  fontsize=13, fontweight="bold")
    ax2.set_xlabel("Sequence Length (bp, log scale)", fontsize=12)
    ax2.set_ylabel("Read Count", fontsize=12)
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.tight_layout()
    _save_figure(fig, output_file, dpi)

    print(f"✅ 单阈值分布图已保存: {os.path.abspath(output_file)}")
    print(
        f"📊 统计: rq≥{min_rq} | Reads={total_reads:,} | Mean={mean_len:.1f} bp | Median={median_len:.1f} bp")


def plot_multi_threshold(
    multi_counts: Dict[float, Dict[int, int]],
    rq_thresholds: List[float],
    output_file: str,
    bins: int = 100,
    figsize: tuple = (12, 8),
    dpi: int = 300
):
    """多阈值分布：单图叠加（对数坐标）"""
    # 按阈值降序排列（高质量在前）
    thresholds = sorted(rq_thresholds, reverse=False)

    # 全局长度范围（用于统一分箱）

    # 创建图形
    sns.set_theme(style="whitegrid", font_scale=1.2)
    fig, axs = plt.subplots(figsize=figsize, nrows=len(thresholds), ncols=1)

    # 配色方案（从高质量到低质量）
    colors = sns.color_palette("husl", len(thresholds))

    # 绘制每个阈值的分布
    stats_summary = []
    for idx, threshold in enumerate(thresholds):

        ax = axs[idx]
        counts_dict = multi_counts.get(f"{threshold:.5f}", {})
        if not counts_dict:
            print(f"⚠️  阈值 rq≥{threshold} 无有效数据，跳过")
            continue

        lengths = np.array(list(counts_dict.keys()))
        counts = np.array(list(counts_dict.values()))

        # ✅ 每个 subplot 自己的 x 轴范围（基于 length）
        cur_min = max(1, lengths.min())
        cur_max = lengths.max()
        bins_linear = np.linspace(cur_min, cur_max, bins + 1)

        # 统计量
        total = counts.sum()
        mean = np.average(lengths, weights=counts)
        median = weighted_median(lengths, counts)
        stats_summary.append((threshold, total, mean, median))

        # 画图
        ax.hist(
            lengths,
            bins=bins_linear,
            weights=counts,
            histtype='step',
            linewidth=2.0,
            color=colors[idx],
            alpha=0.9
        )

        # ===== 字体 & 轴设置 =====
        ax.set_ylabel(
            "Read Count",
            fontsize=9,
            fontweight="bold"
        )

        # 子图标题（小而清晰）
        ax.text(
            0.02, 0.95,
            f"rq ≥ {threshold:.3f} | n = {total:,}",
            transform=ax.transAxes,
            fontsize=5,
            fontweight="bold",
            va="top",
            ha="left"
        )

        ax.set_xlim(cur_min, cur_max)

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=8
        )

        ax.grid(
            True,
            which="both",
            linestyle="--",
            linewidth=0.6,
            alpha=0.6
        )

        # ✅ 只有最后一个 subplot 显示 x label
        if idx == len(thresholds) - 1:
            ax.set_xlabel(
                "Sequence Length (bp)",
                fontsize=10,
                fontweight="bold"
            )
        else:
            ax.set_xlabel("")

    # plt.tight_layout()
    _save_figure(fig, output_file, dpi)

    # 打印统计摘要
    print(f"\n✅ 多阈值分布图已保存: {os.path.abspath(output_file)}")
    print("\n📊 各阈值统计摘要:")
    print(f"{'Threshold':<12} {'Total Reads':<15} {'Mean Length':<15} {'Median Length':<15}")
    print("-" * 60)
    for threshold, total, mean, median in stats_summary:
        print(
            f"rq ≥ {threshold:<5.2f} {total:>14,} {mean:>14.0f} bp {median:>14.0f} bp")


def _save_figure(fig: plt.Figure, output_file: str, dpi: int):
    """保存图形的辅助函数"""
    output_dir = os.path.dirname(os.path.abspath(output_file))
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def generate_output_path(bam_path: str, rq_values: List[float], is_multi: bool) -> str:
    """生成输出文件路径"""
    p = pathlib.Path(bam_path)
    if is_multi:
        # 多阈值：用下划线连接阈值（保留2位小数）
        rq_str = "_".join(
            [f"{v:.2f}" for v in sorted(rq_values, reverse=True)])
        return str(p.parent / f"{p.stem}.rq_{rq_str}.length_dist.png")
    else:
        # 单阈值：保留原命名风格
        rq_str = f"{rq_values[0]:.2f}" if rq_values[0] > 0 else "all"
        return str(p.parent / f"{p.stem}.rq_{rq_str}.length_dist.png")


def main(args):
    # 确保阈值列表非空
    if not args.min_rq:
        args.min_rq = [0.0]

    for bam_path in args.bams:
        # 单次提取所有阈值的分布
        multi_counts = extract_multi_rq_length_counts(
            bam_path,
            args.channel_tag,
            args.min_rq
        )

        # 根据阈值数量选择绘图模式
        is_multi = len(args.min_rq) > 1
        output_path = generate_output_path(bam_path, args.min_rq, is_multi)

        if is_multi:
            plot_multi_threshold(
                multi_counts,
                args.min_rq,
                output_path,
                bins=args.bins,
                figsize=(12, 8),
                dpi=args.dpi
            )
        else:
            # 单阈值：取第一个阈值的数据
            single_counts = multi_counts.get(args.min_rq[0], {})
            plot_single_threshold(
                single_counts,
                bam_path,
                args.min_rq[0],
                output_path,
                bins=args.bins,
                figsize=(10, 8),
                dpi=args.dpi
            )


def main_cli():
    parser = argparse.ArgumentParser(
        prog="bam-length-dist",
        description="BAM file sequence length distribution analysis with rq filtering"
    )
    parser.add_argument("bams", nargs="+", type=str,
                        help="One or more BAM file paths")
    parser.add_argument(
        "--min_rq",
        nargs="+",
        type=float,
        default=[0.0, 0.95, 0.99, 0.999],
        help="One or more rq thresholds (e.g., --min_rq 0.7 0.8 0.9). Default: 0.0 (no filtering)"
    )
    parser.add_argument(
        "--channel_tag",
        type=str,
        default="ch",
        choices=["ch", "zm"],
        help="Channel tag to validate (ch or zm). Default: ch"
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=1000,
        help="Number of histogram bins. Default: 100 (multi-threshold) / 1000 (single-threshold)"
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image DPI. Default: 300"
    )
    args = parser.parse_args()

    # 单阈值时使用更细的分箱
    if len(args.min_rq) == 1 and args.bins == 100:
        args.bins = 1000

    main(args)


if __name__ == "__main__":
    main_cli()
