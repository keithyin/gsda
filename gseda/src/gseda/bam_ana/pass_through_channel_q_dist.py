import pysam
from tqdm import tqdm
import os
import argparse
import pathlib
from multiprocessing import cpu_count
import matplotlib.pyplot as plt


def read_bam_info(
    bam_file: str, np_thr: int
):
    rq_counter = {}

    with pysam.AlignmentFile(
        filename=bam_file, mode="rb", check_sq=False, threads=cpu_count()
    ) as reader:
        for record in tqdm(reader.fetch(until_eof=True), desc=f"reading {bam_file}"):

            np = record.get_tag("np")
            if np > np_thr:
                continue
            rq = float(record.get_tag("rq")) * 100
            rq = int(rq)
            rq_counter.setdefault(rq, 0)
            rq_counter[rq] += 1

    return rq_counter


def plot_rq_histogram(rq_counter, title="Read Quality Distribution", save_path=None):
    """
    绘制RQ值柱状图
    Args:
        rq_counter: 从read_bam_info函数返回的字典
        title: 图表标题
        save_path: 保存图片的路径（可选）
    """
    # 筛选50~100的RQ值
    filtered_rq = {rq: count for rq, count in rq_counter.items()
                   if 50 <= rq <= 100}

    if not filtered_rq:
        print("警告：没有找到50~100范围内的RQ值")
        return

    # 获取所有50~100的RQ值，包括没有计数的
    all_rq_values = list(range(50, 101))
    counts = [filtered_rq.get(rq, 0) for rq in all_rq_values]

    # 创建图表
    plt.figure(figsize=(12, 6))

    # 创建柱状图
    bars = plt.bar(all_rq_values, counts, width=0.8,
                   edgecolor='black', alpha=0.7)

    # 设置颜色渐变（从低到高）
    for i, bar in enumerate(bars):
        # 根据RQ值设置颜色（从蓝到红）
        normalized_rq = (all_rq_values[i] - 50) / 50  # 归一化到0-1
        bar.set_facecolor(plt.cm.viridis(normalized_rq))

    # 添加数值标签（只显示非零值）
    for i, count in enumerate(counts):
        if count > 0:
            plt.text(all_rq_values[i], count + max(counts)*0.01,
                     str(count), ha='center', va='bottom', fontsize=8)

    # 设置图表属性
    plt.xlabel('Read Quality (RQ)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')

    # 设置x轴刻度
    plt.xticks(range(50, 101, 5))
    plt.xlim(49, 101)

    # 添加网格
    plt.grid(True, axis='y', alpha=0.3, linestyle='--')

    # 添加统计信息
    total_reads = sum(counts)

    stats_text = f"Total Reads: {total_reads:,}"
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    # 保存或显示图表
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图表已保存到: {save_path}")

    plt.show()


def main(args):
    rq_counter = read_bam_info(bam_file=args.bam, np_thr=args.np_thr)
    p = pathlib.Path(args.bam)
    save_path = str(p.parent.joinpath(f"{p.stem}.np_thr{args.np_thr}.png"))

    plot_rq_histogram(rq_counter, save_path=save_path)


def main_cli():
    parser = argparse.ArgumentParser(prog="bam basic stat")
    parser.add_argument("bam", type=str)

    parser.add_argument(
        "--np-thr",
        type=float,
        default=0,
        help="only the rq ≥ min-rq will be considered",
        dest="np_thr",
    )
    args = parser.parse_args()
    main(args=args)


if __name__ == "__main__":
    main_cli()
