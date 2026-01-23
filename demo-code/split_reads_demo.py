import mappy
import re
import pysam
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np


def plot_len_dist(lengths, name='length_distribution.jpg'):
    plt.figure(figsize=(10, 6))

    # 绘制直方图（分布图）
    # bins参数控制直方图的柱子数量，可以根据需要调整
    plt.hist(lengths, bins=20, edgecolor='black', alpha=0.7, color='skyblue')

    # 添加标题和标签
    plt.title('Length Distribution', fontsize=16, fontweight='bold')
    plt.xlabel('Length', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)

    # 添加网格（可选）
    plt.grid(True, alpha=0.3, linestyle='--')

    # 如果需要显示统计信息
    mean_length = np.mean(lengths)
    median_length = np.median(lengths)
    plt.axvline(mean_length, color='red', linestyle='--',
                linewidth=2, label=f'Mean: {mean_length:.2f}')
    plt.axvline(median_length, color='green', linestyle='--',
                linewidth=2, label=f'Median: {median_length:.2f}')
    plt.legend()

    # 调整布局
    plt.tight_layout()

    # 保存为JPEG
    plt.savefig(name, dpi=300, format='jpeg')
    print(f"图像已保存为 {name}")


def detailed_cigar_stats(cigar_string):
    """提供详细的CIGAR统计信息"""
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    stats = {'=': 0, 'X': 0, 'I': 0, 'D': 0}

    for length_str, operation in matches:
        length = int(length_str)
        if operation in stats:
            stats[operation] += length

    total_aligned = stats['='] + stats['X'] + stats['I'] + stats['D']
    identity = stats['='] / total_aligned if total_aligned > 0 else 0

    print("\n详细统计信息:")
    print(f"精确匹配(=): {stats['=']} 个碱基")
    print(f"错配(X): {stats['X']} 个碱基")
    print(f"插入(I): {stats['I']} 个碱基")
    print(f"删除(D): {stats['D']} 个碱基")
    print(f"总比对长度: {total_aligned} 个碱基")
    print(f"Identity: {identity:.6f}")

    return stats


def calculate_identity_from_cigar(cigar_string):
    """
    根据CIGAR字符串计算identity

    Parameters:
    cigar_string (str): CIGAR字符串, 如 "370=1I89=1D26=1I208=1D1250=1D578=1D2=2X1=1I578=1I1248=1I210=1D23=1I90=1D373="

    Returns:
    float: identity值
    """
    # 使用正则表达式解析CIGAR字符串
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    # 初始化计数器
    match_count = 0      # 匹配的碱基数（=）
    mismatch_count = 0   # 错配的碱基数（X）
    total_aligned = 0    # 总比对长度（= + X + I + D）

    for length_str, operation in matches:
        length = int(length_str)

        if operation == '=':
            # 精确匹配
            match_count += length
            total_aligned += length
        elif operation == 'X':
            # 错配
            mismatch_count += length
            total_aligned += length
        elif operation == 'I':
            # 插入：参考序列没有，查询序列有
            total_aligned += length
        elif operation == 'D':
            # 删除：参考序列有，查询序列没有
            total_aligned += length

    # 计算identity
    if total_aligned == 0:
        return 0.0

    identity = match_count / total_aligned
    return identity


def try_split_seq(seq: str):
    seq_len = len(seq)
    if seq_len > 20000:
        return None
    aligner = mappy.Aligner(seq=seq, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    for hit in aligner.map(seq):
        if hit.r_st == 0 and hit.r_en == seq_len and hit.strand == 1 and calculate_identity_from_cigar(hit.cigar_str) >= 0.99999:
            continue
        # if hit.r_st == hit.q_st and hit
        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage = (hit.r_en - hit.r_st) / seq_len
        if identity > 0.85 and coverage > 0.85:
            return (hit.r_en + hit.r_st) // 2
    return None


def try_split_seq_2(seq: str):
    seq_len = len(seq)
    if seq_len > 20000:
        return None

    seq1 = seq[:(seq_len // 2)]
    seq2 = seq[(seq_len // 2):]

    aligner = mappy.Aligner(seq=seq1, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    for hit in aligner.map(seq2):

        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage1 = (hit.r_en - hit.r_st) / len(seq1)
        coverage2 = (hit.q_en - hit.q_st) / len(seq2)
        if identity > 0.85 and coverage1 > 0.85 and coverage2 > 0.85:

            print("ctg:{}\tref_start:{}\tref_end:{}\tq_start:{}\tq_end:{}\tstrand:{}\tref_len:{}\tquery_len:{}. ideneity:{:.4f}".format(
                hit.ctg, hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.strand, len(seq1), len(seq2), calculate_identity_from_cigar(hit.cigar_str)))

            return seq_len // 2
    return None


def main():
    run_name = "20251224_250302Y0004_Run0002"
    fpath = f"/data1/ccs_data/20260109-saisuofei-resplit/{run_name}/barcodes_reads_fastq_amplicon/Single-1_Double-257.fastq"
    tot = 0
    splitted = 0

    lengths = []

    ori_lengths = []

    with pysam.FastxFile(fpath) as f:
        for record in tqdm(f, desc="processing"):
            seq = record.sequence
            ori_lengths.append(len(seq))
            position = try_split_seq(seq)
            tot += 1
            if position is not None:
                splitted += 1
                lengths.extend([position, len(seq) - position])
            else:
                lengths.append(len(seq))

    print(f"split_ratio: {splitted / tot * 100:.4f}%")
    plot_len_dist(lengths, name=f"{run_name}-length_distribution.splited.jpg")
    plot_len_dist(ori_lengths, name=f"{run_name}-length_distribution.ori.jpg")


if __name__ == "__main__":
    main()
