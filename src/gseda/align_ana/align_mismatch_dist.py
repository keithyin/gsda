import matplotlib.pyplot as plt
import mappy
import pathlib
import os
import sys
import argparse
import re
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")

cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))

from mappy_ext import revcomp, parse_cigar  # noqa: E402


def calculate_identity_from_cigar(cigar_string):
    """
    根据CIGAR字符串计算identity

    Parameters:
    cigar_string (str): CIGAR字符串, 如 "370=1I89=1D26=1I208=1D1250=1D2=2X1=1I373="

    Returns:
    float: identity值
    """
    pattern = r"(\d+)([=IDX])"
    matches = re.findall(pattern, cigar_string)

    match_count = 0
    total_aligned = 0

    for length_str, operation in matches:
        length = int(length_str)
        if operation == "=":
            match_count += length
            total_aligned += length
        elif operation == "X":
            total_aligned += length
        elif operation in ("I", "D"):
            total_aligned += length

    if total_aligned == 0:
        return 0.0
    return match_count / total_aligned


def read_fasta(fpath: str):
    """读取单序列FASTA文件"""
    headers = []
    seqs = {}
    with open(fpath, mode="r", encoding="utf8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                headers.append(line[1:])
                seqs[line[1:]] = ""
            elif headers:
                seqs[headers[-1]] += line
    return seqs


def extract_mismatch_positions(query, ref, hit):
    """
    提取比对中所有错配位点的位置信息。

    Returns:
    list of (query_pos_0based, target_pos_0based, query_base, target_base):
        query_pos_0based: query上的0-based位置
        target_pos_0based: target上的0-based位置
        query_base: query碱基
        target_base: target碱基
    """
    q_st = hit.q_st
    q_en = hit.q_en
    r_st = hit.r_st
    strand = hit.strand
    cigar_str = hit.cigar_str

    # 处理负链：将query替换为反向互补
    if strand == -1:
        query = revcomp(query)
        strand = 1

    q_pos = q_st
    r_pos = r_st
    mismatches = []

    for length, op in parse_cigar(cigar_str):
        if op in ("=", "X"):
            for _ in range(length):
                qb = query[q_pos]
                rb = ref[r_pos]
                if qb != rb:
                    mismatches.append((q_pos, r_pos, qb, rb))
                q_pos += 1
                r_pos += 1
        elif op == "I":
            q_pos += length
        elif op == "D":
            r_pos += length
        # 忽略 S, H, P 等

    return mismatches


def plot_mismatch_distribution(mismatches, target_len, output_prefix):
    """
    画不一致位点的分布图。

    Parameters:
    mismatches: list of (query_pos_0based, target_pos_0based, query_base, target_base)
    target_len: target序列长度
    output_prefix: 输出文件前缀
    """
    if not mismatches:
        print("No mismatches found in filtered alignments.")
        return

    # 收集target上的位置
    target_positions = [m[1] for m in mismatches]

    step = 5
    bins = list(range(0, target_len + step, step))

    fig, ax = plt.subplots(figsize=(14, 4))
    bin_counts = [sum(1 for p in target_positions if bins[i] <= p < bins[i + 1]) for i in range(len(bins) - 1)]

    centers = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]
    ax.bar(centers, bin_counts, width=step * 0.9, color="steelblue", edgecolor="white")
    ax.set_xlabel("Target Position (0-based)")
    ax.set_ylabel("Mismatch Count")
    ax.set_title(
        f"Mismatch Distribution (total: {len(mismatches)} mismatches)")
    ax.set_xticks(range(0, target_len + 1, 50))

    # 在顶部标出总不一致数
    ax.text(0.02, 0.95, f"Total mismatches: {len(mismatches)}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    out_path = f"{output_prefix}_mismatch_dist.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Mismatch distribution plot saved to: {out_path}")

    # 同时打印概览统计
    print(f"\nTarget length: {target_len}")
    print(f"Total mismatches: {len(mismatches)}")
    if target_len > 0:
        print(
            f"Mean mismatch position: {sum(target_positions)/len(target_positions):.0f}")
        print(
            f"Median mismatch position: {sorted(target_positions)[len(target_positions)//2]}")
        # 前10个和最后10个
        sorted_pos = sorted(target_positions)
        print(f"First mismatch at target position: {sorted_pos[0]}")
        print(f"Last mismatch at target position: {sorted_pos[-1]}")


def main():
    parser = argparse.ArgumentParser(
        description="Align query FASTA to target FASTA using mappy, "
                    "filter by coverage/identity, plot mismatch distribution."
    )
    parser.add_argument("--target", required=True,
                        help="Target FASTA file (single sequence)")
    parser.add_argument("--query", required=True,
                        help="Query FASTA file (one or more sequences)")
    parser.add_argument("--output-prefix", default="mismatch_dist",
                        help="Output file prefix (default: mismatch_dist)")
    parser.add_argument("--min-identity", type=float, default=0.99,
                        help="Minimum identity threshold (default: 0.99)")
    args = parser.parse_args()

    # 读取target (单序列)
    target_seqs = read_fasta(args.target)
    if len(target_seqs) != 1:
        print(
            f"Error: target FASTA should contain exactly 1 sequence, got {len(target_seqs)}")
        sys.exit(1)
    target_name = list(target_seqs.keys())[0]
    target_seq = target_seqs[target_name]
    print(f"Target: {target_name}, length: {len(target_seq)}")

    # 读取query
    query_seqs = read_fasta(args.query)
    print(f"Query sequences: {len(query_seqs)}")

    # 初始化aligner
    aligner = mappy.Aligner(seq=target_seq, extra_flags=67108864,
                            k=11, w=1, best_n=10, n_threads=1)

    # 收集所有符合条件的错配位点
    all_mismatches = []
    aligned_count = 0
    mismatched_cnt_in_aligned = 0

    for qname, qseq in tqdm(query_seqs.items(), desc="processing"):
        for hit in aligner.map(qseq):
            if not hit.is_primary:
                continue

            query_cov = (hit.q_en - hit.q_st) / len(qseq)
            target_cov = (hit.r_en - hit.r_st) / len(target_seq)
            identity = calculate_identity_from_cigar(hit.cigar_str)

            # 过滤条件
            if query_cov != 1.0 or target_cov != 1.0 or identity <= args.min_identity:
                continue

            aligned_count += 1
            # print(f"Aligned [{qname}]: Qcov={query_cov:.4f} Tcov={target_cov:.4f} ident={identity:.6f} "
            #       f"strand={hit.strand} [{hit.q_st}-{hit.q_en}]x[{hit.r_st}-{hit.r_en}]")

            # 提取错配位点
            mismatches = extract_mismatch_positions(qseq, target_seq, hit)
            if len(mismatches) > 0:
                mismatched_cnt_in_aligned += 1
            all_mismatches.extend(mismatches)

    print(
        f"\nTotal aligned (qcov=1, tcov=1, ident>{args.min_identity}): {aligned_count} / {len(query_seqs)} = {aligned_count / len(query_seqs)}")
    print(
        f"mismatched_ratio: {mismatched_cnt_in_aligned} / {aligned_count} = {mismatched_cnt_in_aligned / aligned_count}")
    print(f"Total mismatches collected: {len(all_mismatches)}")

    # 画图
    if all_mismatches:
        plot_mismatch_distribution(
            all_mismatches, len(target_seq), args.output_prefix)
    else:
        "GTAAAACGACGGCCAGTAGAGTTTGATCCTGGCTCAG"
        print("No mismatches to plot (either no alignments passed filter or all perfect matches).")


if __name__ == "__main__":
    main()
