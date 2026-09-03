"""比对错配位置分布 + 各分位点 baseQ 分析工具。

将一组 query reads 比到单条参考序列上，收集所有错配(mismatch)位点，
统计其沿参考序列的位置分布与各分位点的 PHRED baseQ 值，并输出可视化图与
概览统计，用于评估比对在参考序列不同区域的"错配集中度"和"错误置信度"。

处理流程：
1. 读取单序列参考 FASTA（ref）与 query（FASTQ 或 BAM，按 .bam 后缀自动选择读取方式）。
   - FASTQ：read_fastq 逐条读入，qual 解析为 PHRED int 列表。
   - BAM：read_bam 用 pysam 读入，跳过未比对 reads(is_unmapped)，query_sequence
     已含比对方向（负链为反向互补序列），与 FASTQ 路径统一为原始 reads 形式；
     rq 从 BAM tag "rq"（错误概率）换算为 -10*log(1-rq)。
2. 按 read quality(readsq, PHRED) 过滤：rq-thr>0 时丢弃低于阈值的 read。
   FASTQ 无 rq 字段，用 readsq_from_baseq 从 baseq 计算（按错误概率求均值再转
   回 PHRED，而非直接对 baseq 求算术均值）。
3. 用 mappy 将每条 read 比对到参考（负链自动反向互补后处理，统一为正链坐标）。
4. 按 query 覆盖度(qcov)、参考覆盖度(rcov)、序列一致性(identity)三重阈值
   过滤，只保留高质量比对（identity 用 CIGAR 中 = 占比计算）。
5. 对每个保留的 hit，用 extract_mismatch_with_baseq 提取错配位点：
   - 解析 CIGAR，逐列比对 query/ref 碱基，取 baseq（原始 reads 的 PHRED 质量）；
   - 边界收缩：若比对两端落在参考的同聚物(poly)区，把两端 poly run 内的错配
     忽略（poly 长度变异视为噪声），只统计稳定 core 区的错配。
6. plot_mismatch_and_baseq 绘制三联图（同一 figure 的三个 subplot）：
   - 上：错配位置分布直方图（按 5bp 分箱）；
   - 中：各分位点(5/25/50/75/95) baseQ 沿位置变化 + 均值曲线；
   - 下：各分箱 baseQ 的标准差曲线。
7. 打印概览统计：参考长度、错配总数、baseQ 各分位点/均值/取值范围。

输出：
- 一张 PNG 图（{output-prefix}_mismatch_and_baseq.png）
- stdout 概览统计

用法示例：
    python align_mismatch_dist_and_q_dist.py \
        --ref ref.fa --query reads.fq \
        --output-prefix out --min-identity 0.99
"""

import mappy
import matplotlib.pyplot as plt
import argparse
import math
import pathlib
import sys
import os
import multiprocessing as mp

import matplotlib
matplotlib.use("Agg")

from tqdm import tqdm  # noqa: E402

cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))
sys.path.append(str(prev_prev_dir / "src"))

import gseda.align_ana.mappy_ext  # noqa: E402
revcomp = gseda.align_ana.mappy_ext.revcomp  # noqa: E402
parse_cigar = gseda.align_ana.mappy_ext.parse_cigar  # noqa: E402
# from gseda.align_ana.mappy_ext import revcomp, parse_cigar  # noqa: E402


# ---------------------------------------------------------------------------
# FASTQ reader
# ---------------------------------------------------------------------------

def read_fastq(fpath: str):
    """读取 FASTQ 文件，返回 {name: (seq, qual, rq)}，qual 为 PHRED int list。

    FASTQ 无 readsq 字段，rq 置 None，由调用方从 baseq 计算。
    """
    records = {}
    with open(fpath, "r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            name = header[1:].split()[0]
            records[name] = (seq, [ord(c) - 33 for c in qual], None)
    return records


def read_bam(fpath: str):
    """读取 BAM 文件，返回 {name: (seq, qual, rq)}，qual 为 PHRED int list。

    只保留 mapped 的 reads（is_unmapped 跳过）。query_sequence 已含比对方向
    （负链为反向互补序列），与 FASTQ 路径一致，下游按原始 reads 处理。
    rq 从 BAM tag "rq"（错误概率）换算：-10 * log(1 - rq)，用自然对数
    math.log（注意与 readsq_from_baseq 的 log10 口径不一致）。
    """
    import pysam
    records = {}
    with pysam.AlignmentFile(fpath, "rb", check_sq=False, threads=mp.cpu_count() // 2) as bam:
        for read in tqdm(bam.fetch(until_eof=True), desc=f"reading {fpath}"):

            records[read.query_name] = (
                read.query_sequence, list(read.query_qualities), -10 * math.log(1-read.get_tag("rq")))
    return records


def readsq_from_baseq(qual: list[int]) -> float:
    """从 baseq 列表计算 readsq（PHRED）。

    不能直接对 baseq 求算术均值：PHRED 质量是 log10 尺度，正确做法是先把每个
    碱基转成错误概率 10^(-q/10)、对概率求均值、再转回 PHRED，即
    -10 * log10(mean(10^(-q/10)))。
    """
    if not qual:
        return 0.0
    mean_err = sum(10.0 ** (-q / 10.0) for q in qual) / len(qual)
    return -10.0 * math.log10(mean_err)


# ---------------------------------------------------------------------------
# Alignment helpers (mirrors align_mismatch_dist.py)
# ---------------------------------------------------------------------------

def calculate_identity_from_cigar(cigar_string: str) -> float:
    import re
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
    return match_count / total_aligned if total_aligned else 0.0


def _homo_run_from_start(seq: str, pos: int) -> int:
    """统计比对左边界 `pos` 处、朝 core(内)方向的同聚物(poly) run 长度。

    用途(见 extract_mismatch_with_baseq)：比对两端若落在 poly 区(如 "AAAAA")，
    poly 长度变异被视为测序/比对噪声，需要从错配统计中剔除。本函数量出从
    边界向内连续相同碱基的长度 n，下游据此把 [r_st, r_st+n) 收缩掉。

    设计要点：
    - run 只朝"内"(core 方向，即 pos 向右)数；边界"外"(pos-1 及更左)的碱基
      不在比对覆盖范围内，不会被计入。
    - 问题：poly 可能被比对边界从中间切开，poly 主体在外侧，向内只剩 1 个碱基
      （例如 ... A A A A | A G，边界切在倒数第 2 个 A，向内只数到 1 个 A）。
      此时若只看内侧 n=1，会误判为"非 poly"而 return 0，导致这段被切开
      的 poly 没被切除，其中残留的假错配会混入统计。
    - 对策：若内侧不足 2 个，但外侧紧邻碱基 seq[pos-1] 与 seq[pos] 相同，
      即可确认此处并非孤立碱基，而是一个被边界截断的 poly —— 返回 2。

    为何返回 2 而非真实总长：
    - 往外侧延伸多少个碱基本函数无从得知(可能远超 2)，但真实长度无关紧要。
    - 该值有两重下游用途：(a) 二值判据，n>0 表示"落在 poly 区"；
      (b) 用于 new_r_st = r_st + n 实际移动边界切除。用 2 表示"当前碱基 +
      外侧那一个同类碱基"都属 poly 边界区，比返回 1 更保守(多切一个，
      宁可多剔一个 poly 边界碱基，也不留 poly 中段的假错配)，同时不会
      因切太多而吃掉真正的 core 错配。

    返回 0 表示该边界处不构成 poly(相同碱基不足 2 个)。
    """
    base = seq[pos]
    n = 0
    while pos + n < len(seq) and seq[pos + n] == base:
        n += 1
    # 内侧不足 2 但外侧紧邻同类碱基 → 边界切开了一个 poly，至少 2 个，按 poly 处理
    if n < 2 and pos >= 1 and seq[pos - 1] == base:
        n = 2
    return n


def _homo_run_from_end(seq: str, pos_end: int) -> int:
    """统计比对右边界 `pos_end`(exclusive)处、朝 core(内)方向的 poly run 长度。

    与 _homo_run_from_start 完全对称：从 pos_end-1 向左(朝内)数连续相同碱基，
    若不足 2 个但外侧紧邻碱基 seq[pos_end] 与 seq[pos_end-1] 相同，说明边界
    从中间切开了一个 poly(主体在外侧)，返回 2。

    下游用途：返回的 n 用于 new_r_en = r_en - n 切除右端 poly 区，
    同时 n>0 作为"该边界落在 poly 区"的二值判据。返回 0 表示非 poly。
    """
    base = seq[pos_end - 1]
    n = 0
    while pos_end - 1 - n >= 0 and seq[pos_end - 1 - n] == base:
        n += 1
    # 内侧不足 2 但外侧紧邻同类碱基 → 边界切开了一个 poly，至少 2 个，按 poly 处理
    if n < 2 and pos_end < len(seq) and seq[pos_end] == base:
        n = 2
    return n


def extract_mismatch_with_baseq(query: str, qual: list[int], ref: str, hit) -> list[tuple[int, int, str, str, int]]:
    """
    提取比对中所有错配位点的位置 + baseq 信息。

    边界收缩：若比对边界落在 ref 的同聚物(poly)区，则把落在两端 poly run
    内的错配忽略掉（poly 长度变异视为噪声，仅统计稳定 core 区的错配）。
    这里只"过滤"错配，不移动 q/r 比对边界，因此 CIGAR walk 全程保持 1:1 对齐。

    Returns:
        list of (query_pos_0based, target_pos_0based, query_base, target_base, baseq):
            baseq: 该位点在原始 reads 中的 PHRED 质量值
    """
    q_st = hit.q_st
    q_en = hit.q_en
    r_st = hit.r_st
    r_en = hit.r_en
    strand = hit.strand
    cigar_str = hit.cigar_str

    if strand == -1:
        query = revcomp(query)
        qual = qual[::-1]  # reverse qual too
        strand = 1
        # 负链：比对在 rc(query) 上，第一列对齐的位置是 len-q_en
        q_st = len(query) - q_en

    # --- 边界收缩：跳过两端落在 poly 区的一部分 ---------------------------
    # 基于 ref 坐标计算，与 strand 无关。
    n_start = _homo_run_from_start(ref, r_st)
    n_end = _homo_run_from_end(ref, r_en)
    # 防止两端 poly run 收缩到重叠（read 极短或整段 poly 时）
    aln_span = r_en - r_st
    if n_start + n_end > aln_span:
        n_start = max(0, n_start - (n_start + n_end - aln_span))
        n_end = max(0, n_end - (aln_span - n_start))
    new_r_st = r_st + n_start
    new_r_en = r_en - n_end

    q_pos = q_st
    r_pos = r_st
    mismatches = []

    for length, op in parse_cigar(cigar_str):
        if op in ("=", "X"):
            for _ in range(length):
                qb = query[q_pos]
                rb = ref[r_pos]
                # 仅保留 core 区错配；两端 poly 边界区一律丢弃
                if qb != rb and new_r_st <= r_pos < new_r_en:
                    mismatches.append((q_pos, r_pos, qb, rb, qual[q_pos]))
                q_pos += 1
                r_pos += 1
        elif op == "I":
            q_pos += length
        elif op == "D":
            r_pos += length

    return mismatches


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_mismatch_and_baseq(mismatches, target_len, output_prefix: str):
    """
    将 mismatch 位置分布 + 各分位点 baseq 画在同一 figure 的三个 subplot 中
    （上：位置分布直方图；中：各分位点 + 均值；下：各分箱标准差）。

    Parameters:
        mismatches: list of (query_pos, target_pos, q_base, ref_base, baseq)
        target_len: reference 长度
        output_prefix: 输出前缀
    """
    target_positions = [m[1] for m in mismatches]
    baseqs = [m[4] for m in mismatches]

    step = 5
    bins = list(range(0, target_len + step, step))
    bin_counts = [
        sum(1 for p in target_positions if bins[i] <= p < bins[i + 1])
        for i in range(len(bins) - 1)
    ]
    centers = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(
        14, 14), gridspec_kw={"height_ratios": [1, 1, 1]})

    # --- top: mismatch position distribution ---
    ax1.bar(centers, bin_counts, width=step * 0.9,
            color="steelblue", edgecolor="white")
    ax1.set_xlabel("Target Position (0-based)")
    ax1.set_ylabel("Mismatch Count")
    ax1.set_title(
        f"Mismatch Position Distribution  (total: {len(mismatches)} mismatches)")
    ax1.set_xticks(range(0, target_len + 1, 50))
    ax1.text(0.02, 0.95, f"Total mismatches: {len(mismatches)}",
             transform=ax1.transAxes, va="top", ha="left", fontsize=10,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # --- bottom: baseq percentile lines (per-bin, varies with position) ---
    # Bin baseQ values by position so the lines show variation along the target
    bin_baseqs: list[list[int]] = [[] for _ in range(len(bins) - 1)]
    for pos, bq in zip(target_positions, baseqs):
        for i in range(len(bins) - 1):
            if bins[i] <= pos < bins[i + 1]:
                bin_baseqs[i].append(bq)
                break

    quantile_fracs = [0.05, 0.25, 0.50, 0.75, 0.95]
    q_labels = ["5th", "25th", "50th", "75th", "95th"]
    colors = ["purple", "darkorange", "green", "blue", "red"]
    for qi, qval in enumerate(quantile_fracs):
        y_vals = []
        for bin_vals in bin_baseqs:
            if len(bin_vals) >= 2:
                sorted_bin = sorted(bin_vals)
                idx = int(len(sorted_bin) * qval)
                y_vals.append(sorted_bin[min(idx, len(sorted_bin) - 1)])
            else:
                y_vals.append(None)
        # plot with None as NaN so the line breaks at empty bins
        ax2.plot(centers, [float("nan") if v is None else v for v in y_vals],
                 color=colors[qi], linestyle="--", linewidth=1.2,
                 alpha=0.8, label=f"{q_labels[qi]}")

    mean_vals = [sum(b) / len(b) if b else float("nan") for b in bin_baseqs]
    ax2.plot(centers, mean_vals, color="darkred", linestyle="-",
             linewidth=1.5, alpha=0.9, label="mean")

    ax2.set_xlabel("Mismatch Position (0-based)")
    ax2.set_ylabel("Base Quality (PHRED)")
    ax2.set_title("BaseQ Distribution at Mismatch Positions")
    ax2.set_xticks(range(0, target_len + 1, 50))
    ax2.legend(loc="upper right", fontsize=9)

    # --- third: std per bin ---
    def _std(vals):
        if len(vals) < 2:
            return float("nan")
        m = sum(vals) / len(vals)
        return math.sqrt(sum((x - m) ** 2 for x in vals) / len(vals))
    std_vals = [_std(b) for b in bin_baseqs]
    ax3.plot(centers, std_vals, color="gray", linestyle=":", linewidth=1.5,
             label="std")
    ax3.set_xlabel("Mismatch Position (0-based)")
    ax3.set_ylabel("Std Dev")
    ax3.set_title("BaseQ Std Dev at Mismatch Positions")
    ax3.set_xticks(range(0, target_len + 1, 50))
    ax3.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    out_path = f"{output_prefix}_mismatch_and_baseq.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Plot saved to: {out_path}")

    # --- 概览统计 ---
    print(f"\nTarget length: {target_len}")
    print(f"Total mismatches: {len(mismatches)}")
    if baseqs:
        sorted_bq = sorted(baseqs)
        n = len(sorted_bq)
        q05 = sorted_bq[int(n * 0.05)]
        q25 = sorted_bq[int(n * 0.25)]
        q50 = sorted_bq[int(n * 0.50)]
        q75 = sorted_bq[int(n * 0.75)]
        q95 = sorted_bq[int(n * 0.95)]
        mn = sum(baseqs) / n
        print(
            f"BaseQ percentiles: 5th={q05} 25th={q25} 50th={q50} 75th={q75} 95th={q95}")
        print(f"BaseQ mean: {mn:.1f}")
        print(f"BaseQ range: {sorted_bq[0]} – {sorted_bq[-1]}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Align query FASTQ to reference FASTA, collect mismatch positions + baseQ, "
                    "and plot the distribution."
    )
    parser.add_argument("--ref", required=True,
                        help="Reference FASTA file (single sequence)")
    parser.add_argument("--query", required=True,
                        help="Query FASTQ file, or BAM file (detected by .bam suffix)")
    parser.add_argument("--output-prefix",
                        default="mismatch_dist", help="Output file prefix")
    parser.add_argument("--min-identity", type=float,
                        default=0.99, help="Minimum identity (default: 0.99)")
    parser.add_argument("--min-qcov", type=float,
                        default=1.0, help="Minimum query coverage (default: 1.0)")
    parser.add_argument("--min-rcov", type=float,
                        default=1.0, help="Minimum reference coverage (default: 1.0)")
    parser.add_argument("--rq-thr", type=float,
                        default=0.0, help="Minimum read quality (PHRED); reads below this are "
                        "filtered out (default: 0.0, no filter)")
    args = parser.parse_args()

    # Read reference
    ref_seqs = {}
    with open(args.ref, "r") as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].split()[0]
                ref_seqs[header] = ""
            elif header:
                ref_seqs[header] += line
    if len(ref_seqs) != 1:
        print(
            f"Error: reference should contain exactly 1 sequence, got {len(ref_seqs)}")
        sys.exit(1)
    ref_name = list(ref_seqs.keys())[0]
    ref_seq = ref_seqs[ref_name]
    print(
        f"Reference: {ref_name}, length: {len(ref_seq)}, seq={ref_seq[:100]}...")

    # Read query（按后缀自动选择 FASTQ / BAM 读取方式，统一返回 {name: (seq, qual, rq)}；
    # FASTQ 的 rq 为 None，在比对循环中从 baseq 计算）
    if args.query.endswith(".bam"):
        query_records = read_bam(args.query)
    else:
        query_records = read_fastq(args.query)
    print(f"Query records: {len(query_records)}")

    # Align
    # aligner = mappy.Aligner(ref_seq, extra_flags=67108864, k=11, w=1, best_n=10, n_threads=1)
    aligner = mappy.Aligner(seq=ref_seq, extra_flags=67108864,
                            k=11, w=1, best_n=10, n_threads=1)

    all_mismatches = []
    aligned_count = 0
    mismatched_count = 0

    filtered_by_cov_id = 0
    filtered_by_rq = 0

    for _, (qseq, qual, rq) in tqdm(query_records.items(), desc="processing"):
        # 低质量 read 过滤：BAM 直接读 rq 字段；FASTQ 无 rq 字段，从 baseq 计算
        if rq is None:
            rq = readsq_from_baseq(qual)
        if rq < args.rq_thr:
            filtered_by_rq += 1
            continue

        for hit in aligner.map(qseq):
            # print(f"name={qname}, seq={qseq[:100]}...")
            if not hit.is_primary:
                continue

            qcov = (hit.q_en - hit.q_st) / len(qseq)
            rcov = (hit.r_en - hit.r_st) / len(ref_seq)
            identity = calculate_identity_from_cigar(hit.cigar_str)

            if qcov < args.min_qcov or rcov < args.min_rcov or identity <= args.min_identity:
                filtered_by_cov_id += 1
                continue

            aligned_count += 1
            mismatches = extract_mismatch_with_baseq(qseq, qual, ref_seq, hit)
            if mismatches:
                mismatched_count += 1
            all_mismatches.extend(mismatches)

    print(
        f"\nAligned: {aligned_count} / {len(query_records)} "
        f"(qcov>={args.min_qcov}, rcov>={args.min_rcov}, ident>{args.min_identity})")
    print(f"Mismatched reads: {mismatched_count} / {aligned_count}")
    print(f"Total mismatches: {len(all_mismatches)}")
    print(f"filtered_by_cov_id: {filtered_by_cov_id}")
    print(f"filtered_by_rq: {filtered_by_rq}")

    if all_mismatches:
        plot_mismatch_and_baseq(
            all_mismatches, len(ref_seq), args.output_prefix)
    else:
        print("No mismatches to plot.")


if __name__ == "__main__":
    main()
