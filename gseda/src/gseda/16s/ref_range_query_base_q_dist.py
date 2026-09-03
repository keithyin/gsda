"""参考区域 query baseQ 四类分布分析工具。

在 align_mismatch_dist_and_q_dist.py（只统计 mismatch 位点 baseQ）的基础上，
圈中 reference 上的一块区域 [start, end)，统计该区域内四类 CIGAR 操作对应的
query baseQ 分布：
  - eq(=)        匹配位点，取该列 query base 的 baseq；
  - mismatch(X)  错配位点，取该列 query base 的 baseq；
  - insertion(I) 插入的每个 query base，各取其 baseq；
  - deletion(D)  取 del 区域两端 query base 的最小值，重复 del 长度次。

处理流程：
1. 读取单序列参考 FASTA（ref）与 query（FASTQ 或 BAM，按 .bam 后缀自动选择）。
   读取 / readsq 口径与参考脚本一致（FASTQ 无 rq，用 readsq_from_baseq 计算）。
2. rq-thr 过滤低质量 read。
3. mappy 将每条 read 比对到参考（负链自动反向互补后统一为正链坐标）。
4. 按 qcov / rcov / identity 三重阈值过滤，只保留高质量比对。
5. 对每个保留的 hit，用 extract_baseq_by_op 走 CIGAR，把落在 [start, end) 内的
   四类操作分别收集 baseq（deletion 取两端 min，边界缺失则整条排除）。
6. plot_baseq_distributions 绘制 grouped bar chart：x 轴 = baseQ bin，每个 bin
   上 4 根柱子（eq / mismatch / insertion / deletion），每根柱 = 该 op 类型内该
   baseQ bin 的占比（各类内部各自归一，比例之和为 1），保证分布形状可比。
7. 打印概览：每类 count、baseQ 均值 / 取值范围。

输出：
- 一张 PNG 图（{output-prefix}_baseq_distributions.png）
- stdout 概览统计

用法示例：
    python ref_range_query_base_q_dist.py \
        --ref ref.fa --query reads.fq \
        --start 0 --end 1000 \
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


# ---------------------------------------------------------------------------
# FASTQ / BAM readers (与参考脚本一致)
# ---------------------------------------------------------------------------

def read_fastq(fpath: str):
    """读取 FASTQ，返回 {name: (seq, qual, rq)}，qual 为 PHRED int list。"""
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
    """读取 BAM，返回 {name: (seq, qual, rq)}，qual 为 PHRED int list。

    rq 从 BAM tag "rq"（错误概率）换算：-10 * log(1 - rq)。
    """
    import pysam
    records = {}
    with pysam.AlignmentFile(fpath, "rb", check_sq=False, threads=mp.cpu_count() // 2) as bam:
        for read in tqdm(bam.fetch(until_eof=True), desc=f"reading {fpath}"):
            records[read.query_name] = (
                read.query_sequence, list(read.query_qualities), -10 * math.log(1-read.get_tag("rq")))
    return records


def readsq_from_baseq(qual) -> float:
    """从 baseq 列表计算 readsq（PHRED）：对错误概率求均值再转回 PHRED。"""
    if not qual:
        return 0.0
    mean_err = sum(10.0 ** (-q / 10.0) for q in qual) / len(qual)
    return -10.0 * math.log10(mean_err)


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


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def _homo_run_from_start(seq: str, pos: int) -> int:
    """统计比对左边界 `pos` 处、朝 core(内)方向的同聚物(poly) run 长度。

    若边界"外"(pos-1)紧邻碱基与 seq[pos] 相同，说明边界切开了一个 poly，
    内侧即使不足 2 个也按 2 处理（保守多切一个，宁可多剔边界碱基，
    不留 poly 中段的噪声）。返回 0 表示该边界处不构成 poly。
    """
    base = seq[pos]
    n = 0
    while pos + n < len(seq) and seq[pos + n] == base:
        n += 1
    if n < 2 and pos >= 1 and seq[pos - 1] == base:
        n = 2
    return n


def _homo_run_from_end(seq: str, pos_end: int) -> int:
    """统计比对右边界 `pos_end`(exclusive)处、朝 core(内)方向的 poly run 长度。

    与 _homo_run_from_start 对称：从 pos_end-1 向左(朝内)数连续相同碱基，
    若不足 2 个但外侧紧邻碱基 seq[pos_end] 与 seq[pos_end-1] 相同，按 2 处理。
    返回 0 表示非 poly。
    """
    base = seq[pos_end - 1]
    n = 0
    while pos_end - 1 - n >= 0 and seq[pos_end - 1 - n] == base:
        n += 1
    if n < 2 and pos_end < len(seq) and seq[pos_end] == base:
        n = 2
    return n


def extract_baseq_by_op(query, qual, ref, hit,
                        start, end):
    """走 CIGAR，收集落在参考区域 [start, end) 内四类操作的 query baseQ。

    边界收缩（对所有 op 类型统一适用，不限 mismatch）：若比对的 target 边界
    （r_st / r_en）落在 ref 的同聚物(poly)区，则把两端 poly run 收缩进 core，
    只统计 core 区 [core_st, core_en) 内的操作；poly 长度变异视为噪声。
    有效分析窗口 = 用户区域 [start, end) ∩ core 区 [core_st, core_en)。

    返回 dict：
        {"=": [...], "X": [...], "I": [...], "D": [...]}
    其中每个 list 为对应 op 类型在有效窗口内的 baseq 列表：
      - "="  eq：该列 query base 的 baseq；
      - "X"  mismatch：该列 query base 的 baseq；
      - "I"  insertion：插入的每个 query base 各一个 baseq；
      - "D"  deletion：del 两端 query base 的最小值，重复裁剪后的 del 长度次。
            若 del 位于 read 边界导致某侧 query base 缺失，则整条该 deletion 排除。
    """
    q_st = hit.q_st
    q_en = hit.q_en
    r_st = hit.r_st
    r_en = hit.r_en
    strand = hit.strand
    cigar_str = hit.cigar_str

    if strand == -1:
        query = revcomp(query)
        qual = qual[::-1]
        q_st = len(query) - q_en

    # --- 边界收缩：target 边界落在 poly 区则收缩进 core（对所有 op 适用）---
    n_start = _homo_run_from_start(ref, r_st)
    n_end = _homo_run_from_end(ref, r_en)
    # 防止两端 poly run 收缩到重叠（read 极短或整段 poly 时）
    aln_span = r_en - r_st
    if n_start + n_end > aln_span:
        n_start = max(0, n_start - (n_start + n_end - aln_span))
        n_end = max(0, n_end - (aln_span - n_start))
    core_st = r_st + n_start
    core_en = r_en - n_end
    # 有效窗口 = 用户区域 ∩ core 区
    lo = max(start, core_st)
    hi = min(end, core_en)

    out = {"=": [], "X": [], "I": [], "D": []}
    q_pos = q_st
    r_pos = r_st

    for length, op in parse_cigar(cigar_str):
        if op == "=":
            for _ in range(length):
                if lo <= r_pos < hi:
                    out["="].append(qual[q_pos])
                q_pos += 1
                r_pos += 1
        elif op == "X":
            for _ in range(length):
                if lo <= r_pos < hi:
                    out["X"].append(qual[q_pos])
                q_pos += 1
                r_pos += 1
        elif op == "I":
            # 插入发生在 r_pos 处，ref 不动；窗口判据用 r_pos
            if lo <= r_pos < hi:
                for _ in range(length):
                    out["I"].append(qual[q_pos])
                    q_pos += 1
        elif op == "D":
            # del 占 ref [r_pos, r_pos+length)，裁剪到有效窗口
            k = min(r_pos + length, hi) - max(r_pos, lo)
            if k > 0:
                # del 不消耗 query，flanks 为 qual[q_pos-1] 与 qual[q_pos]
                if q_pos >= 1 and q_pos < len(qual):
                    val = min(qual[q_pos - 1], qual[q_pos])
                    out["D"].extend([val] * k)
            r_pos += length
        else:
            # M / N / S / H / P：参考脚本 CIGAR 只出 = X I D，此处兜底不处理
            pass

    return out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_baseq_distributions(all_bq: dict, output_prefix: str):
    """两联图：上 = grouped bar chart（x 轴固定 [0, 60)，step 1），
    下 = 各类 CDF 曲线（同一 x 轴）。

    上：每根柱高度 = 该 op 类型内该 baseQ bin 的占比（各类内部各自归一，sum=1）。
    下：每类一条 CDF 曲线，y = 累积占比，便于比较分布的左/右尾。
    超出 [0, 60) 的 baseQ 计入末 bin 59-60。
    """
    op_order = ["=", "X", "I", "D"]
    op_labels = ["eq", "mismatch", "insertion", "deletion"]
    colors = ["green", "steelblue", "orange", "purple"]

    all_vals = [v for op in op_order for v in all_bq.get(op, [])]
    if not all_vals:
        print("No baseQ values to plot.")
        return

    # 固定 x 轴 [0, 60)，step 1 → 60 个 bin
    bin_width = 1
    axis_max = 60
    n_bins = axis_max // bin_width
    bin_edges = list(range(0, axis_max + bin_width, bin_width))

    # 每类每 bin 的计数与占比
    counts = {}
    prop = {}
    for op in op_order:
        vals = all_bq.get(op, [])
        c = [0] * n_bins
        for v in vals:
            c[min(max(v, 0) // bin_width, n_bins - 1)] += 1
        total = sum(c)
        counts[op] = c
        prop[op] = [x / total if total else 0.0 for x in c]

    n = len(op_order)
    width = 0.8 / n
    bin_centers = [bin_edges[i] + bin_width / 2 for i in range(n_bins)]

    fig, (ax, ax_cdf) = plt.subplots(
        2, 1, figsize=(14, 9), sharex=True,
        gridspec_kw={"height_ratios": [1.2, 1.0]})

    # --- 上：grouped bar ---
    for oi, op in enumerate(op_order):
        xs = [c + (oi - (n - 1) / 2) * width for c in bin_centers]
        ax.bar(xs, prop[op], width=width, label=op_labels[oi],
               color=colors[oi], edgecolor="white", linewidth=0.5)
    ax.set_ylabel("Proportion")
    ax.set_title(
        "BaseQ Distribution by CIGAR Op Type (within selected reference region)")
    ax.set_xticks(range(0, axis_max, 5))
    ax.set_xlim(0, axis_max)
    ax.set_ylim(0, 1)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # --- 下：CDF（同一 [0, 60) 网格）---
    for op, label, color in zip(op_order, op_labels, colors):
        c = counts[op]
        total = sum(c)
        if total == 0:
            continue
        cum = []
        s = 0
        for x in c:
            s += x
            cum.append(s / total)
        ax_cdf.plot(bin_centers, cum, label=label,
                    color=color, linewidth=1.5)
    ax_cdf.set_xlabel("Base Quality (PHRED), step 1")
    ax_cdf.set_ylabel("Cumulative Proportion")
    ax_cdf.set_ylim(0, 1)
    ax_cdf.legend(loc="lower right", fontsize=10)
    ax_cdf.grid(axis="y", linestyle=":", alpha=0.4)

    plt.tight_layout()
    out_path = f"{output_prefix}_baseq_distributions.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Plot saved to: {out_path}")

    # --- 概览统计 ---
    print("\nPer-op baseQ summary:")
    for op, label in zip(op_order, op_labels):
        vals = all_bq.get(op, [])
        if vals:
            sv = sorted(vals)
            print(f"  {label:<10} count={len(vals):>8}  "
                  f"mean={sum(vals) / len(vals):.1f}  range={sv[0]}-{sv[-1]}")
        else:
            print(f"  {label:<10} count={0}")
    return out_path

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Align query to reference, and for a selected reference region "
                    "[start, end), collect query baseQ per CIGAR op type (eq/mismatch/insertion/deletion) "
                    "and plot within-type proportions."
    )
    parser.add_argument("--ref", required=True,
                        help="Reference FASTA file (single sequence)")
    parser.add_argument("--query", required=True,
                        help="Query FASTQ file, or BAM file (detected by .bam suffix)")
    parser.add_argument("--start", type=int, required=True,
                        help="Region start on reference (0-based, inclusive)")
    parser.add_argument("--end", type=int, required=True,
                        help="Region end on reference (0-based, exclusive)")
    parser.add_argument("--output-prefix",
                        default=None,
                        help="Output file prefix (default: query path without extension)")
    parser.add_argument("--min-identity", type=float,
                        default=0.95, help="Minimum identity (default: 0.95)")
    parser.add_argument("--min-qcov", type=float,
                        default=0.5, help="Minimum query coverage (default: 0.5)")
    parser.add_argument("--min-rcov", type=float,
                        default=0.5, help="Minimum reference coverage (default: 0.5)")
    parser.add_argument("--rq-thr", type=float,
                        default=20.0, help="Minimum read quality (PHRED); reads below this are "
                        "filtered out (default: 20.0, no filter)")
    parser.add_argument("--context", type=int,
                        default=100, help="Number of surrounding bases to print around the "
                        "region (default: 100)")
    args = parser.parse_args()

    if args.output_prefix is None:
        args.output_prefix = pathlib.Path(
            args.query).with_suffix("").as_posix()

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
    # print(f"Reference: {ref_name}, length: {len(ref_seq)}, seq={ref_seq[:100]}...")

    # 校验区域
    if not (0 <= args.start < args.end <= len(ref_seq)):
        print(f"Error: region invalid. require 0 <= start < end <= {len(ref_seq)}, "
              f"got start={args.start}, end={args.end}")
        sys.exit(1)
    print(
        f"Selected region: [{args.start}, {args.end}) length {args.end - args.start}")

    # 打印感兴趣区域及其周边碱基（region 高亮，周边各 context 个碱基）
    ctx = args.context
    ctx_st = max(0, args.start - ctx)
    ctx_en = min(len(ref_seq), args.end + ctx)
    def _show(a, b): return ref_seq[a:b]
    print(
        f"Region (with {ctx} surrounding bases each side) [{ctx_st}, {ctx_en}):")
    if ctx_st > 0:
        print(f"  ...{_show(ctx_st, args.start)}", end="")
    print(f"[{_show(args.start, args.end)}]", end="")
    if ctx_en < len(ref_seq):
        print(f"{_show(args.end, ctx_en)}...")

    # Read query（按后缀自动选择 FASTQ / BAM）
    if args.query.endswith(".bam"):
        query_records = read_bam(args.query)
    else:
        query_records = read_fastq(args.query)
    print(f"Query records: {len(query_records)}")

    aligner = mappy.Aligner(seq=ref_seq, extra_flags=67108864,
                            k=11, w=1, best_n=10, n_threads=1)

    all_bq = {"=": [], "X": [], "I": [], "D": []}
    # 记录每类 op 命中过的 query name（同一 query 可同时出现于多类）
    op_query_names = {"X": [], "I": [], "D": []}
    aligned_count = 0
    filtered_by_cov_id = 0
    filtered_by_rq = 0

    for qname, (qseq, qual, rq) in tqdm(query_records.items(), desc="processing"):
        if rq is None:
            rq = readsq_from_baseq(qual)
        if rq < args.rq_thr:
            filtered_by_rq += 1
            continue

        for hit in aligner.map(qseq):
            if not hit.is_primary:
                continue

            qcov = (hit.q_en - hit.q_st) / len(qseq)
            rcov = (hit.r_en - hit.r_st) / len(ref_seq)
            identity = calculate_identity_from_cigar(hit.cigar_str)

            if qcov < args.min_qcov or rcov < args.min_rcov or identity <= args.min_identity:
                filtered_by_cov_id += 1
                continue

            aligned_count += 1
            res = extract_baseq_by_op(
                qseq, qual, ref_seq, hit, args.start, args.end)
            for op in all_bq:
                all_bq[op].extend(res[op])
            for op in ("X", "I", "D"):
                if res[op]:
                    op_query_names[op].append(qname)

    print(f"\nAligned: {aligned_count} / {len(query_records)} "
          f"(qcov>={args.min_qcov}, rcov>={args.min_rcov}, ident>{args.min_identity})")
    print(f"filtered_by_cov_id: {filtered_by_cov_id}")
    print(f"filtered_by_rq: {filtered_by_rq}")

    # 打印各类 op 命中的 query name（同一 query 可同时出现在多类）
    op_labels = {"X": "mismatch", "I": "insertion", "D": "deletion"}
    print("\nQuery names per op type (in selected region):")
    for op in ("X", "I", "D"):
        names = list(dict.fromkeys(op_query_names[op]))  # 去重保序
        print(f"  {op_labels[op]} ({len(names)}):")
        for nm in names:
            print(f"    {nm}")
        print("    {}".format(",".join([nm.split("/")[1] for nm in names])))

    if any(all_bq[op] for op in all_bq):
        import shutil
        image_path = plot_baseq_distributions(all_bq, args.output_prefix)
        gsda_tmp_path = pathlib.Path("/root/projects/gsda/tmp-data-dir")
        # if gsda_tmp_path.exists():

        shutil.copy(image_path, gsda_tmp_path /
                    pathlib.Path(image_path).name)

    else:
        print("No baseQ values in the selected region.")


if __name__ == "__main__":
    main()
