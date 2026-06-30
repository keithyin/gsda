import argparse
import math
import pathlib
import sys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mappy

from tqdm import tqdm  # noqa: E402

cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))
sys.path.append(str(prev_prev_dir / "src"))

import gseda.align_ana.mappy_ext # noqa: E402
revcomp = gseda.align_ana.mappy_ext.revcomp # noqa: E402
parse_cigar = gseda.align_ana.mappy_ext.parse_cigar # noqa: E402
# from gseda.align_ana.mappy_ext import revcomp, parse_cigar  # noqa: E402


# ---------------------------------------------------------------------------
# FASTQ reader
# ---------------------------------------------------------------------------

def read_fastq(fpath: str):
    """读取 FASTQ 文件，返回 {name: (seq, qual)}"""
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
            records[name] = (seq, qual)
    return records


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


def extract_mismatch_with_baseq(query: str, qual: list[int], ref: str, hit) -> list[tuple[int, int, str, str, int]]:
    """
    提取比对中所有错配位点的位置 + baseq 信息。

    Returns:
        list of (query_pos_0based, target_pos_0based, query_base, target_base, baseq):
            baseq: 该位点在原始 reads 中的 PHRED 质量值
    """
    q_st = hit.q_st
    q_en = hit.q_en
    r_st = hit.r_st
    strand = hit.strand
    cigar_str = hit.cigar_str

    if strand == -1:
        query = revcomp(query)
        qual = qual[::-1]  # reverse qual too
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
    将 mismatch 位置分布 + 各分位点 baseq 画在同一 figure 的两个 subplot 中。

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

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 14), gridspec_kw={"height_ratios": [1, 1, 1]})

    # --- top: mismatch position distribution ---
    ax1.bar(centers, bin_counts, width=step * 0.9, color="steelblue", edgecolor="white")
    ax1.set_xlabel("Target Position (0-based)")
    ax1.set_ylabel("Mismatch Count")
    ax1.set_title(f"Mismatch Position Distribution  (total: {len(mismatches)} mismatches)")
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
        print(f"BaseQ percentiles: 5th={q05} 25th={q25} 50th={q50} 75th={q75} 95th={q95}")
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
    parser.add_argument("--ref", required=True, help="Reference FASTA file (single sequence)")
    parser.add_argument("--query", required=True, help="Query FASTQ file")
    parser.add_argument("--output-prefix", default="mismatch_dist", help="Output file prefix")
    parser.add_argument("--min-identity", type=float, default=0.99, help="Minimum identity (default: 0.99)")
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
        print(f"Error: reference should contain exactly 1 sequence, got {len(ref_seqs)}")
        sys.exit(1)
    ref_name = list(ref_seqs.keys())[0]
    ref_seq = ref_seqs[ref_name]
    print(f"Reference: {ref_name}, length: {len(ref_seq)}, seq={ref_seq[:100]}...")

    # Read query
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

    for _, (qseq, qual_str) in tqdm(query_records.items(), desc="processing"):
        qual = [ord(c) - 33 for c in qual_str]  # PHRED+33
        # print(qname)
        for hit in aligner.map(qseq):
            # print(f"name={qname}, seq={qseq[:100]}...")
            if not hit.is_primary:
                continue

            qcov = (hit.q_en - hit.q_st) / len(qseq)
            rcov = (hit.r_en - hit.r_st) / len(ref_seq)
            identity = calculate_identity_from_cigar(hit.cigar_str)

            if qcov != 1.0 or rcov != 1.0 or identity <= args.min_identity:
                filtered_by_cov_id += 1
                continue

            aligned_count += 1
            mismatches = extract_mismatch_with_baseq(qseq, qual, ref_seq, hit)
            if mismatches:
                mismatched_count += 1
            all_mismatches.extend(mismatches)

    print(f"\nAligned: {aligned_count} / {len(query_records)} (qcov=1, rcov=1, ident>{args.min_identity})")
    print(f"Mismatched reads: {mismatched_count} / {aligned_count}")
    print(f"Total mismatches: {len(all_mismatches)}")
    print(f"filtered_by_cov_id: {filtered_by_cov_id}")

    if all_mismatches:
        plot_mismatch_and_baseq(all_mismatches, len(ref_seq), args.output_prefix)
    else:
        print("No mismatches to plot.")


if __name__ == "__main__":
    main()
