import pysam
import mappy
import argparse
import math
import matplotlib.pyplot as plt
import pathlib
from typing import Dict, List, Tuple
import multiprocessing as mp
from tqdm import tqdm


def calculate_identity_from_cigar(cigar_string: str) -> float:
    """
    Calculate identity from a CIGAR string.
    Matches (=) contribute to match_count and total_aligned.
    Mismatches (X) contribute to total_aligned.
    Indels (I, D) contribute to total_aligned.
    """
    import re
    pattern = r'(\d+)([=IDX])'

    matches = re.findall(pattern, cigar_string)

    match_count = 0
    total_aligned = 0
    del_cnt = 0
    ins_cnt = 0
    mismatch_cnt = 0

    for length_str, operation in matches:
        length = int(length_str)
        if operation == '=':
            match_count += length
            total_aligned += length
        elif operation == 'X':
            mismatch_cnt += length
            total_aligned += length
        elif operation == 'I':
            ins_cnt += length
            total_aligned += length
        elif operation == 'D':
            del_cnt += length
            total_aligned += length

    if total_aligned == 0:
        return 0.0
    return match_count / total_aligned, mismatch_cnt, ins_cnt, del_cnt


def load_bam_reads(bam_path: str) -> Dict[int, List[str]]:
    """
    Load reads from a BAM file and index them by the 'ch' tag.
    Returns a dictionary where keys are ch values and values are lists of sequences.
    """
    reads_map = {}
    with pysam.AlignmentFile(bam_path, "rb", threads=mp.cpu_count()) as bam:
        for record in tqdm(bam.fetch(until_eof=True), desc=f"reading {bam_path}"):
            try:
                ch = int(record.get_tag("ch"))
                seq = record.query_sequence
                if ch not in reads_map:
                    reads_map[ch] = []
                reads_map[ch].append(seq)
            except (KeyError, ValueError):
                continue
    return reads_map


def calculate_phred(identity: float) -> float:
    """
    Convert identity to Phred quality score.
    Q = -10 * log10(1 - identity)
    """
    if identity >= 1.0:
        return 60.0  # Cap at 60 for perfect matches
    if identity <= 0.0:
        return 0.0
    return -10 * math.log10(1.0 - identity)


def plot_distribution(data: List[float], title: str, xlabel: str, filename: str):
    """
    Plot a histogram of the given data.
    """
    plt.figure(figsize=(8, 6))
    plt.hist(data, bins=50, edgecolor='black', alpha=0.7)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(filename)
    plt.close()
    print(f"Plot saved to {filename}")


def align_sequences(seq_true: str, seq_false: str) -> Tuple[float, str]:
    """
    Align two sequences using mappy and return identity and CIGAR.
    """
    aligner = mappy.Aligner(seq=seq_true, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    hits = list(aligner.map(seq_false))
    if not hits:
        return 0.0, 0, 0, 0

    best_hit = hits[0]
    identity, mismatch, ins_cnt, del_cnt = calculate_identity_from_cigar(
        best_hit.cigar_str)
    return identity, mismatch, ins_cnt, del_cnt


def process_molecule(ch: int, seqs_true: List[str]) -> Tuple[int, float, int, int, int, float]:
    """
    Worker function to align two sequences for a given molecule and calculate metrics.
    Returns: (ch, identity, mismatch, ins_cnt, del_cnt, phreq)
    """
    identity, mismatch, ins_cnt, del_cnt = align_sequences(
        seq_true=seqs_true[0], seq_false=seqs_true[1])
    phreq = calculate_phred(identity=identity)
    return ch, identity, mismatch, ins_cnt, del_cnt, phreq


def _wrapper(args):
    return process_molecule(*args)


def main():
    """
        分析bystrand=True时, 同一 channel 正反链的结果差异
    """
    parser = argparse.ArgumentParser(
        description="Analyze differences between Combined and Separate SMC consensus modes.")

    parser.add_argument("--bystrand-bam", type=str, required=True,
                        help="BAM file from Separate Mode (byStrand=True)")

    args = parser.parse_args()

    out_prefix = str(pathlib.Path(args.bystrand_bam).with_suffix(''))

    print(f"Loading bystrand  BAM: {args.bystrand_bam}...")
    true_reads = load_bam_reads(args.bystrand_bam)

    # In combined mode, we expect one read per ch.
    # In separate mode, we expect potentially multiple (fwd/rev).
    # We'll simplify true_reads to just one sequence per ch for the purpose of this analysis.

    true_reads_double = {ch: seqs
                         for ch, seqs in true_reads.items() if len(seqs) == 2}

    common_chs = true_reads_double.keys()
    print(f"Found {len(common_chs)} channels (ch).")

    mismatches = []
    insertions = []
    deletions = []
    phred_scores = []
    diff_counts = {i: 0 for i in range(11)}

    fwd_rev_eq = 0

    tasks = [(ch, true_reads_double[ch]) for ch in sorted(common_chs)]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = []
        for r in tqdm(
            pool.imap_unordered(_wrapper, tasks, chunksize=8),
            total=len(tasks),
            desc="processing ..."
        ):
            results.append(r)

    for res in results:
        ch, identity, mismatch, ins_cnt, del_cnt, phreq = res

        tot_diff_cnt = mismatch + ins_cnt + del_cnt
        diff_counts[min(tot_diff_cnt, 10)] += 1

        mismatches.append(mismatch)
        insertions.append(ins_cnt)
        deletions.append(del_cnt)
        phred_scores.append(phreq)

        if identity > 0.9999:
            fwd_rev_eq += 1

    fwd_rev_eq_ratio = fwd_rev_eq / len(common_chs)
    print(f"fwd_rev_eq_ratio:{fwd_rev_eq_ratio:.2f}")

    print("\nTotal Difference Distribution:")
    diff_dist = {}
    for i in range(11):
        ratio = diff_counts[i] / len(common_chs)
        label = f"{i}" if i < 10 else ">=10"
        diff_dist[label] = ratio
        print(f"{label}: {ratio:.4%}")
    print("")

    plot_distribution(mismatches, "Mismatch Count Distribution",
                      "Mismatches", f"{out_prefix}_mismatch_dist.png")
    plot_distribution(insertions, "Insertion Count Distribution",
                      "Insertions", f"{out_prefix}_ins_dist.png")
    plot_distribution(deletions, "Deletion Count Distribution",
                      "Deletions", f"{out_prefix}_del_dist.png")
    plot_distribution(phred_scores, "Phred Quality Score Distribution",
                      "Phred Score", f"{out_prefix}_phred_dist.png")


if __name__ == "__main__":
    main()
