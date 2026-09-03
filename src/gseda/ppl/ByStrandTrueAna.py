import pysam
import mappy
import polars as pl
import argparse
import math
import matplotlib.pyplot as plt
import pathlib
from typing import Dict, List, Tuple
import multiprocessing as mp
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import sys
cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
print(cur_dir)
sys.path.append(cur_dir)
import env_prepare  # noqa: E402


class FwdRev:
    def __init__(self, fwd=None, rev=None):
        self.fwd = fwd
        self.rev = rev

    def has_rev_fwd(self):
        return self.fwd is not None and self.rev is not None

    def __repr__(self):
        return f"FwdRev(\n fwd={self.fwd}, \n rev={self.rev})"


def stat_alignment_info(hit, ref: str, query: str):
    from collections import defaultdict

    cigar = hit.cigar  # [(len, op)]

    # 统计
    match_count = 0
    mismatch_cnt = 0
    ins_cnt = 0
    del_cnt = 0
    total_aligned = 0

    mismatch_detail = defaultdict(int)
    ins_detail = defaultdict(int)
    del_detail = defaultdict(int)

    ref_idx = hit.r_st
    query_idx = hit.q_st

    for length, op in cigar:

        # ===== match =====
        if op == 7:  # '='
            match_count += length
            total_aligned += length
            ref_idx += length
            query_idx += length

        # ===== mismatch =====
        elif op == 8:  # 'X'
            mismatch_cnt += length
            total_aligned += length

            for i in range(length):
                r_base = ref[ref_idx + i]
                q_base = query[query_idx + i]
                mismatch_detail[(r_base, q_base)] += 1

            ref_idx += length
            query_idx += length

        # ===== insertion =====
        elif op == 1:  # 'I'
            ins_cnt += length
            total_aligned += length

            for i in range(length):
                q_base = query[query_idx + i]
                ins_detail[q_base] += 1

            query_idx += length

        # ===== deletion =====
        elif op == 2:  # 'D'
            del_cnt += length
            total_aligned += length

            for i in range(length):
                r_base = ref[ref_idx + i]
                del_detail[r_base] += 1

            ref_idx += length

        # ===== fallback: M（未区分）=====
        elif op == 0:
            raise RuntimeError("eqx needed")

        else:
            # 其他操作（clip等）一般不会出现在 mappy 的 hit.cigar
            pass

    identity = match_count / total_aligned if total_aligned else 0.0

    return {
        "identity": identity,
        "match": match_count,
        "mismatch": mismatch_cnt,
        "ins": ins_cnt,
        "del": del_cnt,
        "mismatch_detail": dict(mismatch_detail),
        "ins_detail": dict(ins_detail),
        "del_detail": dict(del_detail),
    }


def load_bam_reads(bam_path: str, np_thr: int = 5) -> Dict[int, FwdRev]:
    """
    Load reads from a BAM file and index them by the 'ch' tag.
    Returns a dictionary where keys are ch values and values are lists of sequences.
    """
    reads_map = {}
    with pysam.AlignmentFile(bam_path, "rb", threads=mp.cpu_count()) as bam:
        for record in tqdm(bam.fetch(until_eof=True), desc=f"npthr:{np_thr}, reading {bam_path}"):
            try:
                ch = int(record.get_tag("ch"))
                if int(record.get_tag("np")) < np_thr:
                    continue
                seq = record.query_sequence
                if ch not in reads_map:
                    reads_map[ch] = FwdRev()

                if record.query_name.endswith("fwd"):
                    reads_map[ch].fwd = seq
                if record.query_name.endswith("rev"):
                    reads_map[ch].rev = seq

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


def align_sequences(fwd: str, rev: str) -> Tuple[float, str]:
    """
    Align two sequences using mappy and return identity and CIGAR.
    """
    aligner = mappy.Aligner(seq=fwd, extra_flags=67108864 | 1048576,
                            k=9, w=7, best_n=10, n_threads=1)
    hits = list(aligner.map(rev))
    if not hits:
        return None

    best_hit = hits[0]
    assert best_hit.strand == 1

    stat_info = stat_alignment_info(
        best_hit, fwd, rev)
    return stat_info


def process_molecule(ch: int, fwd_rev: FwdRev) -> Tuple[int, float, int, int, int, float, dict, dict, dict]:
    """
    Worker function to align two sequences for a given molecule and calculate metrics.
    Returns: (ch, identity, mismatch, ins_cnt, del_cnt, phreq, mismatch_detail, ins_detail, del_detail)
    """
    stats = align_sequences(
        fwd=fwd_rev.fwd, rev=fwd_rev.rev)
    if stats is None:
        return None

    identity = stats["identity"]
    mismatch = stats["mismatch"]
    ins_cnt = stats["ins"]
    del_cnt = stats["del"]
    phreq = calculate_phred(identity=identity)

    return ch, identity, mismatch, ins_cnt, del_cnt, phreq, stats["mismatch_detail"], stats["ins_detail"], stats["del_detail"]


def _wrapper(args):
    return process_molecule(*args)


def main():
    """
        分析bystrand=True时, 同一 channel 正反链的结果差异
    """
    env_prepare.polars_env_init()

    parser = argparse.ArgumentParser(
        description="Analyze differences between Combined and Separate SMC consensus modes.")

    parser.add_argument("--bystrand-bam", type=str, required=True,
                        help="BAM file from Separate Mode (byStrand=True)")
    parser.add_argument("--np-thr", type=int, default=5)

    args = parser.parse_args()

    out_prefix = str(pathlib.Path(args.bystrand_bam).with_suffix(''))

    print(f"Loading bystrand  BAM: {args.bystrand_bam}...")
    fwd_rev_reads = load_bam_reads(args.bystrand_bam, np_thr=args.np_thr)

    # In combined mode, we expect one read per ch.
    # In separate mode, we expect potentially multiple (fwd/rev).
    # We'll simplify true_reads to just one sequence per ch for the purpose of this analysis.

    # print(fwd_rev_reads[48509])
    # print(fwd_rev_reads[48509].has_rev_fwd())

    true_reads_double = {ch: fwd_rev
                         for ch, fwd_rev in fwd_rev_reads.items() if fwd_rev.has_rev_fwd()}

    common_chs = true_reads_double.keys()
    print(f"Found {len(common_chs)} channels (ch).")

    mismatches = []
    insertions = []
    deletions = []
    phred_scores = []
    diff_counts = {i: 0 for i in range(11)}
    mismatch_counts = {i: 0 for i in range(11)}
    ins_counts = {i: 0 for i in range(11)}
    del_counts = {i: 0 for i in range(11)}

    total_mismatch_detail = {}
    total_ins_detail = {}
    total_del_detail = {}

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

    for res in tqdm(results, desc=f"aggr ..."):
        if res is None:
            continue
        ch, identity, mismatch, ins_cnt, del_cnt, phreq, mis_det, ins_det, del_det = res

        tot_diff_cnt = mismatch + ins_cnt + del_cnt
        diff_counts[min(tot_diff_cnt, 10)] += 1
        mismatch_counts[min(mismatch, 10)] += 1
        ins_counts[min(ins_cnt, 10)] += 1
        del_counts[min(del_cnt, 10)] += 1

        for k, v in mis_det.items():
            total_mismatch_detail[k] = total_mismatch_detail.get(k, 0) + v
        for k, v in ins_det.items():
            total_ins_detail[k] = total_ins_detail.get(k, 0) + v
        for k, v in del_det.items():
            total_del_detail[k] = total_del_detail.get(k, 0) + v

        mismatches.append(mismatch)
        insertions.append(ins_cnt)
        deletions.append(del_cnt)
        phred_scores.append(phreq)

        if identity > 0.9999:
            fwd_rev_eq += 1

    fwd_rev_eq_ratio = fwd_rev_eq / len(common_chs)
    print(f"fwd_rev_eq_ratio:{fwd_rev_eq_ratio:.2f}")

    # Print distributions as a Polars DataFrame
    total_n = len(common_chs)
    labels = [str(i) for i in range(10)] + [">=10"]

    df_dist = pl.DataFrame({
        "Label": labels,
        "Total (%)": [f"{diff_counts[i]/total_n:.4%}" if total_n > 0 else "0.0000%" for i in range(11)],
        "Mismatch (%)": [f"{mismatch_counts[i]/total_n:.4%}" if total_n > 0 else "0.0000%" for i in range(11)],
        "Insertion (%)": [f"{ins_counts[i]/total_n:.4%}" if total_n > 0 else "0.0000%" for i in range(11)],
        "Deletion (%)": [f"{del_counts[i]/total_n:.4%}" if total_n > 0 else "0.0000%" for i in range(11)],
    })

    print("\nDifference Distributions:")
    print(df_dist)

    print("\nBase-specific Error Distributions:")

    # Mismatches
    bases = ['A', 'C', 'G', 'T']
    mismatch_data = []
    for r in bases:
        row = {"Ref": r}
        for q in bases:
            row[q] = total_mismatch_detail.get((r, q), 0)
        mismatch_data.append(row)
    df_mismatch = pl.DataFrame(mismatch_data)
    print("\n--- Mismatches (Ref -> Query) ---")
    print(df_mismatch)

    # Insertions
    ins_total = sum(total_ins_detail.values())
    ins_data = []
    for b in bases:
        count = total_ins_detail.get(b, 0)
        ins_data.append({
            "Base": b,
            "Count": count,
            "Percentage": f"{count/ins_total:.2%}" if ins_total > 0 else "0.00%"
        })
    df_ins = pl.DataFrame(ins_data)
    print("\n--- Insertions (Query base) ---")
    print(df_ins)
    print(f"Total Insertions: {ins_total}")

    # Deletions
    del_total = sum(total_del_detail.values())
    del_data = []
    for b in bases:
        count = total_del_detail.get(b, 0)
        del_data.append({
            "Base": b,
            "Count": count,
            "Percentage": f"{count/del_total:.2%}" if del_total > 0 else "0.00%"
        })
    df_del = pl.DataFrame(del_data)
    print("\n--- Deletions (Ref base) ---")
    print(df_del)
    print(f"Total Deletions: {del_total}")

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
