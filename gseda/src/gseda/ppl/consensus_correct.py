#!/usr/bin/env python3
"""
Determine the correct GGC tandem repeat count in a consensus sequence
by comparing read alignments across candidate variants.

The strategy: the most accurate consensus should produce the fewest
indel + mismatch bases when all reads are aligned against it.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from typing import Tuple, List

import pysam


def swap_ggc_region(fa_path: str, repeat_count: int) -> str:
    """Replace the GGC region in the consensus FASTA with a new repeat count.

    The original GGC region spans positions 240-426 (0-based, exclusive of 426).
    This function replaces that substring with ``repeat_count * 'GGC'`` and writes
    the modified sequence to a temporary file.

    Parameters
    ----------
    fa_path : str
        Path to the original consensus FASTA file.
    repeat_count : int
        The desired number of GGC repeats (e.g. 60, 61, ..., 65).

    Returns
    -------
    str
        Path to the temporary modified FASTA file.
    """
    ref_start = 240
    ref_end = 426

    with pysam.FastxFile(fa_path) as fx:
        name = None
        seq = None
        for entry in fx:
            name = entry.name
            seq = entry.sequence
            break

    if seq is None:
        print("error: no sequence found in reference FASTA", file=sys.stderr)
        sys.exit(1)

    assert len(seq) >= ref_end, (
        f"reference length ({len(seq)}) is shorter than expected region end "
        f"({ref_end}); adjust ref_start/ref_end to match the actual sequence."
    )

    # Verify the original GGC region looks plausible
    original_ggc = seq[ref_start:ref_end]
    expected_ggc = "GGC" * 62
    if original_ggc != expected_ggc:
        print(
            f"warning: region [{ref_start}:{ref_end}] is not exactly 62 GGC repeats. "
            f"Got: {original_ggc[:20]}...{original_ggc[-10:]}",
            file=sys.stderr,
        )

    new_region = "GGC" * repeat_count
    new_seq = seq[:ref_start] + new_region + seq[ref_end:]

    # Write temporary FASTA
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    )
    chunk_size = 60
    tmp.write(f">{name}\n")
    for i in range(0, len(new_seq), chunk_size):
        tmp.write(new_seq[i : i + chunk_size] + "\n")
    tmp.close()
    return tmp.name


def count_errors_in_region(read: pysam.AlignedSegment, consensus_seq: str, ggc_length) -> Tuple[int, int, int]:
    """Count insertion, deletion, and mismatch bases for a single aligned read.

    Only reads that fully span the key region [240, 426) are considered.
    Uses ``read.get_aligned_pairs(matches_only=False)`` which already encodes
    CIGAR semantics: ``qpos is None / rpos is not None`` = deletion,
    ``qpos is not None / rpos is None`` = insertion.

    Parameters
    ----------
    read : pysam.AlignedSegment
        An aligned BAM record (read).
    consensus_seq : str
        The full consensus reference sequence (used to look up reference bases).

    Returns
    -------
    tuple[int, int, int]
        ``(insertion_bases, deletion_bases, mismatch_bases)`` for this read.
    """
    ref_start = 240
    ref_end = ref_start + ggc_length

    # Must fully span the region
    if read.reference_start > ref_start:
        return (0, 0, 0)
    if read.reference_end < ref_end:
        return (0, 0, 0)

    query_seq = read.query_sequence
    insertion_bases = 0
    deletion_bases = 0
    mismatch_bases = 0

    # Track current reference position using cursor-based logic.
    # An insertion "in" the key region if its position falls within
    # [ref_start - 1, ref_end) — i.e. it occurs at or just before the region.
    ref_cursor = -1  # sentinel; never a valid ref position

    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos is not None:
            ref_cursor = rpos

        # Skip reference positions strictly before the key region (including ref_start-1)
        if ref_cursor is not None and ref_cursor < ref_start - 1:
            continue

        # Insertion: query base exists but no corresponding reference position.
        # Count it as "in the key region" if it occurs between ref_cursor and ref_end.
        if qpos is not None and rpos is None:
            if ref_cursor is not None and (ref_cursor == ref_start - 1 or ref_cursor < ref_end):
                insertion_bases += 1
            continue

        # Deletion: reference position exists but no corresponding query base
        if qpos is None and rpos is not None:
            if ref_cursor >= ref_start and ref_cursor < ref_end:
                deletion_bases += 1
            continue

        # Both exist — check for mismatch
        if qpos is not None and rpos is not None:
            if ref_cursor >= ref_start and ref_cursor < ref_end:
                ref_base = consensus_seq[rpos].upper()
                query_base = query_seq[qpos].upper()
                if ref_base != query_base:
                    mismatch_bases += 1

    return (insertion_bases, deletion_bases, mismatch_bases)


def count_all_errors(bam_path: str, consensus_seq: str, ggc_length) -> Tuple[int, int, int]:
    """Count total indel + mismatch bases across all reads covering the key region.

    Parameters
    ----------
    bam_path : str
        Path to a BAM file (may be SAM; pysam will auto-detect).
    consensus_seq : str
        The full consensus reference sequence used during alignment.

    Returns
    -------
    tuple[int, int, int]
        ``(total_insertions, total_deletions, total_mismatches)`` summed across all reads.
    """
    total_ins = 0
    total_dels = 0
    total_mm = 0

    with pysam.AlignmentFile(bam_path, "rb" if bam_path.endswith(".bam") else "r", check_sq=False) as bam:
        for read in bam:
            ins, dels, mms = count_errors_in_region(read, consensus_seq, ggc_length=ggc_length)
            total_ins += ins
            total_dels += dels
            total_mm += mms

    return (total_ins, total_dels, total_mm)


def run_gsmm2(query_fq: str, consensus_fa: str, prefix: str) -> str:
    """Run gsmm2 align and return the path to the output BAM file.

    Parameters
    ----------
    query_fq : str
        Path to the query FASTQ file.
    consensus_fa : str
        Path to the (modified) consensus FASTA file.
    prefix : str
        Output file prefix — gsmm2 will generate ``{prefix}.sam`` and
        ``{prefix}.bam`` among other files.

    Returns
    -------
    str
        Path to the generated BAM file.
    """
    cmd = (
        f"gsmm2 align -q {query_fq} -t {consensus_fa} "
        f"-p {prefix} --noMar "
    ) # --rq-range 0.99:1.1
    print(f"running: {cmd}", file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    bam_path = f"{prefix}.bam"
    sam_path = f"{prefix}.sam"

    for candidate in (bam_path, sam_path):
        if os.path.exists(candidate):
            return candidate

    raise FileNotFoundError(
        f"gsmm2 did not produce output at {bam_path} or {sam_path}"
    )


if __name__ == "__main__":
    """
     有一堆 query 和 一个 consensus seq 序列
     consensus seq 序列中 有 62*(GGC) 的 short tandem repeats 序列

     我现在是不知道 这块区域 是 62, 63 还是 64 x。但是我想到了一个简单的策略来做这个事情。

     最正确的 consenses seq. 对于这块区域的 indel bases + mismatched 数量是最少的。

     请基于上述策略 实现代码。整体的思路如下所示
        遍历 60~65*(GGC)
            生成 consensus_seq.
            调用 gsmm2 align 将 query 和 consensus seq 进行比对。
                 gsmm2 align -q 257-0_STR_filtered.cutBarcodes.fastq  -t barcode257-0-ref-add5gcc.fasta  -p filtered-257-rqQ20--ref-splitted-add5gcc --noMar --rq-range 0.99:1.1
                    -q是  query filename
                    -t 是 consensus seq 的文件名
                    -p 是输出文件的 前缀
                    后面的参数保留
            然后 判断关键区域的 indel bases + mismatch bases 个数 。记录 并输出

        找到 indel bases + mismatch bases 最少的 GCC repeats 数量，并输出


    补充：
    GAAGGACACACGAGGCTGCTTCGCTGCACACCCGAGAAAGTTTCAGCCAAACTTCGGGCGGCGGCTGAGGCGGCGGCCGAGGAGCGGCGGACTCGGGGCGCGGGGAGTCGAGGCATTTGCGCCTGTGCTTCGGACCGTAGCGCCAGGGCCTGAGCCTTTGAAGCAGGAGGAGGGGAGGAGAGAGTGGGGCTCCTCTATCGGGACCCCCTCCCCGAATTCGCCACCATGTGGATCTGCCCAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGAGGAGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGAGGCGGCGGAGGAGGCGGCGACCGAGAAGATGCCCGCCCTGCGCCGCTCTGCTGTGGGCGCTGCTGGCGCTCTGGCTGTGCTGCGCGACCCCCGCGCATGCATTGCAGTGTCGAGATGGCTATGAACCCTGTGGGATCCGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCT
    这个是 原始的 consensus_seq. 其中 240~426 位置是 62 个 GGC. 你要做的就是将这 62 个 GGC 换成 60 个，61个，62个，63个，64个，65 个。
    下面这个代码可以用来参考。用来 计算 GGC repeats 区域的 indel+mismatch base 数量
    def extract_subseq(read, ref_start, ref_end):
        if read.reference_start > ref_start:
            return ""

        if read.reference_end < ref_end:
            return ""

        query = read.query_sequence
        chars = []
        q_cursor = -1
        r_cursor = -1
        for qpos, rpos in read.get_aligned_pairs(matches_only=False):
            if qpos is not None:
                q_cursor = qpos

            if rpos is not None:
                r_cursor = rpos

            if q_cursor is None and r_cursor is None:
                continue

            if r_cursor < (ref_start - 1):
                continue

            if r_cursor == (ref_start - 1) and rpos is not None:
                continue

            if r_cursor < ref_start:
                continue

            if r_cursor >= ref_end:
                break

            # deletion
            if qpos is None:
                continue

            chars.append(query[qpos])

        return "".join(chars)


    实现的过程中，保留需求描述部分，也就是这个注释
    """

    parser = argparse.ArgumentParser(
        prog="consensus_correct",
        description=(
            "Determine the correct GGC repeat count in a consensus sequence. "
            "Runs gsmm2 align against candidate variants (60-65 repeats) and reports "
            "the variant with the fewest indel + mismatch bases."
        ),
    )
    parser.add_argument("query_fastq", help="Query FASTQ file containing aligned reads")
    parser.add_argument("consensus_fasta", help="Consensus FASTA reference with 62 GGC repeats at positions 240-426")

    args = parser.parse_args()

    query_fq: str = args.query_fastq
    ref_fasta: str = args.consensus_fasta

    # Load the original consensus sequence for later base comparison
    with pysam.FastxFile(ref_fasta) as fx:
        orig_seq = None
        for entry in fx:
            orig_seq = entry.sequence
            break

    if orig_seq is None:
        print("error: empty reference FASTA", file=sys.stderr)
        sys.exit(1)

    results: List[Tuple[int, int, int, int]] = []  # (repeat_count, ins, dels, mismatches)
    tmp_dir = tempfile.mkdtemp(prefix="consensus_correct_")

    for count in range(60, 66):
        print(f"\n--- Candidate: {count} GGC repeats ---", file=sys.stderr)
        fasta_path = swap_ggc_region(ref_fasta, count)
        prefix = os.path.join(tmp_dir, f"candidate_{count}")

        bam_path = run_gsmm2(query_fq, fasta_path, prefix)
        ins, dels, mm = count_all_errors(bam_path, orig_seq, ggc_length=count*3)
        results.append((count, ins, dels, mm))
        print(f"  ins={ins}, dels={dels}, mm={mm}", file=sys.stderr)

    # Print summary table
    print("\n" + "=" * 50)
    print(f"{'Repeat':>8} {'Ins':>6} {'Del':>6} {'MM':>6} {'Total':>7}")
    print("-" * 35)
    best_total = min(t for _, _, _, t in results)
    for count, ins, dels, mm in results:
        total = ins + dels + mm
        marker = " <-- best" if total == best_total else ""
        print(f"{count:>8} {ins:>6} {dels:>6} {mm:>6} {total:>7}{marker}")

    best_count, best_ins, best_del, best_mm = min(results, key=lambda r: r[3])
    best_total_val = best_ins + best_del + best_mm
    print(f"\nBest GGC repeat count: {best_count} (indel+mismatch bases: {best_total_val})")
