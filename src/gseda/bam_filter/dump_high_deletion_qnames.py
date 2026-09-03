"""Dump query names with high deletion ratio within a reference range from BAM.

Walks each read's CIGAR operations to compute the fraction of reference bases
covered by the read that are deletions (op code 2), for reads overlapping a
specified [start, end) range on a contig.
"""

import argparse
import os
import sys
from tqdm import tqdm

import pysam


def calc_deletion_ratio(read, start: int, end: int) -> float:
    """
    Calculate deletion ratio within reference interval [start, end).

    deletion_ratio =
        deleted reference bases inside the overlapped region
        ----------------------------------------------------
           reference overlap length between read and region

    Returns
    -------
    float
    """

    ref_start = read.reference_start
    ref_end = _read_ref_end(read)

    # overlap on reference
    ovlp_start = max(ref_start, start)
    ovlp_end = min(ref_end, end)

    if ovlp_start >= ovlp_end:
        return 0.0

    overlap_len = ovlp_end - ovlp_start

    deletions = 0

    # current reference position
    ref_pos = ref_start

    for op, length in read.cigartuples:

        ref_len = _ref_advance(op, length)

        op_ref_start = ref_pos
        op_ref_end = ref_pos + ref_len

        # operation completely before region
        if op_ref_end <= ovlp_start:
            ref_pos = op_ref_end
            continue

        # already beyond region
        if op_ref_start >= ovlp_end:
            break

        # overlap on reference
        covered = min(op_ref_end, ovlp_end) - max(op_ref_start, ovlp_start)

        if covered > 0 and op == 2:      # CIGAR D
            deletions += covered

        ref_pos = op_ref_end

    return deletions / overlap_len

def _ref_advance(op: int, length: int) -> int:
    # M D = X consume reference
    if op in (0, 2, 7, 8):
        return length

    # I S H consume no reference
    return 0


def _read_ref_end(read) -> int:
    """Return the exclusive end coordinate of a read on the reference."""
    offset = 0
    for op_code, length in read.cigartuples:
        if op_code in (0, 2, 3, 7, 8):   # ref-consuming ops
            offset += length
    return read.reference_start + offset


def _ref_advance(op_code: int, length: int) -> int:
    """Return how much this CIGAR op advances the reference coordinate."""
    return length if op_code in (0, 2, 3, 7, 8) else 0


def main(args):
    """Identify reads with deletion_ratio > threshold within [start, end) on contig."""
    n_read = 0
    results = []  # list of (query_name, ratio)

    with pysam.AlignmentFile(
        args.bam, mode="rb", threads=args.threads, check_sq=False
    ) as in_bam:
        for read in tqdm(in_bam.fetch(contig=args.contig, start=args.start, end=args.end, until_eof=True)):
            if read.is_unmapped or read.is_supplementary or read.is_duplicate:
                continue

            n_read += 1
            ratio = calc_deletion_ratio(read, args.start, args.end)
            if ratio > args.threshold:
                results.append((read.query_name, ratio))

    results.sort(key=lambda x: x[1], reverse=True)
    sys.stderr.write(f"total_reads={n_read} filtered={len(results)} ratio:{len(results) / n_read}\n")

    tag = f"{args.contig}:{args.start}-{args.end}"
    print(f"{tag}\tquery_name\tratio", file=sys.stdout)
    for qname, ratio in results:
        print(f"{tag}\t{qname}\t{ratio:.4f}", file=sys.stdout)
        
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--bam", type=str, required=True)
    parser.add_argument("--contig", type=str, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("-t", "--threads", type=int,
                        default=min(40, (os.cpu_count() or 1) // 2))
    main(parser.parse_args())
