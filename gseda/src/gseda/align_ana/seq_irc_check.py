import os
import argparse
import time
import mappy
import re
import pysam
from multiprocessing import Process, Queue, Value, Lock, cpu_count
import signal
import sys

# =========================
# 全局计数器（多进程安全）
# =========================
TOTAL_READS = Value('i', 0)
SPLIT_READS = Value('i', 0)
COUNTER_LOCK = Lock()

# Sentinel for queue termination
_SENTINEL = None


def calculate_identity_from_cigar(cigar_string):
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    match_count = 0
    mismatch_count = 0
    total_aligned = 0

    for length_str, operation in matches:
        length = int(length_str)
        if operation == '=':
            match_count += length
            total_aligned += length
        elif operation == 'X':
            mismatch_count += length
            total_aligned += length
        elif operation in ('I', 'D'):
            total_aligned += length

    if total_aligned == 0:
        return 0.0
    return match_count / total_aligned


def try_split_seq_2(seq: str):
    seq_len = len(seq)
    if seq_len < 200 or seq_len > 20000:
        return False

    mid = seq_len // 2
    seq1 = seq[:mid]
    seq2 = seq[mid:]

    aligner = mappy.Aligner(seq=seq1, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    for hit in aligner.map(seq2):
        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage1 = (hit.r_en - hit.r_st) / len(seq1)
        coverage2 = (hit.q_en - hit.q_st) / len(seq2)
        if identity > 0.80 and coverage1 > 0.85 and coverage2 > 0.85:
            return True
    return False


def main_cli():
    seq = "CCCCGGGGAAAATTTTGAGAAGAGAGCCCCGCACTTCCACCACCAGCTCCTCCATCTTCTCTTCAGCCCTGCTAGCGCCGGGAGCCCGCCCCCGAGAGGTGGGCTGCGGGCGCTCGAGGCCCAGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTCCGCCGCCGCCGCCGCCGCCGCCGCCGCGCTGCCGCACGCCCCCTGGCAGCGGCGCCTCCGTCACCGCCGCCGCCCGCGCTCGCCGTCGGCCCGCCGCCCGCTCAGAGGCGGCCCTCCACCGGAAGTGAAACCGAAACGGAGCTGAGCGCCCCGCACTTCCACCACCAGCTCCTCCATCTTCTCTTCGGCCCTGCTAGCGCCGGGAGCCCGCCCCCGAGAGGTGGGCTGCGGGCGCTCGAGGCCCAGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTCCGCCGCCGCCGCCGCCGCCGCCGCCTCCGCCGCCGCCGCCGCCGCCGCCGCCGCGCTGCCGCACGCCCCCTGGCAGCGGCGCCTCCGTCACCGCCGCCGCCCGCGCTCGCCGTCGGCCCGCCGCCCGCTCAGAGGCGGCCCTCCACCGGAAGTGAAACCGAAACGGAGCTGAGCATCTCTCTCAAAATTTTCCCCGGGG"
    tag = try_split_seq_2(seq=seq)
    print(f"it is ok? = {tag}")


if __name__ == "__main__":
    main_cli()
