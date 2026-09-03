"""Analyze upstream context base distribution of a motif in BAM reads."""

import argparse
import os
import sys
from collections import Counter
from tqdm import tqdm

cur_path = os.path.abspath(__file__)
cur_dir = os.path.dirname(cur_path)
prev_dir = os.path.dirname(cur_dir)
sys.path.insert(0, str(cur_dir))
sys.path.insert(0, str(prev_dir))

import pysam  # noqa: E402


def scan_motif_context(
    bam_path: str,
    motif: str,
    context_len: int,
    top_n: int = 20,
):
    total_reads = 0
    total_matches = 0
    ctx_counter: Counter = Counter()

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=os.cpu_count()) as bam:
        for record in tqdm(bam.fetch(until_eof=True), desc=f"reading"):
            seq = record.query_sequence
            if seq is None:
                continue
            total_reads += 1

            start = 0
            while True:
                pos = seq.find(motif, start)
                if pos == -1:
                    break

                ctx_start = pos - context_len
                if ctx_start >= 0:
                    ctx = seq[ctx_start:pos]
                else:
                    ctx = "N" * (-ctx_start) + seq[0:pos]

                ctx_counter[ctx] += 1
                total_matches += 1
                start = pos + 1

    ctx_total = sum(ctx_counter.values()) if ctx_counter else 1

    print("=" * 40)
    print("=== Motif Context Base Distribution ===")
    print(f"Motif: {motif}")
    print(f"Context length: {context_len}")
    print(f"Total reads scanned: {total_reads}")
    print(f"Total motif matches: {total_matches}")
    print()

    print(f"{'Context':<12} {'Count':>8} {'Frequency':>10}")
    print("-" * 40)
    for ctx, cnt in ctx_counter.most_common(top_n):
        freq = cnt / ctx_total * 100
        print(f"{ctx:<12} {cnt:>8} {freq:>9.1f}%")


def main_cli():
    parser = argparse.ArgumentParser(description="Motif upstream context base distribution in BAM reads")
    parser.add_argument("--bam", required=True, help="Path to BAM file")
    parser.add_argument("--motif", default="T" * 13 + "GAACG", help="Motif sequence (default: TTT TTT TTT TTT TGAACG)")
    parser.add_argument("--context", type=int, default=6, help="Upstream context length (default: 6)")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top contexts to display (default: 20)")
    args = parser.parse_args()

    scan_motif_context(args.bam, args.motif, args.context, args.top_n)


if __name__ == "__main__":
    main_cli()
