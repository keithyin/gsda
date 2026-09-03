#!/usr/bin/env python3

import argparse
import pysam
import os
import multiprocessing as mp
from tqdm import tqdm


def export_bam_subset(bam_path, min_len, max_len):
    # 输出文件名：input.bam -> input.<min>-<max>.bam
    base, ext = os.path.splitext(bam_path)
    out_bam = f"{base}.{min_len}-{max_len}{ext}"

    bam_in = pysam.AlignmentFile(
        bam_path, "rb", check_sq=False, threads=mp.cpu_count())
    bam_out = pysam.AlignmentFile(
        out_bam, "wb", header=bam_in.header, threads=mp.cpu_count())

    total = 0
    selected = 0

    for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {bam_path}"):
        total += 1

        length = read.query_length
        if length is None:
            continue

        if min_len <= length < max_len:
            bam_out.write(read)
            selected += 1

    bam_in.close()
    bam_out.close()

    ratio = selected / total if total > 0 else 0

    print(f"Input BAM: {bam_path}")
    print(f"Output BAM: {out_bam}")
    print(f"Length range: [{min_len}, {max_len})")
    print(f"Reads in range: {selected}")
    print(f"Total reads: {total}")
    print(f"Ratio: {ratio:.4%}")


def main():
    parser = argparse.ArgumentParser(
        description="Export BAM subset by read length range")
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("--min", type=int, required=True,
                        help="Minimum read length")
    parser.add_argument("--max", type=int, required=True,
                        help="Maximum read length")

    args = parser.parse_args()

    export_bam_subset(args.bam, args.min, args.max)


if __name__ == "__main__":
    main()
