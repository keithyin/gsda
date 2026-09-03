#!/usr/bin/env python3

import os
import argparse
import pysam
from collections import defaultdict
from tqdm import tqdm


def get_barcode(path):
    name = os.path.basename(path)
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def collect_fastq_files(input_dirs):
    fastqs = []

    for d in input_dirs:
        demux_path = os.path.join(d, "demuxed")

        if not os.path.exists(demux_path):
            print(f"[WARN] skip {d}, no demuxed/")
            continue

        for f in os.listdir(demux_path):
            if f.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                fastqs.append(os.path.join(demux_path, f))

    return fastqs


def group_by_barcode(fastq_files):
    """
    barcode -> [file1, file2, ...]
    """
    groups = defaultdict(list)

    for f in fastq_files:
        bc = get_barcode(f)
        groups[bc].append(f)

    return groups


def count_reads(paths):
    total = 0
    for p in paths:
        with pysam.FastxFile(p) as f:
            for _ in f:
                total += 1
    return total


def merge_fastq(paths, out_path):
    """
    流式合并 FASTQ
    """
    with open(out_path, "w") as out:
        for p in paths:
            with pysam.FastxFile(p) as f:
                for entry in f:
                    out.write(f"@{entry.name}\n")
                    out.write(f"{entry.sequence}\n")
                    out.write("+\n")
                    out.write(f"{entry.quality}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--min-reads", type=int, default=1000)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    fastq_files = collect_fastq_files(args.inputs)
    print(f"Found {len(fastq_files)} FASTQ files")

    groups = group_by_barcode(fastq_files)
    print(f"Found {len(groups)} barcodes")

    kept = 0
    dropped = 0

    for bc, paths in tqdm(groups.items()):
        total_reads = count_reads(paths)

        if total_reads < args.min_reads:
            dropped += 1
            continue

        out_path = os.path.join(args.output, f"{bc}.fastq")

        merge_fastq(paths, out_path)

        kept += 1

    print("\n=== Summary ===")
    print(f"Kept barcodes   : {kept}")
    print(f"Dropped barcodes: {dropped}")


if __name__ == "__main__":
    main()
