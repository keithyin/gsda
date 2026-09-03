#!/usr/bin/env python3

import os
import argparse
import pysam
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
import math


def get_barcode(path):
    """
    从 fastq 文件名提取 barcode
    """
    name = os.path.basename(path)
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


PHRED_TO_P = [10 ** (-q / 10) for q in range(94)]  # 一般 Q 不超过 93


def calc_read_qscore(qual_str):
    """
    基于 error rate 计算 read 的真实 Q
    """
    if not qual_str:
        return 0

    error_sum = 0.0
    for c in qual_str:
        q = ord(c) - 33
        error_sum += PHRED_TO_P[q]

    mean_error = error_sum / len(qual_str)

    return -10 * math.log10(mean_error)


def process_fastq(path):
    """
    单文件处理
    """
    barcode = get_barcode(path)

    total = 0
    q20 = 0
    q30 = 0

    try:
        with pysam.FastxFile(path) as f:
            for entry in f:
                total += 1

                qscore = calc_read_qscore(entry.quality)

                # 等价于 avg >= 20 / 30（避免浮点）
                if qscore >= 20:
                    q20 += 1
                if qscore >= 30:
                    q30 += 1

    except Exception as e:
        print(f"[ERROR] {path}: {e}")
        return barcode, 0, 0, 0

    return barcode, total, q20, q30


def collect_fastq_files(input_dirs):
    """
    收集多个目录下 demuxed/*.fastq
    """
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


def merge_results(results):
    """
    按 barcode 聚合（支持多个目录合并）
    """
    stat = defaultdict(lambda: [0, 0, 0])

    for barcode, total, q20, q30 in results:
        stat[barcode][0] += total
        stat[barcode][1] += q20
        stat[barcode][2] += q30

    return stat


def pretty_print(stat):
    rows = []
    for bc, (total, q20, q30) in stat.items():
        if total == 0:
            continue

        rows.append({
            "barcode": bc,
            "total": total,
            "q20": q20,
            "q30": q30,
            "q20_ratio": q20 / total,
            "q30_ratio": q30 / total,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("barcode")

    print(df.to_string(index=False))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        required=True,
        help="multiple input dirs (each contains demuxed/)",
    )
    parser.add_argument("-t", "--threads", type=int, default=cpu_count())
    args = parser.parse_args()

    fastq_files = collect_fastq_files(args.inputs)

    print(f"Found {len(fastq_files)} FASTQ files")

    with Pool(args.threads) as pool:
        results = list(
            tqdm(pool.imap(process_fastq, fastq_files), total=len(fastq_files))
        )

    stat = merge_results(results)

    # print("\n=== Barcode Statistics ===")
    # print("barcode\ttotal\tq20\tq30\tq20_ratio\tq30_ratio")

    pretty_print(stat)

    # for bc, (total, q20, q30) in sorted(stat.items()):
    #     if total == 0:
    #         continue

    #     print(
    #         f"{bc}\t{total}\t{q20}\t{q30}\t"
    #         f"{q20/total:.4f}\t{q30/total:.4f}"
    #     )


if __name__ == "__main__":
    main()
