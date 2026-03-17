#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import multiprocessing as mp


def phred_to_accuracy(qual):
    """
    将 FASTQ 的 phred quality 转换为 accuracy (rq)
    """
    if len(qual) == 0:
        return 0

    q = np.array(qual)
    error_prob = 10 ** (-q / 10)
    acc = 1 - np.mean(error_prob)

    return acc


def read_bam_lengths(path, rq_thresholds):
    """
    从 BAM 读取 read length, 并根据 rq 分组
    """
    result = {rq: [] for rq in rq_thresholds}

    with pysam.AlignmentFile(path, check_sq=False, threads=mp.cpu_count()) as bam_file:

        for read in tqdm(bam_file.fetch(until_eof=True), desc=f"reading {path}"):

            length = read.query_length
            try:
                rq = read.get_tag("rq")
            except KeyError:
                continue

            for thr in rq_thresholds:
                if rq >= thr:
                    result[thr].append(length)

    return result


def read_fastq_lengths(path, rq_thresholds):
    """
    从 FASTQ 读取 length 并计算 rq
    """
    result = {rq: [] for rq in rq_thresholds}

    with pysam.FastxFile(filename=path) as f:
        for rec in tqdm(f, desc=f"reading {path}"):
            length = len(rec.seq)
            qual = rec.quality
            rq = phred_to_accuracy(qual)

            for thr in rq_thresholds:
                if rq >= thr:
                    result[thr].append(length)

    return result


def read_fasta_lengths(path):
    """
    FASTA 只统计长度
    """
    lengths = []
    with pysam.FastxFile(filename=path) as f:
        for rec in tqdm(f, desc=f"reading {path}"):
            lengths.append(len(rec.seq))

    return lengths


def plot_distribution(data_dict, output, max_length=99999999999999, length_wrap=99999999999999):

    fig, ax = plt.subplots(figsize=(8, 6))

    all_lengths = []
    for (rq, lengths) in data_dict.items():
        all_lengths.extend(lengths)
    all_lengths = np.array(all_lengths)
    real_max_length = int(np.max(all_lengths))
    max_length = min(max_length, real_max_length)

    for (rq, lengths) in data_dict.items():

        if len(lengths) == 0:
            ax.set_title(f"rq >= {rq} (no reads)")
            continue

        lengths = [l for l in lengths if l < max_length]

        lengths = np.array(lengths)
        lengths = np.clip(lengths, a_min=0, a_max=length_wrap)
        lengths = np.sort(lengths)
        num_sample = len(lengths)
        cum = 0
        for (start, end) in [(0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 2500)]:
            region_count = np.sum((lengths >= start) & (lengths < end))
            cum += region_count
            print(f"rq ≥ {rq}, {start}-{end}: {region_count}, {region_count/num_sample * 100:.2f}%, {cum/num_sample * 100:.2f}%")
        
        print("--------------------------")

        plt.sca(ax)
        plt.yticks(list(map(lambda x: x*0.1,   range(0, 10))))
        y = np.arange(1, len(lengths) + 1) / len(lengths)
        # ax.hist(lengths, bins=100)
        ax.plot(lengths, y, label=f"rq ≥ {rq}")

    ax.set_title(f"CDF")
    ax.grid(True)
    ax.set_ylabel("Prob")
    plt.legend()
    plt.tight_layout()
    plt.xticks(list(range(0, length_wrap, 200)))
    plt.tick_params(axis="x", rotation=45)

    plt.savefig(output, dpi=200)


def main_cli():

    parser = argparse.ArgumentParser()

    parser.add_argument("input")

    parser.add_argument(
        "--min-rq",
        nargs="+",
        type=float,
        default=[0],
        help="rq thresholds",
    )

    parser.add_argument(
        "--max-length-filter",
        type=int,
        default=99999999999999999,
        help="length threshold",
    )

    parser.add_argument(
        "--max-length-wrap",
        type=int,
        default=5000,
        help="length threshold",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="length_distribution.png",
    )

    args = parser.parse_args()

    rq_thresholds = sorted(args.min_rq)

    if args.input.endswith(".bam"):
        data = read_bam_lengths(args.input, rq_thresholds)

    elif args.input.endswith(".fastq") or args.input.endswith(".fq"):
        data = read_fastq_lengths(args.input, rq_thresholds)

    elif args.input.endswith(".fasta") or args.input.endswith(".fa"):
        lengths = read_fasta_lengths(args.input)
        data = {0: lengths}

    else:
        raise ValueError("unsupported file format")

    plot_distribution(data, args.output,
                      args.max_length_filter, args.max_length_wrap)


if __name__ == "__main__":
    main_cli()
