import pysam
import argparse
from tqdm import tqdm
import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)
from range_identity import RangeIdentityCalculator


def calculate_low_identity_windows(
    bam_file, window_size, step_size, identity_threshold
):
    # 打开 BAM 文件
    bam = pysam.AlignmentFile(bam_file, "rb", threads=40)

    total_queries = 0
    queries_with_low_identity = 0
    total_query_length = 0
    total_low_identity_length = 0

    for read in tqdm(bam, desc=f"reading {bam_file}"):
        if read.is_unmapped:
            continue  # 跳过未比对的 reads

        total_queries += 1
        query_length = read.query_length or 0
        total_query_length += query_length

        # 计算每个窗口的 identity
        sequence = read.query_sequence
        low_identity_found = False
        low_identity_length = 0

        query_start = read.query_alignment_start
        query_end = read.query_alignment_end

        # 遍历窗口
        range_ideneity_calc = RangeIdentityCalculator(read.cigartuples)
        for start in range(query_start, query_end, step_size):
            end = min(start + window_size, query_end)

            identity = range_ideneity_calc.range_identity(start=start, end=end)

            if identity < identity_threshold:
                low_identity_found = True
                low_identity_length += end - start

        if low_identity_found:
            queries_with_low_identity += 1
            total_low_identity_length += low_identity_length

    bam.close()

    # 计算统计结果
    proportion_queries_with_low_identity = (
        queries_with_low_identity / total_queries if total_queries > 0 else 0
    )
    proportion_low_identity_length = (
        total_low_identity_length / total_query_length if total_query_length > 0 else 0
    )

    return {
        "total_queries": total_queries,
        "queries_with_low_identity": queries_with_low_identity,
        "proportion_queries_with_low_identity": proportion_queries_with_low_identity,
        "total_query_length": total_query_length,
        "total_low_identity_length": total_low_identity_length,
        "proportion_low_identity_length": proportion_low_identity_length,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BAM file low identity window statistics."
    )
    parser.add_argument(
        "-b", "--bam", required=True, help="Input BAM file (sorted and indexed)."
    )
    parser.add_argument(
        "-w",
        "--window_size",
        type=int,
        default=50,
        help="Size of the sliding window (default: 50).",
    )
    parser.add_argument(
        "-s",
        "--step_size",
        type=int,
        default=10,
        help="Step size for sliding the window (default: 10).",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.8,
        help="Identity threshold for low identity windows (default: 0.8).",
    )
    args = parser.parse_args()

    stats = calculate_low_identity_windows(
        args.bam, args.window_size, args.step_size, args.threshold
    )

    print("BAM File: ", args.bam)
    print("Total Queries: ", stats["total_queries"])
    print("Queries with Low Identity Windows: ", stats["queries_with_low_identity"])
    print(
        "Proportion of Queries with Low Identity Windows: {:.2%}".format(
            stats["proportion_queries_with_low_identity"]
        )
    )
    print("Total Query Length: ", stats["total_query_length"])
    print("Total Low Identity Length: ", stats["total_low_identity_length"])
    print(
        "Proportion of Low Identity Length: {:.2%}".format(
            stats["proportion_low_identity_length"]
        )
    )
