import pysam
import argparse
from tqdm import tqdm


def calculate_long_insertions(bam_file, min_insertion_length):
    # 打开 BAM 文件
    bam = pysam.AlignmentFile(bam_file, "rb", threads=40)

    total_queries = 0
    queries_with_long_insertions = 0
    total_query_length = 0
    total_long_insertion_length = 0

    for read in tqdm(bam.fetch(), desc=f"reading {bam_file}"):
        # 统计总 query 数量
        total_queries += 1

        if read.is_unmapped:
            continue  # 跳过未比对的 reads

        query_length = read.query_length or 0
        total_query_length += query_length

        # 遍历 CIGAR 操作
        has_long_insertion = False
        long_insertion_length = 0

        for op, length in read.cigartuples or []:
            if op == 1 and length >= min_insertion_length:  # 插入 (CIGAR code 1)
                has_long_insertion = True
                long_insertion_length += length

        if has_long_insertion:
            queries_with_long_insertions += 1
            total_long_insertion_length += long_insertion_length

    bam.close()

    # 计算统计结果
    proportion_queries_with_long_insertions = (
        queries_with_long_insertions / total_queries if total_queries > 0 else 0
    )
    proportion_long_insertion_length = (
        total_long_insertion_length / total_query_length
        if total_query_length > 0
        else 0
    )

    return {
        "total_queries": total_queries,
        "queries_with_long_insertions": queries_with_long_insertions,
        "proportion_queries_with_long_insertions": proportion_queries_with_long_insertions,
        "total_query_length": total_query_length,
        "total_long_insertion_length": total_long_insertion_length,
        "proportion_long_insertion_length": proportion_long_insertion_length,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BAM file long insertion statistics.")
    parser.add_argument(
        "-b", "--bam", required=True, help="Input BAM file (sorted and indexed)."
    )
    parser.add_argument(
        "-l",
        "--length",
        type=int,
        default=10,
        help="Minimum length to consider an insertion as long (default: 10).",
    )
    args = parser.parse_args()

    stats = calculate_long_insertions(args.bam, args.length)

    print("BAM File: ", args.bam)
    print("Total Queries: ", stats["total_queries"])
    print("Queries with Long Insertions: ", stats["queries_with_long_insertions"])
    print(
        "Proportion of Queries with Long Insertions: {:.2%}".format(
            stats["proportion_queries_with_long_insertions"]
        )
    )
    print("Total Query Length: ", stats["total_query_length"])
    print("Total Long Insertion Length: ", stats["total_long_insertion_length"])
    print(
        "Proportion of Long Insertion Length: {:.2%}".format(
            stats["proportion_long_insertion_length"]
        )
    )
