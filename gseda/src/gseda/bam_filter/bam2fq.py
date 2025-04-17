import pysam
import argparse
from tqdm import tqdm

def bam_to_fastq(bam_path, fastq_path, rq_threshold):
    tot = 0
    dumped = 0
    with pysam.AlignmentFile(bam_path, "rb", threads=40, check_sq=False) as bam_file, open(fastq_path, "w") as fastq_out:
        for read in tqdm(bam_file.fetch(until_eof=True), desc=f"dumping {bam_path} to {fastq_path}"):

            # 尝试获取 rq 字段
            try:
                rq = read.get_tag("rq")
            except KeyError:
                continue  # 没有 rq 字段，跳过

            if rq < rq_threshold:
                continue

            # 构造 FASTQ 格式
            name = read.query_name
            seq = read.query_sequence
            qual = read.qual  # 转换为 ASCII 的质量字符串

            if seq is None or qual is None:
                continue  # 有可能 read 被软裁剪或缺失，跳过

            fastq_out.write(f"@{name}\n{seq}\n+\n{qual}\n")

    print(f"转换完成，输出文件: {fastq_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BAM to FASTQ with rq filter.")
    parser.add_argument("bam", help="Input BAM file path")
    parser.add_argument("fastq", help="Output FASTQ file path")
    parser.add_argument("--rq", type=float, default=0.0, help="Minimum rq threshold (default: 0.0)")

    args = parser.parse_args()
    bam_to_fastq(args.bam, args.fastq, args.rq)