import pysam
import argparse
from tqdm import tqdm


def bam_to_fastx(bam_path, fastx_path: str, rq_threshold):
    tot = 0
    dumped = 0
    output_fq = fastx_path.endswith("fastq")

    with pysam.AlignmentFile(
        bam_path, "rb", threads=40, check_sq=False
    ) as bam_file, open(fastx_path, "w") as fastx_out:
        for read in tqdm(
            bam_file.fetch(until_eof=True), desc=f"dumping {bam_path} to {fastx_path}"
        ):
            tot += 1

            # 尝试获取 rq 字段
            try:
                rq = read.get_tag("rq")
                if rq < rq_threshold:
                    continue
            except KeyError:
                pass

            # 构造 FASTQ 格式
            name = read.query_name
            seq = read.query_sequence
            qual = read.qual  # 转换为 ASCII 的质量字符串

            dumped += 1

            if seq is None or qual is None:
                continue  # 有可能 read 被软裁剪或缺失，跳过
            if output_fq:
                fastx_out.write(f"@{name}\n{seq}\n+\n{qual}\n")
            else:
                fastx_out.write(f">{name}\n{seq}\n")

    print(f"Tot:{tot}, dumped:{dumped}, ratio:{dumped / tot: .4f}")
    print(f"转换完成，输出文件: {fastx_path}")


def main_cli():
    parser = argparse.ArgumentParser(description="Convert BAM to FASTQ with rq filter.")
    parser.add_argument("bam", help="Input BAM file path")
    parser.add_argument("fastx", help="Output FASTX file path. fasta/fastq")
    parser.add_argument(
        "--rq", type=float, default=0.0, help="Minimum rq threshold (default: 0.0)"
    )

    args = parser.parse_args()
    assert args.fastx.endswith("fastq") or args.fastx.endswith("fasta"), "only fastq/fasta are supported"
    bam_to_fastx(args.bam, args.fastx, args.rq)


if __name__ == "__main__":

    main_cli()
