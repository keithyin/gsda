#!/bin/bash

adapter_file=*_adapter.bam
output_prefix="${adapter_file%_adapter.bam}.bystrand"
echo "processing $adapter_file -> $output_prefix"
IMAGE_BASE="192.168.3.38:5000/algo/adacus"
# 获取所有匹配的镜像，按创建时间倒序，取第一行
LATEST_IMAGE=$(docker images --format "table {{.Repository}}:{{.Tag}}\t{{.CreatedAt}}" | \
               grep "^${IMAGE_BASE}:" | \
               sort -k2 -r | \
               head -n1 | \
               awk '{print $1}')
docker run -it -v `pwd`:/data --gpus all \
    ${LATEST_IMAGE} \
        smc \
        /data/${adapter_file} \
        /data/${output_prefix} \
        --maxLength 10000000 \
        --minLength 10 \
        --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4 \
        --byStrand 1

cat > bam2fx.py << 'EOF'
import pysam
import argparse
from tqdm import tqdm
import multiprocessing as mp


def bam_to_fastx(bam_path, fastx_path: str, rq_threshold):
    tot = 0
    dumped = 0
    output_fq = fastx_path.endswith("fastq")

    with pysam.AlignmentFile(
        bam_path, "rb", threads=mp.cpu_count(), check_sq=False
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
    parser = argparse.ArgumentParser(
        description="Convert BAM to FASTQ with rq filter.")
    parser.add_argument("bam", help="Input BAM file path")
    parser.add_argument("fastx", help="Output FASTX file path. fasta/fastq")
    parser.add_argument(
        "--rq", type=float, default=0.0, help="Minimum rq threshold (default: 0.0)"
    )

    args = parser.parse_args()
    assert args.fastx.endswith("fastq") or args.fastx.endswith(
        "fasta"), "only fastq/fasta are supported"
    bam_to_fastx(args.bam, args.fastx, args.rq)


if __name__ == "__main__":

    main_cli()
EOF

bystrand_filename=${output_prefix}.bam
docker exec -it gseq_product_dna python `pwd`/bam2fx.py $bystrand_filename ${output_prefix}.q8.fastq --rq 0.84151
docker exec -it gseq_product_dna python `pwd`/bam2fx.py $bystrand_filename ${output_prefix}.q10.fastq --rq 0.9
docker exec -it gseq_product_dna python `pwd`/bam2fx.py $bystrand_filename ${output_prefix}.q20.fastq --rq 0.99
docker exec -it gseq_product_dna python `pwd`/bam2fx.py $bystrand_filename ${output_prefix}.q30.fastq --rq 0.999



