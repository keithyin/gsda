#!/usr/bin/env python3
"""
合并多个 BAM 文件为一个 BAM 文件
"""

import pysam
import argparse
import os
import sys
from datetime import datetime
from tqdm import tqdm


def merge_bam_files(input_bams, output_bam):
    """
    合并多个 BAM 文件

    参数:
        input_bams: 输入 BAM 文件路径列表
        output_bam: 输出 BAM 文件路径
        merge_method: 合并方法 ('pysam' 或 'samtools')
    """

    if not input_bams:
        print("错误: 没有提供输入 BAM 文件")
        sys.exit(1)

    # 检查输入文件是否存在
    for bam_file in input_bams:
        if not os.path.exists(bam_file):
            print(f"错误: 文件不存在: {bam_file}")
            sys.exit(1)
        if not os.path.exists(bam_file + '.bai'):
            print(f"警告: 未找到索引文件: {bam_file}.bai")

    print(f"开始合并 {len(input_bams)} 个 BAM 文件...")
    print(f"输入文件: {input_bams}")
    print(f"输出文件: {output_bam}")

    start_time = datetime.now()

    # 方法1: 使用 pysam 合并
    merge_with_pysam(input_bams, output_bam)

    end_time = datetime.now()
    duration = end_time - start_time

    print(f"合并完成! 耗时: {duration}")
    print(f"输出文件: {output_bam}")


def merge_with_pysam(input_bams, output_bam):
    """使用 pysam 合并 BAM 文件"""

    # 打开第一个文件作为模板
    with pysam.AlignmentFile(input_bams[0], "rb") as template:
        header = template.header.to_dict()

    # 创建输出文件
    with pysam.AlignmentFile(output_bam, "wb", header=header, check_sq=False, threads=40) as out_bam:
        for bam_file in input_bams:
            print(f"正在处理: {bam_file}")
            with pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40) as in_bam:
                for read in tqdm(in_bam, desc=f"dumping {bam_file}"):
                    out_bam.write(read)


def main():
    parser = argparse.ArgumentParser(
        description='合并多个 BAM 文件为一个 BAM 文件',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python merge_bam.py -i sample1.bam sample2.bam sample3.bam -o merged.bam
  python merge_bam.py -i *.bam -o merged.bam --method samtools
        '''
    )

    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='输入 BAM 文件路径（支持通配符）')
    parser.add_argument('-o', '--output', required=True,
                        help='输出 BAM 文件路径')

    args = parser.parse_args()

    # 处理通配符
    import glob
    input_files = []
    for pattern in args.input:
        input_files.extend(glob.glob(pattern))

    input_files = list(set(input_files))  # 去重

    if not input_files:
        print("错误: 未找到匹配的 BAM 文件")
        sys.exit(1)

    merge_bam_files(input_files, args.output)


if __name__ == '__main__':
    main()
