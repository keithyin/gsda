#!/usr/bin/env python3
"""
Convert subreads.bam from internal format to PacBio-like format.
Focus on modifying header (@HD, @RG, @PG) and ensuring proper RG tag on reads.
"""

import os
import sys
import argparse
import pysam
from tqdm import tqdm
import multiprocessing as mp
import array


def get_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="将内部格式的 subreads.bam 转换为类 PacBio 格式。"
    )
    parser.add_argument("input_bam", help="输入的内部格式 BAM 文件")
    # parser.add_argument("reference_bam", help="作为 Header 模板的 PacBio 参考 BAM 文件")
    parser.add_argument("output_bam", help="输出的 BAM 文件路径")

    # 如果参数为空，显示帮助信息
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def pb_header():
    hd = {'VN': '1.5', 'SO': 'unknown', 'pb': '3.0.7'}
    rg = [
        {'ID': 'ffffffff',
         'PL': 'PACBIO',
         'DS': 'READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-789-500;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000',
         'PU': 'S',
         'PM': 'SEQUELII'
         }
    ]
    pg = [{'ID': 'samtools', 'PN': 'samtools',
           'VN': '1.21', 'CL': 'samtools view -S -b sd.sam'}]

    return hd, rg, pg


# def extract_pacbio_header_info(reference_bam):
#     """从 PacBio 参考 BAM 中提取 HD, RG, PG 字典信息"""
#     with pysam.AlignmentFile(reference_bam, "rb", check_sq=False) as ref:
#         header = ref.header.to_dict()
#         # 注意：to_dict() 返回的 HD 通常是字典，而 RG 和 PG 是列表
#         hd = header.get('HD', {})
#         rg = header.get('RG', [])
#         pg = header.get('PG', [])
#     return hd, rg, pg


def build_new_header(input_bam, pacbio_hd, pacbio_rg, pacbio_pg):
    """基于输入文件和参考信息构建新的 Header 对象"""
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as inp:
        header_dict = inp.header.to_dict()

        # 替换核心字段
        if pacbio_hd:
            header_dict['HD'] = pacbio_hd
        if pacbio_rg:
            header_dict['RG'] = pacbio_rg
        if pacbio_pg:
            header_dict['PG'] = pacbio_pg

        # 使用 pysam 从字典重新构建 Header 对象，这比手动拼接字符串更安全
        return pysam.AlignmentHeader.from_dict(header_dict)


def convert_bam(input_bam, output_bam):
    """执行转换主逻辑"""
    # print(f"[*] 正在从 {reference_bam} 提取 PacBio Header 信息...")
    hd_info, rg_list, pg_list = pb_header()

    # 获取用于 Read Tag 的 RG ID
    # 通常取参考文件中的第一个 RG ID
    rg_id = "ffffffff"
    if rg_list and 'ID' in rg_list[0]:
        rg_id = rg_list[0]['ID']

    print(f"[*] 目标 RG ID 设为: {rg_id}")

    print(f"[*] 正在构建新 Header 并写入 {output_bam}...")
    new_header = build_new_header(input_bam, hd_info, rg_list, pg_list)

    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, threads=mp.cpu_count() // 2) as infile, \
            pysam.AlignmentFile(output_bam, "wb", header=new_header, threads=mp.cpu_count() // 2) as outfile:

        for read in tqdm(infile, desc=f"processing..."):
            # 更新或设置 RG 标签 (Z 代表字符串类型)
            read.set_tag('RG', rg_id, value_type='Z')

            read.set_tag("ar", None)
            read.set_tag("dw", None)

            read.set_tag("ip", array.array('B', [9]*read.query_length))
            read.set_tag("pw", array.array('B', [9]*read.query_length))

            outfile.write(read)

    print(f"[OK] 转换完成")


def main():
    args = get_args()

    # 检查文件是否存在
    if not os.path.exists(args.input_bam):
        print(f"Error: 输入文件不存在 -> {args.input_bam}")
        sys.exit(1)

    convert_bam(args.input_bam, args.output_bam)


if __name__ == "__main__":
    our_sbr = "/data1/ccs_data/202603-good-ecoli-data/20260311_240601Y0012_Run0003_adapter.bam"
    pb_sbr = "/data1/simulated_data/12k_7_pass.bam"

    # print(extract_pacbio_header_info(pb_sbr))
    main()
