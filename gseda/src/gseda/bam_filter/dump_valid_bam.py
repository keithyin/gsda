import pysam
import re
import argparse
from multiprocessing import cpu_count
from tqdm import tqdm


def extract_channel_ids(fq_file):
    """
    提取 smc.bam 中满足条件的 channel_id 信息
    :param smc_bam: smc.bam 文件路径
    :return: 一个包含符合条件的 channel_id 的集合
    """
    channel_ids = set()
    # 20250604_240901Y0008_Run0001/41790_0_3753

    pattern = r"/(\d+)_"  # 匹配斜杠后的数字，直到下划线

    with pysam.FastxFile(fq_file) as records:
        for record in records:
            channel = re.search(pattern, record.name).group(1)
            channel_ids.add(int(channel))

    return channel_ids


def extract_channel_from_bam(bam_file):
    """
    提取 smc.bam 中满足条件的 channel_id 信息
    :param smc_bam: smc.bam 文件路径
    :return: 一个包含符合条件的 channel_id 的集合
    """
    channel_ids = set()
    with pysam.AlignmentFile(
        bam_file, "rb", check_sq=False, threads=cpu_count()
    ) as in_bam:
        for record in tqdm(in_bam.fetch(until_eof=True), desc=f"reading {bam_file}"):
            channel_id = int(record.get_tag("ch"))
            channel_ids.add(channel_id)
    return channel_ids


def dump_reads_by_channel_id(in_bam, channel_ids, output_bam):
    """
    将 subreads.bam 中符合 channel_id 条件的记录 dump 到新的 BAM 文件
    :param subreads_bam: subreads.bam 文件路径
    :param channel_ids: 从 smc.bam 中提取的 channel_id 集合
    :param output_bam: 输出的新的 BAM 文件路径
    """
    with pysam.AlignmentFile(
        in_bam, "rb", check_sq=False, threads=cpu_count() // 2
    ) as in_bam_reader, pysam.AlignmentFile(
        output_bam, "wb", header=in_bam_reader.header, check_sq=False, threads=cpu_count() // 2
    ) as out_bam_writer:
        for read in tqdm(in_bam_reader.fetch(until_eof=True), desc=f"dumping to {output_bam}"):
            # 如果该 read 的 channel_id 在指定的集合中，则将其写入新的 BAM 文件
            ch = None
            if read.has_tag("ch"):
                ch = int(read.get_tag("ch"))
            else:
                ch = int(read.query_name.split("_")[1])
            if ch in channel_ids:
                out_bam_writer.write(read)


def main():
    """
    主函数：从 smc.bam 中提取符合条件的 channel_id 并从 subreads.bam 中 dump 相应的记录
    :param smc_bam: smc.bam 文件路径
    :param subreads_bam: subreads.bam 文件路径
    :param output_bam: 输出的新的 BAM 文件路径
    """

    parser = argparse.ArgumentParser(prog="dump valid bam")
    parser.add_argument("inp_bam", type=str)
    parser.add_argument("oup_bam", type=str)
    parser.add_argument("anchor_bam", type=str)
    args = parser.parse_args()

    # 提取 smc.bam 中符合条件的 channel_id
    channel_ids = extract_channel_from_bam(args.anchor_bam)

    print(f"提取到的 channel_id 数量: {len(channel_ids)}")

    # 从 subreads.bam 中 dump 符合条件的记录到 output_bam
    dump_reads_by_channel_id(args.inp_bam, channel_ids, args.oup_bam)
    print(f"符合条件的记录已保存至 {args.oup_bam}")


if __name__ == "__main__":
    # 示例：调用主函数

    main()
