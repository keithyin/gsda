import pysam
import re


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


def dump_subreads_by_channel_id(subreads_bam, channel_ids, output_bam):
    """
    将 subreads.bam 中符合 channel_id 条件的记录 dump 到新的 BAM 文件
    :param subreads_bam: subreads.bam 文件路径
    :param channel_ids: 从 smc.bam 中提取的 channel_id 集合
    :param output_bam: 输出的新的 BAM 文件路径
    """
    with pysam.AlignmentFile(
        subreads_bam, "rb", check_sq=False, threads=40
    ) as subreads, pysam.AlignmentFile(
        output_bam, "wb", header=subreads.header, check_sq=False, threads=40
    ) as out_bam:
        for read in subreads.fetch(until_eof=True):
            # 如果该 read 的 channel_id 在指定的集合中，则将其写入新的 BAM 文件

            if read.has_tag("ch") and read.get_tag("ch") in channel_ids:
                out_bam.write(read)


def main():
    """
    主函数：从 smc.bam 中提取符合条件的 channel_id 并从 subreads.bam 中 dump 相应的记录
    :param smc_bam: smc.bam 文件路径
    :param subreads_bam: subreads.bam 文件路径
    :param output_bam: 输出的新的 BAM 文件路径
    """

    anchor_file = "/data/ccs_data/genegle/Barcode-209.fastq"
    target_bam = "/data/ccs_data/genegle/20250604_240901Y0008_Run0001_adapter.smicing.3.smc_all_reads.bam"
    output_bam = "/data/ccs_data/genegle/20250604_240901Y0008_Run0001_adapter.q20.barcode209.smc_all_reads.bam"

    # 提取 smc.bam 中符合条件的 channel_id
    channel_ids = extract_channel_ids(anchor_file)
    print(f"提取到的 channel_id 数量: {len(channel_ids)}")

    # 从 subreads.bam 中 dump 符合条件的记录到 output_bam
    dump_subreads_by_channel_id(target_bam, channel_ids, output_bam)
    print(f"符合条件的记录已保存至 {output_bam}")


if __name__ == "__main__":
    # 示例：调用主函数

    main()
