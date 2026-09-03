import pysam


def extract_channel_ids(smc_bam):
    """
    提取 smc.bam 中满足条件的 channel_id 信息
    :param smc_bam: smc.bam 文件路径
    :return: 一个包含符合条件的 channel_id 的集合
    """
    channel_ids = set()

    with pysam.AlignmentFile(smc_bam, "rb", check_sq=False) as smc:
        for read in smc.fetch(until_eof=True):
            # 检查是否有 'sc' 标签且其值为 1，并且查询长度 > 9000
            if (
                read.has_tag("sc")
                and read.get_tag("sc") == 1
                and read.query_length > 9000
            ):
                # 提取 channel_id 标签（假设 channel_id 存储在 'channel_id' 标签中）
                if read.has_tag("channel_id"):
                    channel_ids.add(read.get_tag("channel_id"))

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
            if read.has_tag("channel_id") and read.get_tag("channel_id") in channel_ids:
                out_bam.write(read)


def main(smc_bam, subreads_bam, output_bam):
    """
    主函数：从 smc.bam 中提取符合条件的 channel_id，并从 subreads.bam 中 dump 相应的记录
    :param smc_bam: smc.bam 文件路径
    :param subreads_bam: subreads.bam 文件路径
    :param output_bam: 输出的新的 BAM 文件路径
    """
    # 提取 smc.bam 中符合条件的 channel_id
    channel_ids = extract_channel_ids(smc_bam)
    print(f"提取到的 channel_id 数量: {len(channel_ids)}")

    # 从 subreads.bam 中 dump 符合条件的记录到 output_bam
    dump_subreads_by_channel_id(subreads_bam, channel_ids, output_bam)
    print(f"符合条件的记录已保存至 {output_bam}")


if __name__ == "__main__":
    # 示例：调用主函数
    smc_bam = "/data/ccs_data/speed-test/13k-5h/4cc/20241220_Sync_Y0701_01_H01_Run0001_called.adapter.smc_all_reads.bam"
    subreads_bam = "/data/ccs_data/speed-test/13k-5h/4cc/20241220_Sync_Y0701_01_H01_Run0001_called.adapter.bam"
    output_bam = "/data/ccs_data/speed-test/13k-5h/4cc/20241220_Sync_Y0701_01_H01_Run0001_called.adapter.filtered.bam"

    main(smc_bam, subreads_bam, output_bam)
