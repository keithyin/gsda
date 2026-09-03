import pysam

if __name__ == "__main__":
    input_bam = "/data/ccs_data/data2025Q1/S_aureus_1h/20250207_Sync_Y0002_02_H01_Run0001_called.adapter.bam"
    output_bam = "/data/ccs_data/data2025Q1/S_aureus_1h/20250207_Sync_Y0002_02_H01_Run0001_called.adapter.subset_2000_channels.bam"

    bam_in = pysam.AlignmentFile(input_bam, "rb", threads=40, check_sq=False)
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header, threads=40, check_sq=False)

    # 用于记录前 1000 个 channel
    channel_set = set()
    max_channels = 2000

    # channel -> list of reads 的 map 也可以加快效率，但不需要时用 set 更省内存
    for read in bam_in:
        try:
            ch = read.get_tag("ch")
        except KeyError:
            continue  # 没有 ch tag 的跳过

        if ch in channel_set or len(channel_set) < max_channels:
            channel_set.add(ch)
            bam_out.write(read)

        # 仅当已收满 1000 个 channel 且当前 read 的 channel 不在其中，才跳过
        elif ch not in channel_set:
            continue

    bam_in.close()
    bam_out.close()
    print(f"✔ 已输出包含前 {max_channels} 个 channel 的 BAM 文件：{output_bam}")