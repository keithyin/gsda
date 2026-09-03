import pysam
from tqdm import tqdm
import argparse


def read_channels(fname: str):
    channels = set()
    with open(fname, mode="r", encoding="utf8") as in_data:
        for idx, line in enumerate(in_data):
            line = line.strip()
            try:
                channel = int(line)
                channels.add(channel)
            except Exception as e:
                print(f"line: {idx}. error {e}")
    return channels


def dump_smc_bam_subset(inp_bam_path: str, out_bam_path: str, channels):
    filtered = 0
    tot = 0
    with pysam.AlignmentFile(
        inp_bam_path, mode="rb", threads=40, check_sq=False
    ) as in_bam:
        with pysam.AlignmentFile(
            out_bam_path, mode="wb", threads=40, check_sq=False, header=in_bam.header
        ) as out_bam:
            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"processing {inp_bam_path}"
            ):
                tot += 1
                ch = int(record.get_tag("ch"))
                if ch in channels:
                    out_bam.write(record)
                else:
                    filtered += 1

    print(f"tot:{tot}, filtered:{filtered}, ratio:{filtered/tot}")


def main(args):
    """
    通过 channel_id 子采样 BAM 文件，只保留指定 channel 的 reads。

    流程：
    1. 读取 channel_filename 中列出的 channel_id（每行一个整数），去重后得到白名单。
    2. 遍历 in_bam 的每一条 alignment record，读取其 "ch" 标签（即 channel_id）。
    3. 若该 channel_id 在白名单中，则写入 out_bam；否则丢弃。
    4. 完成后打印总记录数、丢弃数和丢弃比例。

    典型用途：SMC（Structured Matrix Capture）分析中，从全 BAM 中提取特定 channel 的
    sub-BAM，用于后续定量或对比分析。

    用法：
        python dump_smc_sub_bam_through_channels.py input.bam output.bam channel_list.txt
    """
    # 1. 加载 channel 白名单
    channels = read_channels(args.channel_filename)
    print(f"Loaded {len(channels)} channels from {args.channel_filename}")

    # 2. 按 channel 白名单过滤 BAM，输出 sub-BAM
    dump_smc_bam_subset(args.in_bam, args.out_bam, channels)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("in_bam", type=str)
    parser.add_argument("out_bam", type=str)
    parser.add_argument("channel_filename", type=str)

    main(parser.parse_args())
