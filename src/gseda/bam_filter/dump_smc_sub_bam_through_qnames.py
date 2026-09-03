import pysam
from tqdm import tqdm
import argparse


def read_qnames(fname: str):
    """从文本文件读取 qname 白名单，每行一个 qname（支持前缀匹配用 `:` 分隔）。"""
    qnames = set()
    with open(fname, mode="r", encoding="utf8") as in_data:
        for idx, line in enumerate(in_data):
            line = line.strip()
            if "rev" in line:
                qname = line.split("rev", maxsplit=1)[0] + "rev"
            elif "fwd" in line:
                qname = line.split("fwd", maxsplit=1)[0] + "fwd"
            else:
                continue                
            qnames.add(qname)
    return qnames


def dump_smc_bam_subset_by_qname(inp_bam_path: str, qnames):
    """按 qname 白名单过滤 BAM，只保留 qname 在名单中的 reads。

    对于 SMC 数据，一个 channel 内的 reads 共享相同的前缀 qname（如 channel id）。
    此脚本通过匹配 qname 来提取目标 channel 的 sub-BAM。
    """
    out_bam_path = inp_bam_path.replace(".bam", ".filtered.bam")
    return _dump_bam_by_qname(inp_bam_path, out_bam_path, qnames, match_key="qname")


def dump_sbr_bam_subset(inp_bam_path: str, channels):
    """按 channel_id 过滤 SBR BAM，只保留 ch tag 在 channels 中的 reads。

    SBR（subreads）BAM 中每条 read 的 ch tag 存储了所属 channel id。
    此脚本通过 ch tag 提取目标 channel 的 sub-BAM。
    """
    out_bam_path = inp_bam_path.replace(".bam", ".filtered.bam")
    return _dump_bam_by_qname(inp_bam_path, out_bam_path, channels, match_key="ch_tag")


def _dump_bam_by_qname(inp_bam_path: str, out_bam_path: str, whitelist, match_key: str):
    """通用的 BAM 过滤入口，支持按 qname 或 ch tag 过滤。

    Args:
        inp_bam_path: 输入 BAM 路径
        out_bam_path: 输出 BAM 路径
        whitelist: 白名单（qname 字符串集合 或 channel_id 整数集合）
        match_key: "qname" 按 record.qname 匹配；"ch_tag" 按 record.get_tag("ch") 匹配
    """
    filtered = 0
    tot = 0
    dumped_channels = set()
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
                channel_id = int(record.get_tag("ch"))
                if match_key == "qname":
                    key = record.qname
                else:
                    key = channel_id
                if key in whitelist:
                    out_bam.write(record)

                    dumped_channels.add(channel_id)

                else:
                    filtered += 1

    print(f"tot:{tot}, filtered:{filtered}, ratio:{filtered/tot}")
    return dumped_channels


def main():
    """
    通过 SMC qname 子采样 BAM 文件，只保留指定 qname 的 reads。

    流程：
    1. 读取 qname_filename 中列出的 qname（每行一个），去重后得到白名单。
    2. 遍历 in_bam 的每一条 alignment record，取 record.qname 与白名单匹配。
    3. 若 qname 在白名单中，则写入 out_bam；否则丢弃。
    4. 完成后打印总记录数、丢弃数和丢弃比例。

    典型用途：SMC 分析中，通过 channel 的 qname（即 channel id 本身作为 qname）
    从全 BAM 中提取目标 channel 的 sub-BAM，用于后续定量或对比分析。

    用法：
        python dump_smc_sub_bam_through_qnames.py input.bam output.bam qname_list.txt
    """

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--smc-bam", type=str)
    parser.add_argument("--sbr-bam", type=str)
    parser.add_argument("--qname-fn", type=str)
    parser.add_argument(
        "--sbr", action="store_true",
        help="按 SBR channel_id (ch tag) 过滤；否则按 SMC qname 过滤。",
    )
    args = parser.parse_args()

    qnames = read_qnames(args.qname_fn)
    print(f"Loaded {len(qnames)} qnames from {args.qname_fn}")
    channels = dump_smc_bam_subset_by_qname(args.smc_bam, qnames)
    dump_sbr_bam_subset(args.sbr_bam, channels)


if __name__ == "__main__":
    main()
