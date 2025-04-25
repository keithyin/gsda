import pysam
from tqdm import tqdm
from typing import Set, Mapping, List


def extract_interested_channels(bam_filepath: str, min_length: int):
    channels = set()
    with pysam.AlignmentFile(
        filename=bam_filepath, mode="rb", check_sq=False, threads=40
    ) as bam_file:
        for read in tqdm(
            bam_file.fetch(until_eof=True), desc=f"reading {bam_filepath}"
        ):
            if read.query_name.endswith("sbr"):
                continue
            
            if read.query_length >= min_length:
                channels.add(read.get_tag("ch"))
    return channels


def expand_channels(channels: Set[int], expected_num: int) -> Mapping[int, List[int]]:
    result = {}
    for ch in channels:
        result.setdefault(ch, [ch])

    cur_num = len(channels)

    sorted_channels = sorted(list(channels))
    cursor = 0

    for new_ch in range(100000000):
        if cur_num >= expected_num:
            break
        if new_ch in channels:
            continue

        ch_idx = cursor % len(channels)
        cursor += 1
        cur_num += 1

        anchor_ch = sorted_channels[ch_idx]
        result[anchor_ch].append(new_ch)
    return result


def modify_read_info(read_old: pysam.AlignedSegment, old_ch: int, new_ch: int):
    idx = read_old.query_name.rsplit("/", maxsplit=1)[1]
    read_old.query_name = f"read_{new_ch}/{new_ch}/subread/{idx}"
    read_old.set_tag("ch", new_ch, value_type="i")
    return read_old


def dumping_new_subreads(sbr_bam_filepath: str, channels: Mapping[int, List[int]]):
    o_bam_filename = "{}.{}.bam".format(
        sbr_bam_filepath.rsplit(".", maxsplit=1)[0], "expanded-2"
    )
    with pysam.AlignmentFile(
        sbr_bam_filepath, mode="rb", threads=40, check_sq=False
    ) as in_bam:
        with pysam.AlignmentFile(
            o_bam_filename, mode="wb", threads=40, check_sq=False, header=in_bam.header
        ) as out_bam:

            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"dumping {o_bam_filename}"
            ):
                ch = int(record.get_tag("ch"))
                if ch in channels:
                    old_ch = ch
                    for new_ch in channels[ch]:
                        record = modify_read_info(
                            read_old=record, old_ch=old_ch, new_ch=new_ch
                        )
                        out_bam.write(record)
                        old_ch = new_ch


def main_cli():

    sbr_bam_filepath = "/data/ccs_data/speed-test/8k/haorui_hla_subreads.bam"
    smc_bam_filepath = "/data/ccs_data/speed-test/8k//haorui_hla_subreads.smc_all_reads.bam"

    channels = extract_interested_channels(smc_bam_filepath, min_length=8000)
    channels = expand_channels(channels=channels, expected_num=500000)
    dumping_new_subreads(sbr_bam_filepath=sbr_bam_filepath, channels=channels)
    


if __name__ == "__main__":
    main_cli()
