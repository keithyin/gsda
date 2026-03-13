import pysam
from tqdm import tqdm
from typing import Set, Mapping, List
import multiprocessing as mp


def extract_interested_channels_from_smc_bam(bam_filepath: str, min_length: int, max_length=None, np_lower_bound=3):
    channels = set()
    with pysam.AlignmentFile(
        filename=bam_filepath, mode="rb", check_sq=False, threads=mp.cpu_count()
    ) as bam_file:
        for read in tqdm(
            bam_file.fetch(until_eof=True), desc=f"reading {bam_filepath}"
        ):
            if read.query_name.endswith("sbr"):
                continue

            if read.get_tag("np") < np_lower_bound:
                continue
            if read.query_length < min_length:
                continue
            if max_length is not None and read.query_length > max_length:
                continue
            channels.add(read.get_tag("ch"))
    return channels


def extract_interested_channels_from_called_bam(bam_filepath: str, min_length: int, max_length=None):
    channels = set()
    with pysam.AlignmentFile(
        filename=bam_filepath, mode="rb", check_sq=False, threads=mp.cpu_count()
    ) as bam_file:
        for read in tqdm(
            bam_file.fetch(until_eof=True), desc=f"reading {bam_filepath}"
        ):
            if read.query_length < min_length:
                continue
            if max_length is not None and read.query_length > max_length:
                continue
            channels.add(int(read.query_name.split("_", maxsplit=1)[1]))
    return channels


def expand_channels(channels: Set[int], expected_num: int) -> Mapping[int, List[int]]:
    channels = sorted(list(channels))
    num_ch = len(channels)
    result = {}
    for new_ch in range(expected_num):
        old_ch_idx = new_ch % num_ch
        old_ch = channels[old_ch_idx]
        result.setdefault(old_ch, []).append(new_ch)
    return result


def modify_read_info(read_old: pysam.AlignedSegment, old_ch: int, new_ch: int):
    idx = read_old.query_name.rsplit("/", maxsplit=1)[1]
    read_old.query_name = f"read_{new_ch}/{new_ch}/subread/{idx}"
    read_old.set_tag("ch", new_ch, value_type="i")
    return read_old


def dumping_new_subreads(sbr_bam_filepath: str, channels: Mapping[int, List[int]], min_length, max_length=None):
    o_bam_filename = "{}.{}.{}-{}.bam".format(
        sbr_bam_filepath.rsplit(".", maxsplit=1)[
            0], "expanded-2", min_length, max_length
    )
    o_bam_sorted_filename = "{}.{}.{}-{}.sorted.bam".format(
        sbr_bam_filepath.rsplit(".", maxsplit=1)[
            0], "expanded-2", min_length, max_length
    )
    with pysam.AlignmentFile(
        sbr_bam_filepath, mode="rb", threads=mp.cpu_count() // 2, check_sq=False
    ) as in_bam:
        with pysam.AlignmentFile(
            o_bam_filename, mode="wb", threads=mp.cpu_count() // 2, check_sq=False, header=in_bam.header
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
    print("sorting....")
    pysam.sort("-n", "-t", "ch", "-@", f"{mp.cpu_count()}", "-o",
               o_bam_sorted_filename, o_bam_filename)


def main_cli():
    """
    该脚本 生成 x-Hour y-InsertSize 的数据, 用于后续对于 SMICING 的速度评估
    Input:
        call.bam: 用来识别出 哪些 channel 测序时长达到了目标。是基于 called.bam 中 reads 的长度来近似计算得出
        smc.bam : 用来识别出 哪些 channel 的 Insert-Size 达到了目标。是基于 smc.bam 中 reads 的长度估算的
        sbr.bam : 将满足 条件的 channels 的 dump 出来
    """
    sequencing_speed = 4  # number of bases per secs

    hour_lower_bound = 2
    hour_upper_bound = 4

    insert_size_lower_bound = 2500
    insert_size_upper_bound = 4500

    np_lower_bound = 3

    reads_len_lower_bound = hour_lower_bound * 60 * 60 * sequencing_speed
    reads_len_upper_bound = hour_upper_bound * 60 * 60 * sequencing_speed

    called_bam_filepath = "/data1/ccs_data/20260127-plasmid-speedup/20260127_250302Y0001_Run0001/20260127_250302Y0001_Run0001_called.bam"
    sbr_bam_filepath = "/data1/ccs_data/20260127-plasmid-speedup/20260127_250302Y0001_Run0001/20260127_250302Y0001_Run0001_adapter.bam"
    smc_bam_filepath = "/data1/ccs_data/20260127-plasmid-speedup/20260127_250302Y0001_Run0001/20260127_250302Y0001_Run0001.smc_all_reads.bam"

    smc_channels = extract_interested_channels_from_smc_bam(
        smc_bam_filepath, min_length=insert_size_lower_bound, max_length=insert_size_upper_bound, np_lower_bound=np_lower_bound)
    print(f"SmcChannels:{len(smc_channels)}")

    called_channels = extract_interested_channels_from_called_bam(
        called_bam_filepath, min_length=reads_len_lower_bound, max_length=reads_len_upper_bound)

    print(f"CalledChannels:{len(called_channels)}")

    channels = smc_channels.intersection(called_channels)
    print(f"JoinedChannels:{len(channels)}")

    smc_channels = expand_channels(channels=channels, expected_num=500000)
    dumping_new_subreads(sbr_bam_filepath=sbr_bam_filepath,
                         channels=smc_channels, min_length=insert_size_lower_bound, max_length=insert_size_upper_bound)


if __name__ == "__main__":
    main_cli()
