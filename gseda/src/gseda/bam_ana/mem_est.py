import pysam
from tqdm import tqdm
import numpy as np
from typing import List
import argparse
from numba import njit


@njit
def sliding_topn_square_sum_numba(arr: np.ndarray, k: int, topn: int) -> np.ndarray:
    n = len(arr)
    output = np.empty(n - k + 1)

    for i in range(n - k + 1):
        window = arr[i : i + k]
        # 获取 topN（不排序），使用 np.partition
        idx = np.argpartition(window, -topn)[-topn:]
        square_sum = 0.0
        for j in range(topn):
            val = window[idx[j]]
            square_sum += val * val
        output[i] = square_sum

    return output


class ChannelLength:
    def __init__(self, ch):
        self.ch = ch
        self.lengths = []

    def median(self):
        return float(np.median(np.array(self.lengths)))


def read_channel_lengths(bam_file: str):
    result = []

    with pysam.AlignmentFile(
        filename=bam_file, mode="rb", check_sq=False, threads=40
    ) as reader:
        for record in tqdm(reader.fetch(until_eof=True), desc=f"reading {bam_file}"):
            ch = record.get_tag("ch")
            cx = record.get_tag("cx")
            if cx < 3:
                continue

            if len(result) == 0:
                result.append(ChannelLength(ch=ch))

            if result[-1].ch != ch:
                result.append(ChannelLength(ch=ch))

            assert result[-1].ch == ch, "{}!={}".format(result[-1].ch, ch)
            result[-1].lengths.append(record.query_length)

    return result


def mem_est(
    channel_lengths: List[ChannelLength], min_passes: int, window: int, top_n: int
):
    result = [r for r in channel_lengths if len(r.lengths) >= min_passes]
    medians = [r.median() for r in result]
    medians = np.array(medians)

    mems = sliding_topn_square_sum_numba(medians, window, top_n)
    # 4 for int. 2 for each strand
    print("PeakUsed EST. {:.2f}GB".format(np.max(mems) * 4 * 2 / (1024 * 1024 * 1024)))

def main_cli():
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("sbr_bam", type=str)
    parser.add_argument("--min-passes", required=True, type=int, dest="min_passes")
    parser.add_argument("--win", required=True, type=int)
    parser.add_argument("--topn", required=True, type=int)
    
    args = parser.parse_args()
    channel_lengths = read_channel_lengths(bam_file=args.sbr_bam)
    mem_est(channel_lengths=channel_lengths, min_passes=args.min_passes, window=args.win, top_n=args.topn)
        
    pass


if __name__ == "__main__":
    main_cli()
