import pysam
from tqdm import tqdm
import numpy as np
from typing import List
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib


class ChannelLength:
    def __init__(self, ch):
        self.ch = ch
        self.lengths = []

    def median(self):
        return np.median(np.array(self.lengths))


def read_bam_info(bam_file: str):
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
                # if len(result) > 10000:
                #     break
                result.append(ChannelLength(ch=ch))

            assert result[-1].ch == ch, "{}!={}".format(result[-1].ch, ch)
            result[-1].lengths.append(record.query_length)

    return result


def plot(result: List[ChannelLength], o_path: str):
    result = [r for r in result if len(r.lengths) >= 3]
    medians = [r.median() for r in result]
    # print(medians)
    x = list(range(len(medians)))

    # 设置 seaborn 风格
    sns.set(style="whitegrid")

    # 创建折线图
    plt.figure(figsize=(40, 21))
    sns.lineplot(x=x, y=medians, marker="o")

    plt.xlabel("Channel Index")
    plt.ylabel("Median Read Length")
    plt.title("Median Read Length per Channel")

    # 保存图片
    plt.tight_layout()
    plt.savefig(o_path)
    print(f"Saved plot to {o_path}")
    plt.close()
    pass


def main_cli():
    bam = "/data/ccs_data/speed-test/4.6k5h/20250424_240601Y0005_Run0001_adapter.bam"
    o_path = "{}.incomming_lengths.png".format(pathlib.Path(bam).stem)
    incomming_lengths = read_bam_info(bam)
    plot(incomming_lengths, o_path=o_path)
    
    pass

if __name__ == "__main__":
    main_cli()