import pysam
from tqdm import tqdm
import numpy as np


def count_queries(bam_file):
    lengths = []
    count = 0
    # 打开 BAM 文件
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40) as bam:
        # 遍历 BAM 文件中的每条对齐记录
        for read in tqdm(bam.fetch(until_eof=True), desc=f"reading {bam_file}"):
            # 检查 `sc=1`（这里假设 SC 信息存储在 tags 中，若没有需根据实际情况调整）
            if read.has_tag("sc") and read.get_tag("sc") == 1:
                # 检查 query 长度 > 9000
                if read.query_length > 9000:
                    lengths.append(read.query_length)
                    count += 1
    lengths = sorted(lengths)
    lengths = np.array(lengths)

    num = len(lengths)
    print(
        lengths.min(),
        lengths.max(),
        np.median(lengths),
        lengths[int(0.25 * num)],
        lengths[int(0.50 * num)],
        lengths[int(0.75 * num)],
        lengths[int(0.95 * num)],
        lengths[int(0.99 * num)],
        lengths[int(0.999 * num)],
    )

    return count


if __name__ == "__main__":

    # 示例：调用函数并输出结果
    bam_file = "/data/ccs_data/speed-test/13k-5h/4cc/20241220_Sync_Y0701_01_H01_Run0001_called.adapter.smc_all_reads.bam"
    result = count_queries(bam_file)
    print(f"符合条件的 query 个数: {result}")
