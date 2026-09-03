
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from collections import defaultdict
import numpy as np

# 1. 打开 BAM 文件
smc_bam_path = "/data/ccs_data/subreads-dist-smc-dist-discrepency/20250402_240601Y0012_Run0002.smc_all_reads.bam"  # ← 修改为你的 bam 路径
smc_bamfile = pysam.AlignmentFile(smc_bam_path, "rb", check_sq=False, threads=40)

# 2. 提取每个 read 的长度
channel_ids = set([read.get_tag("ch") for read in tqdm(smc_bamfile, desc=f"reading {smc_bamfile}")])



bam_path = "/data/ccs_data/subreads-dist-smc-dist-discrepency/adapter_bam/20250402_240601Y0012_Run0002_adapter.bam"  # ← 修改为你的 bam 路径
bamfile = pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=40)

channel_lengths = defaultdict(list)

for read in tqdm(bamfile.fetch(until_eof=True), desc=f"reading {bamfile}"):
    
    try:
        channel = read.get_tag("ch")  # channel ID
        if channel not  in channel_ids:
            continue
        cx = int(read.get_tag("cx"))       # cx tag
    except KeyError:
        continue  # skip reads without 'ch' or 'cx'

    if cx == 3:
        channel_lengths[channel].append(read.query_length)

bamfile.close()

# 每个 channel 取中位数
channel_medians = [
    np.median(lengths)
    for lengths in channel_lengths.values()
    if lengths  # 确保不是空的
]
channel_medians = np.array(channel_medians)
max_rl = channel_medians.max()
bins = np.arange(0, max_rl + 10 -1 , 10)


# 绘图
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.histplot(channel_medians, bins=bins, kde=True, color="orange")
plt.xlim(70, 15000)
plt.title("Channel-level Read Length Distribution (cx == 3, median)")
plt.xlabel("Median Read Length per Channel")
plt.ylabel("Channel Count")
plt.tight_layout()
plt.savefig("channel_read_length_distribution.png", dpi=300)