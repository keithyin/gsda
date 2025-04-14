
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from collections import Counter
import numpy as np

# 1. 打开 BAM 文件
bam_path = "/data/ccs_data/subreads-dist-smc-dist-discrepency/20250402_240601Y0012_Run0002.smc_all_reads.bam"  # ← 修改为你的 bam 路径
bamfile = pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=40)

# 2. 提取每个 read 的长度
read_lengths = [read.query_length for read in tqdm(bamfile.fetch(until_eof=True), desc=f"reading {bam_path}")]

read_lengths = np.array(read_lengths)
max_rl = read_lengths.max()

bins = np.arange(0, max_rl + 10 - 1, 10)


# _, bins = np.histogram(read_lengths, bins="auto")
# print("bins:", bins, len(bins))

# print("diffs:", np.diff(bins))
counts = Counter(read_lengths)
counts = sorted(counts.items(), key=lambda x: x[0])
for count in counts:
    if count[0] > 100:
        break
    print(count[0], count[1])
    
# # print(counts)


# print(len(read_lengths))
bamfile.close()

# 3. 可视化分布
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.histplot(read_lengths, bins=bins, kde=True, color="skyblue")
plt.title("Read Length Distribution")
plt.xlabel("Read Length")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("output.png")