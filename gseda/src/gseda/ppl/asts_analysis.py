import subprocess
import argparse
import pathlib
import pysam
from multiprocessing import cpu_count
from tqdm import tqdm
import math
import matplotlib.pyplot as plt
import numpy as np

def asts_bam_identity(asts_bam_file: str):
    q_count = {}
    with pysam.AlignmentFile(asts_bam_file, mode="rb", check_sq=False, threads=cpu_count()) as in_bam:
        for record in tqdm(in_bam.fetch(), desc=f"reading {asts_bam_file}"):
            identity = record.get_tag("iy")
            identity = min(0.9999, identity)
            identity = max(0.0001, identity)
            
            q = int(round((-10 * math.log10(1 - identity))))
            q_count.setdefault(q, 0)
            q_count[q] += 1
    return q_count
    

def plot(q_count, dest_filepath):
    # 创建子图
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10))
    
    all_values = list(range(0, 41))
    total_count = sum( [v for _, v in q_count.items()] )
    frequencies = [q_count.get(value, 0) for value in all_values]
    probabilities = [freq/total_count for freq in frequencies]
    # 频数直方图
    bars1 = ax1.bar(all_values, frequencies, color='skyblue', 
                    edgecolor='navy', alpha=0.7, width=0.8)
    ax1.set_xlabel('(0-40)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Count distribution', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_xticks(range(0, 41, 1))  # 每5个刻度显示一个
    ax1.set_xlim(-0.5, 40.5)

    # 在频数柱子上方显示数值（只显示非零值）
    for bar, freq in zip(bars1, frequencies):
        if freq > 0:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + max(frequencies)*0.01,
                    str(int(freq)), ha='center', va='bottom', fontsize=8)

    # 概率密度直方图
    bars2 = ax2.bar(all_values, probabilities, color='lightcoral', 
                    edgecolor='darkred', alpha=0.7, width=0.8)
    ax2.set_xlabel('(0-40)', fontsize=12)
    ax2.set_ylabel('Prob', fontsize=12)
    ax2.set_title('Prob', fontsize=14, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    ax2.set_xticks(range(0, 41, 1))  # 每5个刻度显示一个
    ax2.set_xlim(-0.5, 40.5)

    # 在概率柱子上方显示数值（只显示非零值）
    for bar, prob in zip(bars2, probabilities):
        if prob > 0:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + max(probabilities)*0.01,
                    f'{prob:.3f}', ha='center', va='bottom', fontsize=8)

    # # 添加统计信息文本框
    # stats_text = f'''统计信息:
    # 总数据点: {total_count}
    # 出现过的数值: {sum(1 for f in frequencies if f > 0)}
    # 最频繁值: {all_values[np.argmax(frequencies)]} (出现{max(frequencies)}次, {max(probabilities):.3f})
    # 平均值: {np.mean(data):.2f}
    # 标准差: {np.std(data):.2f}'''

    # ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
    #         verticalalignment='top', fontsize=10,
    #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    # plt.show()
    plt.savefig(dest_filepath)

    # 可选：打印详细的统计信息
    print("详细统计信息:")
    print(f"数据范围: 0-40")
    print(f"总样本数: {total_count}")
    print(f"覆盖的数值个数: {sum(1 for f in frequencies if f > 0)}")
    # print(f"缺失的数值: {[i for i in range(0, 41) if i not in counts]}")
    print(f"最频繁的值: {all_values[np.argmax(frequencies)]} (频率: {max(probabilities):.4f})")


def main():
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("sbr_bam")
    parser.add_argument("smc_bam")
    parser.add_argument("--rq-range", type=str, default="0.95:1.0", dest="rq_range", help="start:end")
    
    
    args = parser.parse_args()
    
    
    subreads_bam = args.sbr_bam
    smc_bam = args.smc_bam
    
    smc_path = pathlib.Path(smc_bam)
    prefix = smc_path.parent.joinpath(f"{smc_path.stem}.asts")
    
    asts_cmd = f"asts -q {subreads_bam} -t {smc_bam} -p {prefix} --rq-range {args.rq_range}"
    print(f"running: {asts_cmd}")
    
    subprocess.check_call(asts_cmd, shell=True)
    q_count = asts_bam_identity(f"{prefix}.bam")
    dest_img_path = f"{prefix}.q_dist.jpg"
    plot(q_count=q_count, dest_filepath=dest_img_path)
    print(f"dest_image_path: {dest_img_path}")
    
    pass

if __name__ == "__main__":
    main()