from collections import Counter
from typing import List
import json
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

def load_data(file_path: str):
    with open(file_path, "r") as file:
        for line in tqdm(file):
            line = line.strip()
            v = json.loads(line)
            qual = v["msa_seqs"][0]
            assert isinstance(qual, str)
            qual = qual.replace("-", "9")
            qual = qual.replace("#", "9")
            sbrs = v["msa_seqs"][3:]
            yield qual, sbrs

def count_alternating_columns(q_values: str, surbeads: List[str]) -> int:
    # 将 Q 值字符串转换为整数列表
    q_values_list = list(map(int, q_values.strip()))

    # 存储满足条件的列数
    alternating_count = 0

    # 遍历每一列
    for col_index in range(len(q_values_list)):
        if q_values_list[col_index] < 4:  # 判断 Q 值是否小于 4
            # 提取当前列的碱基
            column_bases = [surbead[col_index] for surbead in surbeads if col_index < len(surbead)]

            # 计算碱基的频次
            base_counts = Counter(column_bases)
            total_count = sum(base_counts.values())

            # 如果出现的碱基数小于 2，不可能交替
            if len(base_counts) < 2:
                continue

            # 获取出现频次最高的两个碱基及其比例
            most_common = base_counts.most_common(2)
            top1_base, top1_count = most_common[0]
            top2_base, top2_count = most_common[1]

            top1_ratio = top1_count / total_count
            top2_ratio = top2_count / total_count

            # 检查两者的比例是否都在 40% 以上
            if top1_ratio >= 0.45 and top2_ratio >= 0.45:
                alternating_count += 1

    return alternating_count

def main():
    file_path = "/data/ccs_data/case-study/20250310-lowQ30/sbr2icing.q30.asrtc.txt"
    locus = 0
    channels = 0
    
    in_consis_cnts = []
    
    for q_values, surbeads in load_data(file_path):
        result = count_alternating_columns(q_values, surbeads)
        if result > 0 and result < 50:
            in_consis_cnts.append(min(result, 100))
        locus += result
        if result > 0:
            channels += 1
    print(locus, channels)
    plt.figure(figsize=(10, 6))
    sns.histplot(in_consis_cnts, bins=50, kde=False)
    plt.title('Inconsistent Counts Distribution')
    plt.xlabel('Inconsistent Count')
    plt.ylabel('Frequency')
    plt.xticks(range(0, 50, 1))
    plt.xticks(rotation=45)

    plt.savefig('in_consis_cnts.png')
    plt.close()
            

if __name__ == "__main__":
    main()