import numpy as np
from scipy.stats import permutation_test
from tqdm import tqdm

# ---------- 原始数据 ----------
nA, nB = 400, 700
clicksA, clicksB = 4 * 83, 7 * 85
total_clicks = clicksA + clicksB
total_no_clicks = nA + nB - total_clicks

# 构造二值数组：1=点击，0=未点击
data = np.array([1]*total_clicks + [0]*total_no_clicks)
# 分组标签：前500个为A，后1000个为B
groups = np.array(['A']*nA + ['B']*nB)

# 原始观察到的点击率差（B-A）
pA_obs = clicksA / nA
pB_obs = clicksB / nB
obs_diff = pB_obs - pA_obs
print(f"Observed difference (B-A): {obs_diff:.4f}")

# ---------- 手写置换检验 ----------
np.random.seed(42)
K = 10000
perm_diffs = np.empty(K)

for i in tqdm(range(K), desc="running permutation test"):
    # 随机打乱点击/未点击的序列
    shuffled = np.random.permutation(data)
    # 重新划分伪A和伪B
    perm_pA = np.mean(shuffled[:nA])
    perm_pB = np.mean(shuffled[nA:])
    perm_diffs[i] = perm_pB - perm_pA

# 单侧p值（B > A）
p_value_one_sided = (np.sum(perm_diffs >= obs_diff) + 1) / (K + 1)
# 双侧p值
p_value_two_sided = (np.sum(np.abs(perm_diffs) >= obs_diff) + 1) / (K + 1)

print(f"手写置换检验 - 单侧p值: {p_value_one_sided:.5f}")
print(f"手写置换检验 - 双侧p值: {p_value_two_sided:.5f}")

# ---------- 使用 scipy （更简洁）----------
# 定义统计函数：计算两组点击率的差值（B - A）
def diff_proportion(x, y):
    return np.mean(y) - np.mean(x)

# 准备两个组的原始结果（0/1列表）
A_data = [1]*clicksA + [0]*(nA - clicksA)
B_data = [1]*clicksB + [0]*(nB - clicksB)

res = permutation_test((A_data, B_data), diff_proportion, 
                       n_resamples=10000, alternative='greater')  # 'greater' 表示 B > A
print(f"scipy 单侧p值: {res.pvalue:.5f}")