from scipy.stats import binomtest

def quick_binom_test(k, n, p_0, alternative='two-sided'):
    """
    使用 scipy 库进行二项分布检验
    
    参数:
    k : int - 成功的次数 (阳性样本数)
    n : int - 总试验次数 (总样本数)
    p_0 : float - 零假设中的预期成功概率 (0 到 1 之间)
    alternative : str - 检验方向: 'two-sided' (双边), 'greater' (单边大于), 'less' (单边小于)
    """
    # 执行检验
    result = binomtest(k, n=n, p=p_0, alternative=alternative)
    
    print(f"=== 二项分布检验报告 ({alternative}) ===")
    print(f"总次数 (n): {n}, 成功次数 (k): {k}")
    print(f"观测样本率: {k/n:.4f} vs 预期概率 (p_0): {p_0}")
    print(f"p-value: {result.pvalue:.6f}")
    
    # 通常与 0.05 的显著性水平进行比较
    alpha = 0.05
    if result.pvalue < alpha:
        print(f"结论: 拒绝零假设 (p < {alpha})，差异具有统计学显著性。")
    else:
        print(f"结论: 无法拒绝零假设 (p >= {alpha})，没有足够证据表明有显著差异。")
        
    return result.pvalue


if __name__ == "__main__":
    # 示例：假设一个新注册页面的历史转化率是 10% (p_0 = 0.1)
    # 换了新设计后，100 个人里有 16 个人转化了 (k = 16, n = 100)
    # 我们想知道转化率是否有显著变化（双边检验）
    quick_binom_test(k=2, n=5, p_0=0.90, alternative='two-sided')