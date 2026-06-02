from scipy.stats import fisher_exact

# 2x2 contingency table
table = [
    [4, 396],   # device A
    [35, 665]   # device B
]

odds_ratio, p_value = fisher_exact(table, alternative='two-sided')

print("odds ratio =", odds_ratio)
print("p-value =", p_value)