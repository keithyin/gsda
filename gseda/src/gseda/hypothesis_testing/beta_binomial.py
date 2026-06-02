
import scipy.stats as stats


def beta_binomial_test(
    k: int,
    n: int,
    alpha: float,
    beta: float,
):
    """
    Beta-binomial upper-tail test

    H0:
        variant frequency comes from background error distribution

    Parameters
    ----------
    k : int
        observed variant reads

    n : int
        total reads

    alpha, beta :
        beta prior parameters

    Returns
    -------
    pvalue
    """

    # P(X >= k)
    pvalue = stats.betabinom.sf(k - 1, n, alpha, beta)

    return pvalue

def main():
    k = 35
    n = 700
    alpha = 4
    beta = 400
    pv = beta_binomial_test(k=k, n=n, alpha=alpha, beta=beta)
    print(f"p_value:{pv}")
    

if __name__ == "__main__":
    main()