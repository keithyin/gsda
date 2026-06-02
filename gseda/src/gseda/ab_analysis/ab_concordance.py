#!/usr/bin/env python3

import argparse
import polars as pl


def main():
    parser = argparse.ArgumentParser(
        description="Compare concordance between control and experiment strategies"
    )
    parser.add_argument("--ctrl", required=True, help="control")
    parser.add_argument("--exp", required=True, help="experiment")
    args = parser.parse_args()

    df_x = pl.read_csv(args.ctrl, separator="\t")
    df_y = pl.read_csv(args.exp, separator="\t")

    # join on channel_id, then filter and compare
    merged = df_x.join(df_y, on="channel_id", how="inner", suffix="_y")

    # filter subreadPasses > 5 (use baseline side)
    merged = merged.filter(pl.col("subreadPasses") > 5)

    total = merged.height
    better = merged.filter(pl.col("concordance") > pl.col("concordance_y"))
    n_better = better.height
    worse = merged.filter(pl.col("concordance") < pl.col("concordance_y"))
    n_worse = worse.height
    n_equal = total - n_better - n_worse

    print(f"Total channels (subreadPasses > 5): {total}")
    print(f"ctrl > exp: {n_better}")
    print(f"Proportion: {n_better / total:.4f}")
    print(f"exp > ctrl: {n_worse}")
    print(f"Proportion: {n_worse / total:.4f}")
    print(f"equal: {n_equal}")
    print(f"Proportion: {n_equal / total:.4f}")
    print()
    print("Channel IDs where control concordance is higher:")
    print(better.select("channel_id").to_series())


if __name__ == "__main__":
    main()
