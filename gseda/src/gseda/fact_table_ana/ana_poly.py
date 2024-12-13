import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)
sys.path.insert(0, ".")
import polars as pl
import argparse
import polars_init
import numpy as np


def analysis(fact_poly_info_file: str):

    df = pl.read_csv(fact_poly_info_file, separator="\t")
    df = df.filter(pl.col("qClean")).with_columns(
        [pl.format("{}{}", pl.col("rBase"), pl.col("rRepeats")).alias("key")]
    )

    df = (
        df.group_by(["key", "qRepeats"])
        .agg([pl.len().alias("cnt")])
        .with_columns(
            [
                (pl.col("cnt") / pl.col("cnt").sum().over(pl.col("key"))).alias(
                    "ratio"
                ),
                pl.col("key").str.slice(1).cast(pl.Int32).alias("ref_hp_cnt"),
            ]
        )
    ).sort(by=["key", "qRepeats"])

    df = df.with_columns(
        [
            (pl.col("ref_hp_cnt") + 1).alias("ref_hp_cnt_ADD_1"),
            (pl.col("ref_hp_cnt") - 1).alias("ref_hp_cnt_MINUS_1"),
        ]
    )

    df_eq = df.filter(pl.col("ref_hp_cnt") == pl.col("qRepeats"))
    df_minus_one = df.filter(pl.col("ref_hp_cnt_MINUS_1") == pl.col("qRepeats"))
    df_add_one = df.filter(pl.col("ref_hp_cnt_ADD_1") == pl.col("qRepeats"))

    all = (
        df_eq.join(df_minus_one, on="key", suffix="_minus_one", how="left")
        .join(df_add_one, on="key", suffix="_add_one", how="left")
        .select(
            [
                "key",
                "qRepeats",
                "cnt",
                "ratio",
                "qRepeats_minus_one",
                "cnt_minus_one",
                "ratio_minus_one",
                "qRepeats_add_one",
                "cnt_add_one",
                "ratio_add_one",
            ]
        )
    )
    print(all.head(100))

    print(df.select(["key", "qRepeats", "cnt", "ratio"]).sort(by=["key", "qRepeats"]))
    return df


def map(model: pl.DataFrame, base, cnts):
    results = np.array([-100] * 10)
    for n in range(0, len(results)):

        for cnt in cnts:
            prob = model.filter(
                (pl.col("key") == f"{base}{n}").and_(pl.col("qRepeats") == cnt)
            ).select(pl.col("ratio").log())
            if prob.shape[0] == 0:
                results[n] += -100
            else:
                prob = float(prob.to_pandas()["ratio"])
                results[n] += prob
    results = np.exp(results)
    print(results)
    print(results.argmax())
    pass


def main(args):
    df = analysis(args.fact_poly_info_file)
    return

    examples = [
        ["A", [3, 4, 4, 3, 3, 4, 4, 3, 4]],
        ["A", [4, 4, 3, 3, 3, 3, 3]],
        ["T", [6, 5, 5, 3, 6, 3, 4, 4]],
        ["A", [5, 3, 5, 5, 4, 4, 3, 5]],
        ["T", [2, 2, 2, 2, 2, 4, 5]],
        ["A", [6, 6, 6, 5, 6, 5, 6]],
        ["A", [3, 2, 3, 5, 3, 4, 4]],
        ["T", [5, 3, 5, 4, 4, 6, 5]],
        ["A", [3, 3, 4, 4, 3, 4, 4, 4]],
        ["T", [4, 4, 3, 4, 4, 2, 4]],
        ["A", [2, 2, 3, 2, 2, 3]],
        ["T", [6, 5, 6, 6, 5, 5, 6, 6]],
        ["A", [3, 2, 3, 2, 3, 3, 2]],
        ["G", [2, 2, 1, 3, 3, 2, 3, 2, 2, 2]],
        ["T", [6, 5, 6, 5, 6, 6, 6]],
    ]
    for example in examples:
        map(df, example[0], example[1])

    pass


if __name__ == "__main__":
    polars_init.polars_env_init()
    cli_params = {
        "fact_poly_info_file": "/data/adam/20241203_240901Y0005_Run0001/analysis/fact_poly_info.csv"
    }
    main(argparse.Namespace(**cli_params))
