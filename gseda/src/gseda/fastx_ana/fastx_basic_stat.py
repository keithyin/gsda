import pysam
from tqdm import tqdm
import polars as pl
import os
import argparse
from typing import Tuple
import numpy as np


def phred33_to_rq(qual_str: str) -> float:
    data = [ord(char) - 33 for char in qual_str]
    data = np.array(data)
    data = 1 - np.mean(np.power(10, data / -10))
    return float(data)


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "300"
    os.environ["POLARS_FMT_STR_LEN"] = "100"


def q2phreq_expr(inp_name, oup_name=None):
    oup_name = oup_name if oup_name is not None else inp_name
    return (
        -10.0
        * (
            1
            - pl.when(pl.col(inp_name) > (1 - 1e-10))
            .then(1 - 1e-10)
            .otherwise(pl.col(inp_name))
        ).log10()
    ).alias(oup_name)


def read_fastx_info(fastx_file: str, min_rq: float = None) -> Tuple[pl.DataFrame, bool]:
    nps = []
    seq_lens = []
    rqs = []
    cxs = []

    is_sbr = False

    with pysam.FastxFile(fastx_file) as reader:
        for entry in tqdm(reader, desc=f"reading {fastx_file}"):
            seq_lens.append(len(entry.sequence))
            rq = 0.0
            if entry.quality is not None:
                rq = phred33_to_rq(entry.quality)
            rqs.append(rq)
            nps.append(0)
            cxs.append(3)

    df = pl.DataFrame({"seq_len": seq_lens, "rq": rqs, "np": nps, "cx": cxs})

    df = df.with_columns([q2phreq_expr("rq", "phreq")])
    return (df, is_sbr)


def stat_channel_reads(df: pl.DataFrame):
    res = df.select(
        [
            pl.len().alias("numChannels"),
            pl.col("seq_len")
            .sum()
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String)
            .alias("num_bases"),
            pl.col("seq_len").mean().cast(pl.Int32).alias("seq_len_mean"),
            pl.concat_str(
                pl.col("seq_len").min().cast(pl.Int32),
                pl.quantile("seq_len", quantile=0.05).cast(pl.Int32),
                pl.quantile("seq_len", quantile=0.25).cast(pl.Int32),
                pl.quantile("seq_len", quantile=0.5).cast(pl.Int32),
                pl.quantile("seq_len", quantile=0.75).cast(pl.Int32),
                pl.quantile("seq_len", quantile=0.95).cast(pl.Int32),
                pl.quantile("seq_len", quantile=0.99).cast(pl.Int32),
                pl.col("seq_len").max().cast(pl.Int32),
                separator=", ",
            ).alias("SeqLen0_5_25_50_75_95_99_100"),
        ]
    )
    print(res)

    # q
    res = df.select(
        [
            (pl.col("phreq").ge(pl.lit(8)).sum() / pl.len())
            .alias("≥Q8")
            .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
            (pl.col("phreq").ge(pl.lit(10)).sum() / pl.len())
            .alias("≥Q10")
            .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
            (pl.col("phreq").ge(pl.lit(15)).sum() / pl.len())
            .alias("≥Q15")
            .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
            (pl.col("phreq").ge(pl.lit(20)).sum() / pl.len())
            .alias("≥Q20")
            .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
            (pl.col("phreq").ge(pl.lit(30)).sum() / pl.len())
            .alias("≥Q30")
            .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
            pl.col("phreq").mean().alias("MeanQValue"),
            pl.col("phreq").median().alias("MedianQValue"),
        ]
    )
    print(res)

    # np
    res = (
        df.group_by("np")
        .agg(
            [
                pl.len().alias("numChannels"),
                pl.col("phreq")
                .min()
                .map_elements(lambda x: f"{x: .2f}", return_dtype=pl.String)
                .alias("minQv"),
                pl.quantile("phreq", quantile=0.05)
                .map_elements(lambda x: f"{x: .2f}", return_dtype=pl.String)
                .alias("Qv_5"),
                pl.quantile("phreq", quantile=0.25)
                .map_elements(lambda x: f"{x: .2f}", return_dtype=pl.String)
                .alias("Qv_25"),
                pl.quantile("phreq", quantile=0.50)
                .map_elements(lambda x: f"{x: .2f}", return_dtype=pl.String)
                .alias("Qv_50"),
                pl.col("phreq")
                .max()
                .map_elements(lambda x: f"{x: .2f}", return_dtype=pl.String)
                .alias("maxQv"),
                (pl.col("phreq").ge(pl.lit(8)).sum() / pl.len())
                .alias("≥Q8")
                .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
                (pl.col("phreq").ge(pl.lit(10)).sum() / pl.len())
                .alias("≥Q10")
                .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
                (pl.col("phreq").ge(pl.lit(15)).sum() / pl.len())
                .alias("≥Q15")
                .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
                (pl.col("phreq").ge(pl.lit(20)).sum() / pl.len())
                .alias("≥Q20")
                .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
                (pl.col("phreq").ge(pl.lit(30)).sum() / pl.len())
                .alias("≥Q30")
                .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String),
                pl.col("seq_len").mean().cast(pl.Int32).alias("seq_len_mean"),
                pl.col("seq_len").median().cast(pl.Int32).alias("seq_len_median"),
                pl.quantile("seq_len", quantile=0.99)
                .cast(dtype=pl.Int32)
                .alias("seq_len_99"),
                pl.col("seq_len").max().cast(pl.Int32).alias("seq_len_max"),
            ]
        )
        .sort(by=["np"], descending=[False])
    )

    print("------------------------Channel----------------------------")
    print(res)
    pass


def len_dist(channel_level_info: pl.DataFrame) -> pl.DataFrame:
    res = channel_level_info.select(
        [
            pl.len().alias("numChannels"),
            pl.col("num_bases")
            .sum()
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String)
            .alias("num_bases"),
            pl.col("full_len_bases")
            .sum()
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String)
            .alias("full_len_bases"),
            pl.concat_str(
                pl.col("oriPasses").min(),
                pl.quantile("oriPasses", quantile=0.05).cast(pl.Int32),
                pl.quantile("oriPasses", quantile=0.25).cast(pl.Int32),
                pl.quantile("oriPasses", quantile=0.5).cast(pl.Int32),
                pl.quantile("oriPasses", quantile=0.75).cast(pl.Int32),
                pl.quantile("oriPasses", quantile=0.95).cast(pl.Int32),
                pl.col("oriPasses").max(),
                separator=", ",
            ).alias("oriPasses0_5_25_50_75_95_100"),
            pl.concat_str(
                pl.col("seq_len_median").min().cast(pl.Int32),
                pl.quantile("seq_len_median", quantile=0.05).cast(pl.Int32),
                pl.quantile("seq_len_median", quantile=0.25).cast(pl.Int32),
                pl.quantile("seq_len_median", quantile=0.5).cast(pl.Int32),
                pl.quantile("seq_len_median", quantile=0.75).cast(pl.Int32),
                pl.quantile("seq_len_median", quantile=0.95).cast(pl.Int32),
                pl.quantile("seq_len_median", quantile=0.99).cast(pl.Int32),
                pl.col("seq_len_median").max().cast(pl.Int32),
                separator=", ",
            ).alias("SeqLenMedian0_5_25_50_75_95_99_100"),
            pl.concat_str(
                pl.col("seq_len_mean").min().cast(pl.Int32),
                pl.quantile("seq_len_mean", quantile=0.05).cast(pl.Int32),
                pl.quantile("seq_len_mean", quantile=0.25).cast(pl.Int32),
                pl.quantile("seq_len_mean", quantile=0.5).cast(pl.Int32),
                pl.quantile("seq_len_mean", quantile=0.75).cast(pl.Int32),
                pl.quantile("seq_len_mean", quantile=0.95).cast(pl.Int32),
                pl.quantile("seq_len_mean", quantile=0.99).cast(pl.Int32),
                pl.col("seq_len_mean").max().cast(pl.Int32),
                separator=", ",
            ).alias("SeqLenMean0_5_25_50_75_95_99_100"),
        ]
    )
    return res


def stat_subreads(df: pl.DataFrame):
    df = df.with_columns([(pl.col("cx") == 3).cast(pl.Int32).alias("is_full_len")])

    channel_level_info = df.group_by(["ch"]).agg(
        [
            pl.col("is_full_len").sum().alias("oriPasses"),
            pl.col("seq_len").sum().alias("num_bases"),
            pl.col("seq_len")
            .filter(pl.col("is_full_len") == 1)
            .sum()
            .alias("full_len_bases"),
            pl.col("seq_len")
            .filter(pl.col("is_full_len") == 1)
            .median()
            .alias("seq_len_median"),
            pl.col("seq_len")
            .filter(pl.col("is_full_len") == 1)
            .mean()
            .alias("seq_len_mean"),
        ]
    )
    print("------------------------Adapter.bam----------------------------")

    res = len_dist(channel_level_info=channel_level_info)
    print("------------------------Passes>=0----------------------------")
    print(res)

    res = len_dist(
        channel_level_info=channel_level_info.filter(pl.col("oriPasses") >= 1)
    )
    print("------------------------Passes>=1----------------------------")
    print(res)

    res = len_dist(
        channel_level_info=channel_level_info.filter(pl.col("oriPasses") >= 3)
    )
    print("------------------------Passes>=3----------------------------")
    print(res)


def main(args):
    polars_env_init()
    for fastx in args.fastx:
        print("")
        print("")
        df, is_sbr = read_fastx_info(fastx, min_rq=args.min_rq)
        stat_channel_reads(df=df)


def main_cli():
    parser = argparse.ArgumentParser(prog="bam basic stat")
    parser.add_argument("fastx", nargs="+", type=str)
    parser.add_argument(
        "--min-rq",
        type=float,
        default=None,
        help="only the rq ≥ min-rq will be considered",
        dest="min_rq",
    )
    args = parser.parse_args()
    main(args=args)


if __name__ == "__main__":
    main_cli()
