import pysam
from tqdm import tqdm
import polars as pl
import os
import argparse
from typing import Tuple
import numpy as np


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


def read_bam_info(
    bam_file: str, channel_tag: str, min_rq: float = None
) -> Tuple[pl.DataFrame, bool]:
    assert channel_tag in ("ch", "zm")
    channels = []
    nps = []
    seq_lens = []
    rqs = []
    cxs = []

    is_sbr = False

    with pysam.AlignmentFile(
        filename=bam_file, mode="rb", check_sq=False, threads=40
    ) as reader:
        for record in tqdm(reader.fetch(until_eof=True), desc=f"reading {bam_file}"):

            rq = 0
            if record.has_tag("rq"):
                rq = float(record.get_tag("rq"))
            else:
                rq = 1 - 10 ** (float(record.get_tag("cq")) / -10)

            if min_rq is not None and rq < min_rq:
                continue

            cx = 3
            if record.has_tag("cx"):
                cx = record.get_tag("cx")
                is_sbr = True
            if not record.has_tag(channel_tag):
                ch = int(record.query_name.split("_")[1])
            else:
                ch = int(record.get_tag(channel_tag))
            channels.append(ch)
            seq_lens.append(record.query_length)
            rqs.append(rq)
            cxs.append(cx)
            np = 1
            if record.has_tag("np"):
                np = int(record.get_tag("np"))
            np = 25 if np >= 25 else np

            nps.append(np)

    df = pl.DataFrame(
        {"ch": channels, "seq_len": seq_lens, "rq": rqs, "np": nps, "cx": cxs}
    )

    df = df.with_columns([q2phreq_expr("rq", "phreq")])
    return (df, is_sbr)


def compute_n50(lengths: np.ndarray):
    lengths[::-1].sort()
    total = lengths.sum()
    half_total = total / 2

    cumsum = np.cumsum(lengths)
    n50_index = np.searchsorted(cumsum, half_total)
    return lengths[n50_index]


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

    lengths = df.select([pl.col("seq_len")]).to_pandas()["seq_len"].to_numpy()
    print("N50:{}".format(compute_n50(lengths)))

    print(res)

    # q
    res = df.select(
        [
            pl.col("phreq")
            .ge(pl.lit(8))
            .sum()
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String)
            .alias("≥Q8"),
            pl.col("phreq")
            .ge(pl.lit(10))
            .sum()
            .alias("≥Q10")
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String),
            pl.col("phreq")
            .ge(pl.lit(15))
            .sum()
            .alias("≥Q15")
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String),
            pl.col("phreq")
            .ge(pl.lit(20))
            .sum()
            .alias("≥Q20")
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String),
            pl.col("phreq")
            .ge(pl.lit(30))
            .sum()
            .alias("≥Q30")
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String),
        ]
    )
    print(res)

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
        .with_columns([pl.col("numChannels").sum().alias("TotChannels")])
        .with_columns(
            [
                (pl.col("numChannels") / pl.col("TotChannels"))
                .map_elements(lambda x: f"{x: .2%}", return_dtype=pl.String)
                .alias("ratio")
            ]
        )
        .drop(["TotChannels"])
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
    for bam_path in args.bams:
        print("")
        print("")

        df, is_sbr = read_bam_info(
            bam_path, channel_tag=args.channel_tag, min_rq=args.min_rq
        )

        if is_sbr:
            stat_subreads(df=df)
        else:
            stat_channel_reads(df=df)

    pass


def main_cli():
    parser = argparse.ArgumentParser(prog="bam basic stat")
    parser.add_argument("bams", nargs="+", type=str)
    parser.add_argument(
        "--channel-tag", type=str, default="ch", help="zm/ch", dest="channel_tag"
    )
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
