import pysam
from tqdm import tqdm
import polars as pl
import os
import argparse


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
) -> pl.DataFrame:
    assert channel_tag in ("ch", "zm")
    channels = []
    nps = []
    seq_lens = []
    rqs = []
    with pysam.AlignmentFile(
        filename=bam_file, mode="rb", check_sq=False, threads=40
    ) as reader:
        for record in tqdm(reader.fetch(until_eof=True), desc=f"reading {bam_file}"):
            rq = float(record.get_tag("rq"))
            if min_rq is not None and rq < min_rq:
                continue

            ch = int(record.get_tag(channel_tag))
            channels.append(ch)
            seq_lens.append(record.query_length)
            rqs.append(rq)

            np = int(record.get_tag("np"))
            np = 25 if np >= 25 else np

            nps.append(np)

    df = pl.DataFrame({"ch": channels, "seq_len": seq_lens, "rq": rqs, "np": nps})

    df = df.with_columns([q2phreq_expr("rq", "phreq")])
    return df


def stat_channel_reads(df: pl.DataFrame):
    res = df.select(
        [
            pl.len().alias("numChannels"),
            pl.col("seq_len")
            .sum()
            .map_elements(lambda x: f"{x:,}", return_dtype=pl.String)
            .alias("num_bases"),
            pl.col("seq_len").mean().cast(pl.Int32).alias("seq_len_mean"),
            pl.col("seq_len").median().cast(pl.Int32).alias("seq_len_median"),
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
                pl.col("phreq").min().alias("minQv"),
                pl.quantile("phreq", quantile=0.05).alias("Qv_5"),
                pl.quantile("phreq", quantile=0.25).alias("Qv_25"),
                pl.quantile("phreq", quantile=0.50).alias("Qv_50"),
                pl.col("phreq").max().alias("maxQv"),
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
            ]
        )
        .sort(by=["np"], descending=[False])
    )

    print(res)
    pass


def stat_subreads(df: pl.DataFrame):
    res = (
        df.group_by(["ch"])
        .agg(
            [
                pl.len().alias("oriPasses"),
                pl.col("seq_len").sum().alias("num_bases"),
                pl.col("seq_len").median().alias("seq_len_median"),
                pl.col("seq_len").mean().alias("seq_len_mean"),
            ]
        )
        .select(
            [
                pl.len().alias("numChannels"),
                pl.col("num_bases")
                .sum()
                .map_elements(lambda x: f"{x:,}", return_dtype=pl.String)
                .alias("num_bases"),
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
    )

    print(res)


def bam_basic_stat():
    polars_env_init()

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

    for bam_path in args.bams:
        print("")
        print("")

        df = read_bam_info(bam_path, channel_tag=args.channel_tag, min_rq=args.min_rq)
        max_num_passes = (
            df.group_by(["ch"])
            .agg([pl.len().alias("numPasses")])
            .select([pl.col("numPasses").max()])
            .to_pandas()["numPasses"]
            .values[0]
        )
        if max_num_passes == 1:
            stat_channel_reads(df=df)
        else:
            stat_subreads(df=df)


if __name__ == "__main__":
    bam_basic_stat()
