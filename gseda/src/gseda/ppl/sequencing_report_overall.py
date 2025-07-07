import gen_output_bam
import env_prepare
import subprocess
import pathlib
import os
import logging
import polars as pl
import shutil
import argparse
from multiprocessing import cpu_count
import os
import semver
import sys
import pysam
from tqdm import tqdm
import numpy as np
import math
cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
print(cur_dir)
sys.path.append(cur_dir)


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def compute_n50(lengths: np.ndarray):
    lengths[::-1].sort()
    total = lengths.sum()
    half_total = total / 2

    cumsum = np.cumsum(lengths)
    n50_index = np.searchsorted(cumsum, half_total)
    return lengths[n50_index]


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def bam_basic_ana(bam_file: str, type: str):
    num_passes = []
    rq_values = []
    read_lengths = []
    with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False, threads=os.cpu_count()) as in_bam:
        for record in tqdm(in_bam.fetch(until_eof=True), desc=f"bam_basic_ana, reading {bam_file}"):
            num_pass = 0
            if record.has_tag("np"):
                num_pass = record.get_tag("np")

            if record.has_tag("cq"):  # called.bam or adapter.bam
                rq = record.get_tag('cq')
            else:
                rq = -10 * math.log10(1-record.get_tag("rq"))

            num_passes.append(num_pass)
            rq_values.append(rq)
            read_lengths.append(len(record.query_sequence))

    n50 = compute_n50(np.array(read_lengths))
    df = pl.DataFrame({
        "np": num_passes,
        "rq": rq_values,
        "read_len": read_lengths,
    })

    basic_metric = df.select([
        pl.len().cast(pl.Float64).alias(f"ReadsNum({type})"),
        pl.col("read_len").sum().cast(pl.Float64).alias(f"TotBases({type})"),
        pl.lit(n50).cast(pl.Float64).alias(f"N50({type})"),

        pl.col("read_len").mean().cast(
            pl.Float64).alias(f"ReadLenAvg({type})"),
        pl.col("read_len").min().cast(pl.Float64).alias(f"ReadLenMin({type})"),
        pl.quantile("read_len", 0.25).cast(
            pl.Float64).alias(f"ReadLenP25({type})"),
        pl.quantile("read_len", 0.5).cast(
            pl.Float64).alias(f"ReadLenP50({type})"),
        pl.quantile("read_len", 0.75).cast(
            pl.Float64).alias(f"ReadLenP75({type})"),
        pl.col("read_len").max().cast(pl.Float64).alias(f"ReadLenMax({type})"),

        pl.col("rq").filter(pl.col("rq") >= pl.lit(8)
                            ).count().cast(pl.Float64).alias(f"≥Q8({type})"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(10)
                            ).count().cast(pl.Float64).alias(f"≥Q10({type})"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(15)
                            ).count().cast(pl.Float64).alias(f"≥Q15({type})"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(20)
                            ).count().cast(pl.Float64).alias(f"≥Q20({type})"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(30)
                            ).count().cast(pl.Float64).alias(f"≥Q30({type})"),

        (pl.col("np").filter((pl.col("np") == pl.lit(4)).and_(pl.col("rq") >= pl.lit(20))).count(
        ) / pl.col("np").filter(pl.col("np") == pl.lit(4)).count()).cast(pl.Float64).alias(f"4xQ20({type})")
    ])

    basic_metric = basic_metric.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )
    return basic_metric


def generate_metric_file(
    bam_file: str,
    ref_fasta: str,
    out_filename: str,
    force: bool = False,
    threads=None,
    no_supp=False,
    no_mar=False,
    short_aln=False,
) -> str:

    if no_supp and no_mar:
        raise ValueError("no_supp, no_mar can't be all true")

    if force and os.path.exists(out_filename):
        os.remove(out_filename)

    if os.path.exists(out_filename):
        logging.warning("use the existing metric file : %s", out_filename)
        return out_filename

    threads = cpu_count() if threads is None else threads
    cmd = f"""gsmm2-aligned-metric --threads {threads} \
            -q {bam_file} \
            -t {ref_fasta} \
            --out {out_filename} \
            --kmer 11 \
            --wins 1 """
    if short_aln:
        cmd += " --short-aln"

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    return out_filename


def analysis_aligned(df: pl.DataFrame) -> pl.DataFrame:
    df = (
        df.select(
            [
                pl.when(pl.col("rname").eq(pl.lit("")).or_(
                    pl.col("rname").is_null()))
                .then(pl.lit("notAligned"))
                .otherwise(pl.lit("aligned"))
                .alias("name")
            ]
        )
        .group_by(["name"])
        .agg([pl.len().alias("cnt")])
        .with_columns(
            [(pl.col("cnt") / pl.col("cnt").sum().over(pl.lit(1))).alias("ratio")]
        )
    )

    metric_cnt = df.select(
        [pl.col("name"), pl.col("cnt").cast(pl.Float64).alias("value")]
    ).sort(by=["name"])

    metric_ratio = df.select(
        [pl.format("{}Ratio", pl.col("name")), pl.col("ratio").alias("value")]
    ).sort(by=["name"])

    return pl.concat([metric_cnt, metric_ratio])


def align_stats(metric_filename, filename):
    df = pl.read_csv(
        metric_filename, separator="\t", schema_overrides={"longIndel": pl.String}
    )

    print(df.select(
        [(pl.col("covlen").sum() / pl.col("qlen").sum()).alias("GlobalQueryCoverage")]))

    metric_aligned_not_aligned = analysis_aligned(df=df)
    df = df.filter(pl.col("rname") != "")
    # print(df.head(2))
    # print(
    #     df.filter(pl.col("segs") > 1)
    #     .head(2)
    #     .select(
    #         [
    #             "qname",
    #             "rname",
    #             "qlen",
    #             "segs",
    #             "queryCoverage",
    #             "identity",
    #             "oriAlignInfo",
    #             "oriQGaps",
    #             "qOvlp",
    #             "qOvlpRatio",
    #             "rOvlpRatio",
    #             "mergedQrySpan",
    #         ]
    #     )
    # )

    identity_metric = df.select([
        (pl.col("identity").filter(pl.col("identity") >= pl.lit(0.83)
                                   ).count() / pl.len()).cast(pl.Float64).alias("identity≥0.83"),
        (pl.col("identity").filter(pl.col("identity") >= pl.lit(0.90)
                                   ).count() / pl.len()).cast(pl.Float64).alias("identity≥0.90"),

    ]).transpose(
        include_header=True, header_name="name", column_names=["value"]
    )

    df = df.with_columns(row_align_span())
    aggr_metrics = df.select(aggr_expressions())
    aggr_metrics = aggr_metrics.select(
        [
            pl.col("alignedRatio"),
            pl.col("queryCoverage3"),
            pl.col("identity"),
            pl.col("identity-p50"),
        ]
    )

    aggr_metrics = aggr_metrics.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )

    all_metrics = []

    all_metrics.append(metric_aligned_not_aligned)
    all_metrics.append(aggr_metrics)
    all_metrics.append(identity_metric)

    all_metrics = pl.concat(
        all_metrics
    )

    return all_metrics


def aggr_expressions():

    exprs = [
        pl.col("covlen").sum().cast(pl.Float64).alias("alignedBases"),
        (pl.col("primaryCovlen").sum() / pl.col("qlen").sum()).alias("queryCoverage"),
        (pl.col("miscCovlen").sum() / pl.col("qlen").sum()).alias("queryCoverage2"),
        (pl.col("covlen").sum() / pl.col("qlen").sum()).alias("queryCoverage3"),
        pl.quantile("queryCoverage", 0.25).alias("queryCoverage-p25"),
        pl.col("queryCoverage").median().alias("queryCoverage-p50"),
        pl.quantile("queryCoverage", 0.75).alias("queryCoverage-p75"),
        (pl.col("match").sum() / pl.col("alignSpan").sum()).alias("identity"),
        pl.quantile("identity", 0.25).alias("identity-p25"),
        pl.col("identity").median().alias("identity-p50"),
        pl.quantile("identity", 0.75).alias("identity-p75"),
        (pl.col("misMatch").sum() / pl.col("alignSpan").sum()).alias("mmRate"),
        (pl.col("ins").sum() / pl.col("alignSpan").sum()).alias("NHInsRate"),
        (pl.col("homoIns").sum() / pl.col("alignSpan").sum()).alias("HomoInsRate"),
        (pl.col("del").sum() / pl.col("alignSpan").sum()).alias("NHDelRate"),
        (pl.col("homoDel").sum() / pl.col("alignSpan").sum()).alias("HomoDelRate"),
    ]

    for base in "ACGT":
        exprs.extend(
            [
                (
                    pl.col(f"match-{base}").sum() /
                    pl.col(f"alignSpan-{base}").sum()
                ).alias(f"identity-{base}"),
                (
                    pl.col(f"misMatch-{base}").sum() /
                    pl.col(f"alignSpan-{base}").sum()
                ).alias(f"mmRate-{base}"),
                (pl.col(f"ins-{base}").sum() / pl.col(f"alignSpan-{base}").sum()).alias(
                    f"NHInsRate-{base}"
                ),
                (
                    pl.col(f"homoIns-{base}").sum() /
                    pl.col(f"alignSpan-{base}").sum()
                ).alias(f"HomoInsRate-{base}"),
                (pl.col(f"del-{base}").sum() / pl.col(f"alignSpan-{base}").sum()).alias(
                    f"NHDelRate-{base}"
                ),
                (
                    pl.col(f"homoDel-{base}").sum() /
                    pl.col(f"alignSpan-{base}").sum()
                ).alias(f"HomoDelRate-{base}"),
            ]
        )

    return exprs


def row_align_span():
    exprs = [
        (
            pl.col("match")
            + pl.col("misMatch")
            + pl.col("ins")
            + pl.col("homoIns")
            + pl.col("del")
            + pl.col("homoDel")
        ).alias("alignSpan")
    ]

    for base in "ACGT":
        exprs.append(
            (
                pl.col(f"match-{base}")
                + pl.col(f"misMatch-{base}")
                + pl.col(f"ins-{base}")
                + pl.col(f"homoIns-{base}")
                + pl.col(f"del-{base}")
                + pl.col(f"homoDel-{base}")
            ).alias(f"alignSpan-{base}")
        )
    return exprs


def main(
    called_bam: str,
    adatper_bam: str,
    smc_bam: str,
    ref_fa: str,
    threads=None,
    force=False,
    short_aln=False,
    outdir=None,
    copy_bam_file=False,
    disable_basic_stat=False,
    disable_align_stat=False
) -> str:
    """
        step1: generate detailed metric info
        step2: compute the aggr metric. the result aggr_metric.csv is a '\t' seperated csv file. the header is name\tvalue
            here is a demo.
            ---------aggr_metric.csv
            name    value
            queryCoverage   0.937
            ----------

    requirements:
        mm2: cargo install mm2 (version >= 0.24.0)

    Args:
        bam_file (str): bam file. call.bam/adapter.bam/smc_all_reads.bam
        ref_fa (str): ref genome fa file bam
        threads (int|None): threads for generating detailed metric file
        force (boolean): if force==False, use the existing metric file if exists
        short_aln (boolean): if query or target in [30, 200], set true
        outdir: if None, ${bam_filedir}/${bam_file_stem}-metric as outdir
        copy_bam_file: copy bam file to outdir. Set this parameter to true when the file is on the NAS.

    Return:
        (aggr_metric_filename, fact_metric_filename) (str, str)
    """
    env_prepare.check_and_install(
        "gsmm2-aligned-metric", semver.Version.parse("0.24.0"), "cargo install mm2")
    # mm2_version_check()

    if copy_bam_file:
        assert outdir is not None, "must provide outdir when copy_bam_file=True"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        new_bam_file = os.path.join(outdir, os.path.basename(bam_file))
        if os.path.exists(new_bam_file):
            raise ValueError(f"{new_bam_file} already exists")
        shutil.copy(bam_file, new_bam_file)
        bam_file = new_bam_file

    bam_filedir = os.path.dirname(bam_file)
    stem = extract_filename(bam_file)
    if outdir is None:
        outdir = os.path.join(bam_filedir, f"{stem}-metric")

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    output_bam_path = f"{outdir}/{stem}.output.bam"
    gen_output_bam.filter_called_bam(smc_bam, called_bam, output_bam_path)

    fact_metric_filename = f"{outdir}/{stem}.output.gsmm2_aligned_metric_fact.csv"
    if not disable_align_stat:
        fact_metric_filename = generate_metric_file(
            output_bam_path,
            ref_fa,
            out_filename=fact_metric_filename,
            force=force,
            threads=threads,
            short_aln=short_aln
        )
    aggr_metric_filename = f"{outdir}/{stem}.output.gsmm2_aligned_metric_aggr.csv"
    if force and os.path.exists(aggr_metric_filename):
        os.remove(aggr_metric_filename)

    all_metrics = []
    if not disable_basic_stat:
        all_metrics.append(bam_basic_ana(bam_file=called_bam, type="called"))
        all_metrics.append(bam_basic_ana(bam_file=smc_bam, type="smc"))

    if not disable_align_stat:
        all_metrics.append(
            align_stats(fact_metric_filename))

    if len(all_metrics) > 0:
        all_metrics = pl.concat(all_metrics)
        print(all_metrics)
        all_metrics.write_csv(aggr_metric_filename,
                              include_header=True, separator="\t")

    return (aggr_metric_filename, fact_metric_filename)


def main_cli():
    env_prepare.polars_env_init()

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--called", type=str, required=True, help="called bam")
    parser.add_argument("--adapter", type=str,
                        required=True, help="adapter bam")
    parser.add_argument("--smc", type=str, required=True, help="smc bam")

    parser.add_argument("--refs", type=str, default=None,
                        help="reference fasta")

    parser.add_argument("--short-aln", type=int, default=0,
                        help="for query or target in [30, 200]", dest="short_aln")

    parser.add_argument("--disable-basic-stat",
                        action="store_true", dest="disable_basic_stat")
    parser.add_argument("--disable-align-stat",
                        action="store_true", dest="disable_align_stat")

    parser.add_argument(
        "-f",
        action="store_true",
        default=False,
        help="regenerate the metric file if exists",
    )

    args = parser.parse_args()

    called_bam = args.called
    adapter_bam = args.adapter
    smc_bam = args.smc

    ref = args.ref

    main(
        called_bam=called_bam,
        adatper_bam=adapter_bam,
        smc_bam=smc_bam,
        ref_fa=ref, force=args.f,
        short_aln=args.short_aln == 1,
        disable_basic_stat=args.disable_basic_stat,
        disable_align_stat=args.disable_align_stat)


if __name__ == "__main__":
    main_cli()
