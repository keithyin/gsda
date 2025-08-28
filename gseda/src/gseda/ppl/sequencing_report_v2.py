import subprocess
import pathlib
import os
import logging
import polars as pl
import shutil
import argparse
from multiprocessing import cpu_count
import os
import sys
import pysam
from tqdm import tqdm
import numpy as np
import math
import semver
from typing import List
cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
print(cur_dir)
sys.path.append(cur_dir)
import env_prepare  # noqa: E402


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def is_empty_bam(bam_file: str):
    with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False, threads=os.cpu_count()) as in_bam:
        idx = 0
        for (idx, _) in enumerate(in_bam.fetch(until_eof=True)):
            if idx > 1:
                break
        return idx < 1


def compute_n50(lengths: np.ndarray):
    if len(lengths) == 0:
        return 0
    total = lengths.sum()
    half_total = total / 2

    cumsum = np.cumsum(lengths)
    n50_index = np.searchsorted(cumsum, half_total)
    return lengths[n50_index]


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def dump_default_aggr(aggr_metric_filename: str):
    default_aggr_data = """\
name\tvalue
reads_num\t0
tot_bases\t0
n50\t0
read_len_avg\t0
read_len_p50\t0
≥Q8\t0
≥Q10\t0
≥Q15\t0
≥Q20\t0
≥Q30\t0
4xQ20\t0
alignedRatio\t0
queryCoverage3\t0
identity\t0
identity-p50\t0
identity≥0.83\t0
identity≥0.90\t0
    """
    with open(aggr_metric_filename, mode="w", encoding="utf8") as out:
        out.write(default_aggr_data)


def dump_bam_basic(bam_file: str, output_dir: str) -> str:
    """dump bam basic

    Args:
        bam_file (str): path
        output_dir (str): outdir

    Returns:
        str: outputfile path. if the bam_file is empty. None will be returned
    """
    num_passes = []
    rq_values = []
    read_lengths = []
    qnames = []
    stem = extract_filename(bam_file)
    output_csv = f"{output_dir}/{stem}.basic.csv"
    with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False, threads=os.cpu_count()) as in_bam:
        for record in tqdm(in_bam.fetch(until_eof=True), desc=f"bam_basic_ana, reading {bam_file}"):
            num_pass = 0
            if record.has_tag("np"):
                num_pass = record.get_tag("np")

            if record.has_tag("cq"):  # called.bam or adapter.bam
                rq = record.get_tag('cq')
            else:
                rq = -10 * math.log10(1-record.get_tag("rq"))
            qnames.append(record.query_name)
            num_passes.append(num_pass)
            rq_values.append(rq)
            read_lengths.append(len(record.query_sequence))

    if len(qnames) == 0:
        return None

    df = pl.DataFrame({
        "qname": qnames,
        "np": num_passes,
        "rq": rq_values,
        "read_len": read_lengths,
    })
    df.write_csv(output_csv, separator="\t", include_header=True)
    return output_csv


def bam_basic_ana(basic_csvs: List[str]):
    basic_dataframes = []
    for basic_csv in basic_csvs:
        if basic_csv is None:
            continue
        basic_dataframes.append(pl.read_csv(basic_csv, separator="\t"))
    if len(basic_dataframes) == 0:
        return None

    df = pl.concat(basic_dataframes)
    assert isinstance(df, pl.DataFrame)
    n50 = compute_n50(df.select([pl.col("read_len")]).sort(
        "read_len", descending=True).to_numpy()[:, 0])
    basic_metric = df.select([
        pl.len().cast(pl.Float64).alias("reads_num"),
        pl.col("read_len").sum().cast(pl.Float64).alias("tot_bases"),
        pl.lit(n50).cast(pl.Float64).alias("n50"),

        pl.col("read_len").mean().cast(pl.Float64).alias("read_len_avg"),
        pl.col("read_len").min().cast(pl.Float64).alias("read_len_min"),
        pl.quantile("read_len", 0.25).cast(pl.Float64).alias("read_len_p25"),
        pl.quantile("read_len", 0.5).cast(pl.Float64).alias("read_len_p50"),
        pl.quantile("read_len", 0.75).cast(pl.Float64).alias("read_len_p75"),
        pl.col("read_len").max().cast(pl.Float64).alias("read_len_max"),

        pl.col("rq").filter(pl.col("rq") >= pl.lit(8)
                            ).count().cast(pl.Float64).alias("≥Q8"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(10)
                            ).count().cast(pl.Float64).alias("≥Q10"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(15)
                            ).count().cast(pl.Float64).alias("≥Q15"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(20)
                            ).count().cast(pl.Float64).alias("≥Q20"),
        pl.col("rq").filter(pl.col("rq") >= pl.lit(30)
                            ).count().cast(pl.Float64).alias("≥Q30"),

        (pl.col("np").filter((pl.col("np") == pl.lit(4)).and_(pl.col("rq") >= pl.lit(20))).count(
        ) / pl.col("np").filter(pl.col("np") == pl.lit(4)).count()).cast(pl.Float64).alias("4xQ20")
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


def analisys_long_indel(df: pl.DataFrame) -> pl.DataFrame:
    metric = df.select(
        [
            pl.len().alias("cnt"),
            pl.col("longIndel")
            .filter(pl.col("longIndel").is_not_null())
            .len()
            .alias("longIndelCnt"),
        ]
    ).select(
        [
            pl.lit("longIndelRatio").alias("name"),
            (pl.col("longIndelCnt") / pl.col("cnt")).alias("value"),
        ]
    )
    return metric


def align_stats(metric_filenames: List[str]):

    all_dfs = []
    for metric_fname in metric_filenames:
        if metric_fname is None:
            continue

        all_dfs.append(pl.read_csv(
            metric_fname, separator="\t", schema_overrides={"longIndel": pl.String}
        ))

    if len(all_dfs) == 0:
        return None

    df = pl.concat(all_dfs)

    global_query_coverage_metric = df.select(
        [(pl.col("covlen").sum() / pl.col("qlen").sum()).alias("GlobalQueryCoverage")])
    global_query_coverage_metric = global_query_coverage_metric.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )

    metric_aligned_not_aligned = analysis_aligned(df=df)
    metric_long_indel = analisys_long_indel(df=df)
    df = df.filter(pl.col("rname") != "")

    identity_metric = df.select([
        (pl.col("identity").filter(pl.col("identity") >= pl.lit(0.83)
                                   ).count() / pl.len()).cast(pl.Float64).alias("identity≥0.83"),
        (pl.col("identity").filter(pl.col("identity") >= pl.lit(0.90)
                                   ).count() / pl.len()).cast(pl.Float64).alias("identity≥0.90"),

    ]).transpose(
        include_header=True, header_name="name", column_names=["value"]
    )

    df = df.with_columns(
        [
            ((pl.col("qOvlpRatio") < 0.01).and_(pl.col("rOvlpRatio") < 0.01))
            .or_(
                (pl.col("qOvlpRatio") < 0.01)
                .and_(pl.col("rOvlpRatio") > 0.90)
                .and_(pl.col("oriQGaps").str.split(",").list.get(1).cast(pl.Int32) < 20)
            )
            .alias("valid")
        ]
    )

    df = df.with_columns(
        [
            pl.when(pl.col("valid"))
            .then(pl.col("covlen"))
            .otherwise(pl.col("primaryCovlen"))
            .alias("miscCovlen")
        ]
    )

    df = df.with_columns(row_align_span())
    aggr_metrics = df.select(aggr_expressions())
    aggr_metrics = aggr_metrics.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )

    all_metrics = []

    all_metrics.append(metric_aligned_not_aligned)
    all_metrics.append(aggr_metrics)
    all_metrics.append(identity_metric)
    all_metrics.append(metric_long_indel)
    all_metrics.append(global_query_coverage_metric)

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


def merge_partition(
    output_csv_path: str,
    fact_csvs: List[str] = None,
    basic_csvs: List[str] = None,
):

    all_metrics = []
    if basic_csvs is not None:
        inner = bam_basic_ana(basic_csvs=basic_csvs)
        if inner is not None:
            all_metrics.append(inner)

    if fact_csvs is not None:
        inner = align_stats(fact_csvs)
        if inner is not None:
            all_metrics.append(inner)

    if len(all_metrics) > 0:
        all_metrics = pl.concat(all_metrics)
        print(all_metrics)
        all_metrics.write_csv(output_csv_path,
                              include_header=True, separator="\t")

    return output_csv_path


def main(
    bam_file: str = None,
    ref_fa: str = None,
    threads=None,
    force=False,
    short_aln=False,
    outdir=None,
    copy_bam_file=False,
    disable_basic_stat=False,
    disable_align_stat=False,
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

    aggr_metric_filename = f"{outdir}/{stem}.gsmm2_aligned_metric_aggr.csv"
    if is_empty_bam(bam_file=bam_file):
        dump_default_aggr(aggr_metric_filename)
        return (aggr_metric_filename, None, None)

    basic_metric_filename = None
    if not disable_basic_stat:
        basic_metric_filename = dump_bam_basic(
            bam_file=bam_file, output_dir=outdir)

    fact_metric_filename = None
    if not disable_align_stat:
        fact_metric_filename = f"{outdir}/{stem}.gsmm2_aligned_metric_fact.csv"
        fact_metric_filename = generate_metric_file(
            bam_file,
            ref_fa,
            out_filename=fact_metric_filename,
            force=force,
            threads=threads,
            short_aln=short_aln
        )

    all_metrics = []
    if basic_metric_filename is not None:
        all_metrics.append(bam_basic_ana(
            basic_csvs=[basic_metric_filename]))

    if fact_metric_filename is not None:
        inner = align_stats([fact_metric_filename])
        if inner is not None:
            all_metrics.append(inner)

    if len(all_metrics) > 0:
        all_metrics = pl.concat(all_metrics)
        print(all_metrics)
        if os.path.exists(aggr_metric_filename):
            os.remove(aggr_metric_filename)

        all_metrics.write_csv(aggr_metric_filename,
                              include_header=True, separator="\t")

    return (aggr_metric_filename, fact_metric_filename, basic_metric_filename)


def main_cli():
    env_prepare.polars_env_init()

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--bams", nargs="+", type=str, default=None)
    parser.add_argument("--refs", nargs="+", type=str, default=None)
    parser.add_argument("--fact-csvs", nargs="+", type=str,
                        default=None, dest="fact_csvs")
    parser.add_argument("--basic-csvs", nargs="+", type=str,
                        default=None, dest="basic_csvs")

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

    bam_files = args.bams
    refs = args.refs
    fact_csvs = args.fact_csvs
    basic_csvs = args.basic_csvs

    if fact_csvs is not None or basic_csvs is not None:
        output_dir = None
        stem = None
        if fact_csvs:
            output_dir = os.path.dirname(fact_csvs[0])
            stem = extract_filename(fact_csvs[0])
        else:
            output_dir = os.path.dirname(basic_csvs[0])
            stem = extract_filename(basic_csvs[0])
        out_csv_path = f"{output_dir}/{stem}.all.csv"
        merge_partition(out_csv_path, fact_csvs=fact_csvs,
                        basic_csvs=basic_csvs)
    else:

        if not args.disable_align_stat and refs is not None:
            assert len(refs) > 0
        else:
            refs = [""]

        if len(refs) == 1:
            refs = refs * len(bam_files)

        assert len(bam_files) == len(refs)

        for bam, ref in zip(bam_files, refs):
            print(main(bam_file=bam, ref_fa=ref, force=args.f,
                       short_aln=args.short_aln == 1,
                       disable_basic_stat=args.disable_basic_stat,
                       disable_align_stat=args.disable_align_stat))


if __name__ == "__main__":
    main_cli()
