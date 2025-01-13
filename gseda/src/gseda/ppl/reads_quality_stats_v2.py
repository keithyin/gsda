import subprocess
import pathlib
import os
import logging
import polars as pl
import shutil
import argparse
from multiprocessing import cpu_count

import os


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "300"
    os.environ["POLARS_FMT_STR_LEN"] = "100"


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def generate_metric_file(
    bam_file: str,
    ref_fasta: str,
    out_filename: str,
    force: bool = False,
    threads=None,
    no_supp=False,
    no_mar=False,
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

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    return out_filename


def analysis_segs(df):
    print(
        df.select(
            [
                pl.col("segs").filter(pl.col("segs") == 1).len().alias("segs1"),
                pl.col("segs").filter(pl.col("segs") == 2).len().alias("segs2"),
                pl.col("segs").len().alias("segsAll"),
            ]
        ).with_columns(
            [
                (pl.col("segs1") / pl.col("segsAll")).alias("seg1Ratio"),
                (pl.col("segs2") / pl.col("segsAll")).alias("seg2Ratio"),
            ]
        )
    )


def analysis_segs2(df):
    print(
        df.filter(pl.col("segs") == 2)
        .select(
            [
                pl.col("qOvlpRatio")
                .filter(pl.col("qOvlpRatio") < 0.01)
                .len()
                .alias("nonOvlpQuery"),
                pl.col("qOvlpRatio")
                .filter((pl.col("qOvlpRatio") < 0.01).and_(pl.col("rOvlpRatio") < 0.01))
                .len()
                .alias("nonOvlpQueryAndNonOvlpRef"),
                pl.col("qOvlpRatio")
                .filter((pl.col("qOvlpRatio") < 0.01).and_(pl.col("rOvlpRatio") > 0.90))
                .len()
                .alias("nonOvlpQueryAndOvlpRef"),
                pl.len().alias("seg2"),
            ]
        )
        .with_columns(
            [
                (pl.col("nonOvlpQueryAndNonOvlpRef") / pl.col("seg2")).alias("svRatio"),
                (pl.col("nonOvlpQueryAndOvlpRef") / pl.col("seg2")).alias("noCutRatio"),
            ]
        )
    )


def analysis_segs2_nocut(df):
    print(
        df.filter(pl.col("segs") == 2)
        .filter(((pl.col("qOvlpRatio") < 0.01).and_(pl.col("rOvlpRatio") > 0.90)))
        .with_columns(
            [pl.col("oriQGaps").str.split(",").list.get(1).cast(pl.Int32).alias("gap")]
        )
        .select(
            [
                pl.len().alias("tot"),
                (pl.col("gap") < 20).cast(pl.Int32).sum().alias("gap<20"),
            ]
        )
        .with_columns([(pl.col("gap<20") / pl.col("tot")).alias("RatioOf(gap<20)")])
        # .select([pl.len().alias("cnt")
        #          ])
        # .with_columns(
        #     [
        #         (pl.col("nonOvlpQueryAndNonOvlpRef") / pl.col("seg2")).alias("svRatio"),
        #         (pl.col("nonOvlpQueryAndOvlpRef") / pl.col("seg2")).alias("noCutRatio"),
        #     ]
        # )
    )


def stats(metric_filename, filename):
    df = pl.read_csv(metric_filename, separator="\t").filter(pl.col("rname") != "")
    print(df.head(2))
    print(
        df.filter(pl.col("segs") > 1)
        .head(2)
        .select(
            [
                "qname",
                "rname",
                "qlen",
                "segs",
                "queryCoverage",
                "identity",
                "oriAlignInfo",
                "oriQGaps",
                "qOvlp",
                "qOvlpRatio",
                "rOvlpRatio",
                "mergedQrySpan",
            ]
        )
    )
    analysis_segs(df)
    analysis_segs2(df)
    analysis_segs2_nocut(df)

    df = df.with_columns(
        [
            ((pl.col("qOvlpRatio") < 0.01).and_(pl.col("rOvlpRatio") < 0.01))
            .or_((pl.col("qOvlpRatio") < 0.01).and_(pl.col("rOvlpRatio") > 0.90))
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
    print(aggr_metrics)

    aggr_metrics.write_csv(filename, include_header=True, separator="\t")


def aggr_expressions():

    exprs = [
        (pl.col("covlen").sum() / pl.col("qlen").sum()).alias("queryCoverage"),
        (pl.col("primaryCovlen").sum() / pl.col("qlen").sum()).alias("queryCoverage2"),
        (pl.col("miscCovlen").sum() / pl.col("qlen").sum()).alias("queryCoverage3"),
        (pl.col("match").sum() / pl.col("alignSpan").sum()).alias("identity"),
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
                    pl.col(f"match-{base}").sum() / pl.col(f"alignSpan-{base}").sum()
                ).alias(f"identity-{base}"),
                (
                    pl.col(f"misMatch-{base}").sum() / pl.col(f"alignSpan-{base}").sum()
                ).alias(f"mmRate-{base}"),
                (pl.col(f"ins-{base}").sum() / pl.col(f"alignSpan-{base}").sum()).alias(
                    f"NHInsRate-{base}"
                ),
                (
                    pl.col(f"homoIns-{base}").sum() / pl.col(f"alignSpan-{base}").sum()
                ).alias(f"HomoInsRate-{base}"),
                (pl.col(f"del-{base}").sum() / pl.col(f"alignSpan-{base}").sum()).alias(
                    f"NHDelRate-{base}"
                ),
                (
                    pl.col(f"homoDel-{base}").sum() / pl.col(f"alignSpan-{base}").sum()
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


def main(bam_file: str, ref_fa: str, threads=None, force=False, outdir=None) -> str:
    """
        step1: generate detailed metric info
        step2: compute the aggr metric. the result aggr_metric.csv is a '\t' seperated csv file. the header is name\tvalue
            here is a demo.
            ---------aggr_metric.csv
            name    value
            queryCoverage   0.937
            ----------

    requirements:
        mm2: cargo install mm2 (version >= 0.19.0)

    Args:
        bam_file (str): bam file. only support adapter.bam
        ref_fa (str): ref genome fa file nam
        threads (int|None): threads for generating detailed metric file
        force (boolean): if force==False, use the existing metric file if exists
        outdir: if None, ${bam_filedir}/${bam_file_stem}-metric as outdir

    Return:
        (aggr_metric_filename, fact_metric_filename) (str, str)
    """
    bam_filedir = os.path.dirname(bam_file)
    stem = extract_filename(bam_file)
    if outdir is None:
        outdir = os.path.join(bam_filedir, f"{stem}-metric")

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fact_metric_filename = f"{outdir}/{stem}.gsmm2_aligned_metric_fact.csv"
    fact_metric_filename = generate_metric_file(
        bam_file,
        ref_fa,
        out_filename=fact_metric_filename,
        force=force,
        threads=threads,
    )
    aggr_metric_filename = f"{outdir}/{stem}.gsmm2_aligned_metric_aggr.csv"
    if force and os.path.exists(aggr_metric_filename):
        os.remove(aggr_metric_filename)

    # if not os.path.exists(aggr_metric_filename):
    stats(fact_metric_filename, filename=aggr_metric_filename)
    # else:
    #     logging.warning(
    #         "aggr_metric_file exists, use existing one. %s", aggr_metric_filename
    #     )
    return (aggr_metric_filename, fact_metric_filename)


def test_stat():
    fact_bam_basic = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/fact_aligned_bam_bam_basic.csv"
    aggr_metric_filename = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/aggr_metric.csv"

    with open(aggr_metric_filename, encoding="utf8", mode="w") as file_h:
        file_h.write(f"name\tvalue\n")
        stats(fact_bam_basic, file_h=file_h)


def sv_identification(ori_align_info: str):
    """structural variant identification
    rule:
        small ovlp between aligned segments
        different segments align to the different reference regions
    """
    if ";" not in ori_align_info:
        return False

    align_regions = ori_align_info.split(";")
    align_regions = [align_region[:-1] for align_region in align_regions]
    # align_regions = [align_region.split(":")]
    pass


def adapter_remover_error_identification():
    """adapter remover error identification

    rule:
        small ovlp between aligned segments
        different segments align to the similar reference region

        if gap < 10: treat as adapter_missing
        if gap > 10: treat as adapter_lowq

    """
    pass


if __name__ == "__main__":
    polars_env_init()

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--bams", nargs="+", type=str, required=True)
    parser.add_argument("--refs", nargs="+", type=str, required=True)
    parser.add_argument("-f", action="store_true", default=False)

    args = parser.parse_args()

    bam_files = args.bams
    refs = args.refs
    if len(refs) == 1:
        refs = refs * len(bam_files)

    assert len(bam_files) == len(refs)

    # bam_files = [
    #     "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter.bam",
    #     # "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0002_adapter.bam",
    #     # "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0003_adapter.bam",
    # ]
    # ref = "/data/ccs_data/MG1655.fa"

    for bam, ref in zip(bam_files, refs):
        main(bam_file=bam, ref_fa=ref, force=args.f)
    # test_stat()

    # print(merge_intervals([{"qstart": 11, "qend": 771}]))
