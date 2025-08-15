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

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
print(cur_dir)
sys.path.append(cur_dir)
import env_prepare # noqa: E402


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def mm2_version_check():
    oup = subprocess.getoutput("gsmm2-aligned-metric -V")
    oup = oup.strip()
    version_str = oup.rsplit(" ", maxsplit=1)[1]

    logging.info(f"gsmm2-aligned-metric Version: {version_str}")
    mm2_version = semver.Version.parse(version_str)
    expected_version = "0.24.0"
    assert mm2_version >= semver.Version.parse(
        expected_version
    ), f"current mm2 version:{mm2_version} < {expected_version}, try 'cargo uninstall mm2; cargo install mm2@={expected_version}' "


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "300"
    os.environ["POLARS_FMT_STR_LEN"] = "100"


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
    short_aln=False,
    np_range=None,
    rq_range=None
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
    if np_range is not None:
        cmd += f" --np-range {np_range}"
        
    if rq_range is not None:
        cmd += f" --rq-range {rq_range}"
    
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


def analysis_segs(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.with_columns(
            [
                pl.when(pl.col("segs") >= 20)
                .then(pl.lit(20))
                .otherwise(pl.col("segs"))
                .alias("segs"),
                pl.when(pl.col("segs") >= 20)
                .then(pl.lit("≥20"))
                .otherwise(pl.format("={}", pl.col("segs")))
                .alias("tag"),
            ]
        )
        .group_by(["segs", "tag"])
        .agg([pl.len().alias("cnt")])
        .with_columns(
            [
                (pl.col("cnt") / pl.col("cnt").sum().over(pl.lit("1"))).alias("ratio"),
                pl.col("cnt").sum().over(pl.lit("1")).alias("tot_cnt"),
            ]
        )
        .sort("segs")
        .select(
            [
                pl.format("[SegsRatio]segs{}", pl.col("tag")).alias("name"),
                pl.col("ratio").alias("value"),
            ]
        )
    )


def analysis_segs2(
    df: pl.DataFrame,
    ovlp_low_threshold: float = 0.02,
    ovlp_high_threshold: float = 0.9,
) -> pl.DataFrame:
    return (
        df.filter(pl.col("segs") == pl.lit(2))
        .with_columns(
            [
                (pl.col("qOvlpRatio") < ovlp_low_threshold).alias("nonOvlpQuery"),
            ]
        )
        .with_columns(
            [
                (
                    pl.col("nonOvlpQuery").and_(
                        pl.col("rOvlpRatio") > ovlp_high_threshold
                    )
                ).alias("svCandidate"),
                (
                    pl.col("nonOvlpQuery").and_(
                        pl.col("rOvlpRatio") < ovlp_low_threshold
                    )
                ).alias("noCutCandidate"),
            ]
        )
        .with_columns(
            [pl.col("oriQGaps").str.split(",").list.get(
                1).cast(pl.Int32).alias("gap")]
        )
        .with_columns([pl.col("gap").lt(pl.lit(20)).alias("gap<20")])
        .with_columns(
            [
                pl.col("svCandidate").and_(pl.col("gap<20")).alias("sv"),
                pl.col("noCutCandidate").and_(pl.col("gap<20")).alias("noCut"),
            ]
        )
        .with_columns(
            [
                pl.when(pl.col("sv"))
                .then(pl.lit("sv"))
                .when(pl.col("noCut"))
                .then(pl.lit("noCut"))
                .when(pl.col("svCandidate"))
                .then(pl.lit("svCandidate"))
                .when(pl.col("noCutCandidate"))
                .then(pl.lit("noCutCandidate"))
                .otherwise(pl.lit("badCase"))
                .alias("tag")
            ]
        )
        .with_columns(
            [
                pl.when(
                    pl.col("tag").eq(pl.lit("badCase")).and_(
                        pl.col("identity") < 0.85)
                )
                .then(pl.lit("badCase-lowIdentity"))
                .otherwise(pl.col("tag"))
                .alias("tag")
            ]
        )
        .group_by(["tag"])
        .agg([pl.len().alias("cnt")])
        .select(
            [
                pl.col("tag"),
                pl.col("cnt"),
                pl.col("cnt").sum().over(pl.lit("1")).alias("tot_cnt"),
            ]
        )
        .with_columns([(pl.col("cnt") / pl.col("tot_cnt")).alias("ratio")])
        .sort(["tag"])
        .select(
            [
                pl.format("[segs=2]{}", pl.col("tag")).alias("name"),
                pl.col("ratio").alias("value"),
            ]
        )
    )


def analysis_gaps(df: pl.DataFrame) -> pl.DataFrame:

    segs_ratio = (
        df.select(
            [
                pl.when(pl.col("segs") > 1)
                .then(pl.lit("[segs>1]"))
                .otherwise(pl.lit("[segs=1]"))
                .alias("name")
            ]
        )
        .group_by("name")
        .agg([pl.len().alias("cnt")])
        .with_columns(
            [(pl.col("cnt") / pl.col("cnt").sum().over(pl.lit(1))).alias("ratio")]
        )
    )

    metric_ratio = pl.concat(
        [
            segs_ratio.filter(pl.col("name").eq(pl.lit("[segs>1]"))).select(
                [
                    pl.format("{}Cnt", pl.col("name")),
                    pl.col("cnt").cast(pl.Float64).alias("value"),
                ]
            ),
            segs_ratio.filter(pl.col("name").eq(pl.lit("[segs>1]"))).select(
                [
                    pl.format("{}Ratio", pl.col("name")),
                    pl.col("ratio").cast(pl.Float64).alias("value"),
                ]
            ),
        ]
    )

    top20_gap = (
        df.filter(pl.col("segs") > 1)
        .select(
            [
                pl.col("oriQGaps")
                .str.split(",")
                .list.slice(1, pl.col("segs") - 1)
                .explode()
                .alias("gap")
            ]
        )
        .select([pl.col("gap").cast(pl.Int32)])
        .group_by("gap")
        .agg([pl.len().alias("cnt")])
        .with_columns(
            [(pl.col("cnt") / pl.col("cnt").sum().over(pl.lit(1))).alias("ratio")]
        )
        .sort(["ratio"], descending=[True])
        .head(20)
    )

    top20_gap_metric = top20_gap.select(
        [
            pl.format("[segs>1]gap={}", pl.col("gap")).alias("name"),
            pl.col("ratio").alias("value"),
        ]
    )

    top20_tot = top20_gap.select([pl.col("ratio").sum().alias("value")]).select(
        [pl.lit("[segs>1]gapTop20Ratio").alias("name"), pl.col("value")]
    )

    # .filter(pl.col("ratio") > 0.01)
    # .sort("gap")
    # .select(
    #     [
    #         pl.format("[segs>1]gap={}", pl.col("gap")).alias("name"),
    #         pl.col("ratio").alias("value"),
    #     ]
    # )
    return pl.concat([metric_ratio, top20_tot, top20_gap_metric])


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


def stats(metric_filename, filename):
    
    # check if it is a empty metric file
    with open(metric_filename, mode="r", encoding="utf8") as metric_inp:
        idx = 0
        for (idx, _) in enumerate(metric_inp):
            pass
        if idx < 2:
            logging.warning(f"empty metric file: {metric_filename}")
            return
    
    df = pl.read_csv(
        metric_filename, separator="\t", schema_overrides={"longIndel": pl.String}
    )
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
    metric_long_indel = analisys_long_indel(df=df)
    metric_segs = analysis_segs(df)
    metric_segs2 = analysis_segs2(df)
    metric_gaps = analysis_gaps(df=df)

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

    aggr_metrics = pl.concat(
        [
            metric_aligned_not_aligned,
            metric_long_indel,
            aggr_metrics,
            metric_segs,
            metric_segs2,
            metric_gaps,
        ]
    )

    print(aggr_metrics)

    if os.path.exists(filename):
        os.remove(filename)
    aggr_metrics.write_csv(filename, include_header=True, separator="\t")


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
    bam_file: str,
    ref_fa: str,
    threads=None,
    force=False,
    short_aln=False,
    outdir=None,
    copy_bam_file=False,
    np_range=None,
    rq_range=None
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
        mm2: cargo install mm2 (version >= 0.19.0)

    Args:
        bam_file (str): bam file. only support adapter.bam
        ref_fa (str): ref genome fa file nam
        threads (int|None): threads for generating detailed metric file
        force (boolean): if force==False, use the existing metric file if exists
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

    fact_metric_filename = f"{outdir}/{stem}.gsmm2_aligned_metric_fact.csv"
    fact_metric_filename = generate_metric_file(
        bam_file,
        ref_fa,
        out_filename=fact_metric_filename,
        force=force,
        threads=threads,
        short_aln=short_aln,
        np_range=np_range,
        rq_range=rq_range
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


def main_cli():
    polars_env_init()

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--bams", nargs="+", type=str, required=True)
    parser.add_argument("--refs", nargs="+", type=str, required=True)
    parser.add_argument("--short-aln", type=int, default=0,
                        help="for query or target in [30, 200]", dest="short_aln")
    parser.add_argument(
        "-f",
        action="store_true",
        default=False,
        help="regenerate the metric file if exists",
    )
    parser.add_argument("--np-range", type=str, default=None, dest="np_range", help="1:3,5,7:9 means [[1, 3], [5, 5], [7, 9]]. only valid for bam input that contains np field")
    parser.add_argument("--rq-range", type=str, default=None, dest="rq_range", help="0.9:1.1 means 0.9<=rq<=1.1. only valid for bam input that contains rq field")
    
    args = parser.parse_args()

    bam_files = args.bams
    refs = args.refs
    if len(refs) == 1:
        refs = refs * len(bam_files)

    assert len(bam_files) == len(refs)

    for bam, ref in zip(bam_files, refs):
        main(bam_file=bam, ref_fa=ref, force=args.f,
             short_aln=args.short_aln == 1, np_range=args.np_range, rq_range=args.rq_range)


if __name__ == "__main__":
    main_cli()
