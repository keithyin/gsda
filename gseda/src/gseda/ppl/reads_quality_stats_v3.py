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
import threading
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

# deprecated ...
logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "300"
    os.environ["POLARS_FMT_STR_LEN"] = "100"


def mm2_version_check():
    oup = subprocess.getoutput("gsmm2-aligned-metric -V")
    oup = oup.strip()
    version_str = oup.rsplit(" ", maxsplit=1)[1]

    logging.info(f"gsmm2-aligned-metric Version: {version_str}")
    mm2_version = semver.Version.parse(version_str)
    expected_version = "0.21.0"
    assert mm2_version >= semver.Version.parse(
        expected_version
    ), f"current mm2 version:{mm2_version} < {expected_version}, try 'cargo uninstall mm2; cargo install mm2@={expected_version}' "


def gsetl_version_check():
    oup = subprocess.getoutput("gsetl -V")
    oup = oup.strip()
    version_str = oup.rsplit(" ", maxsplit=1)[1]

    logging.info(f"gsetl Version: {version_str}")
    gsetl_version = semver.Version.parse(version_str)
    expected_version = "0.7.1"
    assert gsetl_version >= semver.Version.parse(
        expected_version
    ), f"current gsetl version:{gsetl_version} < {expected_version}, try 'cargo uninstall gsetl; cargo install gsetl@={expected_version}' "


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def generate_aligned_metric_fact_file(
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


def stats(metric_filename: str, filename: str):
    df = pl.read_csv(
        metric_filename, separator="\t", infer_schema_length=3000, schema_overrides={"longIndel": pl.String}
    )
    metric_aligned_not_aligned = analysis_aligned(df=df)
    df = df.filter(pl.col("rname") != "")
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

    aggr_metrics.write_csv(filename, include_header=True, separator="\t")
    identity_min = 0.6
    identity_max = 1.0
    if df.shape[0] > 10000:
        sampled_identity = df.select([pl.col("identity").clip(
            lower_bound=identity_min, upper_bound=identity_max)]).sample(n=10000, seed=2025).to_pandas()
    else:
        sampled_identity = df.to_pandas()

    identity_hist_filename = "{}.idenity_hist.png".format(
        filename.rsplit(".", maxsplit=1)[0])

    plot_histgram(sampled_identity["identity"], fname=identity_hist_filename, xlabel="Identity",
                  ylabel="Count", title="SampledIdentityHist", xlim=(identity_min, identity_max))


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
        ((pl.col("ins").sum() + pl.col("homoIns").sum()) /
         pl.col("alignSpan").sum()).alias("insRate"),

        (pl.col("del").sum() / pl.col("alignSpan").sum()).alias("NHDelRate"),
        (pl.col("homoDel").sum() / pl.col("alignSpan").sum()).alias("HomoDelRate"),
        ((pl.col("homoDel").sum() + pl.col("del").sum()) /
         pl.col("alignSpan").sum()).alias("delRate"),
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

                (
                    (pl.col(f"homoIns-{base}").sum() + pl.col(f"ins-{base}").sum()
                     ) / pl.col(f"alignSpan-{base}").sum()
                ).alias(f"insRate-{base}"),

                (pl.col(f"del-{base}").sum() / pl.col(f"alignSpan-{base}").sum()).alias(
                    f"NHDelRate-{base}"
                ),
                (
                    pl.col(f"homoDel-{base}").sum() /
                    pl.col(f"alignSpan-{base}").sum()
                ).alias(f"HomoDelRate-{base}"),

                (
                    (pl.col(f"homoDel-{base}").sum() + pl.col(f"del-{base}").sum()
                     ) / pl.col(f"alignSpan-{base}").sum()
                ).alias(f"delRate-{base}"),
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


def aligned_metric_analysis(fact_metric_filename: str, aggr_metric_filename: str, force: bool):
    if force and os.path.exists(aggr_metric_filename):
        os.remove(aggr_metric_filename)

    if os.path.exists(aggr_metric_filename):
        logging.warning(f"{aggr_metric_filename} will be override")

    # if not os.path.exists(aggr_metric_filename):
    stats(fact_metric_filename, filename=aggr_metric_filename)


def generate_non_aligned_metric_fact_file(bam_file: str, out_filepath: str, out_dir: str, force: bool):
    if not force and os.path.exists(out_filepath):
        logging.info(f"{out_filepath} exists, use the existed file")
        return

    cmd = f"gsetl --outdir {out_dir} non-aligned-bam --bam {bam_file} -o {out_filepath}"
    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)


def non_aligned_metric_analysis(fact_metric_filename: str, aggr_metric_filename: str, force: bool, out_dir: str, stem: str):
    if os.path.exists(aggr_metric_filename) and not force:
        logging.warning(f"{aggr_metric_filename} will be override")

    df = pl.read_csv(
        fact_metric_filename, separator="\t",
        infer_schema_length=3000,
        schema_overrides={
            "oe": pl.Float32
        })

    df = df.with_columns([
        (pl.col("dw_sum") * pl.lit(2)).alias("dw_sum"),
        (pl.col("ar_sum") * pl.lit(2)).alias("ar_sum"),
    ])

    whole_aggr = df.select([
        (pl.col("dw_sum").sum() / pl.col("base_cnt").sum()).alias("dw-mean"),
        (pl.col("ar_sum").sum() / pl.col("base_cnt").sum()).alias("ar-mean"),
        pl.col("cq").mean().alias("cq-mean"),
        ((pl.col("base_cnt") * pl.col("cr_mean")).sum() /
         pl.col("base_cnt").sum()).alias("cr-mean"),
        pl.col("oe").median().alias("oe-median"),
        (pl.col("base_cnt").sum() / ((pl.col("dw_sum").sum() +
         pl.col("ar_sum").sum()) * pl.lit(0.001))).alias("speed")
    ])

    base_level_aggr = df.group_by(["base"]).agg([
        (pl.col("dw_sum").sum() / pl.col("base_cnt").sum()).alias("dw-mean"),
        (pl.col("ar_sum").sum() / pl.col("base_cnt").sum()).alias("ar-mean")
    ])

    whole_aggr = whole_aggr.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )

    metrics = [whole_aggr]

    for base in ["A", "C", "G", "T"]:
        tmp = base_level_aggr.filter([pl.col("base") == pl.lit(base)])\
            .select([pl.col("dw-mean").alias(f"dw-mean-{base}"), pl.col("ar-mean").alias(f"ar-mean-{base}")])
        metrics.append(tmp.transpose(
            include_header=True, header_name="name", column_names=["value"]))

    metrics = pl.concat(metrics)

    print(metrics)

    metrics.write_csv(aggr_metric_filename,
                      include_header=True, separator="\t")

    # TODO draw plot
    read_lengths = df.group_by(["qname"])\
        .agg([pl.col("base_cnt").sum()])\
        .select([pl.col("base_cnt")])
    if read_lengths.shape[0] > 10000:
        read_lengths = read_lengths.sample(
            n=10000, seed=2025).to_pandas()["base_cnt"]
    else:
        read_lengths = read_lengths.to_pandas()["base_cnt"]

    read_length_hist_fname = f"{out_dir}/{stem}.readlength_hist.png"
    plot_histgram(read_lengths, read_length_hist_fname,
                  xlabel="ReadLength", ylabel="Count", title="SampledReadLengthHist")

    dw_hist_fname = f"{out_dir}/{stem}.dw_hist.png"
    dw_min = 0
    dw_max = 200
    dw = df.group_by(["qname"])\
        .agg([pl.col("base_cnt").sum(), pl.col("dw_sum").sum()])\
        .select([(pl.col("dw_sum") / pl.col("base_cnt")).alias("dw")])\
        .select([pl.col("dw").clip(lower_bound=dw_min, upper_bound=dw_max)])
    if dw.shape[0] > 10000:
        dw = dw.sample(n=10000, seed=2025).to_pandas()["dw"]
    else:
        dw = dw.to_pandas()["dw"]
    plot_histgram(dw, dw_hist_fname, xlabel="DwellTime", ylabel="Count",
                  title="SampledDwellTimeHist", xlim=(dw_min, dw_max))

    ar_min = 0
    ar_max = 1000
    ar = df.group_by(["qname"])\
        .agg([pl.col("base_cnt").sum(), pl.col("ar_sum").sum()])\
        .select([(pl.col("ar_sum") / pl.col("base_cnt")).alias("ar")])\
        .select([pl.col("ar").clip(lower_bound=ar_min, upper_bound=ar_max)])
    if ar.shape[0] > 10000:
        ar = ar.sample(n=10000, seed=2025).to_pandas()["ar"]
    else:
        ar = ar.to_pandas()["ar"]

    ar_hist_fname = f"{out_dir}/{stem}.ar_hist.png"
    plot_histgram(ar, ar_hist_fname, xlabel="ArrivalTime", ylabel="Count",
                  title="SampledArrivalTimeHist", xlim=(ar_min, ar_max))

    pass


def plot_histgram(data, fname, xlabel, ylabel, title, xlim=None, bins=100):
    plt.figure(figsize=(8, 6))

    # 绘制直方图
    sns.histplot(data, bins=bins, kde=True, color='skyblue', binrange=xlim)

    # 添加标题与标签
    plt.title(title, fontsize=16)
    plt.xlabel(xlabel, fontsize=14)
    if xlim is not None:
        plt.xlim(xlim)
    plt.ylabel(ylabel, fontsize=14)

    # 保存图片到文件
    plt.savefig(fname, dpi=300, bbox_inches='tight')


def main(
    bam_file: str,
    ref_fa: str,
    threads=None,
    force=False,
    outdir=None,
    copy_bam_file=False,
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

    mm2_version_check()
    gsetl_version_check()

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
    aggr_metric_filename = f"{outdir}/{stem}.gsmm2_aligned_metric_aggr.csv"
    no_aligned_fact_filename = f"{outdir}/{stem}.non_aligned_fact.csv"
    no_aligned_aggr_filename = f"{outdir}/{stem}.non_aligned_aggr.csv"

    """
    fact_metric_filename = generate_metric_file(
        bam_file,
        ref_fa,
        out_filename=fact_metric_filename,
        force=force,
        threads=threads,
    )
    """

    processes = []
    if ref_fa != "":
        aligned_fact_thread = threading.Thread(target=generate_aligned_metric_fact_file, args=(
            bam_file, ref_fa, fact_metric_filename, force, threads))
        aligned_fact_thread.start()
        processes.append(aligned_fact_thread)

    non_aligned_fact_thread = threading.Thread(target=generate_non_aligned_metric_fact_file, args=(
        bam_file, no_aligned_fact_filename, outdir, force))
    non_aligned_fact_thread.start()
    processes.append(non_aligned_fact_thread)

    for p in processes:
        p.join()
    if ref_fa != "":
        aligned_metric_analysis(fact_metric_filename,
                                aggr_metric_filename, force=force)
    non_aligned_metric_analysis(
        no_aligned_fact_filename, no_aligned_aggr_filename, force, out_dir=outdir, stem=stem)

    return (aggr_metric_filename, fact_metric_filename, no_aligned_fact_filename)


def test_stat():
    fact_bam_basic = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/fact_aligned_bam_bam_basic.csv"
    aggr_metric_filename = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/aggr_metric.csv"

    with open(aggr_metric_filename, encoding="utf8", mode="w") as file_h:
        file_h.write(f"name\tvalue\n")
        stats(fact_bam_basic, file_h=file_h)


def expand_bam_files(bam_files):
    final_bam_files = []
    for bam_file in bam_files:
        if "*" in bam_file:
            final_bam_files.extend(glob(bam_file))
        else:
            final_bam_files.append(bam_file)
    return final_bam_files


def main_cli():
    """
    aligned bam analysis & origin bam analysis
    在 metric 中使用
    """
    polars_env_init()

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--bams", nargs="+", type=str,
                        required=True, help="wildcard '*' is supported")
    parser.add_argument("--refs", nargs="+", type=str,
                        help="if not provided. the alignment related metric will not output")
    parser.add_argument(
        "-f",
        action="store_true",
        default=False,
        help="regenerate the metric file if exists",
    )
    args = parser.parse_args()

    bam_files = args.bams
    bam_files = expand_bam_files(bam_files)

    refs = args.refs

    if refs is None or len(refs) == 0:
        refs = [""] * len(bam_files)

    if len(refs) == 1:
        refs = refs * len(bam_files)

    assert len(bam_files) == len(refs)

    for bam, ref in zip(bam_files, refs):
        main(bam_file=bam, ref_fa=ref, force=args.f)


if __name__ == "__main__":
    main_cli()
