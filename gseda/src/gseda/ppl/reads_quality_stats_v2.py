import subprocess
import pathlib
import os
import logging
import polars as pl
import shutil
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
    outdir: str,
    force: bool = False,
    threads=None,
    no_supp=False,
    no_mar=False,
) -> str:

    if no_supp and no_mar:
        raise ValueError("no_supp, no_mar can't be all true")

    metric_filename = f"{outdir}/metric.csv"
    if no_supp:
        metric_filename = f"{outdir}/metric_noSupp.csv"
    if no_mar:
        metric_filename = f"{outdir}/metric_noMar.csv"

    if force and os.path.exists(metric_filename):
        os.remove(metric_filename)

    if os.path.exists(metric_filename):
        logging.warn("use the existing metric file : %s", metric_filename)
        return metric_filename

    threads = cpu_count() if threads is None else threads
    cmd = f"""gsmm2-aligned-metric --threads {threads} \
            -q {bam_file} \
            -t {ref_fasta} \
            --oupdir {outdir} \
            --kmer 11 \
            --wins 1 """

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    return metric_filename


def stats(metric_filename, filename):
    df = pl.read_csv(metric_filename, separator="\t").filter(pl.col("rname") != "")
    # print(df.head(2))
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
                "ovlp",
                "mergedQrySpan",
            ]
        )
    )

    aggr_metrics = df.select(
        [
            (pl.col("covlen").sum() / pl.col("qlen").sum()).alias("queryCoverage"),
            (
                pl.col("match").sum()
                / (
                    pl.col("match")
                    + pl.col("misMatch")
                    + pl.col("ins")
                    + pl.col("homoIns")
                    + pl.col("del")
                    + pl.col("homoDel")
                ).sum()
            ).alias("identity"),
        ]
    )
    aggr_metrics = aggr_metrics.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )
    print(aggr_metrics)

    aggr_metrics.write_csv(filename, include_header=True, separator="\t")


def main(bam_file: str, ref_fa: str, threads=None, force=False, outdir=None) -> str:
    """
        step1: do alignment
        step2: generate detailed metric info
        step3: compute the aggr metric. the result aggr_metric.csv is a '\t' seperated csv file. the header is name\tvalue
            here is a demo.
            ---------aggr_metric.csv
            name    value
            queryCoverage   0.937
            ----------

    requirements:
        mm2: cargo install mm2
        gsetl: cargo install gsetl

    Args:
        bam_file (str): bam file. only support adapter.bam
        ref_fa (str): ref genome fa file nam
        force (boolean): if force==False, the outdir must not exists in advance. if force==True, the outdir will be removed if exists
            the proceduer will create a empty outdir for the metric related files
        outdir:
            if outdir provided, read ${outdir}/metric/aggr_metric.csv for metric result
            if not, read ${bam_file_dir}/${bam_file_name}-metric/metric/aggr_metric.csv for metric result

    Return:
        aggr_metric_filename (str): the aggr metric file
    """
    bam_filedir = os.path.dirname(bam_file)
    bam_filename = extract_filename(bam_file)
    if outdir is None:
        outdir = os.path.join(bam_filedir, f"{bam_filename}-metric")
    if force and os.path.exists(outdir):
        shutil.rmtree(outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    metric_filename = generate_metric_file(
        bam_file, ref_fa, outdir=outdir, force=force, threads=threads
    )
    aggr_metric_filename = f"{outdir}/aggr_metric.csv"
    stats(metric_filename, filename=aggr_metric_filename)
    return aggr_metric_filename


def test_stat():
    fact_bam_basic = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/fact_aligned_bam_bam_basic.csv"
    aggr_metric_filename = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/aggr_metric.csv"

    with open(aggr_metric_filename, encoding="utf8", mode="w") as file_h:
        file_h.write(f"name\tvalue\n")
        stats(fact_bam_basic, file_h=file_h)


if __name__ == "__main__":
    polars_env_init()
    bam_files = [
        "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter.bam",
        # "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0002_adapter.bam",
        # "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0003_adapter.bam",
    ]
    ref = "/data/ccs_data/MG1655.fa"

    for bam in bam_files:
        main(bam_file=bam, ref_fa=ref, force=False)
    # test_stat()

    # print(merge_intervals([{"qstart": 11, "qend": 771}]))
