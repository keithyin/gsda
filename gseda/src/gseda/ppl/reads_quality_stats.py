import subprocess
import pathlib
import os
import logging
import polars as pl
import shutil

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def do_alignment(
    bam_file: str, ref_fasta: str, outdir: str, force: bool = False
) -> str:

    res_bam_prefix = "{}/{}.aligned".format(outdir, extract_filename(bam_file))
    result_bam = f"{res_bam_prefix}.bam"

    if force and os.path.exists(result_bam):
        os.remove(result_bam)

    if os.path.exists(result_bam):
        logging.info(f"{result_bam} exists, use it now")
        return result_bam

    cmd = f"""gsmm2 align -q {bam_file} -t {ref_fasta} -p {res_bam_prefix} --noSeco --noSupp"""

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    return result_bam


def generate_fact_table(aligned_bam: str, ref_fasta: str, outdir: str):
    """_summary_

    Args:
        aligned_bam (str): aligned bam
        outdir (str): the outdir must be empty or not in use
    """
    if os.path.exists(outdir):
        raise ValueError(f"{outdir} already exists")
    cmd = f"gsetl -f --outdir {outdir} aligned-bam --bam {aligned_bam} --ref-file {ref_fasta} --factRefLocusInfo 0 --factErrorQueryLocusInfo 0 --factBaseQStat 0 --factPolyInfo 0"
    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    fact_bam_basic = f"{outdir}/fact_aligned_bam_bam_basic.csv"
    fact_bam_record_stat = f"{outdir}/fact_aligned_bam_record_stat.csv"
    return fact_bam_basic, fact_bam_record_stat


def stats(fact_bam_basic, fact_bam_record_stat, file_h):
    fact_bam_basic = pl.read_csv(fact_bam_basic, separator="\t")
    fact_bam_record_stat = pl.read_csv(fact_bam_record_stat, separator="\t")
    # print(fact_bam_basic.head(2))
    # print(fact_bam_record_stat.head(2))

    joined = fact_bam_basic.join(fact_bam_record_stat, on="qname")
    # print(joined.head(2))
    query_coverage = joined.select(
        [
            pl.col("qlen"),
            (
                pl.col("matchBp")
                + pl.col("mismatchBp")
                + pl.col("nonHpInsertionBp")
                + pl.col("hpInsertionBp")
            ).alias("alignedBp"),
        ]
    ).select(
        [(pl.col("alignedBp").sum() / pl.col("qlen").sum()).alias("queryCoverage")]
    )
    qc = query_coverage.to_numpy()[0][0]

    file_h.write(f"queryCoverage\t{qc}\n")
    # print(query_coverage)


def main(bam_file: str, ref_fa: str, force=False, outdir=None) -> str:
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
    if os.path.exists(outdir):
        raise ValueError(f"{outdir} already exists. remove it by manually or force=True")
    
    os.makedirs(outdir)

    metric_outdir = os.path.join(outdir, "metric")

    aligned_bam_file = do_alignment(bam_file, ref_fa, outdir=outdir, force=force)

    fact_bam_basic, fact_bam_record_stat = generate_fact_table(
        aligned_bam_file, ref_fasta=ref_fa, outdir=metric_outdir
    )

    aggr_metric_filename = os.path.join(metric_outdir, "aggr_metric.csv")

    with open(aggr_metric_filename, encoding="utf8", mode="w") as file_h:
        file_h.write(f"name\tvalue\n")
        stats(fact_bam_basic, fact_bam_record_stat, file_h=file_h)
    return aggr_metric_filename


def test_stat():
    fact_bam_basic = "/data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter-metric/metric/fact_aligned_bam_bam_basic.csv"
    fact_bam_record_stat = "/data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter-metric/metric/fact_aligned_bam_record_stat.csv"
    aggr_metric_filename = "/data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter-metric/metric/aggr_metric.csv"

    with open(aggr_metric_filename, encoding="utf8", mode="w") as file_h:
        file_h.write(f"name\tvalue\n")
        stats(fact_bam_basic, fact_bam_record_stat, file_h=file_h)


if __name__ == "__main__":
    bam = "/data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter.bam"
    ref = "/data/ccs_data/MG1655.fa"
    main(bam_file=bam, ref_fa=ref)
    # test_stat()
