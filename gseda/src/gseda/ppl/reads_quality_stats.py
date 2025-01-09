import subprocess
import pathlib
import os
import logging
import polars as pl
import shutil
from multiprocessing import cpu_count

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def do_alignment(
    bam_file: str, ref_fasta: str, outdir: str, force: bool = False, threads=None
) -> str:

    res_bam_prefix = "{}/{}.aligned".format(outdir, extract_filename(bam_file))
    result_bam = f"{res_bam_prefix}.bam"

    if force and os.path.exists(result_bam):
        os.remove(result_bam)

    if os.path.exists(result_bam):
        logging.info(f"{result_bam} exists, use it now")
        return result_bam
    threads = cpu_count() if threads is None else threads
    cmd = f"""gsmm2 --threads {threads} align \
            -q {bam_file} \
            -t {ref_fasta} \
            -p {res_bam_prefix} \
            --kmer 11 \
            --wins 1 \
            --noSeco"""

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
    cmd = f"""gsetl -f --outdir {outdir} \
        aligned-bam \
        --bam {aligned_bam} \
        --ref-file {ref_fasta} \
        --factRefLocusInfo 0 \
        --factErrorQueryLocusInfo 0 \
        --factBaseQStat 0 \
        --factPolyInfo 0 \
        --factRecordStat 0 \
        --useSupp
        """
    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    fact_bam_basic = f"{outdir}/fact_aligned_bam_bam_basic.csv"
    # fact_bam_record_stat = f"{outdir}/fact_aligned_bam_record_stat.csv"
    return fact_bam_basic, None


# def stats(fact_bam_basic, fact_bam_record_stat, file_h):
#     fact_bam_basic = pl.read_csv(fact_bam_basic, separator="\t")
#     fact_bam_record_stat = pl.read_csv(fact_bam_record_stat, separator="\t")
#     # print(fact_bam_basic.head(2))
#     # print(fact_bam_record_stat.head(2))

#     joined = fact_bam_basic.join(fact_bam_record_stat, on="qname")
#     # print(joined.head(2))
#     query_coverage = joined.select(
#         [
#             pl.col("qlen"),
#             (
#                 pl.col("matchBp")
#                 + pl.col("mismatchBp")
#                 + pl.col("nonHpInsertionBp")
#                 + pl.col("hpInsertionBp")
#             ).alias("alignedBp"),
#         ]
#     ).select(
#         [(pl.col("alignedBp").sum() / pl.col("qlen").sum()).alias("queryCoverage")]
#     )
#     qc = query_coverage.to_numpy()[0][0]

#     file_h.write(f"queryCoverage\t{qc}\n")
#     # print(query_coverage)


def stats(fact_bam_basic, file_h):
    fact_bam_basic = pl.read_csv(fact_bam_basic, separator="\t")
    df = (
        fact_bam_basic.select(
            [
                pl.col("qname"),
                pl.col("qlen"),
                pl.when(pl.col("fwd"))
                .then(pl.col("qstart"))
                .otherwise(pl.col("qlen") - pl.col("qend"))
                .alias("qstart"),
                pl.when(pl.col("fwd"))
                .then(pl.col("qend"))
                .otherwise(pl.col("qlen") - pl.col("qstart"))
                .alias("qend"),
            ]
        )
        .select(
            [pl.col("qname"), pl.col("qlen"), pl.struct("qstart", "qend").alias("se")]
        )
        .group_by("qname")
        .agg([pl.col("qlen").max(), pl.col("se"), pl.len().alias("cnt")])
    )

    # print(df.filter(pl.col("cnt") > 1).count())
    # os._exit(0)

    df = df.with_columns(
        [
            pl.col("se")
            .map_elements(
                lambda x: merge_intervals(x),
                return_dtype=pl.UInt64,
                returns_scalar=True,
            )
            .alias("aligned_qlen")
        ]
    )

    # print(df.filter(pl.col("aligned_qlen") < 0).head(2))
    # print(df.head(2))

    df = df.select(
        [
            (pl.col("aligned_qlen").sum() / pl.col("qlen").cast(pl.UInt64).sum()).alias(
                "queryCoverage"
            )
        ]
    )
    print(df)
    qc = df.to_numpy()[0][0]
    file_h.write(f"queryCoverage\t{qc}\n")


def merge_intervals(region_list):
    # print(region_list)
    region_list = [[se["qstart"], se["qend"]] for se in region_list]
    if len(region_list) == 1:
        return region_list[0][1] - region_list[0][0]

    # Step 1: Sort intervals by start position
    sorted_regions = sorted(region_list, key=lambda x: x[0])

    # Step 2: Merge overlapping intervals
    merged_regions = []
    current_start, current_end = sorted_regions[0][0], sorted_regions[0][1]

    for start, end in sorted_regions[1:]:
        if start <= current_end:  # Overlapping or adjacent regions
            current_end = max(current_end, end)  # Extend the current region
        else:
            # No overlap, so save the current region and start a new one
            merged_regions.append((current_start, current_end))
            current_start, current_end = start, end

    # Add the last region
    merged_regions.append((current_start, current_end))

    # Step 3: Calculate the total length of merged intervals
    total_length = sum([(se[1] - se[0]) for se in merged_regions])
    # print(merged_regions, total_length)
    # if total_length < 0:
    #     print(merge_intervals, total_length)

    return total_length


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
    if os.path.exists(outdir):
        raise ValueError(
            f"{outdir} already exists. remove it by manually or force=True"
        )

    os.makedirs(outdir)

    metric_outdir = os.path.join(outdir, "metric")

    aligned_bam_file = do_alignment(
        bam_file, ref_fa, outdir=outdir, force=force, threads=threads
    )

    fact_bam_basic, _ = generate_fact_table(
        aligned_bam_file, ref_fasta=ref_fa, outdir=metric_outdir
    )

    aggr_metric_filename = os.path.join(metric_outdir, "aggr_metric.csv")

    with open(aggr_metric_filename, encoding="utf8", mode="w") as file_h:
        file_h.write(f"name\tvalue\n")
        stats(fact_bam_basic, file_h=file_h)
    return aggr_metric_filename


def test_stat():
    fact_bam_basic = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/fact_aligned_bam_bam_basic.csv"
    aggr_metric_filename = "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter-metric/metric/aggr_metric.csv"

    with open(aggr_metric_filename, encoding="utf8", mode="w") as file_h:
        file_h.write(f"name\tvalue\n")
        stats(fact_bam_basic, file_h=file_h)


if __name__ == "__main__":
    bam_files = [
        "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0001_adapter.bam",
        "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0002_adapter.bam",
        "/data/adapter-query-coverage-valid-data/20250107_240901Y0007_Run0003_adapter.bam",
    ]
    ref = "/data/ccs_data/MG1655.fa"

    for bam in bam_files:
        main(bam_file=bam, ref_fa=ref, force=True)
    # test_stat()

    # print(merge_intervals([{"qstart": 11, "qend": 771}]))
