import reads_quality_stats_v3
import subprocess
import pathlib
import os
import logging
import polars as pl
import shutil
import argparse
import semver
from glob import glob
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
print(cur_dir)
sys.path.append(cur_dir)

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


def gsetl_version_check():
    oup = subprocess.getoutput("gsetl -V")
    oup = oup.strip()
    version_str = oup.rsplit(" ", maxsplit=1)[1]

    logging.info(f"gsetl Version: {version_str}")
    gsetl_version = semver.Version.parse(version_str)
    expected_version = "0.8.0"
    assert gsetl_version >= semver.Version.parse(
        expected_version
    ), f"current gsetl version:{gsetl_version} < {expected_version}, try 'cargo uninstall gsetl; cargo install gsetl@={expected_version}' "


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def generate_non_aligned_metric_fact_file(bam_file: str, n: int, out_filepath: str, out_dir: str, force: bool, length_thr=None, length_percentile_thr=None):
    if not force and os.path.exists(out_filepath):
        logging.info(f"{out_filepath} exists, use the existed file")
        return

    cmd = f"gsetl --outdir {out_dir} non-aligned-bam-seq-n-stats --bam {bam_file} -n {n}"
    if length_thr is not None:
        cmd += f" --length-thr {length_thr}"
    if length_percentile_thr is not None:
        cmd += f" --length-percentile-thr {length_percentile_thr}"

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)


def stat_expr(name: str):
    name_items = name.split("-")
    first = name_items[0]
    name_items[0] = name_items[-1]
    name_items[-1] = first

    out_name_items = [name_items[-1]]
    out_name_items.extend(name_items[:-1])

    out_name = "-".join(out_name_items)
    return [
        pl.col(name).mean().alias(f"MEAN-{out_name}"),
        pl.col(name).median().alias(f"MEDIAN-{out_name}"),
        pl.col(name).std().alias(f"STD-{out_name}"),
    ]


def non_aligned_metric_analysis(fact_metric_filename: str, aggr_metric_filename: str, force: bool):
    if os.path.exists(aggr_metric_filename) and not force:
        logging.warning(f"{aggr_metric_filename} will be override")

    df = pl.read_csv(
        fact_metric_filename, separator="\t",
        infer_schema_length=3000)

    df = df.select([
        pl.col("dw-first-n-median") * 2,
        pl.col("dw-last-n-median") * 2,
        pl.col("ar-first-n-median") * 2,
        pl.col("ar-last-n-median") * 2,
        pl.col("dw-first-n-mean") * 2,
        pl.col("dw-last-n-mean") * 2,
        pl.col("ar-first-n-mean") * 2,
        pl.col("ar-last-n-mean") * 2,
    ])

    exprs = stat_expr("dw-first-n-median")
    exprs.extend(stat_expr("dw-last-n-median"))
    exprs.extend(stat_expr("ar-first-n-median"))
    exprs.extend(stat_expr("ar-last-n-median"))

    exprs.extend(stat_expr("dw-first-n-mean"))
    exprs.extend(stat_expr("dw-last-n-mean"))
    exprs.extend(stat_expr("ar-first-n-mean"))
    exprs.extend(stat_expr("ar-last-n-mean"))

    df = df.select(exprs)

    metrics = df.transpose(
        include_header=True, header_name="name", column_names=["value"]
    )
    metrics = metrics.sort("name")

    print(metrics)

    metrics.write_csv(aggr_metric_filename,
                      include_header=True, separator="\t")


def main(
    bam_file: str,
    n: int,
    length_thr=None,
    length_percentile_thr=None,

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

    seq_n_metric_filename = f"{outdir}/{stem}.seq_n_stats.csv"
    if force and os.path.exists(seq_n_metric_filename):
        os.remove(seq_n_metric_filename)

    generate_non_aligned_metric_fact_file(
        bam_file, n, out_filepath=seq_n_metric_filename, out_dir=outdir,
        force=force, length_thr=length_thr, length_percentile_thr=length_percentile_thr)

    seq_n_metric_aggr_filename = f"{outdir}/{stem}.seq_n_stats.aggr.csv"
    non_aligned_metric_analysis(
        seq_n_metric_filename, aggr_metric_filename=seq_n_metric_aggr_filename, force=True)
    return outdir, stem


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
    parser.add_argument("-n", required=True,
                        type=int)
    parser.add_argument("--length-thr", default=None,
                        type=int, dest="length_thr")
    parser.add_argument("--length-percentile-thr", default=None, type=int,
                        help="[0, 100], compute the length-thr according to the length-percentile-thr", dest="length_percentile_thr")
    parser.add_argument("--ref-fa", default=None, type=str,
                        help="ref fasta", dest="ref_fa")
    parser.add_argument(
        "-f",
        action="store_true",
        default=False,
        help="regenerate the metric file if exists",
    )
    args = parser.parse_args()

    assert args.length_thr is not None or args.length_percentile_thr is not None, "--length-thr and --length-percentile-thr can't all be None"

    bam_files = args.bams
    bam_files = expand_bam_files(bam_files)

    for bam in bam_files:
        main(bam_file=bam, n=args.n, force=args.f, length_thr=args.length_thr,
                            length_percentile_thr=args.length_percentile_thr)


if __name__ == "__main__":
    main_cli()
