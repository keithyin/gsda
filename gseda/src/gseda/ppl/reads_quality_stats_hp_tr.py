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
    rq_range=None,
) -> str:

    if no_supp and no_mar:
        raise ValueError("no_supp, no_mar can't be all true")

    if force and os.path.exists(out_filename):
        os.remove(out_filename)

    if os.path.exists(out_filename):
        logging.warning("use the existing metric file : %s", out_filename)
        return out_filename

    threads = cpu_count() if threads is None else threads
    cmd = f"""gsmm2-metric --mode hp-tr-v2 --threads {threads} \
            -q {bam_file} \
            -t {ref_fasta} \
            --out {out_filename} \
            --kmer 11 \
            --wins 1 """
    if short_aln:
        cmd += " --short-aln"

    if np_range is not None:
        cmd += f" --np-range {np_range}"
    if rq_range is not None:
        cmd += f" --rq-range {rq_range}"

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    return out_filename


def stats(metric_filename, filename):
    pattern = r"\(([^)]+)\)(\d+)"
    # r"\(([^)]+)\)(\d+)"
    df = pl.read_csv(
        metric_filename, separator="\t", schema_overrides={"longIndel": pl.String}
    )
    df = df.filter(pl.col("rname") != "")

    df = (df
          .group_by(["motif"])
          .agg([
              pl.col("eq").sum(),
              pl.col("diff").sum(),
              pl.col("ins").sum(),
              pl.col("del").sum()])
          .with_columns([
              (pl.col("eq") + pl.col("diff") + pl.col("ins") +
               pl.col("del")).alias("aligned_span")
          ])
          .with_columns([
              (pl.col("eq") / pl.col("aligned_span")).alias("eq_rate"),
              (pl.col("diff") / pl.col("aligned_span")).alias("diff_rate"),
              (pl.col("ins") / pl.col("aligned_span")).alias("ins_rate"),
              (pl.col("del") / pl.col("aligned_span")).alias("del_rate"),
          ])
          .with_columns([
              pl.col("motif").str.extract(pattern, 1).alias("key_text"),
              pl.col("motif").str.extract(pattern, 2).cast(pl.Int64).alias("key_cnt")])
          .sort([pl.col("key_text").str.len_chars(), pl.col("key_text"), pl.col("key_cnt")], descending=[False, False, False])
          .drop(["key_text", "key_cnt"])
          )
    print(df)

    df.write_csv(filename, include_header=True, separator="\t")


def main(
    bam_file: str,
    ref_fa: str,
    threads=None,
    force=False,
    short_aln=False,
    outdir=None,
    copy_bam_file=False,
    np_range=None,
    rq_range=None,
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
        "gsmm2-metric", semver.Version.parse("0.4.1"), "cargo install gsmm2-metric")

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

    fact_metric_filename = f"{outdir}/{stem}.gsmm2-hp-tr-fact.csv"
    fact_metric_filename = generate_metric_file(
        bam_file,
        ref_fa,
        out_filename=fact_metric_filename,
        force=force,
        threads=threads,
        short_aln=short_aln,
        np_range=np_range,
        rq_range=rq_range,
    )
    aggr_metric_filename = f"{outdir}/{stem}.gsmm2-hp-tr-aggr.csv"
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


def main_cli():
    polars_env_init()

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--bams", nargs="+", type=str, required=True)
    parser.add_argument("--refs", nargs="+", type=str, required=True)
    parser.add_argument("--np-range", type=str, default=None, dest="np_range")
    parser.add_argument("--rq-range", type=str, default=None, dest="rq_range")
    parser.add_argument("--short-aln", type=int, default=0,
                        help="for query or target in [30, 200]", dest="short_aln")
    parser.add_argument(
        "-f",
        action="store_true",
        default=False,
        help="regenerate the metric file if exists",
    )
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
