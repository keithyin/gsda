import pysam
import pathlib
import os
import logging
import argparse
from glob import glob
import sys
import polars as pl
from multiprocessing import cpu_count
from tqdm import tqdm

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
print(cur_dir)
sys.path.append(cur_dir)
import reads_quality_stats_v3  # noqa: E402


# deprecated ...
logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def is_empty_file(filepath: str):
    valid_line_cnt = 0
    with open(filepath, mode="r", encoding="utf8") as f:
        for line in f:
            if line.strip() != "":
                valid_line_cnt += 1
            if valid_line_cnt > 1:
                return False
    return True


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def compute_percentile_length(bam_path: str, percentile: int) -> int:
    logging.info("computing length thr")
    lengths = []
    with pysam.AlignmentFile(filename=bam_path, mode="rb", check_sq=False, threads=cpu_count()) as bam_in:
        for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {bam_path}"):
            lengths.append(read.query_length)
    lengths = sorted(lengths)
    assert len(lengths) > 0
    if percentile <= 0:
        return lengths[0]
    if percentile >= 100:
        return lengths[-1]

    pos = int(len(lengths) * (percentile / 100))

    return lengths[pos]


def dump_sub_bam(bam_path: str, out_path: str, length_thr: int, n: int, first=True):
    with pysam.AlignmentFile(bam_path, mode="rb", check_sq=False, threads=cpu_count() // 2) as bam_in:
        with pysam.AlignmentFile(
                out_path, mode="wb", check_sq=False, header=bam_in.header, threads=cpu_count() // 2) as bam_out:
            for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {bam_path} for dump to {out_path}"):
                if read.query_length < length_thr:
                    continue
                n = min(n, read.query_length)

                seq_len = read.query_length

                if first:
                    read.query_sequence = read.query_sequence[0:n]
                    dw = read.get_tag("dw")
                    ar = read.get_tag("ar")
                    cr = read.get_tag("cr")
                    read.set_tag("dw", dw[0:n])
                    read.set_tag("ar", ar[0:n])
                    read.set_tag("cr", cr[0:n])
                else:
                    read.query_sequence = read.query_sequence[(seq_len - n):]
                    dw = read.get_tag("dw")
                    ar = read.get_tag("ar")
                    cr = read.get_tag("cr")

                    read.set_tag("dw", dw[(seq_len - n):])
                    read.set_tag("ar", ar[(seq_len - n):])
                    read.set_tag("cr", cr[(seq_len - n):])

                bam_out.write(read=read)


def concat_metrics(non_aligned_metric_filepath, aligned_metric_filepath, new_value_name: str) -> pl.DataFrame:
    concated_metrics = []
    if not is_empty_file(non_aligned_metric_filepath):
        concated_metrics.append(pl.read_csv(non_aligned_metric_filepath, separator="\t", schema_overrides={
            "name": pl.String, "value": pl.Float32}))
    if not is_empty_file(aligned_metric_filepath):
        concated_metrics.append(pl.read_csv(aligned_metric_filepath, separator="\t", schema_overrides={
            "name": pl.String, "value": pl.Float32}))
    concated_metrics = pl.concat(concated_metrics)
    assert isinstance(concated_metrics, pl.DataFrame)
    concated_metrics = concated_metrics.rename({"value": new_value_name})
    return concated_metrics


def main(
    bam_file: str,
    n: int,
    ref_fa: str,
    length_thr=None,
    length_percentile_thr=None,
    force=False,
    outdir=None,
) -> str:

    if length_thr is None:
        length_thr = compute_percentile_length(
            bam_file, percentile=length_percentile_thr)
    logging.info(f"length thr = {length_thr}")
    bam_file_dir = os.path.dirname(bam_file)
    stem = pathlib.Path(bam_file).stem
    if outdir is None:
        outdir = os.path.join(bam_file_dir, f"{stem}-metric")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    first_n_bam = os.path.join(outdir, f"{stem}.first-n.bam")
    last_n_bam = os.path.join(outdir, f"{stem}.last-n.bam")

    final_res_csv = os.path.join(outdir, f"{stem}.seq-n-stats-aggr.csv")

    dump_sub_bam(bam_file, first_n_bam, length_thr=length_thr, n=n, first=True)
    dump_sub_bam(bam_file, last_n_bam, length_thr=length_thr, n=n, first=False)

    (first_n_aligned_aggr, _, first_n_non_aligned_aggr, _) = reads_quality_stats_v3.main(
        bam_file=first_n_bam, ref_fa=ref_fa, force=force)
    (last_n_aligned_aggr, _, last_n_non_aligned_aggr, _) = reads_quality_stats_v3.main(
        bam_file=last_n_bam, ref_fa=ref_fa, force=force)

    first_n_concated_metrics = concat_metrics(
        first_n_non_aligned_aggr, first_n_aligned_aggr, new_value_name="value_first_n")
    last_n_concated_metrics = concat_metrics(
        last_n_non_aligned_aggr, last_n_aligned_aggr, new_value_name="value_last_n")

    metric = first_n_concated_metrics.join(
        last_n_concated_metrics, on="name", how="outer")

    metric = metric.select([
        pl.coalesce(pl.col("name"), pl.col("name_right")).alias("name"),
        pl.col("value_first_n"),
        pl.col("value_last_n"),
    ])

    print(metric)

    metric.write_csv(final_res_csv, separator="\t")

    logging.info(f"result dumped to {final_res_csv}")


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

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--bams", nargs="+", type=str,
                        required=True, help="wildcard '*' is supported")
    parser.add_argument("-n", required=True,
                        type=int)
    parser.add_argument("--length-thr", default=None,
                        type=int, dest="length_thr")
    parser.add_argument("--length-percentile-thr", default=None, type=int,
                        help="[0, 100], compute the length-thr according to the length-percentile-thr", dest="length_percentile_thr")
    parser.add_argument("--ref-fa", default="", type=str,
                        help="ref fasta", dest="ref_fa")
    parser.add_argument(
        "-f",
        action="store_true",
        default=False,
        help="regenerate the metric file if exists",
    )
    args = parser.parse_args()

    assert args.length_thr is not None or args.length_percentile_thr is not None, "--length-thr and --length-percentile-thr can't all be None"

    ref_fa = args.ref_fa

    bam_files = args.bams
    bam_files = expand_bam_files(bam_files)

    for bam in bam_files:
        main(bam_file=bam, n=args.n, force=args.f, length_thr=args.length_thr,
             length_percentile_thr=args.length_percentile_thr, ref_fa=ref_fa)


if __name__ == "__main__":
    main_cli()
