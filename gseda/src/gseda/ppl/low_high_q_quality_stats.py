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


def dump_sub_bam(bam_path: str, out_path: str, q_start: float, q_end: float):
    try:
        with pysam.AlignmentFile(bam_path, mode="rb", check_sq=False, threads=cpu_count() // 2) as bam_in:
            with pysam.AlignmentFile(
                    out_path, mode="wb", check_sq=False, header=bam_in.header, threads=cpu_count() // 2) as bam_out:
                for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {bam_path} for dump to {out_path}"):
                    cq = read.get_tag("cq")
                    if cq >= q_start and cq <= q_end:
                        bam_out.write(read=read)
    except Exception as e:
        logging.error(f"dump sub bam error : {e}")
        


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
    q_ranges: str,
    ref_fa: str,
    force=False,
    outdir=None,
) -> str:

    
    bam_file_dir = os.path.dirname(bam_file)
    stem = pathlib.Path(bam_file).stem
    if outdir is None:
        outdir = os.path.join(bam_file_dir, f"{stem}-metric")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for q_range in q_ranges.split(","):
        q_range_items = q_range.split(":")
        assert len(q_range_items) == 2
        q_start, q_end = float(q_range_items[0]), float(q_range_items[1])
        bam_name  =os.path.join(outdir, f"{stem}.Q{q_start}-{q_end}.bam") 

        dump_sub_bam(bam_file, bam_name, q_start, q_end)

        reads_quality_stats_v3.main(
            bam_file=bam_name, ref_fa=ref_fa, force=force)



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
    parser.add_argument("--q-ranges", required=True,
                        type=str, help="start1:end1,start2:end2", dest="q_ranges")
    parser.add_argument("--ref-fa", default="", type=str,
                        help="ref fasta", dest="ref_fa")
    parser.add_argument(
        "-f",
        action="store_true",
        default=False,
        help="regenerate the metric file if exists",
    )
    args = parser.parse_args()

    ref_fa = args.ref_fa

    bam_files = args.bams
    bam_files = expand_bam_files(bam_files)

    for bam in bam_files:
        main(bam_file=bam, q_ranges=args.q_ranges.strip(), ref_fa=ref_fa, force=args.f)


if __name__ == "__main__":
    main_cli()
