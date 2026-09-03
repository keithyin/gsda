import pysam
from tqdm import tqdm
from argparse import Namespace
import argparse
import subprocess
import os
import logging

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s")


def extracted_wanted_channels(in_bam: str, ref_file: str, target_cov_thr: float = 0.8) -> set:
    """
    1. 调用 gsmm2 命令，生成一个比对 bam 文件， 可以通过 gsmm2 align -h 获取 该命令的调用信息
    2. 调用 bam_stat bam-concordance 得到channel粒度的比对结果, 最终会生成一个 csv 文件
    3. 读取该 csv 文件, 使用 polars，基于  (pl.col("refEnd") - pl.col("refStart")) / ref_len 获取 targetCoverage，然后仅保留 targetCoverage > 0.8 的 channel，将对应的 channel_id 输出
    """
    import polars as pl

    aligned_prefix = in_bam.rsplit(".", maxsplit=1)[0] + ".aligned"
    aligned_bam = f"{aligned_prefix}.bam"

    # if not os.path.exists(aligned_bam):
    cmd = f"gsmm2 --preset map-ont align -q {in_bam} -t {ref_file} -p {aligned_prefix} --noMar"
    logging.info(f"running: {cmd}")
    subprocess.check_call(cmd, shell=True)

    res_csv = f"{aligned_prefix}.metric.csv"
    # if not os.path.exists(res_csv):
    cmd = f"bam_stat bam-concordance {ref_file} {aligned_bam}"
    logging.info(f"running: {cmd}")
    subprocess.check_call(cmd, shell=True)

    if not os.path.exists(res_csv):
        logging.error("bam-concordance output %s not found", res_csv)
        return set()

    # read ref_len from fasta
    ref_len = 0
    with pysam.FastaFile(ref_file) as fa:
        ref_len = fa.lengths[0]

    df = pl.read_csv(res_csv, separator="\t")
    high_cov = (
        df.with_columns(
            ((pl.col("refEnd") - pl.col("refStart")) /
             pl.lit(ref_len)).alias("targetCoverage")
        )
        .filter(pl.col("targetCoverage") > target_cov_thr)
        .select("channel_id")
    )
    result = set(int(x) for x in high_cov.to_numpy()[:, 0])
    logging.info("Wanted %d channels (targetCoverage > %.2f)",
                 len(result), target_cov_thr)
    return result


def dump_sbr_bam(sbr_bam: str, infix: str, wanted_channels: set) -> str:
    """Filter SBR BAM by wanted_channels, write a new BAM and rebuild its index."""
    out_bam = sbr_bam.rsplit(".", maxsplit=1)[0] + f".{infix}.bam"

    with pysam.AlignmentFile(sbr_bam, "rb", check_sq=False, threads=os.cpu_count() // 2) as inf, \
            pysam.AlignmentFile(out_bam, "wb", template=inf, check_sq=False, threads=os.cpu_count() // 2) as outf:
        count = 0
        for rec in tqdm(inf.fetch(until_eof=True), desc=f"dump_sbr_bam -> {infix}"):
            ch = int(rec.query_name.split(":")
                     [-1]) if ":" in rec.query_name else -1
            if ch in wanted_channels:
                outf.write(rec)
                count += 1
    logging.info("Dumped %d reads to %s", count, out_bam)

    return out_bam


def main(args):
    channels = extracted_wanted_channels(
        args.smc_bam, args.ref_fa
    )
    dump_sbr_bam(args.sbr_bam, "filtered", channels)


def main_cli():

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--smc-bam", type=str, required=True)
    parser.add_argument("--sbr-bam", type=str, required=True)
    parser.add_argument("--ref-fa", type=str, required=True)
    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    cli_args = {
        "smc_bam": "/data1/ccs_data/20260428_250302Y0004_Run0002_polyN/20260428_250302Y0004_Run0002.smc_all_reads.bam",
        "sbr_bam": "/data1/ccs_data/20260428_250302Y0004_Run0002_polyN/20260428_250302Y0004_Run0002.smc_all_reads.bam",
        "ref_fa": "/data1/REF_GENOMES/poly-N-20260428_250302Y0004_Run0002.fasta",
        "infix": "filtered"
    }
    cli_args = Namespace(**cli_args)
    main(cli_args)
