import subprocess
import pathlib
import os
import logging

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

    cmd = f"""gsmm2 align -q {bam_file} -t {ref_fasta} -p {res_bam_prefix} --noSeco --noSupp"""

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    return result_bam


def generate_fact_table(aligned_bam: str, outdir: str):
    


def main(bam_file: str, ref_fa: str):
    """
        step1: do alignment
        step2: generate detailed metric info
        step3: compute the aggr metric

    requirements:
        mm2: cargo install mm2
        gsetl: cargo install gsetl

    Args:
        bam_file (str): bam file. only support adapter.bam and smc.bam now
        ref_fa (str): ref genome fa file nam
    """
