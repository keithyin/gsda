import pysam
import matplotlib.pyplot as plt
import os
from multiprocessing import cpu_count
import logging
import subprocess
import pathlib
from tqdm import tqdm
from glob import glob

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def get_ref_seq(fasta_file):
    with pysam.FastxFile(fasta_file) as f:
        for record in f:
            return record.sequence


def smc_bam_basic_stat(bam: str):
    q8 = 0
    q10 = 0
    q20 = 0
    q30 = 0
    with pysam.AlignmentFile(bam, "rb", threads=cpu_count(), check_sq=False) as file:
        for record in tqdm(file.fetch(until_eof=True), desc=f"reading {bam}"):
            rq = float(record.get_tag("rq"))
            if rq > 0.8431:
                q8 += 1
            if rq > 0.9:
                q10 += 1
            if rq > 0.99:
                q20 += 1
            if rq > 0.999:
                q30 += 1

    print(f"q8:{q8}, q10:{q10}, q20:{q20}, q30:{q30}")


def generate_aligned_metric_fact_file(
    bam_file: str,
    ref_fasta: str,
    short_aln=False,
    force: bool = False,
    threads=None
) -> str:

    p = pathlib.Path(bam_file)
    align_out_prefix = p.with_suffix(".align")
    for fpath in glob(f"{align_out_prefix}*"):
        os.remove(fpath)

    threads = cpu_count() if threads is None else threads
    cmd = f"gsmm2 --threads {threads} align -q {bam_file} -t {ref_fasta} -p {str(align_out_prefix)} --kmer 11 --wins 1 --rq-range 0.99:1.1"
    if short_aln:
        cmd += "  --short-aln"

    logging.info("cmd: %s", cmd)
    subprocess.check_call(cmd, shell=True)

    return f"{str(align_out_prefix)}.bam"


def collect_qscores_from_locus(
    bam_path,
    ref_name,
    ref_start,
    ref_end,
    ref_seq
):
    """
    从 BAM 中提取指定 ref 区间的 query base Q 值
    """

    bam = pysam.AlignmentFile(bam_path, "rb", threads=cpu_count())

    qscores_true = []
    qscores_false = []

    # 🚀 只遍历目标区域（关键优化）
    for read in tqdm(bam.fetch(ref_name, ref_start, ref_end), desc=f"processing {bam_path}"):

        if read.is_unmapped:
            continue
        # if  read.is_secondary or read.is_supplementary:
        #     continue

        query_qual = read.query_qualities
        query_seq = read.query_sequence

        # (qpos, rpos)
        for qpos, rpos in read.get_aligned_pairs(matches_only=False):

            # 跳过 indel
            if qpos is None or rpos is None:
                continue

            # 限定在目标 locus
            if rpos < ref_start or rpos >= ref_end:
                continue
            if ref_seq[rpos] == query_seq[qpos]:
                qscores_true.append(query_qual[qpos])
            else:
                qscores_false.append(query_qual[qpos])

    bam.close()
    return qscores_true, qscores_false


def plot_qscore_distribution(qscores, ref_start, ref_end, bins=100):
    plt.figure()
    plt.hist(qscores[0], bins=bins, label="=")
    plt.hist(qscores[1], bins=bins, label="X")
    plt.legend()
    plt.xlabel("Q score")
    plt.ylabel("Frequency")
    plt.title(f"Q score distribution in locus. [{ref_start}, {ref_end})")
    plt.show()


def main():
    """
        主要关注 比对到 参考基因组 某个 位点的 reads 的 Q 值分布
    """
    bam_path = "/data1/ccs_data/20260318-saisuofei-snp-check/barcodes_reads_fastq_amplicon/Group_0_Adaptor-barcode295-2.bam"
    ref_fasta = "/data1/ccs_data/20260318-saisuofei-snp-check/barcodes_reads_plasmid_asm_amplicon/plasmid_assembly/Group_0_Adaptor-barcode295-2_plassembler_plasmids.fasta"
    smc_bam_basic_stat(bam_path)

    aligned_bam_path = generate_aligned_metric_fact_file(
        bam_path, ref_fasta=ref_fasta)
    # ref_name = "Clust_0.1.1|len:1333|num:3824|ratio:0.76|complete:Yes|sus:6"
    ref_name = "1"

    ref_seq = get_ref_seq(ref_fasta)

    # ref_start = 576
    # ref_end = 577

    ref_start = 4994
    ref_end = 4995
    expanded_ref_start = max(ref_start-5, 0)
    expanded_ref_end = min(ref_end + 5, len(ref_seq))
    print(
        f":{ref_seq[expanded_ref_start: ref_start]}[{ref_seq[ref_start: ref_end]}]{ref_seq[ref_end:expanded_ref_end]}")
    # os._exit(0)
    qscores = collect_qscores_from_locus(
        aligned_bam_path,
        ref_name,
        ref_start,
        ref_end,
        ref_seq=ref_seq
    )

    print(f"Collected {len(qscores[0]) + len(qscores[1])} bases")

    plot_qscore_distribution(qscores, ref_start=ref_start, ref_end=ref_end)


if __name__ == "__main__":
    main()
