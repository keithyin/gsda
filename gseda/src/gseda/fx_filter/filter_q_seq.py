from typing import Mapping
import pysam
from pathlib import Path
from Bio import SeqIO


def read_record_rq(bam_file: str) -> Mapping[str, float]:
    """
    从 bam 文件中读取每条 read 的 query_name 和其 rq 值。

    Args:
        bam_file (str): BAM 文件路径

    Returns:
        Mapping[str, float]: qname 到 rq 值的映射
    """
    qname_rq = {}
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40) as bam:
        for read in bam:
            try:
                rq = read.get_tag("rq")
                qname_rq[read.query_name] = float(rq)
            except KeyError:
                continue  # 如果没有 rq tag，跳过
    return qname_rq


def filter_fasta(
    fa_filepath: str, qname_rq_mapping: Mapping[str, float], rq_threshold: float
):
    """
    读取 fasta 文件，筛选出在 qname_rq_mapping 中 rq 值 ≥ 阈值的序列，并写入新的 fasta 文件

    Args:
        fa_filepath (str): FASTA 文件路径
        qname_rq_mapping (Mapping[str, float]): qname 到 rq 的映射
        rq_threshold (float): 阈值，保留 rq ≥ 该值的序列
    """
    out_fpath = "{}.filtered.fasta".format(fa_filepath.rsplit(".", maxsplit=1)[0])

    with open(fa_filepath) as fin, open(out_fpath, "w") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            rq = qname_rq_mapping.get(record.id)
            if rq is not None and rq >= rq_threshold:
                SeqIO.write(record, fout, "fasta")

    print(f"Filtered FASTA saved to: {out_fpath}")


def main_cli():
    bam_file = "/data/ccs_data/ludaopei-false-positive/20250414_240601Y0005_Run0001.smc_all_reads.bam"
    fa_filepath ="/data/ccs_data/ludaopei-false-positive/0414-demux/s2.fasta"
    rq_threshold = 0.999
    qname_rq = read_record_rq(bam_file=bam_file)
    filter_fasta(fa_filepath=fa_filepath, qname_rq_mapping=qname_rq, rq_threshold=rq_threshold)

if __name__ == "__main__":
    main_cli()
    pass