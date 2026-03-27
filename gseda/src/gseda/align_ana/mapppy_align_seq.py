import mappy
import pysam
import re
import pathlib
import os
import sys
cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))


from mappy_ext import blast_like_alignment  # noqa: E402
from mappy_ext import calculate_identity_from_cigar  # noqa: E40


def calculate_identity_from_cigar(cigar_string):
    """
    根据CIGAR字符串计算identity

    Parameters:
    cigar_string (str): CIGAR字符串, 如 "370=1I89=1D26=1I208=1D1250=1D578=1D2=2X1=1I578=1I1248=1I210=1D23=1I90=1D373="

    Returns:
    float: identity值
    """
    # 使用正则表达式解析CIGAR字符串
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    # 初始化计数器
    match_count = 0      # 匹配的碱基数（=）
    mismatch_count = 0   # 错配的碱基数（X）
    total_aligned = 0    # 总比对长度（= + X + I + D）

    for length_str, operation in matches:
        length = int(length_str)

        if operation == '=':
            # 精确匹配
            match_count += length
            total_aligned += length
        elif operation == 'X':
            # 错配
            mismatch_count += length
            total_aligned += length
        elif operation == 'I':
            # 插入：参考序列没有，查询序列有
            total_aligned += length
        elif operation == 'D':
            # 删除：参考序列有，查询序列没有
            total_aligned += length

    # 计算identity
    if total_aligned == 0:
        return 0.0

    identity = match_count / total_aligned
    return identity


def read_fasta(fpath: str):
    with open(fpath, mode="r", encoding="utf8") as f:
        line = f.readline().strip()
        if line.startswith(">"):
            header = line
        seq = ""
        while True:
            new_line = f.readline().strip()
            if new_line == "":
                return header, seq
            seq += new_line


def read_fasta2(fpath):
    with pysam.FastxFile(fpath) as f:
        for record in f:
            return record.name, record.sequence


def main():
    target_seq = "GACTGGCCGGTCTTAGGATGTGTAGGTGAGTCGAACAATCTGACACCCAGGTCCAAGATTTCCCTATTAGCGGTGTCGATCAGTGCGAAGCCGCAAGATGCGATGCCTGGGTCCATACCCAGCACCAGATGTCTAGGCAAACCCAGCTTACCGATTTCCCTGTTGATTGTGATCTCGGCTGCTGGGACTCCGTGGATACCGACCTTCCGCTTCTTCTTTGGGGCCATCTTATCGTCATCGTCTTTGTAATCAATATCATGATCCTTGTAGTCTCCGTCGTGGTCCTTATAGTCCATGGTGGCACCGGTCCAACCTGAAAAAAAGTGATTTCAGGCAGGTGCTCCAGGTAATTAAACATTAATACCCCACCAACCAACCATCCCTTAAACCCTTACCTCTTGCTCAGCTAATTACAGCCCGGAGGAGAAGGGCCGTCCCGCCCGCTCACCTGTGGGAGTAACGCGGTCAGTCAGAGCCGGGGCGGGCGGCGCGAGGCGGCGGCGGAGCGGGGCACGGGGCGAAGGCAGCGCGCAGCGACTCCCGCCCGCCGCGCGCTTCGCTTTTTATAGGGCCGCCGCCGCCGCCGCCTCGCCATAAAAGGAAACTTTCGGAGCGCGCCGCTCTGATTGGCTGCCGCCGCACCTCTCCGCCTCGCCCCGCCCCGCCCCTCGCCCCGCCCCGCCCCGCCTGGCGCGCGCCCCCCCCCCCCCCCCGCCCCCATCGCTGCACAAAATAATTAAAAAATAAATAAATACAAAATTGGGGGTGGGGAGGGGGGGGAGATGGGGAGAGTGAAGCAGAACGTGGGGCTCACCTCGACCATGGTAATAGCGATGACTAATACGTAGATGTACTGCCAAGTAGGAAAGTCCCATAAGGTCATGTACTGGGCACAATGCCAGGCGGGCCATTTACCGTCATTGACGTCAATAGGGGGCGTACTTGGCATATGATACACTTGATGTACTGCCAAGTGGGCAGTTTACCG"
    query_seq = "CGGTAAACTGCCCACTTGGCAGTACATCAAGTGTACATATGCCAAGTACGCCCCTATTTGACGTCCAATGACGGTAAATGGCCCGCCCTGGCATTGTGCCCATACATGACCTTATGGGACTTTCCTACTTGGCAGACATCTACGTATTTAGTCATCGCTATTACCCATGTCGAGGTGAGCCCCACGTTCTGCTTCACTGTCCCCATCTCCACCCCCCTCCCCACCCCCATTTTGTATTTCTTTATTTTTTTAAATTATTTTGTGCAGCCGGAGTGGGGGCGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGAAGAAAAAGAAAGGTGTTGTATA"
    aligner = mappy.Aligner(seq=target_seq, extra_flags=67108864,
                            k=11, w=1, best_n=10, n_threads=1)
    cnt = 0
    for hit in aligner.map(query_seq):
        cnt += 1
        print("ctg:{}\tr_st:{}\tr_end:{}\tq_st:{}\tq_end:{}\tstrand:{}\tref_len:{}\tquery_len:{}\tideneity:{:.4f}\tQueryCoverage:{:.4f}\tTargetCoverage:{:.4f}".format(
            hit.ctg, hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.strand, len(
                target_seq), len(query_seq),
            calculate_identity_from_cigar(hit.cigar_str),
            (hit.q_en - hit.q_st) / len(query_seq),
            (hit.r_en - hit.r_st) / len(target_seq)))
        blast_like_alignment(query_seq, target_seq, hit)

    pass


if __name__ == "__main__":
    main()
