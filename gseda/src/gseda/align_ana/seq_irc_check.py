import os
import mappy
import re
import pathlib
import sys
cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))


from mappy_ext import blast_like_alignment  # noqa: E402
from mappy_ext import calculate_identity_from_cigar  # noqa: E402


def calculate_identity_from_cigar(cigar_string):
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    match_count = 0
    mismatch_count = 0
    total_aligned = 0

    for length_str, operation in matches:
        length = int(length_str)
        if operation == '=':
            match_count += length
            total_aligned += length
        elif operation == 'X':
            mismatch_count += length
            total_aligned += length
        elif operation in ('I', 'D'):
            total_aligned += length

    if total_aligned == 0:
        return 0.0
    return match_count / total_aligned


def try_split_seq_1(seq: str):
    seq_len = len(seq)
    if seq_len < 200 or seq_len > 20000:
        return False

    # mid = seq_len // 2

    seq1 = seq
    seq2 = seq

    # print(seq1[334: 482])
    # print(seq2[2: 151])

    aligner = mappy.Aligner(seq=seq1, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)

    ok = False
    print("###################    try_split_seq_1    ##########################\n")
    for hit in aligner.map(seq2):
        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage1 = (hit.r_en - hit.r_st) / len(seq1)
        coverage2 = (hit.q_en - hit.q_st) / len(seq2)

        print("TSt:{}\tTEn:{}\tQSt:{}\tQEn:{}\tstrand:{}\tTLen:{}\tQLen:{}\tIden:{:.4f}\tQryCov:{:.4f}\tTagtCov:{:.4f}".format(
            hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.strand, len(
                seq1), len(seq2),
            calculate_identity_from_cigar(hit.cigar_str),
            (hit.q_en - hit.q_st) / len(seq2),
            (hit.r_en - hit.r_st) / len(seq1)))

        # blast_like_alignment(seq2, seq1, hit)

        if identity > 0.80 and coverage1 > 0.85 and coverage2 > 0.85:
            ok |= True
    print("\n###################    try_split_seq_1 DONE    ####################")

    return ok


def try_split_seq_2(seq: str):
    seq_len = len(seq)
    if seq_len < 200 or seq_len > 20000:
        return False

    mid = seq_len // 2

    seq1 = seq[:mid]
    seq2 = seq[mid:]

    # print(seq1[334: 482])
    # print(seq2[2: 151])

    aligner = mappy.Aligner(seq=seq1, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    print("###################    try_split_seq_2    ##########################\n")

    ok = False
    for hit in aligner.map(seq2):
        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage1 = (hit.r_en - hit.r_st) / len(seq1)
        coverage2 = (hit.q_en - hit.q_st) / len(seq2)

        print("TSt:{}\tTEn:{}\tQSt:{}\tQEn:{}\tstrand:{}\tTLen:{}\tQLen:{}\tIden:{:.4f}\tQryCov:{:.4f}\tTagtCov:{:.4f}".format(
            hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.strand, len(
                seq1), len(seq2),
            calculate_identity_from_cigar(hit.cigar_str),
            (hit.q_en - hit.q_st) / len(seq2),
            (hit.r_en - hit.r_st) / len(seq1)))

        blast_like_alignment(seq2, seq1, hit)

        if identity > 0.80 and coverage1 > 0.85 and coverage2 > 0.85:
            ok |= True
    print("\n###################    try_split_seq_2    ##########################")

    return ok


def main_cli():
    seq = "CCCCGGGGAAAATTTTGAGAGAGATGCTCAGCTCCGTTTCGGTTTCACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGGAGGCGCCGCTGCCAGGGGGCGTGCGGCAGCGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCGCTAGCAGGGCTGAAGAGGAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTGCTCAGCTCCGTTTCGGTTTCACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGGAGGCGCCGCTGCCAGGGGGCGTGCGGCAGCGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCGCTAGCAGGGCCGAAGAGAAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTCTCAGCTCCGTTTCGGTTTCACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGGAGGCGCCGCTGCCAGGGGGCGTGCGGCAGCGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCGCTAGCAGGGCTGAAGAGAAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTATCTCTCTCAAAATTTTCCCCGGGG"
    try_split_seq_1(seq=seq)

    print("\n\n\n")
    tag = try_split_seq_2(seq=seq)
    print(f"it is ok? = {tag}")


if __name__ == "__main__":
    main_cli()
