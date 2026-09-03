import pysam
import os
import sys
from tqdm import tqdm
import argparse

pre_dir = os.path.abspath(__file__).rsplit("/", maxsplit=2)[0]
print(pre_dir)
sys.path.append(pre_dir)

from fact_table_ana.utils import read_fastx_file


def compress_homopolymer(sequence):
    if not sequence:
        return "", []

    compressed = []
    counts = []

    current_char = sequence[0]
    count = 1

    for char in sequence[1:]:
        if char == current_char:
            count += 1
        else:
            compressed.append(current_char)
            counts.append(count)
            current_char = char
            count = 1

    # 添加最后一个字符和计数
    compressed.append(current_char)
    counts.append(count)

    return "".join(compressed), counts


def bam_hpc(bam_filepath: str):
    out_bam = "{}.hpc.bam".format(bam_filepath.rsplit(".", maxsplit=1)[0])

    with pysam.AlignmentFile(
        bam_filepath, mode="rb", threads=40, check_sq=False
    ) as in_bam_h:
        with pysam.AlignmentFile(
            out_bam, mode="wb", threads=40, check_sq=False, header=in_bam_h.header
        ) as out_bam_h:

            for record in tqdm(
                in_bam_h.fetch(until_eof=True), desc=f"extracting channel ref"
            ):
                record_new = pysam.AlignedSegment()
                record_new.query_name = record.query_name
                hpc_seq, counts = compress_homopolymer(record.query_sequence)
                record_new.query_sequence = hpc_seq
                record_new.set_tags(
                    [
                        ("dw", counts, "I"),
                        ("ch", int(record.get_tag("ch")), "I"),
                    ]
                )
                out_bam_h.write(record_new)


def fa_hpc(filepath: str):
    out_hpc_file = "{}.hpc.fa".format(filepath.rsplit(".", maxsplit=1)[0])
    out_hpc_cnt_file = "{}.hpc.txt".format(filepath.rsplit(".", maxsplit=1)[0])

    out_hpc_file_h = open(out_hpc_file, mode="w", encoding="utf8")
    out_hpc_cnt_file_h = open(out_hpc_cnt_file, mode="w", encoding="utf8")

    fastx_data = read_fastx_file(filepath)

    for qname, seq in fastx_data.items():
        hpc_seq, cnt = compress_homopolymer(seq[0])
        out_hpc_file_h.write(">{}\n{}\n".format(qname, hpc_seq))
        out_hpc_cnt_file_h.write(">{}\n{}\n".format(qname, ",".join(map(str, cnt))))


def main(args):
    for filepath in args.files:
        if filepath.endswith(".bam"):
            bam_hpc(filepath)
            continue
        if (
            filepath.endswith("fa")
            or filepath.endswith("fasta")
            or filepath.endswith("fna")
        ):
            fa_hpc(filepath)
            continue

        raise ValueError("invalid data format")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("files", nargs="+")

    main(args=parser.parse_args())
