import pysam
from typing import List, Set
import argparse
import pysam.samtools
from tqdm import tqdm


def extract_non_multiple_mapping_reads(bam_file: str) -> Set[str]:
    qname_cnts = {}
    with pysam.AlignmentFile(filename=bam_file, mode="rb", threads=40) as bam_h:
        for record in tqdm(
            bam_h.fetch(), desc=f"extract_non_multiple_mapping_reads from {bam_file}"
        ):
            qname = record.query_name
            if qname not in qname_cnts:
                qname_cnts[qname] = 0
            qname_cnts[qname] += 1
    return set([qname for qname, cnt in qname_cnts.items() if cnt == 1])


def dump_non_multiple_mapping_reads(bam_file: str, qname_set):
    o_bam_file = "{}.non_multiple_mapping_reads.bam".format(
        bam_file.rsplit(".", maxsplit=1)[0]
    )

    with pysam.AlignmentFile(filename=bam_file, mode="rb", threads=40) as bam_h:
        with pysam.AlignmentFile(
            filename=o_bam_file, mode="wb", threads=40, header=bam_h.header
        ) as o_bam_h:
            for record in tqdm(
                bam_h.fetch(), desc=f"dump_non_multiple_mapping_reads to {o_bam_file}"
            ):
                if record.query_name in qname_set:
                    o_bam_h.write(record)

    pysam.samtools.index("-@", "40", "-b", o_bam_file)


def main(args):
    qname_set = extract_non_multiple_mapping_reads(args.bam)
    dump_non_multiple_mapping_reads(args.bam, qname_set)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("bam")
    main(args=parser.parse_args())
