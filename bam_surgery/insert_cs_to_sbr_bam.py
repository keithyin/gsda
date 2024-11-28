import pysam
import os
import sys
from tqdm import tqdm
import argparse

pre_dir = os.path.abspath(__file__).rsplit("/", maxsplit=2)[0]
sys.path.append(pre_dir)


def insert_cs_to_sbr_bam(sbr_bam: str, smc_bam: str):
    out_bam = "{}.with-cs.bam".format(sbr_bam.rsplit(".", maxsplit=1)[0])

    with pysam.AlignmentFile(sbr_bam, mode="rb", threads=40, check_sq=False) as in_h:
        with pysam.AlignmentFile(
            out_bam, mode="wb", threads=40, check_sq=False, header=in_h.header
        ) as out_h:
            for record in tqdm(in_h.fetch(until_eof=True), desc=f"dumping1 {out_bam}"):
                out_h.write(record)

            with pysam.AlignmentFile(
                smc_bam, mode="rb", threads=40, check_sq=False
            ) as in_h:
                for in_record in tqdm(
                    in_h.fetch(until_eof=True), desc=f"dumping2 {out_bam}"
                ):
                    qname = in_record.query_name
                    record = pysam.AlignedSegment()
                    record.query_name = f"00_{qname}"
                    record.query_sequence = in_record.query_sequence
                    record.set_tag("ch", int(in_record.get_tag("ch")), value_type="I")
                    out_h.write(record)


def main(args):
    insert_cs_to_sbr_bam(args.sbr_bam, args.smc_bam)
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--sbr-bam", dest="sbr_bam")
    parser.add_argument("--smc-bam", dest="smc_bam")

    main(args=parser.parse_args())
