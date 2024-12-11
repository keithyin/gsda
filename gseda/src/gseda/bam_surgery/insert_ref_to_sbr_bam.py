import pysam
import os
import sys
from tqdm import tqdm
import argparse

pre_dir = os.path.abspath(__file__).rsplit("/", maxsplit=3)[0]
print(pre_dir)
sys.path.append(pre_dir)

import gseda.fact_table_ana.utils as gseda_utils


def insert_ref_to_sbr_bam(sbr_bam: str, smc2ref_bam: str, ref_file: str):
    out_bam = "{}.with-ref.bam".format(sbr_bam.rsplit(".", maxsplit=1)[0])
    ref_data = gseda_utils.read_fastx_file(ref_file)

    channel_refs = {}
    with pysam.AlignmentFile(smc2ref_bam, mode="rb", threads=40) as smc2ref_bam_h:
        for record in tqdm(smc2ref_bam_h.fetch(), desc=f"extracting channel ref"):
            ch = int(record.get_tag("ch"))
            ref_start = record.reference_start
            ref_end = record.reference_end
            ref_name = record.reference_name
            ref_len = len(ref_data[ref_name][0])

            ref_start = max(0, ref_start - 10)
            ref_end = min(ref_len, ref_end + 11)

            ref_sub_seq = ref_data[ref_name][0][ref_start:ref_end]
            assert ch not in channel_refs
            channel_refs[ch] = ref_sub_seq
    with pysam.AlignmentFile(sbr_bam, mode="rb", threads=40, check_sq=False) as in_h:
        with pysam.AlignmentFile(
            out_bam, mode="wb", threads=40, check_sq=False, header=in_h.header
        ) as out_h:
            for record in tqdm(in_h.fetch(until_eof=True), desc=f"dumping1 {out_bam}"):
                out_h.write(record)

            for ch, ref_seq in tqdm(channel_refs.items(), desc=f"dumping2 {out_bam}"):
                record = pysam.AlignedSegment()
                record.query_name = f"09_REF_{ch}"
                record.query_sequence = ref_seq
                record.set_tag("ch", ch, value_type="I")
                out_h.write(record)


def main(args):
    insert_ref_to_sbr_bam(args.sbr_bam, args.smc2ref_bam, args.ref_file)
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--sbr-bam", dest="sbr_bam")
    parser.add_argument("--smc2ref-bam", dest="smc2ref_bam")
    parser.add_argument("--ref-file", dest="ref_file")

    main(args=parser.parse_args())
