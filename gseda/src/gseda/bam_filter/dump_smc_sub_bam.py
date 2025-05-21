import pysam
from tqdm import tqdm
from argparse import Namespace
import argparse


def dump_smc_bam_subset(inp_bam_path: str, out_bam_path: str, rq_thr: float):
    filtered = 0
    tot = 0
    with pysam.AlignmentFile(
        inp_bam_path, mode="rb", threads=40, check_sq=False
    ) as in_bam:
        with pysam.AlignmentFile(
            out_bam_path, mode="wb", threads=40, check_sq=False, header=in_bam.header
        ) as out_bam:
            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"processing {inp_bam_path}"
            ):
                tot += 1
                record_rq = record.get_tag("rq")
                if record_rq < rq_thr:
                    filtered += 1
                    continue
                out_bam.write(record)
    print(f"tot:{tot}, filtered:{filtered}, ratio:{filtered/tot}")


def main(args):

    dump_smc_bam_subset(args.in_bam, args.out_bam, args.rq_low)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("in_bam", type=str)
    parser.add_argument("out_bam", type=str)
    parser.add_argument("--rq-low", type=float, dest="rq_low")

    # main(Namespace(**cli_args))
    main(parser.parse_args())
