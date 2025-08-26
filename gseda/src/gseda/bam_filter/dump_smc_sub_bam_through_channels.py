import pysam
from tqdm import tqdm
import argparse


def read_channels(fname: str):
    channels = set()
    with open(fname, mode="r", encoding="utf8") as in_data:
        for idx, line in enumerate(in_data):
            line = line.strip()
            try:
                channel = int(line)
                channels.add(channel)
            except Exception as e:
                print(f"line: {idx}. error {e}")
    return channels


def dump_smc_bam_subset(inp_bam_path: str, out_bam_path: str, channels):
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
                ch = int(record.get_tag("ch"))
                if ch in channels:
                    out_bam.write(record)
                else:
                    filtered += 1

    print(f"tot:{tot}, filtered:{filtered}, ratio:{filtered/tot}")


def main(args):
    channels = read_channels(args.channel_filename)
    dump_smc_bam_subset(args.in_bam, args.out_bam, channels)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("in_bam", type=str)
    parser.add_argument("out_bam", type=str)
    parser.add_argument("channel_filename", type=str)

    main(parser.parse_args())
