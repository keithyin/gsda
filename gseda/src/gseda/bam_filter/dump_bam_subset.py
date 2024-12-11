import pysam
from tqdm import tqdm
from argparse import Namespace
import argparse


def read_wanted_channels(filename: str):

    channels = set()
    if (
        filename.endswith("fa")
        or filename.endswith("fasta")
        or filename.endswith("fna")
        or filename.endswith("fq")
        or filename.endswith("fastq")
    ):
        with pysam.FastxFile(filename=filename) as fh:
            for entry in fh:
                ch = entry.name.split("/")[1]
                try:
                    ch = int(ch)
                except:
                    continue
                channels.add(ch)
    elif filename.endswith("bam"):
        with pysam.AlignmentFile(
            filename=filename, mode="rb", threads=40, check_sq=False
        ) as bam_h:
            for record in bam_h.fetch(until_eof=True):
                channels.add(int(record.get_tag("ch")))

    return channels


def dump_bam_according_to_fastx_channel(subreads_bam: str, fastx_file: str, infix: str):
    channels = read_wanted_channels(fastx_file)
    o_bam_filename = "{}.{}.bam".format(subreads_bam.rsplit(".", maxsplit=1)[0], infix)
    with pysam.AlignmentFile(
        subreads_bam, mode="rb", threads=40, check_sq=False
    ) as in_bam:
        with pysam.AlignmentFile(
            o_bam_filename, mode="wb", threads=40, check_sq=False, header=in_bam.header
        ) as out_bam:

            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"dumping {o_bam_filename}"
            ):
                ch = int(record.get_tag("ch"))
                if ch in channels:
                    out_bam.write(record)


def main(args):

    dump_bam_according_to_fastx_channel(args.bam, args.fastx_file, args.infix)

    pass


if __name__ == "__main__":
    cli_args = {
        "subreads_bam": "/data/ccs_data/case-study/tr-at-error-kaipu/k2_5_subreads.bam",
        "fastx_file": "/data/ccs_data/case-study/tr-at-error-kaipu/sample5.fastq",
    }

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--bam", type=str, required=True, dest="bam")
    parser.add_argument("--fx-file-or-bam", type=str, required=True, dest="fastx_file")
    parser.add_argument("--infix", type=str, required=True)
    # main(Namespace(**cli_args))
    main(parser.parse_args())
    pass
