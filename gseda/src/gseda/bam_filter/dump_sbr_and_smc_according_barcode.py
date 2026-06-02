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


def dump_bam_according_to_fastx_channel(in_bam: str, fastx_file: str, infix: str):
    channels = read_wanted_channels(fastx_file)
    o_bam_filename = "{}.{}.bam".format(
        in_bam.rsplit(".", maxsplit=1)[0], infix)
    with pysam.AlignmentFile(
        in_bam, mode="rb", threads=40, check_sq=False
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


def main_cli():

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--bam", type=str, required=True)
    parser.add_argument("--demuxed-file", type=str,
                        required=True, help=".fastq")
    parser.add_argument("--barcode-name", type=str, required=True)
    args = parser.parse_args()
    dump_bam_according_to_fastx_channel(
        args.bam, args.demuxed_file, args.barcode_name)

    pass


if __name__ == "__main__":
    cli_args = {
        "smc_bam": "/data1/ccs_data/202605-baiying-500-1M-diff/20260513_250802Y0012_Run0001/smc8.smc_all_reads.bam",
        "adapter_bam": "/data1/ccs_data/202605-baiying-500-1M-diff/20260513_250802Y0012_Run0001/demuxed_q8.bam",
        "demuxed_file": "/data1/ccs_data/202605-baiying-500-1M-diff/Y0012-1M/Adaptor-barcode202-0.fastq",
        "infix": "202"
    }
    cli_args = Namespace(**cli_args)
    dump_bam_according_to_fastx_channel(cli_args.smc_bam, cli_args.demuxed_file, cli_args.infix)
    dump_bam_according_to_fastx_channel(cli_args.adapter_bam, cli_args.demuxed_file, cli_args.infix)
    pass
