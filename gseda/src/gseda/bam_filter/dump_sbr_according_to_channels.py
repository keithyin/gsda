import pysam
from tqdm import tqdm
import argparse


def dump_bam_according_channels(subreads_bam: str, channels, infix: str):
    o_bam_filename = "{}.{}.bam".format(
        subreads_bam.rsplit(".", maxsplit=1)[0], infix)
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


def main_cli():

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--bam", type=str,
                        required=True)
    args = parser.parse_args()
    dump_bam_according_channels(
        args.demuxed_file, args.barcode_name)

    pass


def main():
    bam = "/data1/ccs_data/20260702_250302Y0006_Run0001/20260702_250302Y0006_Run0001_adapter.bam"
    channels = set([1034445])
    infix = "1034445"
    dump_bam_according_channels(bam, channels, infix)


if __name__ == "__main__":
    # cli_args = {
    #     "bam": "/data1/ccs_data/202603-henan-nongda/20260325_240601Y0009_Run0001/20260325_240601Y0009_Run0001_smc_8.bam",
    #     "infix":
    # }
    # cli_args = Namespace(**cli_args)
    # dump_bam_according_channels(
    #     cli_args.bam, cli_args.demuxed_file, cli_args.barcode_name)
    # main(Namespace(**cli_args))
    # main_cli(Namespace(**cli_args))
    main()
    pass
