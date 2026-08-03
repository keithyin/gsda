"""Split a BAM into two files by channel list (keep / drop).

Reads a channel list file (one channel ID per line, extracted from qname
via the "/" delimiter), then splits each read of the input BAM into either
a "kept" BAM (channel in list) or a "dropped" BAM (channel not in list).
"""

import argparse
import os
import sys
from tqdm import tqdm

import pysam


def extract_channel(qname: str) -> int:
    """Extract channel ID from query name.

    qname format: <run_id>/<channel_id>
    e.g. '20260728_240601Y0004_Run0004/21900' → 21900

    Returns -1 on failure (never matches any valid channel).
    """
    try:
        return int(qname.rsplit("/", maxsplit=1)[-1])
    except (ValueError, IndexError):
        return -1


def read_channel_list(filename: str) -> set:
    """Read channel IDs from a text file, one per line."""
    channels = set()
    with open(filename, mode="r", encoding="utf8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                channels.add(int(line))
            except ValueError:
                sys.stderr.write(f"warning: skipping non-integer channel '{line}'\n")
    return channels


def split_bam_by_channel(bam_path: str, channels: set, kept_infix: str = ".kept", dropped_infix: str = ".dropped"):
    """Split input BAM into two output BAMs based on channel membership.

    Reads whose extracted channel is in *channels* are written to the kept BAM;
    all others go to the dropped BAM.

    Returns
    -------
    tuple[str, str]
        (kept_bam_path, dropped_bam_path)
    """
    bam_base = bam_path.rsplit(".", maxsplit=1)[0]
    kept_path = f"{bam_base}{kept_infix}.bam"
    dropped_path = f"{bam_base}{dropped_infix}.bam"

    tot = kept_count = dropped_count = 0

    with pysam.AlignmentFile(bam_path, mode="rb", threads=40, check_sq=False) as in_bam:
        with pysam.AlignmentFile(
            kept_path, mode="wb", threads=40, check_sq=False, header=in_bam.header
        ) as out_kept:
            with pysam.AlignmentFile(
                dropped_path, mode="wb", threads=40, check_sq=False, header=in_bam.header
            ) as out_dropped:
                for record in tqdm(in_bam.fetch(until_eof=True), desc=f"splitting {bam_path}"):
                    tot += 1
                    ch = extract_channel(record.query_name)
                    if ch in channels:
                        out_kept.write(record)
                        kept_count += 1
                    else:
                        out_dropped.write(record)
                        dropped_count += 1

    ratio = kept_count / tot if tot else 0.0
    sys.stderr.write(f"total={tot} kept={kept_count} dropped={dropped_count} ratio={ratio:.4f}\n")

    return kept_path, dropped_path


def main(args):
    channels = read_channel_list(args.channel_list)
    sys.stderr.write(f"loaded {len(channels)} channels from {args.channel_list}\n")
    kept, dropped = split_bam_by_channel(
        args.bam, channels, kept_infix=args.kept, dropped_infix=args.dropped
    )
    sys.stderr.write(f"output: kept={kept}  dropped={dropped}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--bam", type=str, required=True)
    parser.add_argument("--channel-list", type=str, required=True)
    parser.add_argument("--kept", type=str, default=".kept")
    parser.add_argument("--dropped", type=str, default=".dropped")
    main(parser.parse_args())
