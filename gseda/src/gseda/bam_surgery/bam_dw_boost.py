#!/usr/bin/env python3
"""Add a constant to every element of each read's dw tag and write a new BAM.

The 'dw' tag is a per-base auxiliary array (homopolymer run-length / dwell
counts). This script reads a BAM, adds a constant value to each element of the
'dw' tag for every read, and writes the result to a new BAM.

All other fields (sequence, qualities, flags, names, coordinates, other tags)
are preserved.
"""

import argparse
import array
import sys

import pysam
from tqdm import tqdm


def boost_dw(
    input_bam: str,
    output_bam: str,
    boost: int,
    threads: int = 40,
    show_progress: bool = True,
):
    """Read input BAM, add `boost` to every dw element, write to output BAM."""
    n_modified = 0

    with pysam.AlignmentFile(
        input_bam, mode="rb", check_sq=False, threads=threads
    ) as in_f, pysam.AlignmentFile(
        output_bam, mode="wb", check_sq=False, threads=threads, header=in_f.header
    ) as out_f:

        iterator = in_f.fetch(until_eof=True)
        if show_progress:
            iterator = tqdm(iterator, desc=f"boosting dw {input_bam}")

        for record in iterator:
            if record.has_tag("dw"):
                dw = record.get_tag("dw")
                if isinstance(dw, array.array):
                    # Preserve the original typecode (B/S/I) so the stored BAM type is unchanged
                    record.set_tag(
                        "dw", array.array(dw.typecode, (v + boost for v in dw))
                    )
                elif isinstance(dw, (list, tuple)):
                    record.set_tag("dw", [v + boost for v in dw])
                else:  # scalar fallback
                    record.set_tag("dw", dw + boost)
                n_modified += 1

            out_f.write(record)

    print(f"Boosted dw tag of {n_modified} reads")


def main_cli():
    """CLI entry point: parse args and call boost_dw."""
    parser = argparse.ArgumentParser(
        description="Add a constant value to every element of each read's dw tag "
                    "and write the result to a new BAM."
    )
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file path")
    parser.add_argument(
        "--boost", type=int, required=True,
        help="Constant value to add to every element of the dw tag"
    )
    parser.add_argument(
        "--threads", type=int, default=40,
        help="Number of threads for BAM I/O (default: 40)"
    )
    parser.add_argument("--no-progress", action="store_true", help="Disable progress bar")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    boost_dw(
        input_bam=args.input_bam,
        output_bam=args.output_bam,
        boost=args.boost,
        threads=args.threads,
        show_progress=not args.no_progress,
    )


if __name__ == "__main__":
    main_cli()
