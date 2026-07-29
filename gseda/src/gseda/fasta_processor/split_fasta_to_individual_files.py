#!/usr/bin/env python3
"""Split a multi-sequence FASTA file into individual sequence files.

Each sequence is written to its own .fasta file named after the
FASTA header (e.g. >STR1-1 -> STR1-1.fasta).

Example:
    python split_fasta_to_individual_files.py input.fasta --out-dir ./split_out
"""

import argparse
import os
import pathlib
import sys

import pysam


def main(input_fasta: str, out_dir: str) -> int:
    """Split a multi-sequence FASTA into individual per-sequence files.

    Parameters
    ----------
    input_fasta : str
        Path to the input multi-sequence FASTA file.
    out_dir : str
        Path to the output directory (created if it doesn't exist).

    Returns
    -------
    int
        0 on success, 1 on error.
    """
    seqs = {}
    with pysam.FastxFile(input_fasta) as f:
        for record in f:
            seqs[record.name] = record.sequence

    if not seqs:
        print("Error: no sequences found.", file=sys.stderr)
        return 1

    out_path = pathlib.Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    for name, seq in seqs.items():
        out_file = out_path / f"{name}.fasta"
        with open(out_file, "w", encoding="utf8") as fout:
            fout.write(f">{name}\n{seq}\n")

    print(f"Split {len(seqs)} sequences from {input_fasta}", file=sys.stderr)
    return 0


def main_cli():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="split_fasta_to_individual_files",
        description="Split a multi-sequence FASTA into individual sequence files.",
    )
    parser.add_argument("input", help="Input multi-sequence FASTA file")
    parser.add_argument("--out-dir", default="./split_output",
                        help="Output directory (default: ./split_output)")
    args = parser.parse_args()
    return main(args.input, args.out_dir)


if __name__ == "__main__":
    sys.exit(main_cli())
