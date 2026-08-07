#!/usr/bin/env python3
"""Trim polyN (consecutive identical character) runs from BAM read ends.

Detects homopolymer-like runs at the 5' start and/or 3' end of each read's
query_sequence, removes them from both query_sequence and query_qualities,
updates CIGAR soft-clip counts correspondingly, and trims per-base auxiliary
array tags (ar, dw, cr, etc.) by the same slice range.

All other fields (flags, names, reference coordinates, mapping quality) are preserved.
"""

import argparse
import array
import sys

import pysam
from tqdm import tqdm


def detect_polyN_at(seq: str, at_start: bool):
    """Detect longest run of identical characters at one end of a sequence.

    Parameters
    ----------
    seq : str or None
        The query sequence to inspect.
    at_start : bool
        If True, inspect the beginning; if False, inspect the end (via reversal).

    Returns
    -------
    dict or None
        ``{"char": <first char>, "length": <run length>}`` if a run is found, else None.
    """
    if not seq:
        return None

    base = seq[0] if at_start else seq[-1]
    run_len = 0
    for ch in (seq if at_start else reversed(seq)):
        if ch == base:
            run_len += 1
        else:
            break

    return {"char": base, "length": run_len} if run_len > 0 else None


def detect_polyN_both_ends(seq: str):
    """Detect polyN runs at both ends independently.

    Returns (poly5_result_dict or None, poly3_result_dict or None).
    The two detections are performed on the original sequence to avoid interference.
    """
    poly5 = detect_polyN_at(seq, at_start=True)
    # For 3' end, reverse and reuse the same detection logic
    poly3_reversed = detect_polyN_at(seq[::-1], at_start=True)

    return poly5, poly3_reversed


def modify_cigartuples_for_polyN(cigartuples, trim_start_len: int, trim_end_len: int):
    """Adjust soft-clip CIGAR operations to account for polyN removal.

    Only modifies S (soft-clip) operations at the very start or end of the CIGAR tuple.
    Non-S operations (M/I/D/N/X) are never modified for safety.

    Parameters
    ----------
    cigartuples : tuple of (op_code, length) from pysam.AlignedSegment.cigartuples
    trim_start_len : int
        Number of bases removed from the 5' end.
    trim_end_len : int
        Number of bases removed from the 3' end.

    Returns
    -------
    tuple of (op_code, length) or None
        New CIGAR tuple if changed, else None (caller can skip assignment).
    """
    if not cigartuples or trim_start_len == 0 and trim_end_len == 0:
        return None

    cigar_list = list(cigartuples)

    # Trim from 5' end of the CIGAR tuple
    if trim_start_len > 0:
        op, length = cigar_list[0]
        if op == 4:  # S (soft-clip) in SAM spec
            cigar_list[0] = (op, max(0, length - trim_start_len))

    # Trim from 3' end of the CIGAR tuple
    if trim_end_len > 0:
        op, length = cigar_list[-1]
        if op == 4:  # S (soft-clip) in SAM spec
            cigar_list[-1] = (op, max(0, length - trim_end_len))

    # Return None if no changes were made
    if len(cigar_list) != len(cigartuples) or not all(a == b for a, b in zip(cigar_list, cigartuples)):
        return tuple(cigar_list)
    return None


def identify_perbase_tags(tags, query_length: int):
    """Return dict of per-base aux tags whose length matches query_length.

    Scans pysam.AlignedSegment.tags and selects fields that are array.array
    (or list/tuple) with length exactly equal to the original query length.
    These correspond to per-base auxiliary data like 'ar' (base quality),
    'dw' (homopolymer run-lengths), 'cr' (channel offset), etc.
    """
    perbase = {}
    for name, value in tags:
        if isinstance(value, str):
            continue  # string tags (e.g., RG) are not per-base
        if isinstance(value, (int, float)):
            continue  # scalar tags
        if isinstance(value, array.array) and len(value) == query_length:
            perbase[name] = value
        elif isinstance(value, (list, tuple)) and len(value) == query_length:
            perbase[name] = list(value)  # normalize to list for slicing
    return perbase


def trim_polyN(
    input_bam: str,
    output_bam: str,
    min_poly_length: int = 5,
    trim_start: bool = True,
    trim_end: bool = True,
    threads: int = 40,
):
    """Core function: read BAM, trim polyN from each read, write new BAM.

    Parameters
    ----------
    input_bam : str
        Path to the input BAM file.
    output_bam : str
        Path to the output BAM file.
    min_poly_length : int
        Minimum consecutive identical characters to consider as polyN (default: 5).
    trim_start : bool
        If True, trim polyN at the 5' end of each read.
    trim_end : bool
        If True, trim polyN at the 3' end of each read.
    threads : int
        Number of threads for BAM I/O (default: 40).
    """
    n_modified = 0

    with pysam.AlignmentFile(
        input_bam, mode="rb", check_sq=False, threads=threads
    ) as in_f, pysam.AlignmentFile(
        output_bam, mode="wb", threads=threads, check_sq=False, header=in_f.header
    ) as out_f:

        for record in tqdm(in_f.fetch(until_eof=True), desc=f"trimming {input_bam}"):
            seq = record.query_sequence
            if seq is None or len(seq) == 0:
                out_f.write(record)
                continue

            orig_len = len(seq)

            # --- Detect polyN at both ends ---
            poly5, poly3_reversed = detect_polyN_both_ends(seq)
            trim_start_len = (
                poly5["length"] if poly5 and poly5["length"] >= min_poly_length else 0
            )
            trim_end_len = (
                poly3_reversed["length"]
                if poly3_reversed and poly3_reversed["length"] >= min_poly_length
                else 0
            )

            # Nothing to trim? pass through unchanged
            if trim_start_len == 0 and trim_end_len == 0:
                out_f.write(record)
                continue

            # Effective flags based on user --trim-ends option
            do_trim_start = trim_start and trim_start_len > 0
            do_trim_end = trim_end and trim_end_len > 0

            total_trim = (trim_start_len if do_trim_start else 0) + (
                trim_end_len if do_trim_end else 0
            )
            new_len = orig_len - total_trim

            # Save original qualities before CIGAR modification
            # (pysam may clear query_qualities when cigartuples is reassigned)
            saved_qual = record.query_qualities

            # --- Update CIGAR (only for mapped reads with soft-clip at ends) ---
            if not record.is_unmapped and record.cigartuples:
                new_cigar = modify_cigartuples_for_polyN(
                    record.cigartuples,
                    trim_start_len if do_trim_start else 0,
                    trim_end_len if do_trim_end else 0,
                )
                if new_cigar is not None:
                    record.cigartuples = new_cigar

            # --- Update query_sequence ---
            seq_slice_start = trim_start_len if do_trim_start else 0
            seq_slice_end = orig_len - (trim_end_len if do_trim_end else 0)
            record.query_sequence = seq[seq_slice_start:seq_slice_end]

            # --- Update query_qualities ---
            # Use saved_qual because pysam may clear query_qualities after CIGAR modification
            if saved_qual is not None:
                qual_slice = saved_qual[seq_slice_start : seq_slice_end]
                record.query_qualities = array.array("H", qual_slice)
            else:
                # Fallback: create uniform quality array
                record.query_qualities = array.array("H", [50] * new_len)

            # --- Update per-base auxiliary tags ---
            perbase_tags = identify_perbase_tags(record.tags, orig_len)
            for tag_name, tag_value in perbase_tags.items():
                sliced = tag_value[seq_slice_start:seq_slice_end]
                record.set_tag(tag_name, sliced)

            # Write the modified read
            out_f.write(record)
            n_modified += 1

    print(f"Trimmed {n_modified} reads")


def parse_trim_ends(ends_str: str):
    """Convert CLI trim-ends string to (trim_start, trim_end) booleans."""
    return {
        "both": (True, True),
        "5p":   (True, False),
        "3p":   (False, True),
    }[ends_str]


def main_cli():
    """CLI entry point: parse args and call trim_polyN."""
    parser = argparse.ArgumentParser(
        description="Trim polyN (consecutive identical character) runs from BAM read ends "
                    "and synchronously update CIGAR, qualities, and per-base auxiliary tags."
    )
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file path")
    parser.add_argument(
        "--min-poly-length", type=int, default=5,
        help="Minimum consecutive identical characters to consider as polyN (default: 5)"
    )
    parser.add_argument(
        "--trim-ends", default="both", choices=["5p", "3p", "both"],
        help="Which ends to trim: '5p'=start, '3p'=end, 'both'=both (default)"
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
    trim_start, trim_end = parse_trim_ends(args.trim_ends)

    trim_polyN(
        input_bam=args.input_bam,
        output_bam=args.output_bam,
        min_poly_length=args.min_poly_length,
        trim_start=trim_start,
        trim_end=trim_end,
        threads=args.threads,
    )


if __name__ == "__main__":
    main_cli()
