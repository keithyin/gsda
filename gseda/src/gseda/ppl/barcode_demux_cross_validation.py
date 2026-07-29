#!/usr/bin/env python3
"""Cross-validate alignment results against barcode demultiplexing.

Two independent assignment systems for the same set of reads:
  - Alignment CSV maps qname -> rname (reference name) or empty (unaligned)
  - Barcode demux groups reads by barcode in separate fastq/bam files

Assuming alignment is correct, report the barcode demux error rate.
Assuming barcode demux is correct, report the alignment error rate.

Usage:
    python barcode_demux_cross_validation.py \\
        --alignment-csv /path/to/aligned_metric_fact.csv \\
        --barcode-dir /path/to/barcode_assign/
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


# Reference-barcode mapping: expected rname -> barcode name
BARCODE_MAP = {
    "STR1-1": "Barcode01",
    "STR1-2": "Barcode02",
    "STR1-5": "Barcode03",
    "STR1-8": "Barcode04",
    "STR2-1": "Barcode05",
    "STR2-2": "Barcode06",
    "STR2-3": "Barcode07",
    "STR2-4": "Barcode08",
    "STR2-9": "Barcode09",
    "STR3-3": "Barcode10",
}

# Reverse mapping: barcode name -> expected rname
BARCODE_TO_REF = {v: k for k, v in BARCODE_MAP.items()}

HEADER = "\n" + "=" * 72


def extract_fastq_qnames(fastq_path: Path, n_sample: int) -> list[str]:
    """Extract read names from a FASTQ file.

    Each read occupies 4 lines; the name is on line 1 (after @).
    If n_sample > 0, only load that many reads.
    """
    qnames = []
    count = 0
    with open(fastq_path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            f.readline()  # sequence line (skip)
            f.readline()  # '+' line (skip)
            f.readline()  # quality line (skip)
            qname = header[1:].split()[0]  # strip '@' and take first token
            qnames.append(qname)
            count += 1
            if 0 < n_sample <= count:
                break
    return qnames


def load_alignment_csv(csv_path: Path) -> dict[str, str]:
    """Load alignment CSV into {qname: rname} dict.

    Tab-separated with columns: qname, rname, ...
    Empty rname means unaligned.
    """
    mapping = {}
    # Read in two passes for efficiency on large files:
    # First pass: count total lines to report progress
    total_lines = 0
    with open(csv_path, "r") as f:
        _ = f.readline()  # skip header
        for line in f:
            total_lines += 1

    print(f"Loading alignment CSV ({total_lines:,} reads) ...", file=sys.stderr)
    with open(csv_path, "r") as f:
        _ = f.readline()  # skip header
        lineno = 0
        for line in f:
            lineno += 1
            if lineno % 50_000 == 0:
                print(f"  ... {lineno:,}/{total_lines:,}", file=sys.stderr)
            parts = line.rstrip("\n").split("\t")
            qname = parts[0]
            rname = parts[1]
            mapping[qname] = rname

    return mapping


def collect_barcode_qnames(
    barcode_dir: Path, n_sample: int
) -> dict[str, list[str]]:
    """Collect qnames from all barcode fastq files.

    Returns {barcode_name: [qname1, qname2, ...]}.
    Prefers .fastq; falls back to counting bam files if no fastq found.
    """
    result = {}
    for i in range(1, 11):
        barcode = f"Barcode{i:02d}"
        fastq_path = barcode_dir / f"{barcode}.fastq"
        if not fastq_path.exists():
            print(f"WARNING: {fastq_path} not found, skipping", file=sys.stderr)
            continue
        qnames = extract_fastq_qnames(fastq_path, n_sample)
        result[barcode] = qnames
    return result


def cross_validate_alignment_correct(
    barcode_qnames: dict[str, list[str]],
    alignment_map: dict[str, str],
) -> None:
    """Assuming alignment is correct, compute barcode demux error rate."""
    print(HEADER)
    print("  SCENARIO 1: Assume ALIGNMENT is correct")
    print("  --> What fraction of reads in each barcode are on the WRONG reference?")
    print(HEADER)

    # Global stats
    total_barcode_reads = 0
    total_wrong_ref = 0
    total_unaligned_in_barcode = 0

    per_barcode = {}

    for barcode, qnames in sorted(barcode_qnames.items()):
        expected_ref = BARCODE_TO_REF.get(barcode, "UNKNOWN")
        n_total = len(qnames)
        n_aligned_to_expected = 0
        n_wrong_ref = 0
        n_unaligned = 0
        wrong_ref_counts = defaultdict(int)

        for qn in qnames:
            rname = alignment_map.get(qn, "")
            if not rname:
                n_unaligned += 1
                continue
            if rname == expected_ref:
                n_aligned_to_expected += 1
            else:
                n_wrong_ref += 1
                wrong_ref_counts[rname] += 1

        total_barcode_reads += n_total
        total_wrong_ref += n_wrong_ref
        total_unaligned_in_barcode += n_unaligned

        per_barcode[barcode] = {
            "n_total": n_total,
            "n_expected": n_aligned_to_expected,
            "n_wrong_ref": n_wrong_ref,
            "n_unaligned": n_unaligned,
            "wrong_refs": dict(wrong_ref_counts),
        }

        if n_total > 0:
            wrong_pct = 100.0 * n_wrong_ref / n_total
            unalign_pct = 100.0 * n_unaligned / n_total
        else:
            wrong_pct = 0
            unalign_pct = 0

        print(f"\n{barcode} (expected: {expected_ref}):")
        print(
            f"  Total reads in barcode file: {n_total:>8,}"
        )
        print(
            f"  Aligned to expected ref ({expected_ref}): "
            f"{n_aligned_to_expected:>8,} ({100*n_aligned_to_expected/n_total:.2f}%)"
            if n_total > 0 else ""
        )
        if n_wrong_ref > 0:
            print(
                f"  Aligned to WRONG ref:          "
                f"{n_wrong_ref:>8,} ({wrong_pct:.2f}%)",
            )
            for ref, cnt in sorted(wrong_ref_counts.items(), key=lambda x: -x[1]):
                print(f"      -> {ref}: {cnt:,}")
        if n_unaligned > 0:
            print(
                f"  Aligned to nothing (empty rname): "
                f"{n_unaligned:>8,} ({unalign_pct:.2f}%)"
            )

    if total_barcode_reads > 0:
        overall_wrong = 100.0 * total_wrong_ref / total_barcode_reads
        overall_unalign = 100.0 * total_unaligned_in_barcode / total_barcode_reads
    else:
        overall_wrong = overall_unalign = 0

    print(
        f"\n{'─' * 60}",
        file=sys.stderr,
    )
    print(f"\nOVERALL (all barcodes combined):")
    print(f"  Total reads across all barcode files: {total_barcode_reads:,}")
    print(
        f"  Total reads aligned to wrong reference: {total_wrong_ref:,}"
        f" ({overall_wrong:.2f}%)"
    )
    print(
        f"  Total reads with no alignment (empty rname): {total_unaligned_in_barcode:,}"
        f" ({overall_unalign:.2f}%)"
    )
    print(
        f"\n  BARCODE DEMUX ERROR RATE (assuming alignment correct): "
        f"{overall_wrong:.2f}%"
    )

    # Print sample qnames for manual verification
    print("\n\nSample qnames from cross-contaminated reads:")
    for barcode in sorted(barcode_qnames):
        qnames = barcode_qnames[barcode]
        expected_ref = BARCODE_TO_REF.get(barcode, "UNKNOWN")
        wrong_refs_sample = {}
        for qn in qnames:
            rname = alignment_map.get(qn, "")
            if not rname or rname == expected_ref:
                continue
            if len(wrong_refs_sample) < 5:
                wrong_refs_sample.setdefault(rname, []).append(qn)
        for ref, samples in sorted(wrong_refs_sample.items()):
            for sn in samples[:3]:
                print(f"    {barcode} expected={expected_ref}, aligned={ref}: {sn}")


def cross_validate_barcode_correct(
    barcode_qnames: dict[str, list[str]],
    alignment_map: dict[str, str],
) -> None:
    """Assuming barcode demux is correct, compute alignment error rate."""
    print(HEADER)
    print("  SCENARIO 2: Assume BARCODE DEMUX is correct")
    print("  --> What fraction of reads aligned to each reference ended up in a WRONG barcode?")
    print(HEADER)

    # Collect per-reference totals from barcode files
    ref_total = defaultdict(int)

    for barcode, qnames in sorted(barcode_qnames.items()):
        expected_ref = BARCODE_TO_REF.get(barcode, "UNKNOWN")
        for qn in qnames:
            rname = alignment_map.get(qn, "")
            if rname:
                ref_total[rname] += 1

    # For each reference with reads aligned, where do those reads go?
    print(f"\n{'REF':<12} {'Total aligned':>14} "
          f"{'Expected barcode':>18} {'Correct count':>14} "
          f"{'Expected %':>10} {'Wrong ref in other barcodes':>30}")

    total_aligned = 0
    for ref in sorted(ref_total, key=lambda x: -ref_total[x]):
        total = ref_total[ref]
        total_aligned += total
        # Find which barcode should have these reads
        expected_barcode = BARCODE_MAP.get(ref, "UNKNOWN")
        # Reads correctly assigned to expected barcode
        correct_in_expected = 0
        wrong_refs_list = []

        for barcode, qnames in barcode_qnames.items():
            matching = sum(1 for qn in qnames if alignment_map.get(qn, "") == ref)
            if barcode == expected_barcode:
                correct_in_expected = matching
            elif matching > 0:
                wrong_refs_list.append((barcode, matching))

        n_wrong = total - correct_in_expected

        correct_pct = 100.0 * correct_in_expected / total if total > 0 else 0

        wrong_str = ", ".join(f"{b}={c}" for b, c in sorted(wrong_refs_list))
        if not wrong_str:
            wrong_str = "(none)"

        print(
            f"{ref:<12} {total:>14,} "
            f"{expected_barcode:>18} {correct_in_expected:>14,} "
            f"{correct_pct:>9.2f}%  {wrong_str}"
        )

    if total_aligned > 0:
        print(f"\n{'─' * 60}")
        print(f"OVERALL:")
        print(
            f"  Total reads aligned to references: {total_aligned:,}"
        )

    # Also compute per-barcode alignment accuracy
    print(HEADER)
    print("  SCENARIO 2b: For each barcode, what fraction of its reads are correctly aligned?")
    print(HEADER)

    n_total_reads = sum(len(qs) for qs in barcode_qnames.values())
    total_correct_align = 0
    total_wrong_align = 0
    total_noref_align = 0

    for barcode, qnames in sorted(barcode_qnames.items()):
        expected_ref = BARCODE_TO_REF.get(barcode, "UNKNOWN")
        n_total = len(qnames)
        n_correct = 0
        n_wrong = 0
        n_no_ref = 0
        for qn in qnames:
            rname = alignment_map.get(qn, "")
            if not rname:
                n_no_ref += 1
            elif rname == expected_ref:
                n_correct += 1
            else:
                n_wrong += 1

        total_correct_align += n_correct
        total_wrong_align += n_wrong
        total_noref_align += n_no_ref

        if n_total > 0:
            print(
                f"  {barcode}: correct={n_correct:>8,} ({100*n_correct/n_total:.2f}%)"
                f"  wrong-ref={n_wrong:>6,} ({100*n_wrong/n_total:.2f}%)"
                f"  no-ref={n_no_ref:>6,} ({100*n_no_ref/n_total:.2f}%)"
            )

    print(f"\n{'─' * 60}")
    if n_total_reads > 0:
        align_error_rate = (total_wrong_align + total_noref_align) / n_total_reads * 100
        print(
            f"OVERALL: assuming barcode demux is correct, alignment error rate = "
            f"{(total_wrong_align + total_noref_align):,} errors out of {n_total_reads:,} reads ({align_error_rate:.2f}%)"
        )
        print(f"  - Wrong reference:     {total_wrong_align:,}")
        print(f"  - No alignment:        {total_noref_align:,}")


def print_sample_qnames(
    barcode_qnames: dict[str, list[str]],
    alignment_map: dict[str, str],
) -> None:
    """Print sample qnames for manual verification."""
    print(HEADER)
    print("  Sample qnames for manual verification")
    print(HEADER)

    categories = {
        "Aligned to expected ref": [],
        "Cross-contaminated (barcode->wrong ref)": [],
        "In CSV but unaligned (empty rname)": [],
        "Not in CSV (missing alignment record)": [],
    }

    for barcode, qnames in sorted(barcode_qnames.items()):
        expected_ref = BARCODE_TO_REF.get(barcode, "UNKNOWN")
        for qn in qnames:
            if qn not in alignment_map:
                if len(categories["Not in CSV (missing alignment record)"]) < 5:
                    categories["Not in CSV (missing alignment record)"].append(
                        f"{barcode} -> qname NOT in CSV: {qn}"
                    )
                continue

            rname = alignment_map[qn]

            if not rname:
                if len(categories["In CSV but unaligned (empty rname)"]) < 5:
                    categories["In CSV but unaligned (empty rname)"].append(
                        f"{barcode} -> empty rname in CSV: {qn}"
                    )
                continue

            if rname == expected_ref:
                if len(categories["Aligned to expected ref"]) < 5:
                    categories["Aligned to expected ref"].append(
                        f"{barcode} -> {rname}: {qn}"
                    )
            elif len(categories["Cross-contaminated (barcode->wrong ref)"]) < 10:
                categories["Cross-contaminated (barcode->wrong ref)"].append(
                    f"{barcode} (expected {expected_ref}) -> aligned to {rname}: {qn}"
                )

    for cat, samples in categories.items():
        if samples:
            print(f"\n  {cat}:")
            for s in samples:
                print(f"    {s}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cross-validate alignment vs barcode demultiplexing"
    )
    parser.add_argument(
        "--alignment-csv",
        type=Path,
        default="/data1/ccs_data/str-optimization/20260723_250302Y0001_Run0003/"
                "20260723_250302Y0001_Run0003.smc_all_reads-metric/"
                "20260723_250302Y0001_Run0003.smc_all_reads.gsmm2_aligned_metric_fact.csv",
        help="Path to the main alignment metric CSV (tab-separated)",
    )
    parser.add_argument(
        "--barcode-dir",
        type=Path,
        default="/data1/ccs_data/str-optimization/20260723_250302Y0001_Run0003/barcode_assign",
        help="Directory containing barcode fastq/bam files",
    )
    parser.add_argument(
        "--n-sample",
        type=int,
        default=-1,
        help="Limit reads per barcode file (default: all). Set > 0 for sampling.",
    )
    args = parser.parse_args()

    # Validate inputs
    if not args.alignment_csv.exists():
        print(f"ERROR: alignment CSV not found: {args.alignment_csv}", file=sys.stderr)
        sys.exit(1)
    if not args.barcode_dir.is_dir():
        print(f"ERROR: barcode dir not found: {args.barcode_dir}", file=sys.stderr)
        sys.exit(1)

    # Step 1: Load alignment data
    print("=" * 72, file=sys.stderr)
    print("Loading alignment CSV ...", file=sys.stderr)
    alignment_map = load_alignment_csv(args.alignment_csv)
    total_in_csv = len(alignment_map)

    n_aligned = sum(1 for r in alignment_map.values() if r)
    n_unaligned = total_in_csv - n_aligned
    print(f"  Total reads:       {total_in_csv:>10,}", file=sys.stderr)
    print(f"  Aligned:           {n_aligned:>10,}", file=sys.stderr)
    print(f"  Unaligned (empty): {n_unaligned:>10,}", file=sys.stderr)

    # Step 2: Load barcode qnames
    print("\nLoading barcode fastq files ...", file=sys.stderr)
    barcode_qnames = collect_barcode_qnames(args.barcode_dir, args.n_sample)
    total_barcode_reads = sum(len(qs) for qs in barcode_qnames.values())
    print(f"  Total reads across all barcodes: {total_barcode_reads:,}", file=sys.stderr)

    # Step 3: Cross-validation
    cross_validate_alignment_correct(barcode_qnames, alignment_map)
    cross_validate_barcode_correct(barcode_qnames, alignment_map)
    print_sample_qnames(barcode_qnames, alignment_map)

    # Final summary
    print(HEADER)
    print("  SUMMARY")
    print(HEADER)
    n_total_barcode = sum(len(qs) for qs in barcode_qnames.values())
    if n_total_barcode > 0:
        demux_error_rate = 0
        total_wrong_from_align = 0
        for barcode, qnames in barcode_qnames.items():
            expected_ref = BARCODE_TO_REF.get(barcode, "UNKNOWN")
            total_wrong_from_align += sum(
                1 for qn in qnames if alignment_map.get(qn, "") != ""
                               and alignment_map.get(qn, "") != expected_ref
            )
        demux_error_rate = 100.0 * total_wrong_from_align / n_total_barcode

        print(f"  Scenario 1 - Assuming alignment is correct:")
        print(f"    Barcode demux error rate:          {demux_error_rate:.2f}%")
        print(f"\n  Scenario 2 - Assuming barcode demux is correct:")
        align_error = 0
        for barcode, qnames in barcode_qnames.items():
            expected_ref = BARCODE_TO_REF.get(barcode, "UNKNOWN")
            align_error += sum(
                1 for qn in qnames
                if alignment_map.get(qn, "") != ""
                   and alignment_map.get(qn, "") != expected_ref
            ) + sum(
                1 for qn in qnames
                if not alignment_map.get(qn, "")
            )
        print(f"    Alignment error rate:              {100.0 * align_error / n_total_barcode:.2f}%")

    # Read overlap between CSV and barcode files
    print(f"\n  Read overlap:")
    csv_reads = set(alignment_map.keys())
    barcode_reads = set()
    for qnames in barcode_qnames.values():
        barcode_reads.update(qnames)
    both = csv_reads & barcode_reads
    only_csv = csv_reads - barcode_reads
    only_barcode = barcode_reads - csv_reads
    print(f"    Reads in CSV but not in barcodes:  {len(only_csv):,}")
    print(f"    Reads in barcodes but not in CSV:  {len(only_barcode):,}")
    print(f"    Overlapping reads:                 {len(both):,}")


if __name__ == "__main__":
    main()
