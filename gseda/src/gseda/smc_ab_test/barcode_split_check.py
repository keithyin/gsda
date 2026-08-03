#!/usr/bin/env python3
"""AB check: compare query_names between v4 and v5 demuxed FASTQ outputs.

For each filename present in both directories, extract the read names
(query_name from @ lines) and report whether they match exactly.

Usage:
    python barcode_split_check.py <v4_dir> <v5_dir> [--csv output.csv]
"""

import argparse
import csv
from pathlib import Path




def extract_names(fastq_path: Path) -> tuple[set[str], int]:
    """Stream through a FASTQ file, collecting query_name prefix (partition ID).

    Returns names as <run>/<partition>, e.g. '20260706_260602Y0005_Run0001/477153'.
    Returns (names_set, record_count).
    """
    names = set()
    count = 0
    with fastq_path.open("r") as fh:
        header_line = fh.readline()
        while header_line:
            if header_line.startswith("@"):
                name = header_line[1:].strip().split()[0]
                # Strip coordinate suffix: "run/partition/32-2429" -> "run/partition"
                parts = name.rsplit("/", 2)
                key = "/".join(parts[:-1]) if len(parts) >= 3 else name
                names.add(key)
                count += 1
                # skip the next 3 lines (seq, +, qual)
                fh.readline()
                fh.readline()
                fh.readline()
            else:
                fh.readline()
                fh.readline()
                fh.readline()
                fh.readline()
            header_line = fh.readline()
    return names, count


def main():
    parser = argparse.ArgumentParser(description="AB check v4 vs v5 demuxed FASTQ")
    parser.add_argument("dir_v4", help="v4 demuxed FASTQ directory")
    parser.add_argument("dir_v5", help="v5 demuxed FASTQ directory")
    parser.add_argument("--csv", help="Export per-file results to CSV")
    args = parser.parse_args()

    dir_v4 = Path(args.dir_v4)
    dir_v5 = Path(args.dir_v5)

    fastq_files_v4 = {f.name for f in dir_v4.glob("*.fastq")}
    fastq_files_v5 = {f.name for f in dir_v5.glob("*.fastq")}

    common = sorted(fastq_files_v4 & fastq_files_v5)
    only_v4 = sorted(fastq_files_v4 - fastq_files_v5)
    only_v5 = sorted(fastq_files_v5 - fastq_files_v4)

    # Per-file results
    results = []  # list of dicts
    identical_count = 0
    different_count = 0
    total_records_v4 = 0
    total_records_v5 = 0
    total_only_in_v4 = 0
    total_only_in_v5 = 0

    for fname in common:
        path_v4 = dir_v4 / fname
        path_v5 = dir_v5 / fname

        names_v4, count_v4 = extract_names(path_v4)
        names_v5, count_v5 = extract_names(path_v5)
        total_records_v4 += count_v4
        total_records_v5 += count_v5

        if names_v4 == names_v5 and count_v4 == count_v5:
            identical_count += 1
            status = "IDENTICAL"
        else:
            different_count += 1
            only_in_v4 = names_v4 - names_v5
            only_in_v5 = names_v5 - names_v4
            shared = names_v4 & names_v5
            total_only_in_v4 += len(only_in_v4)
            total_only_in_v5 += len(only_in_v5)
            # Check if the same records exist but order differs (same set, diff count)
            if names_v4 == names_v5:
                status = "SAME_NAMES_DIFF_COUNT"
            else:
                status = "DIFFERENT"

            results.append({
                "filename": fname,
                "status": status,
                "v4_records": count_v4,
                "v5_records": count_v5,
                "shared_records": len(shared),
                "only_in_v4": len(only_in_v4),
                "only_in_v5": len(only_in_v5),
            })

    # --- Print summary ---
    print("=" * 70)
    print("AB CHECK SUMMARY: v4 vs v5 demuxed FASTQ")
    print("=" * 70)
    print(f"\nTotal common filenames:       {len(common)}")
    print(f"Only in v4:                   {len(only_v4)}")
    print(f"Only in v5:                   {len(only_v5)}")
    print()
    print(f"Identical query_names:        {identical_count}")
    print(f"Different query_names:        {different_count}")
    print(f"\nTotal records (v4):           {total_records_v4:,}")
    print(f"Total records (v5):           {total_records_v5:,}")
    total_divergent = total_only_in_v4 + total_only_in_v5
    total_all = total_records_v4 + total_records_v5
    if total_divergent > 0:
        pct = total_divergent / total_all * 100
        print(f"Difference in total records:  {total_divergent:,} ({pct:.2f}%)")
        print(f"  (v4-only={total_only_in_v4:,}, v5-only={total_only_in_v5:,}, "
              f"total pool={total_all:,})")

    # --- Files with only one side ---
    if only_v4 or only_v5:
        print()
        print("-" * 70)
        print("FILES ONLY IN ONE SIDE (no pair to compare):")
        for f in only_v4[:10]:
            path = dir_v4 / f
            _, cnt = extract_names(path)
            print(f"  [v4 only] {f} ({cnt:,} records)")
        if len(only_v4) > 10:
            print(f"  ... and {len(only_v4) - 10} more")
        for f in only_v5[:10]:
            path = dir_v5 / f
            _, cnt = extract_names(path)
            print(f"  [v5 only] {f} ({cnt:,} records)")
        if len(only_v5) > 10:
            print(f"  ... and {len(only_v5) - 10} more")

    # --- Details for different files ---
    if results:
        print()
        print("-" * 70)
        print("FILES WITH DIFFERENCES:")
        print(f"\n{'Filename':<40} {'Status':<25} {'v4 recs':>8} {'v5 recs':>8} {'Shared':>8} {'v4-only':>8} {'v5-only':>8}")
        print("-" * 110)
        for r in results:
            print(
                f"{r['filename']:<40} {r['status']:<25} "
                f"{r['v4_records']:>8,} {r['v5_records']:>8,} "
                f"{r['shared_records']:>8,} {r['only_in_v4']:>8,} {r['only_in_v5']:>8,}"
            )

        # Per-file breakdown for files with actual name differences
        differing = [r for r in results if r["status"] == "DIFFERENT"]
        if differing:
            print()
            print("-" * 70)
            print("DETAILED BREAKDOWN FOR FILES WITH DIFFERENCES:")
            for r in differing[:5]:  # top 5 by total diff
                fname = r["filename"]
                path_v4 = dir_v4 / fname
                path_v5 = dir_v5 / fname
                names_v4, count_v4 = extract_names(path_v4)
                names_v5, count_v5 = extract_names(path_v5)
                only_in_v4 = names_v4 - names_v5
                only_in_v5 = names_v5 - names_v4

                print(f"\n  {fname} (v4={count_v4:,}, v5={count_v5:,})")
                if only_in_v4:
                    # Show up to 5 sample names unique to v4
                    samples = sorted(only_in_v4)[:5]
                    print(f"    v4-only ({len(only_in_v4)} reads): first 5 = {samples}")
                if only_in_v5:
                    samples = sorted(only_in_v5)[:5]
                    print(f"    v5-only ({len(only_in_v5)} reads): first 5 = {samples}")

    # --- CSV export ---
    if args.csv:
        csv_path = Path(args.csv)
        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "filename", "status", "v4_records", "v5_records",
                "shared_records", "only_in_v4", "only_in_v5",
            ])
            for r in results:
                writer.writerow([
                    r["filename"], r["status"], r["v4_records"], r["v5_records"],
                    r["shared_records"], r["only_in_v4"], r["only_in_v5"],
                ])
        print(f"\nDetailed results exported to {csv_path}")

    print()
    print("=" * 70)
    if different_count == 0 and not only_v4 and not only_v5:
        print("RESULT: PERFECT MATCH - v4 and v5 outputs are identical!")
    elif different_count == 0:
        print(f"RESULT: All {identical_count} paired files have matching query_names.")
    else:
        print(f"RESULT: Found differences in {different_count}/{len(common)} files.")
    print("=" * 70)


if __name__ == "__main__":
    main()
