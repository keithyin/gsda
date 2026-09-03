#!/usr/bin/env python3
"""Rename barcode-split FASTQ files by a name mapping.

Reads a mapping TSV (col1 = original file base name without .fastq,
col2 = new barcodename) and copies the files that exist in the source
directory into an output directory, renamed to the mapped barcodename.
The source directory is never modified.
"""
import argparse
import os
import shutil
import sys

# Embedded default mapping (barcode-2-barcodename.tsv). Columns are
# tab-separated: col1 = original file base (matches source <name>.fastq),
# col2 = target barcodename, col3 = barcode seq (unused by this script).
# The source TSV's header is mislabeled ("barcodename barcode seq"); the
# data rows above are authoritative — col1 holds the Adaptor-... name.
DEFAULT_MAPPING = """\
Adaptor-barcode251-4\tbarcode-01\tATTGGAGTATCTTTCTGCAGCCCTCCTCCTCCTCCGTTTGTCTGAGCAGAAAGATACTCCAAT
Adaptor-barcode250-4\tbarcode-02\tATAACATTCTACGGTAACAGCCCTCCTCCTCCTCCGTTTGTCTGAGTTACCGTAGAATGTTAT
Adaptor-barcode271-3\tbarcode-03\tATTTCACGTAAGTAATCGGCTTATCCTCCTCCTCCGTTACATACGCGATTACTTACGTGAAAT
Adaptor-barcode266-2\tbarcode-04\tATCTCGGCGGAAAACCTACACCCTCCTCCTCCTCCGTTGCCGCCATAGGTTTTCCGCCGAGAT
Adaptor-barcode229-0\tbarcode-05\tATTGTCTGCGAAATTGGTTTTTTTCCTCCTCCTCCGTTGTTGTTTACCAATTTCGCAGACAAT
Adaptor-barcode255-3\tbarcode-06\tATGAACAACCTAGATGCTGCTTATCCTCCTCCTCCGTTACATACGAGCATCTAGGTTGTTCAT
Adaptor-barcode261-3\tbarcode-07\tATAGAATTAATGCACCACGCTTATCCTCCTCCTCCGTTACATACGGTGGTGCATTAATTCTAT
Adaptor-barcode206-3\tbarcode-08\tATAATTTTCCAAAACCGGGCTTATCCTCCTCCTCCGTTACATACGCCGGTTTTGGAAAATTAT
Adaptor-barcode230-2\tbarcode-09\tATTCCACACTCTCACCACCACCCTCCTCCTCCTCCGTTGCCGCCAGTGGTGAGAGTGTGGAAT
Adaptor-barcode256-3\tbarcode-10\tATGTCTCATGGTGTGAGTGCTTATCCTCCTCCTCCGTTACATACGACTCACACCATGAGACAT
Adaptor-barcode207-4\tbarcode-11\tATCCGGTTGGCCGGAATTAGCCCTCCTCCTCCTCCGTTTGTCTGAAATTCCGGCCAACCGGAT
Adaptor-barcode222-2\tbarcode-12\tATCCCCTTTTGGAATTCCCACCCTCCTCCTCCTCCGTTGCCGCCAGGAATTCCAAAAGGGGAT
Adaptor-barcode270-0\tbarcode-13\tATGGTTCCACACGCAACGTTTTTTCCTCCTCCTCCGTTGTTGTTTCGTTGCGTGTGGAACCAT
Adaptor-barcode201-5\tbarcode-14\tATGGCCTTCCTTTTGGCCCCATTTCCTCCTCCTCCGTTCTGCCAAGGCCAAAAGGAAGGCCAT
Adaptor-barcode253-4\tbarcode-15\tATATGTCTATCCCAAGCAAGCCCTCCTCCTCCTCCGTTTGTCTGATGCTTGGGATAGACATAT
Adaptor-barcode217-2\tbarcode-16\tATCCGGTTCCAAGGGGCCCACCCTCCTCCTCCTCCGTTGCCGCCAGGCCCCTTGGAACCGGAT
Adaptor-barcode214-1\tbarcode-17\tATTTCCAACCGGAATTAAACAACTCCTCCTCCTCCGTTCATGAACTTAATTCCGGTTGGAAAT
Adaptor-barcode221-3\tbarcode-18\tATTTGGCCTTTTCCGGGGGCTTATCCTCCTCCTCCGTTACATACGCCCCGGAAAAGGCCAAAT
Adaptor-barcode216-0\tbarcode-19\tATAACCCCAAGGCCGGTTTTTTTTCCTCCTCCTCCGTTGTTGTTTAACCGGCCTTGGGGTTAT
Adaptor-barcode254-5\tbarcode-20\tATCCAGAGTCAGTCATTCCCATTTCCTCCTCCTCCGTTCTGCCAAGAATGACTGACTCTGGAT
Adaptor-barcode220-0\tbarcode-21\tATAATTGGTTCCTTAATTTTTTTTCCTCCTCCTCCGTTGTTGTTTAATTAAGGAACCAATTAT
Adaptor-barcode225-1\tbarcode-22\tATCAACGTATAAACGTCGACAACTCCTCCTCCTCCGTTCATGAACCGACGTTTATACGTTGAT
Adaptor-barcode233-0\tbarcode-23\tATCTGTTAGGATCGCAGATTTTTTCCTCCTCCTCCGTTGTTGTTTTCTGCGATCCTAACAGAT
Adaptor-barcode239-2\tbarcode-24\tATTGAAGCGTTGCGTGCGCACCCTCCTCCTCCTCCGTTGCCGCCACGCACGCAACGCTTCAAT
"""


def rename_target(name):
    """Normalize a target barcodename to the output form.

    e.g. "barcode-01" -> "Barcode01". Dashes are dropped and the first
    letter is capitalized; the rest is left as-is.
    """
    compact = name.replace("-", "")
    return compact[:1].upper() + compact[1:]


def parse_mapping(text):
    """Return {orig_base: new_name} from mapping text.

    col1 = original file base (no .fastq), col2 = new barcodename.
    Extra columns (e.g. a seq column) and a header line are ignored.
    """
    mapping = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2 or not parts[0]:
            continue
        # skip a header row: its target column has no digits (e.g. "barcode"),
        # whereas real targets like "barcode-01" do.
        if not any(c.isdigit() for c in parts[1]):
            continue
        mapping[parts[0].strip()] = parts[1].strip()
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", required=True,
                    help="directory containing <orig>.fastq files")
    ap.add_argument("--tsv", default=None,
                    help="mapping TSV (col1=original base, col2=new name). "
                         "Omit to use the built-in default mapping.")
    ap.add_argument("--outdir", required=True,
                    help="output directory for renamed copies")
    ap.add_argument("--include-unmapped", action="store_true",
                    help="also copy files with no mapping, keeping their name")
    args = ap.parse_args()

    if args.tsv:
        with open(args.tsv) as f:
            mapping = parse_mapping(f.read())
    else:
        mapping = parse_mapping(DEFAULT_MAPPING)

    src_files = sorted(
        f for f in os.listdir(args.src) if f.endswith(".fastq")
    )
    src_bases = {f[:-len(".fastq")] for f in src_files}

    os.makedirs(args.outdir, exist_ok=True)

    mapped, unmapped = [], []
    for base in sorted(src_bases):
        if base in mapping:
            new_name = rename_target(mapping[base]) + ".fastq"
            shutil.copy2(os.path.join(args.src, base + ".fastq"),
                         os.path.join(args.outdir, new_name))
            mapped.append((base, new_name))
        elif args.include_unmapped:
            shutil.copy2(os.path.join(args.src, base + ".fastq"),
                         os.path.join(args.outdir, base + ".fastq"))
            unmapped.append(base)

    missing = [b for b in mapping if b not in src_bases]

    print(f"source files:            {len(src_bases)}")
    print(f"mapping entries:         {len(mapping)}")
    print(f"renamed (copied):        {len(mapped)}")
    print(f"unmapped (kept as-is):   {len(unmapped)}")
    print(f"mapped-but-missing:      {len(missing)}")
    print(f"total in outdir:         {len(os.listdir(args.outdir))}")
    if missing:
        print("  mapped but no source file:")
        for m in missing:
            print(f"    {m} -> {rename_target(mapping[m])}")
    print(f"outdir: {args.outdir}")


if __name__ == "__main__":
    sys.exit(main())
