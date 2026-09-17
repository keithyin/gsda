#!/usr/bin/env python3
"""
Prepare per-run inputs for the barcode_ref_align pipeline on the 2nd STR batch.

Reads the run-aware mapping (plasmid<TAB>barcode<TAB>RUN号, barcode written as
"24标签-N") and, per run:
  - writes a numeric plasmid<TAB>barcode TSV that run_report.py understands,
  - builds a BarcodeNN.fastq view of that run's query dir. If the query dir is
    already named BarcodeNN.fastq it links those; otherwise it resolves
    "24标签-N" -> barcode-NN -> "<Adaptor-...>.fastq" via barcode-2-barcodename.tsv
    and links the amplicon file. Links, not copies (query dirs are ~1 GB/run).

Prints a per-run manifest and exits non-zero if any mapped sample is missing.
"""
import argparse
import csv
import os
import re
import sys

BARCODE2NAME_TSV = "/data1/ccs_data/str-optimization/barcode-2-barcodename.tsv"


def load_barcode2orig(path):
    """{barcode-NN: Adaptor-...} from the barcodename table (col1=orig, col2=barcode-NN)."""
    m = {}
    with open(path) as f:
        rdr = csv.reader(f, delimiter="\t")
        next(rdr, None)
        for row in rdr:
            if len(row) < 2:
                continue
            orig, name = row[0].strip(), row[1].strip()
            m[name] = orig
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping", required=True, help="plasmid<TAB>barcode<TAB>RUN号 TSV")
    ap.add_argument("--root", required=True, help="scratch root for per-run inputs")
    ap.add_argument("--run-query", action="append", required=True, metavar="RUN=DIR",
                    help="RUN号=query dir, repeatable")
    args = ap.parse_args()

    queries = dict()
    for item in args.run_query:
        run, _, d = item.partition("=")
        queries[run.strip()] = d.strip()

    b2o = load_barcode2orig(BARCODE2NAME_TSV)

    rows = list(csv.DictReader(open(args.mapping), delimiter="\t"))
    per_run = {}
    for r in rows:
        pl, bc, run = r["plasmid"].strip(), r["barcode"].strip(), r["RUN号"].strip()
        if not pl:
            continue
        n = int(re.search(r"(\d+)\s*$", bc).group(1))
        per_run.setdefault(run, []).append((n, pl, bc))

    missing = []
    for run, items in per_run.items():
        if run not in queries:
            print(f"NOTE run={run}: no query dir supplied, skipped")
            continue
        qdir = queries[run]
        rd = os.path.join(args.root, run)
        view = os.path.join(rd, "barcode_view")
        os.makedirs(view, exist_ok=True)
        mpath = os.path.join(rd, "plasmid_2_barcode.tsv")

        with open(mpath, "w") as fh:
            fh.write("plasmid\tbarcode\n")
            for n, pl, bc in sorted(items):
                fh.write(f"{pl}\t{n}\n")

        # already BarcodeNN-named?
        direct = os.path.join(qdir, f"Barcode{1:02d}.fastq")
        print(f"### run={run} samples={len(items)} query={qdir}")
        for n, pl, bc in sorted(items):
            link = os.path.join(view, f"Barcode{n:02d}.fastq")
            if os.path.exists(direct):
                src = os.path.join(qdir, f"Barcode{n:02d}.fastq")
                how = "as-is"
            else:
                orig = b2o.get(f"barcode-{n:02d}")
                if orig is None:
                    print(f"  MISSING barcode-{n:02d} in {BARCODE2NAME_TSV} ({pl})")
                    missing.append((run, pl))
                    continue
                src = os.path.join(qdir, orig + ".fastq")
                how = f"<- {orig}"
            if not os.path.exists(src):
                print(f"  MISSING {src}  ({pl} / {bc})")
                missing.append((run, pl))
                continue
            if os.path.islink(link) or os.path.exists(link):
                os.remove(link)
            os.symlink(os.path.abspath(src), link)
            sz = os.path.getsize(src)
            print(f"  Barcode{n:02d}.fastq {how:32s} {pl:8s} {bc:10s} {sz/1e6:8.1f} MB")
        print(f"  mapping -> {mpath}\n  view    -> {view}")

    if missing:
        print(f"\n{len(missing)} mapped sample(s) could not be resolved:", file=sys.stderr)
        for run, pl in missing:
            print(f"  {run} {pl}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
