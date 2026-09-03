#!/usr/bin/env python3
"""
Barcode-split FASTQ -> reference alignment report.

For each plasmid in a plasmid->barcode mapping that has a matching reference
sequence in a multi-record reference FASTA, this script:
  1. extracts the per-plasmid reference,
  2. converts the barcode-split FASTQ to an unaligned BAM (with an `rq` tag),
  3. runs gseda.ppl.sequencing_report_v2 to align query vs reference and emit
     aggregated metrics,
  4. writes a cross-sample summary TSV.

Everything is run inside a given conda env (default: py38).
"""
import argparse
import csv
import glob
import os
import subprocess
import sys

REPO = "/root/projects/gsda/gseda/src/gseda"


def log(*a):
    print(*a, flush=True)


def build_plasmid_barcode_map(mapping_path):
    """Read the plasmid->barcode TSV (header: plasmid<TAB>barcode)."""
    m = {}
    with open(mapping_path) as f:
        rdr = csv.reader(f, delimiter="\t")
        next(rdr, None)  # header
        for row in rdr:
            if len(row) < 2 or not row[0].strip():
                continue
            m[row[0].strip()] = row[1].strip()
    return m


def list_ref_records(ref_fasta, prefix):
    """Return list of (ref_id, plasmid) for each reference record."""
    out = []
    with open(ref_fasta) as f:
        for line in f:
            if line.startswith(">"):
                rid = line.strip().lstrip(">")
                pl = rid.replace(prefix, "") if prefix else rid
                out.append((rid, pl))
    return out


def extract_ref(ref_fasta, plasmid, prefix, out_fasta):
    """Extract one record (matched by plasmid id) into out_fasta. Returns True on success."""
    snippet = (
        "import sys\n"
        "from Bio import SeqIO\n"
        "want, src, dst, pre = sys.argv[1:5]\n"
        "rs = [r for r in SeqIO.parse(src, 'fasta') if r.id.replace(pre, '') == want]\n"
        "if not rs:\n"
        "    sys.exit(1)\n"
        "SeqIO.write(rs, dst, 'fasta')\n"
    )
    r = subprocess.run(["python3", "-c", snippet, plasmid, ref_fasta, out_fasta, prefix or ""],
                       capture_output=True, text=True)
    return r.returncode == 0


def conda_run(env, cmd):
    return subprocess.run(["conda", "run", "-n", env] + cmd,
                          capture_output=True, text=True)


def main():
    ap = argparse.ArgumentParser(description="Barcode FASTQ vs reference alignment report")
    ap.add_argument("--barcode-dir", required=True,
                    help="directory containing BarcodeNN.fastq files")
    ap.add_argument("--ref-fasta", required=True,
                    help="multi-record reference fasta (e.g. intersect.fa)")
    ap.add_argument("--mapping", required=True,
                    help="plasmid<TAB>barcode TSV with a header line")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--rq-range", default="0.99:1.1",
                    help="rq filter range for report, e.g. 0.99:1.1")
    ap.add_argument("--conda-env", default="py38", help="conda env (default py38)")
    ap.add_argument("--ref-prefix", default="260827STR-",
                    help="prefix stripped from reference ids to get plasmid name")
    ap.add_argument("--short-aln", type=int, default=1,
                    help="1 = short-alignment mode (query or target in [30,200])")
    ap.add_argument("--fastq2bam", default=f"{REPO}/file_format_cvt/fastq2bam.py")
    ap.add_argument("--report", default=f"{REPO}/ppl/sequencing_report_v2.py")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    refd = os.path.join(args.outdir, "refs")
    os.makedirs(refd, exist_ok=True)

    p2b = build_plasmid_barcode_map(args.mapping)
    refs = {pl: rid for rid, pl in list_ref_records(args.ref_fasta, args.ref_prefix)}

    # barcode number -> plasmid (only plasmids that have a reference)
    todo = []
    for pl, bc in p2b.items():
        if pl in refs:
            todo.append((bc, pl))
        else:
            log(f"SKIP plasmid={pl} barcode={bc}: no reference in {args.ref_fasta}")
    todo.sort(key=lambda x: int(x[0]))
    log(f"Running {len(todo)} samples (of {len(p2b)} mapped plasmids).")

    for bc, pl in todo:
        fq = os.path.join(args.barcode_dir, f"Barcode{int(bc):02d}.fastq")
        if not os.path.exists(fq):
            log(f"SKIP plasmid={pl} barcode={bc}: missing {fq}")
            continue
        ref = os.path.join(refd, f"{pl}.fa")
        if not extract_ref(args.ref_fasta, pl, args.ref_prefix, ref):
            log(f"FAIL {pl}: reference extraction failed")
            continue
        bam = os.path.join(args.outdir, f"Barcode{bc}_{pl}.bam")
        log(f"##### plasmid={pl} barcode={bc}")
        r1 = conda_run(args.conda_env,
                       ["python", args.fastq2bam, fq, bam, "--ref", ref])
        if r1.returncode != 0:
            log(f"FAIL {pl}: fastq2bam -> {r1.stderr.strip()}")
            continue
        rpt = ["python", args.report, "--bams", bam, "--refs", ref,
               "--rq-range", args.rq_range]
        if args.short_aln:
            rpt += ["--short-aln", str(args.short_aln)]
        r2 = conda_run(args.conda_env, rpt)
        if r2.returncode != 0:
            log(f"FAIL {pl}: report -> {r2.stderr.strip()}")
            continue
        log(f"done {pl} (Barcode{bc})")

    # summary
    files = sorted(glob.glob(os.path.join(args.outdir, "*-metric",
                                          "*.gsmm2_aligned_metric_aggr.csv")))
    if not files:
        log("No result files produced.")
        return
    sel = ["reads_num", "tot_bases", "n50", "read_len_p50",
           "alignedRatio", "notAlignedRatio", "queryCoverage", "queryCoverage3",
           "identity", "identity-p50", "mmRate", "longIndelRatio",
           "GlobalQueryCoverage", "identity≥0.83", "identity≥0.90", "identity≥0.99"]

    def load(f):
        return {r["name"]: r["value"]
                for r in csv.DictReader(open(f), delimiter="\t")}

    def fnum(d, k):
        try:
            return f"{float(d[k]):.6f}"
        except (KeyError, ValueError, TypeError):
            return ""

    out = os.path.join(args.outdir, "summary.tsv")
    with open(out, "w") as fh:
        fh.write("sample\t" + "\t".join(sel) + "\n")
        for f in files:
            stem = os.path.basename(f).replace(".gsmm2_aligned_metric_aggr.csv", "")
            d = load(f)
            fh.write(stem + "\t" + "\t".join(fnum(d, k) for k in sel) + "\n")
    log(f"Wrote {out}")
    log("ALL DONE")


if __name__ == "__main__":
    main()
