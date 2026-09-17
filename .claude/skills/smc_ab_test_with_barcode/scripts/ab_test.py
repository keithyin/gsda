#!/usr/bin/env python3
"""A/B comparison of two barcode-call outputs via sequencing_report_v2.

Given two per-barcode demux directories (e.g. two barcode-calling models over
the same SMC BAM), the plasmid->barcode mapping, the per-plasmid Sanger
references, and a target RUN号, this:
  1. for every barcode the mapping covers, converts that model's
     BarcodeNN.fastq to an unaligned BAM (rq tag) via fastq2bam.py,
  2. runs gseda.ppl.sequencing_report_v2 against the plasmid reference,
  3. collects the per-barcode aggr metrics for both models,
  4. aggregates all barcodes per model with the report's own
     --fact-csvs/--basic-csvs merge (merge_partition) logic — a read/bases
     weighted global metric, NOT an arithmetic mean of per-barcode metrics.

Outputs (under --outdir):
  <labelA>/  <labelB>/      per-barcode BAMs + *-metric/ dirs + _agg/
  per_barcode.tsv           long: barcode x metric x (A, B, delta)
  aggregate.csv             pooled-metric x (A, B, delta)
  aggregate_<labelA>.csv    full per-model aggregate
  aggregate_<labelB>.csv
  summary.md                human-readable comparison

The two demux dirs must be named BarcodeNN.fastq (NN = the numeric part of
"24标签-NN" in the mapping). Any BarcodeNN.fastq in a demux dir that the
mapping does not cover is silently skipped.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys

GS = "/root/miniconda3/envs/gseda/bin/python"
REPO = "/root/projects/gsda/third_party/gseda/src/gseda"
FASTQ2BAM = f"{REPO}/file_format_cvt/fastq2bam.py"
REPORT = f"{REPO}/ppl/sequencing_report_v2.py"

# metrics shown in the per-barcode / aggregate comparison (read from name\tvalue)
SEL = [
    "reads_num", "tot_bases", "n50", "read_len_p50",
    "alignedRatio", "notAlignedRatio",
    "queryCoverage", "queryCoverage2", "queryCoverage3",
    "identity", "identity-p50", "mmRate",
    "HomoInsRate", "HomoDelRate",
    "NHInsRate", "NHDelRate", "longIndelRatio", "GlobalQueryCoverage",
    "identity≥0.83", "identity≥0.90", "identity≥0.99",
]


def log(*a):
    print(*a, flush=True)


def run(cmd):
    log("$", " ".join(cmd))
    return subprocess.call(cmd)


def load_plasmid_barcode(path, run):
    """{barcode_int: plasmid} for the given RUN号 (barcode column '24标签-N')."""
    m = {}
    with open(path) as f:
        rdr = csv.reader(f, delimiter="\t")
        next(rdr, None)  # header
        for row in rdr:
            if len(row) < 3:
                continue
            pl, bc, r = (x.strip() for x in row[:3])
            if not pl or r != run:
                continue
            try:
                n = int(bc.rsplit("-", 1)[1])
            except (IndexError, ValueError):
                log(f"WARN cannot parse barcode '{bc}'")
                continue
            m[n] = pl
    return m


def read_kv(path):
    d = {}
    with open(path) as f:
        for r in csv.reader(f, delimiter="\t"):
            if len(r) >= 2 and r[0] != "name":
                d[r[0]] = r[1]
    return d


def fval(d, k):
    try:
        return float(d[k])
    except (KeyError, ValueError, TypeError):
        return None


def fmt(v, sign=False):
    if v is None:
        return ""
    return f"{v:+.6g}" if sign else f"{v:.6g}"


def main():
    ap = argparse.ArgumentParser(description="A/B compare two barcode-call outputs via sequencing_report_v2")
    ap.add_argument("--a", required=True, metavar="LABEL=DEMUXDIR",
                    help="model A: label=demuxed dir (containing BarcodeNN.fastq)")
    ap.add_argument("--b", required=True, metavar="LABEL=DEMUXDIR",
                    help="model B: label=demuxed dir (containing BarcodeNN.fastq)")
    ap.add_argument("--mapping", required=True,
                    help="plasmid<TAB>barcode<TAB>RUN号 TSV (barcode like '24标签-N')")
    ap.add_argument("--run", required=True, help="RUN号 to compare")
    ap.add_argument("--ref-dir", required=True,
                    help="dir of single-record STR<plasmid>.fa reference files")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rq-range", default=None,
                    help="optional rq filter, e.g. 0.99:1.1 (default: none, all demuxed reads)")
    ap.add_argument("--np-thr", type=int, default=5,
                    help="minimum np (passes) threshold; passed to sequencing_report_v2 as "
                         f"--np-range {{np-thr}}:10000000 (default: 5; use 0 for no filter)")
    args = ap.parse_args()

    la, da = args.a.split("=", 1)
    lb, db = args.b.split("=", 1)
    models = {la: da, lb: db}
    out = args.outdir
    os.makedirs(out, exist_ok=True)

    p2b = load_plasmid_barcode(args.mapping, args.run)
    log(f"run {args.run}: {len(p2b)} barcodes -> {sorted(p2b)}")
    log(f"A  {la} = {da}")
    log(f"B  {lb} = {db}")
    log(f"ref-dir: {args.ref_dir}")
    log(f"rq-range: {args.rq_range or 'NONE (all reads)'}")
    log(f"np-thr: {args.np_thr} (np >= {args.np_thr})")

    report_extra = []
    if args.rq_range:
        report_extra += ["--rq-range", args.rq_range]
    # np-thr is an "np >= threshold" filter -> sequencing_report_v2 --np-range thr:10000000
    if args.np_thr > 0:
        report_extra += ["--np-range", f"{args.np_thr}:10000000"]

    # ---- per-barcode reports for each model ----
    for model, demux in models.items():
        mdir = os.path.join(out, model)
        os.makedirs(mdir, exist_ok=True)
        for n in sorted(p2b):
            pl = p2b[n]
            fq = os.path.join(demux, f"Barcode{n:02d}.fastq")
            ref = os.path.join(args.ref_dir, f"STR{pl}.fa")
            if not os.path.exists(fq):
                log(f"[{model}] SKIP Barcode{n:02d} ({pl}): missing {fq}")
                continue
            if not os.path.exists(ref):
                log(f"[{model}] SKIP Barcode{n:02d} ({pl}): missing ref {ref}")
                continue
            stem = f"Barcode{n:02d}_{pl}"
            bam = os.path.join(mdir, f"{stem}.bam")
            if run([GS, FASTQ2BAM, fq, bam, "--ref", ref]) != 0:
                log(f"[{model}] FAIL fastq2bam {stem}")
                continue
            cmd = [GS, REPORT, "--bams", bam, "--refs", ref] + report_extra
            if run(cmd) != 0:
                log(f"[{model}] FAIL report {stem}")
                continue
            log(f"[{model}] ok {stem}")

    # ---- aggregate per model via the report's own merge logic ----
    agg_paths = {}
    for model in models:
        mdir = os.path.join(out, model)
        aggdir = os.path.join(mdir, "_agg")
        os.makedirs(aggdir, exist_ok=True)
        fact, basic = [], []
        for d in sorted(os.listdir(mdir)):
            md = os.path.join(mdir, d)
            if not os.path.isdir(md) or d == "_agg":
                continue
            for f in os.listdir(md):
                if f.endswith(".gsmm2_aligned_metric_fact.csv"):
                    dst = os.path.join(aggdir, f"{model}_{d}_{f}")
                    shutil.copy(os.path.join(md, f), dst)
                    fact.append(dst)
                elif f.endswith(".basic.csv"):
                    dst = os.path.join(aggdir, f"{model}_{d}_{f}")
                    shutil.copy(os.path.join(md, f), dst)
                    basic.append(dst)
        if not fact:
            log(f"[{model}] no fact csvs to aggregate")
            continue
        agg_out = os.path.join(aggdir, "aggregate.csv")
        if run([GS, REPORT, "--fact-csvs", *sorted(fact),
                "--basic-csvs", *sorted(basic)]) != 0:
            log(f"[{model}] FAIL aggregate")
            continue
        # merge_partition writes <first-fact-stem>.all.csv next to the fact csvs
        produced = [f for f in os.listdir(aggdir) if f.endswith(".all.csv")]
        if produced:
            shutil.move(os.path.join(aggdir, sorted(produced)[0]), agg_out)
            shutil.copy(agg_out, os.path.join(out, f"aggregate_{model}.csv"))
            agg_paths[model] = agg_out
            log(f"[{model}] aggregate -> {agg_out}")

    # ---- per-barcode comparison (long format) ----
    rows = []
    for n in sorted(p2b):
        pl = p2b[n]
        stem = f"Barcode{n:02d}_{pl}"
        d = {}
        for model in (la, lb):
            c = os.path.join(out, model, f"{stem}-metric",
                             f"{stem}.gsmm2_aligned_metric_aggr.csv")
            d[model] = read_kv(c) if os.path.exists(c) else {}
        for k in SEL:
            a = fval(d.get(la, {}), k)
            b = fval(d.get(lb, {}), k)
            delta = (a - b) if (a is not None and b is not None) else None
            rows.append((n, pl, k, fmt(a), fmt(b), fmt(delta, True)))
    with open(os.path.join(out, "per_barcode.tsv"), "w") as f:
        f.write(f"barcode\tplasmid\tmetric\t{la}\t{lb}\tdelta_{la}_minus_{lb}\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    log(f"wrote per_barcode.tsv ({len(rows)} rows)")

    # ---- aggregate comparison ----
    aggr = {m: read_kv(agg_paths[m]) for m in agg_paths}
    with open(os.path.join(out, "aggregate.csv"), "w") as f:
        f.write(f"metric\t{la}\t{lb}\tdelta_{la}_minus_{lb}\n")
        for k in SEL:
            a = fval(aggr.get(la, {}), k)
            b = fval(aggr.get(lb, {}), k)
            delta = (a - b) if (a is not None and b is not None) else None
            f.write("\t".join([k, fmt(a), fmt(b), fmt(delta, True)]) + "\n")
    log("wrote aggregate.csv")

    # ---- markdown summary ----
    md = [f"# {args.run} {la} vs {lb} — sequencing_report_v2\n"]
    md.append(f"rq-range: {args.rq_range or 'none (all demuxed reads)'}  \n")
    md.append(f"np-thr: >= {args.np_thr}  \n")
    md.append(f"barcodes: {sorted(p2b)}  \n")
    md.append("## Aggregate (pooled across barcodes)\n")
    md.append(f"| metric | {la} | {lb} | Δ({la}−{lb}) |")
    md.append("|---|---:|---:|---:|")
    for k in SEL:
        a = fval(aggr.get(la, {}), k)
        b = fval(aggr.get(lb, {}), k)
        delta = (a - b) if (a is not None and b is not None) else None
        md.append(f"| {k} | {fmt(a)} | {fmt(b)} | {fmt(delta, True)} |")
    md.append(f"\nSee `per_barcode.tsv` for the per-barcode breakdown.\n")
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(md))
    log("wrote summary.md")
    log("ALL DONE")


if __name__ == "__main__":
    main()
