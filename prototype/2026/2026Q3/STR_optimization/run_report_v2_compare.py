#!/usr/bin/env python3
"""Compare the two Run0001 barcode-call models with sequencing_report_v2.

For every barcode that the Run0001 mapping covers, this:
  1. converts the model's BarcodeNN.fastq to an unaligned BAM (rq tag) via fastq2bam.py,
  2. runs gseda.ppl.sequencing_report_v2 against the plasmid's Sanger reference,
  3. collects the per-barcode aggr metrics,
  4. aggregates all barcodes per model with the report's own
     --fact-csvs/--basic-csvs merge (merge_partition) logic.

Outputs (under --outdir):
  new_model/  baseline/   per-barcode BAMs + *-metric/ dirs
  per_barcode.tsv    long-format: sample x metric  (new_model, baseline, delta)
  aggregate.csv      pooled-metric x (new_model, baseline)  for each model
  aggregate_new_model.csv / aggregate_baseline.csv
  summary.md         human-readable comparison
"""
import argparse
import csv
import os
import pathlib
import shutil
import subprocess
import sys

GS = "/root/miniconda3/envs/gseda/bin/python"
REPO = "/root/projects/gsda/third_party/gseda/src/gseda"
FASTQ2BAM = f"{REPO}/file_format_cvt/fastq2bam.py"
REPORT = f"{REPO}/ppl/sequencing_report_v2.py"

BASE = "/data1/ccs_data/str-optimization/second-batch-of-data"
RUN = "20260805_250804Y0004_Run0001"
DEMUX = {
    "new_model": f"{BASE}/{RUN}/{RUN}_called-barcode-v4-2026Q2Model/demuxed",
    "baseline": f"{BASE}/{RUN}/{RUN}_called-barcode-v4-baseline/demuxed",
}
REFD = f"{BASE}/STR第二批一代测序/STR第二批一代测序/merged_output"
MAP = f"{BASE}/plasmid_2_barcode.tsv"

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


def aggr_csv_for(metric_dir, stem):
    return os.path.join(metric_dir, f"{stem}.gsmm2_aligned_metric_aggr.csv")


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rq-range", default=None,
                    help="optional rq filter, e.g. 0.99:1.1 (default: none, all demuxed reads)")
    args = ap.parse_args()
    out = args.outdir
    os.makedirs(out, exist_ok=True)

    p2b = load_plasmid_barcode(MAP, RUN)
    log(f"Run0001 mapping: {len(p2b)} barcodes -> {sorted(p2b)}")
    log(f"rq-range: {args.rq_range or 'NONE (all reads)'}")

    report_extra = ["--rq-range", args.rq_range] if args.rq_range else []

    # ---- per-barcode reports for each model ----
    for model, demux in DEMUX.items():
        mdir = os.path.join(out, model)
        os.makedirs(mdir, exist_ok=True)
        for n in sorted(p2b):
            pl = p2b[n]
            fq = os.path.join(demux, f"Barcode{n:02d}.fastq")
            ref = os.path.join(REFD, f"STR{pl}.fa")
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
    for model in DEMUX:
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
        for model in DEMUX:
            c = aggr_csv_for(os.path.join(out, model, f"{stem}-metric"), stem)
            d[model] = read_kv(c) if os.path.exists(c) else {}
        for k in SEL:
            v_new = fval(d.get("new_model", {}), k)
            v_base = fval(d.get("baseline", {}), k)
            delta = (v_new - v_base) if (v_new is not None and v_base is not None) else None
            rows.append((n, pl, k,
                         f"{v_new:.6g}" if v_new is not None else "",
                         f"{v_base:.6g}" if v_base is not None else "",
                         f"{delta:+.6g}" if delta is not None else ""))
    with open(os.path.join(out, "per_barcode.tsv"), "w") as f:
        f.write("barcode\tplasmid\tmetric\tnew_model\tbaseline\tdelta_new_minus_base\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    log(f"wrote per_barcode.tsv ({len(rows)} rows)")

    # ---- aggregate comparison ----
    aggr = {m: read_kv(agg_paths[m]) for m in agg_paths}
    with open(os.path.join(out, "aggregate.csv"), "w") as f:
        f.write("metric\tnew_model\tbaseline\tdelta_new_minus_base\n")
        for k in SEL:
            a = fval(aggr.get("new_model", {}), k)
            b = fval(aggr.get("baseline", {}), k)
            delta = (a - b) if (a is not None and b is not None) else None
            f.write("\t".join([k,
                               f"{a:.6g}" if a is not None else "",
                               f"{b:.6g}" if b is not None else "",
                               f"{delta:+.6g}" if delta is not None else ""]) + "\n")
    log(f"wrote aggregate.csv")

    # ---- markdown summary ----
    md = []
    md.append(f"# Run0001 new_model vs baseline — sequencing_report_v2\n")
    md.append(f"rq-range: {args.rq_range or 'none (all demuxed reads)'}  \n")
    md.append(f"barcodes: {sorted(p2b)}  \n")
    md.append("## Aggregate (pooled across barcodes)\n")
    md.append("| metric | new_model | baseline | Δ(new−base) |")
    md.append("|---|---:|---:|---:|")
    for k in SEL:
        a = fval(aggr.get("new_model", {}), k)
        b = fval(aggr.get("baseline", {}), k)
        delta = (a - b) if (a is not None and b is not None) else None
        md.append(f"| {k} | {a:.6g} | {b:.6g} | {delta:+.6g} |".replace(".nan", ""))
    md.append("\nSee `per_barcode.tsv` for the per-barcode breakdown.\n")
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(md))
    log("wrote summary.md")
    log("ALL DONE")


if __name__ == "__main__":
    main()
