#!/usr/bin/env python3
"""Per-barcode A/B of per-reference-locus accuracy, via locus_error_rate.

Given two per-barcode demux directories (e.g. a new barcode-calling model vs
the baseline over the same SMC run), the plasmid->barcode mapping, the
per-plasmid Sanger references, and a target RUN号, this:
  1. for every barcode the mapping covers, calls gseda.ab_analysis.
     locus_error_rate.main_cli with:
         control = control/BarcodeNN.fastq
         exp     = exp/BarcodeNN.fastq
         ref     = STR<plasmid>.fa
     producing per-barcode per-locus eq/diff/ins/del/depth tables plus
     per-locus locus_accuracy = eq/(eq+diff+ins+del).
  2. aggregates all barcodes per group (control_all / experiment_all) and
     the joined table (joined_all), each with a leading `barcode` column.
  3. writes a depth-weighted accuracy summary (overall + per-barcode).

IMPORTANT: the query files here are FASTQ, which carry NO CCS read-accuracy
(`rq`) or number-of-passes (`np`) field, so gsmm2's --rq-range/--np-range are
ignored. The metric is therefore alignment identity vs the Sanger reference,
NOT raw read accuracy.

Outputs (under --outdir):
  BarcodeNN/                    per-barcode control_locus_accuracy.csv,
                                experiment_locus_accuracy.csv,
                                locus_accuracy_joined.csv + aligned bams
  control_all.csv               all barcodes, control side, +barcode col
  experiment_all.csv            all barcodes, experiment side, +barcode col
  joined_all.csv                all barcodes, joined, +barcode col
  summary.md                    depth-weighted A/B summary (overall + per-barcode)
"""
import argparse
import os
import sys
from multiprocessing import cpu_count

GS = "/root/miniconda3/envs/gseda/bin/python"
REPO_SRC = "/root/projects/gsda/third_party/gseda/src"
if REPO_SRC not in sys.path:
    sys.path.insert(0, REPO_SRC)

import polars as pl  # noqa: E402
from gseda.ab_analysis.locus_error_rate import main_cli as ler_main  # noqa: E402
from gseda.fact_table_ana.polars_init import polars_env_init  # noqa: E402


def log(*a):
    print(*a, flush=True)


def parse_label_dir(s, name):
    """'LABEL=DIR' -> (label, dir)."""
    if "=" not in s:
        raise SystemExit(f"--{name} must be LABEL=DIR, got: {s}")
    label, d = s.split("=", 1)
    return label, d


def load_pairs(mapping, run):
    """[(plasmid, barcode_int)] for the given RUN号, sorted by barcode int.

    mapping columns: plasmid<TAB>barcode<TAB>RUN号 (barcode like '24标签-N').
    """
    rows = []
    with open(mapping) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("plasmid\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            plasmid, label, r = (x.strip() for x in parts[:3])
            if not plasmid or r != run:
                continue
            try:
                n = int(label.rsplit("-", 1)[1])
            except (IndexError, ValueError):
                log(f"WARN cannot parse barcode '{label}'")
                continue
            rows.append((plasmid, n))
    rows.sort(key=lambda x: x[1])
    return rows


def depth_weighted(sub, tag):
    """Depth-weighted locus accuracy eq/(eq+diff+ins+del) for one side."""
    t = sub["eq_" + tag] + sub["diff_" + tag] + sub["ins_" + tag] + sub["del_" + tag]
    tot = t.sum()
    return sub["eq_" + tag].sum() / tot if tot else None


def main():
    ap = argparse.ArgumentParser(
        description="Per-barcode A/B of per-reference-locus accuracy via "
                    "gseda.ab_analysis.locus_error_rate")
    ap.add_argument("--control", required=True, metavar="LABEL=DEMUXDIR",
                    help="control: label=demuxed dir (BarcodeNN.fastq)")
    ap.add_argument("--exp", required=True, metavar="LABEL=DEMUXDIR",
                    help="experiment: label=demuxed dir (BarcodeNN.fastq)")
    ap.add_argument("--mapping", required=True,
                    help="plasmid<TAB>barcode<TAB>RUN号 TSV (barcode like '24标签-N')")
    ap.add_argument("--run", required=True, help="RUN号 to use (only its rows of the mapping)")
    ap.add_argument("--ref-dir", required=True,
                    help="dir of single-record STR<plasmid>.fa reference files")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rq-thr", default="0",
                    help="read-accuracy thr for locus_error_rate (1 or 2 vals, "
                         "'[ctrl,exp]' or shared). FASTQ has no rq field so this is "
                         "ignored by gsmm2; must still be supplied (required arg). "
                         "Default 0 (non-filtering).")
    ap.add_argument("--np-thr", default="5",
                    help="number-of-passes lower bound (1 or 2 vals). Ignored for "
                         "FASTQ (no np field). Default 5.")
    ap.add_argument("--threads", type=int, default=None,
                    help="gsmm2 thread count (default: CPU count)")
    args = ap.parse_args()

    polars_env_init()

    lc, dctrl = parse_label_dir(args.control, "control")
    lx, dexp = parse_label_dir(args.exp, "exp")
    out = args.outdir
    os.makedirs(out, exist_ok=True)

    pairs = load_pairs(args.mapping, args.run)
    log(f"run {args.run}: {len(pairs)} barcodes -> {[(p, f'{n:02d}') for p, n in pairs]}")
    log(f"control     {lc} = {dctrl}")
    log(f"experiment  {lx} = {dexp}")
    log(f"ref-dir: {args.ref_dir}")
    log(f"rq-thr: {args.rq_thr}   np-thr: {args.np_thr}  (ignored: FASTQ has no rq/np field)")

    # ---- per-barcode locus_error_rate runs ----
    done = []  # (plasmid, n) actually processed
    for plasmid, n in pairs:
        bc = f"Barcode{n:02d}"
        ctrl_fq = os.path.join(dctrl, f"{bc}.fastq")
        exp_fq = os.path.join(dexp, f"{bc}.fastq")
        ref_fa = os.path.join(args.ref_dir, f"STR{plasmid}.fa")
        if not os.path.exists(ctrl_fq):
            log(f"[{bc}] SKIP ({plasmid}): missing control {ctrl_fq}")
            continue
        if not os.path.exists(exp_fq):
            log(f"[{bc}] SKIP ({plasmid}): missing exp {exp_fq}")
            continue
        if not os.path.exists(ref_fa):
            log(f"[{bc}] SKIP ({plasmid}): missing ref {ref_fa}")
            continue
        log(f"\n==== {bc} -> {plasmid} ====")
        ler_main([
            "--control", ctrl_fq,
            "--exp", exp_fq,
            "--ref", ref_fa,
            "--rq-thr", args.rq_thr,
            "--np-thr", args.np_thr,
            "--outdir", os.path.join(out, bc),
        ] + (["--threads", str(args.threads)] if args.threads else []))
        done.append((plasmid, n))

    if not done:
        log("ERROR: no barcodes processed")
        return

    done_ns = sorted(n for _, n in done)

    # ---- aggregate all barcodes per group ----
    def collect(name):
        frames = []
        for n in done_ns:
            p = os.path.join(out, f"Barcode{n:02d}", name)
            if not os.path.exists(p):
                continue
            df = pl.read_csv(p, separator="\t")
            frames.append(df.with_columns(pl.lit(f"Barcode{n:02d}").alias("barcode")))
        outdf = pl.concat(frames)
        return outdf.select(["barcode"] + [c for c in outdf.columns if c != "barcode"])

    for name, tag in [
        ("control_locus_accuracy.csv", "control_all"),
        ("experiment_locus_accuracy.csv", "experiment_all"),
        ("locus_accuracy_joined.csv", "joined_all"),
    ]:
        df = collect(name)
        dest = os.path.join(out, f"{tag}.csv")
        df.write_csv(dest, separator="\t")
        log(f"aggregate {tag}: {df.height} rows, {df.width} cols -> {dest}")

    # ---- depth-weighted accuracy summary ----
    jpath = os.path.join(out, "joined_all.csv")
    if not os.path.exists(jpath):
        log("WARN joined_all.csv missing; skipping summary")
        return
    j = pl.read_csv(jpath, separator="\t")
    c_all = depth_weighted(j, "ctrl")
    e_all = depth_weighted(j, "exp")

    md = [f"# {args.run} {lc} vs {lx} — per-locus accuracy A/B\n"]
    md.append("depth-weighted locus_accuracy = eq/(eq+diff+ins+del) over the "
              "reference, per (barcode, refname, pos) row.\n")
    md.append(f"metric: alignment identity vs the Sanger reference (FASTQ input → "
              f"no rq/np filtering applies).\n")
    md.append("## Overall (all barcodes pooled)\n")
    md.append(f"| side | locus_accuracy |")
    md.append(f"|---|---:|")
    md.append(f"| {lc} (control) | {c_all:.6f} |")
    md.append(f"| {lx} (experiment) | {e_all:.6f} |")
    md.append(f"| Δ (exp − ctrl) | {e_all - c_all:+.6f} |\n")
    md.append("## Per-barcode\n")
    md.append(f"| barcode | locus | {lc} | {lx} | Δ (exp−ctrl) |")
    md.append("|---|---|---:|---:|---:|")
    for n in done_ns:
        bc = f"Barcode{n:02d}"
        sub = j.filter(pl.col("barcode") == bc)
        c = depth_weighted(sub, "ctrl")
        e = depth_weighted(sub, "exp")
        loci = sub["refname"][0]
        d = (e - c) if (c is not None and e is not None) else None
        md.append(f"| {bc} | {loci} | {c:.6f} | {e:.6f} | "
                  f"{d:+.6f} |")
    md.append("")
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(md))
    log(f"overall: {lc}={c_all:.5f}  {lx}={e_all:.5f}  Δ={e_all - c_all:+.5f}")
    log(f"wrote summary.md")
    log("ALL DONE")


if __name__ == "__main__":
    main()
