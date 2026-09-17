#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Homopolymer first-base accuracy (single side), two input modes.

This is the **non-A/B** counterpart of the ``smc_ab_homopolymer_firstbase``
skill. Instead of comparing a control side against an experiment side, it
reports a **single** side's accuracy at the homopolymer run boundaries.

Two ways to give it a single-side per-locus table:

MODE 1 — from a ready table (``--table``):
  Pass the tab-separated per-locus table the
  ``smc_ab_test_with_barcode_ref_locus`` skill writes (``control_all.csv`` /
  ``experiment_all.csv`` with a leading ``barcode`` column, or a per-barcode
  ``control_locus_accuracy.csv``) plus ``--ref-dir`` / ``--mapping`` / ``--run``.
  No gsmm2/gsetl is run; this is a pure table post-process.

MODE 2 — from a raw query (``--query`` + ``--ref``):
  Pass a query file (an **unmapped** BAM, or a FASTQ) and a single-record
  reference FASTA. Every read in the query is treated as a sequencing of that
  one reference. The script gsmm2-aligns the query to the reference and runs
  ``gsetl aligned-bam`` (via ``gseda.ab_analysis.locus_error_rate.process_group``)
  to build the per-locus eq/diff/ins/del/depth table in memory, then runs the
  same homopolymer analysis. ``--mapping`` / ``--run`` / ``--ref-dir`` are NOT
  needed in this mode (there is one reference, no barcodes).

Common to both modes, the script then:
  1. finds every C/G homopolymer run of length >= ``--min-run`` in the
     reference(s) (the ``pos`` in the gsetl table is 0-based),
  2. for each run keeps the requested ``--which`` base (default: the FIRST
     base of the run),
  3. pulls the matching (label, pos) rows,
  4. reports both accuracy metrics —
        locus_accuracy        = eq / (eq+diff+ins+del)
        locus_accuracy_by_depth = eq / depth
     — as a per-locus table (with ``aroundBases`` context) plus a
     depth-weighted pooled number.

Before printing, the script SELF-VERIFIES two invariants on the per-locus
table and aborts with a nonzero exit code if either fails:

  V1. the base bracketed ``[..]`` in ``aroundBases`` equals the reference base
      at the 0-based ``pos`` column, for every row;
  V2. ``locus_accuracy_by_depth`` equals ``eq / depth`` for every row, and
      ``eq <= depth`` everywhere.

Outputs are written to ``--outdir``:
  - homopolymer_firstbase_locus.tsv   (selected rows + added base/runlen cols)
  - homopolymer_firstbase_summary.tsv (pooled metrics, one row per metric)
  - homopolymer_firstbase_report.md   (human-readable report)
  (mode 2 also writes gsmm2/gsetl intermediates under ``--outdir/aligned/``)

Usage (mode 1):
  python homopolymer_firstbase.py \
      --table    .../locus_error_rate/control_all.csv \
      --ref-dir  .../merged_output \
      --mapping  .../plasmid_2_barcode.tsv \
      --run      20260805_250804Y0004_Run0001 \
      --outdir   .../locus_error_rate/homopolymer_firstbase \
      --min-run 4 --which first --char CG

Usage (mode 2):
  python homopolymer_firstbase.py \
      --query    .../reads.bam \
      --ref      .../STR3N-1.fa \
      --outdir   .../homopolymer_firstbase \
      --min-run 4 --which first --char CG
"""

import argparse
import csv
import os
import sys


def read_fasta(path):
    """Single-record FASTA -> uppercase sequence string."""
    seq = ""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq += line.upper()
    return seq


def load_ref(ref_dir, plasmid):
    """``STR<plasmid>.fa`` in ``ref_dir`` -> uppercase sequence string."""
    return read_fasta(os.path.join(ref_dir, f"STR{plasmid}.fa"))


def ref_stem(path):
    """Plasmid key from a reference filename: strip the ext and a leading
    ``STR``. ``STR3N-1.fa`` -> ``3N-1``."""
    stem = os.path.basename(path)
    for ext in (".fa", ".fasta", ".fna", ".faa"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    if stem.startswith("STR"):
        stem = stem[3:]
    return stem or os.path.basename(path)


def read_mapping_run(map_path, run):
    """Return {barcode: plasmid} for the rows of ``run``."""
    bc2pl = {}
    with open(map_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("RUN号") == run:
                label = r["barcode"].split("-")[-1]
                bc2pl[f"Barcode{int(label):02d}"] = r["plasmid"]
    return bc2pl


def build_from_query(args):
    """Mode 2: gsmm2-align query->ref + gsetl locus table, then package the
    rows for the shared analysis. Returns (rows, seqs, bc2pl)."""
    ref_path = args.ref
    seq = read_fasta(ref_path)
    if not seq:
        sys.exit(f"empty reference: {ref_path}")
    plasmid = ref_stem(ref_path)
    label = args.query_label or os.path.splitext(
        os.path.basename(args.query))[0]

    # Lazy import so MODE 1 (table) still needs only the stdlib.
    REPO_SRC = "/root/projects/gsda/third_party/gseda/src"
    if REPO_SRC not in sys.path:
        sys.path.insert(0, REPO_SRC)
    from gseda.fact_table_ana.polars_init import polars_env_init
    from gseda.ab_analysis.locus_error_rate import (
        process_group, parse_rq_thr, parse_np_thr,
        resolve_rq_range, resolve_np_range)
    polars_env_init()

    threads = args.threads or os.cpu_count()
    rq = parse_rq_thr(args.rq_thr)
    npth = parse_np_thr(args.np_thr) if args.np_thr else None
    rq_range = resolve_rq_range(rq, "query")
    np_range = resolve_np_range(npth, "query") if npth is not None else None
    print(f"[query] {args.query}  vs  {ref_path}\n"
          f"        rq-range={rq_range}  np-range={np_range}  threads={threads}")

    align_outdir = os.path.join(args.outdir, "aligned")
    os.makedirs(align_outdir, exist_ok=True)
    df = process_group(args.query, ref_path, rq_range, np_range, threads,
                       "query", align_outdir)
    if df.height == 0:
        sys.exit("query produced no loci (alignment produced nothing?)")

    rows = df.to_dicts()
    for r in rows:
        r["barcode"] = label  # single side, single reference -> one label
    return rows, {plasmid: seq}, {label: plasmid}


def runs_of(seq, chars, min_len):
    """Return list of (start_0based, run_len) for runs of a char in ``chars``
    whose length is >= ``min_len``."""
    out = []
    i, n = 0, len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len and seq[i] in chars:
            out.append((i, j - i))
        i = j
    return out


def which_pos(start, run_len, which):
    if which == "first":
        return start
    if which == "last":
        return start + run_len - 1
    if which == "all":
        return None  # handled by caller
    raise ValueError(f"unknown --which: {which}")


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def bracketed(around):
    """Return the character inside the [..] brackets of an aroundBases string."""
    if around and "[" in around and "]" in around:
        return around[around.index("[") + 1:around.index("]")]
    return None


def verify(rows, seqs, bc2pl):
    """Run invariants V1 and V2. Return (n_rows, n_checked, error_message_or_None)."""
    n = 0
    for r in rows:
        pl = bc2pl.get(r["barcode"])
        if pl is None or pl not in seqs:
            continue
        n += 1
        try:
            pos = int(r["pos"])
        except (TypeError, ValueError):
            return n, 0, f"row with non-integer pos: {r['pos']!r}"
        # V1: bracketed char == ref[0-based pos]
        b = bracketed(r.get("aroundBases", ""))
        refbase = seqs[pl][pos] if 0 <= pos < len(seqs[pl]) else None
        if b is not None and refbase is not None and b != refbase:
            return n, 0, (f"V1 FAIL: barcode={r['barcode']} pos={pos} "
                          f"bracketed={b!r} ref[0-based]={refbase!r} — pos column may "
                          f"not be 0-based, or aroundBases is mis-centered")
        # V2: by_depth == eq/depth and eq <= depth
        eq, dep, byd = r.get("eq"), r.get("depth"), r.get("locus_accuracy_by_depth")
        try:
            eq_i, dep_i = int(eq), int(dep)
        except (TypeError, ValueError):
            continue
        if eq_i > dep_i:
            return n, 0, (f"V2 FAIL: eq>depth in barcode={r['barcode']} "
                          f"pos={pos} (eq={eq_i} depth={dep_i})")
        if dep_i > 0:
            byd_f = num(byd)
            if byd_f is not None and abs(byd_f - eq_i / dep_i) > 1e-6:
                return n, 0, (f"V2 FAIL: by_depth={byd} != eq/depth={eq_i/dep_i:.6f} "
                              f"in barcode={r['barcode']} pos={pos}")
    return n, n, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # --- input source: EITHER --table (mode 1) OR --query + --ref (mode 2) ---
    ap.add_argument("--table", default=None,
                    help="[mode 1] tab-separated single-side per-locus table "
                         "(control_all.csv / experiment_all.csv or a per-barcode "
                         "locus_accuracy.csv)")
    ap.add_argument("--query", default=None,
                    help="[mode 2] query file: an UNMAPPED bam or a fastq. All its "
                         "reads are treated as sequencing the single --ref.")
    ap.add_argument("--ref", default=None,
                    help="[mode 2] single-record reference FASTA (all query reads "
                         "are sequenced against this one reference)")
    ap.add_argument("--query-label", default=None,
                    help="[mode 2] label to stamp on the query's rows (default: "
                         "the query filename without extension)")
    # --- mode 1 only ---
    ap.add_argument("--ref-dir", default=None,
                    help="[mode 1] dir of single-record STR<plasmid>.fa references")
    ap.add_argument("--mapping", default=None,
                    help="[mode 1] plasmid<tab>barcode<tab>RUN号 TSV with header")
    ap.add_argument("--run", default=None, help="[mode 1] which RUN号 rows to use")
    # --- mode 2 only ---
    ap.add_argument("--rq-thr", default="0",
                    help="[mode 2] read-accuracy lower bound forwarded to gsmm2 "
                         "(FASTQ has no rq field so it is ignored). Default 0 "
                         "(non-filtering).")
    ap.add_argument("--np-thr", default="0",
                    help="[mode 2] number-of-passes lower bound forwarded to gsmm2 "
                         "(FASTQ has no np field so it is ignored). Default 0 "
                         "(non-filtering).")
    ap.add_argument("--threads", type=int, default=None,
                    help="[mode 2] gsmm2 thread count (default: CPU count)")
    # --- both modes ---
    ap.add_argument("--outdir", required=True, help="where outputs are written")
    ap.add_argument("--min-run", type=int, default=4,
                    help="minimum homopolymer run length to consider (default 4)")
    ap.add_argument("--which", default="first", choices=["first", "last", "all"],
                    help="which base of each run to measure (default first)")
    ap.add_argument("--char", default="CG",
                    help="chars that count as the homopolymer (default CG)")
    args = ap.parse_args()

    chars = set(args.char.upper())

    # ---- build (rows, seqs, bc2pl) from the chosen input mode ----
    if args.query:
        if not args.ref:
            sys.exit("--query requires --ref")
        rows, seqs, bc2pl = build_from_query(args)
        if not rows:
            sys.exit(f"no rows from query {args.query}")
    else:
        for req in ("table", "ref_dir", "mapping", "run"):
            if getattr(args, req) is None:
                sys.exit(f"missing --{req.replace('_', '-')} (or use "
                         f"--query + --ref)")
        bc2pl = read_mapping_run(args.mapping, args.run)
        if not bc2pl:
            sys.exit(f"no mapping rows for --run {args.run!r}")
        seqs = {pl: load_ref(args.ref_dir, pl) for pl in set(bc2pl.values())}
        with open(args.table, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows:
            sys.exit(f"no rows in {args.table}")
    header = list(rows[0].keys())

    # ---- self-verification (abort on any failure) ----
    n_rows, n_checked, err = verify(rows, seqs, bc2pl)
    if err:
        sys.exit(f"INPUT VERIFICATION FAILED ({n_rows} rows scanned):\n  {err}")
    print(f"[verify] OK: {n_checked}/{n_rows} rows pass V1 (bracket==ref[0-based pos]) "
          f"and V2 (by_depth==eq/depth, eq<=depth)")

    # ---- select target (barcode, pos) pairs ----
    targets = {}  # (barcode, pos) -> (plasmid, run_len)
    for bc, pl in bc2pl.items():
        for start, rl in runs_of(seqs[pl], chars, args.min_run):
            if args.which == "all":
                for k in range(start, start + rl):
                    targets[(bc, k)] = (pl, rl)
            else:
                p = which_pos(start, rl, args.which)
                targets[(bc, p)] = (pl, rl)

    sel = []
    seen = set()
    for r in rows:
        key = (r["barcode"], int(r["pos"]))
        if key in targets and key not in seen:
            seen.add(key)
            d = dict(r)
            d["runlen"] = targets[key][1]
            d["base"] = seqs[targets[key][0]][int(r["pos"])]
            sel.append(d)
    missing = set(targets) - seen
    sel.sort(key=lambda r: (r["barcode"], int(r["pos"])))
    if missing:
        print(f"[warn] {len(missing)} target (barcode,pos) had no row in table "
              f"(e.g. 0-depth position): {sorted(missing)[:10]}")

    if not sel:
        sys.exit("no selected positions — check --min-run / --char / --which")

    # ---- metrics ----
    metrics = []
    # locus_accuracy = eq/(eq+diff+ins+del); pooled depth-weighted = Σeq/Σ(...)
    ec = sum(int(r["eq"]) for r in sel)
    et = sum(int(r["eq"]) + int(r["diff"]) + int(r["ins"]) + int(r["del"]) for r in sel)
    vc = ec / et if et else float("nan")
    metrics.append(("locus_accuracy (eq/(eq+diff+ins+del))", vc))
    # locus_accuracy_by_depth = eq/depth; pooled = Σeq/Σdepth
    ec2 = sum(int(r["eq"]) for r in sel)
    ed = sum(int(r["depth"]) for r in sel)
    vb = ec2 / ed if ed else float("nan")
    metrics.append(("locus_accuracy_by_depth (eq/depth)", vb))

    # ---- write TSV of selected rows ----
    os.makedirs(args.outdir, exist_ok=True)
    tsv_path = os.path.join(args.outdir, "homopolymer_firstbase_locus.tsv")
    with open(tsv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header + ["base", "runlen"], delimiter="\t")
        w.writeheader()
        w.writerows(sel)

    # ---- write summary TSV ----
    sum_path = os.path.join(args.outdir, "homopolymer_firstbase_summary.tsv")
    with open(sum_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "n_positions", "value"])
        for label, v in metrics:
            w.writerow([label, len(sel), f"{v:.6f}"])

    # ---- write markdown report ----
    rep_path = os.path.join(args.outdir, "homopolymer_firstbase_report.md")
    if args.query:
        src = f"query: `{args.query}`  vs  ref: `{args.ref}`"
    else:
        src = f"table: `{args.table}`"
    lines = []
    lines.append(f"# Homopolymer first-base accuracy (char={args.char.upper()}, "
                 f"run>={args.min_run}, which={args.which})")
    lines.append("")
    lines.append(f"- {src}")
    lines.append(f"- n target positions: {len(sel)}  (verified input: {n_checked}/{n_rows} rows)")
    lines.append("")
    lines.append("## Pooled")
    lines.append("")
    lines.append("| metric | value |")
    lines.append("|---|---:|")
    for label, v in metrics:
        lines.append(f"| {label} | {v:.6f} |")
    lines.append("")
    lines.append("## Per-position")
    lines.append("")
    lines.append("| barcode | plasmid | pos(0-based) | run | aroundBases([ ]=target) "
                 "| locus_accuracy | locus_accuracy_by_depth |")
    lines.append("|---|---|---:|---:|---|---:|---:|")
    for r in sel:
        pl = bc2pl[r["barcode"]]
        pos = int(r["pos"])
        a = num(r["locus_accuracy"])
        b = num(r["locus_accuracy_by_depth"])
        lines.append(
            f"| {r['barcode']} | {pl} | {pos} | {r['runlen']} | {r['aroundBases']} "
            f"| {a:.5f} | {b:.5f} |")
    lines.append("")
    with open(rep_path, "w") as f:
        f.write("\n".join(lines))

    print(f"[write] {tsv_path} ({len(sel)} rows)")
    print(f"[write] {sum_path}")
    print(f"[write] {rep_path}")
    print("\n=== pooled ===")
    for label, v in metrics:
        print(f"  {label}:  {v:.6f}")


if __name__ == "__main__":
    main()
