#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Homopolymer first-base A/B comparison on a locus_error_rate joined_all table.

Given the ``joined_all.csv`` (tab-separated) produced by the
``smc_ab_test_with_barcode_ref_locus`` skill, this script:

  1. loads the per-plasmid Sanger reference sequences,
  2. finds every C/G homopolymer run of length >= ``--min-run`` in each
     reference (the ``pos`` in the gsetl table is 0-based),
  3. for each run keeps the requested ``--which`` base (default: the FIRST
     base of the run),
  4. pulls the matching (barcode, pos) rows from the joined table,
  5. reports the A/B comparison for both accuracy metrics —
        locus_accuracy        = eq / (eq+diff+ins+del)
        locus_accuracy_by_depth = eq / depth
     — as a per-locus table (with ``aroundBases`` context) plus a
     depth-weighted pooled number for each side.

Before printing, the script SELF-VERIFIES two invariants on the joined table
and aborts with a nonzero exit code if either fails, so a mis-indexed or
malformed input is caught rather than silently producing a shifted table:

  V1. the base bracketed ``[..]`` in ``aroundBases_exp`` equals the reference
      base at the 0-based ``pos`` column, for every row;
  V2. ``locus_accuracy_by_depth`` equals ``eq / depth`` for every row, and
      ``eq <= depth`` everywhere.

Outputs are written to ``--outdir``:
  - homopolymer_firstbase_locus.tsv   (selected rows + added base/runlen cols)
  - homopolymer_firstbase_summary.tsv (pooled metrics, one row per metric)
  - homopolymer_firstbase_report.md   (human-readable report)

Usage:
  python homopolymer_firstbase_ab.py \
      --joined    .../locus_error_rate/joined_all.csv \
      --ref-dir   .../merged_output \
      --mapping   .../plasmid_2_barcode.tsv \
      --run       20260805_250804Y0004_Run0001 \
      --outdir    .../locus_error_rate/homopolymer_firstbase \
      --min-run 4 --which first --char CG
"""

import argparse
import csv
import os
import sys


def read_mapping_run(map_path, run):
    """Return {barcode: plasmid} for the rows of ``run``."""
    bc2pl = {}
    with open(map_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("RUN号") == run:
                label = r["barcode"].split("-")[-1]
                bc2pl[f"Barcode{int(label):02d}"] = r["plasmid"]
    return bc2pl


def load_ref(ref_dir, plasmid):
    """Single-record FASTA -> uppercase sequence string."""
    seq = ""
    with open(os.path.join(ref_dir, f"STR{plasmid}.fa")) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq += line.upper()
    return seq


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
    if "[" in around and "]" in around:
        return around[around.index("[") + 1:around.index("]")]
    return None


def verify(joined_rows, seqs, bc2pl):
    """Run invariants V1 and V2. Return (n_rows, n_checked, error_message_or_None)."""
    n = 0
    for r in joined_rows:
        pl = bc2pl.get(r["barcode"])
        if pl is None or pl not in seqs:
            continue
        n += 1
        try:
            pos = int(r["pos"])
        except (TypeError, ValueError):
            return n, 0, f"row with non-integer pos: {r['pos']!r}"
        # V1: bracketed char == ref[0-based pos]
        b = bracketed(r.get("aroundBases_exp", ""))
        refbase = seqs[pl][pos] if 0 <= pos < len(seqs[pl]) else None
        if b is not None and refbase is not None and b != refbase:
            return n, 0, (f"V1 FAIL: barcode={r['barcode']} pos={pos} "
                          f"bracketed={b!r} ref[0based]={refbase!r} — pos column may "
                          f"not be 0-based, or aroundBases is mis-centered")
        # V2: by_depth == eq/depth and eq <= depth
        for side in ("ctrl", "exp"):
            eq, dep = r.get(f"eq_{side}"), r.get(f"depth_{side}")
            byd = r.get(f"locus_accuracy_by_depth_{side}")
            try:
                eq_i, dep_i = int(eq), int(dep)
            except (TypeError, ValueError):
                continue
            if eq_i > dep_i:
                return n, 0, (f"V2 FAIL: eq>depth for {side} in "
                              f"barcode={r['barcode']} pos={pos} (eq={eq_i} depth={dep_i})")
            if dep_i > 0:
                byd_f = num(byd)
                if byd_f is not None and abs(byd_f - eq_i / dep_i) > 1e-6:
                    return n, 0, (f"V2 FAIL: by_depth={byd} != eq/depth={eq_i/dep_i:.6f} "
                                  f"for {side} in barcode={r['barcode']} pos={pos}")
    return n, n, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--joined", required=True,
                    help="tab-separated joined_all.csv from the ref-locus skill")
    ap.add_argument("--ref-dir", required=True,
                    help="dir of single-record STR<plasmid>.fa references")
    ap.add_argument("--mapping", required=True,
                    help="plasmid<tab>barcode<tab>RUN号 TSV with header")
    ap.add_argument("--run", required=True, help="which RUN号 rows to use")
    ap.add_argument("--outdir", required=True, help="where outputs are written")
    ap.add_argument("--min-run", type=int, default=4,
                    help="minimum homopolymer run length to consider (default 4)")
    ap.add_argument("--which", default="first", choices=["first", "last", "all"],
                    help="which base of each run to compare (default first)")
    ap.add_argument("--char", default="CG",
                    help="chars that count as the homopolymer (default CG)")
    args = ap.parse_args()

    chars = set(args.char.upper())
    bc2pl = read_mapping_run(args.mapping, args.run)
    if not bc2pl:
        sys.exit(f"no mapping rows for --run {args.run!r}")

    seqs = {pl: load_ref(args.ref_dir, pl) for pl in set(bc2pl.values())}

    with open(args.joined, newline="") as f:
        joined_rows = list(csv.DictReader(f, delimiter="\t"))
    if not joined_rows:
        sys.exit(f"no rows in {args.joined}")
    header = list(joined_rows[0].keys())

    # ---- self-verification (abort on any failure) ----
    n_rows, n_checked, err = verify(joined_rows, seqs, bc2pl)
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
    for r in joined_rows:
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
        print(f"[warn] {len(missing)} target (barcode,pos) had no row in joined table "
              f"(e.g. 0-depth position): {sorted(missing)[:10]}")

    if not sel:
        sys.exit("no selected positions — check --min-run / --char / --which")

    # ---- metrics ----
    metrics = []
    for label, af in (("locus_accuracy (eq/(eq+diff+ins+del))", "locus_accuracy"),
                      ("locus_accuracy_by_depth (eq/depth)", "locus_accuracy_by_depth")):
        if af == "locus_accuracy":
            # pooled depth-weighted = Σeq / Σ(eq+diff+ins+del)
            ec = sum(int(r["eq_ctrl"]) for r in sel)
            et = sum(int(r["eq_ctrl"]) + int(r["diff_ctrl"]) + int(r["ins_ctrl"])
                     + int(r["del_ctrl"]) for r in sel)
            ee = sum(int(r["eq_exp"]) for r in sel)
            ete = sum(int(r["eq_exp"]) + int(r["diff_exp"]) + int(r["ins_exp"])
                      + int(r["del_exp"]) for r in sel)
            vc = ec / et if et else float("nan")
            ve = ee / ete if ete else float("nan")
        else:
            # pooled = Σeq / Σdepth
            ec = sum(int(r["eq_ctrl"]) for r in sel)
            ed = sum(int(r["depth_ctrl"]) for r in sel)
            ee = sum(int(r["eq_exp"]) for r in sel)
            ede = sum(int(r["depth_exp"]) for r in sel)
            vc = ec / ed if ed else float("nan")
            ve = ee / ede if ede else float("nan")
        metrics.append((label, vc, ve, ve - vc if (ve == ve and vc == vc) else float("nan")))

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
        w.writerow(["metric", "n_positions", "control", "experiment", "delta_exp_minus_ctrl"])
        for label, vc, ve, d in metrics:
            w.writerow([label, len(sel), f"{vc:.6f}", f"{ve:.6f}", f"{d:+.6f}"])

    # ---- write markdown report ----
    rep_path = os.path.join(args.outdir, "homopolymer_firstbase_report.md")
    lines = []
    lines.append(f"# Homopolymer first-base A/B (char={args.char.upper()}, "
                 f"run>={args.min_run}, which={args.which})")
    lines.append("")
    lines.append(f"- joined: `{args.joined}`")
    lines.append(f"- n target positions: {len(sel)}  (verified input: {n_checked}/{n_rows} rows)")
    lines.append("")
    lines.append("## Pooled")
    lines.append("")
    lines.append("| metric | control | experiment | Δ (exp−ctrl) |")
    lines.append("|---|---:|---:|---:|")
    for label, vc, ve, d in metrics:
        lines.append(f"| {label} | {vc:.6f} | {ve:.6f} | {d:+.6f} |")
    lines.append("")
    lines.append("## Per-position")
    lines.append("")
    lines.append("| barcode | plasmid | pos(0-based) | run | aroundBases([ ]=target) "
                 "| ctrl_acc | exp_acc | Δ | ctrl_bydepth | exp_bydepth |")
    lines.append("|---|---|---:|---:|---|---:|---:|---:|---:|---:|")
    for r in sel:
        pl = bc2pl[r["barcode"]]
        pos = int(r["pos"])
        ca = num(r["locus_accuracy_ctrl"]); ea = num(r["locus_accuracy_exp"])
        cb = num(r["locus_accuracy_by_depth_ctrl"]); eb = num(r["locus_accuracy_by_depth_exp"])
        d = (ea - ca) if None not in (ca, ea) else float("nan")
        lines.append(
            f"| {r['barcode']} | {pl} | {pos} | {r['runlen']} | {r['aroundBases_exp']} "
            f"| {ca:.5f} | {ea:.5f} | {d:+.5f} | {cb:.5f} | {eb:.5f} |")
    lines.append("")
    with open(rep_path, "w") as f:
        f.write("\n".join(lines))

    print(f"[write] {tsv_path} ({len(sel)} rows)")
    print(f"[write] {sum_path}")
    print(f"[write] {rep_path}")
    print("\n=== pooled ===")
    for label, vc, ve, d in metrics:
        print(f"  {label}:  ctrl={vc:.6f}  exp={ve:.6f}  Δ={d:+.6f}")


if __name__ == "__main__":
    main()
