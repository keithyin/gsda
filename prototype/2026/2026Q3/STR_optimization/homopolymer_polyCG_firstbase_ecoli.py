#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Poly C/G homopolymer (run len 4 or 5) FIRST-base eq/depth A/B on a
single-reference locus_error_rate joined table (E. coli MG1655).

Inputs:
  - locus_accuracy_joined.csv (tab-sep, keyed on refname,pos; _ctrl/_exp cols)
  - MG1655.fa (single contig; gsetl `pos` is 0-based)

For every C or G homopolymer run of length 4 or 5 in the reference, keep the
FIRST base (run start), pull its (refname,pos) row from the joined table, and
report:
  1. per-locus table: refname, locus(pos), aroundBases, base, repeatCount, and
     each metric (eq, depth, eq/depth, diff, ins, del, locus_accuracy) for both
     ctrl and exp.
  2. per (base, repeatCount) group: 25/50/75 percentiles of eq/depth for ctrl
     and exp, plus delta.
"""
import csv
import os

import polars as pl

JOINED = "/data1/ccs_data/202603-good-ecoli-data/Run0002_adapter-v4-ecoliMG1655_locus_ab/locus_accuracy_joined.csv"
REF = "/data1/REF_GENOMES/MG1655.fa"
OUT = "/data1/ccs_data/202603-good-ecoli-data/Run0002_adapter-v4-ecoliMG1655_locus_ab/homopolymer_polyCG_4_5_firstbase"

LEN_SET = {4, 5}
CHARS = set("CG")


def load_ref(path):
    seq = ""
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            continue
        seq += line.upper()
    return seq


def runs_len(seq, chars, lens):
    """[(start_0based, run_len, base)] for runs of chars with len in lens."""
    out = []
    i, n = 0, len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        L = j - i
        if seq[i] in chars and L in lens:
            out.append((i, L, seq[i]))
        i = j
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    seq = load_ref(REF)
    targets = runs_len(seq, CHARS, LEN_SET)  # (start, len, base)
    print(f"reference contig len={len(seq)}; C/G runs len in {sorted(LEN_SET)}: {len(targets)}")

    j = pl.read_csv(JOINED, separator="\t")
    print(f"joined rows={j.height}; cols={len(j.columns)}")

    # Build a lookup pos -> target info
    tmap = {pos: (base, L) for pos, L, base in targets}
    tpos = pl.DataFrame({"pos": [p for p, _, _ in targets]})

    # join on refname,pos (single refname mg1655)
    sel = tpos.join(
        j.filter(pl.col("refname") == "mg1655"), on="pos", how="inner")
    # attach base / repeatCount
    sel = sel.with_columns([
        pl.Series("base", [tmap[p][0] for p in sel["pos"].to_list()]),
        pl.Series("repeatCount", [tmap[p][1] for p in sel["pos"].to_list()]),
    ])
    missing = len(targets) - sel.height
    print(f"selected {sel.height}/{len(targets)} first-base loci (missing={missing})")
    sel = sel.sort("pos")

    # ---- per-locus TSV ----
    per_cols = [
        "refname", "pos", "base", "repeatCount",
        "aroundBases_exp",
        "eq_ctrl", "depth_ctrl", "locus_accuracy_by_depth_ctrl",
        "diff_ctrl", "ins_ctrl", "del_ctrl", "locus_accuracy_ctrl",
        "eq_exp", "depth_exp", "locus_accuracy_by_depth_exp",
        "diff_exp", "ins_exp", "del_exp", "locus_accuracy_exp",
    ]
    per = sel.select(per_cols).rename({"pos": "locus"})
    per_path = os.path.join(OUT, "homopolymer_polyCG_4_5_firstbase_locus.tsv")
    per.write_csv(per_path, separator="\t")

    # ---- percentiles per (base, repeatCount) of eq/depth (locus_accuracy_by_depth) ----
    # eq/depth is undefined (NaN) where depth==0, and polars drop_nulls() keeps NaN.
    # Restrict to loci with finite eq/depth on that side (depth>0). ctrl and exp
    # carry identical per-position depth, so the depth-0 set is the same on both.
    pct_rows = []
    order = [("C", 4), ("C", 5), ("G", 4), ("G", 5)]
    for base, L in order:
        sub = sel.filter((pl.col("base") == base) & (pl.col("repeatCount") == L))
        for col, tag in (("locus_accuracy_by_depth_ctrl", "ctrl"),
                         ("locus_accuracy_by_depth_exp", "exp")):
            vals = [v for v in sub[col].to_list()
                    if v is not None and v == v]  # drop NaN (eq/depth undefined at depth 0)
            vals.sort()
            for q in (25, 50, 75):
                # linear-interpolation percentile (numpy default)
                if not vals:
                    v = float("nan")
                else:
                    idx = (len(vals) - 1) * q / 100.0
                    lo = int(idx); hi = min(lo + 1, len(vals) - 1)
                    frac = idx - lo
                    v = vals[lo] + (vals[hi] - vals[lo]) * frac
                pct_rows.append((base, L, tag, q, v, len(vals)))

    # pivot into a comparison frame
    pctl_cols = ["base", "repeatCount",
                 "eq_over_depth_ctrl_p25", "eq_over_depth_ctrl_p50", "eq_over_depth_ctrl_p75",
                 "eq_over_depth_exp_p25", "eq_over_depth_exp_p50", "eq_over_depth_exp_p75",
                 "delta_p25", "delta_p50", "delta_p75", "n"]
    grid = {(b, L, t, q): v for (b, L, t, q, v, n) in pct_rows}
    ngrid = {(b, L, t): n for (b, L, t, q, v, n) in pct_rows}
    pctl = pl.DataFrame([
        [b, L,
         grid[(b, L, "ctrl", 25)], grid[(b, L, "ctrl", 50)], grid[(b, L, "ctrl", 75)],
         grid[(b, L, "exp", 25)], grid[(b, L, "exp", 50)], grid[(b, L, "exp", 75)],
         grid[(b, L, "exp", 25)] - grid[(b, L, "ctrl", 25)],
         grid[(b, L, "exp", 50)] - grid[(b, L, "ctrl", 50)],
         grid[(b, L, "exp", 75)] - grid[(b, L, "ctrl", 75)],
         ngrid[(b, L, "ctrl")]]
        for (b, L) in order
    ], schema=pctl_cols, orient="row")
    pctl_path = os.path.join(OUT, "homopolymer_polyCG_4_5_firstbase_percentiles.tsv")
    pctl.write_csv(pctl_path, separator="\t")

    # ---- markdown report ----
    md = []
    md.append("# Poly C/G homopolymer (run len 4 / 5) — first-base eq/depth A/B")
    md.append("")
    md.append(f"- joined: `{JOINED}`")
    md.append(f"- reference: `{os.path.basename(REF)}` (contig len {len(seq)})")
    md.append(f"- target: C/G runs of length **4 or 5**, first (5'->3') base of each run")
    md.append(f"- loci selected: **{sel.height}** of {len(targets)} (missing {missing})")
    md.append("")
    md.append("Metric `locus_accuracy_by_depth` = eq/depth per reference base position.")
    md.append("Δ = exp − ctrl.")
    md.append("")
    md.append("## (2) Per (base, repeatCount) percentile comparison of eq/depth")
    md.append("")
    md.append("| group | n | ctrl p25 | ctrl p50 | ctrl p75 | exp p25 | exp p50 | exp p75 | Δ p25 | Δ p50 | Δ p75 |")
    md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in pctl.iter_rows():
        b, L, c25, c50, c75, e25, e50, e75, d25, d50, d75, n = r
        md.append(f"| {b}{L} | {n} | {c25:.4f} | {c50:.4f} | {c75:.4f} "
                  f"| {e25:.4f} | {e50:.4f} | {e75:.4f} | {d25:+.4f} | {d50:+.4f} | {d75:+.4f} |")
    md.append("")
    md.append("## (1) Per-locus first-base table")
    md.append("")
    md.append(f"Full table: `{per_path}` ({per.height} rows).")
    md.append("Columns: refname, locus, base, repeatCount, aroundBases([ ]=first base), "
              "eq/depth/eq÷depth/diff/ins/del/locus_accuracy for _ctrl then _exp.")
    md.append("")
    md.append("| refname | locus | base | rc | aroundBases | eq_c | d_c | eq/dep_c | diff_c | ins_c | del_c | eq_x | d_x | eq/dep_x | diff_x | ins_x | del_x |")
    md.append("|---|---:|:---:|:---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in per.iter_rows():
        (refname, locus, base, rc, around,
         eq_c, d_c, ed_c, diff_c, ins_c, del_c, acc_c,
         eq_x, d_x, ed_x, diff_x, ins_x, del_x, acc_x) = r
        def f(v):
            return "" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v))
        md.append(f"| {refname} | {locus} | {base} | {rc} | {around} "
                  f"| {eq_c} | {d_c} | {f(ed_c)} | {diff_c} | {ins_c} | {del_c} "
                  f"| {eq_x} | {d_x} | {f(ed_x)} | {diff_x} | {ins_x} | {del_x} |")
    md.append("")
    rep_path = os.path.join(OUT, "homopolymer_polyCG_4_5_firstbase_report.md")
    with open(rep_path, "w") as fp:
        fp.write("\n".join(md))

    print(f"[write] {per_path} ({per.height} rows)")
    print(f"[write] {pctl_path}")
    print(f"[write] {rep_path}")
    print("\n=== percentiles (eq/depth) ===")
    with open(pctl_path) as fp:
        print(fp.read())


if __name__ == "__main__":
    main()
