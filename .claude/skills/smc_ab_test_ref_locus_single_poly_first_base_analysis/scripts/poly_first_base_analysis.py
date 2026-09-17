#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Poly-homopolymer first-base A/B on a SINGLE-reference locus_error_rate joined table.

Post-processes the ``locus_accuracy_joined.csv`` (tab-separated, keyed on
``refname,pos``; columns ``…_ctrl`` / ``…_exp``) produced by the sibling skill
``smc_ab_test_ref_locus_single``. That skill gives one per-reference-locus table
per side (eq/diff/ins/del/depth + accuracy at every reference base position);
this one answers the narrower question:

  > At the homopolymer boundaries (default: the FIRST base of every polyC/G run
    of length 4 or 5), how do the two calling models differ, at min / the 5/25/50/75
    percentiles / max of eq/depth?

Why the first base: homopolymer run *interiors* are near-perfect on both models
and carry no discrimination. The first (5′/start) base of a run is where the
indel error concentrates — the model must decide whether the consensus has the
right *number* of repeats, and that decision is scored at the first base.

Steps:
  1. decide the ONE reference sequence the table is indexed against: the table
     must carry a single ``refname``; a single-contig FASTA is used as-is, and
     for a multi-contig FASTA the contig named by that ``refname`` is used
     (anything else is a hard error — run detection on a concatenation would
     invent targets on contigs the table does not cover),
  2. find every C/G (default) homopolymer run of that reference whose length is
     in ``--lens`` (default {4,5}) — ``pos`` 0-based,
  3. for each run keep the requested ``--which`` base (default FIRST),
  4. pull the matching (refname, pos) rows from the joined table,
  5. write a per-locus table, a per-(base, repeatCount) percentile comparison of
     eq/depth over the loci finite on BOTH sides (depth>0 on ctrl and exp), and
     a markdown report, and print both the percentile summary and a per-locus
     comparison (disagreeing-locus count + the most-divergent loci) to stdout.

Before printing, the script SELF-VERIFIES two invariants on the joined table and
aborts with a nonzero exit code if either fails (same invariants as the
``smc_ab_homopolymer_firstbase`` sibling, adapted to one reference):

  V1. the base bracketed ``[..]`` in ``aroundBases_exp`` equals the reference
      base at the 0-based ``pos`` column, for every row;
  V2. ``locus_accuracy_by_depth`` equals ``eq / depth`` where depth > 0, and
      ``eq <= depth`` everywhere.

Both checks are total: a row that cannot be checked (missing column, non-integer
``pos``, ``pos`` outside the reference) aborts rather than being skipped, and the
``[verify] OK`` line reports how many rows each check actually covered.

Outputs (under ``--outdir``), with
``<stem>`` = ``poly<CHAR>_lens<len-lens>_<which>base`` (e.g.
``polyCG_lens4-5_firstbase``):
  - <stem>_locus.tsv        selected rows + added base/repeatCount cols
  - <stem>_percentiles.tsv  per (base, repeatCount) p25/p50/p75 eq/depth A/B
  - <stem>_report.md        human-readable report (percentile table + per-locus)

Usage:
  python poly_first_base_analysis.py \
      --joined    .../locus_accuracy_joined.csv \
      --ref       .../MG1655.fa \
      --outdir    .../poly_first_base \
      --char CG --lens 4,5 --which first
"""

import argparse
import os
import sys

import polars as pl


def load_ref(path):
    """FASTA -> ({name: uppercase_seq}, [names in file order]).

    ``name`` is the first whitespace-delimited field after '>'."""
    seqs, order, cur = {}, [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
                order.append(cur)
                continue
            if cur is None:
                continue
            seqs[cur].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}, order


def pick_refseq(seqs, order, refname):
    """The single reference sequence the joined table's ``pos`` column indexes.

    Run detection must use the coordinates the table was built on. With a
    single-contig FASTA that sequence is the one contig (whatever its name, so
    a table whose refname differs from the FASTA header still works). With a
    multi-contig FASTA the contig has to be the table's ``refname`` — otherwise
    runs would be detected on a concatenation, and runs on the other contigs
    would become phantom targets that no table row can ever match.
    """
    if not order:
        sys.exit("no sequences parsed from the reference FASTA")
    if len(order) == 1:
        return seqs[order[0]]
    if refname in seqs:
        print(f"[info] reference has {len(order)} contigs; using contig "
              f"{refname!r} for run detection (the joined table's refname)")
        return seqs[refname]
    sys.exit(f"reference has {len(order)} contigs {order}, but the joined table's "
             f"refname {refname!r} is not one of them — cannot tell which contig "
             f"the pos column indexes. Pass the single-contig reference the run "
             f"was aligned to.")


def runs_of(seq, chars, lens):
    """Return list of (start_0based, run_len, base) for runs of a char in
    ``chars`` whose length is in ``lens``."""
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


def which_positions(start, run_len, which):
    if which == "first":
        return [start]
    if which == "last":
        return [start + run_len - 1]
    if which == "all":
        return list(range(start, start + run_len))
    raise ValueError(f"unknown --which: {which}")


def bracketed(around):
    """Char inside the [..] brackets of an aroundBases string, else None."""
    if "[" in around and "]" in around:
        return around[around.index("[") + 1:around.index("]")]
    return None


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def finite(v):
    """True for a real float/int that is not None/NaN."""
    return v is not None and v == v


def fmt4(v):
    """Fixed 4-dp, or 'n/a' for a value that is not finite (e.g. an empty group)."""
    return "n/a" if not finite(v) else f"{v:.4f}"


def fmt4s(v):
    """Signed fixed 4-dp, or 'n/a' for a value that is not finite."""
    return "n/a" if not finite(v) else f"{v:+.4f}"


def parse_pcts(spec):
    """'25,50,75' -> [25, 50, 75]. Each must be an integer rank in [0, 100]."""
    out = []
    for tok in spec.replace(" ", "").split(","):
        if not tok:
            continue
        try:
            q = int(tok)
        except ValueError:
            sys.exit(f"--pctl must be comma-separated integers in [0,100], got {spec!r}")
        if not 0 <= q <= 100:
            sys.exit(f"--pctl value out of range [0,100]: {q} (in {spec!r})")
        out.append(q)
    if not out:
        sys.exit(f"--pctl is empty: {spec!r}")
    return out


def verify(df, seq):
    """Run invariants V1 and V2 on the joined table.

    Returns (n_rows, n_v1, n_v2, error_or_None), where ``n_v1`` / ``n_v2`` are
    the row counts each check actually covered. A row that cannot be checked is
    a failure, not a skip: an unverifiable table must not be reported as
    verified.
    """
    n = df.height
    need = ["pos", "aroundBases_exp"]
    for side in ("ctrl", "exp"):
        need += [f"eq_{side}", f"depth_{side}", f"locus_accuracy_by_depth_{side}"]
    absent = [c for c in need if c not in df.columns]
    if absent:
        return n, 0, {}, (f"V1/V2 cannot run: joined table lacks column(s) {absent} — "
                          f"is this a locus_accuracy_joined.csv from "
                          f"smc_ab_test_ref_locus_single?")
    around = df["aroundBases_exp"].to_list()
    pos = df["pos"].to_list()
    n_v1 = 0
    for i in range(n):
        try:
            p = int(pos[i])
        except (TypeError, ValueError):
            return n, n_v1, {}, f"V1 FAIL: row {i} has non-integer pos={pos[i]!r}"
        if not 0 <= p < len(seq):
            return n, n_v1, {}, (f"V1 FAIL: row {i} pos={p} is outside the "
                                 f"reference (len={len(seq)}) — --ref is probably "
                                 f"not the reference the table was built on")
        b = bracketed(around[i])
        if b is None:
            return n, n_v1, {}, (f"V1 FAIL: row {i} pos={p} has no [..] in "
                                 f"aroundBases_exp={around[i]!r}")
        if b != seq[p]:
            return n, n_v1, {}, (f"V1 FAIL: pos={p} bracketed={b!r} ref[0based]={seq[p]!r} — "
                                 f"pos column may not be 0-based, or aroundBases mis-centered")
        n_v1 += 1
    n_v2 = {"ctrl": 0, "exp": 0}
    for side in ("ctrl", "exp"):
        eqcol, depcol, bycol = f"eq_{side}", f"depth_{side}", f"locus_accuracy_by_depth_{side}"
        eq = df[eqcol].to_list()
        dep = df[depcol].to_list()
        byd = df[bycol].to_list()
        for i in range(n):
            eq_i, dep_i = eq[i], dep[i]
            if eq_i is None or dep_i is None:
                continue
            try:
                eq_i = int(eq_i); dep_i = int(dep_i)
            except (TypeError, ValueError):
                return n, n_v1, n_v2, (f"V2 FAIL: non-integer eq/depth for {side} at "
                                       f"pos={pos[i]} (eq={eq[i]!r} depth={dep[i]!r})")
            if eq_i > dep_i:
                return n, n_v1, n_v2, (f"V2 FAIL: eq>depth for {side} at pos={pos[i]} "
                                       f"(eq={eq_i} depth={dep_i})")
            b = byd[i]
            if dep_i > 0 and b is not None and b == b:  # b present & not NaN
                if abs(b - eq_i / dep_i) > 1e-6:
                    return n, n_v1, n_v2, (f"V2 FAIL: by_depth={b} != eq/depth={eq_i/dep_i:.6f} "
                                           f"for {side} at pos={pos[i]}")
            n_v2[side] += 1
    return n, n_v1, n_v2, None


def pct(vals, q):
    """Linear-interpolation percentile (numpy default) of a finite-sorted list."""
    if not vals:
        return float("nan")
    idx = (len(vals) - 1) * q / 100.0
    lo = int(idx); hi = min(lo + 1, len(vals) - 1)
    frac = idx - lo
    return vals[lo] + (vals[hi] - vals[lo]) * frac


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--joined", required=True,
                    help="tab-separated locus_accuracy_joined.csv from "
                         "smc_ab_test_ref_locus_single")
    ap.add_argument("--ref", required=True,
                    help="single-reference FASTA the run was aligned to")
    ap.add_argument("--outdir", required=True, help="where outputs are written")
    ap.add_argument("--char", default="CG",
                    help="chars that count as the homopolymer (default CG)")
    ap.add_argument("--lens", default="4,5",
                    help="comma-separated run lengths to consider (default 4,5)")
    ap.add_argument("--which", default="first", choices=["first", "last", "all"],
                    help="which base of each run to compare (default first)")
    ap.add_argument("--pctl", default="5,25,50,75",
                    help="comma-separated percentile ranks (default 5,25,50,75). "
                         "min and max are ALWAYS reported in addition, as fixed columns")
    args = ap.parse_args()

    chars = set(args.char.upper())
    lens = set(int(x) for x in args.lens.split(",") if x.strip() != "")
    pcts = parse_pcts(args.pctl)
    char_tag = "".join(sorted(chars)) or args.char.upper()
    len_tag = "-".join(str(x) for x in sorted(lens))
    stem = f"poly{char_tag}_lens{len_tag}_{args.which}base"

    df = pl.read_csv(args.joined, separator="\t")
    if df.height == 0:
        sys.exit(f"joined table is empty (header only?): {args.joined}")
    if "refname" not in df.columns:
        sys.exit(f"joined table has no 'refname' column: {args.joined}")
    refnames = df["refname"].unique().to_list()
    if len(refnames) != 1:
        sys.exit(f"joined table has {len(refnames)} refnames: {refnames} — this skill "
                 f"compares ONE reference. For a per-barcode / multi-reference table "
                 f"use smc_ab_homopolymer_firstbase instead.")
    rn = refnames[0]

    seqs, order = load_ref(args.ref)
    seq = pick_refseq(seqs, order, rn)
    print(f"reference contig(s)={order} len={len(seq)}; "
          f"chars={sorted(chars)} lens={sorted(lens)} which={args.which}")
    print(f"joined rows={df.height}; refname={rn}")

    n_rows, n_v1, n_v2, err = verify(df, seq)
    if err:
        sys.exit(f"INPUT VERIFICATION FAILED ({n_rows} rows scanned):\n  {err}")
    print(f"[verify] OK: V1 checked {n_v1}/{n_rows} rows "
          f"(bracket==ref[0-based pos]); V2 checked ctrl {n_v2['ctrl']}/{n_rows} + "
          f"exp {n_v2['exp']}/{n_rows} rows (by_depth==eq/depth, eq<=depth)")

    # ---- select target (pos) set from runs ----
    targets = {}  # pos -> (base, run_len)
    for start, rl, base in runs_of(seq, chars, lens):
        for p in which_positions(start, rl, args.which):
            targets[p] = (base, rl)
    print(f"target {args.which}-base loci: {len(targets)}")

    sub = df.filter(pl.col("refname") == rn).filter(pl.col("pos").is_in(list(targets)))
    sub = sub.with_columns([
        pl.Series("base", [targets[int(p)][0] for p in sub["pos"].to_list()]),
        pl.Series("repeatCount", [targets[int(p)][1] for p in sub["pos"].to_list()]),
    ])
    missing = len(targets) - sub.height
    sub = sub.sort("pos")
    if missing:
        # Depth-0 loci ARE present in the joined table (gsetl emits a row for every
        # reference position, on both sides, and the join is inner on refname,pos),
        # so a nonzero count here is never "missing coverage" — it is a real key
        # mismatch (wrong reference, or a truncated per-side table).
        print(f"[warn] {missing}/{len(targets)} target pos had NO row in the joined "
              f"table. Depth-0 loci are present, so this is a key mismatch, not "
              f"missing coverage: {sorted(set(targets) - set(sub['pos'].to_list()))[:10]}")
    if sub.height == 0:
        sys.exit("no selected positions — check --char / --lens / --which")
    print(f"selected {sub.height}/{len(targets)} loci (missing={missing})")

    # ---- coverage breakdown: what the percentile comparison will drop ----
    def zero_depth(col):
        return pl.col(col).is_null() | (pl.col(col) == 0)

    d0_ctrl, d0_exp = zero_depth("depth_ctrl"), zero_depth("depth_exp")
    n_d0_ctrl = sub.filter(d0_ctrl).height
    n_d0_exp = sub.filter(d0_exp).height
    n_d0_either = sub.filter(d0_ctrl | d0_exp).height
    n_both = sub.height - n_d0_either
    print(f"depth-0 (no coverage at the run edge): ctrl={n_d0_ctrl}, exp={n_d0_exp}, "
          f"either side={n_d0_either} -> finite-both used for the percentiles={n_both}")

    os.makedirs(args.outdir, exist_ok=True)

    # ---- (1) per-locus TSV ----
    per_cols = ["refname", "pos", "base", "repeatCount", "aroundBases_exp",
                "eq_ctrl", "depth_ctrl", "locus_accuracy_by_depth_ctrl",
                "diff_ctrl", "ins_ctrl", "del_ctrl", "locus_accuracy_ctrl",
                "eq_exp", "depth_exp", "locus_accuracy_by_depth_exp",
                "diff_exp", "ins_exp", "del_exp", "locus_accuracy_exp"]
    per = sub.select(per_cols).rename({"pos": "locus"})
    per_path = os.path.join(args.outdir, f"{stem}_locus.tsv")
    per.write_csv(per_path, separator="\t")

    # ---- (2) percentile comparison per (base, repeatCount) of eq/depth ----
    # eq/depth is undefined (NaN) where depth==0. control and experiment are two
    # separate consensus assemblies of the same SMC run, so a run edge may be
    # covered on one side and depth-0 on the other (whether it is — and for how
    # many loci — is data-dependent; the breakdown printed above says). The
    # per-side percentiles are therefore restricted to the loci finite on BOTH
    # sides, so ctrl and exp are distributed over the SAME set of positions and
    # the per-percentile delta is a like-for-like comparison. This is the same
    # finite-both rule the per-locus section below uses. n = the finite-both
    # locus count for the group.
    order = sorted({(b, L) for (b, L) in targets.values()},
                   key=lambda t: (t[0], t[1]))

    # Per-side statistics, in display order: min, the requested percentiles
    # (ascending), then max. min and max are ALWAYS reported, in addition to the
    # --pctl ranks.
    def keytag(key):
        return f"p{key}" if isinstance(key, int) else key  # 'min' / 'max' / 'p{q}'

    stat_keys = ["min"] + sorted(pcts) + ["max"]

    def stat_val(sorted_vals, key):
        # sorted_vals is already finite-both and ascending; '' -> NaN
        if not sorted_vals:
            return float("nan")
        if key == "min":
            return round(sorted_vals[0], 6)
        if key == "max":
            return round(sorted_vals[-1], 6)
        return round(pct(sorted_vals, key), 6)

    pctl_cols = ["base", "repeatCount"]
    for t in ("ctrl", "exp"):
        for key in stat_keys:
            pctl_cols.append(f"eq_over_depth_{t}_{keytag(key)}")
    for key in stat_keys:
        pctl_cols.append(f"delta_{keytag(key)}")
    pctl_cols.append("n")
    pctl_rows = []
    for (b, L) in order:
        gsub = sub.filter((pl.col("base") == b) & (pl.col("repeatCount") == L))
        ctrl_vals = gsub["locus_accuracy_by_depth_ctrl"].to_list()
        exp_vals = gsub["locus_accuracy_by_depth_exp"].to_list()
        both_idx = [i for i in range(len(ctrl_vals))
                    if finite(ctrl_vals[i]) and finite(exp_vals[i])]
        cvals = sorted(ctrl_vals[i] for i in both_idx)
        evals = sorted(exp_vals[i] for i in both_idx)
        cgrid = {key: stat_val(cvals, key) for key in stat_keys}
        egrid = {key: stat_val(evals, key) for key in stat_keys}
        row = [b, L]
        for key in stat_keys:
            row.append(cgrid[key])
        for key in stat_keys:
            row.append(egrid[key])
        # Δ is the difference of the two already-rounded values, so a printed
        # Δ is exactly the difference of the two printed figures (6 dp).
        for key in stat_keys:
            row.append(round(egrid[key] - cgrid[key], 6))
        row.append(len(both_idx))
        pctl_rows.append(row)
    pctl = pl.DataFrame(pctl_rows, schema=pctl_cols, orient="row")
    pctl_path = os.path.join(args.outdir, f"{stem}_percentiles.tsv")
    pctl.write_csv(pctl_path, separator="\t")

    # ---- markdown report ----
    md = []
    md.append(f"# Poly {args.char.upper()} homopolymer (run len {args.lens}) — "
              f"{args.which}-base eq/depth A/B")
    md.append("")
    md.append(f"- joined: `{args.joined}`")
    md.append(f"- reference: `{os.path.basename(args.ref)}` (contig len {len(seq)})")
    md.append(f"- target: {args.char.upper()} runs of length **{args.lens}**, "
              f"{args.which} base of each run")
    md.append(f"- loci selected: **{sub.height}** of {len(targets)} (missing {missing})")
    md.append(f"- depth-0 (no coverage at the run edge): ctrl **{n_d0_ctrl}**, exp "
              f"**{n_d0_exp}**, either side **{n_d0_either}** → **{n_both}** loci are "
              f"finite on both sides, i.e. the set the percentiles below are taken over")
    md.append("")
    md.append("Metric `locus_accuracy_by_depth` = eq/depth per reference base "
              "position, over the loci finite on BOTH sides (depth > 0 on ctrl "
              "and exp); per-percentile Δ is a like-for-like comparison. "
              "Δ = exp − ctrl. **eq/depth ignores insertions** — an extra repeat "
              "unit at the run edge is recorded in the `ins_*` columns and does "
              "not lower eq/depth (see the stdout insertion tally).")
    md.append("")
    md.append("## (2) Per (base, repeatCount) percentile comparison of eq/depth")
    md.append("")
    hdr = ["group", "n"]
    for t in ("ctrl", "exp"):
        hdr += [f"{t[0].upper()}{t[1]} {keytag(key)}" for key in stat_keys]
    hdr += [f"Δ {keytag(key)}" for key in stat_keys]
    md.append("| " + " | ".join(hdr) + " |")
    md.append("|" + "---|" * len(hdr))
    for r in pctl.iter_rows():
        b, L = r[0], r[1]
        cells = [f"{b}{L}", str(r[-1])]            # group, n (last col)
        # order: ctrl stats, exp stats, deltas  (all in stat_keys order)
        k = 2
        cells += [fmt4(r[k+i]) for i in range(len(stat_keys))]          # ctrl
        cells += [fmt4(r[k+len(stat_keys)+i]) for i in range(len(stat_keys))]  # exp
        k = 2 + 2 * len(stat_keys)
        cells += [fmt4s(r[k+i]) for i in range(len(stat_keys))]          # delta
        md.append("| " + " | ".join(cells) + " |")
    md.append("")
    md.append(f"## (1) Per-locus {args.which}-base table")
    md.append("")
    md.append(f"Full table: `{per_path}` ({per.height} rows).")
    md.append("Columns: refname, locus, base, repeatCount, aroundBases([ ]=target "
              "base), eq/depth/eq÷depth/diff/ins/del/locus_accuracy for _ctrl then _exp.")
    md.append("")
    md.append("| refname | locus | base | rc | aroundBases | eq_c | d_c | eq/dep_c "
              "| diff_c | ins_c | del_c | eq_x | d_x | eq/dep_x | diff_x | ins_x | del_x |")
    md.append("|---|---:|:---:|:---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in per.iter_rows():
        (refname, locus, base, rc, around,
         eq_c, d_c, ed_c, diff_c, ins_c, del_c, acc_c,
         eq_x, d_x, ed_x, diff_x, ins_x, del_x, acc_x) = r
        fc = "" if not finite(ed_c) else f"{ed_c:.3f}"
        fx = "" if not finite(ed_x) else f"{ed_x:.3f}"
        md.append(f"| {refname} | {locus} | {base} | {rc} | {around} "
                  f"| {eq_c} | {d_c} | {fc} | {diff_c} | {ins_c} | {del_c} "
                  f"| {eq_x} | {d_x} | {fx} | {diff_x} | {ins_x} | {del_x} |")
    md.append("")
    rep_path = os.path.join(args.outdir, f"{stem}_report.md")
    with open(rep_path, "w") as fp:
        fp.write("\n".join(md))

    print(f"[write] {per_path} ({per.height} rows)")
    print(f"[write] {pctl_path}")
    print(f"[write] {rep_path}")
    print("\n=== percentiles (eq/depth, finite only) ===")
    with open(pctl_path) as fp:
        print(fp.read())

    # ---- per-locus comparison summary (stdout) ----
    # The percentile block above is the aggregate readout; this is the
    # locus-granularity view the user also wants. The full per-locus comparison
    # lives in <stem>_locus.tsv / report.md; here we surface how many loci the
    # two models disagree on plus the most-divergent ones (bounded, so a
    # whole-genome run of ~14k loci stays readable).
    eps = 1e-9
    fin_both = 0
    disagree = []
    better = worse = equal = 0   # exp vs ctrl on eq/depth, finite-both loci
    tot = {t: {'eq': 0, 'diff': 0, 'ins': 0, 'del': 0, 'depth': 0}
           for t in ('ctrl', 'exp')}
    ins_any = 0       # finite-both loci with ins>0 on at least one side
    ins_masked = 0    # ... of those, eq/depth == 1.0 on BOTH (invisible to eq/depth)
    eqs = {}   # locus -> {'c': eq_ctrl, 'e': eq_exp}
    deps = {}  # locus -> {'c': depth_ctrl, 'e': depth_exp}
    for r in per.iter_rows():
        (refname, locus, base, rc, around,
         eq_c, d_c, ed_c, diff_c, ins_c, del_c, acc_c,
         eq_x, d_x, ed_x, diff_x, ins_x, del_x, acc_x) = r
        if not (finite(ed_c) and finite(ed_x)):
            continue
        fin_both += 1
        delta = ed_x - ed_c
        if delta > eps:
            better += 1
        elif delta < -eps:
            worse += 1
        else:
            equal += 1
        eqs[locus] = {'c': eq_c, 'e': eq_x}
        deps[locus] = {'c': d_c, 'e': d_x}
        for t, e, i, d, dep in (('ctrl', eq_c, ins_c, del_c, d_c),
                                ('exp', eq_x, ins_x, del_x, d_x)):
            tot[t]['eq'] += int(e)
            tot[t]['ins'] += int(i)
            tot[t]['del'] += int(d)
            tot[t]['depth'] += int(dep)
        tot['ctrl']['diff'] += int(diff_c)
        tot['exp']['diff'] += int(diff_x)
        if int(ins_c) > 0 or int(ins_x) > 0:
            ins_any += 1
            if ed_c == 1.0 and ed_x == 1.0:
                ins_masked += 1
        if abs(delta) > eps:
            disagree.append((abs(delta), locus, base, rc, around, ed_c, ed_x, delta))
    disagree.sort(reverse=True)
    topn = 20
    print("\n=== per-locus comparison (eq/depth, finite both sides) ===")
    print(f"full per-locus table: {per_path} ({per.height} rows)")
    print(f"finite-both loci: {fin_both}; "
          f"disagreeing (|Δ eq/depth| > 0): {len(disagree)}")
    print("totals over finite-both loci (two lenses at the same loci):")
    for t in ('ctrl', 'exp'):
        tt = tot[t]
        denom = tt['eq'] + tt['diff'] + tt['ins'] + tt['del']
        print(f"  {t}: eq={tt['eq']} depth={tt['depth']} diff={tt['diff']} "
              f"ins={tt['ins']} del={tt['del']} | "
              f"eq/depth={(tt['eq'] / tt['depth']) if tt['depth'] else float('nan'):.4f} "
              f"| locus_accuracy=eq/(eq+diff+ins+del)="
              f"{(tt['eq'] / denom) if denom else float('nan'):.6f}")
    print(f"insertion-carrying loci: {ins_any} of {fin_both} finite-both loci have "
          f"ins>0 on at least one side; {ins_masked} of those have eq/depth == 1.0 "
          f"on BOTH sides, so eq/depth (the percentiles and the tally below) is "
          f"blind to their insertion difference — use the per-locus ins_* columns.")
    print(f"win/lose/tie (finite-both, Δ eq/depth = exp − ctrl): "
          f"exp better={better} "
          f"({(better / fin_both * 100) if fin_both else float('nan'):.2f}%), "
          f"exp worse={worse} "
          f"({(worse / fin_both * 100) if fin_both else float('nan'):.2f}%), "
          f"equal={equal} "
          f"({(equal / fin_both * 100) if fin_both else float('nan'):.2f}%)")
    print(f"top {min(topn, len(disagree))} most-divergent loci "
          f"(ctrl / exp eq and depth, Δ = exp − ctrl):")
    print("  locus  base rc  aroundBases   ctrl eq/dep (eq/depth)   "
          "exp eq/dep (eq/depth)   Δ")
    for _absd, locus, base, rc, around, ed_c, ed_x, delta in disagree[:topn]:
        print(f"  {locus}  {base}  {rc}  {around}   "
              f"c: {eqs[locus]['c']} / {deps[locus]['c']} "
              f"({ed_c:.4f})   "
              f"e: {eqs[locus]['e']} / {deps[locus]['e']} "
              f"({ed_x:.4f})   {delta:+.4f}")
    if not disagree:
        print("  (no loci where the two models differ in eq/depth)")


if __name__ == "__main__":
    main()
