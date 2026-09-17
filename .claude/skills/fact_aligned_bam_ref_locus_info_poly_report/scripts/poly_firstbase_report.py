#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Poly-C/G first-base eq/depth report over a gsetl fact locus table.

Reads a gsetl whole-sample per-locus table (``fact_aligned_bam_ref_locus_info.csv``
style: columns ``refname pos eq diff ins del depth curBase nextBase curIsHomo
nextIsHomo aroundBases diffDetail insDetail``), a multi-contig reference FASTA,
and reports — for every **maximal** poly-C / poly-G run whose length is in
``--run-lengths`` (default ``4,5``) — the **first base** of each run (0-based
``pos``), joined against the table's ``eq/depth``.

Unlike the sibling ``smc_homopolymer_firstbase`` skill this one:
  * takes the gsetl **fact table directly** (``--table``), keyed by
    ``(refname, pos)`` — there is NO ``barcode`` column, because this is the
    whole-sample table, not a per-barcode A/B table;
  * takes a **multi-contig** reference FASTA and joins on the first FASTA
    header field (``>NAME ...`` -> ``NAME``), matching the table's ``refname``;
  * filters runs by **exact** run length (``--run-lengths 4,5``), not a
    ``>= min-run`` floor;
  * reports **eq/depth** (plus the per-contig and distribution roll-ups we want
    to eyeball), and explicitly separates the depth-0 loci (eq/depth = 0/0,
    reported as 0.0000) from the genuinely low-accuracy loci (short depth +
    homopolymer-length indels).

Self-verification (runs before any output; aborts nonzero on failure):
  V1 — for every selected run, the base bracketed ``[..]`` in the table row's
        ``aroundBases`` equals the reference base at the 0-based ``pos``. This
        is the off-by-one guard: gsetl ``pos`` is 0-based (pos 0 = first base);
        a shifted (1-based) join would otherwise silently select the *second*
        base of each run.
  V2 — ``eq <= depth`` for every selected row.

Outputs (written to ``--outdir``):
  - poly_firstbase_locus.tsv   (one row per selected locus: refname pos base
                                depth eq diff eq/depth aroundBases)
  - poly_firstbase_report.md   (per-contig table, eq/depth distribution,
                                lowest eq/depth loci, depth-0 loci)
  (stdout: the same summary so it can be pasted into a chat message.)

Usage:
  python poly_firstbase_report.py \
      --table .../gsetl/fact_aligned_bam_ref_locus_info.csv \
      --ref   .../ref/E_ATCC_25922.fasta \
      --outdir .../gsetl/poly_report \
      --char CG --run-lengths 4,5 --which first
"""

import argparse
import csv
import os
import statistics
from collections import Counter, defaultdict


def read_fasta(path):
    """Multi-record FASTA -> ({name: uppercase_seq}, [names in order]).
    Name = first whitespace-delimited field after '>'."""
    seqs, order, cur = {}, [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
                order.append(cur)
            else:
                if cur is None:
                    continue
                seqs[cur].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}, order


def bracketed(around):
    """The character inside the ``[..]`` brackets of an aroundBases string,
    or None if there is none."""
    if around and "[" in around and "]" in around:
        return around[around.index("[") + 1:around.index("]")]
    return None


def maximal_runs(seq, chars):
    """List of (start_0based, run_len, base) for every maximal run of a single
    character in ``chars`` (a run = longest stretch of one identical char)."""
    out = []
    i, n = 0, len(seq)
    while i < n:
        j = i
        while j + 1 < n and seq[j + 1] == seq[i]:
            j += 1
        if seq[i] in chars:
            out.append((i, j - i + 1, seq[i]))
        i = j + 1
    return out


def parse_run_lengths(spec):
    """``"4,5"`` -> ``{4, 5}``."""
    out = set()
    for tok in str(spec).replace(" ", "").split(","):
        if tok:
            out.add(int(tok))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True,
                    help="gsetl whole-sample per-locus fact table "
                         "(fact_aligned_bam_ref_locus_info.csv), tab-separated")
    ap.add_argument("--ref", required=True,
                    help="multi-contig reference FASTA (header field 1 = refname)")
    ap.add_argument("--outdir", required=True, help="where outputs are written")
    ap.add_argument("--char", default="CG",
                    help="characters that count as the homopolymer (default CG)")
    ap.add_argument("--run-lengths", default="4,5",
                    help="comma list of EXACT maximal run lengths to keep "
                         "(default 4,5)")
    ap.add_argument("--which", default="first", choices=["first", "last"],
                    help="which base of each matching run to take (default first)")
    ap.add_argument("--worst", type=int, default=25,
                    help="how many lowest eq/depth loci to show (default 25)")
    args = ap.parse_args()

    chars = set(args.char.upper())
    lengths = parse_run_lengths(args.run_lengths)
    seqs, order = read_fasta(args.ref)
    if not order:
        sys_exit(f"no sequences parsed from {args.ref}")
    seqs = {k: seqs[k] for k in order}

    # ---- find candidate (refname, pos) pairs on the reference ----
    run_at = {}        # (name, pos) -> (base, run_len)
    per_contig_total = Counter()
    for name in order:
        for start, rl, base in maximal_runs(seqs[name], chars):
            if rl not in lengths:
                continue
            per_contig_total[name] += 1
            p = start if args.which == "first" else start + rl - 1
            run_at[(name, p)] = (base, rl)

    # ---- single streaming pass over the (large) table ----
    matched = {}       # (name, pos) -> row dict
    with open(args.table, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames
        for row in reader:
            key = (row["refname"], int(row["pos"]))
            if key in run_at:
                matched[key] = row
    missing = [k for k in run_at if k not in matched]

    if not matched:
        sys_exit("no candidate loci found in the table — check "
                 "--ref / --char / --run-lengths")

    # ---- self-verification ----
    n_checked = 0
    for (name, pos), row in matched.items():
        base, _ = run_at[(name, pos)]
        n_checked += 1
        # V1: bracketed char == ref[0-based pos]
        b = bracketed(row.get("aroundBases", ""))
        refbase = seqs[name][pos] if 0 <= pos < len(seqs[name]) else None
        if b is not None and refbase is not None and b != refbase:
            sys_exit(f"INPUT VERIFICATION FAILED (V1): {name}:{pos} "
                     f"bracketed={b!r} ref[0-based]={refbase!r} — pos column may "
                     f"not be 0-based, or aroundBases is mis-centred")
        if base is not None and refbase is not None and base != refbase:
            sys_exit(f"INPUT VERIFICATION FAILED (V1): {name}:{pos} run base="
                     f"{base!r} ref[0-based]={refbase!r}")
        # V2: eq <= depth
        eq, dep = int(row["eq"]), int(row["depth"])
        if eq > dep:
            sys_exit(f"INPUT VERIFICATION FAILED (V2): {name}:{pos} "
                     f"eq={eq} > depth={dep}")
    print(f"[verify] OK: {n_checked}/{len(run_at)} candidate loci pass V1 and V2 "
          f"({len(missing)} had no table row)")

    # ---- build selected rows ----
    sel = []
    for (name, pos) in sorted(run_at):
        base, rl = run_at[(name, pos)]
        row = matched.get((name, pos))
        if row is None:
            sel.append(dict(rname=name, pos=pos, pos1=pos + 1, base=base,
                            runlen=rl, depth=0, eq=0, diff=0, ratio=0.0,
                            aroundBases="", note="NOT_IN_CSV"))
            continue
        eq, diff, depth = int(row["eq"]), int(row["diff"]), int(row["depth"])
        ratio = (eq / depth) if depth else 0.0
        sel.append(dict(rname=name, pos=pos, pos1=pos + 1, base=base, runlen=rl,
                        depth=depth, eq=eq, diff=diff, ratio=ratio,
                        aroundBases=row.get("aroundBases", ""), note=""))

    # ---- roll-ups ----
    withdepth = [r for r in sel if r["depth"] > 0]
    depths = [r["depth"] for r in withdepth]
    ratios = [r["ratio"] for r in withdepth]
    diffs = [r["diff"] for r in sel]

    per_in_csv = Counter()     # name -> candidate loci that had a table row
    for (name, _pos) in matched:
        per_in_csv[name] += 1
    per = defaultdict(lambda: [0, 0])   # name -> [loci, diff>0]
    for r in sel:
        per[r["rname"]][0] += 1
        if r["diff"] > 0:
            per[r["rname"]][1] += 1

    buckets = [
        (1.0, 1.0, "== 1.0 (perfect)"),
        (0.95, 0.999999, "0.95-0.99"),
        (0.90, 0.949999, "0.90-0.94"),
        (0.80, 0.899999, "0.80-0.89"),
        (0.0, 0.799999, "< 0.80"),
    ]
    bucket_counts = [(lab, sum(1 for x in ratios if lo <= x <= hi))
                     for lo, hi, lab in buckets]

    zero_depth = [r for r in sel if r["depth"] == 0]
    worst = sorted(withdepth, key=lambda r: (r["ratio"], -r["diff"]))[:args.worst]

    # ---- write TSV ----
    os.makedirs(args.outdir, exist_ok=True)
    tsv_path = os.path.join(args.outdir, "poly_firstbase_locus.tsv")
    with open(tsv_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["refname", "pos", "pos1", "base", "runlen", "depth",
                    "eq", "diff", "eq_div_depth", "aroundBases"])
        for r in sel:
            w.writerow([r["rname"], r["pos"], r["pos1"], r["base"], r["runlen"],
                        r["depth"], r["eq"], r["diff"], f"{r['ratio']:.4f}",
                        r["aroundBases"]])

    # ---- write markdown report ----
    rep_path = os.path.join(args.outdir, "poly_firstbase_report.md")
    L = []
    L.append(f"# Poly-C/G first-base eq/depth (char={args.char.upper()}, "
             f"run-len in {sorted(lengths)}, which={args.which})")
    L.append("")
    L.append(f"- reference: `{args.ref}` ({len(order)} contigs)")
    L.append(f"- table: `{args.table}`")
    L.append(f"- definition: maximal poly-{'/'.join(sorted(chars))} run with "
             f"length in {sorted(lengths)}; take the **{args.which} base** "
             f"(0-based `pos`); join to the table's eq/depth/aroundBases.")
    L.append(f"- candidate loci: **{len(sel)}** "
             f"({len(matched)} matched in CSV, {len(missing)} not present)")
    L.append("")
    L.append("## Per contig")
    L.append("")
    L.append("| contig | candidate | in-CSV | diff>0 | share |")
    L.append("|---|---:|---:|---:|---:|")
    for name in order:
        tot = per_contig_total.get(name, 0)
        in_csv = per_in_csv.get(name, 0)
        c, e = per.get(name, [0, 0])
        L.append(f"| {name} | {tot} | {in_csv} | {e} | {100 * e / c if c else 0:.2f}% |")
    L.append("")
    L.append("## eq/depth distribution (depth > 0 only)")
    L.append("")
    L.append("| eq/depth | loci | share |")
    L.append("|---|---:|---:|")
    nwd = len(withdepth) or 1
    for lab, k in bucket_counts:
        L.append(f"| {lab} | {k} | {100 * k / nwd:.2f}% |")
    L.append("")
    if ratios:
        L.append(f"- overall (depth>0): min={min(ratios):.4f} "
                 f"median={statistics.median(ratios):.4f} "
                 f"mean={statistics.mean(ratios):.4f}; "
                 f"depth min={min(depths)} max={max(depths)} "
                 f"median={statistics.median(depths):.0f}; total diff={sum(diffs)}")
    L.append("")
    if zero_depth:
        L.append("## Depth-0 loci (no reads cover the run edge)")
        L.append("")
        L.append("These report eq/depth = 0/0 = 0.0000; they are **missing "
                 "coverage, not errors**.")
        L.append("")
        L.append("| rname | pos | pos1 | base | runlen | aroundBases |")
        L.append("|---|---:|---:|:-|---:|---|")
        for r in zero_depth:
            L.append(f"| {r['rname']} | {r['pos']} | {r['pos1']} | {r['base']} "
                     f"| {r['runlen']} | `{r['aroundBases']}` |")
        L.append("")
    L.append(f"## Lowest eq/depth loci (depth > 0, top {args.worst})")
    L.append("")
    L.append("| rname | pos | pos1 | base | depth | eq/depth | diff | aroundBases |")
    L.append("|---|---:|---:|:-|---:|---:|---:|---|")
    for r in worst:
        L.append(f"| {r['rname']} | {r['pos']} | {r['pos1']} | {r['base']} "
                 f"| {r['depth']} | {r['ratio']:.4f} | {r['diff']} | "
                 f"`{r['aroundBases']}` |")
    L.append("")
    with open(rep_path, "w") as f:
        f.write("\n".join(L))

    # ---- stdout summary ----
    print(f"\n=== poly-{ '/'.join(sorted(chars)) } run-len in {sorted(lengths)} "
          f"first/last base, which={args.which} ===")
    print(f"candidate loci: {len(sel)}  in-CSV: {len(matched)}  not-present: {len(missing)}")
    print(f"per contig:")
    for name in order:
        c, e = per.get(name, [0, 0])
        print(f"  {name:26s} loci={c:6d}  diff>0={e:5d}  {100 * e / c if c else 0:6.2f}%")
    print(f"eq/depth (depth>0, n={len(withdepth)}):")
    for lab, k in bucket_counts:
        print(f"  {lab:18s}: {k:6d}  ({100 * k / nwd:5.2f}%)")
    if ratios:
        print(f"  min={min(ratios):.4f} median={statistics.median(ratios):.4f} "
              f"mean={statistics.mean(ratios):.4f}; "
              f"depth min={min(depths)} max={max(depths)} median={statistics.median(depths):.0f}")
    print(f"depth-0 loci: {len(zero_depth)}  total diff: {sum(diffs)}")
    print(f"\n[wrote] {tsv_path}  ({len(sel)} rows)")
    print(f"[wrote] {rep_path}")


def sys_exit(msg):
    import sys
    sys.exit(f"ERROR: {msg}")


if __name__ == "__main__":
    main()
