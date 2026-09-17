---
name: fact-aligned-bam-ref-locus-info-poly-report
description: Report the first-base eq/depth for poly-C/G homopolymer runs (exact repeat count, default 4 or 5) by joining a gsetl whole-sample fact locus table (fact_aligned_bam_ref_locus_info.csv) to a multi-contig reference FASTA. Outputs per-contig tables, an eq/depth distribution, the lowest-accuracy loci, and a clean split between depth-0 (no coverage) and short-depth indel loci.
---

# Poly-C/G first-base eq/depth report (gsetl fact table)

## Purpose

Given a **gsetl whole-sample per-locus fact table** (the
`fact_aligned_bam_ref_locus_info.csv` that `gsetl aligned-bam` writes, columns
`refname pos eq diff ins del depth curBase nextBase curIsHomo nextIsHomo
aroundBases diffDetail insDetail`) and the **reference FASTA** it was aligned
to, this skill finds every **maximal poly-C / poly-G run** whose **exact** length
is in `--run-lengths` (default `4,5`), takes the **first base** of each run
(0-based `pos`), and reports that base's `eq / depth` against the table.

It answers the recurring question in the 大肠 / E. coli STR work:

> **At the first base of each 4–5-mer polyC/G homopolymer, how well do the
> reads agree with the reference (eq/depth)?**

The homopolymer run *interior* is usually near-perfect. The **first base** is
the indel-sensitive edge: length-variant reads (calling 3 or 6 copies instead of
the reference 4/5) mis-align there, so the shortfall in `eq/depth` shows up as
`ins`/`del` in the table rather than as base mismatches. This skill isolates
exactly those edge positions and makes that signal visible.

This skill is **pure table post-process** — no gsmm2 / gsetl is run. The table
must already exist (produce it with `gsetl aligned-bam`, or use the sibling
skills that call it).

## How it differs from the sibling `smc_homopolymer_firstbase` skill

| | `smc_homopolymer_firstbase` | this skill |
|---|---|---|
| input | A/B per-barcode table (has a `barcode` col) **or** a raw query+ref it aligns itself | the gsetl **whole-sample fact table** directly — **no `barcode` column**, keyed by `(refname, pos)` |
| reference | single-record plasmid refs | a **multi-contig** FASTA; joins on the header's first field |
| run filter | `>= --min-run` floor | **exact** `--run-lengths 4,5` |
| metric | `locus_accuracy` **and** `locus_accuracy_by_depth` (pooled) | `eq/depth` per locus + per-contig + distribution roll-ups |
| coverage gaps | skipped with a `[warn]` | depth-0 loci are **kept and reported as their own section** (0/0 = 0.0000 = no coverage, not an error) |

## When to use

- you already have a `fact_aligned_bam_ref_locus_info.csv` (or any gsetl fact
  table with the same columns) and want the **polyC/G first-base eq/depth** at
  repeat-count 4 and 5 — the whole-sample homopolymer-edge agreement number.
- you want to see which specific homopolymer edges are dragging the number down
  (the lowest-eq/depth loci) and to tell **missing coverage** (depth 0) apart
  from **short-depth indel** loci.
- you want to re-run this the moment the table is refreshed (it is idempotent
  and single-pass, ~a few seconds for the ~5 M-row E. coli table).

## Required inputs

- `--table` — the gsetl whole-sample per-locus fact table, tab-separated, with
  at least the columns `refname pos eq diff ins del depth aroundBases`. This is
  `fact_aligned_bam_ref_locus_info.csv`. It does **not** need a `barcode`
  column (that's the per-barcode A/B table, a different format).
- `--ref` — the **multi-contig reference FASTA** the table was aligned against.
  Each contig's name is the first whitespace-delimited field after `>`
  (e.g. `>ATCC_25922_contig_4 species=...` → `ATCC_25922_contig_4`), and must
  match the table's `refname` column exactly.
- `--outdir` — where the two output files are written.

## Optional inputs

- `--char` — the characters that count as the homopolymer. **Default `CG`**
  (polyC/G). Pass e.g. `ACGT` for all homopolymers.
- `--run-lengths` — comma list of **exact** maximal run lengths to keep.
  **Default `4,5`** (the polyC/G 4- and 5-mers). A run of length 6 is *not*
  matched (it is its own maximal run), and a length-3 run is not matched either.
- `--which` — `first` (default) | `last`: which base of each matching run to
  take. `first` = the run's 5′ edge (the indel-sensitive position).
- `--worst` — how many lowest-eq/depth loci to show in the report. **Default 25.**

## Execution

```bash
/root/miniconda3/envs/gseda/bin/python \
    /root/projects/gsda/.claude/skills/fact_aligned_bam_ref_locus_info_poly_report/scripts/poly_firstbase_report.py \
    --table    .../大肠/DPN8205/gsetl/fact_aligned_bam_ref_locus_info.csv \
    --ref      /data1/ccs_data/2026型检/ref/E_ATCC_25922.fasta \
    --outdir   .../大肠/DPN8205/gsetl/poly_report \
    --char CG --run-lengths 4,5 --which first
```

> **Use the `gseda` env python** (`/root/miniconda3/envs/gseda/bin/python`) —
> `server/bin` shadows `$PATH` on this machine so `conda activate` / `conda run`
> silently keep the py3.8 `server` python. The script itself only needs the stdlib
> (`argparse`, `csv`, `statistics`, `collections`), so any of these would work —
> but use the env python to be consistent with the rest of the project.

Single pass over the (large) table; a ~5.2 M-row E. coli genome table finishes in
a few seconds. No backgrounding needed, though you may run it as a background
task and poll `--outdir` if the table is huge.

## Output (under `--outdir`)

- `poly_firstbase_locus.tsv` — one row per selected locus:
  `refname  pos  pos1  base  runlen  depth  eq  diff  eq_div_depth  aroundBases`.
  `pos` is 0-based (gsetl convention); `pos1` = pos+1 for convenience.
- `poly_firstbase_report.md` — human-readable report:
  1. the definition and candidate/matched counts,
  2. a **per-contig** table (candidate / in-CSV / diff>0 / share),
  3. the **eq/depth distribution** (perfect / 0.95-0.99 / 0.90-0.94 / 0.80-0.89 /
     <0.80), over depth>0 loci only,
  4. a **depth-0 loci** section (missing coverage — reported explicitly, not
     mixed into the accuracy distribution),
  5. the **lowest eq/depth loci** (depth>0), with `aroundBases` (the bracketed
     letter is the target base).

The script also prints the same summary to stdout for easy pasting into a chat
message.

## Built-in verification (runs before any output — do not skip)

Before it writes anything, the script **self-checks the joined rows** and exits
nonzero if either invariant fails:

- **V1 — indexing.** For every selected run, the base bracketed `[..]` in the
  table row's `aroundBases` must equal the reference base at the **0-based**
  `pos`. This is the off-by-one guard: gsetl `pos` is **0-based** (pos 0 = first
  base). If a 1-based coordinate ever got joined into the 0-based `pos` column,
  V1 fails immediately instead of silently selecting the *second* base of each
  run.
- **V2 — metric integrity.** `eq <= depth` for every selected row.

On success it prints `[verify] OK: N/M candidate loci pass V1 and V2
(K had no table row)`. If it prints `ERROR: INPUT VERIFICATION FAILED …` and
exits nonzero, **stop and fix the input** — do not trust a table with a shifted
`pos`.

## Presenting the final result

When reporting to the user:

- State the **run filter** explicitly — e.g. "first base of each maximal
  polyC/G run of **exact** length 4 or 5" — so it's clear a length-6 run was
  *not* counted and a length-3 run was not counted.
- Lead with the **candidate count** and how many matched in the table.
- Give the **eq/depth distribution** (the perfect / 0.95-0.99 / … buckets) and
  the overall median/mean — this is the headline "how well do the 4/5-mers agree"
  number.
- Show the **lowest eq/depth loci** with `aroundBases`.
- **Separate the depth-0 loci explicitly.** They are *missing coverage* (eq/depth
  = 0/0, reported 0.0000), **not** accuracy errors. Do not present them as
  "0% accuracy" — say they are uncovered run edges.
- Note the interpretation: the genuinely-low (depth>0) loci are typically
  **indel-driven** (homopolymer length variants: reads calling 3 or 6 copies),
  not base mismatches — check the `diff`/`ins`/`del` columns to confirm.
- The **absolute path** of `poly_firstbase_report.md` and
  `poly_firstbase_locus.tsv` must appear in the final message.

### Final message requirement

Include the **absolute path** of `poly_firstbase_report.md`, the candidate
count, and the eq/depth distribution headline (perfect-count and median).

## Notes / gotchas

- **`pos` is 0-based** in the gsetl fact table (pos 0 = first reference base).
  All run-boundary math uses 0-based indices; V1 enforces it.
- **Whole-sample table = no `barcode` column.** Do not try to run the sibling
  `smc_homopolymer_firstbase` skill's mode-1 path on this table — it expects a
  `barcode` column and a `--mapping`. This skill keys on `(refname, pos)` only.
- **Exact length, not a floor.** `--run-lengths 4,5` keeps *exactly* 4- and
  5-mers. A 6-mer is one maximal run of length 6 and is excluded. If you actually
  want "runs of length ≥ 4", that is a different query — say so and the filter
  changes (the `smc_homopolymer_firstbase` skill's `--min-run` is the floor
  variant).
- **The FASTA name must equal the table's `refname`** (first header field). A
  mismatch means the loci won't match and the report will show everything
  not-present. If the contig names have extra prefixes/suffixes, fix the join key.
- **First base is the indel-sensitive edge.** Interior / trailing bases of a
  run are typically `eq/depth ≈ 1.0000` and never discriminate; if the pooled
  number looks suspiciously high you may have selected the wrong base — check
  `--which`.
- **Depth-0 rows are real and expected** on low-coverage contigs (e.g. the
  short `contig_4`/`contig_5` edges). They are reported in their own section and
  excluded from the accuracy distribution (dividing by 0).
- **Idempotent.** Re-running with the same `--outdir` overwrites the two output
  files. The table is read once in a single streaming pass (low memory).
- **Refresh-friendly.** When the CSV is re-generated, just re-run the same
  command — the candidate set is recomputed from the reference each time, and
  the per-loci values track the new table.
