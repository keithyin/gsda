---
name: smc-homopolymer-firstbase
description: Report one calling model's homopolymer (polyC/G, repeat >= N) first-base (or last / all) boundary accuracy — locus_accuracy and locus_accuracy_by_depth with built-in input verification. Input is either a pre-built per-locus table (--table) OR a raw query (unmapped BAM or FASTQ) + a single reference (--query + --ref), which the skill aligns via gsmm2+gsetl first.
---

# SMC — homopolymer first-base accuracy (single side)

## Purpose

This skill is the **non-A/B** counterpart of the sibling skill
`smc_ab_homopolymer_firstbase`. That one reads the two-sided `joined_all.csv`
and compares `control` vs `experiment` at homopolymer boundaries. **This one
measures a single side** — one set of reads vs one reference — and reports how
accurate it is at the homopolymer run boundaries, with **no second side and no Δ**.

It answers:

> **At the homopolymer boundaries — the base positions where indel error
> actually concentrates — how well does *this* calling model / read set
> resolve the repeat length?**

Homopolymer *run interiors* (e.g. the 4th of five G's in the plasmid backbone
`…GGGG[A]TCCTCTAGAG`) are usually near-perfect and carry no information. The
**first base of each run** (the run's 5′/start edge) is the indel-sensitive
position: the model must decide whether the consensus has the right *number* of
repeats, and that decision is scored at the first base. This skill isolates
exactly those positions.

It accepts the per-locus table **in either of two ways**, then does the same
analysis on both:

- **Mode 1 — from a ready table** (`--table`): pass a per-locus accuracy table
  that already exists (from the sibling
  `smc_ab_test_with_barcode_ref_locus` skill: `control_all.csv`,
  `experiment_all.csv`, or a per-barcode `locus_accuracy.csv`). Pure post-process,
  no gsmm2/gsetl.
- **Mode 2 — from a raw query** (`--query` + `--ref`): pass a **query file
  (an unmapped BAM, or a FASTQ)** and **one single-record reference FASTA**.
  Every read in the query is treated as a sequencing of *that one reference*.
  The skill gsmm2-aligns the query to the reference, runs `gsetl aligned-bam`
  to build the per-locus eq/diff/ins/del/depth table in memory, then runs the
  same homopolymer analysis. **No `--table`, `--mapping`, `--run`, or
  `--ref-dir` needed** — there is one reference and no barcode split.

For either mode, the script:
1. loads the Sanger reference `STR<plasmid>.fa` and finds every C/G homopolymer
   run of length `>= --min-run` (run detection on the raw reference sequence),
2. for each run selects the requested base — `--which first` (default, the
   run's first base), `last`, or `all` (every base in the run),
3. pulls the matching `(label, pos)` rows from the table,
4. reports per-position and depth-weighted-pooled `locus_accuracy` **and**
   `locus_accuracy_by_depth`.

## When to use

- **Mode 2 (query + ref)**: you have **reads** (an unmapped BAM or a FASTQ)
  that were all generated against **one reference**, and you want the
  **repeat-boundary** (first-base indel) accuracy of those reads against it —
  no pre-existing table. This is the direct "these reads vs this reference"
  case.
- **Mode 1 (table)**: you already have a **per-locus accuracy table** and want
  to re-read its **repeat-boundary** slice (or one side of an A/B run).
- you want a repeat-length / homopolymer-specific accuracy (first-base indel
  accuracy) rather than whole-plasmid or per-locus identity
- you want both metrics side by side — `eq/(eq+diff+ins+del)` (strict) and
  `eq/depth` (lenient) — to see whether the number is stable across the two

## Required inputs

### Mode 2 — query + reference (aligns on the fly)

- `--query` — the query file. Either an **unmapped BAM** (not an aligned one —
  gsmm2 does the alignment to `--ref` itself) or a **FASTQ**. All its reads are
  sequenced against the single `--ref`.
- `--ref` — the **single-record reference FASTA** every query read is aligned
  to. The plasmid name is recovered from the filename (leading `STR` + ext
  stripped, e.g. `STR3N-1.fa` → `3N-1`) and is only used to label the output.
- `--outdir` — where the three output files (and `aligned/` intermediates) go.

### Mode 1 — pre-built table

- `--table` — the tab-separated single-side per-locus table. Columns:
  `barcode refname pos eq diff ins del depth aroundBases locus_accuracy
  locus_accuracy_by_depth`. This is exactly the format of
  `control_all.csv` / `experiment_all.csv` (with a leading `barcode` column) or
  a per-barcode `control_locus_accuracy.csv`.
- `--ref-dir` — directory of single-record references `STR<plasmid>.fa`.
- `--mapping` — `plasmid<TAB>barcode<TAB>RUN号` TSV with a header (same file
  as the ref-locus skill).
- `--run RUN号` — which run's barcodes to use.
- `--outdir` — where the three output files are written.

The script errors if a required input for the detected mode is missing (and
tells you to use `--query + --ref` if you only gave `--table` args).

## Optional inputs

- `--min-run` — minimum homopolymer run length. **Default `4`.**
- `--which` — `first` (default) | `last` | `all`: which base(s) of each run to
  measure.
- `--char` — the characters that count as the homopolymer. **Default `CG`**
  (polyC/G). Pass e.g. `ACGT` for all homopolymers, or `A` for polyA only.

Mode 2 only:

- `--query-label` — the label to stamp on the query's rows in the output
  (default: the query filename without its extension, e.g. `reads.bam` →
  `reads`).
- `--rq-thr` — read-accuracy lower bound forwarded to gsmm2. **Default `0`
  (non-filtering).** Only bites if the query BAM actually carries an `rq` field;
  a FASTQ (or a plain unmapped BAM) has none, so gsmm2 ignores it.
- `--np-thr` — number-of-passes lower bound forwarded to gsmm2. **Default `0`
  (non-filtering).** Same story as `--rq-thr`.
- `--threads` — gsmm2 thread count (default: CPU count).

## Execution

```bash
# MODE 2 — query (unmapped bam or fastq) + one reference
/root/miniconda3/envs/gseda/bin/python \
    /root/projects/gsda/.claude/skills/smc_homopolymer_firstbase/scripts/homopolymer_firstbase.py \
    --query    .../reads.bam \
    --ref      .../merged_output/STR3N-1.fa \
    --outdir   .../homopolymer_firstbase \
    --min-run 4 --which first --char CG

# MODE 1 — pre-built per-locus table
/root/miniconda3/envs/gseda/bin/python \
    /root/projects/gsda/.claude/skills/smc_homopolymer_firstbase/scripts/homopolymer_firstbase.py \
    --table    .../locus_error_rate/control_all.csv \
    --ref-dir  .../merged_output \
    --mapping  .../plasmid_2_barcode.tsv \
    --run      20260805_250804Y0004_Run0001 \
    --outdir   .../locus_error_rate/homopolymer_firstbase \
    --min-run 4 --which first --char CG
```

> **Use the `gseda` env python** (`/root/miniconda3/envs/gseda/bin/python`) —
> `server/bin` shadows `$PATH` on this machine so `conda activate`/`conda run`
> silently keep the py3.8 `server` python. `gsmm2` and `gsetl` are on `PATH`
> (`/usr/bin/`); the script adds the `gseda` `src/` dir to `sys.path` itself.

- **Mode 1** is a table-only step — it reads the input table once and is
  **instant** (no gsmm2/gsetl). No backgrounding needed.
- **Mode 2** runs gsmm2 align + gsetl. Cost is a few seconds to ~1 min for a
  single query file. Run it as a background task and poll `--outdir` for
  `homopolymer_firstbase_report.md` if the read count is large.

## Output (under `--outdir`)

- `homopolymer_firstbase_locus.tsv` — the selected positions, all input columns
  plus two added: `base` (the reference base at the target, 0-based) and
  `runlen` (the run length it belongs to).
- `homopolymer_firstbase_summary.tsv` — one row per metric: `metric`,
  `n_positions`, `value` (the depth-weighted pooled number for that side).
- `homopolymer_firstbase_report.md` — human-readable report: the two pooled
  metrics, then a per-position table showing `aroundBases` (the `[ ]` bracket
  marks the target base) and both `locus_accuracy` and
  `locus_accuracy_by_depth`.
- (mode 2 only) `aligned/` — the `query.aligned.bam` and `query-gsetl/`
  intermediates from the alignment + gsetl step.

## The metrics (both reported)

- **`locus_accuracy = eq/(eq+diff+ins+del)`** — the strict form; the denominator
  counts only bases the gsetl table classified. Pooled = `Σeq/Σ(eq+diff+ins+del)`
  over the selected (label × position) rows (depth-weighted, **not** an
  arithmetic mean).
- **`locus_accuracy_by_depth = eq/depth`** — the lenient form; it also credits
  observed depth that was *not* counted as eq/diff/ins/del (e.g. discarded or
  unclassified depth) as read-correct. Pooled = `Σeq/Σdepth`. Because
  `eq ≤ depth` always, `by_depth ≥ locus_accuracy` at every position; the two
  differ only where some depth fell into none of the four buckets.

## Built-in verification (run before reporting — do not skip)

Before it prints anything, the script **self-checks the per-locus table** and
exits nonzero if either invariant fails:

- **V1 — indexing.** The base bracketed `[..]` in `aroundBases` must equal the
  reference base at the **0-based** `pos` column, for every row. This is the
  guard against the classic off-by-one: gsetl's `pos` is **0-based** (pos 0 =
  first base). If you (or an earlier step) joined a 1-based coordinate into the
  0-based `pos` column, V1 fails immediately instead of silently showing the
  *second* base of each run.
- **V2 — metric integrity.** `locus_accuracy_by_depth` must equal `eq/depth`
  (to 1e-6) and `eq ≤ depth` must hold for every row.

The script prints `[verify] OK: N/N rows pass V1 and V2` on success. If it
prints `INPUT VERIFICATION FAILED …` and exits nonzero, **stop and fix the
input** — do not trust a table produced with a shifted `pos`.

## Presenting the final result

When reporting to the user:

- State **which mode** was used (a pre-built `--table`, or `--query` vs `--ref`)
  and **which base** was measured (`--which`) and the run filter (`--char`,
  `--min-run`) — e.g. "first base of each polyC/G run ≥ 4".
- Lead with the **two pooled metrics** from `homopolymer_firstbase_summary.tsv`
  — give `locus_accuracy` and `locus_accuracy_by_depth` separately; they can
  differ in magnitude (by_depth is systematically higher) but should agree.
- Give the **per-position** table from the report, with `aroundBases` so the
  target base is visible (the bracketed letter is the target). Note the
  positions with the lowest accuracy.
- Say explicitly that this is **alignment identity vs the Sanger reference** at
  repeat boundaries (mode 2: FASTQ / plain unmapped BAM → no rq/np filtering),
  measuring how well the reads/model resolve the *repeat length* at the run's
  edge, not raw per-base quality.
- The **absolute path** of `homopolymer_firstbase_report.md` (and
  `homopolymer_firstbase_locus.tsv`) must appear in the final message.

### Final message requirement

Include the **absolute path** of `homopolymer_firstbase_report.md`, the
`n_positions` count, and the pooled value of **both** metrics.

## Notes / gotchas

- **`pos` is 0-based** in the gsetl `fact_aligned_bam_ref_locus_info` table
  (pos 0 = first reference base). Run boundaries and any position math must use
  0-based indices. V1 enforces this — if a table was ever produced with 1-based
  positions, this skill refuses it.
- **First base ≈ the indel-sensitive position.** Interior and trailing bases of
  a run (esp. the backbone polyG `…GGGG[A]T…`) are typically `locus_accuracy ≈
  1.0000` and never discriminate. If your pooled number looks suspiciously high
  (~0.997–1.0) you may have selected interior bases rather than the first —
  check `--which`.
- **by_depth ≥ locus_accuracy always.** If you see by_depth *lower* than
  locus_accuracy somewhere, V2 should have already failed — treat the table as
  untrustworthy.
- **Unmapped / zero-depth positions.** A target `(label, pos)` with no row in
  the table (e.g. a 0-depth run edge) is skipped with a `[warn]` line and
  counted — the report's `n_positions` reflects what was actually matched.
- **Mode 2 — query must be UNMAPPED if it's a BAM.** gsmm2 aligns the query to
  `--ref` itself; feed it the raw unmapped BAM (as the SMC all-reads BAM), not
  a BAM already aligned to something else. A FASTQ works the same way.
- **Mode 2 — one reference per run.** All reads in the query are aligned to the
  single `--ref`. If your data spans several plasmids, run the skill once per
  (query, reference) pair — there is no mapping/demux step here.
- **Mode 2 — `gsetl` prints `fetch error target_name:…-fwd/_rev`** — these are
  **benign** (fwd/rev sub-region lookups); the locus table still populates
  fully. Do not treat them as failures.
- **Mode 2 — very short references give depth 0.** gsmm2 may refuse to align
  reads against a very short target, leaving the locus table all zero-depth (the
  pooled metric then reads `nan`). This is an alignment artifact of the tiny
  reference, not a skill bug.
- **Mode 1 — same run, same mapping** as the ref-locus skill — the barcode set
  is taken from the rows of `--mapping` whose `RUN号` equals `--run`, not from
  whatever FASTQs happen to exist. The input table must have a `barcode` column
  (`BarcodeNN`) for the mapping to match it.
- **Idempotent.** Re-running with the same `--outdir` just overwrites the three
  output files (mode 2 also re-runs the alignment into `aligned/`).
