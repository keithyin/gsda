---
name: smc-ab-homopolymer-firstbase
description: Post-process a locus_error_rate joined_all table to A/B compare the two calling models at homopolymer (polyC/G, repeat >= N) boundaries — the first (or last / all) base of each run — reporting both locus_accuracy and locus_accuracy_by_depth with built-in input verification.
---

# SMC A/B Test — homopolymer first-base accuracy

## Purpose

This skill is the **homopolymer-focused post-processing** step on top of the
per-reference-locus A/B that the sibling skill
`smc_ab_test_with_barcode_ref_locus` produces. That skill gives you a whole
`joined_all.csv` (every reference base position, both sides, ctrl/exp counts);
this one answers the narrower question:

> **At the homopolymer boundaries — the base positions where indel error
> actually concentrates — how do the two calling models differ?**

Homopolymer *run interiors* (e.g. the 4th of five G's in the plasmid backbone
`…GGGG[A]TCCTCTAGAG`) are usually near-perfect on both models and carry no
discrimination. The **first base of each run** (the run's 5′/start edge) is the
indel-sensitive position: the model must decide whether the consensus has the
right *number* of repeats, and that decision is scored at the first base. This
skill isolates exactly those positions and compares `control` vs `experiment`
on both accuracy metrics.

For every reference in the run's mapping it:
1. loads the Sanger reference `STR<plasmid>.fa` and finds every C/G homopolymer
   run of length `>= --min-run` (run detection on the raw reference sequence),
2. for each run selects the requested base — `--which first` (default, the
   run's first base), `last`, or `all` (every base in the run),
3. pulls the matching `(barcode, pos)` rows from `joined_all.csv`,
4. reports per-position and depth-weighted-pooled `locus_accuracy` **and**
   `locus_accuracy_by_depth` for both sides.

## When to use

- the overall / per-locus A/B (sibling skill) is done, and you want to know
  **where the difference is** specifically at **repeat boundaries**
- you want a repeat-length / homopolymer-specific readout (first-base
  indel accuracy) rather than whole-plasmid or per-locus identity
- you want both metrics side by side — `eq/(eq+diff+ins+del)` (strict, the
  "counted-outcome" form) and `eq/depth` (lenient, treats uncounted depth as
  read-correct) — to see whether the difference is real or a denominator artifact

> This skill **does not** re-run gsmm2/gsetl and does not re-demux. It reads
> the `joined_all.csv` that `smc_ab_test_with_barcode_ref_locus` already wrote.
> Run that skill first if you don't have `joined_all.csv` yet.

## Required inputs

- `--joined` — the tab-separated `joined_all.csv` from the ref-locus skill
  (columns `barcode refname pos eq_ctrl diff_ctrl ins_ctrl del_ctrl depth_ctrl
  aroundBases_ctrl locus_accuracy_ctrl locus_accuracy_by_depth_ctrl …_exp`).
- `--ref-dir` — directory of single-record references `STR<plasmid>.fa`.
- `--mapping` — `plasmid<TAB>barcode<TAB>RUN号` TSV with a header (same file
  as the ref-locus skill).
- `--run RUN号` — which run's rows to use.
- `--outdir` — where the three output files are written.

## Optional inputs

- `--min-run` — minimum homopolymer run length. **Default `4`.**
- `--which` — `first` (default) | `last` | `all`: which base(s) of each run to
  compare.
- `--char` — the characters that count as the homopolymer. **Default `CG`**
  (polyC/G). Pass e.g. `ACGT` for all homopolymers, or `A` for polyA only.

## Execution

```bash
/root/miniconda3/envs/gseda/bin/python \
    /root/projects/gsda/.claude/skills/smc_ab_homopolymer_firstbase/scripts/homopolymer_firstbase_ab.py \
    --joined    .../locus_error_rate/joined_all.csv \
    --ref-dir   .../merged_output \
    --mapping   .../plasmid_2_barcode.tsv \
    --run       20260805_250804Y0004_Run0001 \
    --outdir    .../locus_error_rate/homopolymer_firstbase \
    --min-run 4 --which first --char CG
```

> **Use the `gseda` env python** (`/root/miniconda3/envs/gseda/bin/python`) —
> `server/bin` shadows `$PATH` on this machine so `conda activate`/`conda run`
> silently keep the py3.8 `server` python. This script is pure stdlib (csv/os/sys)
> and only needs any 3.8+ interpreter, but use the gseda one for consistency.

This is a table-only step — it reads the joined table once and is **instant**
(no gsmm2/gsetl). No backgrounding needed.

## Output (under `--outdir`)

- `homopolymer_firstbase_locus.tsv` — the selected positions, all `joined_all.csv`
  columns plus two added: `base` (the reference base at the target, 0-based) and
  `runlen` (the run length it belongs to).
- `homopolymer_firstbase_summary.tsv` — one row per metric: `metric`,
  `n_positions`, `control`, `experiment`, `delta_exp_minus_ctrl` (depth-weighted
  pooled number for that side).
- `homopolymer_firstbase_report.md` — human-readable report: the two pooled
  metrics, then a per-position table showing `aroundBases` (the `[ ]` bracket
  marks the target base), both `locus_accuracy` and `locus_accuracy_by_depth`
  for each side, and Δ.

## The metrics (both reported)

- **`locus_accuracy = eq/(eq+diff+ins+del)`** — the strict form; the denominator
  counts only bases the gsetl table classified. Pooled = `Σeq/Σ(eq+diff+ins+del)`
  over the selected (barcode × position) rows (depth-weighted, **not** an
  arithmetic mean).
- **`locus_accuracy_by_depth = eq/depth`** — the lenient form; it also credits
  observed depth that was *not* counted as eq/diff/ins/del (e.g. discarded or
  unclassified depth) as read-correct. Pooled = `Σeq/Σdepth`. Because
  `eq ≤ depth` always, `by_depth ≥ locus_accuracy` at every position; the two
  differ only where some depth fell into none of the four buckets.

## Built-in verification (run before reporting — do not skip)

Before it prints anything, the script **self-checks the input** and exits
nonzero if either invariant fails:

- **V1 — indexing.** The base bracketed `[..]` in `aroundBases_exp` must equal
  the reference base at the **0-based** `pos` column, for every row. This is the
  guard against the classic off-by-one: gsetl's `pos` is **0-based** (pos 0 =
  first base). If you (or an earlier step) joined a 1-based coordinate into the
  0-based `pos` column, V1 fails immediately instead of silently showing the
  *second* base of each run.
- **V2 — metric integrity.** `locus_accuracy_by_depth` must equal `eq/depth`
  (to 1e-6) and `eq ≤ depth` must hold, for both sides, every row.

The script prints `[verify] OK: N/N rows pass V1 and V2` on success. If it
prints `INPUT VERIFICATION FAILED …` and exits nonzero, **stop and fix the
input** — do not trust a table produced with a shifted `pos`.

## Presenting the final result

When reporting to the user:

- State **which base** was compared (`--which`) and the run filter
  (`--char`, `--min-run`) — e.g. "first base of each polyC/G run ≥ 4".
- Lead with the **two pooled metrics** (both sides, Δ) from
  `homopolymer_firstbase_summary.tsv` — give `locus_accuracy` and
  `locus_accuracy_by_depth` separately; they can differ in magnitude (by_depth
  is systematically higher) but should agree in **direction**.
- Give the **per-position** table from the report, with `aroundBases` so the
  target base is visible (the bracketed letter is the target). Note the
  positions that moved most, and in which direction.
- Say explicitly that this is **alignment identity vs the Sanger reference** at
  repeat boundaries (FASTQ input → no rq/np filtering), measuring how well each
  model resolves the *repeat length* at the run's edge, not raw per-base quality.
- The **absolute path** of `homopolymer_firstbase_report.md` (and
  `homopolymer_firstbase_locus.tsv`) must appear in the final message.

### Final message requirement

Include the **absolute path** of `homopolymer_firstbase_report.md`, the
`n_positions` count, and the pooled number for **both** metrics (ctrl, exp, Δ)
for each.

## Notes / gotchas

- **`pos` is 0-based** in the gsetl `fact_aligned_bam_ref_locus_info` table
  (pos 0 = first reference base). Run boundaries and any position math must use
  0-based indices. V1 enforces this — if a table was ever produced with 1-based
  positions, this skill refuses it.
- **First base ≈ the indel-sensitive position.** Interior and trailing bases of
  a run (esp. the backbone polyG `…GGGG[A]T…`) are typically `locus_accuracy ≈
  1.0000` on both sides and never discriminate. If your pooled number looks
  suspiciously high (~0.997–1.0) you may have selected interior bases rather
  than the first — check `--which`.
- **by_depth ≥ locus_accuracy always.** If you see by_depth *lower* than
  locus_accuracy somewhere, V2 should have already failed — treat the table as
  untrustworthy.
- **Direction can shift between the two metrics.** A position can look negative
  under `eq/(eq+diff+ins+del)` (denominator includes classified diffs) but
  neutral under `eq/depth` (if both sides' eq and depth are equal and only the
  classified diff/ins/del counts differ). Report the metric you judge by and note
  the other.
- **Unmapped / zero-depth positions.** A target `(barcode, pos)` with no row in
  the joined table (e.g. a 0-depth run edge) is skipped with a `[warn]` line and
  counted — the report's `n_positions` reflects what was actually matched.
- **Same run, same mapping** as the ref-locus skill — only rows whose `RUN号`
  equals `--run` are used; the barcode set is taken from the mapping, not from
  whatever FASTQs happen to exist.
- **Idempotent.** Re-running with the same `--outdir` just overwrites the three
  output files.
