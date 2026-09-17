---
name: smc-ab-test-with-barcode-ref-locus
description: Per-barcode A/B compare two barcode-call outputs' per-reference-locus accuracy (eq/diff/ins/del per locus + locus_accuracy), via gsmm2 + gsetl through gseda.ab_analysis.locus_error_rate.
---

# SMC A/B Test with Barcode — per-reference-locus accuracy

## Purpose

Compare **two** barcode-split outputs that differ only in how the barcodes
were called (e.g. a new calling model vs. the baseline over the same SMC run),
but at a **finer grain** than whole-sample identity: accuracy is computed
**per reference locus** (each base position of the plasmid's Sanger reference).
Each side is a directory of `BarcodeNN.fastq` demuxed from the **same** SMC BAM
with a different calling model, so the comparison isolates the calling model.

For every barcode the plasmid→barcode mapping covers:
1. the control side's `BarcodeNN.fastq` is gsmm2-aligned to the plasmid's
   Sanger reference `STR<plasmid>.fa`,
2. the experiment side's `BarcodeNN.fastq` is aligned the same way,
3. `gsetl aligned-bam` emits a per-locus table (eq/diff/ins/del/depth per
   reference position),
4. a per-locus **`locus_accuracy` = eq/(eq+diff+ins+del)** is added.

The two sides' per-locus tables are joined on `(refname, pos)`, then aggregated
across all barcodes into per-group and joined tables, and a depth-weighted
accuracy summary is written.

## When to use

- you need **per-reference-locus** (base-position) accuracy, not a single
  whole-sample identity — e.g. to see *where* on the plasmid one model errs
  more (which regions, which homopolymer stretches)
- comparing two calling models / two demux passes at locus resolution
- producing a per-barcode + pooled per-locus A/B comparison

> For a single whole-sample identity / mismatch / indel / `identity≥0.99`
> comparison instead, use the sibling skill `smc_ab_test_with_barcode`.

## Required inputs

- `--control LABEL=DIR` and `--exp LABEL=DIR` — the two demuxed dirs, each
  containing `BarcodeNN.fastq`. The label becomes the column header (e.g.
  `baseline`, `new_model`). Use labels without spaces/`=`. `--control` is the
  "ctrl" side and `--exp` is the "exp" side; Δ is reported as exp − ctrl.
- `--mapping` — `plasmid<TAB>barcode<TAB>RUN号` TSV with a header; barcode is
  written like `24标签-N`.
- `--run RUN号` — which run's rows of the mapping to use (the file holds
  several runs; only the matching rows are selected).
- `--ref-dir` — directory of single-record references named `STR<plasmid>.fa`
  (the plasmid id in the mapping, e.g. `3N-1` → `STR3N-1.fa`).
- `--outdir` — where per-barcode dirs, aggregate tables, and `summary.md` go.

## Optional inputs

- `--rq-thr` — read-accuracy threshold forwarded to
  `locus_error_rate.main_cli` (1 value = shared, or 2 = `[ctrl,exp]`).
  **Default `0` (non-filtering).** Because the input is FASTQ (no CCS `rq`
  field), gsmm2 **ignores** it; it is still a required argument of the
  underlying script so a value must be supplied.
- `--np-thr` — number-of-passes lower bound (1 or 2 values). **Default `5`.**
  Ignored for FASTQ (no `np` field).
- `--threads` — gsmm2 thread count (default: CPU count).

## Execution

```bash
/root/miniconda3/envs/gseda/bin/python \
    /root/projects/gsda/.claude/skills/smc_ab_test_with_barcode_ref_locus/scripts/locus_error_rate_ab.py \
    --control baseline=.../called-barcode-v4-baseline/demuxed \
    --exp     new_model=.../called-barcode-v4-2026Q2Model/demuxed \
    --mapping .../plasmid_2_barcode.tsv \
    --run     20260805_250804Y0004_Run0001 \
    --ref-dir .../merged_output \
    --outdir  .../locus_error_rate
```

> **Use the `gseda` env python, not `conda run -n ...`.** On this machine
> `server/bin` shadows `$PATH` so `conda activate`/`conda run` silently keep the
> py3.8 `server` python — call the absolute interpreter
> `/root/miniconda3/envs/gseda/bin/python`. `gsmm2` and `gsetl` are on `PATH`
> (`/usr/bin/`). The script adds the `gseda` `src/` dir to `sys.path` itself.

Per-barcode cost is a few seconds per side (gsmm2 align + gsetl locus table).
For ~15 barcodes × 2 sides that is ~1–2 min. Run it as a background task and
poll the outdir for `summary.md`.

## Output (under `--outdir`)

- `Barcode01/` … `BarcodeNN/` — per barcode: `control_locus_accuracy.csv`,
  `experiment_locus_accuracy.csv` (columns `refname pos eq diff ins del depth
  aroundBases locus_accuracy locus_accuracy_by_depth group`),
  `locus_accuracy_joined.csv`
  (per-locus ctrl/exp columns with `_ctrl`/`_exp` suffixes, keyed on
  `refname`,`pos`), plus the `control.aligned.bam` / `experiment.aligned.bam`
  and `*-gsetl/` dirs.
- `control_all.csv` / `experiment_all.csv` — all barcodes stacked, with a
  leading `barcode` column.
- `joined_all.csv` — all barcodes' joined per-locus tables stacked, leading
  `barcode` column.
- `summary.md` — depth-weighted A/B summary (overall + per-barcode).

## The accuracy metric

`locus_accuracy = eq/(eq+diff+ins+del)` is computed **per reference base
position** (each row of the gsetl `fact_aligned_bam_ref_locus_info` table).
The pooled / per-barcode number in `summary.md` is a **depth-weighted** sum:
`Σeq / Σ(eq+diff+ins+del)` over all (barcode × locus) rows for that side — **not**
an arithmetic mean of per-locus accuracies.

## Presenting the final result

When reporting to the user:

- Lead with the **overall** depth-weighted `locus_accuracy` for both sides and
  Δ (exp − ctrl) from `summary.md`.
- Give the **per-barcode** table (barcode, locus, ctrl, exp, Δ) — note which
  barcodes moved most and in which direction.
- State explicitly that this is **alignment identity vs the Sanger reference**
  (FASTQ input → the `rq`/`np` thresholds are inert), so it measures how close
  each model's consensus is to the first-gen reference at each base, not raw
  per-base read quality.
- For *where* the difference is, point at `joined_all.csv`: filter to rows with
  `locus_accuracy_ctrl != locus_accuracy_exp` (or a large per-locus gap) to see
  the specific positions.

### Final message requirement

The message sent back to the user **must** include the **absolute path** of the
final `joined_all.csv` (e.g. `ls -d $OUT/joined_all.csv`), in addition to the
overall depth-weighted `locus_accuracy` (both sides) and Δ (exp − ctrl).

## Notes / gotchas

- **FASTQ input → `--rq-thr` / `--np-thr` are ignored.** gsmm2 only applies
  `--rq-range`/`--np-range` when the query BAM actually carries those fields.
  These demuxed FASTQs do not, so the metric is plain alignment identity. (If
  you ever feed an rq-tagged BAM, the thresholds would then filter reads.)
- **`gsetl` prints `fetch error target_name:STRxxx-fwd/_rev, err=failed to
  fetch region`** — these are **benign** (fwd/rev sub-region lookups); the
  locus table still populates fully. Do not treat them as failures.
- **The two sides must share the same run.** Only rows of `--mapping` whose
  `RUN号` equals `--run` are used. A `BarcodeNN.fastq` present in a demux dir
  but not covered by that run's mapping (e.g. an unmapped extra barcode) is
  skipped with a log line.
- **Barcode naming must already be `BarcodeNN.fastq`.** This skill does not do
  the `24标签-N ↔ Adaptor-* ↔ barcode-NN` remapping. If the demux dirs are still
  in `Adaptor-barcodeNNN-M.fastq` form, build the `BarcodeNN.fastq` view first
  (see the `barcode24renaming` skill).
- **Short / incomplete references.** If the first-gen reference is shorter than
  the target (e.g. some STR33 series), low `locus_accuracy`/coverage at the
  uncovered positions is an **incomplete reference**, not a model problem —
  treat those barcodes separately when judging the A/B.
- **Idempotency / rerun.** `locus_error_rate` deletes and regenerates its
  `<side>.aligned.bam` each call, and the aggregate tables + `summary.md` are
  rewritten from the per-barcode CSVs every run, so re-running on the same
  `--outdir` is safe and idempotent.
- **`locus_accuracy_by_depth`** is also emitted per locus (`eq/depth`); it differs
  from `locus_accuracy` at positions where `depth` (observed) exceeds the
  eq+diff+ins+del count. The summary uses `locus_accuracy` (the
  eq/(eq+diff+ins+del) form).
