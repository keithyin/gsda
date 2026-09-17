---
name: smc-ab-test-ref-locus-single
description: Single-reference A/B compare two call outputs' per-reference-locus accuracy (eq/diff/ins/del per locus + locus_accuracy) over one pooled query per side, via gsmm2 + gsetl through gseda.ab_analysis.locus_error_rate.
---

# SMC A/B Test — single-reference per-locus accuracy

## Purpose

Compare **two** call outputs that differ only in how the reads were called
(e.g. a new calling model vs. the baseline over the same SMC run) at a
**per-reference-locus** grain, for the case where **the whole run maps to ONE
reference** — instead of the usual "each barcode → its own plasmid reference".
Each side is a **single query** (one FASTQ/BAM, or a demux dir whose
`*.fastq*` are pooled into one FASTQ), aligned to the **same** reference, so the
comparison isolates the calling model.

For each side the script:
1. resolves the query — a file is used as-is, a directory is pooled
   (sorted `*.fastq*` concatenated) into one FASTQ,
2. gsmm2-aligns it to the single `--ref` Sanger reference,
3. `gsetl aligned-bam` emits a per-locus table (eq/diff/ins/del/depth per
   reference position),
4. a per-locus **`locus_accuracy` = eq/(eq+diff+ins+del)** is added.

The two sides' per-locus tables are joined on `(refname, pos)` and a
depth-weighted accuracy summary (overall Δ) is written.

## When to use

- **one run → one reference**: a whole SMC run's reads (or one sample's reads)
  compared against a single Sanger reference, at per-locus resolution — you want
  to see *where* on the reference one model errs more (which regions, which
  homopolymer stretches).
- comparing two calling models / two demux passes at locus resolution, with a
  single shared reference (no per-barcode / per-plasmid mapping).
- you already have the reads as one FASTQ/BAM, or as a demux dir to pool.

> **Pick the sibling skill by the reference layout:**
> - **this skill** (`smc_ab_test_ref_locus_single`) — one run, **one** reference,
>   pooled query per side, A/B.
> - `smc_ab_test_with_barcode_ref_locus` — one run, **many barcodes each with its
>   own plasmid reference**; per-barcode + pooled per-locus A/B (needs a
>   `plasmid↔barcode↔RUN` mapping + a `--ref-dir`).
> - `smc_ab_test_with_barcode` — same per-barcode A/B but at **whole-sample
>   identity / mismatch / indel / `identity≥0.99`** grain via
>   `sequencing_report_v2` (not per-locus).

## Required inputs

- `--control LABEL=QUERY` and `--exp LABEL=QUERY` — the two query sources. `QUERY`
  is either a **file** (`.fastq`/`.fq`/`.fastq.gz`, or an **unaligned/unmapped**
  `.bam`) used as-is, **or a directory** whose `*.fastq*` are concatenated (sorted
  by filename) into one pooled FASTQ. The label becomes the summary header (e.g.
  `baseline`, `new_model`) — use labels without spaces/`=`. `--control` is the
  "ctrl" side and `--exp` is the "exp" side; Δ is reported as exp − ctrl.
  - **BAM queries must be fully unmapped.** gsmm2 does the alignment to `--ref`,
    so a BAM side is expected to carry *unaligned* reads — **every** read flagged
    UNMAP (`0x4`). The script checks this up front and **exits with an error** if
    any read is already mapped, rather than silently re-aligning an already-aligned
    BAM and producing a bogus identity. If you have an aligned BAM, supply the
    unaligned source (or a FASTQ) instead.
- `--ref` — the single reference FASTA (the run's Sanger reference, e.g.
  `STR3N-1.fa`).
- `--outdir` — where the pooled FASTQs, per-locus tables, aligned bams, and
  `summary.md` go.

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
    /root/projects/gsda/.claude/skills/smc_ab_test_ref_locus_single/scripts/ref_locus_single.py \
    --control baseline=.../called-barcode-v4-baseline/demuxed \
    --exp     new_model=.../called-barcode-v4-2026Q2Model/demuxed \  
    --ref     .../merged_output/STR3N-1.fa \
    --outdir  .../locus_error_rate_single
```

> **Use the `gseda` env python, not `conda run -n ...`.** On this machine
> `server/bin` shadows `$PATH` so `conda activate`/`conda run` silently keep the
> py3.8 `server` python — call the absolute interpreter
> `/root/miniconda3/envs/gseda/bin/python`. `gsmm2` and `gsetl` are on `PATH`
> (`/usr/bin/`). The script adds the `gseda` `src/` dir to `sys.path` itself.

Cost is a few seconds per side (gsmm2 align + gsetl locus table) — fast. Run it
as a background task and poll the outdir for `summary.md`.

## Output (under `--outdir`)

- `.pool_<label>.fastq` — the pooled FASTQ for each side that was a directory
  (a file query is used as-is, no pool file).
- `control_locus_accuracy.csv`, `experiment_locus_accuracy.csv` — columns
  `refname pos eq diff ins del depth aroundBases locus_accuracy
  locus_accuracy_by_depth group` (tab-separated).
- `locus_accuracy_joined.csv` — per-locus ctrl/exp columns with `_ctrl`/`_exp`
  suffixes, keyed on `refname`,`pos` (tab-separated).
- `summary.md` — depth-weighted A/B summary (overall both sides + Δ).
- `control.aligned.bam` / `experiment.aligned.bam` + `*-gsetl/` dirs.

## The accuracy metric

`locus_accuracy = eq/(eq+diff+ins+del)` is computed **per reference base
position** (each row of the gsetl `fact_aligned_bam_ref_locus_info` table).
The number in `summary.md` is a **depth-weighted** sum: `Σeq /
Σ(eq+diff+ins+del)` over all (refname, pos) rows for that side — **not** an
arithmetic mean of per-locus accuracies.

## Presenting the final result

When reporting to the user:

- Lead with the **overall** depth-weighted `locus_accuracy` for both sides and
  Δ (exp − ctrl) from `summary.md`.
- State explicitly that this is **alignment identity vs the Sanger reference**
  (FASTQ input → the `rq`/`np` thresholds are inert), so it measures how close
  each model's consensus is to the first-gen reference at each base, not raw
  per-base read quality.
- For *where* the difference is, point at `locus_accuracy_joined.csv`: filter to
  rows with `locus_accuracy_ctrl != locus_accuracy_exp` (or a large per-locus gap)
  to see the specific positions.

### Final message requirement

The message sent back to the user **must** include the **absolute path** of the
final `locus_accuracy_joined.csv` (e.g. `ls -d $OUT/locus_accuracy_joined.csv`),
in addition to the overall depth-weighted `locus_accuracy` (both sides) and Δ
(exp − ctrl).

## Notes / gotchas

- **BAM query must be unmapped (checked).** A BAM side is only valid if **all**
  its reads are flagged UNMAP (`0x4`) — the query feeds gsmm2, which performs the
  alignment to `--ref`. The script walks the BAM and **hard-exits** on the first
  mapped read (it will not re-align an already-aligned BAM). Feed the unaligned
  source BAM or a FASTQ; do not pass an aligned BAM.
- **FASTQ input → `--rq-thr` / `--np-thr` are ignored.** gsmm2 only applies
  `--rq-range`/`--np-range` when the query BAM actually carries those fields.
  These FASTQs do not, so the metric is plain alignment identity. (If you ever
  feed an rq-tagged BAM, the thresholds would then filter reads.)
- **Pooling order.** A directory query concatenates `*.fastq*` in **sorted
  filename** order (so `Barcode01..NN` pool in a stable, barcode-preserving
  order); read headers keep their original names. Pooling is read-order
  independent for the per-locus counts, so order only affects determinism, not
  the numbers. If you need a per-barcode breakdown instead of pooled, use the
  `smc_ab_test_with_barcode_ref_locus` sibling.
- **`gsetl` prints `fetch error target_name:STRxxx-fwd/_rev, err=failed to
  fetch region`** — these are **benign** (fwd/rev sub-region lookups); the
  locus table still populates fully. Do not treat them as failures.
- **Short / incomplete references.** If the first-gen reference is shorter than
  the target (e.g. some STR33 series), low `locus_accuracy`/coverage at the
  uncovered positions is an **incomplete reference**, not a model problem —
  note it when judging the A/B.
- **Idempotency / rerun.** `locus_error_rate` deletes and regenerates its
  `<side>.aligned.bam` each call, the pooled FASTQs are rewritten each run, and
  the per-locus tables + `summary.md` are regenerated from the joined table
  every run, so re-running on the same `--outdir` is safe and idempotent.
- **`locus_accuracy_by_depth`** is also emitted per locus (`eq/depth`); it differs
  from `locus_accuracy` at positions where `depth` (observed) exceeds the
  eq+diff+ins+del count. The summary uses `locus_accuracy` (the
  eq/(eq+diff+ins+del) form).
