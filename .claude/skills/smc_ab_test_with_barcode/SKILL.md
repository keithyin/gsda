---
name: smc-ab-test-with-barcode
description: A/B compare two barcode-call outputs (e.g. two calling models over the same SMC BAM) per-barcode and pooled, via sequencing_report_v2.
---

# SMC A/B Test with Barcode

## Purpose

Compare **two** barcode-split outputs that differ only in how the barcodes were
called (e.g. a new barcode-calling model vs. the baseline), and quantify the
difference in sequencing-quality metrics. Each side is a directory of
`BarcodeNN.fastq` produced by demuxing the **same** SMC BAM with a different
barcode-calling model — so the comparison isolates the effect of the calling
model, not of the underlying reads.

For every barcode the plasmid→barcode mapping covers, per model:
1. convert `BarcodeNN.fastq` → unaligned BAM (with an `rq` tag) via `fastq2bam.py`,
2. run `gseda.ppl.sequencing_report_v2` against the plasmid's Sanger reference,
3. emit per-barcode aggregated metrics.

Then it compares the two models **per-barcode** and **pooled**.

## When to use

- comparing two barcode-calling models / two demux passes over the same run
- checking whether a new model changes identity / mismatch / indel / homopolymer
  error rates or the `identity≥0.99` fraction
- producing a per-barcode + aggregate A/B comparison table

## Required inputs

- `--a LABEL=DIR` and `--b LABEL=DIR` — the two demuxed dirs, each containing
  `BarcodeNN.fastq`. The label becomes the column header (e.g. `new_model`,
  `baseline`). Use labels without spaces/`=`.
- `--mapping` — `plasmid<TAB>barcode<TAB>RUN号` TSV with a header; barcode is
  written like `24标签-N`.
- `--run RUN号` — which run's rows of the mapping to use (the file holds several
  runs; only the matching rows are selected).
- `--ref-dir` — directory of single-record references named `STR<plasmid>.fa`
  (the plasmid id in the mapping, e.g. `3N-1` → `STR3N-1.fa`).
- `--outdir` — where BAMs, per-sample metric dirs, and the comparison tables go.

## Optional inputs

- `--rq-range` — per-read accuracy filter, e.g. `0.99:1.1`. Default **none**
  (all demuxed reads are reported).
- `--np-thr` — minimum `np` (passes) threshold; only channels with `np >= np-thr`
  are reported. **Default `5`.** It is forwarded to `sequencing_report_v2` as
  `--np-range <np-thr>:10000000`. Pass `0` to disable the filter (all channels).

## Execution

```bash
/root/miniconda3/envs/gseda/bin/python \
    /root/projects/gsda/.claude/skills/smc_ab_test_with_barcode/scripts/ab_test.py \
    --a new_model=.../called-barcode-v4-2026Q2Model/demuxed \
    --b baseline=.../called-barcode-v4-baseline/demuxed \
    --mapping  .../plasmid_2_barcode.tsv \
    --run      20260805_250804Y0004_Run0001 \
    --ref-dir  .../merged_output \
    --outdir   .../report_v2_new_model_vs_baseline
```

> **Use the `gseda` env python, not `conda run -n ...`.** On this machine
> `server/bin` shadows `$PATH` so `conda activate`/`conda run` silently keep the
> py3.8 `server` python — call the absolute interpreter
> `/root/miniconda3/envs/gseda/bin/python`. `gsmm2-aligned-metric` (mm2) is on
> PATH and satisfies `sequencing_report_v2`'s env check.

Per-barcode cost is a few seconds (fastq2bam) + a few seconds to ~1 min
(sequencing_report_v2) depending on read count. For ~15 barcodes × 2 models that
is ~5–10 min; run it as a background task and poll the outdir.

## Output (under `--outdir`)

- `<labelA>/`, `<labelB>/` — per-barcode `BarcodeNN_<plasmid>.bam` +
  `BarcodeNN_<plasmid>-metric/` (`..._aggr.csv`, `..._fact.csv`, `..._basic.csv`)
  and an `_agg/` dir.
- `per_barcode.tsv` — long format: `barcode  plasmid  metric  A  B  delta_A_minus_B`.
- `aggregate.csv` — pooled `metric  A  B  delta`.
- `aggregate_<labelA>.csv`, `aggregate_<labelB>.csv` — full per-model aggregates.
- `summary.md` — the pooled comparison table.

## How the "aggregate" (汇总) is computed

It is **not** an arithmetic mean of the per-barcode metrics. The script uses
`sequencing_report_v2.py`'s own `--fact-csvs`/`--basic-csvs` merge path
(`merge_partition` → `align_stats` + `bam_basic_ana`): it **concatenates all
per-barcode per-read fact rows into one dataset** and recomputes, so
`identity`, `queryCoverage`, `mmRate`, the indel/homopolymer rates and the
`identity≥0.x` fractions are **bases-/read-weighted global values**. Sanity
check: the aggregate `reads_num` must equal the sum of the per-barcode
`reads_num`.

## Key metrics

- `reads_num`, `tot_bases`, `n50`, `read_len_p50` — read stats (should be nearly
  identical across the two models; only demux filtering differs).
- `alignedRatio` / `notAlignedRatio` — fraction aligned (see gotcha below).
- `queryCoverage`, `queryCoverage2`, `queryCoverage3` — read coverage of ref.
- `identity`, `identity-p50`, `mmRate` — match/mismatch over aligned span.
- `NHInsRate`/`NHDelRate` (non-homopolymer) and `HomoInsRate`/`HomoDelRate`
  (homopolymer) — per-error-type rates.
- `longIndelRatio`, `GlobalQueryCoverage`.
- `identity≥0.83` / `≥0.90` / `≥0.99` / `≥0.999` — fraction of reads meeting the threshold.
  `identity≥0.99` is usually the most sensitive A/B signal.

## Presenting the final result

When reporting the A/B outcome to the user, **list every metric in the aggregate
(`aggr`) output — not just the `SEL` subset in `summary.md`**. Read
`aggregate_<labelA>.csv` and `aggregate_<labelB>.csv` and present all rows, not only
the pooled table in `summary.md` (which shows only the `SEL` list). Format:

- One table per metric group, columns `metric | <labelA> | <labelB> | Δ(A−B)`,
  where Δ = A − B (the first/labeled model minus the second). Group the rows as:
  1. **Read / base stats** — `reads_num`, `tot_bases`, `n50`, `read_len_avg`,
     `read_len_min/p25/p50/p75/max`, `≥Q8/≥Q10/≥Q15/≥Q20/≥Q30`, `4xQ20`.
  2. **Alignment / coverage** — `aligned`, `notAligned`, `alignedRatio`,
     `notAlignedRatio`, `alignedBases`, `queryCoverage`/`2`/`3`/`-p25/p50/p75`,
     `GlobalQueryCoverage`.
  3. **Match / error rates (aligned span)** — `identity` + `-p25/p50/p75`, `mmRate`,
     `NHInsRate`/`NHDelRate`, `HomoInsRate`/`HomoDelRate`, `longIndelRatio`.
  4. **Per-base identity / mmRate** — `identity-{A,C,G,T}` and
     `mmRate-{A,C,G,T}`, `NHInsRate-{A,C,G,T}`, `HomoInsRate-{A,C,G,T}`,
     `NHDelRate-{A,C,G,T}`, `HomoDelRate-{A,C,G,T}`.
  5. **Identity-threshold fractions** — `identity≥0.83 / ≥0.90 / ≥0.99 /
     ≥0.999 / ≥0.9999 / ≥0.99999`.
- Numeric formatting: ratios/rates to ~6 sig figs; `Δ` signed; count metrics
  (reads_num, base counts, ≥Q*) as plain integers.
- Also give a **per-barcode** breakdown of the key A/B metrics
  (`reads_num`, `identity`, `mmRate`, `identity≥0.99`) from `per_barcode.tsv`,
  and a one-line **sanity check** that aggregate `reads_num` = Σ per-barcode
  `reads_num` per model.
- End with a short read of the result: judge the A/B on `identity`, `mmRate`,
  and the `identity≥0.x` fractions (esp. `≥0.99`); state that `alignedRatio` is
  inflated (see gotcha) and is not the deciding metric; and list the barcodes that
  were **skipped (no first-gen reference)** on both sides, since they are excluded
  from the pooled aggregates.

## Notes / gotchas

- **`alignedRatio` is inflated on full-plasmid CCS data.** Each well mixes the
  target amplicon with ~2.7 kb circular plasmid background that matches only a
  ~48 bp shared vector at the reference's 3' end; gsmm2 counts that as aligned.
  `alignedRatio` lands ~0.94–0.99 almost regardless of quality, so judge the A/B
  on `identity` and the `identity≥0.9x` fractions, not `alignedRatio`.
- **The two models must share the same run.** Only rows of `--mapping` whose
  `RUN号` equals `--run` are used. A demux dir is expected to hold `Barcode01..N`
  for that run; a `BarcodeNN.fastq` present but not covered by the mapping (e.g.
  an unmapped extra barcode) is skipped with a log line.
- **Barcode naming must already be `BarcodeNN.fastq`.** This skill does not do
  the `24标签-N ↔ Adaptor-* ↔ barcode-NN` remapping. If the demux dirs are still
  in `Adaptor-barcodeNNN-M.fastq` form, build the `BarcodeNN.fastq` view first
  (see the `barcode24renaming` skill / `prepare_second_batch_runs.py`).
- **Short/complete references.** If the 一代 reference is shorter than the target
  (e.g. STR33 series 27-1/27-9, refs ~377 bp vs ~750 bp target), low
  `queryCoverage`/`identity` is an **incomplete reference**, not a sequencing or
  calling-model problem — treat those barcodes separately when judging the A/B.
- **Idempotency.** Per-barcode `*-metric` dirs from a prior run are reused by
  `sequencing_report_v2` unless `--rq-range` or `--np-thr` changed (the mm2 metric
  file is keyed by path, not by filter). To force a full recompute,
  delete `--outdir`. The comparison tables are regenerated from the aggr CSVs
  each run.
- **`--short-aln` is intentionally OFF.** Reads are full plasmids and refs are
  ~300–900 bp, so the short-alignment mode (query or target in [30,200]) does not
  apply. `fastq2bam` writes the `rq` tag; `--rq-range` (if given) filters on it.
