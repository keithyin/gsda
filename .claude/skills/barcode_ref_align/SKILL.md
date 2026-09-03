---
name: gsda-barcode-ref-align
description: Align barcode-split FASTQ reads against per-plasmid references and generate a sequencing report (identity/coverage/mismatch metrics) via sequencing_report_v2.
---

# Barcode FASTQ vs Reference Alignment Report

## Purpose

Compare barcode-split FASTQ reads (PacBio full-plasmid CCS reads, one FASTQ
per barcode) against a set of reference sequences, and produce alignment
metrics (identity, coverage, mismatch rate, indel rate) per sample plus a
cross-sample summary.

The relationship between the reads and the references is given by a
`plasmid -> barcode` mapping file.

Pipeline per sample:

1. Extract the per-plasmid reference from a multi-record reference FASTA.
2. Convert the barcode FASTQ to an unaligned BAM (with an `rq` tag) via
   `fastq2bam.py`.
3. Run `gseda.ppl.sequencing_report_v2` to align query vs reference and write
   aggregated metrics.
4. Collect all per-sample aggr metrics into one `summary.tsv`.

## When to use

Use this skill when the user wants to:

- align / compare barcode-split FASTQ data against reference sequences
- generate a sequencing report (`sequencing_report_v2`) for multiple barcodes
- produce identity / coverage / mismatch metrics per plasmid or per barcode
- build a summary table across many barcodes/samples

## Required inputs

- `--barcode-dir`: directory containing `BarcodeNN.fastq` files
- `--ref-fasta`: multi-record reference FASTA (e.g. `intersect.fa`). The sample
  to compare against is the one whose header matches the plasmid name.
- `--mapping`: `plasmid<TAB>barcode` TSV (with a header line) linking each
  barcode's reads to a plasmid, and (via plasmid name) to a reference record.
- `--outdir`: where BAMs, per-sample metric dirs, and `summary.tsv` go.

## Optional inputs

- `--rq-range`: rq (per-read accuracy) filter for the report. Default `0.99:1.1`.
  Means "keep reads whose rq is in [low, high]".
- `--conda-env`: conda env used to run the tools. Default `py38` (gseda is
  installed there; `gsmm2-aligned-metric` is on PATH).
- `--ref-prefix`: prefix stripped from reference ids to recover the plasmid
  name. Default `260827STR-` (so `260827STR-93-1` -> `93-1`). Set to empty
  string if reference ids already equal plasmid names.
- `--short-aln`: 1 to enable short-alignment mode (use when query or target
  length is in [30, 200]). Default 1.
- `--fastq2bam`, `--report`: override the tool script paths.

## Execution

Run the bundled script (it handles mapping, reference extraction, the
per-sample loop, and the summary in one go, inside the conda env):

```bash
python3 /root/projects/gsda/.claude/skills/barcode_ref_align/scripts/run_report.py \
    --barcode-dir <dir with BarcodeNN.fastq> \
    --ref-fasta   <reference intersect.fa> \
    --mapping     <plasmid_name_2_barcode.tsv> \
    --outdir      <output dir> \
    --rq-range    0.99:1.1
```

The job is CPU-bound but fast per sample (a few seconds to ~1 min depending on
read count). For many samples it is fine to run it as a background task and
poll the output.

## Output

Under `--outdir`:

- `refs/<plasmid>.fa` — the extracted per-plasmid reference
- `BarcodeNN_<plasmid>.bam` — unaligned BAM (with `rq` tag)
- `BarcodeNN_<plasmid>-metric/`
  - `...gsmm2_aligned_metric_aggr.csv` — aggregated metrics (`name\tvalue`)
  - `...gsmm2_aligned_metric_fact.csv` — per-read fact metrics
  - `...basic.csv` — per-read basic info
- `summary.tsv` — one row per sample with the key metrics

## Key metrics (in summary.tsv / aggr csv)

- `reads_num`, `tot_bases`, `n50`, `read_len_p50` — read stats
- `alignedRatio` / `notAlignedRatio` — fraction of reads that aligned
- `queryCoverage`, `queryCoverage3` — how much of each read was covered
- `identity`, `identity-p50` — match rate over aligned span
- `mmRate` — mismatch rate
- `longIndelRatio`, `GlobalQueryCoverage`
- `identity≥0.83` / `≥0.90` / `≥0.99` — fraction of reads meeting the threshold

## Notes / gotchas

- **Reference coverage is often partial.** The reference is usually a short
  "intersect" region (a few hundred bp) while the reads are full plasmids
  (~1.5 kb). Low `alignedRatio` / `queryCoverage3` is therefore expected and
  does not indicate a problem.
- **Not every plasmid has a reference.** Samples whose plasmid has no matching
  record in `--ref-fasta` are skipped (the script logs each SKIP). Report to
  the user which plasmids/barcodes were skipped.
- **Plasmid names may not match reference ids exactly.** Adjust `--ref-prefix`
  (or pre-map) if the reference header convention differs from the plasmid
  naming in the mapping file. Verify the mapping before trusting results.
- Runs are per-sample independent and idempotent per output path, but existing
  `*-metric` dirs from a previous run with a different rq-range will be
  reused. Delete `--outdir` to force a full recompute.
