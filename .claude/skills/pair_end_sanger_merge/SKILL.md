---
name: gsda-pair-end-sanger-merge
description: Merge Sanger (一代) double-ended .seq reads into a full-length / overlap consensus per sample using mappy, and emit a detailed per-sample merge stats table.
---

# Pair-End Sanger Result Merge

## Purpose

Merge double-ended Sanger (一代测序) reads into a single sequence per sample
and report per-sample alignment metrics.

Input is a directory of `.seq` files (plain sequence, no FASTA header). Each
sample is sequenced with an M13 forward primer and an M13 reverse primer; the
two reads are matched to each other with `mappy` and merged.

Two output modes:

- **intersect (default)** — emit only the region where the forward and reverse
  reads overlap (the high-confidence consensus region).
- **union (`--union`)** — emit the full-length merged sequence (5'/3' ends
  taken from whichever read extends farther).

This skill reuses the merge logic in
`gseda.ppl.pair_end_sequencing_result_merge` but adds structured stats output
(`merge_stats.tsv`) that the bare module does not write.

## When to use

Use this skill when the user wants to:

- merge / combine double-ended Sanger `.seq` reads into one sequence
- merge STR / plasmid first-generation sequencing results
- check forward/reverse read consistency (identity, coverage) per sample
- produce a merge stats table across many samples

## Required inputs

- `input_dir` (positional): directory containing the `.seq` files. Filenames
  must encode the strand in the second `_`-separated token, e.g.
  `260731STR35-2_M13F-47_...D05.seq` (forward) and
  `260731STR35-2_M13R-48_...D06.seq` (reverse). The sample id is the part
  before the first `_`.

## Optional inputs

- `--outdir`: output directory. Default `<input_dir>/merged_output`.
- `--union`: emit the full-length union sequence instead of the overlap
  intersect region.

## Execution

Run the bundled script (needs the `gseda` package + `mappy` importable, e.g.
the `py38` conda env where gseda is installed):

```bash
python3 /root/projects/gsda/.claude/skills/pair_end_sanger_merge/scripts/run_merge.py \
    <input_dir> \
    --outdir <output_dir>        # optional
    --union                      # optional, for full-length output
```

Job is CPU-bound and fast (a few seconds for tens of samples). No need to
background it unless there are thousands of samples.

## Output

Under `--outdir`:

- `<sample>.intersect.fa` (or `<sample>.union.fa`) — per-sample merged sequence
- `intersect.fa` / `union.fa` — aggregated multi-record FASTA
- `merge_stats.tsv` — one row per sample:

| column | meaning |
|---|---|
| `sample` | sample id (prefix before first `_`) |
| `f_len` / `r_len` | forward / reverse raw read length |
| `overlap_len` | F–R overlap (intersect) length |
| `output_len` | length of the emitted sequence |
| `f_r_identity` | identity of the F vs Rrc merge alignment |
| `f_realign_identity` / `f_query_cov` / `f_target_cov` | forward read realigned to output: identity / queryCoverage / targetCoverage |
| `r_realign_identity` / `r_query_cov` / `r_target_cov` | reverse read (revcomp'd) realigned to output |
| `status` | `OK` / `LOW_IDENTITY` / `NO_HIT` / `MISSING_PAIR` |

## Interpreting results

- **`NO_HIT`** — no primary mappy hit between forward and reverse reads; no
  sequence produced. Usually a poor/short/contaminated read. Inspect the raw
  `.seq` files for that sample.
- **`LOW_IDENTITY`** — merged, but `f_r_identity < 98%`. Look at which strand's
  realign identity is low: if `f_realign_identity` is far below
  `r_realign_identity` (or vice versa), that strand is the likely culprit
  (sequencing quality, or an STR repeat region producing stacked indels).
- **Short reads** — a read much shorter than its partner (e.g. 300 bp vs
  ~1000 bp) shrinks the overlap region; identity is still valid but the
  intersect output will be short. Use `--union` to recover full length.
- `targetCoverage` of the realigned read is ~1.0 when that read spans the whole
  emitted region; a low value means the emitted region extends beyond that read.

## Notes / gotchas

- Strand is inferred from the **second** filename token: `M13F*` → forward,
  `M13R*` → reverse. Any other primer prefix is skipped. If the data uses a
  different primer naming, the files will be dropped silently — check the
  `.seq` filenames first.
- `MISSING_PAIR` means a sample has only one strand present.
- Existing files in `--outdir` are overwritten for samples in the input; stale
  files from a prior run with more samples are not removed. Delete `--outdir`
  to force a clean recompute.
