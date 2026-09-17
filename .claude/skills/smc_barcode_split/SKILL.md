---
name: smc_barcode_split
description: Split an SMC (consensus) all-reads BAM into per-barcode BAMs, grouping by subread id taken from the barcode-split FASTQ headers, and optionally split the source adapter (subread) BAM the same way by joining on the channel (`ch`) field.
---

# SMC BAM → per-barcode BAM split (+ matching adapter BAM split)

## Purpose

Given the SMC (consensus) all-reads BAM for a run and the barcode-split
FASTQ directory (one `BarcodeNN.fastq` per barcode), dump the SMC BAM into
one BAM per barcode. Each SMC read is routed to the barcode whose FASTQ
contains that subread, matched on the read's queryname.

With `--adapter-bam`, the same partitioning is then applied to
`<Run>_adapter.bam` — the subread BAM that smicing consumed to build the SMC
BAM — so each barcode gets its raw subreads alongside its consensus reads.
The two files are joined on the **channel** (`ch`): the SMC read carries the
channel it was consensused from, and every adapter record carries the channel
it came out of. The adapter split is derived from the SMC result, never from
the FASTQs, so the two output sets agree by construction on which channel
belongs to which barcode.

The source BAMs and the FASTQs are never modified.

## When to use

Use this skill when the user wants to:

- split an SMC all-reads BAM into per-barcode BAMs
- associate SMC consensus reads with their demux barcode
- produce a `barcoded_smc/` set of BAMs mirroring a `barcode_assign/` set of FASTQs
- split the run's `*_adapter.bam` (subreads) by barcode, keeping raw
  subread passes for the same molecules the consensus set covers — e.g. to
  re-polish, re-basecall, look at per-pass error, or feed a subread-level
  analysis that needs adapter/barcode sequence still in the read

## Inputs

- `--barcode-dir`: directory with the per-barcode FASTQs (e.g. `barcode_assign/`)
- `--smc-bam`: the SMC all-reads BAM to split
- `--outdir`: output directory for the per-barcode BAMs (created if missing)
- `--adapter-bam` (optional): the run's `_adapter.bam` to split the same way
- `--adapter-outdir` (optional): where its outputs go (default: `--outdir`; use
  a separate dir such as `barcoded_adapter/` so the two sets stay apart)

The SMC split is usually already done when the adapter need comes up, so the
channel map can be taken from its output instead of re-reading the big SMC
BAM: pass `--smc-dir <barcoded_smc>` (per-barcode SMC BAMs) in place of
`--smc-bam`/`--barcode-dir`, together with `--adapter-bam`. That is the
adapter-only mode.

## Execution

Both passes in one command (SMC split, then adapter split using the channels
it just saw):

```bash
python3 /root/projects/gsda/.claude/skills/smc_barcode_split/scripts/split_smc_by_barcode.py \
    --barcode-dir     <dir with BarcodeNN.fastq> \
    --smc-bam         <Run>.smc_all_reads.bam \
    --outdir          <dir>/barcoded_smc \
    --adapter-bam     <Run>_adapter.bam \
    --adapter-outdir  <dir>/barcoded_adapter
```

Adapter only, off an existing SMC split:

```bash
python3 /root/projects/gsda/.claude/skills/smc_barcode_split/scripts/split_smc_by_barcode.py \
    --smc-dir        <dir>/barcoded_smc \
    --adapter-bam    <Run>_adapter.bam \
    --outdir         <dir>/barcoded_adapter
```

Add `--index` to run `samtools index` on each output BAM (both sets). Add
`--unassigned` to also dump SMC reads that match no barcode into
`<outdir>/<smc-stem>.unassigned.bam`, and `--unassigned-adapters` for the
adapter equivalent (that one can be tens of GB — see gotchas).

Cost: the SMC pass is one streaming pass over the SMC BAM (plus reading the
FASTQ headers once up front), roughly 1–2 min for a ~360k-read BAM. The
adapter pass is one full pass over a 15–50 GB file and writes roughly its own
copy of everything it matched, so budget disk before starting. Measured on
a 18.9 GB / 4.39M-record `adapter.bam` with 16 barcodes: **385 s** end to end
(~11k records/s) writing 17 GB, `--index` included, at the default
`--threads 4`. Writing is the bottleneck — `--threads 1` runs ~10x slower —
and `--threads` is per output BAM, so on a small box lower it (peak threads
≈ `--threads` × barcodes touched). Check free disk on the target volume first.

On a many-core box the scaling is near-linear, so those numbers are pessimistic:
`--threads 8` on a 256-core host sustains **85–100k records/s**, i.e. the whole
two-pass split finishes in **~1.5–4.5 min** rather than needing backgrounding.
Two measured runs at `--threads 8`: a 31.4 GB / 7.88M-record adapter BAM with 21
barcodes → **89 s** writing 29 GB (5.9 TB free on the target volume); a 22.4 GB /
5.83M-record one with 15 barcodes → **266 s** writing 21 GB (`--index` included).
A third on the same 18.9 GB / 4.39M-record file as the threads=4 measurement above,
also `--threads 8` and no `--index`, took **42 s** at 117k records/s — ~9x the
threads=4 wall time, which cannot be the thread count alone; a warm page cache on
the input (360 GB of buff/cache on that host) is the likely difference. Treat
threads=8 as "a few minutes", not as a 2x-on-385s improvement.
At that speed run it as a plain foreground command with a generous timeout and
read the log afterwards; backgrounding only pays off above ~50 GB. The `.bai`
files `--index` leaves on these all-unmapped BAMs are 32-byte no-coordinate
artifacts that get deleted downstream anyway — skip the flag unless asked.

## qname mapping (SMC pass, important)

The FASTQ header (`@` stripped) and the BAM queryname are matched to decide
which barcode a read belongs to. These two strings are **not guaranteed to be
byte-identical**, so the script normalizes **both** sides with the same rule
before matching — this is done unconditionally, there is no flag to turn it
off. The rule keeps only the first two `/`-separated fields:
`"/".join(qname.split("/")[:2])`. It is applied equally to the FASTQ headers
and the BAM querynames, so a header and a queryname that differ only after the
second `/` still join correctly (and if they are already ≤2 fields the mapping
is a no-op).

## channel mapping (adapter pass)

`ch` is the join key, and it is a per-run integer well/channel id:

| file | queryname | `ch` tag | records per channel |
|---|---|---|---|
| `<Run>.smc_all_reads.bam` | `<Run>/1076757` | `1076757` | exactly 1 (1:1 on the runs checked) |
| `<Run>_adapter.bam` | `read_1076757/1076757/subread/3` | `1076757` | median ~17, 1 to ~140 observed |
| `BarcodeNN.fastq` header | `@<Run>/1076757[/32-926]` | — | — |

So `ch` is the second `/`-field on both sides, and the tag agrees with it
whenever it is present. `--channel-tag` changes the tag; `--channel-qname-field`
picks the fallback field index used when a record has no tag (an adapter BAM
stripped of `ch` still splits correctly via field 1 — same channels, same
outputs).

Channels missing from the adapter BAM are *not* an error: `--smc-dir`/`--smc-bam`
defines the barcode sets, and each barcode's `adapter_ch` column then just falls
short of its `smc_reads` column.

**Channel ids are only unique within a run** — two runs both have a channel
`612233`. Since adapter querynames carry no run name, the script compares the
read-group run (`RG` `rn` in both BAM headers) before streaming the adapter BAM
and aborts if they differ; `--allow-run-mismatch` overrides. In `--smc-dir` mode
it additionally refuses up front if the directory holds splits of two different
source BAMs that repeat barcode names — narrow `--smc-glob` to one stem
(e.g. `Run0001.smc_all_reads.*.bam`) if you meant to mix.

## Output

Under `--outdir` (and `--adapter-outdir`), one BAM per barcode, named after the
**source** BAM's stem with the barcode name inserted as an infix before the extension:

```
<smc-bam-stem>.<BarcodeNN>.bam
<adapter-bam-stem>.<BarcodeNN>.bam
Run0001.smc_all_reads.bam + Barcode07.fastq -> Run0001.smc_all_reads.Barcode07.bam
Run0001_adapter.bam       + Barcode07.fastq -> Run0001_adapter.Barcode07.bam
```

The barcode part is the FASTQ base name without `.fastq`. Carrying the source
stem in the name lets BAMs from several runs share one `--outdir` without
colliding. Each BAM preserves its own source BAM header (the adapter BAMs keep
`dw`/`ar`/`be`/`cx` and the adapter/barcode bases, since records are copied
verbatim); with `--index`, each also gets a `.bam.bai`. The script prints a
reconciliation summary:

- `total SMC records` — records read from the SMC BAM
- `unmatched` — records whose queryname appeared in no barcode FASTQ
- `written records` — `total - unmatched`
- `channels mapped` — the size of the `ch -> barcode` map handed to the adapter pass
- per-output-file counts, keyed by the full output file name
- `total / written / unmatched adapter records`, `channels in map`,
  `channels covered` (how many of the SMC channels were found in the adapter BAM)
- a `per-barcode reconciliation` table: `smc_reads`, `adapter_ch`, `adapter_recs`,
  `recs/ch` (mean subread depth) per barcode

Two known-good runs on second-batch `20260805_250804Y0004_Run0001` (16 barcodes
each) — **from two different generations of input file living in the same run dir**.
They differ only in the 4th digit of every total, so you cannot tell afterwards
which pair you consumed; the PG line is the only discriminator.

Aug-7 generation (`Run0001.smc_all_reads.bam` + `Run0001_adapter.bam`), SMC pass
driven by `barcode_assign-renamed/`:

```
total SMC records: 347970    written: 323081    unmatched: 24889    channels mapped: 323081
total adapter records: 4389147   written: 4353156   unmatched: 35991
channels covered: 323081 (100.0% of the SMC channels)
adapter_ch == smc_reads on every row; per-barcode recs/ch 6.5-17.4
```

Sep-4 generation (`Run0001_called-demuxed-v4-.smc_all_reads.bam` +
`Run0001_called-demuxed.bam`), SMC pass driven by `called-barcode-v4-/demuxed/`:

```
total SMC records: 348238    written: 325652    unmatched: 22586    channels mapped: 325652
total adapter records: 4387215   written: 4353521   unmatched: 33694
channels covered: 325652 (100.0%); adapter_ch == smc_reads on every row; recs/ch 11.6-17.3
```

Same generation, second-batch `..._Run0002`, 45 s at `--threads 8`, 16G out:
SMC 270595 → 255277 (unmatched 15318); adapter 4069102 → 4042808 (unmatched 26294);
100.0% covered, `adapter_ch == smc_reads` on all 16 rows, recs/ch 12.1-19.4.

`channels covered: 100.0%` plus `adapter_ch == smc_reads` for every barcode is
the signature of a correct pairing — it says every consensus read's channel was
found among the subreads, and nothing stray joined in. Anything meaningfully
below that is worth reporting to the user rather than shipping.

## Notes / gotchas

- **More SMC reads than FASTQ subreads is normal.** The SMC BAM often
  contains subreads that failed consensus / demux, so `unmatched` is usually
  > 0. Report this gap to the user; pass `--unassigned` if they want those
  reads captured.
- **More adapter channels than SMC reads is also normal** — a channel whose
  subreads never reached consensus has no SMC read, so its subreads are
  dropped by the adapter pass. That is the "SMC reads as the baseline"
  behaviour the user asked for; if they want every subread of every channel
  that *did* demux, that's a different key (the demux `barcode.tsv`), not this
  skill.
- **A run can hold several `BarcodeNN.fastq` dirs, and they are not
  interchangeable.** Demux output (`*_called-barcode-v4-/demuxed/`), a renamed
  copy, and the set a downstream report was built on commonly differ by ~1% of
  read membership (measured on second-batch `Run0001`: 325652 vs 323081 mapped
  ids, and 10287 vs 10143 on `Barcode01` alone). `--barcode-dir` silently picks
  one, so if these BAMs must agree read-for-read with an existing report, confirm
  which dir that report consumed before splitting. Likewise a barcode holding only
  a handful of reads (1–3) is the FASTQ's own content, not the split dropping
  reads — check `awk 'NR%4==1' BarcodeNN.fastq | wc -l` before chasing it. Second
  batch `Run0001` has `Barcode16`=2 and `Run0002` has `Barcode15`=3 / `Barcode16`=2;
  those runs' real samples are Barcode01–15 and Barcode01–14 respectively.
- **The adapter BAM must be the one smicing consumed.** `<Run>_adapter.bam`
  next to `<Run>.smc_all_reads.bam` is the right pairing; a `_called.bam` or a
  `_called-demuxed*.bam` from a different pipeline revision is not (the run-name
  guard catches the obvious cases, not a right-run-wrong-file mixup). Check the
  SMC BAM's PG line — `smicing consensus ... <Run>_adapter.bam` names its input.
  On the STR second batch this is not a hypothetical: `..._Run0001/` holds **both**
  generations side by side — Aug-7 `<Run>_adapter.bam` (18963035082 B) +
  `<Run>.smc_all_reads.bam`, and Sep-4 `<Run>_called-demuxed.bam` (18942260791 B) +
  `<Run>_called-demuxed-v4-.smc_all_reads.bam`. The two adapter BAMs are within 0.1%
  of each other's size, so nothing but the PG line distinguishes them; `-v4-` marks
  the *output prefix* of the newer smicing run, not a version of the input.
- **Output volume.** Per-barcode adapter BAMs together are close to the size of
  the input `--adapter-bam` minus the dropped channels, and each barcode is its
  own bgzf stream, so small barcodes compress worse than one big file. Don't
  use `--unassigned-adapters` unless asked; it re-adds everything.
- **Record order.** Outputs keep the source file's order: the adapter BAM is not
  sorted by channel, but a channel's subreads are contiguous, so each output BAM
  has each channel's passes together. `--index` on these all-unmapped BAMs yields
  a 32-byte no-coordinate `.bai` (samtools accepts it, region queries are still
  useless); for name-keyed access use `samtools sort -n -@ 8`.
- **qname mismatch yields near-zero output.** If the FASTQ header and BAM
  queryname don't line up on their first two `/`-fields, every read is
  `unmatched` and the script exits non-zero rather than leaving an empty output
  set. Since the two-field mapping is forced, that failure means the convention
  is *not* "same first two fields" — check one FASTQ header against one BAM
  queryname by eye (strip the `@`) and adjust the `map_qname` rule in the script
  if the run uses a different key layout.
- **Zero adapter matches = wrong pairing**, and the script says so instead of
  quietly writing nothing.
- **Idempotent per output path**, but re-running overwrites files in
  `--outdir` / `--adapter-outdir`. Delete them to force a clean recompute.
- Requires `pysam` (and `samtools` only for `--index`). `pysam` here is built
  without `compression_level`/`ThreadPool` support, so writer compression is
  htslib's default and only `threads=` is tunable.
