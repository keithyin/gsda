---
name: smc-ab-test-ref-locus-single-poly-first-base-analysis
description: Post-process a single-reference locus_error_rate joined table (locus_accuracy_joined.csv) to A/B compare the two calling models at homopolymer (polyC/G, run length 4 or 5) boundaries — the first base of each run — reporting eq/depth min / max and at the 5/25/50/75 percentiles per (base, repeatCount) plus a per-locus table, with built-in input verification.
---

# SMC A/B Test — single-ref poly-homopolymer first-base eq/depth

## Purpose

This skill is the **homopolymer-focused post-processing** step on top of the
single-reference per-locus A/B that the sibling skill
`smc_ab_test_ref_locus_single` produces. That skill gives you one
`locus_accuracy_joined.csv` — every reference base position, both sides,
ctrl/exp counts — for a run that maps to **one** reference (no barcodes).
This one answers the narrower question:

> **At the homopolymer boundaries (default: the first base of every polyC/G run
> of length 4 or 5), how do the two calling models differ, at min / the
> 5/25/50/75 percentiles / max of eq/depth?**

Homopolymer *run interiors* are near-perfect on both models and carry no
discrimination. The **first base of each run** (the run's 5′/start edge) is the
indel-sensitive position: the model must decide whether the consensus has the
right *number* of repeats, and that decision is scored at the first base. This
skill isolates exactly those positions and compares `control` vs `experiment`.

It:
1. picks the one reference sequence the joined table's `pos` column indexes
   (the single contig, or the contig named by the table's `refname`) and finds
   every C/G (default) homopolymer run whose length is in `--lens` (default
   `{4,5}`) — run detection on the raw reference sequence, `pos` 0-based,
2. for each run keeps the requested base — `--which first` (default), `last`,
   or `all`,
3. pulls the matching `(refname, pos)` rows from the joined table,
4. writes a per-locus table, a per-`(base, repeatCount)` comparison of **eq/depth**
   (min, the `--pctl` percentiles — default 5/25/50/75 — and max, plus Δ), and a
   markdown report.

## When to use

- the single-ref per-locus A/B (`smc_ab_test_ref_locus_single`) is done, and you
  want to know **where the difference is** at **repeat boundaries** —
  specifically polyC/G runs of a given length.
- you want a **repeat-length-stratified** readout: first-base eq/depth (min,
  5/25/50/75 percentiles, max), grouped by `C`/`G` and by run length.
- whole-genome (e.g. an E. coli/MG1655) single-reference comparison — not a
  barcode-per-plasmid set. For the barcode case use the sibling
  `smc_ab_homopolymer_firstbase` instead.
- a joined table that carries **exactly one `refname`**. A multi-`refname` table
  is refused (`smc_ab_test_ref_locus_single` always produces a single one) —
  that configuration belongs to the barcode sibling.

> **Pick the sibling skill by the reference layout:**
> - **this skill** — one run, **one** reference, poly-homopolymer first-base
>   eq/depth (min / 5/25/50/75 pctl / max) (single-ref `locus_accuracy_joined.csv`).
> - `smc_ab_homopolymer_firstbase` — one run, **many barcodes each with its own
>   plasmid reference**; first/last/all-base A/B with pooled `locus_accuracy`
>   **and** `locus_accuracy_by_depth` (barcode `joined_all.csv` + a
>   plasmid↔barcode↔RUN mapping).
> - `smc_ab_test_ref_locus_single` — the upstream step that builds the single-ref
>   `locus_accuracy_joined.csv` this skill consumes.

> This skill **does not** re-run gsmm2/gsetl. It reads the
> `locus_accuracy_joined.csv` that `smc_ab_test_ref_locus_single` already wrote.
> Run that skill first if you don't have it yet.

## Required inputs

- `--joined` — the tab-separated `locus_accuracy_joined.csv` from
  `smc_ab_test_ref_locus_single` (columns `refname pos eq_ctrl diff_ctrl ins_ctrl
  del_ctrl depth_ctrl aroundBases_ctrl locus_accuracy_ctrl
  locus_accuracy_by_depth_ctrl …_exp`). It must carry a **single `refname`**.
- `--ref` — the reference FASTA the run was aligned to (the same `--ref` used by
  `smc_ab_test_ref_locus_single`). A **single-contig** FASTA is used as-is (its
  header name need not match the table's `refname`). For a **multi-contig**
  FASTA, run detection uses the contig named by the table's `refname`; if that
  name is not among the contigs the run aborts rather than guessing.
- `--outdir` — where the three output files are written.

## Optional inputs

- `--char` — the characters that count as the homopolymer. **Default `CG`**
  (polyC/G). Pass e.g. `ACGT` for all homopolymers, or `A` for polyA only.
- `--lens` — comma-separated run **lengths** to consider. **Default `4,5`**
  (exactly runs of length 4 or 5, not "≥").
- `--which` — `first` (default) | `last` | `all`: which base(s) of each run to
  compare.
- `--pctl` — comma-separated percentile ranks, each an **integer in [0, 100]**.
  **Default `5,25,50,75`.** Anything else (a non-integer, or an out-of-range rank)
  is rejected rather than silently extrapolated. **min and max are always
  reported** as their own columns, in addition to the `--pctl` ranks — the
  `--pctl` list only controls which interior percentiles are added (they are
  sorted ascending; duplicates are ignored).

## Execution

```bash
/root/miniconda3/envs/gseda/bin/python \
    /root/projects/gsda/.claude/skills/smc_ab_test_ref_locus_single_poly_first_base_analysis/scripts/poly_first_base_analysis.py \
    --joined .../Run0002_adapter-v4-ecoliMG1655_locus_ab/locus_accuracy_joined.csv \
    --ref    /data1/REF_GENOMES/MG1655.fa \
    --outdir .../Run0002_adapter-v4-ecoliMG1655_locus_ab/poly_first_base \
    --char CG --lens 4,5 --which first
```

> **Use the `gseda` env python** (`/root/miniconda3/envs/gseda/bin/python`) —
> `server/bin` shadows `$PATH` on this machine so `conda activate`/`conda run`
> silently keep the py3.8 `server` python. This script needs **polars** (which
> the `gseda` env has); use the gseda interpreter.

This is a table-only step — it reads the joined table once and does the run
detection in memory. On a whole-genome 4.6 Mb reference (4.6 M loci) it runs in
a few seconds. Run it in the foreground; no backgrounding needed.

## Output (under `--outdir`)

`<stem>` = `poly<CHAR>_lens<len-lens>_<which>base`, e.g. `polyCG_lens4-5_firstbase`
for the default `--which first`, `polyCG_lens4-5_lastbase` for `--which last`,
`polyCG_lens4-5_allbase` for `--which all`. `<CHAR>` is the `--char` set sorted
(`--char GC` and `--char CG` give the same stem), and the run lengths are sorted
and `-`-joined (`--lens 5,4` → `lens4-5`). Because `--which` is part of the
stem, running two `--which` modes into the same `--outdir` does **not** overwrite
one with the other.

- `<stem>_locus.tsv` — the selected positions, columns `refname locus base
  repeatCount aroundBases eq/depth/eq÷depth/diff/ins/del/locus_accuracy` for
  `_ctrl` then `_exp`. `locus` is the 0-based reference `pos`; `aroundBases`
  brackets the target base in `[ ]`.
- `<stem>_percentiles.tsv` — one row per `(base, repeatCount)` group:
  `eq_over_depth_{ctrl,exp}_{min,p5,p25,p50,p75,max}` (default `--pctl`; the
  interior percentiles follow `--pctl`, min and max are always present),
  `delta_{min,p5,p25,p50,p75,max}` (exp − ctrl), and `n` (the number of loci
  **finite on both sides** in that group).
- `<stem>_report.md` — human-readable report: the count/coverage header, the
  percentile comparison table, then the full per-locus table.

## The metric

- **`locus_accuracy_by_depth = eq/depth`** per reference base position. It is
  **undefined (NaN) where `depth == 0`** (no coverage at that base). The
  percentile comparison uses **loci that are finite on BOTH sides** — depth>0 on
  ctrl *and* exp. control and experiment are two separate consensus assemblies of
  the same SMC run, so an edge may be covered on one side and depth-0 on the
  other (how often that actually happens is data-dependent — the printed
  coverage breakdown tells you); restricting to the finite-both set means ctrl
  and exp are distributed over the *same* positions and the per-percentile Δ is
  a like-for-like comparison. `n` is that finite-both count.
  (This is the same finite-both rule the per-locus section uses.)
- The per-locus table also carries `locus_accuracy = eq/(eq+diff+ins+del)` for
  reference, but the **percentile comparison is on `eq/depth`** (the lenient
  form, which is `≥ locus_accuracy` at every position).
- **`eq/depth` is blind to insertions.** gsetl's `depth` counts the reads
  spanning a position; an insertion is recorded in a separate `ins` column (at
  the anchor position after the inserted bases) and never lowers `depth`. So a
  consensus that resolves the repeat length with an **extra** repeat unit at the
  edge shows `eq/depth == 1.0` there, while `locus_accuracy` (and the raw `ins`
  count) drops. This matters because "did the model get the number of repeats
  right" is exactly an insertion/deletion question: `del` and mismatches *do*
  move `eq/depth`, insertions do not. The script therefore prints an
  **insertion tally** alongside every `eq/depth` figure (pooled
  `locus_accuracy` for both sides, plus how many finite-both loci carry `ins>0`
  and how many of those are `eq/depth == 1.0` on both sides — i.e. invisible to
  the percentile/tally numbers). When the insertion tally is non-trivial, say so
  and read the per-locus `ins_*` columns rather than presenting `eq/depth` alone.
- **Δ = exp − ctrl** per percentile; positive means the experiment (new) model
  has higher first-base eq/depth at that percentile. Each percentile is rounded
  to 6 dp (the joined table's own precision) before the Δ is taken, so a printed
  Δ is exactly the difference of the two printed percentiles.

## Built-in verification (run before reporting — do not skip)

Before it prints anything, the script **self-checks the joined table** and exits
nonzero if either invariant fails (same invariants as the
`smc_ab_homopolymer_firstbase` sibling, adapted to one reference):

- **V1 — indexing.** The base bracketed `[..]` in `aroundBases_exp` must equal
  the reference base at the **0-based** `pos` column, for every row. gsetl's
  `pos` is **0-based** (pos 0 = first base). This guards against the classic
  off-by-one that would silently show the *second* base of each run.
- **V2 — metric integrity.** `locus_accuracy_by_depth` must equal `eq/depth`
  (to 1e-6, where depth > 0) and `eq ≤ depth` must hold, for both sides, every
  row.

The script prints
`[verify] OK: V1 checked <n>/<N> rows …; V2 checked ctrl <n>/<N> + exp <n>/<N> rows …`
on success. If it prints `INPUT VERIFICATION FAILED …` and exits nonzero, **stop
and fix the input**.

Both checks are **total — a row that cannot be checked is a failure, not a
skip**, so the `OK` line can never claim coverage it did not have. A missing
required column (`aroundBases_exp`, `eq_*`, `depth_*`,
`locus_accuracy_by_depth_*`), a non-integer `pos`, or a `pos` outside the
reference all abort with `INPUT VERIFICATION FAILED`. A `pos` outside the
reference is the signature of a `--ref` that is not the reference the table was
built on.

## Data accuracy is mandatory — reporting wrong numbers is unacceptable

The results in this skill are scientific conclusions about which calling model
performs better at homopolymer boundaries. **A wrong number in the report is
not a minor slip — it is an unacceptable failure.** Before you present any
result, treat every figure as an untrusted claim you must confirm against its
source:

- **Only report numbers the script actually emitted.** Pull every count, total,
  eq/depth value, proportion, and Δ straight from the script's stdout and the
  three output files. Do not paraphrase, round-then-invent, "remember", or
  estimate a value. If a number is not in the output, it does not exist.
- **Cross-check before you report.** Reconcile the per-locus and aggregate
  figures against one another:
  - `exp better + exp worse + equal` **must equal** the `finite-both` locus
    count (and their proportions must sum to 100%).
  - `disagreeing` **must equal** `exp better + exp worse`.
  - The sum of the per-group `n` values **must equal** the `finite-both` count,
    and **must equal** the `finite-both used for the percentiles` figure the
    script prints in the coverage breakdown (every finite-both locus falls in
    exactly one `(base, repeatCount)` group).
  - `selected` **must equal** `target loci − missing`, and `finite-both` **must
    equal** `selected − depth-0 (either side)`.
  - The pooled `eq/depth` **must equal** that side's `Σeq / Σdepth` from the
    totals line, and the pooled `locus_accuracy` **must equal** `Σeq /
    Σ(eq+diff+ins+del)`; each printed per-locus `eq/depth` must equal that
    locus's own `eq / depth`.
  - A percentile must sit between the group's min and max eq/depth; a p50 that
    is not between p25 and p75 is a bug, not a finding.
  - Every top-N divergent locus must actually be a run boundary base in the
    target set (a locus you can't place is a mis-mapping).
- **Never present a failed or empty result as success.** If `[verify]` fails,
  if `sub.height` is 0, if `missing` is nonzero, or if any cross-check above
  does not balance, **say so explicitly and stop** — do not fill the gap with a
  plausible-looking number.
- **Flag the caveats the data carries.**
  - **Depth-0 loci** (NaN, excluded from the percentiles) — the script prints
    the breakdown `ctrl=N, exp=N, either side=N -> finite-both=N`; quote it.
  - **`missing` is not "no coverage".** Depth-0 loci *are* rows in the joined
    table (gsetl emits every reference position on both sides), so a nonzero
    `missing` is a **key mismatch** — a different reference, or a truncated
    per-side table — and must be investigated, never reported as depth-0. The
    healthy value is `missing=0`.
  - The **insertion tally** — how many finite-both loci carry `ins>0`, and how
    many of those are invisible to `eq/depth`. Surface it: it bounds how much of
    the model difference the headline numbers cannot see.
  - For a multi-contig `--ref`, state which contig the run was detected on.
- **When unsure, re-run.** Re-running is cheap and idempotent; fabricating a
  consistent number is not. If you cannot trace a reported value to its source
  line, regenerate it rather than guess.

If any of the above cannot be satisfied, do **not** present the numbers —
report the failure and its cause instead.

## Presenting the final result

When reporting to the user:

- State **which base** (`--which`) and the run filter (`--char`, `--lens`) —
  e.g. "first base of each polyC/G run of length 4 or 5".
- Lead with the **per-`(base, repeatCount)` percentile table** (ctrl and exp
  min / p5 / p25 / p50 / p75 / max of eq/depth, plus Δ). Note which groups
  differ and in which direction — typically p50/p75/max are 1.0 on both sides
  and only the **low end (min and p5, the worst tail)** carries signal; report
  the full min→max spread, not just the quartiles.
- **Also present the per-locus (位点粒度) comparison — not just the summary.**
  The percentile table is the aggregate view; the user also wants to see which
  individual loci the two models treat differently. From the script's stdout
  "per-locus comparison" section (or `<stem>_locus.tsv`), show:
  - the **count of loci where the two models disagree** on eq/depth
    (finite-both, |Δ| > 0), and
  - the **group totals over finite-both loci** for each side: total **eq**, total
    **depth**, their ratio (eq/depth), and — from the same stdout line — the
    **diff / ins / del** totals and the pooled **locus_accuracy =
    eq/(eq+diff+ins+del)**, i.e. ctrl vs exp absolute counts, not just
    per-locus, and
  - the **overall win/lose/tie** over finite-both loci (Δ eq/depth = exp −
    ctrl): the count **and proportion** of loci where **exp is better**, where
    **exp is worse**, and where the two are **equal** (the three sum to the
    finite-both count). Lead the per-locus view with this one-line tally, and
  - the **most-divergent loci** (top-N by |Δ|): locus, base, repeatCount,
    aroundBases (target base in `[ ]`), and — **separately for ctrl and exp** —
    the raw **eq** count and **depth** count plus the derived eq/depth, and Δ.
    The script's stdout prints these as
    `c: <eq> / <depth> (<eq/depth>)` and `e: <eq> / <depth> (<eq/depth>)`.
  For a whole-genome run there can be thousands of loci — present the top-N
  divergent subset inline and point to `<stem>_locus.tsv` / `report.md` for the
  full per-locus table rather than dumping every row.
- Give the **selected-locus count** (e.g. 14,236 of 14,236), the printed
  **coverage breakdown** (`depth-0: ctrl=…, exp=…, either side=… →
  finite-both=…`), and, if useful, the further split of eq/depth == 1.0 vs < 1.0
  per side, plus how many loci the two models **disagree** on.
- State the **insertion tally** in one line: how many finite-both loci carry
  `ins>0` on either side, and how many of those the eq/depth numbers cannot see
  (eq/depth == 1.0 on both sides). If that masked count is a sizeable share of
  the disagreeing loci, say plainly that the eq/depth tally understates the
  difference and that the `Δ locus_accuracy` (and the `ins_*` columns) are the
  other lens on the same loci.
- Say explicitly that this is **alignment identity vs the Sanger reference** at
  repeat boundaries (FASTQ/unmapped input → no rq/np filtering), measuring how
  well each model resolves the *repeat length* at the run's edge, not raw
  per-base quality.
- The **absolute paths** of `<stem>_report.md` and `<stem>_locus.tsv` **must**
  appear in the final message.

### Final message requirement

Include the **absolute paths** of `<stem>_report.md`, `<stem>_percentiles.tsv`,
and `<stem>_locus.tsv`; the selected-locus count and the **coverage breakdown**
(depth-0 per side and finite-both); the **insertion tally** (loci with `ins>0`,
and how many of those are invisible to eq/depth); the per-group percentile
comparison (both sides, Δ); **and** the per-locus comparison — the overall
win/lose/tie tally (count **and** proportion of loci where exp is better, worse,
or equal), the count of loci where the two models disagree, plus the top-N
most-divergent loci. For each
top-N locus, report the **absolute eq and depth counts for BOTH the control and
experiment groups** (not just the eq/depth ratio): locus, base, repeatCount,
aroundBases, `eq_ctrl / depth_ctrl (eq/depth)`, `eq_exp / depth_exp
(eq/depth)`, and Δ.

## Notes / gotchas

- **`pos` is 0-based** in the gsetl `fact_aligned_bam_ref_locus_info` table
  (pos 0 = first reference base). Run boundaries and any position math use
  0-based indices. V1 enforces this — a 1-based table is refused.
- **eq/depth is NaN at depth 0 — `drop_nulls()` does NOT drop NaN.** Polars
  `Float64` nulls and NaN are different: a 0-depth position is written as NaN,
  and `Series.drop_nulls()` keeps it. If you percentile over NaN values you get
  **non-monotonic, wrong percentiles** (NaN sorts first). Always filter to
  finite values (`v is not None and v == v`) before taking percentiles. This
  script does that.
- **`--lens` is exact, not "≥".** `--lens 4,5` selects runs of length exactly 4
  *or* 5, not runs of length ≥ 4. Pass `4,5,6,7` to widen.
- **First base ≈ the indel-sensitive position.** Interior/trailing bases of a
  run are typically eq/depth ≈ 1.000 on both sides and never discriminate. If
  your percentiles look flat at 1.0 everywhere, confirm you selected `first`,
  not an interior base. Note the "first" base is the leftmost base on the
  **reference forward strand**: for a read aligned to the reverse strand it is
  that read's 3′ end, so both orientations are pooled. (On the E. coli/MG1655
  Run0002 data, `--which first` discriminates ≈19× more loci than `--which
  last`: 719 vs 38.)
- **Depth-0 loci are common at run edges** (the read starts/ends inside the run).
  They show as NaN in the per-locus table and are excluded from percentiles. They
  are **rows in the joined table**, so they never show up as `missing` — the
  coverage breakdown the script prints is what tells you how much was absent.
- **Multi-contig reference.** Run detection is done on the single contig the
  joined table's `refname` names (a single-contig FASTA is used as-is). If the
  `refname` is not among the FASTA's contigs, the run aborts rather than
  concatenating — a concatenation would invent targets on contigs the table does
  not cover. Prefer passing just that contig.
- **Idempotent within one `--which`.** Re-running with the same `--outdir` and
  the same `--which` overwrites the three output files; a different `--which`
  writes a different stem (`…_lastbase`, `…_allbase`) instead of clobbering it.
