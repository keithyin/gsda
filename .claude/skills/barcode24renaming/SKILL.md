---
name: barcode24renaming
description: Rename barcode-split FASTQ files into a new directory using a name-mapping TSV (original filename -> barcodename), copying only the mapped files.
---

# Barcode FASTQ Renaming by Mapping

## Purpose

Rename a directory of barcode-split FASTQ files (`<original-name>.fastq`)
using a name-mapping table, producing a new output directory where each file
is renamed to its mapped `barcodename`. The source directory is never
modified — files are **copied**, not moved.

Typical use: a run's `barcode_assign/` directory holds files named by
instrument barcode index (e.g. `Adaptor-barcode251-4.fastq`). A mapping TSV
relates each such file to a friendly `barcodename` (e.g. `barcode-01`).
This skill rebuilds the file set under the friendly names.

## When to use

Use this skill when the user wants to:

- rename barcode FASTQ files to a set of barcodenames
- remap `barcode_assign/` output names to a `barcodename` mapping
- build a renamed/cleaned directory of barcode FASTQ files from a mapping table

## Required inputs

- `--src`: directory containing the original `<name>.fastq` files.
- `--outdir`: where the renamed copies go (created if missing).

## Optional inputs

- `--tsv`: a mapping TSV to use instead of the built-in one. Column layout:
  col1 = original file base name **without** `.fastq` (matches a file in
  `--src`), col2 = new barcodename, any further columns (e.g. a `seq` column)
  ignored. If omitted, the **built-in default mapping** (the
  `barcode-2-barcodename.tsv` table, embedded in the script as `DEFAULT_MAPPING`)
  is used.
- `--include-unmapped`: also copy files that have no mapping, keeping their
  original name. Without this flag, only mapped files are copied (default).

The built-in default mapping (the `Adaptor-barcode... -> barcode-NN` table of
24 entries) is embedded directly in the script, so the common case needs no
external file. The source TSV's header is mislabeled (`barcodename barcode
seq`); the data rows are authoritative — col1 holds the `Adaptor-...` name
(the key), col2 the `barcode-NN` target.

## Execution

```bash
python3 /root/projects/gsda/.claude/skills/barcode24renaming/scripts/rename_by_mapping.py \
    --src    <dir with Adaptor-barcodeNNN-M.fastq> \
    --outdir <dir>barcode_assign-renamed
```

Omitting `--tsv` uses the built-in `Adaptor-barcode... -> barcode-NN` mapping.
Pass `--tsv <file>` to use a different mapping instead. Add
`--include-unmapped` if the user wants unmapped files carried over too.

## Output

A new directory (`--outdir`) containing `<Barcodename>.fastq` for every source
file that had a mapping (plus unmapped originals if `--include-unmapped`).
The target name from the mapping is normalized: the dash is dropped and the
first letter capitalized, so a mapping value `barcode-01` becomes
`Barcode01.fastq` (see `rename_target`). The script prints a reconciliation summary:

- `source files` — total `.fastq` files found in `--src`
- `mapping entries` — rows in the TSV
- `renamed (copied)` — files actually copied and renamed
- `unmapped (kept as-is)` — files with no mapping (only nonzero with the flag)
- `mapped-but-missing` — TSV rows whose source file was not found (listed)

## Notes / gotchas

- **Mapping rarely covers every file.** In practice the TSV holds only a
  subset of the run's barcodes. The default behavior (copy mapped only) means
  the output count will be far smaller than the source. Report the
  reconciliation numbers to the user so the gap is explicit.
- **The script is a copy, not a move.** `--src` is left intact; re-running into
  the same `--outdir` overwrites previously written files with identical names.
- **Key vs. target column.** The mapping keys on col1 and renames to col2.
  If a given TSV stores these in a different order, confirm which column
  matches the source filenames before trusting the result.
- **Timestamps preserved** via `shutil.copy2`.
