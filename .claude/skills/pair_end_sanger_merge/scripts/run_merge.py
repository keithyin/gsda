#!/usr/bin/env python3
"""Merge Sanger (一代) double-ended .seq reads per sample and emit a detailed
stats table.

Thin driver over ``gseda.ppl.pair_end_sequencing_result_merge`` — reuses its
``read_sanger_seqs`` / ``merge_pair`` / ``align_to_reference`` so the merge
logic stays in one place, but runs its own per-sample loop so it can collect
structured metrics and write:

  <outdir>/
    <sample>.intersect.fa   (or .union.fa)   per-sample merged sequence
    intersect.fa / union.fa aggregated multi-record FASTA
    merge_stats.tsv         one row per sample with all metrics
"""
import argparse
import os

from gseda.ppl.pair_end_sequencing_result_merge import (
    MIN_IDENTITY,
    align_to_reference,
    merge_pair,
    read_sanger_seqs,
    revcomp,
)

COLS = [
    "sample", "f_len", "r_len", "overlap_len", "output_len",
    "f_r_identity", "f_realign_identity", "f_query_cov", "f_target_cov",
    "r_realign_identity", "r_query_cov", "r_target_cov", "status",
]


def run(input_dir, outdir, union):
    mode = "union" if union else "intersect"
    os.makedirs(outdir, exist_ok=True)
    suffix = ".union.fa" if union else ".intersect.fa"
    sum_fa = os.path.join(outdir, f"{mode}.fa")

    by_sample = read_sanger_seqs(input_dir)
    rows = []
    merged_lines = []

    for sample in sorted(by_sample):
        seqs = by_sample[sample]
        if "F" not in seqs or "R" not in seqs:
            rows.append({"sample": sample, "status": "MISSING_PAIR"})
            continue

        f_len, r_len = len(seqs["F"]), len(seqs["R"])
        result = merge_pair(seqs["F"], seqs["R"])
        if result is None:
            rows.append({"sample": sample, "f_len": f_len, "r_len": r_len,
                         "status": "NO_HIT"})
            continue

        out_seq = result["merged"] if union else result["overlap"]
        status = "LOW_IDENTITY" if result["identity"] < MIN_IDENTITY else "OK"

        row = {
            "sample": sample,
            "f_len": f_len,
            "r_len": r_len,
            "overlap_len": result["overlap_len"],
            "output_len": len(out_seq),
            "f_r_identity": result["identity"],
            "status": status,
        }
        # realign each raw read back to the output sequence
        for label, raw in (("f", seqs["F"]), ("r", revcomp(seqs["R"]))):
            st = align_to_reference(raw, out_seq)
            if st is None:
                continue
            row[f"{label}_realign_identity"] = st["identity"]
            row[f"{label}_query_cov"] = st["query_coverage"]
            row[f"{label}_target_cov"] = st["target_coverage"]
        rows.append(row)

        with open(os.path.join(outdir, f"{sample}{suffix}"), "w", encoding="utf8") as fh:
            fh.write(f">{sample}\n{out_seq}\n")
        merged_lines.append(f">{sample}\n{out_seq}\n")

    with open(sum_fa, "w", encoding="utf8") as fh:
        fh.writelines(merged_lines)

    return by_sample, rows, mode


def fmt(v):
    if isinstance(v, float):
        return f"{v:.4f}"
    return "" if v is None else str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_dir", help="directory containing .seq files")
    ap.add_argument("--outdir", default=None,
                    help="output dir (default <input_dir>/merged_output)")
    ap.add_argument("--union", action="store_true",
                    help="emit full-length union instead of overlap intersect")
    args = ap.parse_args()

    outdir = args.outdir or os.path.join(args.input_dir, "merged_output")
    by_sample, rows, mode = run(args.input_dir, outdir, args.union)

    stats_path = os.path.join(outdir, "merge_stats.tsv")
    with open(stats_path, "w", encoding="utf8") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c)) for c in COLS) + "\n")

    n = len(by_sample)
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    summary = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"mode={mode} samples={n} {summary}")
    print(f"stats: {stats_path}")
    for r in rows:
        if r["status"] != "OK":
            print(f"  {r['status']}: {r['sample']}")


if __name__ == "__main__":
    main()
