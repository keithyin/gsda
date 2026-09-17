#!/usr/bin/env python
"""Aggregate per-pair rq_iy_analysis.py outputs into one summary table + markdown report."""
import argparse
import math
import pathlib
import subprocess
import sys

import polars as pl

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from run_rq_iy_batch import label_of  # noqa: E402

ROOT = pathlib.Path("/root/projects/gsda")
TUPLES = ROOT / "prototype/2026/2026Q3/STR_optimization/str_datas_tuples.txt"
OUT_DIR = ROOT / "prototype/2026/2026Q3/STR_optimization"


def q2phreq(v):
    return -10 * math.log10(max(1 - v, 1e-12))


def count_bam(path, extra=()):
    r = subprocess.run(["samtools", "view", "-c", *extra, path],
                       capture_output=True, text=True)
    return int(r.stdout) if r.returncode == 0 else None


def fasta_max_len(path):
    best = cur = 0
    for line in open(path):
        if line.startswith(">"):
            best, cur = max(best, cur), 0
        else:
            cur += len(line.strip())
    return max(best, cur)


def stat_pair(csv_path, ref_path):
    df = pl.read_csv(csv_path, separator="\t", schema_overrides={"iy": pl.Float64})
    df = df.with_columns(
        pl.when(pl.col("rq") > 0.99999).then(0.99999).otherwise(pl.col("rq")).alias("rq"),
        pl.when(pl.col("iy") > 0.99999).then(0.99999).otherwise(pl.col("iy")).alias("iy"),
    ).with_columns(
        pl.col("rq").map_elements(q2phreq, return_dtype=pl.Float64).alias("prq"),
        pl.col("iy").map_elements(q2phreq, return_dtype=pl.Float64).alias("piy"),
    )

    def both(t):
        th = 1 - 10 ** (t / -10)
        return df.filter((pl.col("rq") >= th) & (pl.col("iy") >= th)).height

    def acc(t):  # of reads predicted >=Q(t), how many really are
        th = 1 - 10 ** (t / -10)
        n = df.filter(pl.col("rq") >= th).height
        return both(t) / n if n else None

    def rec(t):  # of reads really >=Q(t), how many predicted so
        th = 1 - 10 ** (t / -10)
        n = df.filter(pl.col("iy") >= th).height
        return both(t) / n if n else None

    return dict(
        n_rows=df.height,
        n_reads=df.select(pl.col("qname").n_unique()).item(),
        iy_mean=df["iy"].mean(),
        piy_mean=df["piy"].mean(),
        prq_mean=df["prq"].mean(),
        MAE=(df["rq"] - df["iy"]).abs().mean(),
        reflen=fasta_max_len(ref_path),
        cov_ref=df.select((pl.col("rend") - pl.col("rstart") + 1).max()).item(),
        aln_len=df.select((pl.col("qend") - pl.col("qstart") + 1).mean()).item(),
        qlen=df["qlen"].mean(),
        A20=acc(20), R20=rec(20), A30=acc(30), R30=rec(30),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-total", action="store_true",
                    help="skip samtools counting of query BAMs")
    ap.add_argument("--tuples", default=str(TUPLES))
    ap.add_argument("--out-prefix", default=str(OUT_DIR / "rq_iy_summary"))
    args = ap.parse_args()

    ns = {}
    exec(open(args.tuples).read(), ns)

    rows = []
    for t in ns["datas"]:
        smc, ref = t[1], t[2]
        batch, run, bc, sample = label_of(smc, ref).split("__")
        gsetl = f"{pathlib.Path(smc).parent}/{pathlib.Path(smc).stem}.aligned-gsetl"
        csv_path = f"{gsetl}/fact_aligned_bam_bam_basic.csv"
        row = dict(batch=batch, run=run, barcode=bc, sample=sample,
                   ref=pathlib.Path(ref).name, status="ok")
        if not pathlib.Path(csv_path).exists():
            row.update(status="NO_OUTPUT", n_reads=0)
            rows.append(row)
            continue
        if not args.skip_total:
            row["query_total"] = count_bam(smc)
        row.update(stat_pair(csv_path, ref))
        rows.append(row)

    cols = ["batch", "R", "BC", "sample", "query_total", "n_reads", "kept",
            "iy_mean", "piy_mean", "prq_mean", "A20", "R20", "A30", "R30",
            "aln_len", "qlen", "cov_ref", "reflen", "refcov", "run", "ref"]
    df0 = (
        pl.DataFrame(rows)
        .with_columns([
            pl.col("run").str.extract(r"Run(\d+)$").alias("R"),
            pl.col("barcode").str.extract(r"Barcode(\d+)$").alias("BC"),
        ])
        .with_columns([
            (pl.col("n_reads") / pl.col("query_total")).alias("kept"),
            (pl.col("cov_ref") / pl.col("reflen")).alias("refcov"),
        ])
    )
    out = (
        df0.with_columns([pl.col(c).round(4) for c, t in df0.schema.items()
                          if t == pl.Float64])
        .sort(["batch", "R", "BC"])
    )
    out.select(cols).write_csv(pathlib.Path(args.out_prefix).with_suffix(".tsv"),
                               separator="\t")

    # ---- markdown report
    lines = ["| batch | Run | BC | sample | reads(kept/total) | %kept | mean iy | iy(Q) | "
             "cov/len(ref) | A@Q20 | R@Q20 | A@Q30 | R@Q30 |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]

    def f(v, n=3):
        return "-" if v is None else f"{v:.{n}f}"

    for r in out.iter_rows(named=True):
        cov = "-" if r["cov_ref"] is None or r["reflen"] is None else \
            f"{int(r['cov_ref'])}/{int(r['reflen'])} ({r['refcov']*100:.1f}%)"
        lines.append(
            f"| {r['batch'].replace('-batch-of-data','')} | {r['R']} | {r['BC']} "
            f"| {r['sample']} | {r['n_reads']}/{r['query_total']} "
            f"| {f(r['kept']*100 if r['kept'] is not None else None, 1)} "
            f"| {f(r['iy_mean'])} | {f(r['piy_mean'], 1)} | {cov} "
            f"| {f(r['A20'])} | {f(r['R20'])} | {f(r['A30'])} | {f(r['R30'])} |")
    pathlib.Path(args.out_prefix).with_suffix(".md").write_text("\n".join(lines) + "\n")

    pl.Config.set_tbl_formatting("UTF8_FULL_CONDENSED")
    pl.Config.set_tbl_rows(100)
    pl.Config.set_tbl_cols(40)
    print(out.select(cols))


if __name__ == "__main__":
    main()
