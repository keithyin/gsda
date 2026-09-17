#!/usr/bin/env python
"""Run rq_iy_analysis.py over every (query smc-bam, ref) pair in str_datas_tuples.txt.

tuple layout: (demux_bam, smc_bam[query], ref_fa, x, y)
"""
import argparse
import concurrent.futures
import os
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path("/root/projects/gsda")
PY38 = "/root/miniconda3/envs/py38/bin/python"
SCRIPT = str(ROOT / "third_party/gseda/src/gseda/ppl/rq_iy_analysis.py")
TUPLES = ROOT / "prototype/2026/2026Q3/STR_optimization/str_datas_tuples.txt"
LOG_DIR = ROOT / "prototype/2026/2026Q3/STR_optimization/rq_iy_logs"


def label_of(smc_bam, ref):
    parts = pathlib.Path(smc_bam).parts
    batch = next(p for p in parts if p.endswith("-batch-of-data"))
    run = parts[parts.index(batch) + 1]
    bc = re.search(r"Barcode\d+", pathlib.Path(smc_bam).stem).group(0)
    sample = pathlib.Path(ref).name
    sample = re.sub(r"(\.intersect)?\.fa(sta)?$", "", sample)
    sample = re.sub(r"^\d{6}STR-?", "", sample)
    sample = re.sub(r"^STR", "", sample)
    if ".with_di.q20." in smc_bam:      # 同一 barcode 的 di 补基变体，标出来以免两行撞名
        sample += "~q20"
    return f"{batch}__{run}__{bc}__{sample}"


def run_one(label, smc_bam, ref, rq_range, log_dir=LOG_DIR):
    env = dict(os.environ, MPLBACKEND="Agg")
    log = log_dir / f"{label}.log"
    cmd = [PY38, SCRIPT, "--smc-bam", smc_bam, "--ref", ref]
    if rq_range:
        cmd += ["--rq-range", rq_range]
    with open(log, "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    return label, smc_bam, ref, rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tuples", default=str(TUPLES))
    ap.add_argument("--log-dir", default=str(LOG_DIR))
    ap.add_argument("--rq-range", default="0.99:1.1")
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--only", default=None, help="substring filter on label")
    args = ap.parse_args()

    ns = {}
    exec(open(args.tuples).read(), ns)
    jobs = []
    for t in ns["datas"]:
        lab = label_of(t[1], t[2])
        if args.only and args.only not in lab:
            continue
        jobs.append((lab, t[1], t[2]))

    log_dir = pathlib.Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    print(f"{len(jobs)} pairs, {args.jobs} workers", flush=True)
    fails = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_one, *j, args.rq_range, log_dir) for j in jobs]
        for n, f in enumerate(concurrent.futures.as_completed(futs), 1):
            lab, _, _, rc = f.result()
            print(f"[{n}/{len(jobs)}] rc={rc} {lab}", flush=True)
            if rc != 0:
                fails.append(lab)
    print("FAILED:", fails or "none")


if __name__ == "__main__":
    main()
