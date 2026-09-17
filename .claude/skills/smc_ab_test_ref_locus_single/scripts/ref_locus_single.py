#!/usr/bin/env python3
"""单 reference 的 per-locus A/B 准确率分析 (无 barcode 循环)。

针对 "一个 run 只有一个 reference" 的场景: 两侧的 query 各自汇总成**单个**
查询 (一个 FASTQ/BAM, 或把一个 demux 目录里的 *.fastq* 拼接成一个 pooled FASTQ),
对齐到**同一个** reference, 逐位点比较两组。

底层复用 gseda.ab_analysis.locus_error_rate.main_cli: 它对 (一组 control query,
一组 exp query, 一个 ref) 做 一次 gsmm2 比对 + gsetl 位点表 + 位点准确率列 +
按 (refname,pos) join, 本脚本不重写任何比对逻辑, 只做两件小事:
  1. resolve_query: 把每侧的 QUERY (文件或 demux 目录) 归一成单个 query 路径;
  2. 在产出的 locus_accuracy_joined.csv 上写一个 depth-weighted A/B summary。

IMPORTANT: query 是 FASTQ 时不带 CCS 的 rq/np 字段, gsmm2 的 --rq-range/--np-range
被忽略, 指标是**相对 Sanger reference 的对齐 identity**, 不是原始 read 准确率。
"""
import argparse
import glob
import os
import sys

GS = "/root/miniconda3/envs/gseda/bin/python"
REPO_SRC = "/root/projects/gsda/third_party/gseda/src"
if REPO_SRC not in sys.path:
    sys.path.insert(0, REPO_SRC)

import polars as pl  # noqa: E402
from gseda.ab_analysis.locus_error_rate import main_cli as ler_main  # noqa: E402
from gseda.fact_table_ana.polars_init import polars_env_init  # noqa: E402


def log(*a):
    print(*a, flush=True)


def parse_label_query(s, name):
    """'LABEL=QUERY' -> (label, query)."""
    if "=" not in s:
        raise SystemExit(f"--{name} must be LABEL=QUERY, got: {s}")
    label, q = s.split("=", 1)
    if not label:
        raise SystemExit(f"--{name} label is empty: {s}")
    return label, q


def _check_bam_unmapped(bam: str):
    """校验 BAM query 的**所有 read 都是 unmapped**。

    gsmm2 负责把 query 对齐到 --ref; 因此喂进来的 BAM 必须是**未比对**的
    (所有 read 都带 UNMAP/0x4 标志)。若已有 read 是 mapped, 说明输入是已经
    对齐过的 BAM, 再对齐会拿到错误的 identity。这里抽样检查, 一旦发现
    mapped read 就直接报错退出, 避免静默产出错误指标。
    """
    import pysam
    n = 0
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(until_eof=True):
            if not (r.flag & 0x4):  # 0x4 = UNMAP bit
                raise SystemExit(
                    f"BAM query {bam} 含 mapped read (第 {n+1} 条, qname={r.query_name}, "
                    f"flag={r.flag}); 该 skill 要求 query BAM 内所有 read 都是 unmapped。"
                    f"请提供未比对的 BAM。")
            n += 1
    log(f"checked {bam}: all {n} reads unmapped OK")


def resolve_query(label: str, query: str, outdir: str) -> str:
    """把一侧的 QUERY 归成单个 query 路径。

    - 文件 (fastq/fq/fastq.gz/bam...) -> 原样返回。
    - 目录 -> 把其中所有 *.fastq* 按文件名排序拼接成一个 pooled FASTQ,
      写到 <outdir>/.pool_<label>.fastq, 返回该 pooled 路径。
    """
    if os.path.isfile(query):
        if query.endswith(".bam") and not query.endswith(".bam.bai"):
            _check_bam_unmapped(query)
        return query
    if os.path.isdir(query):
        files = sorted(glob.glob(os.path.join(query, "*.fastq*")))
        if not files:
            raise SystemExit(f"query dir has no *.fastq* files: {query}")
        pooled = os.path.join(outdir, f".pool_{label}.fastq")
        n = 0
        with open(pooled, "wb") as out:
            for f in files:
                with open(f, "rb") as fh:
                    for line in fh:
                        out.write(line)
                        if line.startswith(b"@") or line.startswith(b">"):
                            n += 1
        log(f"pooled {len(files)} files from {query} -> {pooled} ({n} reads)")
        return pooled
    raise SystemExit(f"query is neither a file nor a dir: {query}")


def depth_weighted(sub: pl.DataFrame, tag: str):
    """Depth-weighted locus accuracy eq/(eq+diff+ins+del) for one side."""
    t = sub["eq_" + tag] + sub["diff_" + tag] + sub["ins_" + tag] + sub["del_" + tag]
    tot = t.sum()
    return sub["eq_" + tag].sum() / tot if tot else None


def write_summary(joined_path: str, outdir: str, lc: str, lx: str, ref: str):
    j = pl.read_csv(joined_path, separator="\t")
    c_all = depth_weighted(j, "ctrl")
    e_all = depth_weighted(j, "exp")

    md = [f"# {lc} vs {lx} — single-reference per-locus accuracy A/B\n"]
    md.append(f"reference: `{os.path.basename(ref)}`  ({j.height} loci)\n")
    md.append("depth-weighted locus_accuracy = eq/(eq+diff+ins+del) over the "
              "reference, per (refname, pos) row.\n")
    md.append("metric: alignment identity vs the Sanger reference (FASTQ input → "
              "no rq/np filtering applies).\n")
    md.append("## Overall (all loci pooled)\n")
    md.append("| side | locus_accuracy |")
    md.append("|---|---:|")
    md.append(f"| {lc} (control) | {c_all:.6f} |")
    md.append(f"| {lx} (experiment) | {e_all:.6f} |")
    md.append(f"| Δ (exp − ctrl) | {e_all - c_all:+.6f} |\n")
    md.append("")
    with open(os.path.join(outdir, "summary.md"), "w") as f:
        f.write("\n".join(md))
    log(f"overall: {lc}={c_all:.5f}  {lx}={e_all:.5f}  Δ={e_all - c_all:+.5f}")
    log(f"wrote summary.md")


def main():
    ap = argparse.ArgumentParser(
        description="Single-reference per-locus A/B of alignment accuracy via "
                    "gseda.ab_analysis.locus_error_rate (no barcode loop). Each "
                    "side is one query (a FASTQ/BAM, or a demux dir pooled into "
                    "one FASTQ); both align to the same --ref.")
    ap.add_argument("--control", required=True, metavar="LABEL=QUERY",
                    help="control: label=query (file or dir-of-*.fastq* to pool)")
    ap.add_argument("--exp", required=True, metavar="LABEL=QUERY",
                    help="experiment: label=query (file or dir-of-*.fastq* to pool)")
    ap.add_argument("--ref", required=True, help="single reference fasta")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rq-thr", default="0",
                    help="read-accuracy thr for locus_error_rate (1 val; 2 vals "
                         "'[ctrl,exp]' also accepted). FASTQ has no rq field so "
                         "this is ignored by gsmm2; must still be supplied. "
                         "Default 0 (non-filtering).")
    ap.add_argument("--np-thr", default="5",
                    help="number-of-passes lower bound. Ignored for FASTQ (no np "
                         "field). Default 5.")
    ap.add_argument("--threads", type=int, default=None,
                    help="gsmm2 thread count (default: CPU count)")
    args = ap.parse_args()

    polars_env_init()

    lc, qc = parse_label_query(args.control, "control")
    lx, qx = parse_label_query(args.exp, "exp")
    out = args.outdir
    os.makedirs(out, exist_ok=True)

    log(f"control     {lc} = {qc}")
    log(f"experiment  {lx} = {qx}")
    log(f"ref: {args.ref}")
    log(f"rq-thr: {args.rq_thr}   np-thr: {args.np_thr}  (ignored: FASTQ has no rq/np field)")

    cq = resolve_query(lc, qc, out)
    eq = resolve_query(lx, qx, out)

    ler_main([
        "--control", cq,
        "--exp", eq,
        "--ref", args.ref,
        "--rq-thr", args.rq_thr,
        "--np-thr", args.np_thr,
        "--outdir", out,
    ] + (["--threads", str(args.threads)] if args.threads else []))

    joined_path = os.path.join(out, "locus_accuracy_joined.csv")
    if not os.path.exists(joined_path):
        log("WARN locus_accuracy_joined.csv missing; skipping summary")
        return
    write_summary(joined_path, out, lc, lx, args.ref)
    log("ALL DONE")


if __name__ == "__main__":
    main()
