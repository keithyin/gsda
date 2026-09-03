"""两个 barcode-basecalling 模型（A/B）相对 consensus 的 identity 评估。

对比两套 barcode 拆分后的 CCS reads（例如 baseline 模型 vs 优化后模型），
把每个 barcode 的 reads 比对到该 barcode 对应的 consensus 序列，用
gseda.ppl.sequencing_report_v2 的口径（gsmm2-aligned-metric）统计 identity，
从而评估哪个模型的 basecalling 更准。

配对规则（默认）：
    run 目录下  Adaptor-barcodeXXX-Y.fastq
    consensus 目录  {prefix}Adaptor-barcodeXXX-Y{suffix}
    即 consensus 文件 = prefix + <fastq 去掉 .fastq> + suffix
    （本仓库场景 prefix=Group_0_，suffix=.consensus.fasta）

只有「A、B 两个 run 都有该 barcode」且「consensus 序列非空」的 barcode 才纳入比较，
保证严格的逐 barcode 一一对应（like-for-like）。

流程（对每个纳入的 barcode）：
    1. fastq -> bam（gseda.file_format_cvt.fastq2bam，自动写入 rq tag，
       ref = 该 barcode 的 consensus，供 --rq-range 过滤使用）
    2. sequencing_report_v2.main() 跑 gsmm2-aligned-metric（带 --rq-range）
    3. 汇总：
       - 池化指标：把该 run 所有 barcode 的 fact csv 交给 align_stats()
         （复用 report 自身的聚合逻辑，口径与单 barcode 完全一致）
       - 逐 barcode identity：读各 aggr csv 的 identity 行
       - 逐 barcode 匹配比较（A vs B），统计 B 胜出的 barcode 数、平均/中位差

输出（--outdir）：
    - per_barcode_identity.csv   barcode, a_identity, b_identity, diff(B-A)
    - summary.txt                池化指标 + 逐 barcode 比较的文本报告
    - stdout                     同 summary.txt

用法示例：
    python ab_consensus_identity.py \
        --run-a /path/to/baseline/demuxed --label-a baseline \
        --run-b /path/to/new_model/demuxed --label-b new_model \
        --consensus /path/to/Consensus/Sequences \
        --rq-range 0.99:1.1 \
        --outdir /tmp/ab_eval
"""

import os
import sys
import glob
import argparse
import statistics
import multiprocessing as mp

import polars as pl

cur_path = os.path.abspath(__file__)
cur_dir = os.path.dirname(cur_path)
# cur_dir = .../gseda/src/gseda/16s
gseda_pkg_dir = os.path.dirname(cur_dir)          # .../gseda/src/gseda
src_dir = os.path.dirname(gseda_pkg_dir)          # .../gseda/src
sys.path.insert(0, src_dir)

# 复用仓库内已有的两个模块，保证口径与人工流程一致
from gseda.file_format_cvt.fastq2bam import fastx_to_bam            # noqa: E402
from gseda.ppl import sequencing_report_v2 as seq_rep               # noqa: E402


def consensus_has_seq(consensus_path: str) -> bool:
    """consensus 文件是否含至少一条序列（去掉 header 行后有内容）。"""
    if not os.path.exists(consensus_path) or os.path.getsize(consensus_path) == 0:
        return False
    with open(consensus_path) as fh:
        for line in fh:
            if not line.startswith(">") and line.strip():
                return True
    return False


def discover_barcodes(run_a_dir, run_b_dir, cons_dir, prefix, suffix):
    """返回 list of dict: {stem, fastq_a, fastq_b, consensus}。

    只保留 A、B 都有 fastq 且 consensus 序列非空的 barcode。
    """
    def stems(d):
        return {os.path.basename(f)[:-len(".fastq")]
                for f in glob.glob(os.path.join(d, "*.fastq"))}

    common = stems(run_a_dir) & stems(run_b_dir)
    pairs = []
    for stem in sorted(common):
        cons = os.path.join(cons_dir, prefix + stem + suffix)
        if not consensus_has_seq(cons):
            continue
        pairs.append({
            "stem": stem,
            "fastq_a": os.path.join(run_a_dir, stem + ".fastq"),
            "fastq_b": os.path.join(run_b_dir, stem + ".fastq"),
            "consensus": cons,
        })
    return pairs


def _convert(pair, workdir, label, threads):
    """fastq -> bam（幂等：已存在则跳过）。pair 为单个 barcode dict。"""
    os.makedirs(workdir, exist_ok=True)
    fastq = pair["fastq_a"] if label == "a" else pair["fastq_b"]
    out = os.path.join(workdir, pair["stem"] + ".bam")
    if not (os.path.exists(out) and os.path.getsize(out) > 0):
        fastx_to_bam(fastq, out, ref_fasta=pair["consensus"], threads=threads)
    return out


def _report(bam, cons, outdir, rq_range, threads):
    """对单个 bam 跑 sequencing_report_v2，返回 aggr csv 路径。"""
    os.makedirs(outdir, exist_ok=True)
    aggr, fact, _basic = seq_rep.main(
        bam_file=bam, ref_fa=cons, outdir=outdir,
        force=True, threads=threads, rq_range=rq_range,
    )
    return aggr, fact


def pooled_metrics(fact_files):
    """把多个 fact csv 交给 report 自身的 align_stats，取池化指标。

    返回 dict：identity / identity-p50 / mmRate / NHInsRate / homoInsRate /
               NHDelRate / homoDelRate / queryCoverage3 / alignedRatio /
               aligned / notAligned
    """
    if not fact_files:
        return {}
    df = seq_rep.align_stats(fact_files)  # DataFrame[name, value]
    return dict(zip(df["name"].to_list(), df["value"].to_list()))


def per_barcode_identity(aggr_files):
    """从各 aggr csv 提取 identity 行 -> {stem: identity}。"""
    out = {}
    for f in aggr_files:
        stem = os.path.basename(f)[: -len(".gsmm2_aligned_metric_aggr.csv")]
        d = pl.read_csv(f, separator="\t")
        row = d.filter(pl.col("name") == "identity")
        if row.height:
            out[stem] = row["value"][0]
    return out


def _f(v, fmt="{:.6f}"):
    try:
        return fmt.format(v)
    except (TypeError, ValueError):
        return str(v)


def evaluate(run_a_dir, run_b_dir, cons_dir, outdir, prefix, suffix,
             rq_range, label_a, label_b, parallel):
    pairs = discover_barcodes(run_a_dir, run_b_dir, cons_dir, prefix, suffix)
    if not pairs:
        print("ERROR: 没有找到可配对的 barcode（A、B 都有 fastq 且 consensus 非空）。")
        print(f"       run_a={run_a_dir}")
        print(f"       run_b={run_b_dir}")
        print(f"       consensus={cons_dir}  prefix={prefix!r} suffix={suffix!r}")
        return 1
    print(f"纳入比较的 barcode: {len(pairs)}")

    work_a = os.path.join(outdir, f"{label_a}_bams")
    work_b = os.path.join(outdir, f"{label_b}_bams")

    # 每进程 report 线程数：并行时均分 CPU，避免 gsmm2 过度超订
    per_threads = max(1, mp.cpu_count() // max(1, parallel))

    # --- step1: fastq -> bam（并行）---
    with mp.Pool(parallel) as pool:
        pool.starmap(_convert, [(p, work_b, "b", 0) for p in pairs], chunksize=1)
        pool.starmap(_convert, [(p, work_a, "a", 0) for p in pairs], chunksize=1)

    # --- step2: 逐 barcode 跑 report（并行，幂等）---
    def jobs(label, workdir):
        return [(os.path.join(workdir, p["stem"] + ".bam"),
                 p["consensus"],
                 os.path.join(workdir, p["stem"] + "-metric"))
                for p in pairs]

    with mp.Pool(parallel) as pool:
        res_a = pool.starmap(
            _report,
            [(j[0], j[1], j[2], rq_range, per_threads) for j in jobs(label_a, work_a)],
            chunksize=1)
        res_b = pool.starmap(
            _report,
            [(j[0], j[1], j[2], rq_range, per_threads) for j in jobs(label_b, work_b)],
            chunksize=1)

    fact_a = [r[1] for r in res_a if r[1]]
    fact_b = [r[1] for r in res_b if r[1]]
    aggr_a_files = [r[0] for r in res_a if r[0]]
    aggr_b_files = [r[0] for r in res_b if r[0]]

    # --- step3: 汇总 ---
    pool_a = pooled_metrics(fact_a)
    pool_b = pooled_metrics(fact_b)
    pb_a = per_barcode_identity(aggr_a_files)
    pb_b = per_barcode_identity(aggr_b_files)
    common = sorted(set(pb_a) & set(pb_b))

    # 逐 barcode 差异
    diffs = [(s, pb_a[s], pb_b[s], pb_b[s] - pb_a[s]) for s in common]
    win_b = sum(1 for d in diffs if d[3] > 0)
    win_a = sum(1 for d in diffs if d[3] < 0)
    tie = len(diffs) - win_b - win_a
    dvals = [d[3] for d in diffs]

    # 写 per-barcode CSV
    csv_path = os.path.join(outdir, "per_barcode_identity.csv")
    with open(csv_path, "w") as fh:
        fh.write(f"barcode\t{label_a}_identity\t{label_b}_identity\tdiff_{label_b_minus(label_a, label_b)}\n")
        for s, ia, ib, dd in diffs:
            fh.write(f"{s}\t{_f(ia)}\t{_f(ib)}\t{_f(dd, '{:+.6f}')}\n")

    # 组 summary
    def pool_line(name):
        va, vb = pool_a.get(name, 0), pool_b.get(name, 0)
        try:
            d = vb - va
            return f"  {name:<22} {label_a}={_f(va)}  {label_b}={_f(vb)}  diff={_f(d, '{:+.6f}')}"
        except Exception:
            return f"  {name:<22} {label_a}={va}  {label_b}={vb}"

    lines = []
    lines.append("=" * 72)
    lines.append(f"A/B consensus identity 评估  (rq-range={rq_range})")
    lines.append(f"barcode 数: {len(pairs)}   可比较: {len(common)}")
    lines.append(f"  {label_a}: {run_a_dir}")
    lines.append(f"  {label_b}: {run_b_dir}")
    lines.append("=" * 72)
    lines.append("[池化指标 —— 各 run 所有比对上 reads 汇总]")
    for k in ["aligned", "alignedRatio", "identity", "identity-p50",
              "mmRate", "NHInsRate", "homoInsRate", "NHDelRate", "homoDelRate",
              "queryCoverage3"]:
        lines.append(pool_line(k))
    # insertion / deletion 合计
    for k, ia, ib in [("insRate", "NHInsRate", "homoInsRate"),
                      ("delRate", "NHDelRate", "homoDelRate")]:
        va = pool_a.get(ia, 0) + pool_a.get(ib, 0)
        vb = pool_b.get(ia, 0) + pool_b.get(ib, 0)
        lines.append(f"  {k:<22} {label_a}={_f(va)}  {label_b}={_f(vb)}  diff={_f(vb - va, '{:+.6f}')}")
    lines.append("")
    lines.append(f"[逐 barcode 比较 —— {label_b} vs {label_a}]")
    lines.append(f"  {label_b} 更高: {win_b}   {label_a} 更高: {win_a}   相等: {tie}")
    if dvals:
        lines.append(f"  平均 diff({label_b}-{label_a}) = {_f(statistics.mean(dvals), '{:+.6f}')}")
        lines.append(f"  中位 diff({label_b}-{label_a}) = {_f(statistics.median(dvals), '{:+.6f}')}")
    summary = "\n".join(lines)

    print(summary)
    with open(os.path.join(outdir, "summary.txt"), "w") as fh:
        fh.write(summary + "\n")
    print(f"\n[输出] {csv_path}")
    print(f"[输出] {os.path.join(outdir, 'summary.txt')}")
    return 0


def label_b_minus(a, b):
    return f"{b}_minus_{a}"


def main():
    ap = argparse.ArgumentParser(description="A/B barcode basecalling 模型对 consensus 的 identity 评估")
    ap.add_argument("--run-a", required=True, help="模型 A 的 demuxed fastq 目录")
    ap.add_argument("--run-b", required=True, help="模型 B 的 demuxed fastq 目录")
    ap.add_argument("--consensus", required=True, help="consensus 序列目录")
    ap.add_argument("--label-a", default="A", help="A 的标签（默认 A）")
    ap.add_argument("--label-b", default="B", help="B 的标签（默认 B）")
    ap.add_argument("--prefix", default="Group_0_",
                    help="consensus 文件名前缀（默认 Group_0_）")
    ap.add_argument("--suffix", default=".consensus.fasta",
                    help="consensus 文件名后缀（默认 .consensus.fasta）")
    ap.add_argument("--rq-range", default=None,
                    help="gsmm2 --rq-range，如 0.99:1.1")
    ap.add_argument("--outdir", required=True, help="输出目录")
    ap.add_argument("--parallel", type=int, default=16,
                    help="并行度（默认 16）")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sys.exit(evaluate(
        args.run_a, args.run_b, args.consensus, args.outdir,
        args.prefix, args.suffix, args.rq_range,
        args.label_a, args.label_b, args.parallel,
    ))


if __name__ == "__main__":
    main()
