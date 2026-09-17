#!/usr/bin/env python3
"""
Merge per-run barcode_ref_align summary.tsv files into one cross-run table.

Each run dir produced by the barcode_ref_align skill holds a summary.tsv
(sample = "BarcodeNN_<plasmid>") plus refs/<plasmid>.fa. This joins those with
the run-aware mapping (plasmid<TAB>barcode<TAB>RUN号, barcode as "24标签-N") and
the barcode-2-barcodename table, so a row reads:

    run  plasmid  24标签  barcodename  barcode_orig  fastq  ref  ref_len_bp  <metrics...>

so every alignment number can be traced back to a concrete FASTQ and reference.
"""
import argparse
import csv
import glob
import os
import re
import sys

BARCODE2NAME_TSV = "/data1/ccs_data/str-optimization/barcode-2-barcodename.tsv"
META = ["run", "plasmid", "24标签", "barcodename", "barcode_orig", "fastq",
        "ref", "ref_len_bp"]
RUN_SHORT = {"20260805_250302Y0001_Run0007": "250302Y0001/Run0007",
             "20260805_250804Y0004_Run0001": "250804Y0004/Run0001",
             "20260805_250804Y0004_Run0002": "250804Y0004/Run0002"}
MD_HEAD = ("| 样本 | 24标签 | barcode | 参考长度bp | reads | 比对率% | query覆盖 | "
           "identity% | identity-p50 | mmRate | identity≥0.99% | 自身参考最优% |")
MD_SEP = "|---|---|---|---|---|---|---|---|---|---|---|---|"


def pct(row, key, digits=1):
    try:
        return f"{100 * float(row[key]):.{digits}f}"
    except (KeyError, TypeError, ValueError):
        return ""


def frac(row, key, digits=3):
    try:
        return f"{float(row[key]):.{digits}f}"
    except (KeyError, TypeError, ValueError):
        return ""


def write_markdown(rows, selfpref, path):
    """Per-run metric tables plus the cross-map self-reference preference."""
    def sort_key(r):
        return int(float(r["reads_num"])) if r["reads_num"] else 0
    with open(path, "w") as fh:
        fh.write("# 第二批 barcode reads × 一代参考比对（gsmm2 sequencing_report_v2）\n\n")
        fh.write("- 输入：三个 run 的 query 统一取 `barcode_assign/`（Run0007 已是 "
                 "`BarcodeNN.fastq`；Run0001/Run0002 是 `Adaptor-barcodeNNN-M.fastq`，"
                 "经 `barcode-2-barcodename.tsv` 解析为 `BarcodeNN.fastq` 软链）；"
                 "ref = `STR第二批一代测序/STR第二批一代测序/merged_output/merged.fa`"
                 "（`--ref-prefix STR`）；对应关系 = `plasmid_2_barcode.tsv`"
                 "（`24标签-N` → `barcode-NN` → `Adaptor-*`）\n")
        fh.write("- 参数：`--rq-range 0.99:1.1 --short-aln 1`，conda env `py38`；44 组全部成功，"
                 "无一 SKIP（29 条参考覆盖全部质粒）\n")
        fh.write("- 校验：`barcode_assign/` 即早先 `barcodes_reads_fastq_amplicon/` 的重命名落位；"
                 "`merged.fa` 的 29 条记录与同目录单样本 `STRxx-x.fa` 序列全等，"
                 "`refs/` 抽出的 29 条参考与 `merged.fa` 逐条全等，44 组无 SKIP；"
                 "各组 reads 数与 identity 均值与早先记录（9/4、9/7 两次）一致，"
                 "说明原始 query 数据未变。"
                 "Run0002 未被点名，但 `plasmid_2_barcode.tsv` 覆盖它且 query 目录同在，一并算出\n")
        fh.write("- 明细：`combined_summary.tsv`（44 行，含 fastq / ref 绝对路径）；"
                 "每组指标目录 `<run>/BarcodeNN_<plasmid>-metric/`\n")
        fh.write("- `自身参考最优%`：read 同时比对 29 条参考时，最佳命中仍为「自己那组参考」的"
                 "比例，见 `cross_map_diagnostic.txt`"
                 "（`--max-reads 1500`，与 9/7 03:50 那次同参数）；"
                 "另一套（minimap2、按目标/背景 read 分开）的统计在 "
                 "`prototype/2026/2026Q3/STR_optimization/ccs比对统计报告.md`\n")
        for run in sorted({r["run"] for r in rows}):
            rs = [r for r in rows if r["run"] == run]
            fh.write(f"\n## {RUN_SHORT.get(run, run)}（{len(rs)} 组）\n\n")
            fh.write(MD_HEAD + "\n" + MD_SEP + "\n")
            for r in sorted(rs, key=sort_key, reverse=True):
                sp = selfpref.get((run, r["plasmid"]), "")
                fh.write("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n".format(
                    r["plasmid"], r["24标签"], r["barcodename"], r["ref_len_bp"],
                    f"{int(float(r['reads_num'])):,}" if r["reads_num"] else "",
                    pct(r, "alignedRatio"), frac(r, "queryCoverage"),
                    pct(r, "identity", 3), frac(r, "identity-p50", 4),
                    frac(r, "mmRate", 5), pct(r, "identity≥0.99"),
                    f"{sp:.1f}" if sp != "" else ""))
            ids = [float(r["identity"]) for r in rs if r["identity"]]
            al = [float(r["alignedRatio"]) for r in rs if r["alignedRatio"]]
            lo = min(rs, key=lambda r: float(r["identity"] or 1))
            hi = max(rs, key=lambda r: float(r["identity"] or 0))
            fh.write(f"\n- 本 run 汇总：identity mean={100*sum(ids)/len(ids):.3f}% "
                     f"min={100*min(ids):.3f}% max={100*max(ids):.3f}%；"
                     f"比对率 mean={100*sum(al)/len(al):.1f}% min={100*min(al):.1f}%\n")
            fh.write(f"- identity 最低：{lo['plasmid']}（{100*float(lo['identity']):.3f}%）；"
                     f"最高：{hi['plasmid']}（{100*float(hi['identity']):.3f}%）\n")
        fh.write(ANALYSIS)


ANALYSIS = """
## 读数须知：alignedRatio / queryCoverage 不能直接当质量看

- 每组数据混合两个物种（见 `ccs比对统计报告.md`）：短 read（~600-1100bp）= 目标扩增子，完整比对参考；
  长 read（~2.7kb）= 环状质粒全长背景分子，仅与参考 3' 端 ~48bp 公共载体 MCS 区匹配。
- gsmm2 把这 48bp 也算「比对上」，所以 alignedRatio 普遍 0.94-0.99 并不代表数据好；
  同理 queryCoverage 被长 read 稀释（长 read 只覆盖自身 ~1.8%，各组 queryCoverage-p25≈0.0176）。
  判读优先级：`identity` / `identity≥0.99%` / `自身参考最优%` > `alignedRatio`。
- 「自身参考最优%」低不一定是分样错：同一位点的参考（5-x、27-x、28-x、33-x、29-x、38-x）
  互为近似序列，重复次数不同会让 read 更偏爱「兄弟」参考。5-12（52.5-55.8%）、27-1（2.9-3.9%）、
  27-9（0-4.5%）等都属于这一类，且两个 run 表现一致。
- 真正需要区分的是**参考长度 < 目标 read 长度**的组：这些一代参考缺侧翼/不完整，
  identity 会被系统性低估。第二批里命中该条件的是 27-1（384bp）、27-9（377bp）
  与 Run0002 的全部 33 系列（782-975bp vs 目标 read p50 1021-1051bp）。

## 两 run 同样本一致性（Run0007 72服务器 vs Run0001 153服务器-PRO，同 15 组质粒）

- identity 逐组差异 ≤0.20 个百分点（最大 27-4）；identity 均值 99.556% vs 99.605%
- mmRate 不同量级：Run0007 中位 0.058%（最高 5-10 为 0.198%），Run0001 中位 0.019%（最高 0.032%）

## 异常组（需复核）

| 样本 | run | 现象 | 判断 |
|---|---|---|---|
| STR29-10 | Run0002 | 比对率 55.9%，queryCoverage 0.025，identity 92.88%（全批最低），mmRate 1.66%，identity≥0.83 只有 86.8%；≥800bp read 中 **0%** 以 STR29-10 为最佳命中，81.4% 最佳命中 STR33-6、13.1% STR33-5、3.4% STR33-3 | 该孔里没有 29-10 的插入序列，实际是 33 位点的物质（污染或分样/孔位命名错位）。与 8/18 `ccs比对统计.tsv` 的「目标分子基本缺失」（目标 read p50=48bp、目标比对率 1.3%）互相印证。建议复核 24标签-7 / barcode-07 的孔位对应关系并重做该孔 |
| STR33-3 | Run0002 | identity 97.80%（除 29-10 外全批最低），identity≥0.99% = 0（没有任何 read 达到 99%），mmRate 0.94% 与 longIndelRatio 1.3% 同时偏高，比对率 58.1%；read 最佳命中 82.6% 给 STR33-6，自身仅 3.8% | 差异不是零散测序错误：参考（782bp）短于目标 read（p50 1044bp），一代参考不完整／等位不合。需重新确认 33-3 的一代 merged 序列 |
| STR33-5 / 33-6 / 33-11 | Run0002 | identity 98.38-98.98%（33-11/33-6/33-5），mmRate 0.29-0.67%，longIndelRatio 2.5-5.4%，比对率 50.8-59.7%，identity≥0.99% 49-80%；三孔 read 之间互为最佳命中（33-11 孔 79.3% 命中 STR33-5，33-6 孔 77.4% 命中 STR33-5） | 同 33 位点不同重复次数的参考彼此无法区分，且目标 read 长于参考。属参考问题，不是测序问题；若要用 CCS 定重复数，应以 CCS 自组装序列为参考 |
| STR5-10 / STR25-1 | Run0007 | mmRate 0.198% / 0.106%，明显高于 Run0001 同一样本（0.016% / 0.021%）；STR25-1 是本 run identity 最低（99.343%），STR5-10 倒数第三（99.414%） | 不是 read 质量标签造成的（≥Q20 两 run 基本相同：87.7% vs 89.2%），是比对后错配率确实更高；量级仍小（identity ≥99.3%），可继续使用，统计时不要与 Run0001 混算 |

其余 37 组：identity 99.39-99.88%，mmRate ≤0.090%，比对率 ≥94.5%。
"""


def ref_length(ref_fa):
    """Reference length in bp: .fai if present, else sum the sequence lines."""
    fai = ref_fa + ".fai"
    if os.path.exists(fai):
        parts = open(fai).read().split("\t")
        if len(parts) > 1 and parts[1].strip().isdigit():
            return parts[1].strip()
    if not os.path.exists(ref_fa):
        return ""
    total = 0
    for line in open(ref_fa):
        if not line.startswith(">"):
            total += len(line.strip())
    return str(total)


def parse_selfpref(path):
    """[(plasmid, own-ref win %)] in the order cross_map_diagnostic.py printed them.

    The log carries no run label (plasmid names repeat across runs), so blocks
    are matched positionally against the combined table's row order; main()
    verifies the names line up before trusting that.
    """
    out = []
    for line in open(path, errors="ignore"):
        m = re.match(r"^(\S+)\s+reads>=.*own-ref wins:\s+\d+ \(\s*([\d.]+)%\)", line)
        if m:
            out.append((m.group(1), float(m.group(2))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="dir holding <run>/summary.tsv")
    ap.add_argument("--mapping", required=True, help="plasmid<TAB>barcode<TAB>RUN号 TSV")
    ap.add_argument("--out", required=True, help="combined TSV to write")
    ap.add_argument("--report", help="also write a per-run markdown report here")
    ap.add_argument("--diag", help="cross_map_diagnostic.py output, for the 自身参考最优 column")
    args = ap.parse_args()

    # barcode-NN -> Adaptor-barcodeNNN-M
    b2o = {}
    for row in csv.reader(open(BARCODE2NAME_TSV), delimiter="\t"):
        if len(row) >= 2 and row[0] != "barcodename":
            b2o[row[1].strip()] = row[0].strip()

    # (run, plasmid) -> ("24标签-N", N)
    lab = {}
    for r in csv.DictReader(open(args.mapping), delimiter="\t"):
        pl = (r.get("plasmid") or "").strip()
        bc = (r.get("barcode") or "").strip()
        run = (r.get("RUN号") or "").strip()
        if not pl or not bc:
            continue
        lab[(run, pl)] = (bc, int(re.search(r"(\d+)\s*$", bc).group(1)))

    rows, metric_cols, missing = [], None, []
    for sfile in sorted(glob.glob(os.path.join(args.root, "*/summary.tsv"))):
        rund = os.path.dirname(sfile)
        run = os.path.basename(rund)
        rdr = csv.DictReader(open(sfile), delimiter="\t")
        if metric_cols is None:
            metric_cols = rdr.fieldnames[1:]
        elif rdr.fieldnames[1:] != metric_cols:
            print(f"WARN {run}: metric columns differ; extra columns dropped",
                  file=sys.stderr)
        for row in rdr:
            m = re.match(r"Barcode(\d+)_(.+)$", row["sample"])
            n, pl = int(m.group(1)), m.group(2)
            name = f"barcode-{n:02d}"
            ref_fa = os.path.join(rund, "refs", f"{pl}.fa")
            if (run, pl) not in lab:
                missing.append((run, pl))
            rows.append({
                "run": run,
                "plasmid": f"STR{pl}",
                "24标签": lab.get((run, pl), ("", n))[0],
                "barcodename": name,
                "barcode_orig": b2o.get(name, ""),
                "fastq": os.path.join(rund, "barcode_view", f"Barcode{n:02d}.fastq"),
                "ref": ref_fa,
                "ref_len_bp": ref_length(ref_fa),
                **{k: row.get(k, "") for k in metric_cols},
            })

    with open(args.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=META + (metric_cols or []),
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        ordered = sorted(rows, key=lambda r: (r["run"], r["barcodename"]))
        for r in ordered:
            w.writerow(r)
    print(f"{len(rows)} rows -> {args.out}")
    if missing:
        print("WARN rows with no mapping entry: " +
              ", ".join(f"{r}/{p}" for r, p in missing), file=sys.stderr)

    if args.report:
        pref = parse_selfpref(args.diag) if args.diag and os.path.exists(args.diag) else []
        if pref:
            if len(pref) != len(ordered):
                sys.exit(f"{args.diag}: {len(pref)} blocks but {len(ordered)} rows — "
                         "positional join is unsafe, regenerate the diagnostic")
            for (pl, _), r in zip(pref, ordered):
                if pl != r["plasmid"].replace("STR", "", 1):
                    sys.exit(f"{args.diag}: block {pl} != row {r['plasmid']} — "
                             "diagnostic order no longer matches the combined table")
            selfpref = {(r["run"], r["plasmid"]): p
                        for (_, p), r in zip(pref, ordered)}
        else:
            if args.diag:
                print(f"WARN {args.diag}: no blocks parsed", file=sys.stderr)
            selfpref = {}
        write_markdown(ordered, selfpref, args.report)
        print(f"report -> {args.report}")


if __name__ == "__main__":
    main()
