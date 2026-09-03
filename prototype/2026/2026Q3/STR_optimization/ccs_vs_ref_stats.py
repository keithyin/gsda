#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基于「样本-run-barcode对应关系.tsv」，把每组 barcode fastq 与其参考序列（一代测序 merged 结果）
做 minimap2 比对，按 read 物种分开统计准确率 / 覆盖率。

关键背景（经序列检查确认）：
  - 每组数据里混合两个物种：
    短 read（<=2kb，~600-1100bp）：完整比对参考序列 → 目标分子（STR 扩增子）
    长 read（>2kb，2.4-5.9kb）：环状质粒全长分子（首尾相同=闭环），
      仅与参考 3' 端公共载体 MCS 区（~48bp，XbaI/XhoI/LacZ）匹配 → 背景物种
  - 因此「准确率/覆盖率」基于短 read 计算；长 read 只报告占比（纯度指标）。

指标（每组 = 样本 × run）：
  目标(短read)：reads 数 / 占比 / 比对率 / 准确率(1-NM/(M+D)) / QV / indel率
                平均/中位/最低覆盖，>=1x/5x/10x 参考碱基占比
  背景(长read)：占比、长度 p50
  异常标记：目标 read 中位长度 > 参考长度*1.1（STR 扩增）；目标比对率 < 50%（目标缺失）

输出：ccs比对统计.tsv / ccs比对统计报告.md / ref_per_sample/（单样本参考）
"""
import csv
import math
import re
import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pysam

HERE = Path(__file__).resolve().parent
MAPPING_TSV = HERE / "样本-run-barcode对应关系.tsv"
REF_DIR = HERE / "ref_per_sample"
BAM_DIR = HERE / "ccs_vs_ref_bams"
N_WORKERS = 16
N_THREADS = 8
LONG_CUTOFF = 2000  # bp，短/长 read 分界（实测双峰：~600-1100 vs ~2400-6000）


def sample_token(sample):
    return re.sub(r"^STR", "", sample)


def build_single_ref(sample, src_path):
    """从（可能多记录的）参考文件里抽出该样本的单条记录，写到 REF_DIR。"""
    REF_DIR.mkdir(exist_ok=True)
    out = REF_DIR / f"{sample}.fa"
    if out.exists():
        return out
    tok = sample_token(sample)
    hdr_pat = re.compile(r"^.*STR" + re.escape(tok) + r"(\s|$)")
    lines_out = []
    with open(src_path, encoding="utf-8") as f:
        in_rec, started = False, False
        for line in f:
            if line.startswith(">"):
                if started:
                    break
                in_rec = bool(hdr_pat.match(line[1:].strip()))
                if in_rec:
                    started = True
                    lines_out.append(">" + sample + "\n")
            elif in_rec:
                lines_out.append(line)
    if not started:
        raise ValueError(f"参考文件 {src_path} 中找不到样本 {sample}")
    out.write_text("".join(lines_out), encoding="utf-8")
    return out


def run_one(row):
    """row: 映射表一行。返回统计 dict。"""
    batch, sample, run = row["批次"], row["样本名"], row["RUN号"]
    fastq = Path(row["fastq文件"])
    ref_src = Path(row["参考序列"])
    ref = build_single_ref(sample, ref_src)

    BAM_DIR.mkdir(exist_ok=True)
    bam = BAM_DIR / f"{sample}__{run}.bam"
    if not bam.exists():
        cmd = (f"minimap2 -ax sr -t {N_THREADS} {ref} {fastq} 2>/dev/null"
               f" | samtools view -b - 2>/dev/null | samtools sort -O bam -o {bam} - 2>/dev/null")
        p = subprocess.run(cmd, shell=True, executable="/bin/bash")
        if p.returncode != 0 or not bam.exists():
            return {"批次": batch, "样本": sample, "RUN号": run, "错误": "比对失败"}

    aln = pysam.AlignmentFile(str(bam))
    ref_len = aln.get_reference_length(aln.references[0])
    total = 0
    long_lens = []
    s_total = s_mapped = s_nm = s_mref = s_ins = s_dels = 0
    s_qlens = []
    starts, ends = [], []
    for r in aln:
        total += 1
        if r.query_length > LONG_CUTOFF:  # 长 read：背景质粒物种，只计占比
            long_lens.append(r.query_length)
            continue
        s_total += 1
        if r.flag & 0x4:
            continue
        ct = r.cigartuples
        m = sum(c[1] for c in ct if c[0] == 0)
        d = sum(c[1] for c in ct if c[0] == 2)
        i = sum(c[1] for c in ct if c[0] == 1)
        if m + d == 0:
            continue
        s_mapped += 1
        s_mref += m + d
        s_ins += i
        s_dels += d
        s_nm += r.get_tag("NM", 0)
        s_qlens.append(r.query_length)
        s = min(r.reference_start, ref_len)
        e = min(s + m + d, ref_len)
        if e > s:
            starts.append(s)
            ends.append(e)
    aln.close()

    res = {
        "批次": batch, "样本": sample, "RUN号": run, "参考长度bp": ref_len,
        "总reads": total,
        "目标reads": s_total,
        "目标占比%": round(s_total / total * 100, 1) if total else 0,
        "背景reads": len(long_lens),
        "背景占比%": round(len(long_lens) / total * 100, 1) if total else 0,
    }
    if long_lens:
        res["背景长度p50"] = int(statistics.median(long_lens))
    if s_total:
        res["目标比对上reads"] = s_mapped
        res["目标比对率%"] = round(s_mapped / s_total * 100, 1) if s_total else 0
        res["目标read长度p50"] = int(statistics.median(s_qlens)) if s_qlens else 0
    if s_mref:
        err = s_nm / s_mref
        res["准确率%"] = round((1 - err) * 100, 3)
        res["QV"] = round(-10 * math.log10(err), 2) if err > 0 else ">99"
        res["indel率%"] = round((s_ins + s_dels) / s_mref * 100, 3)
    if starts:
        st = np.array(starts, dtype=np.int64)
        en = np.array(ends, dtype=np.int64)
        delta = np.bincount(st, minlength=ref_len + 2) - np.bincount(en, minlength=ref_len + 2)
        cov = np.cumsum(delta)[:ref_len]
        res["平均覆盖"] = round(float(cov.mean()), 0)
        res["中位覆盖"] = round(float(np.median(cov)), 0)
        res["最低覆盖"] = int(cov.min())
        res["≥1x覆盖%"] = round(float((cov >= 1).mean() * 100), 1)
        res["≥5x覆盖%"] = round(float((cov >= 5).mean() * 100), 1)
        res["≥10x覆盖%"] = round(float((cov >= 10).mean() * 100), 1)
    # 异常标记
    flags = []
    if res.get("目标比对率%") is not None and res["目标比对率%"] < 50:
        flags.append("目标分子基本缺失")
    if res.get("目标read长度p50", 0) > ref_len * 1.1:
        flags.append("目标reads显著长于参考(STR扩增?)")
    if flags:
        res["标记"] = "；".join(flags)
    return res


def main():
    rows = list(csv.DictReader(open(MAPPING_TSV, encoding="utf-8"), delimiter="\t"))
    print(f"共 {len(rows)} 组待比对")
    results = [None] * len(rows)
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        for i, res in enumerate(ex.map(run_one, rows)):
            results[i] = res
            print(f"[{i+1}/{len(rows)}] {res['批次']} {res['样本']} "
                  f"map%={res.get('目标比对率%')} acc={res.get('准确率%')} "
                  f"cov={res.get('平均覆盖')} 背景%={res.get('背景占比%')} {res.get('标记','')}")

    cols = ["批次", "样本", "RUN号", "参考长度bp", "总reads",
            "目标reads", "目标占比%", "目标read长度p50", "目标比对上reads", "目标比对率%",
            "平均覆盖", "中位覆盖", "最低覆盖", "≥1x覆盖%", "≥5x覆盖%", "≥10x覆盖%",
            "准确率%", "QV", "indel率%",
            "背景reads", "背景占比%", "背景长度p50", "标记", "错误"]
    tsv = HERE / "ccs比对统计.tsv"
    with open(tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in results:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    def agg(vals):
        a = np.array([v for v in vals if isinstance(v, (int, float))])
        if len(a) == 0:
            return "-"
        return f"mean={a.mean():.1f} min={a.min():.1f} max={a.max():.1f}"

    report = [
        "# CCS reads 比对参考序列：准确率 / 覆盖率统计",
        "",
        "- 流程：`minimap2 -ax sr` 比对单样本参考（一代测序 merged 序列），"
        "NM/(M+D) 计碱基错误率；覆盖直方图由比对区间累加",
        "- **read 物种分界**：短 read（≤2kb）= 目标分子（完整比对参考）；"
        "长 read（>2kb）= 环状质粒全长背景分子（仅与参考 3' 端 ~48bp 公共载体 MCS 匹配），"
        "准确率/覆盖率仅按目标 read 计算，背景 read 报告占比",
        "- 明细：`ccs比对统计.tsv`",
        "",
    ]
    batches = OrderedDict()
    for r in results:
        batches.setdefault(r["批次"], []).append(r)
    for bname, rs in batches.items():
        ok = [r for r in rs if "错误" not in r]
        report.append(f"## {bname}（{len(ok)} 组）")
        report.append("")
        report.append("| 样本 | run | 目标reads% | 目标比对率% | 准确率% | QV | 平均覆盖 | ≥10x% | 背景reads% | 标记 |")
        report.append("|---|---|---|---|---|---|---|---|---|---|")
        for r in sorted(ok, key=lambda x: x["样本"]):
            run_short = "/".join(r["RUN号"].split("_")[1:])
            report.append(
                f"| {r['样本']} | {run_short} | {r.get('目标占比%','-')} | {r.get('目标比对率%','-')} | "
                f"{r.get('准确率%','-')} | {r.get('QV','-')} | {r.get('平均覆盖','-')} | "
                f"{r.get('≥10x覆盖%','-')} | {r.get('背景占比%','-')} | {r.get('标记','')} |")
        report.append("")
        report.append(f"- 批次汇总：目标比对率 {agg([r.get('目标比对率%') for r in ok])} | "
                      f"准确率% {agg([r.get('准确率%') for r in ok])} | "
                      f"平均覆盖 {agg([r.get('平均覆盖') for r in ok])} | "
                      f"背景reads占比% {agg([r.get('背景占比%') for r in ok])}")
        accs = sorted((r.get("准确率%"), r["样本"]) for r in ok if isinstance(r.get("准确率%"), (int, float)))
        bg = sorted((r.get("背景占比%"), r["样本"]) for r in ok if isinstance(r.get("背景占比%"), (int, float)))
        if accs:
            report.append(f"- 准确率最低：{accs[0][1]}（{accs[0][0]}%）；最高：{accs[-1][1]}（{accs[-1][0]}%）")
        if bg:
            report.append(f"- 背景reads占比最高：{bg[-1][1]}（{bg[-1][0]}%）；最低：{bg[0][1]}（{bg[0][0]}%）")
        flagged = [r for r in ok if r.get("标记")]
        if flagged:
            report.append("- 异常组：" + "；".join(f"{r['样本']}（{r['标记']}）" for r in flagged))
        report.append("")
    (HERE / "ccs比对统计报告.md").write_text("\n".join(report), encoding="utf-8")
    print("DONE: ccs比对统计.tsv / ccs比对统计报告.md")


if __name__ == "__main__":
    main()
