#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pick_representative.py —— 为每个 (位点×分量) 挑选"代表性样本"并打印展示

目的：最终/收敛报告按"位点粒度"统计，多个样品测同一个位点，看聚合统计不直观。
      本脚本在每个分量内挑出能"把这个位点讲清楚"的原型样本（archetype），
      让报告能配 1~3 个具体样品 + 它们的 CN 分布，把该位点的行为落地。

分类依据（与《13-收敛性专项报告》10.1/10.2/10.3 完全一致，直接沿用其结论）：
  clean    = 10.1 收敛良好分量          —— 每个分量挑 1 个干净主峰样品即可
  hetz     = 10.2 杂合双态分量 (3 个)   —— 挑 1 纯合样品 + 1 杂合双峰样品
  stutter  = 10.3 收敛差分量 (11 个)    —— 挑 1 相对好的 + 1 最差发散样品

为什么不再自己发明阈值：报告 13 已经按分量把 43 个分量分好类（10.1/10.2/10.3），
直接沿用可保证"展示"与"结论"口径一致，且避免把长重复位点的"鬼峰/嵌合"误判成杂合。

深度只作平手打破（同档里挑 depth 足够、图看得清楚的），因为报告 13 已证明深度不影响收敛。

所有指标直接读 data/cn/cn_summary.tsv，输出可回溯。
用法：
    python3 pick_representative.py               # 打印全部 43 分量
    python3 pick_representative.py --only PM12 PM25
    python3 pick_representative.py --out data/cn/representative_samples.tsv
"""
import argparse
import csv
import os
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DEF_SUMMARY = os.path.normpath(os.path.join(HERE, "..", "data", "cn", "cn_summary.tsv"))

# ---------------------------------------------------------------------------
# 分类：键 = (locus, unit)。值 = "clean" | "hetz" | "stutter"
# 来源：13 号报告 第 10 章
#   10.1 收敛良好（29 个分量）→ clean
#   10.2 杂合双态（3 个分量） → hetz
#   10.3 收敛差（11 个分量）  → stutter
# ---------------------------------------------------------------------------
CLASS = {}
_CLEAN = [
    "PM06:GTAT", "PM07:AT", "PM07:TGTT", "PM10:CATC", "PM12:CT",
    "PM13:AAGA", "PM14:AG", "PM14:TCGC", "PM16:GGAGA", "PM17:AGTGA",
    "PM18:AGC", "PM20:TTTG", "PM21:CA", "PM23:AGGT", "PM24:GGTGT",
    "PM25:TA", "PM26:GCC", "PM26:GCG", "PM27:TA", "PM28:CTG",
    "PM29:TGCCA", "PM31:AG", "PM33:AG", "PM34:TG", "PM36:GGC",
    "PM36:TAT", "PM37:CTAG", "PM38:GA", "PM39:TGG",
]
_HETZ = ["PM25:CT", "PM32:CTG", "PM19:TC"]
_STUTTER = [
    "PM01:TC", "PM04:CT", "PM05:AG", "PM08:AG", "PM09:TC", "PM10:GA",
    "PM15:GA", "PM22:CT", "PM25:AC", "PM30:GT", "PM30:GCA",
]
for k in _CLEAN:
    loc, unit = k.split(":"); CLASS[(loc, unit)] = "clean"
for k in _HETZ:
    loc, unit = k.split(":"); CLASS[(loc, unit)] = "hetz"
for k in _STUTTER:
    loc, unit = k.split(":"); CLASS[(loc, unit)] = "stutter"

# 报告 13 未明确归类的分量（10-17-40 深度不足 + PM23/PM24/PM31/PM33/PM34/PM38/PM39 等）
# 用组件级统计兜底判类（中位 P≥60 且无强杂合 → clean，否则 stutter）。
def _fallback_class(rs):
    P = [f(r["P_pct"]) for r in rs]; P = [p for p in P if p is not None]
    if not P:
        return "stutter"
    medP = sorted(P)[len(P)//2]
    return "clean" if medP >= 60 else "stutter"


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# archetype 挑选
# ---------------------------------------------------------------------------
def _row(r):
    return {
        "sample": r["sample"], "P": f(r["P_pct"]), "s1": f(r["s1_pct"]),
        "s2": f(r["s2_pct"]), "mode": r["mode_cn"], "second": r["second_cn"],
        "second_pct": f(r["second_pct"]) or 0.0, "depth": int(f(r["depth"]) or 0),
        "cn_dist": r["cn_dist"],
        "shape": "two-peak" if (f(r["second_pct"]) or 0) >= 30 else "spike",
    }


def pick_clean(rs):
    """干净分量：众数主 CN 里 P 最高的一个。"""
    c = Counter(r["mode_cn"] for r in rs)
    top = c.most_common(1)[0][0]
    cands = [r for r in rs if r["mode_cn"] == top]
    best = max(cands, key=lambda r: (f(r["P_pct"]) or 0, f(r["depth"]) or 0))
    return [{**_row(best), "archetype": "clean"}]


def pick_hetz(rs):
    """杂合分量：1 个纯合（单峰干净）+ 1 个杂合双峰（second_pct 最大）。

    注意：对**已分类为 hetz 的分量**（报告 13 10.2 已认定杂合双态），
    直接取 second_pct 最大样品作为"杂合代表"，不再用间隔整除约束——
    因为 2bp 基序的真等位基因本就可能相差 1 个单元（如 PM19 TC6/TC7），
    这是**真实双等位**而非 stutter。
    """
    c = Counter(r["mode_cn"] for r in rs)
    # 主等位 = 分量众数（在 hetz 位点通常是较大 CN，如 CT9 / CTG6 / TC7）
    # 但从纯合角度，取"该等位里 P 最高"的样品
    homo = []
    homo_cands = [r for r in rs if (f(r["second_pct"]) or 0) < 30 and (f(r["P_pct"]) or 0) >= 60]
    if homo_cands:
        hbest = max(homo_cands, key=lambda r: (f(r["P_pct"]) or 0, f(r["depth"]) or 0))
        homo = [{**_row(hbest), "archetype": "homo"}]
    # 杂合代表：second_pct 最大
    hetz_cands = [r for r in rs if (f(r["second_pct"]) or 0) >= 30]
    hetz = []
    if hetz_cands:
        best = max(hetz_cands, key=lambda r: (f(r["second_pct"]) or 0, f(r["depth"]) or 0))
        hetz = [{**_row(best), "archetype": "hetero"}]
    return homo + hetz


def pick_stutter(rs):
    """stutter 分量：1 个相对好的（P 最高的低收敛）+ 1 个最差发散的（P 最低、depth 足够）。"""
    good = max(rs, key=lambda r: (f(r["P_pct"]) or 0, f(r["depth"]) or 0))
    # 最烂：P 最低，但 depth 至少 >=10（图看得见）
    bad_cands = [r for r in rs if (f(r["depth"]) or 0) >= 10]
    if not bad_cands:
        bad_cands = rs
    bad = min(bad_cands, key=lambda r: (f(r["P_pct"]) or 99,))
    g, b = _row(good), _row(bad)
    g["archetype"] = "best-of-stutter"
    b["archetype"] = "worst-stutter"
    if g["sample"] == b["sample"]:
        return [b]
    return [g, b]


def pick(rs, kind):
    if kind == "clean":
        return pick_clean(rs)
    if kind == "hetz":
        return pick_hetz(rs)
    return pick_stutter(rs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=DEF_SUMMARY)
    ap.add_argument("--only", nargs="*", default=None, help="只看指定位点")
    ap.add_argument("--out", default=None, help="写 TSV")
    args = ap.parse_args()

    with open(args.summary) as fh:
        rs = list(csv.DictReader(fh, delimiter="\t"))
    if args.only:
        rs = [r for r in rs if r["locus"] in args.only]

    comp = defaultdict(list)
    for r in rs:
        comp[(r["locus"], r["unit"])].append(r)

    out_rows, lines = [], []
    for (locus, unit) in sorted(comp):
        kind = CLASS.get((locus, unit), _fallback_class(comp[(locus, unit)]))
        picks = pick(comp[(locus, unit)], kind)
        tag = {"clean": "干净收敛", "hetz": "杂合双态", "stutter": "stutter 发散"}[kind]
        lines.append(f"\n### {locus} / {unit}  ——  {tag}  共 {len(comp[(locus,unit)])} 样品")
        for e in picks:
            if e.get("sample") in (None, "-"):
                lines.append(f"  {e['archetype']:<20} （{e['archetype']}）")
                continue
            lines.append(
                f"  {e['archetype']:<20} {e['sample']:<14} P={e['P']:>6.1f}%  "
                f"s1={e['s1']:>5.1f}% s2={e['s2']:>6.1f}%  mode={e['mode']}  "
                f"2nd={e['second']}({e['second_pct']:.0f}%)  depth={e['depth']}"
            )
            lines.append(f"     cn_dist: {e['cn_dist']}")
            out_rows.append([locus, unit, kind, e["archetype"], e["sample"],
                             e["P"], e["s1"], e["s2"], e["mode"], e["second"],
                             e["second_pct"], e["depth"], e["cn_dist"], e["shape"]])

    print("\n".join(lines))
    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["locus", "unit", "class", "archetype", "sample",
                        "P_pct", "s1_pct", "s2_pct", "mode_cn", "second_cn",
                        "second_pct", "depth", "cn_dist", "shape"])
            w.writerows(out_rows)
        print(f"\n[saved] {args.out}  ({len(out_rows)} 个代表样本)")


if __name__ == "__main__":
    main()
