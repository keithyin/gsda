"""SUSPICIOUS_SITES 新旧模型对比。

给定同一源数据的两组 consensus suspicious_sites 输出（旧模型 vs 新模型），
量化新模型相对旧模型在 **降噪**（CLEANED）与 **敏感性**（NEW_SUS）上的差异。

需求文档见同目录 suspicious_sites_model_comparison.md。

设计要点
--------
- 输入是两个目录，各含一批 `<sample>.consensus.suspicious_sites.{txt,csv}`。
  `.txt` 为 tab 分隔、`.csv` 为逗号分隔，是同一份数据的两种格式 → **每组只用一种**
  （--format 指定），不叠加。
- 按表头名解析字段（真实文件为 20 列，含 Name / #Clust / ... / Mut_context / insert_seq，
  文档写的是 18 列；按名解析可兼容两种）。
- 以 (样本, Position) 为键做差异分类：
    CLEANED  : OLD 有、NEW 无 —— 降噪
    NEW_SUS  : NEW 有、OLD 无 —— 新爆
    SHIFTED  : 同一样本内 indel 坐标整体重定位（见 detect_shifted），不计净增减
- 对每个差异位点给归因标签（低比例 SNP / 超低深度 indel / 重复同聚区 indel / 其它）。

用法示例
--------
    conda run -n py38 python suspicious_sites_compare.py \
        --new /data1/.../called_barcode_v4_new_model/Consensus/Consensus/Suspicious_sites \
        --old /data1/.../called_barcode_v4/Consensus/Consensus/Suspicious_sites \
        --format txt \
        --outdir /tmp/susp_cmp
"""

import os
import sys
import glob
import argparse
from collections import OrderedDict


SUSP_SUFFIX = ".consensus.suspicious_sites"

# 解析所需字段（按表头名匹配，缺列时为 None）
NEED_FIELDS = [
    "Position", "Consensus", "Depth", "Top1_base", "Freq_of_top1_base",
    "Top2_base", "Freq_of_top2_base", "Homo_num", "Sus_types", "insert_seq",
]


def _fnum(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _fint(v, default=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _sample_name(path, fmt):
    base = os.path.basename(path)
    # 去掉扩展名与后缀 → 样本名
    ext = ".txt" if fmt == "txt" else ".csv"
    if base.endswith(ext):
        base = base[: -len(ext)]
    if base.endswith(SUSP_SUFFIX):
        base = base[: -len(SUSP_SUFFIX)]
    return base


def load_group(dirpath, fmt):
    """读一个目录，返回 OrderedDict: sample -> list[site dict]（按 Position 排序）。

    site dict 含 Position(int) 及 NEED_FIELDS 字段。空行 / 非数据行忽略。
    """
    sep = "\t" if fmt == "txt" else ","
    pat = os.path.join(dirpath, "*" + SUSP_SUFFIX + ("." + fmt))
    files = sorted(glob.glob(pat))
    out = OrderedDict()
    for f in files:
        sample = _sample_name(f, fmt)
        sites = []
        with open(f) as fh:
            header = fh.readline().rstrip("\n").split(sep)
            idx = {name: i for i, name in enumerate(header)}
            if "Position" not in idx:
                continue
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                parts = line.split(sep)
                if len(parts) < len(header):
                    # 短行：不足补齐，避免越界
                    parts = parts + [""] * (len(header) - len(parts))

                def g(col):
                    return parts[idx[col]] if col in idx else None

                pos = _fint(g("Position"), None)
                if pos is None:
                    continue
                sites.append({
                    "Position": pos,
                    "Consensus": g("Consensus"),
                    "Depth": _fint(g("Depth")),
                    "Top1_base": g("Top1_base"),
                    "Freq_of_top1_base": _fnum(g("Freq_of_top1_base")),
                    "Top2_base": g("Top2_base"),
                    "Freq_of_top2_base": _fnum(g("Freq_of_top2_base")),
                    "Homo_num": _fint(g("Homo_num")),
                    "Sus_types": g("Sus_types"),
                    "insert_seq": g("insert_seq"),
                })
        sites.sort(key=lambda s: s["Position"])
        out[sample] = sites
    return out


def attr_label(site, low_vaf_thresh, ultra_low_depth, homo_thresh):
    """D4 归因标签。"""
    t = site["Sus_types"]
    top2 = (site["Top2_base"] or "")
    freq2 = site["Freq_of_top2_base"]
    depth = site["Depth"]
    homo = site["Homo_num"]
    if t == "SNP":
        if top2 in ("", "-") or freq2 >= low_vaf_thresh:
            return "SNP"
        return "low-VAF-SNP"
    if t == "Indel":
        if depth <= ultra_low_depth:
            return "ultra-low-depth-indel"
        if homo >= homo_thresh:
            return "repeat/homo-indel"
        return "indel"
    return t or "other"


def fmt_site(s):
    """逐位点一行画像（D4）。"""
    ins = (s["insert_seq"] or "").strip()
    if ins in ("", "NA"):
        ins = ""
    else:
        ins = ins[:24]
    t1 = "{}({:.3f})".format(s["Top1_base"], s["Freq_of_top1_base"])
    t2 = "{}({:.3f})".format(s["Top2_base"], s["Freq_of_top2_base"])
    bits = "pos={} Cons={} Depth={} Top1={} Top2={} Homo={} Types={}".format(
        s["Position"], s["Consensus"], s["Depth"], t1, t2, s["Homo_num"], s["Sus_types"])
    if ins:
        bits += "  ins={}".format(ins)
    return bits


def detect_shifted(cleaned, newsus, pos_gap, depth_ratio):
    """把同一样本内 CLEANED / NEW_SUS 的位点聚类；若一个簇内两类并存、
    全为 Indel、深度量级一致，则整个簇判定为 SHIFTED（重定位）。

    返回 (shifted_cleaned_positions:set, shifted_newsus_positions:set, shifts:list[str])
    shifts 元素为 "OLD区间 a–b (n=..) → NEW区间 c–d (n=..)" 描述。
    """
    # 所有差异位点（带来源），按位置排
    pts = []
    for p, s in cleaned.items():
        pts.append((p, "old", s))
    for p, s in newsus.items():
        pts.append((p, "new", s))
    if not pts:
        return set(), set(), []
    pts.sort(key=lambda x: x[0])

    # 贪心按位置间距聚簇
    clusters = []
    cur = [pts[0]]
    for item in pts[1:]:
        if item[0] - cur[-1][0] <= pos_gap:
            cur.append(item)
        else:
            clusters.append(cur)
            cur = [item]
    clusters.append(cur)

    shifted_old, shifted_new = set(), set()
    shifts = []
    for cl in clusters:
        olds = [it for it in cl if it[1] == "old"]
        news = [it for it in cl if it[1] == "new"]
        if not olds or not news:
            continue
        # 必须全是 indel（snp 不参与重定位判定）
        if not all(it[2]["Sus_types"] == "Indel" for it in cl):
            continue
        # 深度量级一致：最大/最小 Depth 比值
        depths = [it[2]["Depth"] for it in cl if it[2]["Depth"] > 0]
        if len(depths) >= 2:
            if max(depths) / min(depths) > depth_ratio:
                continue
        shifted_old.update(p for p, _, _ in olds)
        shifted_new.update(p for p, _, _ in news)
        lo_o, hi_o = min(p for p, _, _ in olds), max(p for p, _, _ in olds)
        lo_n, hi_n = min(p for p, _, _ in news), max(p for p, _, _ in news)
        shifts.append("OLD区间 {}–{} (n={}) → NEW区间 {}–{} (n={})".format(
            lo_o, hi_o, len(olds), lo_n, hi_n, len(news)))
    return shifted_old, shifted_new, shifts


def compare(new_dir, old_dir, fmt, outdir, pos_gap, depth_ratio,
            low_vaf_thresh, ultra_low_depth, homo_thresh):
    new = load_group(new_dir, fmt)
    old = load_group(old_dir, fmt)

    all_samples = sorted(set(new) | set(old))
    new_only_samples = sorted(set(new) - set(old))
    old_only_samples = sorted(set(old) - set(new))

    d1_new = sum(len(v) for v in new.values())
    d1_old = sum(len(v) for v in old.values())
    d2_new = sum(1 for v in new.values() if v)
    d2_old = sum(1 for v in old.values() if v)

    # 逐位点分类 + 归因标签分布
    label_dist = {}
    sample_blocks = []
    totals = {"cleaned": 0, "new_sus": 0, "shifted": 0}

    for sample in all_samples:
        n_sites = old.get(sample) or []
        m_sites = new.get(sample) or []
        old_by_pos = OrderedDict((s["Position"], s) for s in n_sites)
        new_by_pos = OrderedDict((s["Position"], s) for s in m_sites)

        cleaned = OrderedDict((p, s) for p, s in old_by_pos.items()
                              if p not in new_by_pos)
        newsus = OrderedDict((p, s) for p, s in new_by_pos.items()
                             if p not in old_by_pos)

        sh_old, sh_new, shifts = detect_shifted(cleaned, newsus, pos_gap, depth_ratio)

        # 剔除 shifted，得到净增减
        net_cleaned = OrderedDict((p, s) for p, s in cleaned.items() if p not in sh_old)
        net_newsus = OrderedDict((p, s) for p, s in newsus.items() if p not in sh_new)

        # 空样本（两组都为 0）省略明细
        if not cleaned and not newsus:
            continue

        lines = []
        lines.append("### {}   OLD={} NEW={}  cleaned={} newly_susp={} shifted={}".format(
            sample, len(n_sites), len(m_sites),
            len(net_cleaned), len(net_newsus), len(sh_old)))
        for p, s in net_cleaned.items():
            lb = attr_label(s, low_vaf_thresh, ultra_low_depth, homo_thresh)
            label_dist[lb] = label_dist.get(lb, 0) + 1
            lines.append("  [CLEANED]    {}   #{}".format(fmt_site(s), lb))
        for p, s in net_newsus.items():
            lb = attr_label(s, low_vaf_thresh, ultra_low_depth, homo_thresh)
            label_dist[lb] = label_dist.get(lb, 0) + 1
            lines.append("  [NEW_SUS]    {}   #{}".format(fmt_site(s), lb))
        for sh in shifts:
            lines.append("  [SHIFTED]    {}".format(sh))
        sample_blocks.append("\n".join(lines))

        totals["cleaned"] += len(net_cleaned)
        totals["new_sus"] += len(net_newsus)
        totals["shifted"] += len(sh_old)

    # ---- 汇总 ----
    L = []
    L.append("=" * 72)
    L.append("SUSPICIOUS_SITES 新旧模型对比   (format={})".format(fmt))
    L.append("  NEW: {}".format(new_dir))
    L.append("  OLD: {}".format(old_dir))
    L.append("=" * 72)
    L.append("")
    L.append("[D1] 位点总量（数据行数之和）")
    L.append("  OLD={:<8d} NEW={:<8d} diff(NEW-OLD)={:+d}".format(d1_old, d1_new, d1_new - d1_old))
    L.append("")
    L.append("[D2] 命中样本数（位点数>0 的样本）")
    L.append("  OLD={:<8d} NEW={:<8d} diff={:+d}".format(d2_old, d2_new, d2_new - d2_old))
    L.append("  样本总数: {}   NEW-only: {}   OLD-only: {}".format(
        len(all_samples), len(new_only_samples), len(old_only_samples)))
    if new_only_samples:
        L.append("    NEW-only: " + ", ".join(new_only_samples))
    if old_only_samples:
        L.append("    OLD-only: " + ", ".join(old_only_samples))
    L.append("")
    L.append("[D3] 逐位点差异（剔除 SHIFTED 后的净增减）")
    L.append("  CLEANED (OLD有NEW无, 降噪) = {}".format(totals["cleaned"]))
    L.append("  NEW_SUS (NEW有OLD无, 新爆) = {}".format(totals["new_sus"]))
    L.append("  SHIFTED (重定位, 不计净增减) = {}".format(totals["shifted"]))
    L.append("")
    L.append("[D4] 差异位点归因标签分布（CLEANED + NEW_SUS）")
    for lb in sorted(label_dist, key=lambda k: -label_dist[k]):
        L.append("  {:<24} {}".format(lb, label_dist[lb]))
    L.append("")
    L.append("-" * 72)
    L.append("逐样本明细")
    L.append("-" * 72)
    if sample_blocks:
        L.append("\n\n".join(sample_blocks))
    else:
        L.append("（两组无差异位点）")

    report = "\n".join(L)
    os.makedirs(outdir, exist_ok=True)
    rpath = os.path.join(outdir, "suspicious_sites_comparison.txt")
    with open(rpath, "w") as fh:
        fh.write(report + "\n")
    print(report)
    print("\n[输出] {}".format(rpath))
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="SUSPICIOUS_SITES 新旧模型对比（降噪 vs 敏感性）")
    ap.add_argument("--new", required=True, help="新模型 Suspicious_sites 目录")
    ap.add_argument("--old", required=True, help="旧模型 Suspicious_sites 目录")
    ap.add_argument("--format", choices=["txt", "csv"], default="txt",
                    help="用哪种分隔格式解析（默认 txt）")
    ap.add_argument("--outdir", required=True, help="输出目录")
    ap.add_argument("--pos-gap", type=int, default=20,
                    help="SHIFTED 聚类时相邻位点最大间距（默认 20）")
    ap.add_argument("--depth-ratio", type=float, default=3.0,
                    help="SHIFTED 判定最大/最小深度比阈值（默认 3.0）")
    ap.add_argument("--low-vaf-thresh", type=float, default=0.20,
                    help="Top2 频率低于此值判为低比例 SNP（默认 0.20）")
    ap.add_argument("--ultra-low-depth", type=int, default=5,
                    help="Indel 深度 <= 此值判为超低深度（默认 5）")
    ap.add_argument("--homo-thresh", type=int, default=2,
                    help="Homo_num >= 此值判为重复/同聚区 indel（默认 2）")
    args = ap.parse_args()

    for d, tag in ((args.new, "new"), (args.old, "old")):
        if not os.path.isdir(d):
            print("ERROR: {} 目录不存在: {}".format(tag, d))
            return 1

    sys.exit(compare(args.new, args.old, args.format, args.outdir,
                     args.pos_gap, args.depth_ratio, args.low_vaf_thresh,
                     args.ultra_low_depth, args.homo_thresh))


if __name__ == "__main__":
    main()
