#!/usr/bin/env python3
"""检查 reference 上指定区域 (substr) 下方是否有 subreads 支撑。

在每条 ASRTC record 的 reference 行（去 gap 后）中精确搜索给定序列字符串，
对每个出现位置，统计其对应 MSA 列下方的 subread 覆盖（非 gap 碱基，
不要求与 ref 碱基一致）。

每个出现位点的统计:
    - per-subread 覆盖率 (cov_pct, 基于完整 MSA 跨度含插入列)、完整跨越数
      (n_spanning, span 全列非 gap)
    - supporting subreads: cov_pct >= --support-cov
    - 逐列 depth: 每个 ref 碱基列上非 gap subread 数
    - supported 判定: n_supporting >= --support-reads

汇总统计:
    - 出现数、supported 比例、平均覆盖率直方图
    - 按区域内 offset 聚合的 depth profile (同一 substr 各出现位置可比)

排除条件:
    record 内所有出现位点的 SMC-REF 一致率均 > --smc-match-thr (默认 90%) 时,
    该 record 不纳入分析 (一致率只比较 ref 区域与对应 SMC 区域, 双 gap 列不计入)。

用法:
    python region_analysis.py input.asrtc.txt --region-seq ACGTACGT
    python region_analysis.py input.asrtc.txt --region-seq ACGTACGT --support-cov 80 --support-reads 3
    python region_analysis.py input.asrtc.txt --region-seq ACGTACGT --output json --out-file result.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from statistics import median


class RegionSupportAnalyzer:
    """检查 msa_seqs 中 reference 指定区域下方的 subread 覆盖。

    msa_seqs 布局: [0]=qual, [1]=SMC consensus, [2]=reference, [3:]=subreads。
    """

    def __init__(
        self,
        filepath: str,
        region_seq: str,
        identity_thr: float = 1.01,
        identity_min_thr: float = 0.0,
        support_cov: float = 50.0,
        support_reads: int = 1,
        smc_match_thr: float = 90.0,
        sample_n: int = 20,
    ):
        self.filepath = filepath
        self.region_seq = region_seq.upper()
        self.identity_thr = identity_thr
        self.identity_min_thr = identity_min_thr
        self.support_cov = support_cov
        self.support_reads = support_reads
        self.smc_match_thr = smc_match_thr
        self.sample_n = sample_n
        self.records: list[dict] = []
        self._load()

    # ---- loading ---------------------------------------------------------

    def _load(self) -> None:
        """加载 JSONL, 仅保留 identity_min_thr <= identity < identity_thr 的记录。"""
        self.records = []
        line_no = 0
        line = ""
        try:
            with open(self.filepath, "r") as fh:
                for line in fh:
                    line = line.strip()
                    line_no += 1
                    if not line:
                        continue
                    rec = json.loads(line)
                    rec_id = rec.get("identity", 1.0)
                    if rec_id < self.identity_min_thr or rec_id >= self.identity_thr:
                        continue
                    self.records.append(rec)

        except Exception:
            raise ValueError(f"line_no:{line_no} line={line}")

    # ---- analysis --------------------------------------------------------

    @staticmethod
    def _find_all(seq: str, sub: str) -> list[int]:
        """查找 sub 在 seq 中的所有出现位置 (步进 1, 含重叠)。"""
        hits = []
        pos = seq.find(sub)
        while pos != -1:
            hits.append(pos)
            pos = seq.find(sub, pos + 1)
        return hits

    @staticmethod
    def _smc_agreement(smc_s: str, ref_s: str, c0: int, c1: int) -> float:
        """区间 [c0, c1] (双闭) 上 SMC-REF 一致率 (%)。

        双 gap 列 (SMC 与 REF 均为 '-') 不计入; 其余列计入分母,
        SMC == REF 计入分子。计入列数为 0 时返回 0.0。
        """
        match = 0
        n = 0
        for cj in range(c0, c1 + 1):
            sc = smc_s[cj] if cj < len(smc_s) else "-"
            rc = ref_s[cj]
            if sc == "-" and rc == "-":
                continue
            n += 1
            if sc == rc:
                match += 1
        return 100.0 * match / n if n else 0.0

    def _occurrence_stats(self, subreads: list[str], cols: list[int], u0: int) -> dict:
        """单个出现位点的覆盖统计。cols 为该区域对应的 ref 碱基列。

        cov 基于 MSA 完整跨度 [cols[0], cols[-1]]（含插入列）;
        col_depth 仍按 ref 碱基列, 以与 region offset 对齐。
        """
        region_len = len(cols)
        c0, c1 = cols[0], cols[-1]
        span_len = c1 - c0 + 1
        n_sub = len(subreads)
        cov_counts = []
        col_depth = [0] * region_len
        for sr in subreads:
            cov = 0
            for i, cj in enumerate(cols):
                if cj < len(sr) and sr[cj] != "-":
                    cov += 1
                    col_depth[i] += 1
            # 完整 MSA 跨度上的覆盖: 含插入列 (subread 插入碱基计入)
            cov += sum(1 for cj in range(c0, c1 + 1)
                       if cj not in cols and cj < len(sr) and sr[cj] != "-")
            cov_counts.append(cov)

        cov_pcts = [100.0 * c / span_len for c in cov_counts]
        n_spanning = sum(1 for c in cov_counts if c == span_len)
        n_supporting = sum(
            1 for c, p in zip(cov_counts, cov_pcts)
            if c > 0 and p >= self.support_cov)

        return {
            "u_start": u0,
            "u_end": u0 + region_len,
            "msa_col_start": c0,
            "msa_col_end": c1,
            "msa_span_len": span_len,
            "n_subreads": n_sub,
            "n_spanning": n_spanning,
            "n_supporting": n_supporting,
            "supported": n_supporting >= self.support_reads,
            "cov_pct_min": round(min(cov_pcts), 1),
            "cov_pct_mean": round(sum(cov_pcts) / n_sub, 1),
            "cov_pct_max": round(max(cov_pcts), 1),
            "sub_cov_pct": [round(p, 1) for p in cov_pcts],
            "_sub_spans": [c == span_len for c in cov_counts],
            "col_depth": col_depth,
            "depth_min": min(col_depth),
            "depth_mean": round(sum(col_depth) / region_len, 1),
            "depth_max": max(col_depth),
        }

    def analyze_all(self) -> dict:
        """遍历所有记录, 搜索 region 并统计每个出现位点的 subread 支撑。"""
        region_len = len(self.region_seq)
        per_record = []
        all_occ = []
        n_records_no_subreads = 0
        n_records_excluded_smc_ref = 0

        for rec_idx, rec in enumerate(self.records):
            seqs = rec["msa_seqs"]
            if len(seqs) < 3:
                n_records_no_subreads += 1
                continue
            ref_s = seqs[2]
            subreads = [s for s in seqs[3:] if isinstance(s, str)]
            if not subreads:
                n_records_no_subreads += 1
                continue

            # ungapped ref index → MSA column
            u2col = [j for j, c in enumerate(ref_s) if c != "-"]
            ungapped = "".join(ref_s[j] for j in u2col)
            starts = self._find_all(ungapped, self.region_seq)
            if not starts:
                continue

            # 排除: 所有出现位点的 SMC-REF 一致率均 > 阈值 (SMC 支持 REF)
            smc_s = seqs[1]
            agrees = []
            for u0 in starts:
                cols = u2col[u0:u0 + region_len]
                agrees.append(
                    self._smc_agreement(smc_s, ref_s, cols[0], cols[-1]))
            if all(a > self.smc_match_thr for a in agrees):
                n_records_excluded_smc_ref += 1
                continue

            names = rec.get("names") or []
            rec_entry = {
                "rec_idx": rec_idx,
                "name": names[1] if len(names) > 1 else "?",
                "identity": rec.get("identity", 0.0),
                "occurrences": [],
            }
            for u0, agree in zip(starts, agrees):
                cols = u2col[u0:u0 + region_len]
                if len(cols) < region_len:
                    continue  # 防御: 不应发生
                occ = self._occurrence_stats(subreads, cols, u0)
                occ["smc_agree"] = round(agree, 1)
                rec_entry["occurrences"].append(occ)
                occ["_cols"] = cols  # 仅供样例展示, 报告中会移除
                all_occ.append(occ)
            per_record.append(rec_entry)

        return self._build_report(
            per_record, all_occ, n_records_no_subreads, n_records_excluded_smc_ref)

    # ---- report building -------------------------------------------------

    def _get_samples(self, per_record: list[dict], n: int = 20, flank: int = 3) -> list[dict]:
        """获取前 n 个出现位点的 MSA 片段样例。"""
        samples = []
        for rec_entry in per_record:
            if len(samples) >= n:
                break
            rec = self.records[rec_entry["rec_idx"]]
            seqs = rec["msa_seqs"]
            ref_s, smc_s = seqs[2], seqs[1]
            subreads = [s for s in seqs[3:] if isinstance(s, str)]

            for occ_idx, occ in enumerate(rec_entry["occurrences"]):
                if len(samples) >= n:
                    break
                cols = occ["_cols"]
                col_set = set(cols)
                c0, c1 = cols[0], cols[-1]
                s0 = max(0, c0 - flank)
                s1 = min(len(ref_s), c1 + flank + 1)
                marker = "".join(
                    "^" if j in col_set else " " for j in range(s0, s1))

                sub_rows = []
                for k, sr in enumerate(subreads[:5]):
                    sub_rows.append({
                        "name": f"sub_{k}",
                        "seq": sr[s0:s1] if s0 < len(sr) else "",
                        "cov_pct": occ["sub_cov_pct"][k],
                        "spans": occ["_sub_spans"][k],
                    })

                samples.append({
                    "rec_idx": rec_entry["rec_idx"],
                    "name": rec_entry["name"],
                    "identity": rec_entry["identity"],
                    "occ_idx": occ_idx,
                    "u_start": occ["u_start"],
                    "u_end": occ["u_end"],
                    "msa_col_start": occ["msa_col_start"],
                    "msa_col_end": occ["msa_col_end"],
                    "span_len": occ["msa_span_len"],
                    "n_subreads": occ["n_subreads"],
                    "n_spanning": occ["n_spanning"],
                    "n_supporting": occ["n_supporting"],
                    "supported": occ["supported"],
                    "cov_pct_mean": occ["cov_pct_mean"],
                    "smc_agree": occ["smc_agree"],
                    "marker": marker,
                    "ref_seg": ref_s[s0:s1],
                    "smc_seg": smc_s[s0:s1],
                    "subread_segs": sub_rows,
                })
        return samples

    def _build_report(
        self, per_record: list[dict], all_occ: list[dict],
        n_records_no_subreads: int, n_records_excluded_smc_ref: int,
    ) -> dict:
        region_len = len(self.region_seq)
        n_occ = len(all_occ)

        occ_per_rec = Counter(len(r["occurrences"]) for r in per_record)
        n_supported = sum(1 for o in all_occ if o["supported"])
        n_occ_any_span = sum(1 for o in all_occ if o["n_spanning"] > 0)
        spanning_dist = Counter(o["n_spanning"] for o in all_occ)

        # mean cov_pct histogram (per occurrence)
        n_bins, bin_width = 20, 5.0
        hist = [0] * n_bins
        for o in all_occ:
            hist[min(int(o["cov_pct_mean"] / bin_width), n_bins - 1)] += 1

        # depth profile aggregated by offset within region
        profile_values = [[] for _ in range(region_len)]
        for o in all_occ:
            for i, d in enumerate(o["col_depth"]):
                profile_values[i].append(d)
        offset_profile = []
        for i in range(region_len):
            vals = profile_values[i]
            offset_profile.append({
                "offset": i,
                "base": self.region_seq[i],
                "mean_depth": round(sum(vals) / n_occ, 2) if n_occ else 0.0,
                "max_depth": max(vals) if vals else 0,
                "min_depth": min(vals) if vals else 0,
                "median_depth": round(median(vals), 2) if vals else 0.0,
            })

        samples = self._get_samples(per_record, self.sample_n)
        for o in all_occ:
            o.pop("_cols", None)  # 内部字段, 不进报告
            o.pop("_sub_spans", None)

        return {
            "file": self.filepath,
            "region_seq": self.region_seq,
            "region_len": region_len,
            "support_cov": self.support_cov,
            "support_reads": self.support_reads,
            "smc_match_thr": self.smc_match_thr,
            "sample_n": self.sample_n,
            "identity_min_thr": self.identity_min_thr,
            "identity_thr": self.identity_thr,
            "n_records": len(self.records),
            "n_records_no_subreads": n_records_no_subreads,
            "n_records_excluded_smc_ref": n_records_excluded_smc_ref,
            "n_records_with_match": len(per_record),
            "total_occurrences": n_occ,
            "occurrences_per_record": dict(sorted(occ_per_rec.items())),
            "n_supported_occurrences": n_supported,
            "supported_pct": round(100 * n_supported / max(n_occ, 1), 1),
            "n_occurrences_any_spanning": n_occ_any_span,
            "spanning_distribution": dict(sorted(spanning_dist.items())),
            "mean_cov_histogram": hist,
            "offset_depth_profile": offset_profile,
            "records": per_record,
            "samples": samples,
        }

    # ---- text report -----------------------------------------------------

    @staticmethod
    def _ascii_histogram(bins: list[int], n_bins: int, title: str, max_bar_width: int = 50) -> str:
        """ASCII histogram (bar + count + cum_pct)。"""
        total = sum(bins)
        L = []
        L.append(f"\n{'=' * 70}")
        L.append(title)
        L.append(f"共 {total} 个出现位点, bin_width={100 // n_bins}%")
        max_count = max(bins) if bins else 1
        cum = 0.0
        for i in range(n_bins):
            lo = i * (100 // n_bins)
            hi = lo + (100 // n_bins)
            bar_len = round(bins[i] / max_count *
                            max_bar_width) if max_count > 0 else 0
            bar = "#" * bar_len
            cum += round(100 * bins[i] / max(total, 1), 1)
            L.append(
                f"  {lo:>3d}-{hi:<3d}% | {bar} ({bins[i]:>5d}, 累积={cum:5.1f}%)")
        return "\n".join(L)

    def text_report(self) -> str:
        """生成文本报告。"""
        r = self.analyze_all()
        n_occ = r["total_occurrences"]
        L = []

        # Header
        L.append("=" * 70)
        L.append("Reference region subread support — 分析")
        L.append("=" * 70)
        L.append(f"\n文件:       {self.filepath}")
        L.append(f"region:     {r['region_seq']} (len={r['region_len']})")
        L.append(
            f"support 判定: cov_pct >= {self.support_cov}% 的 subread >= {self.support_reads} 条")
        L.append(
            f"排除条件: 所有出现位点 SMC-REF 一致率均 > {self.smc_match_thr}% 的 record")
        L.append(
            f"identity 范围: [{self.identity_min_thr}, {self.identity_thr})")
        L.append(f"记录数:     {r['n_records']} (无 subread 跳过: {r['n_records_no_subreads']})")

        # Match summary
        L.append("\n" + "=" * 70)
        L.append("匹配情况")
        L.append("=" * 70)
        L.append(f"SMC 支持 REF 排除:  {r['n_records_excluded_smc_ref']}")
        L.append(f"含匹配的 record 数: {r['n_records_with_match']}")
        L.append(f"总出现位点数:       {n_occ}")
        L.append(f"每 record 出现数分布: {r['occurrences_per_record']}")

        if n_occ == 0:
            L.append("\n未在任何 record 的 reference 中找到该区域。")
            return "\n".join(L)

        # Support summary
        L.append("\n" + "=" * 70)
        L.append("支撑统计 (出现位点级)")
        L.append("=" * 70)
        L.append(
            f"supported 位点:      {r['n_supported_occurrences']} ({r['supported_pct']}%)")
        L.append(f"有完整跨越的位点:    {r['n_occurrences_any_spanning']}")
        L.append(f"n_spanning 分布:     {r['spanning_distribution']}")
        L.append(self._ascii_histogram(
            r["mean_cov_histogram"], 20, "出现位点平均覆盖率 (mean_cov_pct) 分布"))

        # Offset depth profile
        L.append("\n" + "=" * 70)
        L.append("区域内 offset 聚合 depth profile (所有出现位点平均)")
        L.append("=" * 70)
        prof = r["offset_depth_profile"]
        worst = min(prof, key=lambda e: e["mean_depth"])
        L.append(f"\n最低 depth: offset {worst['offset']} "
                 f"(base={worst['base']}, mean_depth={worst['mean_depth']})")
        L.append(f"\n{'offset':>6} {'base':>4} {'mean_depth':>11} "
                 f"{'max_depth':>9} {'min_depth':>9} {'median_depth':>12}")
        for e in prof:
            L.append(f"{e['offset']:>6} {e['base']:>4} {e['mean_depth']:>11.2f} "
                     f"{e['max_depth']:>9} {e['min_depth']:>9} "
                     f"{e['median_depth']:>12.2f}")

        # Samples
        if r["samples"]:
            L.append("\n" + "=" * 70)
            L.append(f"样例 (前 {self.sample_n} 个出现位点)")
            L.append("=" * 70)
            for s in r["samples"]:
                L.append(
                    f"\n  [rec {s['rec_idx']}] {s['name']} id={s['identity']:.6f} "
                    f"occ#{s['occ_idx']} u=[{s['u_start']},{s['u_end']}) "
                    f"msa_cols=[{s['msa_col_start']},{s['msa_col_end']}] "
                    f"span={s['span_len']}")
                L.append(
                    f"    n_subreads={s['n_subreads']} n_spanning={s['n_spanning']} "
                    f"n_supporting={s['n_supporting']} supported={s['supported']} "
                    f"mean_cov={s['cov_pct_mean']}% smc_agree={s['smc_agree']}%")
                L.append(f"{'':>7}  {s['marker']}")
                L.append(f"{'ref':>7}: {s['ref_seg']}")
                L.append(f"{'smc':>7}: {s['smc_seg']}")
                for row in s["subread_segs"]:
                    tag = "span" if row["spans"] else ""
                    L.append(
                        f"{row['name']:>7}: {row['seq']}  cov={row['cov_pct']}% {tag}")

        # Summary
        L.append("\n" + "=" * 70)
        L.append("总结")
        L.append("=" * 70)
        L.append(f"\nregion '{r['region_seq']}' 共出现 {n_occ} 次 "
                 f"({r['n_records_with_match']} 条 record)")
        L.append(
            f"supported (cov>={self.support_cov}% subread>={self.support_reads} 条): "
            f"{r['n_supported_occurrences']} ({r['supported_pct']}%)")
        L.append(f"有 subread 完整跨越的位点: {r['n_occurrences_any_spanning']}")
        L.append(f"区域内最低平均 depth: offset {worst['offset']} "
                 f"mean_depth={worst['mean_depth']}")

        return "\n".join(L)

    # ---- json report -----------------------------------------------------

    def json_report(self) -> str:
        """生成 JSON 报告。"""
        r = self.analyze_all()
        return json.dumps(r, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reference region subread support — 检查 ref substr 下方是否有 subreads 支撑",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", help="ASRTC JSONL 文件路径 (.asrtc.txt)")
    parser.add_argument("--region-seq", required=True,
                        help="要在 reference 中搜索的序列字符串 (不区分大小写, 精确匹配)")
    parser.add_argument("--support-cov", type=float, default=50.0,
                        help="per-subread 区域覆盖率%% >= 此值才算 supporting (default: 50)")
    parser.add_argument("--support-reads", type=int, default=1,
                        help="supporting subreads >= N 时, 该出现位点判为 supported (default: 1)")
    parser.add_argument("--smc-match-thr", type=float, default=90.0,
                        help="所有出现位点 SMC-REF 一致率均 > 此值(%%) 的 record 被排除 "
                             "(default: 90; 设为 101 可关闭)")
    parser.add_argument("--sample-n", type=int, default=20,
                        help="样例出现位点个数 (default: 20)")
    parser.add_argument("--output", choices=["text", "json"], default="text",
                        help="输出格式 (default: text)")
    parser.add_argument("--out-file", metavar="PATH", default=None,
                        help="输出文件路径 (default: stdout)")
    parser.add_argument("--identity-thr", type=float, default=1.01,
                        help="identity 上限, >=此值的记录被排除 (default: 1.01, 即不过滤)")
    parser.add_argument("--identity-min", type=float, default=0.0,
                        help="identity 下限, <此值的记录被排除 (default: 0.0)")
    return parser


def main(args=None) -> None:
    opts = build_parser().parse_args(args)

    region = opts.region_seq.strip().upper()
    if not region:
        print("Error: --region-seq 不能为空", file=sys.stderr)
        sys.exit(1)

    analyzer = RegionSupportAnalyzer(
        opts.input, region,
        identity_thr=opts.identity_thr, identity_min_thr=opts.identity_min,
        support_cov=opts.support_cov, support_reads=opts.support_reads,
        smc_match_thr=opts.smc_match_thr, sample_n=opts.sample_n)

    report = analyzer.json_report() if opts.output == "json" else analyzer.text_report()

    if opts.out_file:
        with open(opts.out_file, "w") as fh:
            fh.write(report + "\n")
        print(f"Report written to {opts.out_file}")
    else:
        print(report)


if __name__ == "__main__":
    main()
