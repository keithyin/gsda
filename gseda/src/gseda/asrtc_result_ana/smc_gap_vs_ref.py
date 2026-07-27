#!/usr/bin/env python3
"""SMC gap vs reference base 差异分析 — smc=`-` 且 ref 有 base 的场景。

分类逻辑:
  1. 真 deletion (≥70% subreads 在此位置也有 gap): SMC gap 正确, ref 多碱基
  2. SMC gap artifact (<70% subreads 有 gap):
     a) shift 模式 — 前 SHIFT_WINDOW bp 内 ≥3 个相邻列 ref base / smc gap
     b) local gap interference — 距参考序列已有 gap ≤5bp
     c) 残余异常 — 以上两者都不是

分布图:
    python smc_gap_vs_ref.py input.asrtc.txt                       # ASCII histogram in text report (overall + per-HP)
    python smc_gap_vs_ref.py input.asrtc.txt --plot-format png     # save PNG multi-subplot to current dir

按 homopolymer 长度分组的分布:
    python smc_gap_vs_ref.py input.asrtc.txt --output json -o out.json  # "hp_grouped_distribution" in JSON

用法:
    python smc_gap_vs_ref.py input.asrtc.txt
    python smc_gap_vs_ref.py input.asrtc.txt --output json -o result.json
    python smc_gap_vs_ref.py input.asrtc.txt --identity-thr 0.98 --identity-min 0.95
"""
from __future__ import annotations

import json
import sys
import argparse
from collections import Counter as C
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# bp: +2 buffer → total threshold = _GAP_WINDOW+2=5bp, local interference
_GAP_WINDOW = 3
_HP_THRESHOLD = 2        # homopolymer run ≥ this → HP context
_STRONG_PCT = 70         # ≥70% subreads have gap → true deletion
_SHIFT_WINDOW = 10       # bp: 偏移证据检测窗口 (前 N bp 内 ref base / smc gap)
_MIN_SHIFT_EVIDENCE = 3  # ≥3 个相邻位置 → offset 确认


class SmcGapAnalyzer:
    """分析 msa_seqs 中 smc_seq=`-` 但 reference_seq 有 base 的位点。"""

    def __init__(self, filepath: str, identity_thr: float = 1.0, identity_min_thr: float = 0.95):
        self.filepath = filepath
        self.identity_thr = identity_thr
        self.identity_min_thr = identity_min_thr
        self.records: list[dict] = []
        # sub_gap_pct (0-100) for all analyzed positions
        self._gap_ratios: list[float] = []
        # "no_HP" or HP_length (int) → [gap_ratios]
        self._gap_by_hp: dict[str | int, list[float]] = {}
        self._load()

    def _reset(self) -> None:
        """Reset runtime state between classify_all calls."""
        self._gap_ratios.clear()
        self._gap_by_hp.clear()

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

        except Exception as e:
            raise ValueError(f"line_no:{line_no} line={line}")

    # ---- helpers ---------------------------------------------------------

    @staticmethod
    def _count_gaps_before(seq: str, stop: int) -> int:
        """seq[:stop] 中 gap 的数量。"""
        return sum(1 for i in range(stop) if seq[i] == "-")

    @staticmethod
    def _max_hp_in_context(ctx: str) -> tuple[int, str]:
        """在一段序列上下文中查找最大 homopolymer 长度 (忽略 gap)。"""
        best = 0
        best_base = ""
        run = 0
        prev = ""
        for ch in ctx:
            if ch == "-":
                continue
            if ch == prev:
                run += 1
            else:
                if run > best:
                    best = run
                    best_base = prev
                run = 1
                prev = ch
        if run > best:
            best = run
            best_base = prev
        return (best, best_base) if best >= _HP_THRESHOLD else (0, "")

    @staticmethod
    def _local_hp_at(seq: str, pos: int) -> int:
        """计算 seq[pos] 所在的最大同聚物长度（忽略 gap）。"""
        target = seq[pos]
        if target == "-":
            return 0
        left = pos - 1
        cnt = 1
        while left >= 0 and (seq[left] in [target, '-']):
            if seq[left] == target:
                cnt += 1
            left -= 1
            
        right = pos + 1
        while right < len(seq) and (seq[right] in [target, '-']):
            if seq[right] == target:
                cnt += 1
            right += 1
            
        return cnt

    # ---- classification --------------------------------------------------

    def classify_all(self) -> dict:
        """对每条 smc_gap/ref_base 位点分类, 返回汇总+详情。"""
        self._reset()
        total = 0
        true_deletion = []       # cat 1: ≥70% subreads also have gap
        artifact_offset = []     # cat 2a
        artifact_local_gap = []  # cat 2b
        residual_anomaly = []    # cat 2c

        for rec_idx, rec in enumerate(self.records):
            seqs = rec["msa_seqs"]
            if len(seqs) < 3:
                continue
            smc_s, ref_s = seqs[1], seqs[2]
            sbr_seq = [s for s in seqs[3:] if isinstance(s, str)]
            n_sbr = len(sbr_seq)
            if not n_sbr:
                continue
            min_len = min(len(smc_s), len(ref_s))

            # Precompute gap info for this record
            ref_gap_set = {j for j, c in enumerate(ref_s) if c == "-"}

            for pos in range(min_len):
                if smc_s[pos] != "-" or ref_s[pos] == "-":
                    continue  # not our target case

                total += 1
                ref_base = ref_s[pos]

                # Count subread gaps and bases at this position
                sub_gaps = 0
                for sr in sbr_seq:
                    if pos < len(sr) and sr[pos] == "-":
                        sub_gaps += 1
                pct_gap = sub_gaps / n_sbr * 100

                # Collect gap ratio for global + HP-grouped distributions
                self._gap_ratios.append(pct_gap)
                hp = self._local_hp_at(ref_s, pos)
                key = "no_HP" if hp < _HP_THRESHOLD else hp
                key = f"{ref_s[pos]}({key})"
                self._gap_by_hp.setdefault(key, []).append(pct_gap)

                # ---- Category 1: true deletion ----
                if pct_gap >= _STRONG_PCT:
                    hp_len, hp_base = self._max_hp_in_context(
                        ref_s[max(0, pos - 8):pos + 9])
                    true_deletion.append({
                        "rec_idx": rec_idx,
                        "pos": pos,
                        "ref_base": ref_base,
                        "sub_gap_pct": round(pct_gap, 1),
                        "n_sub_gap": sub_gaps,
                        "hp_length": hp_len,
                        "hp_base": hp_base,
                    })
                    continue

                # ---- Category 2: SMC gap artifact ----

                # Sub-check a: alignment offset? Detect from adjacent MSA column content.
                # Count positions in preceding SHIFT_WINDOW where ref_s has base and smc_s has gap.
                shift_evidence = 0
                for j in range(max(0, pos - _SHIFT_WINDOW), pos):
                    if ref_s[j] != "-" and smc_s[j] == "-":
                        shift_evidence += 1

                # Sub-check b: distance to nearest local ref gap
                dist_to_ref_gap = min((abs(pos - g)
                                      for g in ref_gap_set), default=999)

                if shift_evidence >= _MIN_SHIFT_EVIDENCE:
                    artifact_offset.append({
                        "rec_idx": rec_idx,
                        "pos": pos,
                        "ref_base": ref_base,
                        "sub_gap_pct": round(pct_gap, 1),
                        "n_sub_gap": sub_gaps,
                        "shift_evidence": shift_evidence,
                    })
                elif dist_to_ref_gap <= _GAP_WINDOW + 2:
                    # Close to a local ref gap → local interference (no offset)
                    artifact_local_gap.append({
                        "rec_idx": rec_idx,
                        "pos": pos,
                        "ref_base": ref_base,
                        "sub_gap_pct": round(pct_gap, 1),
                        "n_sub_gap": sub_gaps,
                        "dist_to_ref_gap": dist_to_ref_gap,
                    })
                else:
                    # Sub-check c: HP context?
                    hp_len, hp_base = self._max_hp_in_context(
                        ref_s[max(0, pos - 8):pos + 9])
                    residual_anomaly.append({
                        "rec_idx": rec_idx,
                        "pos": pos,
                        "ref_base": ref_base,
                        "sub_gap_pct": round(pct_gap, 1),
                        "n_sub_gap": sub_gaps,
                        "hp_length": hp_len,
                        "hp_base": hp_base,
                    })

        return self._build_report(total, true_deletion, artifact_offset, artifact_local_gap, residual_anomaly)

    # ---- report building -------------------------------------------------

    def _build_hp_grouped(self, n_bins: int, total: int) -> dict:
        """构建按 homopolymer 长度分组的 gap_ratio 分布。"""
        result: dict[str | int, dict] = {}
        for key in sorted(
            self._gap_by_hp.keys(),
            # key=lambda k: (0, k) if k == "no_HP" else (1, int(k)),
        ):
            ratios = self._gap_by_hp[key]
            if not ratios:
                continue
            count_dist = [0] * n_bins
            for v in ratios:
                b = min(int(v / 5.0), n_bins - 1)
                count_dist[b] += 1
            count_pct = [round(100 * c / max(len(ratios), 1), 1)
                         for c in count_dist]
            cum = 0.0
            cum_pct = []
            for p in count_pct:
                cum += p
                cum_pct.append(round(cum, 1))
            result[key] = {
                "n_positions": len(ratios),
                "pct_total": round(100 * len(ratios) / max(total, 1), 1),
                "min": round(min(ratios), 1),
                "max": round(max(ratios), 1),
                "mean": round(sum(ratios) / len(ratios), 1),
                "bins": count_dist,
                "bins_pct": count_pct,
                "cum_pct": cum_pct,
            }
        return result

    def _build_report(
        self, total: int,
        true_del: list, offset_art: list, local_gap_art: list, residual: list,
    ) -> dict:
        """构建报告数据结构。"""

        # gap_ratio distribution (all analyzed positions)
        n_bins = 20
        bin_width = 5.0
        gap_dist = [0] * n_bins
        for v in self._gap_ratios:
            b = min(int(v / bin_width), n_bins - 1)
            gap_dist[b] += 1
        gap_dist_pct = [round(100 * c / max(total, 1), 1) for c in gap_dist]

        # Cumulative probability (running total of percentages)
        cum = 0.0
        gap_dist_cum_pct = []
        for p in gap_dist_pct:
            cum += p
            gap_dist_cum_pct.append(round(cum, 1))

        # True deletion stats by ref base
        del_by_base = C(d["ref_base"] for d in true_del)
        del_hp_count = sum(1 for d in true_del if d["hp_length"])
        del_no_hp_count = len(true_del) - del_hp_count

        # Offset artifact: breakdown by shift_evidence (相邻位置 ref base / smc gap 数)
        offset_by_shift = C(a["shift_evidence"] for a in offset_art)

        # Residual: by ref base and HP status
        res_by_base = C(r["ref_base"] for r in residual)
        res_hp_count = sum(1 for r in residual if r["hp_length"])
        res_no_hp_count = len(residual) - res_hp_count

        return {
            "file": self.filepath,
            "identity_min_thr": self.identity_min_thr,
            "identity_thr": self.identity_thr,
            "n_records": len(self.records),
            "total_smc_gap_ref_base": total,
            # Category 1: True deletion
            "true_deletion": {
                "count": len(true_del),
                "pct": round(100 * len(true_del) / max(total, 1), 1),
                "by_ref_base": dict(del_by_base.most_common()),
                "in_hp_context": del_hp_count,
                "not_in_hp": del_no_hp_count,
                "samples": self._get_samples(true_del, 5, show_subreads=True),
            },
            # Category 2a: Alignment offset
            "artifact_offset": {
                "count": len(offset_art),
                "pct": round(100 * len(offset_art) / max(total, 1), 1),
                "by_shift_evidence": dict(sorted(offset_by_shift.items())),
                "samples": self._get_samples(offset_art, 5),
            },
            # Category 2b: Local gap interference
            "artifact_local_gap": {
                "count": len(local_gap_art),
                "pct": round(100 * len(local_gap_art) / max(total, 1), 1),
                "by_ref_base": dict(C(a["ref_base"] for a in local_gap_art).most_common()),
                "samples": self._get_samples(local_gap_art, 5),
            },
            # Category 2c: Residual anomalies
            "residual_anomalies": {
                "count": len(residual),
                "pct": round(100 * len(residual) / max(total, 1), 1),
                "by_ref_base": dict(res_by_base.most_common()),
                "in_hp_context": res_hp_count,
                "not_in_hp": res_no_hp_count,
                "samples": self._get_samples(residual, 5),
            },
            # gap_ratio distribution across all analyzed positions
            "gap_ratio_distribution": {
                "n_positions": len(self._gap_ratios),
                "min": round(min(self._gap_ratios), 1) if self._gap_ratios else 0,
                "max": round(max(self._gap_ratios), 1) if self._gap_ratios else 0,
                "mean": round(sum(self._gap_ratios) / len(self._gap_ratios), 1) if self._gap_ratios else 0,
                "bins": gap_dist,
                "bins_pct": gap_dist_pct,
                "cum_pct": gap_dist_cum_pct,  # 累积概率分布 (running total)
            },
            # Per-homopolymer-length gap_ratio distribution
            "hp_grouped_distribution": self._build_hp_grouped(n_bins, total),
        }

    def _get_samples(self, items: list, n: int = 5, show_subreads: bool = False) -> list[dict]:
        """获取样例详情。"""
        samples = []
        shown = 0
        for item in items:
            if shown >= n:
                break
            rec = self.records[item["rec_idx"]]
            seqs = rec["msa_seqs"]
            smc_s, ref_s = seqs[1], seqs[2]
            pos = item["pos"]

            ctx_start = max(0, pos - 8)
            ctx_end = min(len(ref_s), pos + 9)

            sample: dict = {
                "rec_idx": item["rec_idx"],
                "name": rec["names"][1],
                "identity": rec["identity"],
                "pos": pos,
                "ref_base": item["ref_base"],
                "ref_ctx": f"{ref_s[ctx_start:pos]}[{ref_s[pos]}]{ref_s[pos+1:ctx_end]}",
                "smc_ctx": f"{smc_s[ctx_start:pos]}[-]{smc_s[pos+1:ctx_end]}",
            }

            if show_subreads:
                sbr_seq = [s for s in seqs[3:] if isinstance(s, str)]
                sub_bases = []
                for sr in sbr_seq:
                    if pos < len(sr):
                        sub_bases.append(sr[pos])
                sample["sub_gap_pct"] = item.get("sub_gap_pct", 0)
                sample["subreads_at_pos"] = dict(C(sub_bases).most_common())
                sample["n_subreads"] = len(sbr_seq)
            else:
                sample["sub_gap_pct"] = item.get("sub_gap_pct", 0)
                if "shift_evidence" in item:
                    sample["shift_evidence"] = item["shift_evidence"]
                if "dist_to_ref_gap" in item:
                    sample["dist_to_ref_gap"] = item["dist_to_ref_gap"]
                if "hp_length" in item and item.get("hp_length"):
                    sample["hp_length"] = item["hp_length"]
                    sample["hp_base"] = item["hp_base"]

            samples.append(sample)
            shown += 1
        return samples

    # ---- gap_ratio distribution ------------------------------------------

    @staticmethod
    def _ascii_histogram(bins: list[int], n_bins: int, title: str = "gap_ratio distribution", max_bar_width: int = 50) -> str:
        """Generate ASCII histogram (bar + count + cum_pct)."""
        # Compute cum_pct inline for display
        total = sum(bins)
        L = []
        L.append(f"\n{'=' * 70}")
        L.append(title)
        L.append(f"共 {total} 个位点, bin_width={100 // n_bins}%")
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

    def write_gap_ratio_plot(self, filepath: str) -> None:
        """Generate and save gap ratio distribution as PNG with HP-grouped subplots."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: PLC0414

        n_bins = 20
        bin_width = 5.0
        N_COLS = 3

        overall_counts = [0] * n_bins
        for v in self._gap_ratios:
            b = min(int(v / bin_width), n_bins - 1)
            overall_counts[b] += 1

        keys = sorted(self._gap_by_hp.keys(), 
                    #   key=lambda k: (0, k) if k == "no_HP" else (1, int(k))
                      )
        N_TOTAL = N_COLS + len(keys)
        N_ROWS = (N_TOTAL - 1 + N_COLS) // N_COLS 

        fig_width = N_COLS * 5.0       # wider columns → more space for tick labels
        fig_height = 4.0 + N_ROWS * 3.0  # taller rows → clearer labels & titles
        fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(
            fig_width, fig_height), squeeze=False)
        bin_labels = [f"{i * 5}-{(i + 1) * 5}%" for i in range(n_bins)]
        # Subsample ticks: show every 4th bin (6 labels) to avoid crowding
        _tick_idx = list(range(0, n_bins, 2))   # every 2nd bin → 10% steps: [0%, 10%, ..., 90%]
        _tick_lbl = [f"{_tick_idx[i] * 5}%" for i in range(len(_tick_idx))]

        # Overall (first subplot, left panel)
        ax_main = axes[0, 0]
        ax_cum = ax_main.twinx()
        ax_main.bar(bin_labels, overall_counts,
                    color="steelblue", edgecolor="white")
        ax_main.set_ylabel("Count", fontsize=9)
        ax_main.axvline(x=6.5, color="green", linestyle="--",
                        linewidth=1.2, label=f"≥{_STRONG_PCT}% 阈值")
        ax_main.legend(fontsize=8)

        cum_overall = [round(sum(overall_counts[:i + 1]) /
                             max(sum(overall_counts), 1) * 100, 1) for i in range(n_bins)]
        ax_cum.plot(bin_labels, cum_overall, color="crimson",
                    marker="o", markersize=3, label="Cumulative %")
        ax_cum.set_ylabel("Cumulative %", fontsize=9)
        ax_cum.legend(fontsize=8)

        # Rotate x-ticks for readability
        for _ax in (ax_main, ax_cum):
            _ax.set_xticks(_tick_idx)
            _ax.set_xticklabels(_tick_lbl, rotation=30, ha='right', fontsize=8)

        overall = len(self._gap_ratios)
        min_val = round(min(self._gap_ratios), 1) if self._gap_ratios else 0
        mean_val = round(sum(self._gap_ratios) /
                         max(overall, 1), 1) if overall else 0
        max_val = round(max(self._gap_ratios), 1) if self._gap_ratios else 0
        ax_main.set_title(
            f"Overall (n={overall}, min={min_val}, mean={mean_val}, max={max_val})", fontsize=10)

        # Per HP group subplots
        for idx, key in enumerate(keys):
            row = -(-(idx + 1) // N_COLS)  # ceil division from 1-based
            col = (idx) % N_COLS
            ax_main = axes[row, col] if row * N_COLS + col < N_ROWS * \
                N_COLS else fig.add_subplot(N_ROWS, N_COLS, row * N_COLS + col + 1)

            ratios = self._gap_by_hp[key]
            counts = [0] * n_bins
            for v in ratios:
                b = min(int(v / bin_width), n_bins - 1)
                counts[b] += 1

            ax_cum = ax_main.twinx()
            ax_main.bar(bin_labels, counts, color="teal", edgecolor="white")
            ax_main.set_ylabel("Count", fontsize=9)

            cum_val = [round(sum(counts[:i + 1]) / max(sum(counts), 1) * 100, 1)
                       for i in range(n_bins)]
            ax_cum.plot(bin_labels, cum_val, color="crimson",
                        marker="o", markersize=3)
            ax_cum.set_ylabel("Cumulative %", fontsize=9)

            # Rotate x-ticks for readability
            for _ax in (ax_main, ax_cum):
                _ax.set_xticks(_tick_idx)
                _ax.set_xticklabels(_tick_lbl, rotation=30, ha='right', fontsize=8)

            mean_g = sum(ratios) / max(len(ratios), 1)
            ax_main.set_title(
                f"HP {str(key)}bp (n={len(ratios)}, mean_gap_ratio={mean_g:.1f}%)", fontsize=10)

        # Prevent empty trailing space from autoscaling
        for row in axes:
            for ax in row:
                ax.set_xlim(-0.5, n_bins - 0.5)

        fig.suptitle("Gap Ratio Distribution — Homopolymer Grouped",
                     fontsize=12, y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(filepath, dpi=150)
        plt.close(fig)

    # ---- text report -----------------------------------------------------

    def text_report(self, top_k: int = 5) -> str:
        """生成文本报告。"""
        r = self.classify_all()
        total = r["total_smc_gap_ref_base"]
        L = []

        def hdr(title):
            L.append("\n" + "=" * 70)
            L.append(title)
            L.append("=" * 70)

        def sub(label, val):
            L.append(f"{label}: {val}")

        # Header
        L.append("=" * 70)
        L.append("SMC gap vs reference base — 综合分析")
        L.append("=" * 70)
        L.append(f"\n文件:    {self.filepath}")
        L.append(
            f"identity 范围: [{self.identity_min_thr}, {self.identity_thr})")
        L.append(f"记录数:  {r['n_records']}")
        L.append(f"\nsmc=`-` / ref有base 总位点数: {total}\n")

        # Category 1: True deletion
        hdr(f"CATEGORY A: 真 deletion (≥{_STRONG_PCT}% subreads 也有 gap — SMC gap 正确)")
        td = r["true_deletion"]
        L.append(f"\n数量: {td['count']} ({td['pct']}%)")
        L.append(f"按 ref base 分布: {td['by_ref_base']}")
        L.append(
            f"同聚物区域 (HP≥5): {td['in_hp_context']} | 非 HP: {td['not_in_hp']}")

        if td["samples"]:
            L.append("\nSample:")
            for s in td["samples"]:
                sub_reads_str = ""
                if "subreads_at_pos" in s:
                    sub_reads_str = f", subreads={s['subreads_at_pos']} ({s['n_subreads']}总)"
                L.append(
                    f"  [{s['rec_idx']}] {s['name']} id={s['identity']:.6f} pos={s['pos']}: ref_base={s['ref_base']}, sub_gap_pct={s['sub_gap_pct']}%{sub_reads_str}")
                L.append(f"    ref: {s['ref_ctx']}")
                L.append(f"    smc: {s['smc_ctx']}")

        # Category 2a: Alignment offset
        hdr("CATEGORY B: SMC gap artifact — shift 模式 (相邻列 ref base / smc gap ≥3)")
        ao = r["artifact_offset"]
        L.append(f"\n数量: {ao['count']} ({ao['pct']}%)")
        L.append(f"shift_evidence 分布: {ao['by_shift_evidence']}")

        if ao["samples"]:
            L.append("\nSample:")
            for s in ao["samples"]:
                L.append(f"  [{s['rec_idx']}] {s['name']} id={s['identity']:.6f} pos={s['pos']}: ref_base={s['ref_base']}, sub_gap_pct={s['sub_gap_pct']}%, shift_evidence={s.get('shift_evidence','?')}")
                L.append(f"    ref: {s['ref_ctx']}")
                L.append(f"    smc: {s['smc_ctx']}")

        # Category 2b: Local gap interference
        hdr("CATEGORY C: SMC gap artifact — local gap interference (距参考已有 gap ≤5bp)")
        lg = r["artifact_local_gap"]
        L.append(f"\n数量: {lg['count']} ({lg['pct']}%)")
        L.append(f"按 ref base 分布: {lg['by_ref_base']}")

        if lg["samples"]:
            L.append("\nSample:")
            for s in lg["samples"]:
                L.append(f"  [{s['rec_idx']}] {s['name']} id={s['identity']:.6f} pos={s['pos']}: ref_base={s['ref_base']}, sub_gap_pct={s['sub_gap_pct']}%, dist_ref_gap={s['dist_to_ref_gap']}")
                L.append(f"    ref: {s['ref_ctx']}")
                L.append(f"    smc: {s['smc_ctx']}")

        # Category 2c: Residual anomalies
        hdr(f"CATEGORY D: 残余异常 (<{_STRONG_PCT}% sub gap, 非 offset, 非 local gap)")
        ra = r["residual_anomalies"]
        L.append(f"\n数量: {ra['count']} ({ra['pct']}%)")
        L.append(f"按 ref base 分布: {ra['by_ref_base']}")
        L.append(
            f"同聚物区域 (HP≥5): {ra['in_hp_context']} | 非 HP: {ra['not_in_hp']}")

        if ra["samples"]:
            L.append("\nSample:")
            for s in ra["samples"]:
                hp_str = f", HP={s.get('hp_length', 0)}" if s.get(
                    "hp_length") else ""
                L.append(
                    f"  [{s['rec_idx']}] {s['name']} id={s['identity']:.6f} pos={s['pos']}: ref_base={s['ref_base']}, sub_gap_pct={s['sub_gap_pct']}%{hp_str}")
                L.append(f"    ref: {s['ref_ctx']}")
                L.append(f"    smc: {s['smc_ctx']}")

        # Gap ratio distribution (overall)
        L.append(self._ascii_histogram(
            r["gap_ratio_distribution"]["bins"], 20, "Gap ratio 分布 (所有 smc=`-`/ref有base 位点)"))

        # Per-homopolymer-length gap_ratio distribution
        hp_dist = r.get("hp_grouped_distribution", {})
        if hp_dist:
            L.append("\n" + "=" * 70)
            L.append(f"Gap ratio 分布 (按 Reference Homopolymer 长度分组)")
            L.append("=" * 70)
            # Summary table
            total_pos = r["total_smc_gap_ref_base"]
            L.append(
                f"\n{'Group':<12} {'Count':>6} | {'% of Total':>10} | {'Min':>5} {'Mean':>5} {'Max':>5}")
            for key in sorted(hp_dist.keys(),
                            #   key=lambda k: (0, k) if k == "no_HP" else (1, int(k))
                              ):
                d = hp_dist[key]
                pct = round(100 * d['n_positions'] / max(total_pos, 1), 1)
                L.append(
                    f"{str(key):<12} {d['n_positions']:>6} | {pct:>10.1f} | {d['min']:>5.1f} {d['mean']:>5.1f} {d['max']:>5.1f}")

            # Detailed histograms per group
            for key in sorted(hp_dist.keys(), 
                            #   key=lambda k: (0, k) if k == "no_HP" else (1, int(k))
                              ):
                d = hp_dist[key]
                L.append(f"\n{'-' * 70}")
                L.append(self._ascii_histogram(
                    d["bins"], 20, f"Group: {str(key)}bp (n={d['n_positions']})"))

        # Summary
        L.append("\n\n" + "=" * 70)
        L.append("总结")
        L.append("=" * 70)
        L.append(f"\n总位点数:          {total}")
        L.append(
            f"A. 真 deletion:    {td['count']} ({td['pct']:>5.1f}%) ← SMC gap 可信, ref 多碱基")
        L.append(
            f"B. shift 模式:     {ao['count']} ({ao['pct']:>5.1f}%) ← 相邻列 ref base / smc gap ≥3, 比对偏移")
        L.append(
            f"C. local gap:      {lg['count']} ({lg['pct']:>5.1f}%) ← 参考序列已有 gap 的局部效应")
        L.append(f"D. 残余异常:       {ra['count']} ({ra['pct']:>5.1f}%) ← 需手动审查")

        return "\n".join(L)

    # ---- json report -----------------------------------------------------

    def json_report(self) -> str:
        """生成 JSON 报告。"""
        r = self.classify_all()
        return json.dumps(r, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="SMC gap vs reference base — smc=`-` / ref 有 base 场景分析",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", help="ASRTC JSONL 文件路径 (.asrtc.txt)")
    parser.add_argument("--output", "-o", choices=["text", "json"], default="text",
                        help="输出格式 (default: text)")
    parser.add_argument("--out-file", metavar="PATH", default=None,
                        help="输出文件路径 (default: stdout)")
    parser.add_argument("--identity-thr", type=float, default=0.99999,
                        help="identity 上限, >=此值的记录被排除 (default: 0.99999)")
    parser.add_argument("--identity-min", type=float, default=0.95,
                        help="identity 下限, <此值的记录被排除 (default: 0.95)")
    parser.add_argument("--plot-format", choices=["none", "ascii", "png"], default="ascii",
                        help="gap_ratio 分布图格式 (default: ascii)")
    return parser


def main(args=None) -> None:
    parser = build_parser()
    opts = parser.parse_args(args)

    analyzer = SmcGapAnalyzer(
        opts.input, identity_thr=opts.identity_thr, identity_min_thr=opts.identity_min)

    if opts.output == "json":
        report = analyzer.json_report()
    else:
        report = analyzer.text_report()

    # Handle PNG plot
    if opts.plot_format == "png":
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            print("Warning: matplotlib not installed, skipping PNG plot",
                  file=sys.stderr)
            opts.plot_format = "ascii"
        else:
            out_dir = str(Path(opts.out_file).parent) if opts.out_file else "."
            plot_path = Path(out_dir) / "gap_ratio_distribution.png"
            analyzer.write_gap_ratio_plot(str(plot_path))
            print(f"Plot saved to {plot_path}")

    if opts.out_file:
        with open(opts.out_file, "w") as fh:
            fh.write(report + "\n")
        print(f"Report written to {opts.out_file}")
    else:
        print(report)


if __name__ == "__main__":
    main()
