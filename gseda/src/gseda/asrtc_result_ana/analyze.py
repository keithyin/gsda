#!/usr/bin/env python3
"""ASRTC JSONL 综合分析 — SMC consensus vs reference 差异全维度分析。

涵盖6个分析角度 + SMC 错误溯源（比对伪影 vs 真实共识错误）：

1. 全局 mismatch 分类统计（subreads 投票 ref / smc / ambiguous）
2. SMC error bias + transition/transversion
3. Indel 分析（gap 位点分类、SMC gap bias、del_in_ref 分析）
4. smc_support 深度挖掘（潜在生物学差异候选）
5. Identity 分层对比（高/中/低 identity 各层独立统计）
6. SMC 错误溯源：比对伪影 vs 真实共识错误

用法:
    python analyze.py input.asrtc.txt [options]
    python analyze.py input.asrtc.txt --threshold 90
    python analyze.py input.asrtc.txt --output markdown -o report.md
"""
from __future__ import annotations

import json
import sys
import argparse
from collections import Counter as C


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_TRANSITIONS = frozenset(["GA", "AG", "CT", "TC"])
_GAP_WINDOW = 5  # bp: mismatch within ±this of any gap → likely artifact
_HP_THRESHOLD = 5  # homopolymer run length ≥ this → problematic context


# ---------------------------------------------------------------------------
# Core analysis class
# ---------------------------------------------------------------------------

class AsrtcAnalyzer:
    """对 ASRTC JSONL 文件执行全维度 SMC vs reference 差异分析。"""

    def __init__(self, filepath: str, identity_thr: float = 1.0, threshold: int = 75):
        """
        Args:
            filepath: JSONL 文件路径
            identity_thr: identity >= 此值的记录被排除 (default: 1.0)
            threshold: subreads 强支持判定阈值 (default: 75 %)
        """
        self.filepath = filepath
        self.identity_thr = identity_thr
        self.threshold = threshold
        self.records: list[dict] = []
        # Analysis results (populated by run())
        self._cat_ref: list[dict] = []
        self._cat_smc: list[dict] = []
        self._cat_ambig: list[dict] = []
        self._indel_smc: list[dict] = []
        self._indel_ref: list[dict] = []
        self._indel_ambig: list[dict] = []

    # ---- loading ---------------------------------------------------------

    def _load(self) -> None:
        """读取 JSONL 文件，排除 identity >= threshold 的记录。"""
        self.records = []
        with open(self.filepath, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                if rec.get("identity", 1.0) >= self.identity_thr:
                    continue
                self.records.append(rec)

    # ---- helpers ---------------------------------------------------------

    @staticmethod
    def _gap_positions(s: str) -> set[int]:
        """返回序列中所有 gap 位点的集合。"""
        return {i for i, ch in enumerate(s) if ch == "-"}

    @staticmethod
    def _min_distance_to_gap(pos: int, gaps: set[int]) -> int:
        """pos 到最近 gap 的距离；gaps 为空返回 inf。"""
        if not gaps:
            return 9999
        # 只检查最近的几个 gap 即可（有序集合的前后）
        min_d = 9999
        for g in gaps:
            d = abs(pos - g)
            if d < min_d:
                min_d = d
        return min_d

    def _check_hp_at(self, seqs: list[str], pos: int) -> tuple[int, str]:
        """检查所有序列在 pos 处是否有 homopolymer ≥ HP_THRESHOLD。

        Returns:
            (max_hp_len, base_char): 最大 HP 长度和对应碱基；若 < threshold 返回 (0, "")
        """
        best = 0
        best_base = ""
        for seq in seqs:
            if pos >= len(seq):
                continue
            base = seq[pos]
            # count run containing pos
            start = pos
            while start > 0 and seq[start - 1] == base:
                start -= 1
            end = pos
            while end + 1 < len(seq) and seq[end + 1] == base:
                end += 1
            hp_len = end - start + 1
            if hp_len > best:
                best = hp_len
                best_base = base
        return (best, best_base) if best >= _HP_THRESHOLD else (0, "")

    # ---- analysis dimensions 1-3: mismatch & indel classification --------

    def classify_mismatches_and_indels(self) -> None:
        """执行角度 1-3：mismatch 分类 + Indel 分类。

        结果写入 self._cat_ref / _cat_smc / _cat_ambig / _indel_*。
        """
        self._cat_ref = []
        self._cat_smc = []
        self._cat_ambig = []
        self._indel_smc = []
        self._indel_ref = []
        self._indel_ambig = []

        for rec_idx, rec in enumerate(self.records):
            seqs = rec["msa_seqs"]
            if len(seqs) < 3:
                continue
            smc_s = seqs[1]
            ref_s = seqs[2]
            sbr_seq = [s for s in seqs[3:] if isinstance(s, str)]
            min_len = min(len(smc_s), len(ref_s))

            # Collect gap positions once per record
            gaps_ref = self._gap_positions(ref_s)
            gaps_smc = self._gap_positions(smc_s)
            all_gaps = gaps_ref | gaps_smc

            for pos in range(min_len):
                rc = ref_s[pos]
                sc = smc_s[pos]

                # --- base mismatch (both are bases) ---
                if rc != "-" and sc != "-" and rc != sc:
                    non_gap = []
                    for sr in sbr_seq:
                        ch = sr[pos] if pos < len(sr) else None
                        if ch is not None and ch != "-":
                            non_gap.append(ch)

                    if non_gap:
                        n_ref_c = sum(1 for c in non_gap if c == rc)
                        n_smc_c = sum(1 for c in non_gap if c == sc)
                        pct_ref = n_ref_c / len(non_gap) * 100
                        pct_smc_v = n_smc_c / len(non_gap) * 100

                        entry: dict = {
                            "rec_idx": rec_idx,
                            "pos": pos,
                            "ref": rc,
                            "smc": sc,
                            "identity": rec.get("identity", 0),
                            "query_aln_len": rec.get("query_aln_len", 0),
                            "sub_counts": C(non_gap).most_common(),
                            "n_subreads": len(sbr_seq),
                        }

                        if pct_ref >= self.threshold and pct_smc_v < (100 - self.threshold):
                            entry["pct_ref"] = pct_ref
                            entry["category"] = "ref_support"
                            self._cat_ref.append(entry)
                        elif pct_smc_v >= self.threshold and pct_ref < (100 - self.threshold):
                            entry["pct_smc"] = pct_smc_v
                            entry["category"] = "smc_support"
                            self._cat_smc.append(entry)
                        else:
                            entry["pct_ref"] = pct_ref
                            entry["pct_smc"] = pct_smc_v
                            entry["category"] = "ambiguous"
                            self._cat_ambig.append(entry)

                # --- indel (one side is gap) ---
                elif (rc == "-") ^ (sc == "-"):  # exactly one is gap
                    sub_bases = []
                    for sr in sbr_seq:
                        ch = sr[pos] if pos < len(sr) else None
                        if ch is not None and ch != "-":
                            sub_bases.append(ch)
                    n_base = len(sub_bases)
                    pct = n_base / len(sbr_seq) * 100 if sbr_seq else 0

                    indel_entry: dict = {
                        "rec_idx": rec_idx,
                        "pos": pos,
                        "n_subreads": len(sbr_seq),
                        "identity": rec.get("identity", 0),
                        "sub_counts": C(sub_bases).most_common() if sub_bases else [],
                        "n_base_in_sbr": n_base,
                    }

                    if rc == "-" and sc != "-":
                        # del_in_ref: ref gap, SMC calls base
                        indel_entry["indel_type"] = "del_in_ref"
                        indel_entry["smc_base"] = sc
                        strong = pct >= self.threshold
                        weak = pct < (100 - self.threshold)
                        if strong and weak:
                            self._indel_smc.append(indel_entry)
                        elif not strong and weak:
                            self._indel_ref.append(indel_entry)
                        else:
                            self._indel_ambig.append(indel_entry)

                    elif sc == "-" and rc != "-":
                        # ins_in_smc: SMC gap, ref has base
                        indel_entry["indel_type"] = "ins_in_smc"
                        indel_entry["ref_base"] = rc
                        strong = pct >= self.threshold
                        weak = pct < (100 - self.threshold)
                        if strong and weak:
                            self._indel_smc.append(indel_entry)
                        elif not strong and weak:
                            self._indel_ref.append(indel_entry)
                        else:
                            self._indel_ambig.append(indel_entry)

    # ---- analysis dimension 4: smc_support deep dive ---------------------

    def _analyse_smc_support(self) -> dict:
        """角度 4：smc_support（subreads 支持 SMC）深度挖掘。

        返回包含 top-N records 及详细上下文的字典。
        """
        if not self._cat_smc:
            return {"total": 0, "by_record": {}, "top_records": []}

        # Group by record
        by_record: dict[int, list[dict]] = {}
        for e in self._cat_smc:
            by_record.setdefault(e["rec_idx"], []).append(e)

        # Top-N records by smc_support count
        ranked = sorted(by_record.items(), key=lambda x: len(x[1]), reverse=True)

        top_records: list[dict] = []
        for rec_idx, entries in ranked[:20]:
            rec = self.records[rec_idx]
            seqs = rec["msa_seqs"]
            smc_s, ref_s = seqs[1], seqs[2]
            sbr_seq = [s for s in seqs[3:] if isinstance(s, str)]

            # Error matrix for this record
            err_matrix: dict[str, int] = {}
            for e in entries:
                key = f"{e['ref']}->{e['smc']}"
                err_matrix[key] = err_matrix.get(key, 0) + 1

            # Transition / Transversion
            transitions = sum(1 for e in entries if (e["ref"] + e["smc"]) in _TRANSITIONS)
            transversions = len(entries) - transitions

            # Average subread confidence
            avg_pct = sum(e.get("pct_smc", 0) for e in entries) / len(entries) if entries else 0

            top_records.append({
                "rec_idx": rec_idx,
                "identity": rec["identity"],
                "n_mismatches": len(entries),
                "aln_len": rec.get("query_aln_len", 0),
                "avg_subread_confidence": round(avg_pct, 1),
                "transitions": transitions,
                "transversions": transversions,
                "error_matrix": err_matrix,
            })

        return {"total": len(self._cat_smc), "by_record_count": len(by_record), "top_records": top_records}

    # ---- analysis dimension 5: identity stratified comparison ------------

    def _stratify_by_identity(self) -> dict:
        """角度 5：按 identity 分层对比。

        Returns:
            {low: {...}, mid: {...}, high: {...}} where low < 0.94, mid 0.94-0.98, high >= 0.98
        """
        strata = {"low": [], "mid": [], "high": []}
        # Map each entry to its stratum by looking up identity from the record
        for e in self._cat_ref + self._cat_smc + self._cat_ambig:
            rid = e["rec_idx"]
            ident = self.records[rid]["identity"]
            if ident < 0.94:
                strata["low"].append(e)
            elif ident < 0.98:
                strata["mid"].append(e)
            else:
                strata["high"].append(e)

        def _stratum_report(entries: list[dict]) -> dict:
            if not entries:
                return {"n_records": 0, "total": 0}
            # Unique record indices
            rec_indices = set(e["rec_idx"] for e in entries)
            n_records = len(rec_indices)
            avg_identity = sum(self.records[i]["identity"] for i in rec_indices) / n_records if n_records else 0

            total = len(entries)
            n_ref = sum(1 for e in entries if e.get("category") == "ref_support")
            n_smc = sum(1 for e in entries if e.get("category") == "smc_support")
            n_ambig = sum(1 for e in entries if e.get("category") == "ambiguous")

            # Mismatch type distribution
            types: dict[str, int] = {}
            for e in entries:
                key = f"{e['ref']}->{e['smc']}"
                types[key] = types.get(key, 0) + 1

            return {
                "n_records": n_records,
                "identity_range": (min(self.records[i]["identity"] for i in rec_indices),
                                   max(self.records[i]["identity"] for i in rec_indices)),
                "avg_identity": round(avg_identity, 6),
                "total_mismatches": total,
                "ref_support": n_ref,
                "smc_support": n_smc,
                "ambiguous": n_ambig,
                "top_mismatch_types": C(types).most_common(10),
            }

        return {s: _stratum_report(strata[s]) for s in ("low", "mid", "high")}

    # ---- analysis dimension 6: SMC error provenance ----------------------

    def classify_error_provenance(self) -> None:
        """角度 6：SMC 错误溯源。

        对每个 ref_support 的 mismatch（subreads 支持 ref，SMC 错了），判断：
          - artifact_near_gap: 距离任一序列 gap ±_GAP_WINDOW bp → 比对伪影
          - artifact_hp_context: HP >= 5 同聚物区域 → 同聚物相关错误
          - true_error: 远离 gap 且不在高 HP 区域 → 真实共识算法错误

        结果写入每条 ref_support entry 的 'provenance' 和 'hp_info' 字段。
        """
        for e in self._cat_ref:
            rec = self.records[e["rec_idx"]]
            seqs = rec["msa_seqs"]
            smc_s, ref_s = seqs[1], seqs[2]
            pos = e["pos"]

            # Check distance to nearest gap in either sequence
            gaps_ref = self._gap_positions(ref_s)
            gaps_smc = self._gap_positions(smc_s)
            dist_to_gap = self._min_distance_to_gap(pos, gaps_ref | gaps_smc)

            # Check homopolymer context across all sequences
            all_str_seqs = [s for s in seqs if isinstance(s, str)]
            hp_len, hp_base = self._check_hp_at(all_str_seqs, pos)

            # Classify provenance
            if dist_to_gap <= _GAP_WINDOW:
                prov = "artifact_near_gap"
            elif hp_len >= _HP_THRESHOLD:
                prov = "artifact_hp_context"
            else:
                prov = "true_error"

            e["provenance"] = prov
            e["hp_length"] = hp_len
            e["hp_base"] = hp_base
            e["dist_to_gap"] = dist_to_gap

    def provenance_summary(self) -> dict:
        """SMC 错误溯源汇总。"""
        prov_counts: dict[str, int] = C()
        prov_by_hp: dict[str, list[dict]] = {}
        provenance_details: dict[str, list[dict]] = {"true_error": [], "artifact_near_gap": [], "artifact_hp_context": []}

        for e in self._cat_ref:
            prov = e["provenance"]
            prov_counts[prov] += 1
            if len(provenance_details.get(prov, [])) < 50:  # keep up to 50 examples each
                provenance_details[prov].append({
                    "rec_idx": e["rec_idx"],
                    "pos": e["pos"],
                    "identity": e["identity"],
                    "ref": e["ref"],
                    "smc": e["smc"],
                    "pct_ref": e.get("pct_ref", 0),
                    "hp_length": e.get("hp_length", 0),
                    "dist_to_gap": e.get("dist_to_gap", 9999),
                })

        total = len(self._cat_ref)
        return {
            "total_smc_errors": total,
            "provenance_counts": dict(prov_counts),
            "provenance_pct": {k: round(100 * v / total, 1) if total else 0 for k, v in prov_counts.items()},
            "examples": provenance_details,
        }

    # ---- indel helpers ---------------------------------------------------

    def _indel_summary(self) -> dict:
        """Indels 三类统计摘要。"""
        total = len(self._indel_smc) + len(self._indel_ref) + len(self._indel_ambig)
        by_type: dict[str, int] = {}
        for e in self._indel_smc + self._indel_ref + self._indel_ambig:
            t = e.get("indel_type", "unknown")
            by_type[t] = by_type.get(t, 0) + 1

        return {
            "total_indel_positions": total,
            "cat_smc_support": len(self._indel_smc),
            "pct_smc": round(100 * len(self._indel_smc) / total, 1) if total else 0,
            "cat_ref_support": len(self._indel_ref),
            "pct_ref": round(100 * len(self._indel_ref) / total, 1) if total else 0,
            "cat_ambiguous": len(self._indel_ambig),
            "pct_ambiguous": round(100 * len(self._indel_ambig) / total, 1) if total else 0,
            "by_type": by_type,
        }

    def _indel_gap_bias(self) -> dict:
        """分析 SMC gap bias — SMC 在哪些 ref base 上倾向于引入 gap。"""
        ins_wrong = [e for e in self._indel_smc if e.get("indel_type") == "ins_in_smc"]
        if not ins_wrong:
            return {"gap_bias": {}, "total": 0}

        bias: dict[str, dict] = {}
        ref_base_counts = C(e["ref_base"] for e in ins_wrong)
        for base, cnt in ref_base_counts.most_common():
            hits = [hit for hit in ins_wrong if hit["ref_base"] == base]
            strong = sum(1 for hit in hits if hit.get("n_base_in_sbr", 0) >= self.threshold * hit["n_subreads"] / 100)
            bias[base] = {"count": cnt, "strong_support": strong, "weak_support": cnt - strong}
        return {"gap_bias": bias, "total": len(ins_wrong)}

    # ---- transition stats ------------------------------------------------

    def _transition_stats(self) -> dict:
        """Transition vs transversion 统计（仅对 smc_support 组）。"""
        transitions = sum(1 for e in self._cat_smc if (e["ref"] + e["smc"]) in _TRANSITIONS)
        total = len(self._cat_smc)
        return {
            "transitions": transitions,
            "transversions": total - transitions,
            "total": total,
            "transition_pct": round(100 * transitions / total, 1) if total else 0,
        }

    # ---- unified report --------------------------------------------------

    def _report_summary(self) -> dict:
        """Mismatch 全局统计摘要。"""
        total = len(self._cat_ref) + len(self._cat_smc) + len(self._cat_ambig)
        return {
            "total_mismatches": total,
            "cat_ref_support": len(self._cat_ref),
            "pct_ref": round(100 * len(self._cat_ref) / total, 1) if total else 0,
            "cat_smc_support": len(self._cat_smc),
            "pct_smc": round(100 * len(self._cat_smc) / total, 1) if total else 0,
            "cat_ambiguous": len(self._cat_ambig),
            "pct_ambiguous": round(100 * len(self._cat_ambig) / total, 1) if total else 0,
        }

    def _error_matrix_smc(self) -> list[tuple[str, int]]:
        """ref -> SMC 错误矩阵 (仅 smc_support)。"""
        matrix: dict[str, int] = {}
        for e in self._cat_smc:
            key = f"{e['ref']}->{e['smc']}"
            matrix[key] = matrix.get(key, 0) + 1
        return C(matrix).most_common()

    def _bias_analysis(self) -> dict:
        """SMC error bias — SMC 错误调用碱基的分布。"""
        smc_wrong = self._cat_ref
        if not smc_wrong:
            return {"bias": {}, "total": 0}

        calls = C(e["smc"] for e in smc_wrong)
        bias: dict[str, dict] = {}
        for target, cnt in calls.most_common():
            hits = [e for e in smc_wrong if e["smc"] == target]
            from_ref = C(e["ref"] for e in hits)
            avg_pct = sum(e.get("pct_ref", 0) for e in hits) / len(hits)
            bias[target] = {
                "count": cnt,
                "from_ref": dict(from_ref.most_common()),
                "avg_subread_confidence": round(avg_pct, 1),
            }
        return {"bias": bias, "total": len(smc_wrong)}

    def _smc_error_details(self) -> list[dict]:
        """SMC errors (ref_support) 的详细列表，含 provenance。"""
        results: list[dict] = []
        for e in self._cat_ref:
            results.append({
                "rec_idx": e["rec_idx"],
                "pos": e["pos"],
                "ref": e["ref"],
                "smc": e["smc"],
                "pct_ref": e.get("pct_ref", 0),
                "n_subreads": e.get("n_subreads", 0),
                "provenance": e.get("provenance", "unknown"),
                "dist_to_gap": e.get("dist_to_gap", 9999),
                "hp_length": e.get("hp_length", 0),
                "identity": self.records[e["rec_idx"]]["identity"],
            })
        return results

    # ---- full pipeline ---------------------------------------------------

    def run(self) -> None:
        """执行完整分析流水线（6个角度 + provenance）。"""
        if not self.records:
            self._load()

        # Dimensions 1-3
        self.classify_mismatches_and_indels()

        # Dimension 4
        self._smc_support_info = self._analyse_smc_support()

        # Dimension 5
        self._stratified = self._stratify_by_identity()

        # Dimension 6: provenance for SMC errors
        self.classify_error_provenance()
        self._provenance_info = self.provenance_summary()

    @property
    def results(self) -> dict:
        """返回完整分析结果（确保 run() 已执行）。"""
        if not hasattr(self, "_smc_support_info"):
            self.run()

        return {
            "file": self.filepath,
            "identity_thr": self.identity_thr,
            "threshold": self.threshold,
            "n_records": len(self.records),
            "summary": self._report_summary(),
            "error_matrix_smc_to_ref": [list(k) + [v] for k, v in self._error_matrix_smc()],
            "bias_analysis": self._bias_analysis(),
            "transition_stats": self._transition_stats(),
            "indel_summary": self._indel_summary(),
            "indel_bias": self._indel_gap_bias(),
            "smc_support_deep_dive": self._smc_support_info,
            "identity_stratified": self._stratified,
            "error_provenance": self._provenance_info,
        }

    # ---- output formatters -----------------------------------------------

    def json_report(self) -> str:
        """生成 JSON 报告。"""
        r = self.results
        data = {k: v for k, v in r.items() if k != "examples"}
        return json.dumps(data, indent=2, ensure_ascii=False)

    def text_summary(self) -> str:
        """生成简明文本摘要（stdout 用）。"""
        r = self.results
        s = r["summary"]
        L = []

        L.append("=" * 70)
        L.append("ASRTC JSONL 综合分析 — Summary")
        L.append("=" * 70)
        L.append(f"文件:      {self.filepath}")
        L.append(f"参数:      identity_thr={self.identity_thr}, threshold={self.threshold}")
        L.append(f"记录数:    {r['n_records']}")

        L.append(f"\nMismatch 总计:       {s['total_mismatches']}")
        L.append(f"  Subreads → REF:    {s['cat_ref_support']:>6d} ({s['pct_ref']:>5.1f}%)")
        L.append(f"  Subreads → SMC:    {s['cat_smc_support']:>6d} ({s['pct_smc']:>5.1f}%)")
        L.append(f"  Ambiguous:          {s['cat_ambiguous']:>6d} ({s['pct_ambiguous']:>5.1f}%)")

        # Bias
        bias = r["bias_analysis"]
        if bias["total"]:
            L.append(f"\nSMC Error Bias (subreads 支持 ref): {bias['total']} 位")
            for target, info in list(bias["bias"].items())[:8]:
                from_ref_str = ", ".join(f"{k}:{v}" for k, v in info["from_ref"].items())
                L.append(f"  smc='{target}' ({info['count']}x): from={from_ref_str}, avg_conf={info['avg_subread_confidence']}%")

        # Transition stats
        tr = r["transition_stats"]
        if tr["total"]:
            L.append(f"\nTransition / Transversion: {tr['transitions']} / {tr['transversions']} ({tr['transition_pct']}% transition)")

        # Error matrix SMC
        em = r["error_matrix_smc_to_ref"]
        if em:
            L.append("\nError Matrix (ref -> SMC, subreads support SMC):")
            for key in em[:12]:
                tag = " [trans]" if len(key) == 3 and (key[0] + key[2]) in _TRANSITIONS else ""
                L.append(f"  {key[0]}->{key[2]}: {key[3]}{tag}")

        # Indel summary
        indel = r["indel_summary"]
        if indel["total_indel_positions"]:
            L.append(f"\nIndel 总计: {indel['total_indel_positions']}")
            L.append(f"  Subreads → SMC (gap 更可能是伪影): {indel['cat_smc_support']} ({indel['pct_smc']}%)")
            L.append(f"  Subreads → ref   (gap 更可能是真实缺失): {indel['cat_ref_support']} ({indel['pct_ref']}%)")
            L.append(f"  Ambiguous:       {indel['cat_ambiguous']} ({indel['pct_ambiguous']}%)")
            L.append(f"  按类型: {indel['by_type']}")

        # SMC support deep dive
        ss = r["smc_support_deep_dive"]
        if ss["total"]:
            L.append(f"\nsmc_support 候选 (潜在生物学差异):")
            L.append(f"  总计位点数: {ss['total']}, 涉及 records: {ss['by_record_count']}")
            L.append("\n  Top-5 records:")
            for rec_info in ss["top_records"][:5]:
                t_tag = f"T={rec_info['transitions']}, tv={rec_info['transversions']}"
                L.append(
                    f"    Rec[{rec_info['rec_idx']}] id={rec_info['identity']:.6f} "
                    f"mismatches={rec_info['n_mismatches']} avg_conf={rec_info['avg_subread_confidence']}% {t_tag}"
                )

        # Identity stratified
        st = r["identity_stratified"]
        L.append("\nIdentity 分层:")
        for name in ("low", "mid", "high"):
            info = st[name]
            if info["n_records"]:
                lo, hi = info["identity_range"]
                L.append(f"\n  [{name.upper()}] {lo:.4f}-{hi:.4f}: records={info['n_records']} avg_id={info['avg_identity']:.6f} mm={info['total_mismatches']}")
                L.append(f"    ref_sup={info['ref_support']} smc_sup={info['smc_support']} ambig={info['ambiguous']}")

        # Error provenance
        ep = r["error_provenance"]
        if ep["total_smc_errors"]:
            L.append("\nSMC 错误溯源 (provenance):")
            for prov, pct in ep["provenance_pct"].items():
                L.append(f"  {prov}: {ep['provenance_counts'].get(prov, 0):>5d} ({pct}%)")

        # Gap bias
        gb = r["indel_bias"]
        if gb["total"]:
            L.append(f"\nSMC Gap Bias (SMC 错误引入 gap): {gb['total']} 次")
            for base, info in list(gb["gap_bias"].items())[:5]:
                L.append(f"  ref='{base}': count={info['count']}, strong={info['strong_support']}, weak={info['weak_support']}")

        return "\n".join(L)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="ASRTC JSONL 综合分析 — SMC vs reference 差异全维度分析",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", help="ASRTC JSONL 文件路径 (.asrtc.txt)")
    parser.add_argument("--output", "-o", choices=["text", "json"], default="text",
                        help="输出格式 (default: text)")
    parser.add_argument("--out-file", metavar="PATH", default=None,
                        help="输出文件路径 (default: stdout)")
    parser.add_argument("--identity-thr", type=float, default=1.0,
                        help="identity 阈值，>=此值的记录被排除 (default: 1.0)")
    parser.add_argument("--threshold", type=int, default=75,
                        help="subreads 强支持判定阈值 % (default: 75)")
    return parser


def main(args=None) -> None:
    parser = build_parser()
    opts = parser.parse_args(args)

    analyzer = AsrtcAnalyzer(opts.input, identity_thr=opts.identity_thr, threshold=opts.threshold)
    analyzer.run()

    if opts.output == "json":
        report = analyzer.json_report()
    else:
        report = analyzer.text_summary()

    if opts.out_file:
        with open(opts.out_file, "w") as fh:
            fh.write(report + "\n")
        print(f"Report written to {opts.out_file}")
    else:
        print(report)


if __name__ == "__main__":
    main()
