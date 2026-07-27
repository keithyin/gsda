"""通用 ASRTC JSONL 分析工具 — SMC consensus vs reference 差异分析。

用法:
    python asrtc_analyzer.py input.asrtc.txt
    python asrtc_analyzer.py input.asrtc.txt --output json --out-file result.json
    python asrtc_analyzer.py input.asrtc.txt --threshold 90 --top-k 3

分析每个 record 的 msa_seqs 字段:
    - 第0个: qual (忽略)
    - 第1个: smc_seq (subreads consensus)
    - 第2个: reference_seq
    - 第3+: subreads seq
关注 ref != SMC 的位置，用 subreads 投票判断谁更可信。
"""
import json
import sys
import argparse
from collections import Counter as C

# ---------------------------------------------------------------------------
# Core analysis class
# ---------------------------------------------------------------------------

class ASRTCAnalyzer:
    """分析 ASRTC JSONL 文件中 SMC consensus 与 reference 的差异。"""

    _TRANSITIONS = frozenset(['GA', 'AG', 'CT', 'TC'])

    def __init__(self, filepath, identity_thr=1.0, threshold=75):
        self.filepath = filepath
        self.identity_thr = identity_thr
        self.threshold = threshold          # strong consensus %
        self.records = []
        self._load()

    # ---- loading ---------------------------------------------------------

    def _load(self):
        """加载 JSONL 文件，过滤 identity >= threshold 的记录。"""
        self.records = []
        with open(self.filepath, 'r') as fh:
            for line in fh:
                rec = json.loads(line)
                if rec.get('identity', 1.0) >= self.identity_thr:
                    continue          # skip "perfect" by default
                self.records.append(rec)

    # ---- classification --------------------------------------------------

    def classify_mismatches(self):
        """遍历所有记录，将 mismatch 分为三类。

        Returns
        -------
        tuple[list, list, list]
            (cat_ref_support, cat_smc_support, cat_ambig)
        """
        cat_ref = []
        cat_smc = []
        cat_ambig = []

        for rec_idx, rec in enumerate(self.records):
            seqs = rec['msa_seqs']
            if len(seqs) < 3:
                continue
            smc_s = seqs[1]
            ref_s = seqs[2]
            sbr_seq = seqs[3:]

            for pos in range(len(ref_s)):
                rc = ref_s[pos]
                sc = smc_s[pos]
                if rc == '-' or sc == '-' or rc == sc:
                    continue

                non_gap = []
                for sr in sbr_seq:
                    ch = sr[pos] if pos < len(sr) else None
                    if ch is not None and ch != '-':
                        non_gap.append(ch)
                if not non_gap:
                    continue

                n_ref_c = sum(1 for c in non_gap if c == rc)
                n_smc_c = sum(1 for c in non_gap if c == sc)
                pct_ref = n_ref_c / len(non_gap) * 100
                pct_smc_v = n_smc_c / len(non_gap) * 100

                entry = {
                    'rec_idx': rec_idx,
                    'pos': pos,
                    'ref': rc,
                    'smc': sc,
                    'identity': rec.get('identity', 0),
                    'query_aln_len': rec.get('query_aln_len', 0),
                }

                if pct_ref >= self.threshold and pct_smc_v < (100 - self.threshold):
                    entry['sub_counts'] = C(non_gap).most_common()
                    entry['pct_ref'] = pct_ref
                    entry['n_subreads'] = len(sbr_seq)
                    cat_ref.append(entry)
                elif pct_smc_v >= self.threshold and pct_ref < (100 - self.threshold):
                    entry['sub_counts'] = C(non_gap).most_common()
                    entry['pct_smc'] = pct_smc_v
                    entry['n_subreads'] = len(sbr_seq)
                    cat_smc.append(entry)
                else:
                    entry['sub_counts'] = C(non_gap).most_common()
                    entry['pct_ref'] = pct_ref
                    entry['pct_smc'] = pct_smc_v
                    entry['n_subreads'] = len(sbr_seq)
                    cat_ambig.append(entry)

        self._cat_ref = cat_ref
        self._cat_smc = cat_smc
        self._cat_ambig = cat_ambig
        return cat_ref, cat_smc, cat_ambig

    # ---- indel reporting helper (shared by text/markdown formatters) ------

    def _add_indel_sections(self, lines, r):
        """为 text 格式追加 indel 分析章节。"""
        indel = r['indel_summary']
        bias = r['indel_bias']
        gap_ref = r['indel_gap_in_ref']

        lines.append('\n' + '=' * 70)
        lines.append('INDEL ANALYSIS (gap 位点分类)')
        lines.append('=' * 70)
        total_indel = indel['total_indel_positions']
        lines.append(f'\nIndel 总计: {total_indel}')
        sub = lambda hdr_text, val: lines.append(f'{hdr_text}: {val}')
        sub('Subreads 支持 SMC (gap 可能错误)', f'{indel["cat_smc_support"]} ({indel["pct_smc"]:.1f}%)')
        sub('Subreads 支持 ref (gap 可能真实)',  f'{indel["cat_ref_support"]} ({indel["pct_ref"]:.1f}%)')
        sub('Ambiguous',                           f'{indel["cat_ambig"]} ({indel["pct_ambig"]:.1f}%)')

        lines.append(f'\n按类型分布: {indel["by_type"]}')

        # --- SMC gap bias (SMC 不该有 gap 但引入了) -----------------
        lines.append('\n--- SMC gap bias: SMC 错误引入 gap (subreads 支持 ref base) ---')
        if bias['total']:
            for base, info in bias['gap_bias'].items():
                strong = info['strong_support']
                weak = info['weak_support']
                lines.append(
                    f"  ref='{base}' -> SMC gap: {info['count']} 次"
                    f" (subreads强支持{strong}, 弱支持{weak})"
                )
        else:
            lines.append('  无')

        # --- del_in_ref ---
        lines.append('\n--- del_in_ref: ref gap 但 SMC/Subreads 有 base ---')
        if gap_ref['total']:
            lines.append(f"总计 {gap_ref['total']} 位")
            lines.append(f'SMC/子序列在此处 CALL的碱基分布: {gap_ref["by_smc_base"]}')
            # Show examples
            shown = 0
            for e in self._indel_smc:
                if e['indel_type'] != 'del_in_ref':
                    continue
                if shown >= 3:
                    break
                rec = self.records[e['rec_idx']]
                seqs = rec['msa_seqs']
                smc_s, ref_s = seqs[1], seqs[2]
                sbr_seq = seqs[3:]
                pos = e['pos']
                ctx_start = max(0, pos - 6)
                ctx_end = min(len(ref_s), pos + 7)
                rctx = ref_s[ctx_start:pos] + '[' + ref_s[pos] + ']' + ref_s[pos+1:ctx_end]
                sctx = smc_s[ctx_start:pos] + '[' + smc_s[pos] + ']' + smc_s[pos+1:ctx_end]

                sub_str = f"{e['n_base_in_sbr']}/{e['n_subreads']} subreads 有 base"
                lines.append(f"\n  Rec[{e['rec_idx']}] id={rec['identity']:.6f} pos={pos}: ref=gap, smc='{e['smc_base']}'")
                lines.append(f"    ref: {rctx}")
                lines.append(f"    smc: {sctx}")
                lines.append(f"    subreads: {sub_str}")

                for i, sr in enumerate(sbr_seq[:6]):
                    ch = sr[pos] if pos < len(sr) else '?'
                    lines.append(f"      [{i}] '{ch}'")
                shown += 1
        else:
            lines.append('  无')

        # --- indel ambiguous ---
        lines.append('\n--- Indel Ambiguous (subreads 不强烈支持任一) ---')
        if self._indel_ambig:
            lines.append(f'数量: {len(self._indel_ambig)}')
            shown = 0
            for e in self._indel_ambig[:3]:
                rec = self.records[e['rec_idx']]
                seqs = rec['msa_seqs']
                smc_s, ref_s = seqs[1], seqs[2]
                sbr_seq = seqs[3:]
                pos = e['pos']

                t_label = f"{'ref' if e['indel_type']=='del_in_ref' else 'smc'} gap"
                lines.append(f"\n  Rec[{e['rec_idx']}] id={rec['identity']:.6f} pos={pos}: {t_label}")
                pct = e.get('pct_sub_has_base', 0)
                lines.append(f"    subreads base ratio: {pct:.1f}% ({e['n_base_in_sbr']}/{e['n_subreads']})")
                shown += 1

        return lines

    # ---- indel classification --------------------------------------------

    def classify_indels(self):
        """分类 gap 位点：判断 ref 或 SMC 引入的 gap 是否合理。

        Returns
        -------
        tuple[list, list, list]
            (cat_smc_support, cat_ref_support, cat_ambig)
        三种 indel 类型:
          - del_in_ref:   ref=gap / smc=base
          - ins_in_smc:   smc=gap / ref=base
          - gap_both:     ref=gap / smc=gap (无分析价值，跳过)

        subreads 投票规则与 base mismatch 一致:
            strong >= threshold, weak < 100-threshold, 其余 ambiguous.
        """
        cat_smc = []   # subreads 支持 SMC (gap 更可能是算法错误引入)
        cat_ref = []   # subreads 支持 ref (gap 更可能是真实缺失)
        cat_ambig = []

        for rec_idx, rec in enumerate(self.records):
            seqs = rec['msa_seqs']
            if len(seqs) < 3:
                continue
            smc_s = seqs[1]
            ref_s = seqs[2]
            sbr_seq = seqs[3:]

            for pos in range(len(ref_s)):
                rc = ref_s[pos]
                sc = smc_s[pos]

                is_gap_ref = (rc == '-')
                is_gap_smc = (sc == '-')

                if not is_gap_ref and not is_gap_smc:
                    continue  # non-gap columns already handled by classify_mismatches

                # Collect subread bases at this position
                sub_bases = []
                for sr in sbr_seq:
                    ch = sr[pos] if pos < len(sr) else None
                    if ch is not None and ch != '-':
                        sub_bases.append(ch)
                n_base = len(sub_bases)

                entry = {
                    'rec_idx': rec_idx,
                    'pos': pos,
                    'n_subreads': len(sbr_seq),
                    'identity': rec.get('identity', 0),
                    'sub_counts': C(sub_bases).most_common() if sub_bases else [],
                    'n_base_in_sbr': n_base,
                }

                pct = n_base / len(sbr_seq) * 100 if sbr_seq else 0
                entry['pct_sub_has_base'] = pct

                if is_gap_ref and not is_gap_smc:
                    # del_in_ref: ref gap, SMC calls base
                    entry['indel_type'] = 'del_in_ref'
                    entry['smc_base'] = sc
                    entry['ref_ctx'] = 'gap'

                    strong = (pct >= self.threshold)
                    weak   = (pct < (100 - self.threshold))

                    if strong and weak:
                        cat_smc.append(entry)  # subreads have base → SMC correct
                    elif not strong and weak:
                        cat_ref.append(entry)  # few subreads have base → ref gap likely real
                    else:
                        cat_ambig.append(entry)

                elif is_gap_smc and not is_gap_ref:
                    # ins_in_smc: SMC gap, ref has base
                    entry['indel_type'] = 'ins_in_smc'
                    entry['ref_base'] = rc
                    entry['smc_ctx'] = 'gap'

                    strong = (pct >= self.threshold)
                    weak   = (pct < (100 - self.threshold))

                    if strong and weak:
                        cat_smc.append(entry)  # subreads have base → SMC gap is wrong
                    elif not strong and weak:
                        cat_ref.append(entry)  # few subreads have base → SMC gap likely real
                    else:
                        cat_ambig.append(entry)

                # gap_both: skip (alignment artifact, no informative signal)

        self._indel_smc = cat_smc
        self._indel_ref = cat_ref
        self._indel_ambig = cat_ambig
        return cat_smc, cat_ref, cat_ambig

    # ---- indel helpers ---------------------------------------------------

    def indel_summary(self):
        """indels 三类统计摘要。"""
        total = len(self._indel_smc) + len(self._indel_ref) + len(self._indel_ambig)
        by_type = {}
        for e in self._indel_smc + self._indel_ref + self._indel_ambig:
            t = e['indel_type']
            by_type[t] = by_type.get(t, 0) + 1

        return {
            'total_indel_positions': total,
            'cat_smc_support': len(self._indel_smc),
            'cat_ref_support': len(self._indel_ref),
            'cat_ambig': len(self._indel_ambig),
            'pct_smc': 100 * len(self._indel_smc) / total if total else 0,
            'pct_ref': 100 * len(self._indel_ref) / total if total else 0,
            'pct_ambig': 100 * len(self._indel_ambig) / total if total else 0,
            'by_type': by_type,
        }

    def indel_bias_analysis(self):
        """分析 SMC gap bias — SMC 在哪些 ref base 上倾向于引入 gap。"""
        # Only ins_in_smc where subreads disagree with ref (SMC gap is likely wrong)
        ins_wrong = [e for e in self._indel_smc if e['indel_type'] == 'ins_in_smc']
        bias = {}
        if not ins_wrong:
            return {'gap_bias': {}, 'total': 0}

        # Which ref bases does SMC most often turn into gaps?
        ref_bases = C(e['ref_base'] for e in ins_wrong)
        for base, cnt in ref_bases.most_common():
            hits = [e for e in ins_wrong if e['ref_base'] == base]
            # Among these, what fraction have strong subread support?
            strong = sum(1 for e in hits if e.get('pct_sub_has_base', 0) >= self.threshold)
            bias[base] = {
                'count': cnt,
                'strong_support': strong,
                'weak_support': cnt - strong,
            }
        return {'gap_bias': bias, 'total': len(ins_wrong)}

    def indel_gap_in_ref_analysis(self):
        """分析 del_in_ref — ref 有 gap 但 SMC/Subreads 有 base 的情况。"""
        # Where ref has gap but subreads have base (SMC likely correct)
        entries = [e for e in self._indel_smc if e['indel_type'] == 'del_in_ref']
        return {
            'total': len(entries),
            'by_smc_base': dict(C(e['smc_base'] for e in entries).most_common()),
        }

    # ---- results aggregation -----------------------------------------------

    @property
    def results(self):
        """Ensure classify_mismatches and classify_indels have been called."""
        if not hasattr(self, '_cat_ref'):
            self.classify_mismatches()
        if not hasattr(self, '_indel_smc'):
            self.classify_indels()
        return {
            'summary': self.report_summary(),
            'bias': self.bias_analysis(),
            'transition_stats': self.transition_stats(),
            'error_matrix_smc': self.error_matrix('smc'),
            'error_matrix_ref': self.error_matrix('ref'),
            'indel_summary': self.indel_summary(),
            'indel_bias': self.indel_bias_analysis(),
            'indel_gap_in_ref': self.indel_gap_in_ref_analysis(),
        }

    def to_dict(self):
        """完整结果序列化为 dict。"""
        r = self.results
        return {
            'file': self.filepath,
            'identity_thr': self.identity_thr,
            'threshold': self.threshold,
            'n_records': len(self.records),
            'summary': r['summary'],
            'bias_analysis': r['bias'],
            'transition_stats': r['transition_stats'],
            'error_matrix_smc_to_ref': [list(k) + [v] for k, v in r['error_matrix_smc']],
            'error_matrix_ref_to_smc': [list(k) + [v] for k, v in r['error_matrix_ref']],
            'indel_summary': r['indel_summary'],
            'indel_bias': r['indel_bias'],
            'indel_gap_in_ref': r['indel_gap_in_ref'],
        }

    # ---- summary ---------------------------------------------------------

    def report_summary(self):
        """返回三类统计摘要。"""
        total = len(self._cat_ref) + len(self._cat_smc) + len(self._cat_ambig)
        return {
            'total_mismatches': total,
            'cat_ref_support': len(self._cat_ref),
            'cat_smc_support': len(self._cat_smc),
            'cat_ambig': len(self._cat_ambig),
            'pct_ref': 100 * len(self._cat_ref) / total if total else 0,
            'pct_smc': 100 * len(self._cat_smc) / total if total else 0,
            'pct_ambig': 100 * len(self._cat_ambig) / total if total else 0,
        }

    def error_matrix(self, category='smc'):
        """ref -> SMC 错误矩阵 (只统计 cat_smc 或 cat_ref)。"""
        cats = {'smc': self._cat_smc, 'ref': self._cat_ref}[category]
        matrix = {}
        for entry in cats:
            key = (entry['ref'], entry['smc'])
            matrix[key] = matrix.get(key, 0) + 1
        return C(matrix).most_common()

    def bias_analysis(self):
        """分析 SMC consensus 的碱基偏倚。"""
        # SMC wrong: subreads support ref
        smc_wrong = self._cat_ref
        if not smc_wrong:
            return {'bias': {}, 'total': 0}

        # What does the SMC call most when it is wrong?
        calls = C(e['smc'] for e in smc_wrong)
        bias = {}
        for target, cnt in calls.most_common():
            hits = [e for e in smc_wrong if e['smc'] == target]
            from_ref = C(e['ref'] for e in hits)
            avg_pct = sum(e['pct_ref'] for e in hits) / len(hits)
            bias[target] = {
                'count': cnt,
                'from_ref': dict(from_ref.most_common()),
                'avg_subread_confidence': round(avg_pct, 1),
            }
        return {'bias': bias, 'total': len(smc_wrong)}

    def transition_stats(self):
        """Transition vs transversion 统计。"""
        transitions = sum(1 for e in self._cat_smc if (e['ref'] + e['smc']) in self._TRANSITIONS)
        transversions = len(self._cat_smc) - transitions
        total = len(self._cat_smc)
        return {
            'transitions': transitions,
            'transversions': transversions,
            'total': total,
            'transition_pct': round(100 * transitions / total, 1) if total else 0,
        }

    # ---- detailed example lookup -----------------------------------------

    def _get_context(self, entry):
        """根据 entry dict 获取 ref/smc context 和 subreads 详情。"""
        rec = self.records[entry['rec_idx']]
        seqs = rec['msa_seqs']
        smc_s, ref_s = seqs[1], seqs[2]
        sbr_seq = seqs[3:]
        pos = entry['pos']

        cs = max(0, pos - 8)
        ce = min(len(ref_s), pos + 9)
        rctx = ref_s[cs:pos] + '[' + ref_s[pos] + ']' + ref_s[pos+1:ce]
        sctx = smc_s[cs:pos] + '[' + smc_s[pos] + ']' + smc_s[pos+1:ce]

        sub_details = []
        for i, sr in enumerate(sbr_seq[:8]):
            ch = sr[pos] if pos < len(sr) else '?'
            strand = rec['strands'][i + 2] if i + 2 < len(rec.get('strands', [])) else '?'
            sname = rec['names'][i + 2].split('/')[-1][:40] if i + 2 < len(rec['names']) else '?'
            sub_details.append({'strand': strand, 'name': sname, 'call': ch})

        nc = entry.get('sub_counts', [])
        return {
            'rctx': rctx,
            'sctx': sctx,
            'sub_reads': sub_details,
            'top2': (nc[0][0], nc[0][1], nc[1][0], nc[1][1]) if len(nc) >= 2 else (nc[0][0], nc[0][1], None, None),
        }

    # ---- report formatters -----------------------------------------------

    def text_report(self, top_k=5):
        """生成纯文本报告。"""
        r = self.results
        s = r['summary']
        lines = []

        def hdr(title):
            lines.append('\n' + '=' * 70)
            lines.append(title)
            lines.append('=' * 70)

        def sub(hdr_text, val):
            lines.append(f'{hdr_text}: {val}')

        # ---- summary ---------------------------------------------------
        lines.append('=' * 70)
        lines.append('SUMMARY')
        lines.append('=' * 70)
        lines.append(f'文件:    {self.filepath}')
        lines.append(f'identity_thr={self.identity_thr} threshold={self.threshold}')
        lines.append(f'记录数:  {len(self.records)}')
        lines.append(f'\nMismatch 总计: {s["total_mismatches"]}')
        sub('Subreads support REFERENCE', f'{s["cat_ref_support"]} ({s["pct_ref"]:.1f}%)')
        sub('Subreads support SMC',      f'{s["cat_smc_support"]} ({s["pct_smc"]:.1f}%)')
        sub('Ambiguous (no clear)',     f'{s["cat_ambig"]} ({s["pct_ambig"]:.1f}%)')

        # ---- bias ------------------------------------------------------
        hdr('CATEGORY A: SMC ERRORS (subreads support reference)')
        bias = r['bias']
        if bias['total']:
            lines.append(f'\nSMC 错误 CALL碱基分布: {len(bias["bias"])} 种类型, 共 {bias["total"]} 位')
            for target, info in bias['bias'].items():
                lines.append(f"\n  SMC  CALL '{target}' ({info['count']} 次):")
                for ref_base, cnt in info['from_ref'].items():
                    lines.append(
                        f"    ref='{ref_base}' -> smc='{target}': {cnt} 次,"
                        f" subread 平均置信度 {info['avg_subread_confidence']}%"
                    )
        else:
            lines.append('\n无 (SMC 与 reference 在所有位置一致)')

        # ---- error matrix SMC ------------------------------------------
        hdr('CATEGORY B: POTENTIAL BIOLOGICAL DIFFERENCES (subreads support SMC)')
        sub('SMC  CALL碱基分布', str(C(e['smc'] for e in self._cat_smc)))

        lines.append('\n错误矩阵 (ref base -> what SMC calls):')
        refs = sorted(set(e['ref'] for e in self._cat_smc))
        for ref_base in refs:
            hits = [e for e in self._cat_smc if e['ref'] == ref_base]
            smc_for_ref = C(e['smc'] for e in hits)
            avg_pct = sum(e.get('pct_smc', 0) for e in hits) / len(hits)
            lines.append(f'  ref=\'{ref_base}\' (total {len(hits)} positions):')
            for sc, cnt in smc_for_ref.most_common():
                tag = ' [trans]' if (ref_base + sc) in self._TRANSITIONS else ''
                lines.append(f"    -> '{sc}': {cnt}{tag}")
            lines.append(f'    Avg subread confidence: {avg_pct:.1f}%')

        tr = r['transition_stats']
        lines.append(f"\nTransition errors: {tr['transitions']} ({tr['transition_pct']}%)")
        lines.append(f"Transversion errors: {tr['transversions']} ({100 - tr['transition_pct']}%)")

        # ---- indel analysis ----------------------------------------------
        self._add_indel_sections(lines, r)

        # ---- ambiguous -------------------------------------------------
        hdr('CATEGORY C: AMBIGUOUS (subreads don\'t strongly support either)')
        if self._cat_ambig:
            rec_ids = set(e['rec_idx'] for e in self._cat_ambig)
            avg_id = sum(self.records[i]['identity'] for i in rec_ids) / len(rec_ids) if rec_ids else 0
            lines.append(f'\n数量: {len(self._cat_ambig)}')
            lines.append(f'受影响记录平均 identity: {avg_id:.6f}')
            lines.append('\nSample ambiguous positions:')
            shown = 0
            for e in self._cat_ambig:
                if shown >= top_k:
                    break
                rec = self.records[e['rec_idx']]
                ctx = self._get_context(e)
                t2 = ctx['top2']
                lines.append(f"\n  Rec[{e['rec_idx']}] id={rec['identity']:.6f} pos={e['pos']} ref={e['ref']} smc={e['smc']}")
                lines.append(f"    ref: {ctx['rctx']}")
                lines.append(f"    smc: {ctx['sctx']}")
                n = t2[1]
                denom = str(n) + '/' + str(sum(v for v in t2[1:] if isinstance(v, int)))
                lines.append(f"    subreads: {t2[0]}({n}/{denom}) / {t2[2]}({t2[3] if t2[3] is not None else '?'})")
                for sb in ctx['sub_reads'][:6]:
                    lines.append(f"      [{sb['strand']}] {sb['name']}: '{sb['call']}'")
                shown += 1

        return '\n'.join(lines)

    def json_report(self):
        """生成 JSON 输出 (可直接保存为文件)。"""
        import collections

        data = self.to_dict()
        # add examples
        examples = {'cat_ref': [], 'cat_smc': [], 'cat_ambig': []}
        for e in self._cat_ref[:10]:
            ctx = self._get_context(e)
            examples['cat_ref'].append({**e, **{'context': {'rctx': ctx['rctx'], 'sctx': ctx['sctx'], 'subreads': ctx['sub_reads']}}})
        for e in self._cat_smc[:10]:
            ctx = self._get_context(e)
            examples['cat_smc'].append({**e, **{'context': {'rctx': ctx['rctx'], 'sctx': ctx['sctx'], 'subreads': ctx['sub_reads']}}})
        for e in self._cat_ambig[:10]:
            ctx = self._get_context(e)
            examples['cat_ambig'].append({**e, **{'context': {'rctx': ctx['rctx'], 'sctx': ctx['sctx'], 'subreads': ctx['sub_reads']}}})

        # indel examples
        examples['indel_smc_support'] = []
        for e in self._indel_smc[:10]:
            rec = self.records[e['rec_idx']]
            seqs = rec['msa_seqs']
            smc_s, ref_s = seqs[1], seqs[2]
            sbr_seq = seqs[3:]
            pos = e['pos']
            cs = max(0, pos - 6)
            ce = min(len(ref_s), pos + 7)
            rctx = ref_s[cs:pos] + '[' + ref_s[pos] + ']' + ref_s[pos+1:ce]
            sctx = smc_s[cs:pos] + '[' + smc_s[pos] + ']' + smc_s[pos+1:ce]
            sub_details = []
            for i, sr in enumerate(sbr_seq[:8]):
                ch = sr[pos] if pos < len(sr) else '?'
                strand = rec['strands'][i+2] if i+2 < len(rec.get('strands',[])) else '?'
                sname = rec['names'][i+2].split('/')[-1][:30] if i+2 < len(rec['names']) else '?'
                sub_details.append({'strand': strand, 'call': ch, 'name': sname})
            examples['indel_smc_support'].append({
                **e, 'context': {'rctx': rctx, 'sctx': sctx}, 'subreads': sub_details,
            })
        examples['indel_ref_support'] = [
            {**e} for e in self._indel_ref[:10]
        ]
        examples['indel_ambig'] = [
            {**e} for e in self._indel_ambig[:10]
        ]

        data['examples'] = examples
        return json.dumps(data, indent=2, ensure_ascii=False)

    def markdown_report(self, top_k=5):
        """生成 Markdown 报告。"""
        r = self.results
        s = r['summary']
        L = []

        def h1(title):
            L.append(f'\n## {title}\n')

        def h2(title):
            L.append(f'\n### {title}\n')

        def tbl(headers, rows):
            L.append('| ' + ' | '.join(headers) + ' |')
            L.append('|' + '|'.join(['---' for _ in headers]) + '|')
            for row in rows:
                L.append('| ' + ' | '.join(str(x) for x in row) + ' |')
            L.append('')

        # ---- header ----------------------------------------------------
        L.append(f'# ASRTC Mismatch Analysis Report')
        L.append(f'\n**文件**: {self.filepath}\n\n')
        L.append(f'**参数**: identity_thr={self.identity_thr}, threshold={self.threshold}\n\n')
        L.append(f'**记录数**: {len(self.records)}\n')

        # ---- summary ---------------------------------------------------
        h1('Summary')
        tbl(
            ['类别', '数量', '占比'],
            [
                ['Subreads support REFERENCE', s['cat_ref_support'], f'{s["pct_ref"]:.1f}%'],
                ['Subreads support SMC',      s['cat_smc_support'], f'{s["pct_smc"]:.1f}%'],
                ['Ambiguous (no clear)',     s['cat_ambig'],       f'{s["pct_ambig"]:.1f}%'],
                ['Mismatch 总计',            s['total_mismatches'], '100%'],
            ],
        )

        # ---- bias ------------------------------------------------------
        h2('SMC Error Bias (Category A)')
        bias = r['bias']
        if bias['total']:
            tbl(
                ['SMC CALL', '次数', '来源分布', '平均置信度'],
                [
                    (target, info['count'], str(info['from_ref']), f"{info['avg_subread_confidence']}%")
                    for target, info in bias['bias'].items()
                ],
            )

        # ---- error matrix SMC ------------------------------------------
        h2('Error Matrix — SMC Supported (Category B)')
        rows = []
        refs = sorted(set(e['ref'] for e in self._cat_smc))
        for ref_base in refs:
            hits = [e for e in self._cat_smc if e['ref'] == ref_base]
            smc_for_ref = C(e['smc'] for e in hits)
            for sc, cnt in smc_for_ref.most_common():
                tag = 'T' if (ref_base + sc) in self._TRANSITIONS else 'tv'
                rows.append((ref_base, sc, cnt, tag))
        tbl(['ref', 'SMC calls', 'count', 'type'], rows)

        tr = r['transition_stats']
        L.append(f'**Transition / Transversion**: {tr["transitions"]} / {tr["transversions"]} ({tr["transition_pct"]}% transition)\n')

        # ---- indel analysis ----------------------------------------------
        h2('Indel Analysis (gap 位点分类)')
        indel = r['indel_summary']
        tbl(
            ['类别', '数量', '占比'],
            [
                ['Subreads 支持 SMC',      indel['cat_smc_support'], f'{indel["pct_smc"]:.1f}%'],
                ['Subreads 支持 ref',       indel['cat_ref_support'], f'{indel["pct_ref"]:.1f}%'],
                ['Ambiguous',              indel['cat_ambig'],       f'{indel["pct_ambig"]:.1f}%'],
                ['Indel 总计',             indel['total_indel_positions'], '100%'],
            ],
        )
        L.append(f'**按类型**: {json.dumps(indel["by_type"])}\n')

        # SMC gap bias
        ib = r['indel_bias']
        if ib['total']:
            h2('SMC Gap Bias (subreads 支持 ref base)')
            tbl(
                ['ref', 'SMC引入gap次数', '强支持', '弱支持'],
                [(base, info['count'], info['strong_support'], info['weak_support'])
                 for base, info in ib['gap_bias'].items()],
            )

        # del_in_ref
        gr = r['indel_gap_in_ref']
        if gr['total']:
            L.append(f'**del_in_ref (ref gap 但 SMC/Subreads 有 base)**: {gr["total"]} 位\n')
            L.append(f'SMC/子序列在此处 CALL的碱基分布: `{json.dumps(gr["by_smc_base"])}`\n')
            shown = 0
            for e in self._indel_smc:
                if e['indel_type'] != 'del_in_ref':
                    continue
                if shown >= 3:
                    break
                rec = self.records[e['rec_idx']]
                seqs = rec['msa_seqs']
                smc_s, ref_s = seqs[1], seqs[2]
                pos = e['pos']
                ctx_start = max(0, pos - 6)
                ctx_end = min(len(ref_s), pos + 7)
                rctx = ref_s[ctx_start:pos] + '[' + ref_s[pos] + ']' + ref_s[pos+1:ctx_end]
                sctx = smc_s[ctx_start:pos] + '[' + smc_s[pos] + ']' + smc_s[pos+1:ctx_end]
                sub_str = f"{e['n_base_in_sbr']}/{e['n_subreads']} subreads 有 base"
                L.append(f'- Rec[{e["rec_idx"]}] id={rec["identity"]:.6f} pos={pos}: ref=gap, smc=`{e["smc_base"]}`')
                L.append(f'  - ref: `{rctx}` | smc: `{sctx}` | subreads: {sub_str}\n')
                shown += 1

        if self._indel_ambig:
            h2('Indel Ambiguous (subreads 不强烈支持任一)')
            L.append(f'数量: {len(self._indel_ambig)}\nSample:\n')
            for e in self._indel_ambig[:3]:
                rec = self.records[e['rec_idx']]
                t_label = 'ref gap' if e['indel_type'] == 'del_in_ref' else 'smc gap'
                pct = e.get('pct_sub_has_base', 0)
                L.append(f'- Rec[{e["rec_idx"]}] id={rec["identity"]:.6f} pos={e["pos"]}: {t_label}, subreads base ratio: {pct:.1f}% ({e["n_base_in_sbr"]}/{e["n_subreads"]})\n')

        # ---- ambiguous -------------------------------------------------
        h2('Ambiguous Positions (Category C)')
        if self._cat_ambig:
            rec_ids = set(e['rec_idx'] for e in self._cat_ambig)
            avg_id = sum(self.records[i]['identity'] for i in rec_ids) / len(rec_ids) if rec_ids else 0
            L.append(f'数量: {len(self._cat_ambig)} | 平均 identity: {avg_id:.6f}\n')
            L.append('Sample:')
            shown = 0
            for e in self._cat_ambig:
                if shown >= top_k:
                    break
                rec = self.records[e['rec_idx']]
                ctx = self._get_context(e)
                t2 = ctx['top2']
                L.append(f'  - Rec[{e["rec_idx"]}] pos={e["pos"]} ref={e["ref"]} smc={e["smc"]} id={rec["identity"]:.6f} sub={t2[0]}({t2[1]})/{t2[2]}({t2[3] if t2[3] is not None else "?"})')
                shown += 1

        return '\n'.join(L)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _jsonl_iter(fh):
    """迭代 JSONL 文件中的每条记录。"""
    for line in fh:
        line = line.strip()
        if line:
            yield json.loads(line)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description='通用 ASRTC JSONL 分析工具 — SMC consensus vs reference 差异分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('input', help='ASRTC JSONL 文件路径 (.asrtc.txt)')
    parser.add_argument('--output', '-o', choices=['text', 'json', 'markdown'], default='text',
                        help='输出格式 (default: text)')
    parser.add_argument('--out-file', metavar='PATH', default=None,
                        help='输出文件路径 (default: stdout)')
    parser.add_argument('--identity-thr', type=float, default=1.0,
                        help='identity 阈值，>=此值的记录被排除 (default: 1.0)')
    parser.add_argument('--threshold', type=float, default=75,
                        help='subread consensus 强支持阈值 (default: 75)')
    parser.add_argument('--top-k', type=int, default=5,
                        help='每个 category 显示 sample 数 (default: 5)')
    return parser


def main(args=None):
    parser = build_parser()
    opts = parser.parse_args(args)

    analyzer = ASRTCAnalyzer(opts.input, identity_thr=opts.identity_thr, threshold=opts.threshold)
    analyzer.classify_mismatches()

    if opts.output == 'text':
        report = analyzer.text_report(top_k=opts.top_k)
    elif opts.output == 'json':
        report = analyzer.json_report()
    elif opts.output == 'markdown':
        report = analyzer.markdown_report(top_k=opts.top_k)

    if opts.out_file:
        with open(opts.out_file, 'w') as fh:
            fh.write(report + '\n')
        print(f'Report written to {opts.out_file}')
    else:
        print(report)


if __name__ == '__main__':
    main()
