#!/usr/bin/env python3
"""prepare_data.py — 玉米 SSR 数据准备：把 docx 参考序列 + xlsx 引物 + 国标信息
整理成标准生信文件（FASTA / BED / TSV）。

═══ 引物溯源（2026-08 确认）═══
docx 中每个位点标注【三对引物】：原对(黄/绿) + 新1对(青/紫) + 新2对(波浪/粗下划线)。
但【真实使用】的引物 = 合成单 xlsx 中的 F1/R1/F2/R2（40 位点 × 4 = 160 条）：
  - 「【海南】合成-海南测序部.xlsx」 含 PM1-40 全部 160 条（最全，基准）
  - 「【海南】合成PM21-PM40.xlsx」   PM21-40 80 条（与测序部一致，交叉验证）
  - 「SSR玉米引物汇总表PM1-20.xlsx」 PM1-20 80 条汇总（与测序部一致，包含于其内）
三份文件重叠区序列 0 差异；合成单中无 GB 原对（原对未合成，仅为参考标注）。
→ all_primers.tsv / all_primers.fa = 真实使用引物集（160 条）。
→ gb_primers.tsv = docx 原对标注 + 国标表1（供对标国标用，非平台使用）。

产出目录：<repo>/gseda/src/gseda/玉米SSR分析/data/
  refs/pm_all.fa           40 位点参考序列
  primers/all_primers.tsv  160 条新引物（含结合位点坐标）
  primers/all_primers.fa   引物 FASTA
  primers/gb_primers.tsv   国标表1 40 对原始引物
  annot/ssr_targets.bed    SSR 靶点（0-based）
  annot/primer_sites.bed   所有引物结合位点
  annot/locus_info.tsv     位点总表（含国标 Bin / 参照样品）
  README.md                说明

引物结合位点定位策略：
  以「引物序列匹配参考序列」为主（正向引物正向匹配；反向引物以反向互补匹配），
  docx 颜色标注作为辅助/交叉验证。
SSR 靶点：算法检测（最长 2-6bp 串联重复）+ docx 红色标注验证。
"""
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import edlib

from ssr_common import classify_region

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

SRC = Path("/data1/ccs_data/20260813-corn-ssr-analysis")
OUT = Path(__file__).resolve().parent.parent / "data"

RC_MAP = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(RC_MAP)[::-1]


# ---------------------------------------------------------------------------
# 1. docx 解析
# ---------------------------------------------------------------------------
def docx_paras(path):
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read('word/document.xml'))
    paras = []
    for p in root.iter(f'{{{W}}}p'):
        runs = []
        for r in p.iter(f'{{{W}}}r'):
            t_el = r.find(f'{{{W}}}t')
            if t_el is None or not t_el.text:
                continue
            attrs = {}
            rpr = r.find(f'{{{W}}}rPr')
            if rpr is not None:
                hl = rpr.find(f'{{{W}}}highlight')
                if hl is not None:
                    attrs['hl'] = hl.get(f'{{{W}}}val')
                u = rpr.find(f'{{{W}}}u')
                if u is not None:
                    attrs['u'] = u.get(f'{{{W}}}val')
                color = rpr.find(f'{{{W}}}color')
                if color is not None:
                    attrs['color'] = color.get(f'{{{W}}}val')
                if rpr.find(f'{{{W}}}b') is not None:
                    attrs['bold'] = True
            runs.append((t_el.text, attrs))
        paras.append(runs)
    return paras


def classify(attrs):
    # SSR 标注颜色：PM1-20 docx 用 EE0000，PM21-40 docx 用 FF0000（两种都要认）
    if attrs.get('bold') and attrs.get('color', '').upper() in ('FF0000', 'EE0000'):
        return 'ssr_target'
    hl = attrs.get('hl')
    if hl == 'yellow': return 'gb_fwd'
    if hl == 'green': return 'gb_rev'
    if hl == 'cyan': return 'new_f1'
    if hl == 'darkMagenta': return 'new_r1'
    u = attrs.get('u')
    if u == 'wave': return 'new_f2'
    if u == 'thick': return 'new_r2'
    return None


def parse_docx_refs(docx_path, locus_names):
    """Return dict locus -> {marker, seq, annos:[(cat,start,end,text)]}.
    Handles PM20 whose reference sequence is glued to the title line.
    """
    paras = docx_paras(docx_path)
    locus_lines = []
    for i, runs in enumerate(paras):
        text = ''.join(t for t, a in runs)
        m = re.match(r'PM(\d+)[：:]\s*(\S*)', text)
        if m:
            num = int(m.group(1))
            locus = f"PM{int(num):02d}"
            rest = m.group(2)
            known = locus_names.get(locus)
            if known:
                if rest.startswith(known):
                    rest = rest[len(known):]
                elif len(rest) >= len(known):
                    mm = sum(1 for a, b in zip(rest[:len(known)], known) if a != b)
                    if mm <= 2:
                        rest = rest[len(known):]
            locus_lines.append((i, locus, rest))

    out = {}
    for k, (title_idx, locus, title_rest) in enumerate(locus_lines):
        name = locus_names.get(locus, locus)
        end_idx = locus_lines[k + 1][0] if k + 1 < len(locus_lines) else len(paras)
        seq_parts = [title_rest]
        annos, pos = [], len(title_rest)
        # 标题行内（title_rest 部分）的标注也要收——PM20 参考序列整段都粘在标题行上，
        # 它的 SSR/引物标注都在标题段内，旧逻辑（只扫标题之后的段落）会全部漏掉。
        ttext = ''.join(t for t, a in paras[title_idx])
        m0 = re.match(r'PM(\d+)[：:]\s*', ttext)
        mark_end = m0.end() if m0 else 0
        seg = ttext[mark_end:]
        removed = 0
        if seg.startswith(name):
            removed = len(name)
        elif len(seg) >= len(name):
            mm = sum(1 for a, b in zip(seg[:len(name)], name) if a != b)
            if mm <= 2:
                removed = len(name)
        title_off = mark_end + removed      # title_rest 在标题段文本中的起点
        pos0 = 0
        for txt, attrs in paras[title_idx]:
            cat = classify(attrs)
            if cat and pos0 >= title_off:
                rel = pos0 - title_off
                annos.append((cat, rel, rel + len(txt), txt))
            pos0 += len(txt)
        for pi in range(title_idx + 1, end_idx):
            for txt, attrs in paras[pi]:
                seq_parts.append(txt)
                cat = classify(attrs)
                if cat:
                    annos.append((cat, pos, pos + len(txt), txt))
                pos += len(txt)
        seq = ''.join(seq_parts)
        merged = []
        for cat, s, e, t in annos:
            if merged and merged[-1][0] == cat and merged[-1][2] == s:
                merged[-1] = (cat, merged[-1][1], e, merged[-1][3] + t)
            else:
                merged.append((cat, s, e, t))
        out[locus] = {'marker': name, 'seq': seq, 'annos': merged}
    return out


# ---------------------------------------------------------------------------
# 2. xlsx 解析
# ---------------------------------------------------------------------------
XLSX_NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'


def xlsx_shared_strings(z):
    if 'xl/sharedStrings.xml' not in z.namelist():
        return []
    root = ET.fromstring(z.read('xl/sharedStrings.xml'))
    return [''.join(t.text or '' for t in si.iter(f'{{{XLSX_NS}}}t'))
            for si in root.iter(f'{{{XLSX_NS}}}si')]


def xlsx_sheet_rows(z, shared, sheet='xl/worksheets/sheet1.xml'):
    root = ET.fromstring(z.read(sheet))
    rows = []
    for row in root.iter(f'{{{XLSX_NS}}}row'):
        cells = {}
        for c in row.iter(f'{{{XLSX_NS}}}c'):
            col = re.match(r'[A-Z]+', c.get('r')).group()
            v = c.find(f'{{{XLSX_NS}}}v')
            is_ = c.find(f'{{{XLSX_NS}}}is')
            if v is not None and v.text:
                cells[col] = shared[int(v.text)] if c.get('t') == 's' else v.text
            elif is_ is not None:
                cells[col] = ''.join(t.text or '' for t in is_.iter(f'{{{XLSX_NS}}}t'))
            else:
                cells[col] = ''
        if cells:
            rows.append(cells)
    return rows


def parse_xlsx_primers():
    """Extract the 160 ACTUALLY-USED primers (F1/R1/F2/R2) from the synthesis
    order sheets. 测序部 sheet covers all 40 loci; the PM21-40 sheet and the
    PM1-20 summary are redundant copies (0 sequence diff) kept as cross-check."""
    primers = {}  # name -> seq
    for fname in ['【海南】合成-海南测序部.xlsx', 'SSR玉米引物汇总表PM1-20.xlsx', '【海南】合成PM21-PM40.xlsx']:
        p = SRC / fname
        if not p.exists():
            continue
        with zipfile.ZipFile(p) as z:
            shared = xlsx_shared_strings(z)
            rows = xlsx_sheet_rows(z, shared)
        for r in rows:
            name = r.get('B', '')
            seq = (r.get('C', '') or '').upper()
            if re.match(r'^PM\d+[FR][12]$', name) and re.fullmatch(r'[ACGTN]+', seq):
                primers[name] = seq
    return primers


# ---------------------------------------------------------------------------
# 3. 标准 md 解析
# ---------------------------------------------------------------------------
def parse_standard():
    std = (SRC / "corn-ssr-standard.md").read_text(encoding='utf-8')
    # 表1: | PM01 | bnlg439w1 | 1.03 | 正向: ...<br>反向: ... | NED |
    gb_primers = {}
    for m in re.finditer(
            r'\|\s*(PM\d+)\s*\|\s*([\w.]+)\s*\|\s*([\d.]+)\s*\|\s*正向:\s*([ACGT]+)\s*<br>\s*反向:\s*([ACGT]+)\s*\|\s*(\w+)\s*\|',
            std):
        gb_primers[m.group(1)] = {
            'marker': m.group(2), 'chr': m.group(3),
            'fwd': m.group(4), 'rev': m.group(5), 'dye': m.group(6),
        }
    # 附录B: | PM01 | bnlg439w1 | 320～368 | 320 | 0.007 | 绵单1号 | 320/350 |
    bin_info = {}
    for m in re.finditer(r'\|\s*(PM\d+)\s*\|\s*([\w.]+)\s*\|\s*(\d+～\d+)\s*\|', std):
        locus = m.group(1)
        bin_info.setdefault(locus, {'range': m.group(3), 'refs': set()})
    # reference samples: rows with 6th col (name) and 7th col (alleles)
    for m in re.finditer(r'\|\s*(PM\d+)\s*\|\s*\|\s*\|?\s*\|\s*\|?\s*\|\s*([一-鿿\w]+)\s*\|\s*([\d/]+)\s*\|', std):
        locus = m.group(1)
        if locus in bin_info:
            bin_info[locus]['refs'].add(f"{m.group(2)}:{m.group(3)}")
    return gb_primers, bin_info


# ---------------------------------------------------------------------------
# 4. primer locating by sequence match
# ---------------------------------------------------------------------------
def locate_primer(seq, ref, is_reverse):
    """Return (start, end, strand) of primer binding site on reference, or None.
    Forward primer: match seq on ref (+). Reverse primer: match rc(seq) on ref (-).
    Allow up to 2 mismatches via edlib HW.
    """
    probe = rc(seq) if is_reverse else seq
    r = edlib.align(probe, ref, mode="HW", task="path")
    if r['editDistance'] is None:
        return None
    # HW: probe fully aligns inside ref. location on ref = locations[0]
    st, en = r['locations'][0]
    # end is exclusive-ish; edlib end is last matched index inclusive? verify
    return st, en + 1, '-' if is_reverse else '+'


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    std_text = (SRC / "corn-ssr-standard.md").read_text(encoding='utf-8')
    # 表1格式: | PM01 | bnlg439w1 | 1.03 | ...
    locus_names = {m.group(1): m.group(2)
                   for m in re.finditer(r'\|\s*(PM\d+)\s*\|\s*([\w.]+)\s*\|', std_text)}

    # 1. docx refs
    refs = {}
    for f in ['玉米SSR参考序列及新引物PM1-20.docx', '玉米SSR序列及引物(PM21-40).docx']:
        refs.update(parse_docx_refs(SRC / f, locus_names))
    clean = lambda s: re.sub(r'[^ACGT]', '', s.upper())
    for locus in refs:
        refs[locus]['seq'] = clean(refs[locus]['seq'])
    refs = {k: v for k, v in refs.items() if 700 <= len(v['seq']) <= 2500}
    print(f"docx 参考序列: {len(refs)} 个")

    # 2. xlsx primers
    primers = parse_xlsx_primers()
    print(f"xlsx 新引物: {len(primers)} 条")
    # map primer name -> locus/type
    primer_meta = {}
    for name in primers:
        m = re.match(r'PM(\d+)([FR])([12])', name)
        primer_meta[name] = (f"PM{int(m.group(1)):02d}", m.group(2), m.group(3))

    # 3. standard
    gb_primers, bin_info = parse_standard()
    print(f"国标表1引物: {len(gb_primers)} 对; 附录B Bin: {len(bin_info)} 个")

    # 4. locate new primers on references
    primer_sites = []  # (locus, marker, primer, type, start, end, strand, seq)
    unmatched = []
    for name, seq in primers.items():
        locus, fwdrev, num = primer_meta[name]
        ref = refs.get(locus, {}).get('seq')
        if not ref:
            unmatched.append((name, 'no_ref'))
            continue
        is_rev = (fwdrev == 'R')
        loc = locate_primer(seq, ref, is_rev)
        if loc is None:
            unmatched.append((name, 'no_match'))
            continue
        st, en, strand = loc
        primer_sites.append((locus, refs[locus]['marker'], name, fwdrev + num, st, en, strand, seq))

    print(f"新引物定位成功: {len(primer_sites)}/{len(primers)}")
    for u in unmatched[:20]:
        print(f"  ⚠️ 未定位: {u}")

    # 5. SSR targets: docx red annotation (authoritative) + algorithm fallback.
    #    注意：
    #    - 红色标注颜色有两种（PM1-20=EE0000，PM21-40=FF0000），classify 已都识别。
    #    - 部分位点有【多个】红色标注（不连续的独立 SSR 区），全部提取并按出现顺序
    #      编号 target（1 = 主靶点/国标对标）。
    #    - 区域内部用 classify_region 区分「复合」（PM25 CT+AC）与「简单/断裂简单」
    #      （PM05 AG×13+AA+AG×3 → 长度法 CN=17）；ssr_targets 每区域一行，
    #      comps 保留各分量坐标供 locus_info 每分量一行输出。
    ssr_targets = []   # (locus, marker, start, end, motif_display, ref_cn_display, target, comps)
    for locus in sorted(refs, key=lambda x: int(x[2:])):
        seq = refs[locus]['seq']
        marker = refs[locus]['marker']
        red = [a for a in refs[locus].get('annos', []) if a[0] == 'ssr_target']
        if red:
            for ti, (_cat, s, e, t_orig) in enumerate(red, 1):
                kind, comps = classify_region(seq[s:e], s, e)
                if kind == 'compound':
                    motif = '+'.join(f"{u}{c}" for u, c, _, _ in comps)
                    dom = max(comps, key=lambda x: x[1] * len(x[0]))
                    ref_cn = round((e - s) / len(dom[0]))
                else:
                    motif, ref_cn = comps[0][0], comps[0][1]
                ssr_targets.append((locus, marker, s, e, motif, ref_cn, ti, comps))
        else:
            best = None
            for u in range(2, 7):
                for start in range(len(seq) - u + 1):
                    unit = seq[start:start + u]
                    j, c = start, 0
                    while j + u <= len(seq) and seq[j:j + u] == unit:
                        c += 1; j += u
                    if c >= 4 and (best is None or c > best[1]):
                        best = (unit, c, start)
            if best:
                unit, cnt, st = best
                en = st + len(unit) * cnt
                ssr_targets.append((locus, marker, st, en, unit, cnt, 1, [(unit, cnt, st, en)]))

    # 6. locate GB primers: docx yellow/green annotations (authoritative) +
    #    seq-match fallback for PM20 (no annotation).
    gb_sites = []
    for locus in sorted(gb_primers, key=lambda x: int(x[2:])):
        g = gb_primers[locus]
        ref = refs.get(locus, {}).get('seq')
        if not ref:
            continue
        annos = refs[locus].get('annos', [])
        fwd_anno = [a for a in annos if a[0] == 'gb_fwd']
        rev_anno = [a for a in annos if a[0] == 'gb_rev']

        if fwd_anno:
            fs, fe = fwd_anno[0][1], fwd_anno[0][2]
            ref_fwd = ref[fs:fe]
        else:
            f_loc = locate_primer(g['fwd'], ref, False)
            fs, fe = (f_loc[0], f_loc[1]) if f_loc else ('', '')
            ref_fwd = ref[fs:fe] if f_loc else ''

        if rev_anno:
            rs, re_ = rev_anno[0][1], rev_anno[0][2]
            ref_rev = rc(ref[rs:re_])          # reference holds rc(rev primer)
        else:
            r_loc = locate_primer(g['rev'], ref, True)
            rs, re_ = (r_loc[0], r_loc[1]) if r_loc else ('', '')
            ref_rev = rc(ref[rs:re_]) if r_loc else ''

        gb_sites.append((locus, g['marker'], g['chr'], g['fwd'], g['rev'], g['dye'],
                         fs, fe, rs, re_, ref_fwd, ref_rev))

    # ---- write outputs ----
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'refs').mkdir(exist_ok=True)
    (OUT / 'primers').mkdir(exist_ok=True)
    (OUT / 'annot').mkdir(exist_ok=True)

    # pm_all.fa
    with open(OUT / 'refs' / 'pm_all.fa', 'w') as fh:
        for locus in sorted(refs, key=lambda x: int(x[2:])):
            r = refs[locus]
            chrm = gb_primers.get(locus, {}).get('chr', '?')
            fh.write(f">{locus}|{r['marker']}|chr{chrm}|len{len(r['seq'])}\n")
            for i in range(0, len(r['seq']), 80):
                fh.write(r['seq'][i:i + 80] + '\n')

    # all_primers.tsv
    with open(OUT / 'primers' / 'all_primers.tsv', 'w') as fh:
        fh.write('locus\tmarker\tprimer_name\ttype\tseq\tlen\tref_start\tref_end\tstrand\n')
        for locus, marker, name, ptype, st, en, strand, seq in sorted(
                primer_sites, key=lambda x: (int(x[0][2:]), x[2])):
            fh.write(f"{locus}\t{marker}\t{name}\t{ptype}\t{seq}\t{len(seq)}\t{st}\t{en}\t{strand}\n")

    # all_primers.fa
    with open(OUT / 'primers' / 'all_primers.fa', 'w') as fh:
        for locus, marker, name, ptype, st, en, strand, seq in sorted(
                primer_sites, key=lambda x: (int(x[0][2:]), x[2])):
            fh.write(f">{name}|{locus}|{marker}|{ptype}|{strand}\n{seq}\n")

    # gb_primers.tsv
    with open(OUT / 'primers' / 'gb_primers.tsv', 'w') as fh:
        fh.write('locus\tmarker\tchr\tgb_fwd\tgb_rev\tfwd_start\tfwd_end\trev_start\trev_end\t'
                 'ref_fwd\tref_rev\tdye\tbin_range\tnote\n')
        for locus, marker, chrm, fwd, rev, dye, fs, fe, rs, re_, ref_fwd, ref_rev in gb_sites:
            br = bin_info.get(locus, {}).get('range', '')
            note = ''
            if not (fs and fe):
                note = 'gb_fwd_未定位'
            if not (rs and re_):
                note += ';gb_rev_未定位'
            fh.write(f"{locus}\t{marker}\t{chrm}\t{fwd}\t{rev}\t{fs}\t{fe}\t{rs}\t{re_}\t"
                     f"{ref_fwd}\t{ref_rev}\t{dye}\t{br}\t{note}\n")

    # ssr_targets.bed
    with open(OUT / 'annot' / 'ssr_targets.bed', 'w') as fh:
        fh.write('#locus\tmarker\tstart\tend\tmotif\tref_cn\ttarget\n')
        for locus, marker, st, en, motif, ref_cn, ti, comps in sorted(
                ssr_targets, key=lambda x: (int(x[0][2:]), x[6])):
            fh.write(f"{locus}\t{marker}\t{st}\t{en}\t{motif}\t{ref_cn}\t{ti}\n")

    # primer_sites.bed (new + GB)
    with open(OUT / 'annot' / 'primer_sites.bed', 'w') as fh:
        fh.write('#locus\tstart\tend\tprimer_name\ttype\tstrand\tseq\n')
        for locus, marker, name, ptype, st, en, strand, seq in sorted(
                primer_sites, key=lambda x: (int(x[0][2:]), x[2])):
            fh.write(f"{locus}\t{st}\t{en}\t{name}\t{ptype}\t{strand}\t{seq}\n")
        for locus, marker, chrm, fwd, rev, dye, fs, fe, rs, re_, ref_fwd, ref_rev in gb_sites:
            if fs and fe:
                fh.write(f"{locus}\t{fs}\t{fe}\tGB_{locus}_F\tGB_fwd\t+\t{ref_fwd}\n")
            if rs and re_:
                fh.write(f"{locus}\t{rs}\t{re_}\tGB_{locus}_R\tGB_rev\t-\t{rev}\n")

    # locus_info.tsv（每 SSR 分量一行：复合拆开、断裂简单合并；target = 区域编号）
    with open(OUT / 'annot' / 'locus_info.tsv', 'w') as fh:
        fh.write('locus\tmarker\tchr\ttarget\tssr_motif\tref_cn\tssr_start\tssr_end\tgb_bin_range\tgb_ref_samples\n')
        for locus in sorted(refs, key=lambda x: int(x[2:])):
            marker = refs[locus]['marker']
            chrm = gb_primers.get(locus, {}).get('chr', '?')
            br = bin_info.get(locus, {}).get('range', '')
            refs_ = ';'.join(sorted(bin_info.get(locus, {}).get('refs', [])))
            for tgt in ssr_targets:
                if tgt[0] != locus:
                    continue
                for u, cn, cs, ce in tgt[7]:
                    fh.write(f"{locus}\t{marker}\t{chrm}\t{tgt[6]}\t{u}\t{cn}\t{cs}\t{ce}\t{br}\t{refs_}\n")

    print("\n✅ 输出已生成:")
    for f in sorted(OUT.rglob('*')):
        if f.is_file() and 'README' not in f.name:
            print(f"  {f.relative_to(OUT)}  ({f.stat().st_size}B)")


if __name__ == "__main__":
    main()
