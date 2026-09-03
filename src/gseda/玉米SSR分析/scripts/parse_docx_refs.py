#!/usr/bin/env python3
"""Parse the two docx files: extract per-locus reference sequence + annotation
runs (primer binding sites by color/underline, SSR target by red bold).

Run-level annotation legend (from the docx header note):
  yellow        (hl=yellow)       = 原始正向引物 (GB original forward)
  green         (hl=green)        = 原始反向引物 (GB original reverse)
  cyan          (hl=cyan)         = 新正向引物1  (new F1)
  darkMagenta   (hl=darkMagenta)  = 新反向引物1  (new R1)
  wave          (u=wave)          = 新正向引物2  (new F2)
  thick         (u=thick)         = 新反向引物2  (new R2)
  red bold      (w:color FF0000 + w:b) = SSR 靶点
"""
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'


def parse_docx(path):
    """Return list of paragraphs from a .docx file. Each paragraph = list of runs.
    run = (text, attrs_dict) where attrs has keys: hl, u, color, bold.
    """
    with zipfile.ZipFile(path) as z:
        xml_data = z.read('word/document.xml')
    root = ET.fromstring(xml_data)
    paras = []
    for p in root.iter(f'{{{W}}}p'):
        runs = []
        for r in p.iter(f'{{{W}}}r'):
            t_el = r.find(f'{{{W}}}t')
            if t_el is None or not t_el.text:
                continue
            txt = t_el.text
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
            runs.append((txt, attrs))
        paras.append(runs)
    return paras


def classify(attrs):
    """Map run attrs -> annotation category."""
    if attrs.get('bold') and attrs.get('color') == 'FF0000':
        return 'ssr_target'
    hl = attrs.get('hl')
    if hl == 'yellow':
        return 'gb_fwd'
    if hl == 'green':
        return 'gb_rev'
    if hl == 'cyan':
        return 'new_f1'
    if hl == 'darkMagenta':
        return 'new_r1'
    u = attrs.get('u')
    if u == 'wave':
        return 'new_f2'
    if u == 'thick':
        return 'new_r2'
    return None


def parse_locus_refs(docx_path, locus_names=None):
    """Yield (locus, marker_name, seq, annotations) for each PM locus.
    annotations: list of (category, start, end, text)
    Handles PM20 whose reference sequence is glued to the title line.
    """
    paras = parse_docx(docx_path)
    # First pass: find locus title paragraphs and the text after "PMxx：name".
    locus_lines = []  # (para_index, locus, title_rest)
    for i, runs in enumerate(paras):
        text = ''.join(t for t, a in runs)
        m = re.match(r'PM(\d+)[：:]\s*(\S*)', text)
        if m:
            num = int(m.group(1))
            locus = f"PM{int(num):02d}"
            rest = m.group(2)
            known = (locus_names or {}).get(locus)
            if known:
                if rest.startswith(known):
                    rest = rest[len(known):]
                elif len(rest) >= len(known):
                    mm = sum(1 for a, b in zip(rest[:len(known)], known) if a != b)
                    if mm <= 2:
                        rest = rest[len(known):]
            locus_lines.append((i, locus, rest))

    results = []
    for k, (title_idx, locus, title_rest) in enumerate(locus_lines):
        # reference = text on title line (after name) + paragraphs after title until next title
        end_idx = locus_lines[k + 1][0] if k + 1 < len(locus_lines) else len(paras)
        seq_parts = [title_rest]
        annos = []  # (category, start, end, text)
        pos = len(title_rest)
        for pi in range(title_idx + 1, end_idx):
            for txt, attrs in paras[pi]:
                seq_parts.append(txt)
                cat = classify(attrs)
                if cat:
                    annos.append((cat, pos, pos + len(txt), txt))
                pos += len(txt)
        seq = ''.join(seq_parts)
        # merge contiguous runs of the same category
        merged = []
        for cat, s, e, t in annos:
            if merged and merged[-1][0] == cat and merged[-1][2] == s:
                merged[-1] = (cat, merged[-1][1], e, merged[-1][3] + t)
            else:
                merged.append((cat, s, e, t))
        marker = (locus_names or {}).get(locus, locus)
        results.append({'locus': locus, 'marker': marker, 'seq': seq, 'annos': merged})
    return results


def main():
    DOCX = [
        "/data1/ccs_data/20260813-corn-ssr-analysis/玉米SSR参考序列及新引物PM1-20.docx",
        "/data1/ccs_data/20260813-corn-ssr-analysis/玉米SSR序列及引物(PM21-40).docx",
    ]
    all_loci = {}
    for d in DOCX:
        for rec in parse_locus_refs(d):
            all_loci[rec['locus']] = rec
    print(f"解析位点数: {len(all_loci)}")
    for locus in sorted(all_loci, key=lambda x: int(x[2:])):
        rec = all_loci[locus]
        seq = re.sub(r'[^ACGT]', '', rec['seq'].upper())
        # merge contiguous same-category annotations
        merged = {}
        for cat, s, e, t in rec['annos']:
            merged.setdefault(cat, []).append((s, e, t))
        print(f"\n{locus} {rec['marker']:12s} seq_len={len(seq)}")
        for cat in ['gb_fwd', 'gb_rev', 'new_f1', 'new_r1', 'new_f2', 'new_r2', 'ssr_target']:
            if cat in merged:
                for s, e, t in merged[cat]:
                    t_seq = re.sub(r'[^ACGT]', '', t.upper())
                    print(f"    {cat:10s} pos{s}-{e} : {t_seq[:45]}")
            else:
                print(f"    {cat:10s} (无)")


if __name__ == "__main__":
    main()
