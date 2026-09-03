#!/usr/bin/env python3
"""bridge_cn2bp.py — 桥A：CN → 片段长度(bp) → 国标Bin 换算 + 附录B验证

背景：CE（毛细管电泳）测的是 GB 引物间 PCR 产物片段长度(bp)，国标等位基因 Bin 也是 bp。
本脚本把平台的 reads CN（重复数）换算成片段长度，映射到附录B Bin，并用附录B 验证。

两套独立换算（相互验证）：
  1. reads 实测法（主）：在每条 read 里找 GB-F / GB-R 结合位点锚 → 物理片段长度（= CE 测的东西）。
     对嵌套位点（GB 引物在新引物扩增子内，37/40）有效，绕开参考坐标误差。
  2. CN 法（交叉验证）：fragment = flanking + Σ(GB扩增子内所有SSR靶点的 CN×unit)，
     其中 flanking = GB扩增子长 - Σ(ref_cn×unit)。多靶点位点必须把所有 SSR 求和。

用法：
  python3 bridge_cn2bp.py                     # 全量
  python3 bridge_cn2bp.py --only barcode201-0 # 单样品调试
  python3 bridge_cn2bp.py --threads 16
"""
import argparse
import re
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import edlib

_RC = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(_RC)[::-1]


# ---------------------------------------------------------------------------
# 数据加载（worker 初始化）
# ---------------------------------------------------------------------------
_REFSEQ = {}
_GB = {}            # locus -> dict(fs,fe,rs,re, gb_fwd, gb_rev)
_SSR = defaultdict(list)   # locus -> [(target,start,end,unit,ref_cn)]
_BINS = {}          # locus -> dict(range, lo, hi, alleles)
_CN = {}            # (sample,locus,target,unit) -> mode_cn
_FQ_DIR = ''
_ONLY = ''
_THREADS = 4


def _load(data_dir, cn_dir):
    # 参考序列
    cur = None
    with open(data_dir / 'refs' / 'pm_all.fa') as fh:
        for line in fh:
            if line.startswith('>'):
                if cur:
                    _REFSEQ[cur] = ''.join(seq)
                cur = line[1:].split('|')[0]
                seq = []
            else:
                seq.append(line.strip())
        if cur:
            _REFSEQ[cur] = ''.join(seq)

    # GB 引物（docx 黄/绿标注坐标 + 标准序列 + 参考实际序列）
    with open(data_dir / 'primers' / 'gb_primers.tsv') as fh:
        for line in fh:
            if line.startswith('locus'):
                continue
            p = line.rstrip('\n').split('\t')
            if not (p[5] and p[8]):
                continue
            _GB[p[0]] = dict(fs=int(p[5]), fe=int(p[6]), rs=int(p[7]), re=int(p[8]),
                             gb_fwd=p[3], gb_rev=p[4],
                             ref_fwd=p[9], ref_rev=p[10])

    # SSR 靶点（全部，含多靶点）
    with open(data_dir / 'annot' / 'ssr_targets.bed') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            _SSR[p[0]].append((int(p[6]), int(p[2]), int(p[3]), p[4], int(p[5])))

    # 附录B：范围 + 等位基因片段
    std = open(Path('/data1/ccs_data/20260813-corn-ssr-analysis/corn-ssr-standard.md'),
               encoding='utf-8').read()
    for m in re.finditer(r'\|\s*(PM\d+)\s*\|\s*[\w.]+\s*\|\s*(\d+～\d+)\s*\|\s*(\d+)\s*\|\s*([\d.]+)',
                         std):
        pm, rng, frag = m.group(1), m.group(2), int(m.group(3))
        lo, hi = map(int, rng.split('～'))
        b = _BINS.setdefault(pm, {'lo': lo, 'hi': hi, 'alleles': set()})
        b['alleles'].add(frag)

    # cn_summary → (sample, locus, target, unit) -> mode_cn × unit长度(bp)
    with open(cn_dir / 'cn_summary.tsv') as fh:
        for line in fh:
            if line.startswith('sample'):
                continue
            p = line.rstrip('\n').split('\t')
            if p[6] != '.':
                _CN[(p[0], p[1], int(p[2]), p[4])] = int(p[10]) * len(p[4])


def _worker_init(data_dir, cn_dir, fq_dir, only):
    _load(Path(data_dir), Path(cn_dir))
    global _FQ_DIR, _ONLY
    _FQ_DIR = fq_dir
    _ONLY = only


def _hw(probe, text, maxd=2):
    r = edlib.align(probe, text, mode='HW', task='path')
    if r['editDistance'] is None or r['editDistance'] > maxd:
        return None
    return r['locations'][0]


# ---------------------------------------------------------------------------
# 第1步：修复 GB 结合位点坐标
# ---------------------------------------------------------------------------
def fix_gb_sites():
    """对 docx 片段不在附录B 范围的位点，用标准全长引物在参考中重新定位。
    返回 locus -> dict(fs,fe,rs,re, fix_note)。"""
    out = {}
    for locus, g in _GB.items():
        ref = _REFSEQ.get(locus, '')
        fs, fe, rs, re_ = g['fs'], g['fe'], g['rs'], g['re']
        ref_frag = re_ - fs
        b = _BINS.get(locus)
        in_range = b is not None and b['lo'] <= ref_frag <= b['hi']
        note = ''
        if b and not in_range:
            # 尝试标准引物定位
            f_loc = _hw(g['gb_fwd'], ref, maxd=3)
            r_loc = _hw(rc(g['gb_rev']), ref, maxd=3)
            new_frag = ref_frag
            if f_loc:
                fs, fe = f_loc[0], f_loc[1] + 1
            if r_loc:
                rs, re_ = r_loc[0], r_loc[1] + 1
            new_frag = re_ - fs
            if b['lo'] <= new_frag <= b['hi']:
                note = f'标准引物修正 {ref_frag}→{new_frag}'
            else:
                note = f'坐标仍超范围({new_frag}, 范围{b["lo"]}-{b["hi"]})'
        out[locus] = dict(fs=fs, fe=fe, rs=rs, re=re_, fix=note, ref_frag=re_ - fs,
                          gb_fwd=g['gb_fwd'], gb_rev=g['gb_rev'],
                          ref_fwd=ref[fs:fe], ref_rev=rc(ref[rs:re_]))
    return out


# ---------------------------------------------------------------------------
# 第2步：reads 实测片段
# ---------------------------------------------------------------------------
def read_fragment(read, g):
    """用 GB 引物结合位点实际序列（docx ref_fwd/ref_rev，Sanger 验证）在 read 里定位，
    量物理片段长度 bp = 测序结果中的 CE 片段。
    正向引物结合位点=ref_fwd 在 5'；反向引物结合位点=rc(ref_rev) 在 3'。"""
    f = _hw(g['ref_fwd'], read, maxd=3)
    if f is None:
        return None
    r = _hw(rc(g['ref_rev']), read, maxd=3)
    if r is None or r[0] <= f[0]:
        return None
    return (r[0] + len(g['ref_rev'])) - f[0]


def process_sample(sample):
    if _ONLY and sample != _ONLY:
        return None
    sdir = Path(_FQ_DIR) / sample
    if not sdir.is_dir():
        return None
    sites = fix_gb_sites()
    rows = []
    for fq in sorted(sdir.glob('*.fastq')):
        locus = fq.stem
        g = sites.get(locus)
        ref = _REFSEQ.get(locus)
        if not g or not ref:
            continue
        frags = []
        lines = fq.read_text().splitlines()
        for i in range(len(lines) // 4):
            fr = read_fragment(lines[i * 4 + 1], g)
            if fr:
                frags.append(fr)
        if not frags:
            continue
        c = Counter(frags)
        mode_frag, mode_n = c.most_common(1)[0]
        n = sum(c.values())
        pct = mode_n / n * 100
        rows.append((locus, n, mode_frag, pct, dict(sorted(c.items()))))
    return sample, rows


def cn_fragment(sample, locus, g):
    """CN 法片段：fragment = flanking + Σ(GB扩增子内所有SSR靶点的重复总长bp)。
    复合靶点（PM25 target1）在 cn_summary 里拆成多个分量，需各分量 CN×unit 求和。
    返回 (flanking, ssr_bp, frag, all_real)。
    all_real=False 表示有靶点无该样品的 CN 数据（用了参考默认长度），
    此时 frag 仅供参考，**不应计入一致性验证**（低深度过滤后常见）。"""
    in_ssr = [t for t in _SSR[locus] if t[1] >= g['fs'] and t[2] <= g['re']]
    if not in_ssr:
        return None, None, None, False
    ref_ssr = sum(e - s for _, s, e, _, _ in in_ssr)   # 靶点坐标长度之和（复合也准）
    flanking = (g['re'] - g['fs']) - ref_ssr
    ssr_bp = 0
    all_real = True
    for tg, s, e, unit, refcn in in_ssr:
        rows = [cn for (sm, lc, t, u), cn in _CN.items()
                if sm == sample and lc == locus and t == tg]
        if rows:
            ssr_bp += sum(rows)                 # 各分量 CN×unit 已在加载时乘好
        else:
            ssr_bp += (e - s)                    # 该靶点无数据 → 用参考长度（默认参考等位基因）
            all_real = False
    return flanking, ssr_bp, flanking + ssr_bp, all_real


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    global _FQ_DIR, _ONLY, _THREADS
    ap = argparse.ArgumentParser(description='CN→bp→国标Bin 换算（桥A）')
    ap.add_argument('--data-dir', default=str(Path(__file__).resolve().parent.parent / 'data'))
    ap.add_argument('--cn-dir', default=None, help='cn_summary 所在目录（默认 <data-dir>/cn）')
    ap.add_argument('--split-dir', default=None, help='拆分 fastq 目录（默认 <data-dir>/assign/split）')
    ap.add_argument('--out', default=None, help='输出目录（默认 <data-dir>/bridge）')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--only', default=None)
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    cn_dir = Path(args.cn_dir) if args.cn_dir else data_dir / 'cn'
    _FQ_DIR = str(Path(args.split_dir) if args.split_dir else data_dir / 'assign' / 'split')
    _ONLY = args.only
    out_dir = Path(args.out) if args.out else data_dir / 'bridge'

    _load(data_dir, cn_dir)   # 主进程也要加载（fix_gb_sites / 汇总用）

    samples = sorted(d.name for d in Path(_FQ_DIR).iterdir() if d.is_dir())
    if args.only:
        samples = [s for s in samples if args.only in s]
    print(f'待处理样品: {len(samples)}')

    with Pool(args.threads, initializer=_worker_init,
              initargs=(str(data_dir), str(cn_dir), _FQ_DIR, args.only)) as pool:
        results = pool.map(process_sample, samples, chunksize=1)
    results = [r for r in results if r]

    # 主进程：读法片段汇总
    readmode = {}   # (sample,locus) -> (n, mode_frag, pct)
    for sample, rows in results:
        for locus, n, mf, pct, dist in rows:
            readmode[(sample, locus)] = (n, mf, pct)

    out_dir.mkdir(parents=True, exist_ok=True)
    sites = fix_gb_sites()

    # ---- fragment_summary.tsv ----
    with open(out_dir / 'fragment_summary.tsv', 'w') as fh:
        fh.write('sample\tlocus\tread_depth\tread_frag_bp\tread_pct\t'
                 'cn_frag_bp\tflanking\tssr_total_bp\tdiff_bp\tgb_range\tin_range\tcn_real\tnote\n')
        for sample in sorted({s for s, _ in readmode} | {s for s, l, t, u in _CN}):
            for locus in sorted(_GB, key=lambda x: int(x[2:])):
                key = (sample, locus)
                b = _BINS.get(locus)
                rng = f"{b['lo']}-{b['hi']}" if b else '?'
                flanking, ssr_bp, cn_frag, all_real = cn_fragment(sample, locus, sites[locus])
                creal = 'Y' if all_real else ('N' if cn_frag is not None else '')
                rd = readmode.get(key)
                if rd:
                    n, mf, pct = rd
                    diff = abs(mf - cn_frag) if cn_frag is not None else ''
                    in_range = b is not None and (b['lo'] <= mf <= b['hi']) if rd else ''
                    fh.write(f"{sample}\t{locus}\t{n}\t{mf}\t{pct:.0f}\t"
                             f"{cn_frag if cn_frag is not None else ''}\t{flanking if flanking is not None else ''}"
                             f"\t{ssr_bp if ssr_bp is not None else ''}\t{diff}\t{rng}\t"
                             f"{'Y' if in_range else 'N'}\t{creal}\t{sites[locus]['fix']}\n")
                elif cn_frag is not None:
                    in_range = b is not None and (b['lo'] <= cn_frag <= b['hi'])
                    fh.write(f"{sample}\t{locus}\t0\t\t\t{cn_frag}\t{flanking}\t{ssr_bp}\t\t"
                             f"{rng}\t{'Y' if in_range else 'N'}\t{creal}\t{sites[locus]['fix']}\n")

    # ---- gb_sites.tsv ----
    with open(out_dir / 'gb_sites.tsv', 'w') as fh:
        fh.write('locus\tfwd_start\tfwd_end\trev_start\trev_end\tref_frag_bp\tgb_range\tfix_note\n')
        for locus in sorted(sites, key=lambda x: int(x[2:])):
            g = sites[locus]
            b = _BINS.get(locus)
            rng = f"{b['lo']}-{b['hi']}" if b else '?'
            fh.write(f"{locus}\t{g['fs']}\t{g['fe']}\t{g['rs']}\t{g['re']}\t{g['ref_frag']}\t"
                     f"{rng}\t{g['fix']}\n")

    # ---- validation.tsv ----
    with open(out_dir / 'validation.tsv', 'w') as fh:
        fh.write('locus\tunit\tgb_range\tref_frag_bp\tref_in_range\t'
                 'read_frag_mode\tread_in_range\tcn_vs_read_ok\tallele_grid_ok\tnote\n')
        for locus in sorted(_GB, key=lambda x: int(x[2:])):
            g = sites[locus]
            b = _BINS.get(locus)
            rng = f"{b['lo']}-{b['hi']}" if b else '?'
            ref_in = b is not None and b['lo'] <= g['ref_frag'] <= b['hi']
            # reads 众数片段（跨样品取中位）
            rd_frags = [mf for (s, l), (n, mf, p) in readmode.items() if l == locus]
            if rd_frags:
                rd_mode = sorted(rd_frags)[len(rd_frags) // 2]
                rd_in = b is not None and b['lo'] <= rd_mode <= b['hi']
            else:
                rd_mode, rd_in = '', ''
            # 两法一致性（样品级，多数 |diff|<=3；只统计该样品有真实 CN 数据的，排除默认回退）
            n_ok = n_tot = 0
            for (s, l), (n, mf, p) in readmode.items():
                if l != locus:
                    continue
                _, _, cn_frag, all_real = cn_fragment(s, locus, g)
                if cn_frag is None or not all_real:
                    continue
                n_tot += 1
                if abs(mf - cn_frag) <= 3:
                    n_ok += 1
            cn_ok = f"{n_ok}/{n_tot}" if n_tot else ''
            # 等位基因网格：reads 片段是否落在附录B 等位基因 ±2bp
            alleles = b['alleles'] if b else set()
            grid_ok = ''
            if rd_frags and alleles:
                hits = sum(1 for f in rd_frags if any(abs(f - a) <= 2 for a in alleles))
                grid_ok = f"{hits}/{len(rd_frags)}"
            fh.write(f"{locus}\t\t{rng}\t{g['ref_frag']}\t{'Y' if ref_in else 'N'}\t"
                     f"{rd_mode}\t{rd_in}\t{cn_ok}\t{grid_ok}\t{g['fix']}\n")

    # 汇总打印
    n_read_loci = sum(1 for k in readmode)
    print(f'\n✅ 完成: {len(results)} 样品, {n_read_loci} (样品×位点) 有 reads 实测片段')
    print(f'输出目录: {out_dir}')


if __name__ == '__main__':
    main()
