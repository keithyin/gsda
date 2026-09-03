#!/usr/bin/env python3
"""count_cn.py — 第三步：copy number 计数（reads 一致性口径）

输入：data/assign/split/ 下第二步拆出的逐位点 fastq（正向归一、接头已剪、header 含 read_id|样品|位点|conf）。
方法：SSR 靶点（ssr_targets.bed，docx 红色标注权威）→ 参考取 22bp 侧翼锚 → read 中 edlib 找左翼（≤2 错配）
      → SSR 起点 → 窗口内数最长纯重复运行（基序全部旋转）；复合重复自动分解分量。
输出：data/cn/ 下 cn_summary.tsv（样品×位点×分量：主CN/P/s1/s2/次峰/判定）、
      cn_per_read.tsv（逐 read 明细）、cn_summary_stats.tsv（每样品汇总）。

用法：
  python3 count_cn.py                 # 全量
  python3 count_cn.py --threads 16
  python3 count_cn.py --only barcode201-0   # 只处理某样品
"""
import argparse
import re
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import edlib

from ssr_common import classify_region

# ---------------------------------------------------------------------------
# 全局（worker 初始化时加载）
# ---------------------------------------------------------------------------
_REFSEQ = {}        # locus -> 参考序列
_LOCUS = {}         # locus -> [region, ...]，每个 region 一个 SSR 靶点（含复合分量）
_MIN_DEPTH = 10     # 每分量最低计数深度（depth<此值的分量不写入输出）


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


# ---------------------------------------------------------------------------
# 数据加载
# ---------------------------------------------------------------------------
def load_data(data_dir: Path):
    global _REFSEQ, _LOCUS
    with open(data_dir / 'refs' / 'pm_all.fa') as fh:
        cur = None
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

    with open(data_dir / 'annot' / 'ssr_targets.bed') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            locus, st, en = p[0], int(p[2]), int(p[3])
            target = int(p[6]) if len(p) > 6 else 1
            ref = _REFSEQ[locus]
            L = 22
            # classify_region：复合（PM25 CT+AC）保留分量；断裂简单（PM05）归简单 + 长度法 CN
            kind, comps = classify_region(ref[st:en], st, en)
            _LOCUS.setdefault(locus, []).append({
                'target': target, 'start': st, 'end': en,
                'left': ref[st - L:st],
                'right': ref[en:en + L],
                'ref_ssr_len': en - st,
                'kind': kind,
                'comps': [(u, c) for u, c, _, _ in comps],
            })


def _worker_init(data_dir, only, min_depth):
    load_data(Path(data_dir))
    global _ONLY, _MIN_DEPTH
    _ONLY = only
    _MIN_DEPTH = min_depth


# ---------------------------------------------------------------------------
# 计数核心
# ---------------------------------------------------------------------------
def _hw(probe, text, maxd=2):
    r = edlib.align(probe, text, mode='HW', task='path')
    if r['editDistance'] is None or r['editDistance'] > maxd:
        return None
    return r['locations'][0]


def count_run(seg, motif):
    """seg 中最长的连续 motif 拷贝数（基序全部旋转）。O(n)。"""
    m = len(motif)
    best = 0
    for rot in range(m):
        mm = (motif + motif)[rot:rot + m]
        c = 0
        i = 0
        while i <= len(seg) - m:
            if seg[i:i + m] == mm:
                c += 1
                i += m
                if c > best:
                    best = c
            else:
                c = 0
                i += 1
    return best


def count_read_cn(read, region):
    """返回 (len_cns, run_cns, found) 或 (None, None, False)。
    len_cns：左右翼都有时侧翼间长度法 CN（每分量；简单位点才有效，鲁棒于内部错误）。
    run_cns：左翼窗口内最长纯运行 CN（每分量）。
    found=False：左翼未命中 → 无法可靠定位（多为伪影/侧翼差异），调用方应排除。"""
    loc = _hw(region['left'], read)
    if loc is None:
        return None, None, False
    L = loc[1] + 1
    win = read[L:min(len(read), L + region['ref_ssr_len'] + 150)]
    run = [count_run(win, unit) for unit, _ in region['comps']]
    rr = _hw(region['right'], read)
    if rr is not None and rr[0] > L:
        R = rr[0]
        lens = [round((R - L) / len(unit)) for unit, _ in region['comps']]
        return lens, run, True
    return None, run, True


# ---------------------------------------------------------------------------
# 单样品处理（worker 内执行）
# ---------------------------------------------------------------------------
def process_sample(sample):
    """处理一个样品目录，返回 (per_read, summaries, sample_stats)。"""
    if _ONLY and sample != _ONLY:
        return None
    sdir = Path(_SPLIT_DIR) / sample
    if not sdir.is_dir():
        return None
    per_read = []          # 逐 read 行
    summaries = []         # 样品×位点×分量 汇总行
    n_loci_ge5 = 0
    p_sums = []
    tot_depth = 0

    for fq in sorted(sdir.glob('*.fastq')):
        locus = fq.stem
        regions = _LOCUS.get(locus)
        if not regions:
            continue
        lines = fq.read_text().splitlines()

        # 每个 SSR 靶点区域分别计数
        for region in regions:
            recs = []            # (read_id, conf, len_cns, run_cns)
            no_flank = 0
            for i in range(0, len(lines), 4):
                header = lines[i].strip()
                seq = lines[i + 1].strip()
                parts = header.lstrip('@').split('|')
                read_id = parts[0]
                conf = parts[3] if len(parts) >= 4 else ''
                lens, run, found = count_read_cn(seq, region)
                if not found:
                    no_flank += 1        # 伪影/侧翼差异，不强行计数
                    continue
                recs.append((read_id, conf, lens, run))

            if not recs:
                continue

            # 逐分量
            for ci, (unit, refcn) in enumerate(region['comps']):
                # 方法选择：真复合→run（PM25 型，按分量）；简单/断裂简单→len 为主
                # （PM05 型断裂重复用长度法，鲁棒于内部 AA 插入），len/run 众数差>2 且 run≈ref_cn 才切 run
                if region['kind'] == 'compound':
                    method = 'run'
                else:
                    len_vals = [r[2][ci] for r in recs if r[2] and r[2][ci] is not None]
                    run_vals = [r[3][ci] for r in recs]
                    len_mode = Counter(len_vals).most_common(1)[0][0] if len_vals else None
                    run_mode = Counter(run_vals).most_common(1)[0][0] if run_vals else None
                    if (len_mode is not None and run_mode is not None
                            and abs(len_mode - run_mode) > 2
                            and abs(run_mode - refcn) <= 1):
                        method = 'run'
                    else:
                        method = 'len'

                dist = Counter()
                both_n = 0
                for read_id, conf, lens, run in recs:
                    if method == 'len' and lens and lens[ci] is not None:
                        cn = lens[ci]
                    else:
                        cn = run[ci]
                    dist[cn] += 1
                    if conf == 'both':
                        both_n += 1
                depth = sum(dist.values())
                if depth < _MIN_DEPTH:
                    # 深度不足的分量不写入输出（summary + per_read 都跳过）
                    continue
                mode_cn, mode_n = dist.most_common(1)[0]
                p = mode_n / depth
                s1 = (dist.get(mode_cn + 1, 0) + dist.get(mode_cn - 1, 0)) / depth
                s2 = 1 - p - s1
                second = dist.most_common(2)
                sec_cn = second[1][0] if len(second) > 1 else '.'
                sec_pct = second[1][1] / depth * 100 if len(second) > 1 else 0
                fonly_n = depth - both_n
                verdict = '✅' if p >= 0.6 else ('⚠️' if p >= 0.4 else '❌')
                dist_str = ','.join(f"{k}:{v}" for k, v in sorted(dist.items()))
                flag = f'no_flank{no_flank}' if no_flank else ''
                ssr_pos = f"{region['start']}-{region['end']}"
                summaries.append(
                    f"{sample}\t{locus}\t{region['target']}\t{ssr_pos}\t{unit}\t{refcn}\t"
                    f"{method}\t{depth}\t{both_n}\t{fonly_n}\t"
                    f"{mode_cn}\t{p * 100:.1f}\t{s1 * 100:.1f}\t{s2 * 100:.1f}\t"
                    f"{sec_cn}\t{sec_pct:.1f}\t{flag}\t{verdict}\t{dist_str}")
                if depth >= 5:
                    n_loci_ge5 += 1
                    p_sums.append(p)
                    tot_depth += depth

                for read_id, conf, lens, run in recs:
                    cn = (lens[ci] if (method == 'len' and lens and lens[ci] is not None) else run[ci])
                    per_read.append(f"{sample}\t{locus}\t{region['target']}\t{unit}\t{cn}\t{conf}")

    mean_p = f"{sum(p_sums) / len(p_sums) * 100:.1f}" if p_sums else '.'
    sample_stats = f"{sample}\t{n_loci_ge5}\t{mean_p}\t{tot_depth}"
    return per_read, summaries, sample_stats


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    global _SPLIT_DIR
    ap = argparse.ArgumentParser(description='数 SSR copy number（reads 一致性口径）')
    ap.add_argument('--split-dir', default=None,
                    help='第二步拆分输出目录（默认 <data-dir>/assign/split）')
    ap.add_argument('--data-dir', default=str(
        Path(__file__).resolve().parent.parent / 'data'))
    ap.add_argument('--out', default=None,
                    help='输出目录（默认 <data-dir>/cn）')
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--only', default=None, help='只处理某样品（调试用）')
    ap.add_argument('--min-depth', type=int, default=10,
                    help='每分量最低计数深度，低于此值的分量不输出（默认10）')
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    _SPLIT_DIR = Path(args.split_dir) if args.split_dir else data_dir / 'assign' / 'split'
    out_dir = Path(args.out) if args.out else data_dir / 'cn'

    samples = sorted(d.name for d in _SPLIT_DIR.iterdir() if d.is_dir())
    if args.only:
        samples = [s for s in samples if args.only in s]
    if not samples:
        print(f'⚠️  未找到样品目录: {_SPLIT_DIR}', file=__import__('sys').stderr)
        return
    print(f'待处理样品: {len(samples)} 个（{_SPLIT_DIR}，min-depth={args.min_depth}）')

    threads = args.threads or 4
    with Pool(threads, initializer=_worker_init,
              initargs=(str(data_dir), args.only, args.min_depth)) as pool:
        results = pool.map(process_sample, samples, chunksize=1)
    results = [r for r in results if r]

    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / 'cn_per_read.tsv', 'w') as fh:
        fh.write('sample\tlocus\ttarget\tunit\tcn\tconf\n')
        for per_read, summaries, ss in results:
            for r in per_read:
                fh.write(r + '\n')

    with open(out_dir / 'cn_summary.tsv', 'w') as fh:
        fh.write('sample\tlocus\ttarget\tssr_pos\tunit\tref_cn\tmethod\tdepth\tdepth_both\tdepth_fonly\t'
                 'mode_cn\tP_pct\ts1_pct\ts2_pct\tsecond_cn\tsecond_pct\t'
                 'flag\tverdict\tcn_dist\n')
        for per_read, summaries, ss in results:
            for s in summaries:
                fh.write(s + '\n')

    with open(out_dir / 'cn_summary_stats.tsv', 'w') as fh:
        fh.write('sample\tn_loci_ge5\tmean_P_pct\ttot_depth_ge5\n')
        for per_read, summaries, ss in results:
            fh.write(ss + '\n')

    # 汇总打印
    n_sum = sum(len(r[1]) for r in results)
    n_read = sum(len(r[0]) for r in results)
    print(f'\n✅ 完成: {len(results)} 样品, {n_sum} 行(样品×位点×分量), {n_read} 条 read 计数')
    print(f'输出目录: {out_dir}')


if __name__ == '__main__':
    main()
