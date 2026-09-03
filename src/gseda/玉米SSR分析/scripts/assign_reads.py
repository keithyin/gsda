#!/usr/bin/env python3
"""assign_reads.py — 第二步：reads 归属（用真实引物拆数据）

输入：barcodes_reads_fastq_amplicon 下的 PacBio CCS fastq（96 主文件 + 17 复孔），
      每条 read 是 40 位点多重的混合扩增子。
方法：用 data/primers/all_primers.tsv 的 160 条真实使用引物（F1/R1/F2/R2）做双向锚定
      —— 5' 找正向引物、3' 找反向引物 rc，同点 F+rcR 都命中为 both（高置信），
      仅 F 为 fonly（可能截短）；无引物的 read 可选 --ref-fallback 用参考局部比对回收。
产出：data/assign/ 下逐 read 明细、样品×位点汇总、统计、（可选）拆分 fastq。

用法：
  python3 assign_reads.py                          # 默认全量
  python3 assign_reads.py --ref-fallback           # 开启参考回退（回收无引物 read）
  python3 assign_reads.py --write-split            # 同时写拆分 fastq
  python3 assign_reads.py --threads 8
"""
import argparse
import math
import re
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import edlib

# ---------------------------------------------------------------------------
# 全局（worker 初始化时加载）
# ---------------------------------------------------------------------------
_AD = 'GGCCAAAAGGAAGGCCAT'      # 通用 5' 接头；3' 为其反向互补
_AD_RC = _AD.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

# Phred 质量字符 → 误差率 e = 10^(-Q/10)（Q = ord(ch) - 33）。用于读质量计算。
_PHRED_ERR = {chr(c): 10 ** (-(c - 33) / 10.0) for c in range(33, 127)}

_PRIMERS = []      # (locus, type, seq)
_FLIST = []        # (locus, seq)  正向引物 F1/F2
_RLIST = []        # (locus, seq)  反向引物 R1/R2
_REFS = []         # [(locus, seq), ...]
_USE_FALLBACK = False
_FQ_DIR = '.'
_WRITE_SPLIT = False
_MIN_Q = 30.0


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def read_q(qual):
    """碱基质量 → 读质量（PacBio RQ 式概率转换，非 Q 值平均）。

    readsQ = -10·log10(mean(10^(-Q_i/10)))。
    即把每条碱基的误差率求平均，再转回 Phred 刻度。
    注意：直接对 Q 值平均是错的（Jensen 不等式）；也不是整条读无错概率
    （450bp 读下几乎不可能到 Q30）。
    空/异常返回 None（调用方对 None 不过滤，保持兼容）。
    """
    n = len(qual)
    if n == 0:
        return None
    total_err = 0.0
    for ch in qual:
        e = _PHRED_ERR.get(ch)
        if e is None:
            return None          # 非 Phred 字符 → 判不了，不过滤
        total_err += e
    return -10.0 * math.log10(total_err / n)


# ---------------------------------------------------------------------------
# 数据加载
# ---------------------------------------------------------------------------
def load_data(data_dir: Path):
    global _PRIMERS, _FLIST, _RLIST, _REFS
    with open(data_dir / 'primers' / 'all_primers.tsv') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('locus'):
                continue
            p = line.split('\t')
            locus, ptype, seq = p[0], p[3], p[4]
            _PRIMERS.append((locus, ptype, seq))
            if ptype in ('F1', 'F2'):
                _FLIST.append((locus, seq))
            else:
                _RLIST.append((locus, seq))
    with open(data_dir / 'refs' / 'pm_all.fa') as fh:
        cur, seq = None, []
        for line in fh:
            if line.startswith('>'):
                if cur:
                    _REFS.append((cur, ''.join(seq)))
                cur = line[1:].split('|')[0]
                seq = []
            else:
                seq.append(line.strip())
        if cur:
            _REFS.append((cur, ''.join(seq)))


def _worker_init(data_dir, use_fallback, fq_dir, write_split, min_q):
    load_data(Path(data_dir))
    global _USE_FALLBACK, _FQ_DIR, _WRITE_SPLIT, _MIN_Q
    _USE_FALLBACK = use_fallback
    _FQ_DIR = fq_dir
    _WRITE_SPLIT = write_split
    _MIN_Q = min_q


# ---------------------------------------------------------------------------
# 核心：单条 read 归属
# ---------------------------------------------------------------------------
def _find(probe, text, maxd=2):
    """edlib HW：probe 完整匹配在 text 内，返回 (start, end) 或 None。≤maxd 错配。"""
    r = edlib.align(probe, text, mode='HW', task='path')
    if r['editDistance'] is None or r['editDistance'] > maxd:
        return None
    return r['locations'][0]


def trim_adapter(s):
    """修剪通用接头：5' 剪 AD（容错≤2），3' 剪 rc(AD)。reads 两端都带该接头（方向无关）。
    返回 (修剪后序列, 5'修剪数, 3'修剪数)。"""
    n5 = n3 = 0
    if len(s) >= 20:
        r = edlib.align(_AD, s[:20], mode='HW', task='path')
        if r['editDistance'] is not None and r['editDistance'] <= 2:
            st, en = r['locations'][0]
            if st <= 2:
                n5 = en + 1
                s = s[n5:]
    if len(s) >= 20:
        r = edlib.align(_AD_RC, s[-20:], mode='HW', task='path')
        if r['editDistance'] is not None and r['editDistance'] <= 2:
            st, en = r['locations'][0]
            n3 = 20 - st
            if 0 < n3 < len(s) - 10:
                s = s[:-n3]
            else:
                n3 = 0
    return s, n5, n3


def scan(seq):
    """双向引物扫描（一个方向的序列）。返回 locus -> {'F':[(pos,name)], 'R':[(pos,name)]}。"""
    res = {}
    for locus, pseq in _FLIST:
        loc = _find(pseq, seq)
        if loc:
            res.setdefault(locus, {}).setdefault('F', []).append((loc[0], pseq))
    for locus, pseq in _RLIST:
        loc = _find(rc(pseq), seq)          # 反向引物以 rc 形式出现在正链 read 3' 端
        if loc:
            res.setdefault(locus, {}).setdefault('R', []).append((loc[0], pseq))
    return res


def _best_assignment(scanned):
    """从单个方向的 scan 结果里选最优 (locus, conf, f_pos, r_pos, n_loci)，无则 None。"""
    cands = []
    loci = set()
    for locus, dd in scanned.items():
        nf = len(dd.get('F', []))
        nr = len(dd.get('R', []))
        loci.add(locus)
        if nf and nr:
            fpos = min(f[0] for f in dd['F'])
            rpos = min(r[0] for r in dd['R'])
            if fpos < rpos:                    # 保证 F 在 rcR 之前（正确方向）
                cands.append((locus, 'both', fpos, rpos))
        elif nf:
            cands.append((locus, 'fonly', min(f[0] for f in dd['F']), None))
    if not cands:
        return None, len(loci)
    # 优先 both，其次 fonly；同置信级取 F 位置最靠前
    cands.sort(key=lambda c: (0 if c[1] == 'both' else 1, c[2]))
    locus, conf, fpos, rpos = cands[0]
    return (locus, conf, fpos, rpos), len(loci)


def ref_fallback(seq):
    """无引物时对 40 参考×双向做局部比对，返回 (locus, ed, strand) 或 None。"""
    best = None
    for locus, rseq in _REFS:
        for strand, rseq2 in (('+', rseq), ('-', rc(rseq))):
            r = edlib.align(seq, rseq2, mode='HW', task='path')
            if r['editDistance'] is None:
                continue
            st, en = r['locations'][0]
            cov = en - st + 1
            if cov < 0.5 * len(seq):           # 覆盖不足一半则不算
                continue
            key = r['editDistance']
            if best is None or key < best[0]:
                best = (key, locus, strand, cov)
    if best is None:
        return None
    return best[1], best[0], best[2]           # locus, ed, strand


def assign_read(read_seq):
    """返回 (locus, conf, orient, f_pos, r_pos, n_loci, ed)。
    orient: '+' = read 即正向；'-' = read 需取 rc。ed 仅 reffallback 用。"""
    t, _, _ = trim_adapter(read_seq)
    if len(t) < 30:
        return '.', 'none', '+', None, None, 0, None

    best = None
    best_nloci = 0
    for orient, seq in (('+', t), ('-', rc(t))):
        scanned = scan(seq)
        ass, nloci = _best_assignment(scanned)
        best_nloci = max(best_nloci, nloci)
        if ass is None:
            continue
        locus, conf, fpos, rpos = ass
        # 跨方向比较：both > fonly，再比 F 位置
        key = (0 if conf == 'both' else 1, fpos)
        if best is None or key < best[0]:
            best = (key, (locus, conf, orient, fpos, rpos))

    if best:
        _, (locus, conf, orient, fpos, rpos) = best
        return locus, conf, orient, fpos, rpos, best_nloci, None

    # 无引物 → 可选参考回退
    if _USE_FALLBACK:
        fb = ref_fallback(t)
        if fb:
            locus, ed, strand = fb
            return locus, 'reffallback', strand, None, None, 0, ed
    return '.', 'none', '+', None, None, best_nloci, None


# ---------------------------------------------------------------------------
# 单文件处理（worker 内执行）
# ---------------------------------------------------------------------------
def process_file(fname):
    """处理一个 fastq，返回 (sample_key, records, stats, split_dict)。"""
    m = re.match(r'Group_0_Adaptor-barcode(\d+)-(\d+)_(.+?)\.fastq$', fname)
    if not m:
        return None
    barcode, suffix, rest = m.group(1), m.group(2), m.group(3)
    sample = f"barcode{barcode}-{suffix}"
    sample_id = rest if rest != 'None' else ''

    records = []
    split = {}          # locus -> list of (header, seq, qual)
    stats = Counter()   # conf -> n
    reads = 0
    path = Path(_FQ_DIR) / fname
    with open(path) as fh:
        lines = fh.readlines()
    for i in range(0, len(lines), 4):
        if i + 3 >= len(lines):
            break
        header = lines[i].strip()
        seq = lines[i + 1].strip()
        qual = lines[i + 3].strip()          # fastq: @id / seq / + / qual
        reads += 1
        read_id = header[1:] if header.startswith('@') else header

        # Q30 过滤：碱基质量 → 读质量（概率法），低于阈值不进入后续分析
        rq = read_q(qual)
        if rq is not None and rq < _MIN_Q:
            stats['qfiltered'] += 1
            records.append('\t'.join([
                sample, sample_id, read_id, str(len(seq)),
                '.', 'qfiltered', '.', '.', '.', '0', '.', f"{rq:.1f}",
            ]))
            continue

        locus, conf, orient, fpos, rpos, nloci, ed = assign_read(seq)

        if conf != 'none':
            stats[conf] += 1
            # 归一化为正向序列（用于拆分输出）
            fwd_seq = rc(seq) if orient == '-' else seq
            fwd_qual = qual[::-1] if orient == '-' else qual
            fwd_seq, n5, n3 = trim_adapter(fwd_seq)
            fwd_qual = fwd_qual[n5:-n3] if n3 else fwd_qual[n5:]
            if _WRITE_SPLIT:
                split.setdefault(locus, []).append(
                    (f"{read_id}|{sample}|{locus}|{conf}", fwd_seq, fwd_qual))
        else:
            stats['none'] += 1

        records.append('\t'.join([
            sample, sample_id, read_id, str(len(seq)),
            locus, conf, orient,
            str(fpos) if fpos is not None else '.',
            str(rpos) if rpos is not None else '.',
            str(nloci),
            str(ed) if ed is not None else '.',
            f"{rq:.1f}",
        ]))
    return sample, sample_id, reads, records, dict(stats), split


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    global _FQ_DIR, _WRITE_SPLIT
    ap = argparse.ArgumentParser(description='用真实引物把 CCS reads 归属到位点')
    ap.add_argument('--fastq-dir', default=(
        '/data1/ccs_data/20260813-corn-ssr-analysis/'
        '20260702_251201Y0001_Run0003/Group_0/barcodes_reads_fastq_amplicon'))
    ap.add_argument('--data-dir', default=str(
        Path(__file__).resolve().parent.parent / 'data'))
    ap.add_argument('--out', default=None,
                    help='输出目录（默认 <data-dir>/assign）')
    ap.add_argument('--ref-fallback', action='store_true',
                    help='对无引物 read 用参考局部比对回收（默认关）')
    ap.add_argument('--write-split', action='store_true',
                    help='写拆分 fastq（每样品每位置一个文件）')
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--min-q', type=float, default=30.0,
                    help='读质量阈值：碱基质量→读质量(概率法)低于此值的 read 被过滤（默认30）')
    ap.add_argument('--only', default=None,
                    help='只处理文件名包含该串的文件（调试用）')
    args = ap.parse_args()

    _FQ_DIR = args.fastq_dir
    _WRITE_SPLIT = args.write_split
    _MIN_Q = args.min_q
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out) if args.out else data_dir / 'assign'

    fq_files = sorted(Path(_FQ_DIR).glob('*.fastq'))
    if args.only:
        fq_files = [f for f in fq_files if args.only in f.name]
    if not fq_files:
        print(f'⚠️  未找到 fastq 文件: {_FQ_DIR}', file=sys.stderr)
        sys.exit(1)
    print(f'待处理文件: {len(fq_files)} 个（读质量阈值 Q{args.min_q:g}）')

    threads = args.threads or 4
    with Pool(threads, initializer=_worker_init,
              initargs=(str(data_dir), args.ref_fallback, args.fastq_dir,
                        args.write_split, args.min_q)) as pool:
        results = pool.map(process_file, [f.name for f in fq_files], chunksize=1)

    results = [r for r in results if r]

    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- read_assignments.tsv ----
    with open(out_dir / 'read_assignments.tsv', 'w') as fh:
        fh.write('sample\tsample_id\tread_id\tread_len\tlocus\tconf\torient\tf_pos\tr_pos\tn_loci\ted\treadsQ\n')
        for sample, sample_id, reads, records, stats, split in results:
            for rec in records:
                fh.write(rec + '\n')

    # ---- summary.tsv（样品×位点）----
    summ = {}
    for sample, sample_id, reads, records, stats, split in results:
        summ[sample] = {'sample_id': sample_id, 'reads': reads, 'loci': Counter()}
    # 重新统计每个 read 记录
    with open(out_dir / 'read_assignments.tsv') as fh:
        next(fh)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            sample, locus, conf = p[0], p[4], p[5]
            if conf not in ('none', 'qfiltered'):
                summ[sample]['loci'][(locus, conf)] += 1
    with open(out_dir / 'summary.tsv', 'w') as fh:
        fh.write('sample\tsample_id\ttotal_reads\tlocus\tboth\tfonly\treffallback\n')
        for sample in sorted(summ, key=lambda s: int(re.search(r'(\d+)', s).group(1))):
            s = summ[sample]
            for locus in sorted({k[0] for k in s['loci']}, key=lambda x: int(x[2:])):
                both = s['loci'].get((locus, 'both'), 0)
                fo = s['loci'].get((locus, 'fonly'), 0)
                fb = s['loci'].get((locus, 'reffallback'), 0)
                if not (both + fo + fb):
                    continue
                fh.write(f"{sample}\t{s['sample_id']}\t{s['reads']}\t{locus}\t{both}\t{fo}\t{fb}\n")

    # ---- stats.tsv（每文件归属统计）----
    with open(out_dir / 'stats.tsv', 'w') as fh:
        fh.write('sample\tsample_id\ttotal_reads\tboth\tfonly\tnone\tqfiltered\tassigned\tassigned_pct\tmulti_loci\n')
        for sample, sample_id, reads, records, stats, split in results:
            b = stats.get('both', 0)
            fo = stats.get('fonly', 0)
            no = stats.get('none', 0)
            qf = stats.get('qfiltered', 0)
            assigned = b + fo
            pct = assigned / reads * 100 if reads else 0
            multi = sum(1 for rec in records if int(rec.split('\t')[9]) > 1)
            fh.write(f"{sample}\t{sample_id}\t{reads}\t{b}\t{fo}\t{no}\t{qf}"
                     f"\t{assigned}\t{pct:.1f}\t{multi}\n")

    # ---- split fastq（可选）----
    if args.write_split:
        spdir = out_dir / 'split'
        spdir.mkdir(exist_ok=True)
        for sample, sample_id, reads, records, stats, split in results:
            for locus, recs in split.items():
                sub = spdir / sample
                sub.mkdir(exist_ok=True)
                with open(sub / f"{locus}.fastq", 'w') as fh:
                    for header, seq, qual in recs:
                        fh.write(f"@{header}\n{seq}\n+\n{qual}\n")

    # ---- 汇总打印 ----
    tot_reads = sum(r[2] for r in results)
    tot_both = sum(r[4].get('both', 0) for r in results)
    tot_fo = sum(r[4].get('fonly', 0) for r in results)
    tot_none = sum(r[4].get('none', 0) for r in results)
    tot_qf = sum(r[4].get('qfiltered', 0) for r in results)
    print(f'\n✅ 完成: {len(results)} 文件, {tot_reads} reads')
    print(f'   both   {tot_both} ({tot_both / tot_reads * 100:.1f}%)')
    print(f'   fonly  {tot_fo} ({tot_fo / tot_reads * 100:.1f}%)')
    print(f'   合计   {tot_both + tot_fo} ({(tot_both + tot_fo) / tot_reads * 100:.1f}%)')
    print(f'   none   {tot_none} ({tot_none / tot_reads * 100:.1f}%)')
    print(f'   Q过滤   {tot_qf} ({tot_qf / tot_reads * 100:.1f}%)')
    print(f'\n输出目录: {out_dir}')


if __name__ == '__main__':
    main()
