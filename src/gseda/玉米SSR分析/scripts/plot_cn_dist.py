#!/usr/bin/env python3
"""plot_cn_dist.py — CN 分布图（按位点拆页，样品×位点×基序 粒度）

输入：data/cn/cn_summary.tsv（每行 = 样品×位点×target×基序 的 CN 分布 cn_dist）。
输出：data/cn/plots/cn_dist_<locus>[_pK].jpeg —— 每个位点一/多页，页内网格每个子图 =
      该位点一个 (样品, 基序) 的 CN 柱状分布。

子图内容：
  x = CN（重复单元数，每子图自适应范围），y = read 数
  柱色按距主峰距离：|CN - mode_cn|：0 深蓝（主峰）、1 浅蓝、≥2 灰
  红色虚线 = ref_cn（参考重复数）
  标题 = 样品 基序(ref_cn) m主CN P% d深度，颜色按判定 ✅绿/⚠️橙/❌红

用法：
  python3 plot_cn_dist.py                          # 全量
  python3 plot_cn_dist.py --only-locus PM05 PM25   # 只画指定位点
  python3 plot_cn_dist.py --only-sample barcode201-0   # 只画某样品的子图（页仍按位点）
  python3 plot_cn_dist.py --min-depth 5            # 只画深度≥5 的子图
  python3 plot_cn_dist.py --max-subplots 45        # 每页最多子图数（大位点自动拆页）
  python3 plot_cn_dist.py --ncols 8 --dpi 150      # 布局/分辨率
"""
import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

VERDICT_COLOR = {'✅': 'green', '⚠️': 'orange', '❌': 'red'}


def parse_cn_dist(s):
    """'5:1,6:1,15:4,16:30,17:64,18:23' -> [(5,1),(6,1),...]"""
    out = []
    for pair in s.split(','):
        pair = pair.strip()
        if not pair:
            continue
        cn, _, n = pair.partition(':')
        out.append((int(cn), int(n)))
    return out


def _unit_xlim(subplots, pad=3.0):
    """同基序组的共享 x 范围（整数边界，调用方 ±0.5 留白）。

    基于各组子图 mode 的 IQR 稳健窗口（Tukey fence）：
      1. fence = [q1-1.5·iqr, q3+1.5·iqr]，落在 fence 内的 mode 视为"主体"；
      2. 窗口 = 主体 mode 的 [min, max] ± pad（避免 fence 对宽等位基因位点过度外扩）。
      - 紧簇位点（PM05 全在 13，3 个离群 24/26）→ 窗口 [mode±pad]，离群被排除。
      - 真实宽等位基因位点（PM22 mode 10~40）→ 窗口覆盖全部真实 mode。
      - 双峰真实等位基因（一半 6 一半 10）→ fence 覆盖两者，均保留。
    """
    modes = sorted(int(row['mode_cn']) for row, _ in subplots)
    n = len(modes)
    q1 = modes[n // 4]
    q3 = modes[(3 * n) // 4]
    iqr = q3 - q1
    lo = q1 - 1.5 * iqr
    hi = q3 + 1.5 * iqr
    bulk = [m for m in modes if lo <= m <= hi]
    xlo = max(1, min(bulk) - pad)
    xhi = max(bulk) + pad
    return xlo, xhi


def plot_page(locus, rows, out_path, page_label, ncols, dpi):
    """把一批 (样品, 基序) CN 分布画成一页。rows 已是过滤排序后的。"""
    n = len(rows)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.2, nrows * 2.4))
    axes = axes.flatten() if nrows * ncols > 1 else [axes]

    # 同一图内 x 轴一致：按基序组共享稳健范围（多基序位点各基序各自一致，
    # 避免 PM40 的 CA(≈4) 被 AG(≈29) 拉没；单基序位点即整页一致）
    by_unit = defaultdict(list)
    for row, dist in rows:
        by_unit[row['unit']].append((row, dist))
    unit_xlims = {u: _unit_xlim(subs) for u, subs in by_unit.items()}

    # y 轴统一用比例（0~1）：每个子图按自身总 reads 归一化，深度差异不影响读图
    for ax, (row, dist) in zip(axes, rows):
        cns = [c for c, _ in dist]
        total = sum(n for _, n in dist)
        props = [n / total for _, n in dist]
        mode = int(row['mode_cn'])
        ref = int(row['ref_cn'])
        # 柱色：主峰深蓝、±1 浅蓝、更远灰
        colors = ['#4C72B0' if c == mode else '#A6C3E5' if abs(c - mode) == 1
                  else '#C9C9C9' for c in cns]
        ax.bar(cns, props, width=0.8, color=colors)
        ax.axvline(ref, color='red', ls='--', lw=1.2)
        # x 轴：同基序组一致（IQR 稳健窗口）；刻度强制整数（CN 是整数）
        xlo, xhi = unit_xlims[row['unit']]
        ax.set_xlim(xlo - 0.5, xhi + 0.5)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylim(0, 1.0)
        ax.tick_params(labelsize=8)
        ax.set_xlabel('CN', fontsize=8)
        ax.set_ylabel('proportion', fontsize=8)
        tc = VERDICT_COLOR.get(row['verdict'], 'black')
        ax.set_title(
            f"{row['sample']} {row['unit']}(ref{ref}) m{mode} "
            f"P{float(row['P_pct']):.0f}% d{row['depth']}",
            fontsize=9, color=tc)

    for ax in axes[n:]:
        ax.axis('off')
    title = f"{locus} — CN distribution (each subplot = sample × motif; "
    if page_label:
        title += f"  {page_label}; "
    title += "red dashed = ref_cn; dark blue = mode)"
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='CN 分布图（按位点拆页）')
    ap.add_argument('--cn-dir', default=None, help='cn_summary 所在目录（默认 <data>/cn）')
    ap.add_argument('--out', default=None, help='输出目录（默认 <cn-dir>/plots）')
    ap.add_argument('--only-locus', nargs='*', help='只画指定位点（如 PM05 PM25）')
    ap.add_argument('--only-sample', nargs='*', help='只画指定样品的子图')
    ap.add_argument('--min-depth', type=int, default=0, help='只画深度≥此值的子图')
    ap.add_argument('--max-subplots', type=int, default=45,
                    help='每页最多子图数（超出按排序拆页；0=不拆）')
    ap.add_argument('--ncols', type=int, default=8, help='每页列数')
    ap.add_argument('--dpi', type=int, default=150)
    args = ap.parse_args()

    cn_dir = Path(args.cn_dir) if args.cn_dir else \
        Path(__file__).resolve().parent.parent / 'data' / 'cn'
    out_dir = Path(args.out) if args.out else cn_dir / 'plots'

    rows = []
    with open(cn_dir / 'cn_summary.tsv') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < len(header):
                continue
            rows.append(dict(zip(header, p)))

    # 过滤
    only_l = set(args.only_locus) if args.only_locus else None
    only_s = set(args.only_sample) if args.only_sample else None
    kept = []
    for r in rows:
        if only_l and r['locus'] not in only_l:
            continue
        if only_s and r['sample'] not in only_s:
            continue
        if args.min_depth and int(r['depth']) < args.min_depth:
            continue
        dist = parse_cn_dist(r['cn_dist'])
        if not dist:
            continue
        kept.append((r, dist))

    by_locus = defaultdict(list)
    for r, dist in kept:
        by_locus[r['locus']].append((r, dist))

    print(f'共 {len(by_locus)} 个位点、{len(kept)} 个子图（样品×位点×基序）')
    n_page = 0
    for locus in sorted(by_locus, key=lambda x: int(x[2:])):
        loc_rows = by_locus[locus]
        # 排序：按基序聚在一起（多靶点 CT/AC/TA 分组），再按主CN
        loc_rows.sort(key=lambda rd: (rd[0]['unit'], int(rd[0]['mode_cn'])))
        # 拆页：均衡分页，避免最后一页只有几个子图（0=不拆）
        if args.max_subplots and len(loc_rows) > args.max_subplots:
            n_pg = (len(loc_rows) + args.max_subplots - 1) // args.max_subplots
            chunk_size = (len(loc_rows) + n_pg - 1) // n_pg
            chunks = [loc_rows[i:i + chunk_size]
                      for i in range(0, len(loc_rows), chunk_size)]
        else:
            chunks = [loc_rows]
        n_page_total = len(chunks)
        for k, chunk in enumerate(chunks, 1):
            page_label = f'page {k}/{n_page_total}' if n_page_total > 1 else ''
            fname = f'cn_dist_{locus}.jpeg' if n_page_total == 1 \
                else f'cn_dist_{locus}_p{k}.jpeg'
            plot_page(locus, chunk, out_dir / fname, page_label, args.ncols, args.dpi)
            n_page += 1

    print(f'✅ 生成 {n_page} 页 jpeg → {out_dir}')
    print('  用 --only-locus/--only-sample/--min-depth 过滤；'
          '--max-subplots/--ncols/--dpi 调布局')


if __name__ == '__main__':
    main()
