#!/usr/bin/env python3
"""ssr_common.py — SSR 区域分解与分类（prepare_data.py / count_cn.py 共享）

职责：
  decompose_ssr  把一段序列拆成连续重复分量 [(unit, cn), ...]
  classify_region 把一个 SSR 区域分为「简单」或「复合」，并给出分量与坐标

背景：docx 红色标注是 SSR 区域的权威来源。但区域内部可能是：
  - 简单重复     (PM12 CT×10, PM16 GGAGA×5)         → 简单，用长度法 CN
  - 断裂简单重复 (PM05 AG×13+AA+AG×3，参考有个 AA 插入) → 简单，用长度法 CN（覆盖断裂）
  - 真复合重复   (PM25 CT×9+AC×6，两种基序各 ≥2 份)  → 复合，按分量 run 计数
classify_region 用「显著分量 ≥2」区分复合与简单：显著 = c×len(unit) ≥ max(8, 0.3×区域长)。
"""


def decompose_ssr(region):
    """把 SSR 区拆成分量 [(unit, cn), ...]。简单重复=1 分量；复合（PM25 CT+AC）=多分量。"""
    runs = []
    i = 0
    while i < len(region):
        best = None
        for u in range(2, 7):
            unit = region[i:i + u]
            j, c = i, 0
            while j + u <= len(region) and region[j:j + u] == unit:
                c += 1
                j += u
            if c >= 2 and (best is None or c > best[1]):
                best = (unit, c)
        if best is None:
            runs.append((region[i:], 1))     # 残余非重复序列
            break
        unit, c = best
        runs.append((unit, c))
        i += len(unit) * c
    return runs


def classify_region(region, st, en):
    """把 SSR 区域分类为 'simple' 或 'compound'，返回 (kind, comps)。

    comps = [(unit, cn, abs_start, abs_end), ...]
    - simple：主分量 + 长度法 CN（round(区域长/基序长)），断裂简单（PM05 AG×13+AA+AG×3）归此类。
    - compound：≥2 个显著分量（各 c×len(u) ≥ max(8, 0.3×区域长)），保留分解坐标（PM25 CT+AC）。
    """
    L = en - st
    comps = decompose_ssr(region)
    sub = [(u, c) for u, c in comps if c * len(u) >= max(8, 0.3 * L)]
    if len(sub) >= 2:
        out = []
        pos = st
        for u, c in comps:
            out.append((u, c, pos, pos + len(u) * c))
            pos += len(u) * c
        return 'compound', out
    dom = max(comps, key=lambda uc: uc[1] * len(uc[0]))
    cn = round(L / len(dom[0]))
    return 'simple', [(dom[0], cn, st, en)]
