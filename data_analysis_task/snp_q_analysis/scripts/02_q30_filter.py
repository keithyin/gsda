#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 2a：Q30 过滤。

将 FASTQ 中每个碱基的 baseQ（Phred+33）换算为整条 read 的 readsQ，
保留 readsQ >= min-q 的 reads。

换算方式（与 gseda/ppl/fastq_quality_filter.py::filter_by_avg_quality、
gseda/file_format_cvt/fastq2bam.py::calc_rq 一致）：
    单碱基错误率  p_i = 10^(-baseQ_i/10)
    read 平均错误率 e = mean(p_i)
    readsQ = -10 * log10(e)
即先转碱基正确率/错误率后取平均，再转回 Phred，不能对 baseQ 直接求平均。

用法:
    python3 02_q30_filter.py results/merged/Barcode01.fastq \
        -o results/filtered/Barcode01.q30.fastq --min-q 30
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)


def reads_q_from_base_quals(qual_str: str) -> float:
    """由 FASTQ 质量字符串计算 readsQ（不直接对 baseQ 求平均）。"""
    if not qual_str:
        return 0.0
    err_sum = 0.0
    for c in qual_str:
        err_sum += 10.0 ** (-(ord(c) - 33) / 10.0)
    mean_err = err_sum / len(qual_str)
    if mean_err <= 0.0:
        return float("inf")
    return -10.0 * math.log10(mean_err)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="按 readsQ(由 baseQ 换算) 过滤 FASTQ。",
    )
    parser.add_argument("fastq", help="输入 FASTQ")
    parser.add_argument("-o", "--output", required=True, help="输出 FASTQ 路径")
    parser.add_argument("--min-q", type=float, default=30.0,
                        help="readsQ 阈值，保留 readsQ>=min-q (默认 30)")
    args = parser.parse_args()

    if not os.path.exists(args.fastq):
        log.error("输入文件不存在: %s", args.fastq)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    total = 0
    kept = 0
    kept_bases = 0
    total_bases = 0
    with open(args.fastq, "r") as fin, open(args.output, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            total += 1
            total_bases += len(seq.strip())
            if reads_q_from_base_quals(qual.rstrip("\n")) >= args.min_q:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
                kept += 1
                kept_bases += len(seq.strip())

    log.info("total reads=%d, kept reads=%d (%.1f%%)",
             total, kept, 100.0 * kept / total if total else 0.0)
    log.info("total bases=%d, kept bases=%d (%.1f%%)",
             total_bases, kept_bases, 100.0 * kept_bases / total_bases if total_bases else 0.0)
    log.info("输出: %s", args.output)


if __name__ == "__main__":
    main()
