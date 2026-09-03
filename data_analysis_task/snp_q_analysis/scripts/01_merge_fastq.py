#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 1：合并两次运行 (run1/run2) 中相同 barcode 的 FASTQ 文件。

将 run1/BarcodeXX.fastq 与 run2/BarcodeXX.fastq 按顺序拼接为单个 FASTQ。
流式处理，不将全部序列载入内存；并校验两个 run 的 read 名无重复。

用法:
    python3 01_merge_fastq.py Barcode01 \
        --run1-dir /data1/ccs_data/20260805-xunming/run1 \
        --run2-dir /data1/ccs_data/20260805-xunming/run2 \
        -o results/merged/Barcode01.fastq
"""

from __future__ import annotations

import argparse
import logging
import os
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)


def _count_and_copy(path: str, out_fh, name_seen: set[str] | None) -> int:
    """将单个 FASTQ 逐行拷入 out_fh，返回 read 数。

    name_seen: 若不为 None，收集该文件中的 read 名（用于后续去重校验）。
    """
    count = 0
    with open(path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if name_seen is not None:
                name_seen.add(header[1:].split()[0])
            out_fh.write(header)
            out_fh.write(seq)
            out_fh.write(plus)
            out_fh.write(qual)
            count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="合并 run1/run2 中相同 barcode 的 FASTQ 文件。",
    )
    parser.add_argument("barcode", help="barcode 名，如 Barcode01")
    parser.add_argument("--run1-dir", required=True, help="run1 目录")
    parser.add_argument("--run2-dir", required=True, help="run2 目录")
    parser.add_argument("-o", "--output", required=True, help="输出 FASTQ 路径")
    args = parser.parse_args()

    run1_path = os.path.join(args.run1_dir, f"{args.barcode}.fastq")
    run2_path = os.path.join(args.run2_dir, f"{args.barcode}.fastq")
    for p in (run1_path, run2_path):
        if not os.path.exists(p):
            log.error("输入文件不存在: %s", p)
            sys.exit(1)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    run1_names: set[str] = set()
    with open(args.output, "w") as out_fh:
        n1 = _count_and_copy(run1_path, out_fh, run1_names)
        n2 = _count_and_copy(run2_path, out_fh, None)

    # 校验 run2 与 run1 无同名 read（若重复需二次扫描 run2 才能检测，
    # 这里对 run2 再扫一遍只取其名字，O(N) 一次）
    dup = 0
    if run1_names:
        with open(run2_path, "r") as f:
            while True:
                header = f.readline()
                if not header:
                    break
                f.readline()
                f.readline()
                f.readline()
                if header[1:].split()[0] in run1_names:
                    dup += 1

    log.info("%s: run1=%d reads, run2=%d reads, merged=%d reads", args.barcode, n1, n2, n1 + n2)
    if dup:
        log.warning("%s: run1/run2 存在 %d 个同名 read", args.barcode, dup)
    else:
        log.info("%s: run1/run2 无同名 read", args.barcode)
    log.info("输出: %s", args.output)


if __name__ == "__main__":
    main()
