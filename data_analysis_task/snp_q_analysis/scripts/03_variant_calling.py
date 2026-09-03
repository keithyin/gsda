#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 2b：将 Q30 过滤后的 reads 比对到参考序列，输出 SNP / Insertion / Deletion 位点。

对每条 read 取 primary 比对，要求 full-length（qcov=1.0 且 rcov=1.0）且
identity >= min-identity，然后按 CIGAR 提取所有变体事件：
    - SNP:    碱基错配，baseQ = read 碱基质量
    - INS:    read 相对参考的插入，baseQ = 插入碱基的最小 Q
    - DEL:    read 相对参考的删除，baseQ = 删除位点两侧紧邻 read 碱基的最小 Q
              （即离被删参考片段最近的两个 read 碱基的 min Q；仅单侧存在时用单侧）

比对与提取逻辑参考 gseda/16s/align_mismatch_dist_and_q_dist.py；
负链 read 使用 gseda/align_ana/mappy_ext.py 的 revcomp / parse_cigar。

输出 TSV 列（前 8 列与旧 SNP 版一致，末尾追加 variant_type）：
    barcode  read_name  ref_pos(1-based)  ref_base  alt_base  base_q  strand  read_len  variant_type

ref_pos 约定（统一 rp_0based + 1）：
    SNP: 变体位点；INS: 插入片段后第一个参考碱基；DEL: 第一个被删参考碱基
ref_base/alt_base：SNP=参考碱基/read 碱基；INS=-/插入序列；DEL=被删参考序列/-

用法:
    python3 03_variant_calling.py --fastq results/filtered/Barcode01.q30.fastq \
        --ref /data1/ccs_data/20260805-xunming/Barcode01_300_draft.fasta \
        --barcode Barcode01 -o results/variant/Barcode01.variants.tsv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys

import mappy
from tqdm import tqdm

from gseda.align_ana.mappy_ext import revcomp, parse_cigar

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)

# 与 gseda/16s/align_mismatch_dist_and_q_dist.py 一致的成熟比对配置
ALIGNER_EXTRA_FLAGS = 67108864
ALIGNER_K = 11
ALIGNER_W = 1
ALIGNER_BEST_N = 10


def read_reference(ref_path: str) -> tuple[str, str]:
    """读取参考 FASTA，返回 (名称, 序列)。要求只含一条序列。"""
    name = None
    seq_parts = []
    with open(ref_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    raise ValueError("参考 FASTA 需仅包含一条序列")
                name = line[1:].split()[0]
            elif name is not None:
                seq_parts.append(line)
    if name is None:
        raise ValueError("参考 FASTA 为空")
    return name, "".join(seq_parts)


def calculate_identity_from_cigar(cigar_string: str) -> float:
    """根据 CIGAR 计算 identity（同 gseda/align_ana/align_mismatch_dist.py）。"""
    matches = re.findall(r"(\d+)([=IDX])", cigar_string)
    match_count = 0
    total_aligned = 0
    for length_str, operation in matches:
        length = int(length_str)
        if operation == "=":
            match_count += length
            total_aligned += length
        elif operation == "X":
            total_aligned += length
        elif operation in ("I", "D"):
            total_aligned += length
    return match_count / total_aligned if total_aligned else 0.0


def extract_variants(query: str, qual: list[int], ref: str, hit) -> list[tuple]:
    """提取比对中所有变体事件（SNP / INS / DEL）。

    Returns:
        list of (variant_type, ref_pos_1based, ref_base, alt_base, baseq)

    链方向与坐标处理：
        hit.q_st/q_en 基于原始 query 序列坐标。反向比对（strand=-1）时，
        先将 query revcomp、qual 反转，再按 mappy_ext 的约定重映射坐标：
            q_start = len(qseq) - hit.q_en
        CIGAR 在 revcomp 后的序列上遍历。

    baseq 规则：
        SNP: 该位点 read 碱基质量
        INS: 插入碱基的最小 Q
        DEL: 删除位点两侧紧邻 read 碱基的最小 Q（仅单侧存在时用单侧）
    """
    qseq = query
    qqual = qual
    q_start, q_end = hit.q_st, hit.q_en
    if hit.strand == -1:
        # 原始坐标 q_st/q_en 映射到 revcomp 后新序列上
        qseq = revcomp(query)
        qqual = qual[::-1]
        q_start = len(qseq) - hit.q_en
        q_end = len(qseq) - hit.q_st

    q_pos = q_start
    r_pos = hit.r_st
    variants = []
    for length, op in parse_cigar(hit.cigar_str):
        if op in ("=", "X"):
            for _ in range(length):
                qb = qseq[q_pos]
                rb = ref[r_pos]
                if qb != rb:
                    variants.append(("SNP", r_pos + 1, rb, qb, qqual[q_pos]))
                q_pos += 1
                r_pos += 1
        elif op == "I":
            ins_seq = qseq[q_pos:q_pos + length]
            base_q = min(qqual[q_pos:q_pos + length])
            variants.append(("INS", r_pos + 1, "-", ins_seq, base_q))
            q_pos += length
        elif op == "D":
            del_seq = ref[r_pos:r_pos + length]
            flank_qs = []
            if q_pos - 1 >= 0:
                flank_qs.append(qqual[q_pos - 1])
            if q_pos < len(qqual):
                flank_qs.append(qqual[q_pos])
            base_q = min(flank_qs) if flank_qs else 0
            variants.append(("DEL", r_pos + 1, del_seq, "-", base_q))
            r_pos += length
    return variants


def main() -> None:
    parser = argparse.ArgumentParser(description="比对参考并输出 SNP/INS/DEL 位点。")
    parser.add_argument("--fastq", required=True, help="Q30 过滤后的 FASTQ")
    parser.add_argument("--ref", required=True, help="参考 FASTA (单条序列)")
    parser.add_argument("--barcode", required=True, help="barcode 名，如 Barcode01")
    parser.add_argument("-o", "--output", required=True, help="输出变体 TSV 路径")
    parser.add_argument("--min-identity", type=float, default=0.9,
                        help="最低 identity (默认 0.9)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="比对线程数")
    args = parser.parse_args()

    for p in (args.fastq, args.ref):
        if not os.path.exists(p):
            log.error("输入文件不存在: %s", p)
            sys.exit(1)

    ref_name, ref_seq = read_reference(args.ref)
    log.info("参考: %s, 长度: %d", ref_name, len(ref_seq))

    aligner = mappy.Aligner(
        seq=ref_seq,
        extra_flags=ALIGNER_EXTRA_FLAGS,
        k=ALIGNER_K,
        w=ALIGNER_W,
        best_n=ALIGNER_BEST_N,
        n_threads=args.threads,
    )
    if aligner is None:
        log.error("mappy 无法构建 aligner")
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    total = 0
    n_aligned = 0
    n_full = 0
    n_var_reads = 0
    n_events = {"SNP": 0, "INS": 0, "DEL": 0}

    with open(args.fastq, "r") as fin, open(args.output, "w") as fout:
        fout.write("barcode\tread_name\tref_pos\tref_base\talt_base\tbase_q\tstrand\tread_len\tvariant_type\n")
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().strip()
            fin.readline()  # +
            qual_str = fin.readline().strip()
            total += 1
            read_name = header[1:].split()[0]

            read_primary = False
            for hit in aligner.map(seq):
                if not hit.is_primary:
                    continue
                read_primary = True

                qcov = (hit.q_en - hit.q_st) / len(seq)
                rcov = (hit.r_en - hit.r_st) / len(ref_seq)
                identity = calculate_identity_from_cigar(hit.cigar_str)
                if qcov != 1.0 or rcov != 1.0 or identity < args.min_identity:
                    continue  # 试下一个 primary；仍非 full-length 则该 read 不产出变体
                n_full += 1

                strand = hit.strand
                qual = [ord(c) - 33 for c in qual_str]

                variants = extract_variants(seq, qual, ref_seq, hit)
                if variants:
                    n_var_reads += 1
                for vtype, r_pos, r_base, q_base, bq in variants:
                    n_events[vtype] += 1
                    fout.write("\t".join([
                        args.barcode, read_name, str(r_pos),
                        r_base, q_base, str(bq), str(strand), str(len(seq)), vtype,
                    ]) + "\n")
                break  # 每个 read 只取一个通过的比对
            if read_primary:
                n_aligned += 1

    log.info("total reads=%d", total)
    log.info("primary aligned=%d", n_aligned)
    log.info("full-length (qcov=1,rcov=1,ident>=%.2f)=%d", args.min_identity, n_full)
    log.info("reads with variant=%d", n_var_reads)
    log.info("events: SNP=%d INS=%d DEL=%d", n_events["SNP"], n_events["INS"], n_events["DEL"])
    log.info("输出: %s", args.output)


if __name__ == "__main__":
    main()
