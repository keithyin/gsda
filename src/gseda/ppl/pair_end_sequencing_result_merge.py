import argparse
import logging
import os

import mappy

from gseda.align_ana.mappy_ext import calculate_identity_from_cigar, parse_cigar, revcomp

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)

SEQ_EXT = ".seq"
MIN_IDENTITY = 0.98


def sample_name(fname: str) -> str:
    """第一个下划线前为样本编号。"""
    return fname.split("_")[0]


def read_sanger_seqs(input_dir: str) -> dict:
    """把 sanger_result 下的 .seq 文件按样本聚合为 {'sample': {'F': seq, 'R': seq}}。

    .seq 文件内容为无 FASTA 头的纯序列。引物前缀为 M13F / M13R。
    """
    by_sample = {}
    for f in sorted(os.listdir(input_dir)):
        if not f.endswith(SEQ_EXT):
            continue
        primer = f.split("_")[1] if "_" in f else ""
        if primer.startswith("M13F"):
            strand = "F"
        elif primer.startswith("M13R"):
            strand = "R"
        else:
            continue
        with open(os.path.join(input_dir, f), encoding="utf8") as fh:
            seq = fh.read().strip()
        by_sample.setdefault(sample_name(f), {})[strand] = seq
    return by_sample


def build_aligner(f_seq: str) -> mappy.Aligner:
    return mappy.Aligner(
        seq=f_seq,
        extra_flags=67108864,
        k=7,
        w=5,
        best_n=10,
        n_threads=1,
        min_cnt=2,
        min_chain_score=4,
    )


def align_to_reference(query: str, ref: str):
    """将 query 比对到 ref，返回 identity / query_coverage / target_coverage。

    query_coverage = 比对到的 query 区间长度 / query 全长
    target_coverage = 比对到的 ref 区间长度 / ref 全长
    """
    aligner = build_aligner(ref)
    hits = list(aligner.map(query))
    primary_hits = [h for h in hits if h.is_primary]
    if not primary_hits:
        return None
    hit = primary_hits[0]
    identity = calculate_identity_from_cigar(hit.cigar_str)
    query_cov = (hit.q_en - hit.q_st) / len(query) if query else 0.0
    target_cov = (hit.r_en - hit.r_st) / len(ref) if ref else 0.0
    return {
        "identity": identity,
        "query_coverage": query_cov,
        "target_coverage": target_cov,
    }


def merge_pair(f_seq: str, r_seq: str):
    """将正向 F 与反向 R 比对合并为一条全长序列。

    R 读段从片段 3' 端反向测序，因此将 R 反向互补为 Rrc 后与 F 同向比对。
    返回 dict：merged / identity / overlap_len / lead / tail / r_st / r_en / q_st / q_en。
    """
    aligner = build_aligner(f_seq)
    hits = list(aligner.map(r_seq))
    primary_hits = [h for h in hits if h.is_primary]
    if not primary_hits:
        return None

    hit = primary_hits[0]

    # 将比对坐标统一到 (F, Rrc) 正链坐标系
    r_st, r_en = hit.r_st, hit.r_en
    q_st, q_en = len(r_seq) - hit.q_en, len(r_seq) - hit.q_st
    identity = calculate_identity_from_cigar(hit.cigar_str)

    rrc = revcomp(r_seq)
    cols = []
    f_pos, r_pos = r_st, q_st
    for length, op in parse_cigar(hit.cigar_str):
        if op in ("=", "X"):
            for _ in range(length):
                cols.append((f_seq[f_pos], rrc[r_pos]))
                f_pos += 1
                r_pos += 1
        elif op == "I":  # 插入：F 无、Rrc 有
            for _ in range(length):
                cols.append(("-", rrc[r_pos]))
                r_pos += 1
        elif op == "D":  # 缺失：F 有、Rrc 无
            for _ in range(length):
                cols.append((f_seq[f_pos], "-"))
                f_pos += 1
        elif op == "S":  # soft clip，仅推进 Rrc 坐标
            r_pos += length

    def is_mismatch(col):
        return col[0] != "-" and col[1] != "-" and col[0] != col[1]

    lead = 0
    while lead < len(cols) and is_mismatch(cols[lead]):
        lead += 1
    tail = 0
    while tail < len(cols) - lead and is_mismatch(cols[len(cols) - 1 - tail]):
        tail += 1

    kept = cols[lead:len(cols) - tail]
    overlap = "".join(c[0] if c[0] != "-" else c[1] for c in kept)

    # 5' 端取测得更远的 read，3' 端同样取测得更远的 read
    head = rrc[:q_st] if q_st >= r_st else f_seq[:r_st]
    tail_seq = f_seq[r_en:] if len(f_seq) - r_en >= len(r_seq) - q_en else rrc[q_en:]

    merged = head + overlap + tail_seq
    return {
        "merged": merged,
        "overlap": overlap,
        "identity": identity,
        "overlap_len": len(kept),
        "lead": lead,
        "tail": tail,
        "r_st": r_st,
        "r_en": r_en,
        "q_st": q_st,
        "q_en": q_en,
    }


def main_cli():
    parser = argparse.ArgumentParser(
        prog="pair_end_sequencing_result_merge",
        description="将一代测序（Sanger）双端结果合并为一条全长序列",
    )
    parser.add_argument("input_dir", type=str, help="包含 .seq 文件的目录")
    parser.add_argument("--outdir", type=str, default=None,
                        help="输出目录，默认 <input_dir>/merged_output")
    parser.add_argument("--union", action="store_true",
                        help="输出拼接后的全长序列（默认只输出两条 read 比对的重叠区/交集）")
    args = parser.parse_args()

    intersect_only = not args.union

    input_dir = args.input_dir
    outdir = args.outdir or os.path.join(input_dir, "merged_output")
    os.makedirs(outdir, exist_ok=True)

    mode = "intersect" if intersect_only else "union"
    sum_fa = os.path.join(outdir, "intersect.fa" if intersect_only else "union.fa")
    per_sample_suffix = ".intersect.fa" if intersect_only else ".union.fa"

    by_sample = read_sanger_seqs(input_dir)
    low_identity = []
    n_ok = 0
    merged_lines = []

    for sample in sorted(by_sample):
        seqs = by_sample[sample]
        if "F" not in seqs or "R" not in seqs:
            logging.warning("sample %s missing paired read: %s", sample, sorted(seqs.keys()))
            continue
        result = merge_pair(seqs["F"], seqs["R"])
        if result is None:
            logging.warning("sample %s: no mappy hit between forward/reverse read", sample)
            continue
        n_ok += 1
        out_seq = result["overlap"] if intersect_only else result["merged"]
        identity = result["identity"]
        logging.info(
            "sample %s: forward len=%d, reverse len=%d, overlap len=%d, output len=%d, identity=%.4f",
            sample, len(seqs["F"]), len(seqs["R"]), result["overlap_len"], len(out_seq), identity,
        )

        # 将原始两条 read 比对回实际输出的序列，记录 identity / queryCoverage / targetCoverage
        merged_seq = out_seq
        for label, raw_seq in (("forward", seqs["F"]), ("reverse", revcomp(seqs["R"]))):
            stats = align_to_reference(raw_seq, merged_seq)
            if stats is None:
                logging.warning("sample %s: %s read realigned to output but no hit", sample, label)
                continue
            logging.info(
                "sample %s: %s realigned to output identity=%.4f queryCoverage=%.4f targetCoverage=%.4f",
                sample, label, stats["identity"],
                stats["query_coverage"], stats["target_coverage"],
            )

        if identity < MIN_IDENTITY:
            low_identity.append((sample, identity))
            logging.warning("sample %s identity %.2f%% < 98%%", sample, identity * 100)

        header = f">{sample}\n"
        body = out_seq + "\n"
        with open(os.path.join(outdir, f"{sample}{per_sample_suffix}"), "w", encoding="utf8") as fh:
            fh.write(header)
            fh.write(body)
        merged_lines.append(header)
        merged_lines.append(body)

    with open(sum_fa, "w", encoding="utf8") as fh:
        fh.writelines(merged_lines)

    print(f"samples: {len(by_sample)}, merged ok: {n_ok}, low identity (<98%): {len(low_identity)}, mode: {mode}")
    if low_identity:
        print("low identity samples: " + ", ".join(
            f"{s}({i * 100:.2f}%)" for s, i in low_identity))


if __name__ == "__main__":
    main_cli()
