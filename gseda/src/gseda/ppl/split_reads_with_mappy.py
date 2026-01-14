import os
import sys
from multiprocessing import Pool, Value, Lock
from glob import glob
import argparse
import time
import mappy
import re

# =========================
# 全局计数器（多进程安全）
# =========================
TOTAL_READS = Value('i', 0)
SPLIT_READS = Value('i', 0)

TOTAL_BASES_BEFORE = Value('l', 0)
TOTAL_BASES_AFTER = Value('l', 0)

TOTAL_READS_AFTER = Value('i', 0)

COUNTER_LOCK = Lock()


# =========================
# 生信逻辑
# =========================
def reverse_complement(seq):
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'
    }
    return "".join(complement.get(b, b) for b in reversed(seq))


def calculate_identity_from_cigar(cigar_string):
    """
    根据CIGAR字符串计算identity

    Parameters:
    cigar_string (str): CIGAR字符串, 如 "370=1I89=1D26=1I208=1D1250=1D578=1D2=2X1=1I578=1I1248=1I210=1D23=1I90=1D373="

    Returns:
    float: identity值
    """
    # 使用正则表达式解析CIGAR字符串
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    # 初始化计数器
    match_count = 0      # 匹配的碱基数（=）
    mismatch_count = 0   # 错配的碱基数（X）
    total_aligned = 0    # 总比对长度（= + X + I + D）

    for length_str, operation in matches:
        length = int(length_str)

        if operation == '=':
            # 精确匹配
            match_count += length
            total_aligned += length
        elif operation == 'X':
            # 错配
            mismatch_count += length
            total_aligned += length
        elif operation == 'I':
            # 插入：参考序列没有，查询序列有
            total_aligned += length
        elif operation == 'D':
            # 删除：参考序列有，查询序列没有
            total_aligned += length

    # 计算identity
    if total_aligned == 0:
        return 0.0

    identity = match_count / total_aligned
    return identity


def try_split_seq(seq: str):
    seq_len = len(seq)
    if seq_len > 20000:
        return None
    aligner = mappy.Aligner(seq=seq, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    for hit in aligner.map(seq):
        if hit.r_st == 0 and hit.r_en == seq_len and hit.strand == 1 and calculate_identity_from_cigar(hit.cigar_str) >= 0.99999:
            continue
        # if hit.r_st == hit.q_st and hit
        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage = (hit.r_en - hit.r_st) / seq_len
        if identity > 0.85 and coverage > 0.85:
            return (hit.r_en + hit.r_st) // 2
    return None

def try_split_seq_2(seq: str):
    seq_len = len(seq)
    if seq_len > 20000:
        return None
    
    seq1 = seq[:(seq_len // 2)]
    seq2 = seq[(seq_len // 2):]
    
    aligner = mappy.Aligner(seq=seq1, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    for hit in aligner.map(seq2):
        
        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage1 = (hit.r_en - hit.r_st) / len(seq1)
        coverage2 = (hit.q_en - hit.q_st) / len(seq2)
        if identity > 0.85 and coverage1 > 0.85 and coverage2 > 0.85:
            
            # print("ctg:{}\tref_start:{}\tref_end:{}\tq_start:{}\tq_end:{}\tstrand:{}\tref_len:{}\tquery_len:{}. ideneity:{:.4f}".format(
            #     hit.ctg, hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.strand, len(seq1), len(seq2), calculate_identity_from_cigar(hit.cigar_str)))
            
            return (hit.r_en + hit.r_st) // 2
    return None


# =========================
# 单个 FASTQ 文件处理
# =========================
def process_fastq(fq_path, out_dir, threshold):
    basename = os.path.basename(fq_path)
    out_path = os.path.join(
        out_dir, basename.replace(".fastq", ".split.fastq"))

    local_total = 0
    local_split = 0

    local_bases_before = 0
    local_bases_after = 0
    local_reads_after = 0

    with open(fq_path) as fin, open(out_path, "w") as fout:
        while True:
            header = fin.readline().rstrip()
            if not header:
                break
            seq = fin.readline().rstrip()
            plus = fin.readline().rstrip()
            qual = fin.readline().rstrip()

            local_total += 1

            read_len = len(seq)
            local_bases_before += read_len

            if len(seq) < 200:
                fout.write(f"{header}\n{seq}\n+\n{qual}\n")
                continue

            split_position = try_split_seq_2(seq=seq)

            if split_position is not None:
                local_split += 1
                

                fout.write(
                    f"{header}_split_L\n"
                    f"{seq[:split_position]}\n+\n"
                    f"{qual[:split_position]}\n"
                )
                fout.write(
                    f"{header}_split_R\n"
                    f"{seq[split_position:]}\n+\n"
                    f"{qual[split_position:]}\n"
                )

                local_bases_after += read_len
                local_reads_after += 2
            else:
                fout.write(f"{header}\n{seq}\n+\n{qual}\n")
                local_bases_after += read_len
                local_reads_after += 1

    # === 安全更新全局计数 ===
    with COUNTER_LOCK:
        TOTAL_READS.value += local_total
        SPLIT_READS.value += local_split
        TOTAL_BASES_BEFORE.value += local_bases_before
        TOTAL_BASES_AFTER.value += local_bases_after
        TOTAL_READS_AFTER.value += local_reads_after


# =========================
# 多进程 worker
# =========================
def worker(args):
    fq_files, out_dir, threshold = args
    for fq in fq_files:
        assert isinstance(fq, str)
        if fq.endswith("unidentified.fastq"):
            continue
        process_fastq(fq, out_dir, threshold)


# =========================
# 主函数
# =========================
def main(fq_dir, out_dir, threshold=0.85, threads=24):

    start = time.time()

    os.makedirs(out_dir, exist_ok=True)

    fq_files = sorted(glob(os.path.join(fq_dir, "*.fastq")))
    if not fq_files:
        print("No fastq files found!")
        return

    # 按文件拆分任务
    chunk_size = max(1, len(fq_files) // threads)
    tasks = [
        (fq_files[i:i+chunk_size], out_dir, threshold)
        for i in range(0, len(fq_files), chunk_size)
    ]

    with Pool(threads) as pool:
        pool.map(worker, tasks)

    print("===== SUMMARY =====")
    print(f"Total reads : {TOTAL_READS.value}")
    print(f"Split reads : {SPLIT_READS.value}")
    print(f"Split ratio : {SPLIT_READS.value / max(1, TOTAL_READS.value):.4f}")
    before_avg = TOTAL_BASES_BEFORE.value / max(1, TOTAL_READS.value)
    after_avg = TOTAL_BASES_AFTER.value / max(1, TOTAL_READS_AFTER.value)

    print(f"Avg read length before split : {before_avg:.2f}")
    print(f"Avg read length after split  : {after_avg:.2f}")
    print(f"Elapsed Time: {time.time() - start}s")


def main_cli():
    parser = argparse.ArgumentParser(
        description="Split chimeric reads in FASTQ files using k-mer similarity"
    )

    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Directory containing input FASTQ files"
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Directory to write output FASTQ files"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.85,
        help="Similarity threshold for splitting (default: 0.85)"
    )
    parser.add_argument(
        "-p", "--threads",
        type=int,
        default=24,
        help="Number of worker processes (default: 24)"
    )

    args = parser.parse_args()
    main(args.input_dir, args.output_dir,
         threshold=args.threshold, threads=args.threads)


if __name__ == "__main__":
    # 用法:
    # python script.py fastq_dir output_dir 0.85 32
    main_cli()
