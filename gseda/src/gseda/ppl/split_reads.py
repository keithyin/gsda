import os
import sys
from multiprocessing import Pool, Value, Lock
from glob import glob
import argparse
import time

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


def check_kmer_similarity(s1, s2, k=13, step=3):
    if len(s1) < k or len(s2) < k:
        return 0
    set1 = {s1[i:i+k] for i in range(0, len(s1)-k+1, step)}
    set2 = {s2[i:i+k] for i in range(0, len(s2)-k+1, step)}
    if not set1 or not set2:
        return 0
    return len(set1 & set2) / min(len(set1), len(set2))


# =========================
# 单个 FASTQ 文件处理
# =========================
def process_fastq(fq_path, out_dir, threshold):
    basename = os.path.basename(fq_path)
    out_path = os.path.join(out_dir, basename.replace(".fastq", ".split.fastq"))

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

            mid = len(seq) // 2
            l_seq, r_seq = seq[:mid], seq[mid:]
            l_qual, r_qual = qual[:mid], qual[mid:]

            sim_rc = check_kmer_similarity(l_seq, reverse_complement(r_seq))
            sim_fwd = check_kmer_similarity(l_seq, r_seq)

            if max(sim_rc, sim_fwd) >= threshold:
                local_split += 1
                tag = "RC" if sim_rc >= sim_fwd else "FWD"
                
                left_len = mid + 20
                right_len = len(seq) - (mid - 20)

                fout.write(
                    f"{header}_split_L_{tag}\n"
                    f"{seq[:mid+20]}\n+\n"
                    f"{qual[:mid+20]}\n"
                )
                fout.write(
                    f"{header}_split_R_{tag}\n"
                    f"{seq[mid-20:]}\n+\n"
                    f"{qual[mid-20:]}\n"
                )
                
                local_bases_after += left_len + right_len
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
    main(args.input_dir, args.output_dir, threshold=args.threshold, threads=args.threads)
    

if __name__ == "__main__":
    # 用法:
    # python script.py fastq_dir output_dir 0.85 32
    main_cli()
