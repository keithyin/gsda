import os
import argparse
import time
import mappy
import re
import pysam
from multiprocessing import Process, Queue, Value, Lock, cpu_count
import signal
import sys

# =========================
# 全局计数器（多进程安全）
# =========================
TOTAL_READS = Value('i', 0)
SPLIT_READS = Value('i', 0)
COUNTER_LOCK = Lock()

# Sentinel for queue termination
_SENTINEL = None


def calculate_identity_from_cigar(cigar_string):
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    match_count = 0
    mismatch_count = 0
    total_aligned = 0

    for length_str, operation in matches:
        length = int(length_str)
        if operation == '=':
            match_count += length
            total_aligned += length
        elif operation == 'X':
            mismatch_count += length
            total_aligned += length
        elif operation in ('I', 'D'):
            total_aligned += length

    if total_aligned == 0:
        return 0.0
    return match_count / total_aligned


def try_split_seq_2(seq: str):
    seq_len = len(seq)
    if seq_len < 200 or seq_len > 20000:
        return False

    mid = seq_len // 2
    seq1 = seq[:mid]
    seq2 = seq[mid:]

    aligner = mappy.Aligner(seq=seq1, extra_flags=67108864,
                            k=9, w=7, best_n=10, n_threads=1)
    for hit in aligner.map(seq2):
        identity = calculate_identity_from_cigar(hit.cigar_str)
        coverage1 = (hit.r_en - hit.r_st) / len(seq1)
        coverage2 = (hit.q_en - hit.q_st) / len(seq2)
        if identity > 0.80 and coverage1 > 0.85 and coverage2 > 0.85:
            return True
    return False


# =========================
# 消费者 Worker
# =========================
def consumer(queue):
    local_total = 0
    local_split = 0

    while True:
        seq = queue.get()
        if seq is _SENTINEL:
            # Put sentinel back for other consumers
            queue.put(_SENTINEL)
            break

        if seq and len(seq) > 0:
            local_total += 1
            if try_split_seq_2(seq):
                local_split += 1

    # Update global counters
    with COUNTER_LOCK:
        TOTAL_READS.value += local_total
        SPLIT_READS.value += local_split


# =========================
# 生产者：从 BAM 读取并放入队列
# =========================
def producer(bam_path, queue):
    count = 0
    with pysam.AlignmentFile(bam_path, "rb", threads=10, check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            seq = read.query_sequence
            if seq:
                queue.put(seq)
                count += 1
                if count % 10000 == 0:
                    print(f"\rProduced {count} reads...",
                          end='', flush=True)
    print(f"\nProducer finished. Total reads enqueued: {count}")


# =========================
# 主函数
# =========================
def main(bam_path, threads):
    start = time.time()

    # 创建队列（适当大小避免爆内存）
    queue = Queue(maxsize=threads * 10)

    # 启动消费者进程
    consumers = []
    for _ in range(threads):
        p = Process(target=consumer, args=(queue,))
        p.start()
        consumers.append(p)

    # 启动生产者（在主进程中运行，简化控制）
    print("Starting producer...")
    producer(bam_path, queue)

    # 发送终止信号
    queue.put(_SENTINEL)

    # 等待所有消费者结束
    for p in consumers:
        p.join()

    # 输出结果
    total = TOTAL_READS.value
    split = SPLIT_READS.value
    ratio = split / total if total > 0 else 0.0
    report_str = f"""
================= Macebell Report =================
- Reads Processed : {total}
- Macebell Reads  : {split}
- Macebell Ratio  : {ratio*100: .2f}%
- Elapsed Time    : {time.time() - start:.2f}s
=======================================================
"""
    print(report_str)
    return report_str


def main_cli():
    parser = argparse.ArgumentParser(
        description="Estimate chimeric read split ratio from BAM using producer-consumer model."
    )
    parser.add_argument(
        "-i", "--input-bam",
        required=True,
        help="Input BAM file"
    )
    parser.add_argument(
        "-p", "--threads",
        type=int,
        default=cpu_count(),
        help="Number of consumer processes (default: all CPUs)"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input_bam):
        raise FileNotFoundError(f"BAM file not found: {args.input_bam}")

    # Handle Ctrl+C gracefully
    def signal_handler(sig, frame):
        print("\nInterrupted by user. Exiting...")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)

    main(args.input_bam, args.threads)


if __name__ == "__main__":
    main_cli()
