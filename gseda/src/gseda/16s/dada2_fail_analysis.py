import math
import sys
import argparse
from collections import defaultdict
import numpy as np
from tqdm import tqdm


# -----------------------------
# Q model (correct)
# -----------------------------
_Q_MAX = 93
_Q_TO_ER = np.array([10 ** (-q / 10) for q in range(_Q_MAX + 1)], dtype=np.float64)


def er_to_q(er: np.ndarray) -> np.ndarray:
    """vectorized error-rate -> Q"""
    er = np.maximum(er, 1e-300)
    q = -10.0 * np.log10(er)
    return np.clip(q, 0, _Q_MAX).astype(np.int32)


# -----------------------------
# FASTQ streaming parser
# -----------------------------
def fastq_iter(fp):
    while True:
        h = fp.readline()
        if not h:
            break
        if not h.startswith(b"@"):
            continue

        seq = fp.readline().strip()
        fp.readline()
        qual = fp.readline().strip()

        seq = seq.decode()
        q = np.frombuffer(qual, dtype=np.uint8) - 33

        yield seq, q


# -----------------------------
# STREAMING DEREP ENGINE
# -----------------------------
class StreamingDerep:
    """
    O(1) per sequence memory:
    only stores:
      - count
      - sum_error_rate vector
    """

    def __init__(self):
        self.counts = {}
        self.er_sums = {}

    def add(self, seq, q_array):
        er = _Q_TO_ER[q_array]

        if seq not in self.counts:
            self.counts[seq] = 0
            self.er_sums[seq] = np.zeros(len(q_array), dtype=np.float64)

        self.counts[seq] += 1
        self.er_sums[seq] += er

    def finalize(self):
        result = []
        for seq, cnt in tqdm(self.counts.items(), desc=f"compute mean Q"):
            er_mean = self.er_sums[seq] / cnt
            q_mean = er_to_q(er_mean)
            mean_er = -10 * math.log10(float(er_mean.mean()))

            result.append((seq, cnt, q_mean, mean_er))

        result.sort(key=lambda x: (-x[1], x[0]))
        return result


# -----------------------------
# FASTQ writer
# -----------------------------
def write_fastq(records, out):
    for i, (seq, cnt, q, reads_q) in tqdm(enumerate(records), desc="writing ... "):
        qual = "".join(chr(int(x) + 33) for x in q)
        out.write(f"@read_{i} count={cnt} len={len(qual)} Q={reads_q}\n")
        out.write(seq + "\n+\n")
        out.write(qual + "\n")


def write_tsv(records, out):
    out.write("seq\tcount\tmean_ER\tlen\n")
    for seq, cnt, q, er in records:
        out.write(f"{seq}\t{cnt}\t{er:.8f}\t{len(seq)}\n")


# -----------------------------
# MAIN PIPELINE
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastq")
    ap.add_argument("--out-fastq")
    ap.add_argument("--out-tsv")
    args = ap.parse_args()

    derep = StreamingDerep()

    with open(args.fastq, "rb") as f:
        for seq, q in tqdm(fastq_iter(f), desc=f"reading ..."):
            derep.add(seq, q)

    results = derep.finalize()

    if args.out_fastq:
        print(f"write to {args.out_fastq}")
        with open(args.out_fastq, "w") as f:
            write_fastq(results, f)
    else:
        write_fastq(results, sys.stdout)

    if args.out_tsv:
        with open(args.out_tsv, "w") as f:
            write_tsv(results, f)


if __name__ == "__main__":
    main()