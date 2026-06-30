import argparse
import os
import pathlib
import sys
import re

import mappy
from tqdm import tqdm  # noqa: E402

cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))
sys.path.append(str(prev_prev_dir / "src"))


# ---------------------------------------------------------------------------
# FASTQ reader
# ---------------------------------------------------------------------------

def read_fastq(fpath: str):
    """读取 FASTQ 文件，返回 {name: (seq, qual)}"""
    records = {}
    with open(fpath, "r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            name = header[1:].split()[0]
            records[name] = (seq, qual)
    return records


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def calculate_identity_from_cigar(cigar_string: str) -> float:
    pattern = r"(\d+)([=IDX])"
    matches = re.findall(pattern, cigar_string)
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Align query FASTQ to reference FASTA, filter by qcov=1, rcov=1, identity>min_identity, "
                    "and export matching reads to a FASTQ file."
    )
    parser.add_argument("--ref", required=True, help="Reference FASTA file (single sequence)")
    parser.add_argument("--query", required=True, help="Query FASTQ file")
    parser.add_argument("--min-identity", type=float, default=0.99, help="Minimum identity (default: 0.99)")
    parser.add_argument("--output-matched", help="Output FASTQ file for reads with qcov==1, rcov==1, identity>min_identity")
    args = parser.parse_args()

    # Read reference
    ref_seqs = {}
    with open(args.ref, "r") as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].split()[0]
                ref_seqs[header] = ""
            elif header:
                ref_seqs[header] += line
    if len(ref_seqs) != 1:
        print(f"Error: reference should contain exactly 1 sequence, got {len(ref_seqs)}")
        sys.exit(1)
    ref_name = list(ref_seqs.keys())[0]
    ref_seq = ref_seqs[ref_name]
    print(f"Reference: {ref_name}, length: {len(ref_seq)}, seq={ref_seq[:100]}...")

    # Read query
    query_records = read_fastq(args.query)
    print(f"Query records: {len(query_records)}")

    # Align & filter
    aligner = mappy.Aligner(seq=ref_seq, extra_flags=67108864,
                            k=11, w=1, best_n=10, n_threads=1)

    matched_fh = open(args.output_matched, "w") if args.output_matched else None
    matched_count = 0
    not_matched_count = 0

    for qname, (qseq, qual_str) in tqdm(query_records.items(), desc="processing"):
        for hit in aligner.map(qseq):
            if not hit.is_primary:
                continue

            qcov = (hit.q_en - hit.q_st) / len(qseq)
            rcov = (hit.r_en - hit.r_st) / len(ref_seq)
            identity = calculate_identity_from_cigar(hit.cigar_str)

            if qcov == 1.0 and rcov == 1.0 and identity > args.min_identity:
                if matched_fh:
                    matched_fh.write(f"@{qname}\n{qseq}\n+\n{qual_str}\n")
                matched_count += 1
            else:
                not_matched_count += 1

    if matched_fh:
        matched_fh.close()

    print(f"\nTotal: {len(query_records)}")
    print(f"Matched (qcov=1, rcov=1, ident>{args.min_identity}): {matched_count}")
    print(f"Not matched: {not_matched_count}")
    if args.output_matched:
        print(f"Matched reads written to: {args.output_matched}")


if __name__ == "__main__":
    main()
