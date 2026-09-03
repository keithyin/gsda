#!/usr/bin/env python3
"""Pairwise sequence identity and coverage from a multi-sequence FASTA file.

Uses mappy to index all sequences in the input FASTA as targets, then maps
each sequence against all others. Outputs identity and coverage matrices as
TSV files (one row/column per sequence).

Example:
    python fasta_pairwise_identity_coverage.py sequences.fasta
    python fasta_pairwise_identity_coverage.py sequences.fasta --no-self -o pw_out
"""

import mappy
import pysam
import argparse
import sys
import os
import pathlib
from collections import defaultdict

cur_path = pathlib.Path(os.path.abspath(__file__))
sys.path.append(str(cur_path.parent))

from mappy_ext import calculate_identity_from_cigar  # noqa: E402


def read_fasta(fpath):
    """Read a multi-sequence FASTA file, return {name: sequence} dict."""
    seqs = {}
    with pysam.FastxFile(fpath) as f:
        for record in f:
            seqs[record.name] = record.sequence
    return seqs


def build_matrices(seq_names, pairwise_results):
    """Build identity and coverage matrices from pairwise results.

    Parameters
    ----------
    seq_names : list[str]
        Ordered list of sequence names (row/column order).
    pairwise_results : dict[str, dict[str, dict]]
        Nested dict: results[query_name][target_name] = {identity, q_cov, t_cov, n_hits}

    Returns
    -------
    tuple[dict, dict]
        (identity_matrix, coverage_matrix) each is {name: {name: float}}
    """
    identity_matrix = {n: {m: 0.0 for m in seq_names} for n in seq_names}
    coverage_matrix = {n: {m: 0.0 for m in seq_names} for n in seq_names}

    for q_name, targets in pairwise_results.items():
        if q_name not in identity_matrix:
            continue
        for t_name, hit_info in targets.items():
            if t_name not in identity_matrix.get(q_name, {}):
                continue
            identity_matrix[q_name][t_name] = hit_info["identity"]
            coverage_matrix[q_name][t_name] = (hit_info["q_cov"] + hit_info["t_cov"]) / 2.0

    return identity_matrix, coverage_matrix


def write_matrix_tsv(seq_names, matrix, out_path):
    """Write a seq-by-seq matrix as TSV."""
    with open(out_path, "w", encoding="utf8") as f:
        # Header row
        f.write("Name\t" + "\t".join(f"{n:>30}" for n in seq_names) + "\n")
        # Data rows
        for row_name in seq_names:
            vals = "\t".join(f"{matrix[row_name][col]:.4f}" for col in seq_names)
            f.write(f"{row_name:>30}\t{vals}\n")


def pairwise_align(seq_a_name, seq_a_seq, seq_b_name, seq_b_seq, kmer, word_size):
    """Align seq_a against seq_b (indexed as reference).

    Returns dict with identity, q_cov, t_cov or None if no alignment.
    """
    aligner = mappy.Aligner(
        seq=seq_b_seq,
        extra_flags=67108864,
        k=kmer,
        w=word_size,
        best_n=1,
        n_threads=1,
    )
    best = None
    for hit in aligner.map(seq_a_seq):
        identity = calculate_identity_from_cigar(hit.cigar_str)
        entry = {
            "identity": identity,
            "q_cov": (hit.q_en - hit.q_st) / len(seq_a_seq),
            "t_cov": (hit.r_en - hit.r_st) / len(seq_b_seq),
            "r_st": hit.r_st,
            "r_en": hit.r_en,
            "q_st": hit.q_st,
            "q_en": hit.q_en,
        }
        if best is None or identity > best["identity"]:
            best = entry
    return best


def main(args):
    """Run pairwise identity/coverage alignment."""
    # 1. Read all sequences
    seqs = read_fasta(args.input)
    seq_names = list(seqs.keys())

    if len(seqs) < 2:
        print("Error: input FASTA must contain at least 2 sequences.", file=sys.stderr)
        sys.exit(1)

    n_seqs = len(seqs)
    # Unique unordered pairs (i < j)
    n_pairs = n_seqs * (n_seqs - 1) // 2
    if not args.no_self:
        n_pairs += n_seqs
    print(f"Loaded {n_seqs} sequences, computing pairwise alignments ...", file=sys.stderr)

    # 2. Compute pairwise alignments for each unique pair (i < j)
    # results[a][b] and results[b][a] are populated symmetrically
    results = defaultdict(lambda: defaultdict(list))

    pairs_processed = 0
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            name_a, name_b = seq_names[i], seq_names[j]
            seq_a, seq_b = seqs[name_a], seqs[name_b]

            # Align A -> B (B indexed as reference)
            hit_ab = pairwise_align(name_a, seq_a, name_b, seq_b, args.kmer, args.word_size)
            if hit_ab is not None:
                results[name_a][name_b].append(hit_ab)
                results[name_b][name_a].append({
                    "identity": hit_ab["identity"],
                    "q_cov": hit_ab["t_cov"],  # swap q/t for reverse direction
                    "t_cov": hit_ab["q_cov"],
                    "r_st": hit_ab["q_st"],
                    "r_en": hit_ab["q_en"],
                    "q_st": hit_ab["r_st"],
                    "q_en": hit_ab["r_en"],
                })
            pairs_processed += 1

            # Self-alignments if requested
            if not args.no_self:
                for name in (name_a, name_b):
                    seq = seqs[name]
                    results[name][name].append({
                        "identity": 1.0,
                        "q_cov": 1.0,
                        "t_cov": 1.0,
                        "n_hits": 0,
                    })

            # Progress indicator for large inputs
            if (pairs_processed + 1) % 50 == 0:
                print(f"  ... {pairs_processed}/{n_pairs} pairs processed", file=sys.stderr)

    # Take best-hit per directed pair (highest identity)
    best_results = {}
    for q_name in results:
        best_results[q_name] = {}
        for t_name in results[q_name]:
            hits = results[q_name][t_name]
            best_results[q_name][t_name] = max(hits, key=lambda h: h["identity"])

    # 3. Build matrices
    identity_matrix = {n: {m: 0.0 for m in seq_names} for n in seq_names}
    coverage_matrix = {n: {m: 0.0 for m in seq_names} for n in seq_names}
    for q_name, targets in best_results.items():
        if q_name not in identity_matrix:
            continue
        for t_name, hit_info in targets.items():
            if t_name not in identity_matrix.get(q_name, {}):
                continue
            identity_matrix[q_name][t_name] = hit_info["identity"]
            coverage_matrix[q_name][t_name] = (hit_info["q_cov"] + hit_info["t_cov"]) / 2.0

    # 4. Write output files
    prefix = args.out_prefix
    id_path = f"{prefix}_identity.tsv"
    cov_path = f"{prefix}_coverage.tsv"
    write_matrix_tsv(seq_names, identity_matrix, id_path)
    write_matrix_tsv(seq_names, coverage_matrix, cov_path)

    # 5. Summary to stderr
    print(f"\n=== Pairwise alignment summary ===", file=sys.stderr)
    print(f"Sequences: {n_seqs}", file=sys.stderr)
    aligned_pairs = sum(1 for q in seq_names for t in seq_names if identity_matrix[q][t] > 0.0 and (not args.no_self or q != t))
    print(f"Aligned pairs: {aligned_pairs} / {n_pairs}", file=sys.stderr)

    # Compute stats over aligned non-self pairs
    identities = []
    covs = []
    for q in seq_names:
        for t in seq_names:
            if args.no_self and q == t:
                continue
            val_id = identity_matrix[q][t]
            val_cov = coverage_matrix[q][t]
            if val_id > 0.0:
                identities.append(val_id)
                covs.append(val_cov)

    if identities:
        print(f"Identity   : min={min(identities):.4f}  max={max(identities):.4f}  mean={sum(identities)/len(identities):.4f}", file=sys.stderr)
        print(f"Coverage   : min={min(covs):.4f}  max={max(covs):.4f}  mean={sum(covs)/len(covs):.4f}", file=sys.stderr)
    print(f"\nOutput files:", file=sys.stderr)
    print(f"  Identity   : {id_path}", file=sys.stderr)
    print(f"  Coverage   : {cov_path}", file=sys.stderr)

    return id_path, cov_path


def main_cli():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="fasta_pairwise_identity_coverage",
        description="Compute pairwise sequence identity and coverage from a multi-sequence FASTA file.",
    )
    parser.add_argument("input", help="Input FASTA file (multiple sequences)")
    parser.add_argument("--no-self", action="store_true",
                        help="Exclude self-alignments from output")
    parser.add_argument("-o", "--out-prefix", default="pw",
                        help="Output filename prefix (default: pw)")
    parser.add_argument("--kmer", type=int, default=7,
                        help="k-mer size for indexing (default: 7)")
    parser.add_argument("-w", "--word-size", type=int, default=5,
                        help="Word size for seeding (default: 5)")
    parser.add_argument("--best-n", type=int, default=10,
                        help="Max hits per query to report (default: 10)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (default: 1)")

    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    main_cli()
