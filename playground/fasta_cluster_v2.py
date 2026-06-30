#!/usr/bin/env python3
"""
fasta_cluster_v2.py — all-vs-all clustering via mappy AVA.

Architecture (key optimization: zero preparation overhead):
  Main process:
    - Loads FASTA once
    - Splits into N batches of query sequences
    - Launches nproc workers

  Each worker (process):
    - Loads ALL sequences from FASTA on disk  (one pass per worker, not per pair)
    - Builds ONE mappy index                   (one index per worker, not per pair)
    - Maps only its assigned batch of queries  (N pairs per worker)

Total work = nproc × (1 FASTA read + 1 index build) + N×(N-1)/2 × mappy map
"""
import argparse
import mappy
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm


# ===== Union-Find (path compression) — copied verbatim from fasta_cluster.py =====
class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    @staticmethod
    def from_edges(edges):
        """Build a UnionFind from a list of (a, b) edges."""
        uf = UnionFind()
        for a, b in edges:
            uf.union(a, b)
        return uf


# ===== FASTA reader =====
def read_fasta(fpath: str):
    """Yield (name, sequence) tuples from a FASTA file."""
    name, seq = None, []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq)


# ===== Identity from CIGAR =====
_cigar_re = re.compile(r'(\d+)([=IDX])')


def _identity_from_cigar(cigar_str: str):
    """Compute alignment identity from a mappy compact CIGAR string."""
    m = 0  # matches
    aln = 0  # total aligned length
    for length_str, op in _cigar_re.findall(cigar_str):
        length = int(length_str)
        aln += length
        if op == '=':
            m += length
    return m / aln if aln else 0.0


# ===== Coverage helper =====
def _hit_cov(q_st, q_en, t_st, t_en, qlen, tlen):
    """Compute the stricter of query/target coverage."""
    qcov = (q_en - q_st) / qlen if qlen else 0
    tcov = (t_en - t_st) / tlen if tlen else 0
    return min(qcov, tcov)


# ===== Worker =====
def _worker(queries, fasta, k, w, min_cov, min_id, n_threads):
    """Worker process: load FASTA, build index, map queries.

    Each worker reads FASTA once and builds one mappy index,
    then maps only its batch of queries against it.
    """
    # Load all sequences from FASTA and build index (once per worker)
    all_seqs = list(read_fasta(fasta))
    seq_len = {name: len(seq) for name, seq in all_seqs}

    # Build mappy index
    ref_seqs = [seq for _, seq in all_seqs]
    ref_idx = mappy.index(ref_seqs, k=k, w=w, mask=-1)
    aligner = mappy.Aligner(ref_idx, n_threads=n_threads)

    edges = []
    for qname, qseq in queries:
        for hit in aligner.map(qseq):
            tname = hit.ctg
            if tname == qname:
                continue
            if not hit.is_primary:
                continue
            tlen = seq_len.get(tname, 0)
            if tlen == 0:
                continue
            cov = _hit_cov(hit.q_st, hit.q_en, hit.r_st, hit.r_en,
                           hit.q_en - hit.q_st, tlen)
            if cov < min_cov:
                continue
            identity = _identity_from_cigar(hit.cigar_str)
            if identity < min_id:
                continue
            edges.append((qname, tname))
    return edges


# ===== Main clustering pipeline =====
def cluster_with_mappy(fasta: str, min_cov: float = 0.99, min_id: float = 0.999,
                       k: int = 15, w: int = 10, n_threads: int = 1, nproc: int = 1):
    """All-vs-all cluster sequences using mappy.

    Each sequence compared exactly once against every other → N*(N-1)/2 pairs.

    Parameters
    ----------
    fasta : str
        Path to input FASTA.
    min_cov : float
        Minimum alignment coverage (stricter of query/target).
    min_id : float
        Minimum alignment identity.
    k, w : int
        mappy k-mer size and seed interval.
    n_threads : int
        mappy internal threads (per process).
    nproc : int
        Number of parallel processes.

    Returns
    ----- --
    UnionFind
        Clustering result.
    """
    # 1. load all sequences
    seqs = list(read_fasta(fasta))
    n = len(seqs)
    print(f"Loaded {n} sequences from {fasta}")

    if n <= 1:
        return UnionFind()

    total_pairs = n * (n - 1) // 2
    print(f"Will compare {total_pairs:,} pairs")

    # 2. split queries into batches (O(n) only, no preparation of targets)
    #    This is just indexing — no tuple/string allocation
    batches = []
    batch_size = max(1, (n - 1) // nproc)
    for i in range(0, n - 1, batch_size):
        end = min(i + batch_size, n - 1)
        batches.append([(name, seq) for name, seq in seqs[i:end]])

    # 3. parallel mapping: each worker loads FASTA + builds index independently
    #    Worker receives only its queries (O(batch_size) data per worker)
    all_edges = []
    with ProcessPoolExecutor(max_workers=nproc) as executor:
        with tqdm(total=len(batches), desc="AVA mapping", unit="batch") as pbar:
            futures = {
                executor.submit(_worker, batch, fasta, k, w, min_cov, min_id, n_threads): idx
                for idx, batch in enumerate(batches)
            }
            for fut in as_completed(futures):
                all_edges.extend(fut.result())
                pbar.update(1)

    return UnionFind.from_edges(all_edges)


# ===== Build clusters from UnionFind =====
def build_clusters(uf, fasta):
    """Group all sequence names by their cluster root."""
    clusters = defaultdict(list)
    for name, _ in read_fasta(fasta):
        root = uf.find(name)
        clusters[root].append(name)
    return dict(clusters)


# ===== Write cluster TSV =====
def write_clusters(clusters, out):
    with open(out, "w") as f:
        f.write("sequence\tcluster\trepresentative\n")
        for cid, members in enumerate(clusters.values()):
            rep = members[0]
            for m in members:
                f.write(f"{m}\tcluster_{cid}\t{rep}\n")


# ===== CLI entry point =====
def main():
    parser = argparse.ArgumentParser(
        description="All-vs-all sequence clustering via mappy AVA")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("-o", "--output", default="clusters_v2.tsv",
                        help="Output TSV (default: clusters_v2.tsv)")
    parser.add_argument("--min-cov", type=float, default=0.99,
                        help="Minimum coverage (default: 0.99)")
    parser.add_argument("--min-id", type=float, default=1.0,
                        help="Minimum identity (default: 0.999)")
    parser.add_argument("-k", type=int, default=15,
                        help="mappy k-mer size (default: 15)")
    parser.add_argument("-w", type=int, default=10,
                        help="mappy seed interval (default: 10)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="mappy threads per process (default: 1)")
    parser.add_argument("-p", "--nproc", type=int, default=1,
                        help="Number of parallel processes (default: 1)")
    args = parser.parse_args()

    uf = cluster_with_mappy(args.fasta, min_cov=args.min_cov, min_id=args.min_id,
                            k=args.k, w=args.w, n_threads=args.threads,
                            nproc=args.nproc)
    clusters = build_clusters(uf, args.fasta)
    write_clusters(clusters, args.output)

    sizes = [len(m) for m in clusters.values()]
    print(f"Done: {len(clusters)} clusters ({min(sizes)}–{max(sizes)} members)")


if __name__ == "__main__":
    main()
