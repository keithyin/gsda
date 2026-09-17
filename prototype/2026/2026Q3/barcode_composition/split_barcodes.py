#!/usr/bin/env python
"""Split an smc_all_reads file (FASTQ or BAM) into per-sample barcodeNN.fastq.

Inputs
  --pool   64-barcode pool TSV (id group sequence len gc_pct n_samples n_reads)
  --pairs  768-pair TSV       (sample_id outer_id outer_seq inner_id inner_seq n_reads)
  positional input            .fq / .fastq (auto-detected by extension) or .bam

Read grammar (see extract_barcodes.py)
    (AAGGTTAA + 24nt + CAGCACCA)x2  ...  x2(TGGTGCTG + 24nt + TTAACCTT)
Each end's barcode is read from its two tandem copies (consensus), then
error-corrected to the pool at Hamming <= 3. The pair is treated as an
UNORDERED set (consensus strand is random):

    {A, B}, A != B  ->  the matching row of --pairs  ->  barcodeNN.fastq
    {A, A}          ->  self_pairs.fastq            (index-hopping leakage)
    unassignable    ->  unmapped.fastq

Quality strings are carried through unchanged; BAM qual is decoded from its
binary encoding to ASCII and re-emitted as a FASTQ record.
"""
import argparse
import os
import sys
from collections import Counter

import pysam

ANCHOR = "AAGGTTAA"
SEP = "CAGCACCA"
BC_LEN = 24
TAG_LEN = 8 + BC_LEN + 8
SCAN = 200
MM_TOL = 2     # for locating a tag (anchor+separator fingerprint)
MMD = 3        # for assigning a 24nt barcode to the pool


def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def mm(a, b):
    if len(a) != len(b):
        return 99
    return sum(1 for x, y in zip(a, b) if x != y)


def extract_tag(seq):
    """(barcode, score) of the first tag in seq[:SCAN], else (None, None)."""
    best_i, best_s = -1, 99
    lim = min(SCAN, len(seq) - TAG_LEN) + 1
    for i in range(lim):
        s = mm(seq[i:i + 8], ANCHOR)
        if s > MM_TOL:
            continue
        s += mm(seq[i + 8 + BC_LEN:i + TAG_LEN], SEP)
        if s < best_s:
            best_i, best_s = i, s
            if s == 0:
                break
    if best_i < 0 or best_s > MM_TOL:
        return None, None
    return seq[best_i + 8:best_i + 8 + BC_LEN], best_s


def consensus2(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return "".join(x if x == y else x for x, y in zip(a, b))


def end_tag(seq):
    """One end of a molecule: two tandem copies -> per-read consensus."""
    b1, _ = extract_tag(seq)
    if b1 is None:
        return None
    b2, _ = extract_tag(seq[40:44 + TAG_LEN + 2])
    return consensus2(b1, b2)


def load_pool(path):
    pool = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            pool.append((p[0], p[2]))
    return pool


def load_pairs(path):
    """frozenset(outer_seq, inner_seq) -> sample id, plus the full row count.

    Pairs file is header-less TSV:  sample_id  outer_seq  inner_seq
    """
    m = {}
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            n += 1
            m[frozenset((p[1], p[2]))] = p[0]
    return m, n


def build_lookup(pool):
    """raw 24nt string -> pool member (or None) within Hamming <= MMD, cached."""
    cache = {}

    def lookup(b):
        if b in cache:
            return cache[b]
        hit, best = None, MMD + 1
        for _, seq in pool:
            d = mm(b, seq)
            if d < best:
                best, hit = d, seq
                if d == 0:
                    break
        cache[b] = hit
        return hit

    return lookup


def iter_reads(fa):
    """Yield (name, seq, qual, is_bam, rq). rq is None for FASTQ input."""
    if fa.endswith(".bam") or fa.endswith(".bam.bai"):
        with pysam.AlignmentFile(fa) as sam:
            for r in sam:
                q = r.qual
                if isinstance(q, bytes):
                    q = q.decode("latin-1")
                yield r.qname, r.seq, (q or ""), True, r.get_tag("rq")
    else:
        with open(fa) as fh:
            while True:
                name = fh.readline().lstrip("@").rstrip()
                if not name:
                    break
                seq = fh.readline().rstrip()
                fh.readline()
                qual = fh.readline().rstrip()
                yield name, seq, qual, False, None


class Worker:
    def __init__(self, pool, pairmap, args):
        self.lookup = build_lookup(pool)
        self.pairmap = pairmap
        self.args = args
        self.counts = Counter()
        self.outs = {}
        self.min_rq = args.min_rq
        os.makedirs(args.outdir, exist_ok=True)

    def out(self, key):
        fh = self.outs.get(key)
        if fh is None:
            if key.startswith("S"):
                path = os.path.join(self.args.outdir, f"barcode{key[1:]}.fastq")
            else:
                path = os.path.join(self.args.outdir, key + ".fastq")
            fh = self.outs[key] = open(path, "w")
        return fh

    def handle(self, rec):
        name, seq, qual, is_bam, rq = rec
        if rq is None:
            self.counts["no_rq"] += 1
        elif rq < self.min_rq:
            self.counts["rq_filtered"] += 1
            return
        if not seq:
            self.counts["NOSKIP"] += 1
            return
        L = end_tag(seq)
        Rr = end_tag(rc(seq))
        if L is None or Rr is None:
            self.write("unmapped", name, seq, qual, is_bam)
            return
        A = self.lookup(L)
        B = self.lookup(Rr)
        if A is None or B is None:
            self.write("unmapped", name, seq, qual, is_bam)
            return
        if A == B:
            self.write("self_pairs", name, seq, qual, is_bam)
            return
        sid = self.pairmap.get(frozenset((A, B)))
        if sid is None:
            self.write("unmapped", name, seq, qual, is_bam)
            return
        self.write(sid, name, seq, qual, is_bam)

    def write(self, key, name, seq, qual, is_bam):
        fh = self.out(key)
        fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
        self.counts[key] += 1

    def flush(self):
        for fh in self.outs.values():
            fh.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="smc_all_reads .fq or .bam")
    ap.add_argument("--pool", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max-reads", type=int, default=0, help="0 = all (debug)")
    ap.add_argument("--min-rq", type=float, default=0.0,
                    help="drop BAM reads with rq < this (BAM only; 0 = no filter)")
    a = ap.parse_args()

    pool = load_pool(a.pool)
    pairmap, npairs = load_pairs(a.pairs)
    print(f"pool={len(pool)}  pairs={npairs}", file=sys.stderr)
    assert len(pairmap) == npairs, "duplicate sample id in pairs file"

    w = Worker(pool, pairmap, a)
    n = 0
    for rec in iter_reads(a.input):
        if a.max_reads and n >= a.max_reads:
            break
        n += 1
        w.handle(rec)
        if n % 100000 == 0:
            print(f"  {n} reads processed", file=sys.stderr)
    w.flush()

    n_input = n                      # all reads seen in the input
    n_rq_drop = w.counts["rq_filtered"]
    kept = n_input - n_rq_drop       # reads that passed the rq filter
    assigned = kept - w.counts["unmapped"] - w.counts["self_pairs"] - w.counts["NOSKIP"]
    print(f"input reads:   {n_input}")
    if n_rq_drop:
        print(f"rq<{a.min_rq} dropped: {n_rq_drop:>10} ({n_rq_drop/max(n_input,1):.1%})")
        print(f"kept reads:    {kept:>10}")
    print(f"assigned:      {assigned:>10} ({assigned/max(kept,1):.1%} of kept)")
    print(f"self_pairs:    {w.counts['self_pairs']:>10} ({w.counts['self_pairs']/max(kept,1):.1%} of kept)")
    print(f"unmapped:      {w.counts['unmapped']:>10} ({w.counts['unmapped']/max(kept,1):.1%} of kept)")
    if w.counts["NOSKIP"]:
        print(f"no-seq skipped: {w.counts['NOSKIP']:>8}")

    # per-sample summary, sorted by sample id
    samples = sorted(w.pairmap.values())
    with open(os.path.join(a.outdir, "summary.tsv"), "w") as fh:
        fh.write("sample_id\tn_reads\n")
        for sid in samples:
            fh.write(f"{sid}\t{w.counts.get(sid, 0)}\n")
    print(f"wrote {os.path.join(a.outdir, 'summary.tsv')}")


if __name__ == "__main__":
    main()
