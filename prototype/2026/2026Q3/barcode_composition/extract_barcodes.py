#!/usr/bin/env python
"""Extract the two-end (paired) barcodes from an smc_all_reads FASTA.

Read grammar observed in
  20260824_250302Y0001_Run0001_called_demuxed.smc_all_reads.q20.fa

    (AAGGTTAA + BC_left + CAGCACCA) x2   ... insert ...   x2 (TGGTGCTG + BC_right + TTAACCTT)

so every read carries a *pair* of 24-nt end barcodes, each present in two
tandem copies (one for each strand of the duplex). The right tag appears on
the read's own strand as revcomp of its canonical form, so we recover it by
running the same head-grammar scan over the reverse complement.

Both tandem copies of a tag are used to build a per-read consensus, which
absorbs most remaining sequencing error.
"""
import sys
import argparse
from collections import Counter, defaultdict

ANCHOR = "AAGGTTAA"          # 5' anchor of a tag
SEP = "CAGCACCA"             # 3' separator of a tag
BC_LEN = 24
TAG_LEN = 8 + BC_LEN + 8     # 40
SCAN = 200                   # how far into the read to look for a tag
MAX_MM = 2                   # total mismatches (anchor + separator) to accept


def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def mm(a, b):
    """Hamming distance; unequal length => reject."""
    if len(a) != len(b):
        return 99
    return sum(1 for x, y in zip(a, b) if x != y)


def extract_tag(seq):
    """Return (barcode, score) for the first tag in seq[:SCAN], or (None, None).

    Picks the offset minimising mismatches against the anchor+separator
    fingerprint. The barcode is the 24 nt sitting between them.
    """
    best_i, best_s = -1, 99
    lim = min(SCAN, len(seq) - TAG_LEN) + 1
    for i in range(lim):
        s = mm(seq[i:i + 8], ANCHOR)
        if s > MAX_MM:
            continue
        s += mm(seq[i + 8 + BC_LEN:i + TAG_LEN], SEP)
        if s < best_s:
            best_i, best_s = i, s
            if s == 0:
                break
    if best_i < 0 or best_s > MAX_MM:
        return None, None
    return seq[best_i + 8:best_i + 8 + BC_LEN], best_s


def iter_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].rstrip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def consensus(a, b):
    """Majority per position of the two tandem copies; ties -> the first copy."""
    if a is None:
        return b
    if b is None:
        return a
    return "".join(x if x == y else x for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("--max-reads", type=int, default=0, help="0 = all")
    ap.add_argument("--out-prefix", default="/tmp/jinsirui")
    a = ap.parse_args()

    left = Counter()
    right = Counter()
    pair = Counter()
    n = n_both = n_l = n_r = 0

    for name, seq in iter_fasta(a.fasta):
        if a.max_reads and n >= a.max_reads:
            break
        n += 1

        # --- copy 1 = tandem copy starting near the head, copy 2 follows it
        b1, s1 = extract_tag(seq)
        b1c = b1
        if b1 is not None:
            # the second tandem copy starts 40 nt later (allow +-2 for indels)
            b2, s2 = extract_tag(seq[40:44 + TAG_LEN + 2])
            b1c = consensus(b1, b2)
        else:
            b2, s2 = None, None

        # --- right tag: same grammar on the reverse complement
        r = rc(seq)
        c1, t1 = extract_tag(r)
        c1c = c1
        if c1 is not None:
            c2, t2 = extract_tag(r[40:44 + TAG_LEN + 2])
            c1c = consensus(c1, c2)
        else:
            c2, t2 = None, None

        if b1c:
            n_l += 1
            left[b1c] += 1
        if c1c:
            n_r += 1
            right[c1c] += 1
        if b1c and c1c:
            n_both += 1
            pair[(b1c, c1c)] += 1

    print(f"reads={n}  left_ok={n_l} ({n_l/n:.1%})  right_ok={n_r} ({n_r/n:.1%})  "
          f"both_ok={n_both} ({n_both/n:.1%})", file=sys.stderr)
    print(f"distinct_left={len(left)}  distinct_right={len(right)}  "
          f"distinct_pairs={len(pair)}", file=sys.stderr)

    for tag, ctr in (("left", left), ("right", right)):
        with open(f"{a.out_prefix}_{tag}_raw.tsv", "w") as fh:
            fh.write("count\tbarcode\n")
            for bc, c in ctr.most_common():
                fh.write(f"{c}\t{bc}\n")
    with open(f"{a.out_prefix}_pair_raw.tsv", "w") as fh:
        fh.write("count\tbarcode_left\tbarcode_right\n")
        for (l, r_), c in pair.most_common():
            fh.write(f"{c}\t{l}\t{r_}\n")


if __name__ == "__main__":
    main()
