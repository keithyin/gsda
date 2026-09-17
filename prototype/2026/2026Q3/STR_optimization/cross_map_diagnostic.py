#!/usr/bin/env python3
"""
Which reference does each read actually prefer? (cross-plasmid diagnostic)

The barcode_ref_align report only ever shows a read against ITS OWN assigned
reference, so a well that got the wrong plasmid — or a reference that is the
wrong allele — just looks like "low identity / low alignedRatio". This re-maps
a sample's reads against ALL reference records at once and picks each read's
best target, so mis-assignment becomes visible.

Best-target rule: among alignments covering at least --min-ref-cov of the
target, take the highest per-base identity (ties broken by more matched bases).
Ranking on identity rather than raw match count matters here: the STR33-series
references are the same locus at different repeat counts, so a length-biased
score would always hand the win to the longest record.

Reads shorter than --min-len are skipped: they are the amplicon species and the
background full-plasmid species differ sharply (see ccs_vs_ref_stats.py), and
only the former carry enough target sequence to discriminate.

  python3 cross_map_diagnostic.py --assigned <plasmid>=<file.fastq> --refs merged.fa
"""
import argparse
import collections
import os
import subprocess
import sys
import tempfile


def sample_reads(path, min_len, max_n):
    out, n, seen = [], 0, 0
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline(); fh.readline()
            seen += 1
            if len(s) >= min_len:
                out.append((h[1:].split()[0], s))
                n += 1
                if n >= max_n:
                    break
    return out, seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assigned", action="append", required=True, metavar="PLASMID=FASTQ",
                    help="plasmid name (as a reference id) = its barcode fastq, repeatable")
    ap.add_argument("--refs", required=True, help="multi-record reference fasta")
    ap.add_argument("--prefix", default="STR", help="prefix stripped to match --assigned names")
    ap.add_argument("--min-len", type=int, default=800)
    ap.add_argument("--min-ref-cov", type=float, default=0.5)
    ap.add_argument("--max-reads", type=int, default=3000)
    ap.add_argument("--preset", default="map-hifi")
    args = ap.parse_args()

    for spec in args.assigned:
        pl, _, fq = spec.partition("=")
        reads, seen = sample_reads(fq, args.min_len, args.max_reads)
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
            q = tf.name
            for name, seq in reads:
                tf.write(f">{name}\n{seq}\n")
        r = subprocess.run(["minimap2", "-x", args.preset, "-c", "-t", "8",
                            args.refs, q], capture_output=True, text=True)
        os.unlink(q)
        if r.returncode != 0:
            print(f"{pl}: minimap2 failed -> {r.stderr.strip()[:200]}", file=sys.stderr)
            continue
        best = {}
        for line in r.stdout.splitlines():
            c = line.split("\t")
            qn, tn, tlen = c[0], c[5], int(c[6])
            mlen, alen = int(c[9]), int(c[10])
            if alen == 0 or alen / tlen < args.min_ref_cov:
                continue
            key = (mlen / alen, mlen)
            if qn not in best or key > best[qn][1]:
                best[qn] = (tn, key)
        cnt = collections.Counter(v[0] for v in best.values())
        assigned = args.prefix + pl
        wins = cnt.most_common(5)
        self_pct = 100 * cnt.get(assigned, 0) / len(best) if best else 0.0
        own_id = [v[1][0] for v in best.values() if v[0] == assigned]
        print(f"{pl:8s} reads>={args.min_len}: {len(reads):5d}  with>=50%-ref-hit: {len(best):5d}"
              f"  own-ref wins: {cnt.get(assigned,0):5d} ({self_pct:5.1f}%)"
              f"  own-ref identity p50: {sorted(own_id)[len(own_id)//2] if own_id else float('nan'):.4f}")
        for t, c in wins:
            print(f"           best-hit {t:12s} {c:5d} ({100*c/len(best):5.1f}%)")


if __name__ == "__main__":
    main()
