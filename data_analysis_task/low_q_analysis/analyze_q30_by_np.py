#!/usr/bin/env python3
"""
Analyze why Q30 is low for np>=10 reads in SMC all-reads consensus BAM.

Reads `samtools view` TSV, groups by the `np` (number of passes) tag, and for
each np bucket reports base counts and Q30/Q20 fractions over the sequenced
(call) region. Also emits a per-position Q30 profile and a rq-vs-Q30 cross
tabulation to separate "read-level model low-confidence" from "read-level good
but position-driven 3'/5' degradation".
"""
import sys, collections, subprocess, statistics, json, argparse

Q30 = 30

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--min-np", type=int, default=10, help="only include reads with np >= this")
    ap.add_argument("--max-lines", type=int, default=0, help="debug: only process first N reads (0=all)")
    return ap.parse_args()

def base_q(ch):
    o = ord(ch)
    if ch in "=-":
        return None      # N / gap
    if o < 33:           # below '!' (33) -> not a real quality
        return None
    return o - 33

def main():
    a = parse_args()
    cmd = ["samtools", "view", a.bam]
    n_iter = a.max_lines if a.max_lines > 0 else None

    # per-np aggregates
    np_stats = collections.defaultdict(lambda: dict(n=0, len=0, Q30=0, Q20=0, bases=0, rq_list=[]))
    pos_profile = collections.defaultdict(lambda: [0,0])  # pos -> [bases, q30bases]
    len_profile = collections.defaultdict(lambda: [0,0])   # len->[bases,q30]
    rq_buckets = collections.defaultdict(lambda: [0,0])    # rqbucket->[bases,q30]
    maxpos = 0
    total_reads = 0
    included_reads = 0

    for i, line in enumerate(subprocess.run(cmd, stdout=subprocess.PIPE).stdout):
        if n_iter is not None and i >= n_iter:
            break
        try:
            fields = line.rstrip(b"\n").decode().split("\t")
        except Exception:
            continue
        if len(fields) < 11:
            continue
        total_reads += 1
        seq = fields[9]
        qual = fields[10]
        tags = fields[11:]

        np_val = None
        rq_val = None
        for t in tags:
            if t.startswith("np:"):
                np_val = int(t.split(":")[2])
            elif t.startswith("rq:"):
                rq_val = float(t.split(":")[2])
        if a.min_np > 0 and (np_val is None or np_val < a.min_np):
            continue
        included_reads += 1

        s = np_stats[np_val if np_val is not None else -1]
        s["n"] += 1
        L = len(qual)
        s["len"] += L
        maxpos = max(maxpos, L)
        if rq_val is not None:
            s["rq_list"].append(rq_val)
        for ci, ch in enumerate(qual):
            q = base_q(ch)
            if q is None:
                continue
            s["bases"] += 1
            if q >= Q30:
                s["Q30"] += 1
            if q >= 20:
                s["Q20"] += 1
            # position profile (cap to avoid huge arrays, but SMC reads are long)
            pos_profile[ci][0] += 1
            if q >= Q30:
                pos_profile[ci][1] += 1
            len_profile[L][0] += 1
            if q >= Q30:
                len_profile[L][1] += 1
            # rq bucket (5% bins)
            if rq_val is not None:
                b = int(rq_val * 100 // 5)
                rq_buckets[b][0] += 1
                if q >= Q30:
                    rq_buckets[b][1] += 1

    # ---- report ----
    print("\n==== SUMMARY (np >= %d) ====" % a.min_np)
    print("total reads in BAM: %d   included np>=%d reads: %d" % (total_reads, a.min_np, included_reads))
    if included_reads == 0:
        print("No reads matched np>=%d — check that the np tag exists." % a.min_np)
        return

    totbases = sum(s["bases"] for s in np_stats.values())
    totq30 = sum(s["Q30"] for s in np_stats.values())
    all_len = [l for s in np_stats.values() for _ in range(s['n']) for l in [s['len']/s['n']]]
    medlen = statistics.median([s['len']/s['n'] for s in np_stats.values() if s['n']])
    print("OVERALL: bases=%d  Q30=%.2f%%  (median mean read len ~ %.0f bp)" %
          (totbases, 100.0*totq30/max(totbases,1), statistics.median([s['len']/s['n'] for s in np_stats.values()])))

    print("\n==== Q30 vs np ====")
    print("np     reads     mean_len   bases      Q20%     Q30%     mean_rq")
    for npv in sorted(np_stats, key=lambda x: (x is None, x)):
        s = np_stats[npv]
        if s["bases"] == 0:
            continue
        mean_rq = (sum(s["rq_list"])/len(s["rq_list"])) if s["rq_list"] else float("nan")
        print("%-4s %8d %10.1f %10d %7.2f %7.2f %9.4f" % (
            npv, s["n"], s["len"]/s["n"], s["bases"],
            100.0*s["Q20"]/s["bases"], 100.0*s["Q30"]/s["bases"], mean_rq))

    print("\n==== Q30 vs mean read length (np>=%d pooled) ====" % a.min_np)
    print("len        bases     Q30%")
    for L in sorted(len_profile, key=lambda x:-len_profile[x][0])[:15]:
        b, q = len_profile[L]
        print("%6d %10d %7.2f" % (L, b, 100.0*q/b))

    print("\n==== Q30 vs rq bucket (5%% bins, np>=%d pooled) ====" % a.min_np)
    print("rq_bucket  lo-hi    bases     Q30%")
    for b in sorted(rq_buckets):
        lo = b*0.05
        hi = lo + 0.05
        bb, q = rq_buckets[b]
        if bb == 0: continue
        print("%8.3f   bases=%9d  Q30=%.2f" % (lo, bb, 100.0*q/bb))

    # ---- per-position profile (downsampled, first 40 positions + tail) ----
    print("\n==== Q30 profile by position (first 30 & last 20 positions, sampled) ====")
    pos_sorted = [p for p in pos_profile if pos_profile[p][0] > 0]
    if pos_sorted:
        def bucketed(idx, step):
            return [i*step+i2 for i2 in [0]]
        pos_sorted.sort()
        sample = []
        import math
        if maxpos > 80:
            step = max(1, maxpos//80)
            sample = list(range(0, maxpos, step))
        else:
            sample = list(range(0, maxpos))
        print("pos    bases     Q30%")
        for p in sample[:80]:
            b, q = pos_profile.get(p, (0,0))
            if b == 0:
                continue
            print("%6d %9d %7.2f" % (p, b, 100.0*q/b))
        # tail positions
        print("  ... (tail) ...")
        for p in sorted(pos_sorted, reverse=True)[:12]:
            b, q = pos_profile[p]
            print("%6d %9d %7.2f" % (p, b, 100.0*q/b))

if __name__ == "__main__":
    main()
