#!/usr/bin/env python3
"""
Analysis: why Q30 is "low" for np>=10 reads
Target run : /data1/ccs_data/low-q-analysis/20260819_250214YJ006_Run0002
Reference  : /data1/ccs_data/16s-low-q30/20251126_250214YJ006_Run0004
             (same YJ006 sample class, same GSEQ500 instrument, Nov 2025)
"""
import subprocess, numpy as np, collections, argparse

def analyze(bam, min_np=10, sample_cap=0):
    per = collections.defaultdict(lambda: dict(n=0, bases=0, q30=0, q20=0,
                                                all30=0, m95=0, m90=0))
    rq_ct = collections.Counter()
    ends5 = collections.Counter(); mid = collections.Counter(); ends3 = collections.Counter()
    tot = 0; included = 0

    p = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE)
    for line in p.stdout:
        tot += 1
        if sample_cap and tot > sample_cap: break
        parts = line.rstrip(b"\n").split(b"\t")
        if len(parts) < 11: continue
        npv = -1; rq = 0.0
        for t in parts[11:]:
            if t[:5] == b"np:i:": npv = int(t[5:])
            elif t[:5] == b"rq:f:": rq = float(t[5:])
        if npv < min_np: continue
        q = parts[10]; L = len(q)
        if L == 0: continue
        code = np.frombuffer(q[:min(L, 100000)], dtype=np.uint8).astype(np.int16) - 33
        code = code[(code >= 0) & (code <= 99)]
        if code.size == 0: continue
        included += 1
        s = per[npv]
        s['n'] += 1
        s['bases'] += code.size
        s['q30'] += (code >= 30).sum()
        s['q20'] += (code >= 20).sum()
        frac30 = (code >= 30).mean()
        if frac30 == 1.0: s['all30'] += 1
        if frac30 >= 0.95: s['m95'] += 1
        if frac30 >= 0.90: s['m90'] += 1
        rq_ct[int(rq * 100)] += 1
        if L >= 600:
            def rfq30(s_, e_):
                seg = code[s_:e_]; seg = seg[(seg >= 0) & (seg <= 99)]
                if seg.size: return float((seg >= 30).mean())
                return None
            r5 = rfq30(0, 100); rm = rfq30(L // 2 - 200, L // 2 + 200); r3 = rfq30(L - 100, L)
            if r5 is not None: ends5[r5] += 1
            if rm is not None: mid[rm] += 1
            if r3 is not None: ends3[r3] += 1
    p.stdout.close(); p.wait()

    print("=" * 78)
    print(f"  TARGET run analysis (np>={min_np})")
    print(f"  Total reads scanned: {tot};  included: {included}")
    print("=" * 78)

    # --- base-level ---
    totb = sum(s['bases'] for s in per.values())
    totq30 = sum(s['q30'] for s in per.values())
    totq20 = sum(s['q20'] for s in per.values())
    tot_all30 = sum(s['all30'] for s in per.values())
    tot_incl = sum(s['n'] for s in per.values())
    print(f"\n  ## Base-level (per base, all reads pooled)")
    print(f"     base Q20 : {100*totq20/totb:.2f}%")
    print(f"     base Q30 : {100*totq30/totb:.2f}%")
    print(f"\n  ## Read-level QC (per read, all-or-nothing metrics)")
    print(f"     strict (every base Q30+): {100*tot_all30/tot_incl:6.2f}%  of reads pass")
    m95 = sum(s['m95'] for s in per.values())
    m90 = sum(s['m90'] for s in per.values())
    print(f"     lenient (>=95% Q30+)     : {100*m95/tot_incl:6.2f}%  of reads pass")
    print(f"     lenient (>=90% Q30+)     : {100*m90/tot_incl:6.2f}%  of reads pass")

    # --- per-np base Q30 ---
    print(f"\n  ## Base-level Q30% per np (target run)")
    print(f"  {'np':>4}  {'reads':>9}  {'bases':>11}  {'Q20%':>6}  {'Q30%':>6}  "
          f"{'allQ30':>7}  {'>=95%':>6}  {'>=90%':>6}")
    for n in sorted(per):
        s = per[n]
        print(f"  {n:>4}  {s['n']:>9d}  {s['bases']:>d}  "
              f"{100*s['q20']/s['bases']:>5.2f}  {100*s['q30']/s['bases']:>5.2f}  "
              f"{100*s['all30']/s['n']:>5.2f}  {100*s['m95']/s['n']:>5.2f}  "
              f"{100*s['m90']/s['n']:>5.2f}")

    # --- region profile ---
    if included and (sum(ends5.values()) + sum(mid.values()) + sum(ends3.values())) > 0:
        print(f"\n  ## Region Q30 distribution (reads with len>=600, "
              f"n={sum(mid.values())})")
        def reg_stats(name, d):
            if not d: return
            total = sum(d.values())
            mean = sum(k * v for k, v in d.items()) / total
            print(f"     {name:12s}: n={total:8d}   mean-Q30-fraction={mean:.4f}")
        reg_stats("5' 100bp", ends5)
        reg_stats("middle 400bp", mid)
        reg_stats("3' 100bp", ends3)

    # --- rq (read-quality model output) ---
    print(f"\n  ## Consensus-read quality score `rq` (target run)")
    tot_r = sum(rq_ct.values())
    if tot_r:
        for k in sorted(rq_ct):
            print(f"     rq={k / 100:.2f}   reads={rq_ct[k]:>8d}   "
                  f"{100 * rq_ct[k] / tot_r:6.2f}%")
    return per

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--baseline", default="")
    ap.add_argument("--min-np", type=int, default=10)
    a = ap.parse_args()
    analyze(a.target, a.min_np)
    if a.baseline:
        print("\n" + "=" * 78)
        print(f"  BASELINE run analysis (np>={a.min_np}) — for comparison")
        analyze(a.baseline, a.min_np, sample_cap=300000)
