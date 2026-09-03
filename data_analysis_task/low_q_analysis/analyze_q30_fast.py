#!/usr/bin/env python3
"""
Why is Q30 low for np>=10 reads in the SMC all-reads consensus BAM?

Stream `samtools view` TSV; keep only reads whose `np` (number of passes) tag is
>= threshold; then compute Q30/Q20, a per-position Q30 profile, length-vs-Q30 and
read-quality(rq)-vs-Q30, vectorized with numpy over the qual byte array.

Qscore mapping: phred = qualchar - 33 (ASCII '!' = 33 = Q0). Q20 iff code>=20,
Q30 iff code>=30. N/'-' and any code outside [0,99] are excluded (not a real base).
"""
import subprocess, numpy as np, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--min-np", type=int, default=10)
    ap.add_argument("--maxlen", type=int, default=3000, help="truncate qual analysis to N bp")
    a = ap.parse_args()

    per_np = {}            # np -> [reads, bases, q30, q20]
    rq_hist = {}           # rq 5% bin -> [bases, q30]
    len_hist = {}          # len -> [bases, q30]
    posbase = np.zeros(a.maxlen, dtype=np.int64)
    pos30   = np.zeros(a.maxlen, dtype=np.int64)
    bases_total = q30sum = q20sum = 0
    nreads = total = 0

    p = subprocess.Popen(["samtools", "view", a.bam], stdout=subprocess.PIPE)
    for line in p.stdout:
        total += 1
        parts = line.rstrip(b"\n").split(b"\t")
        if len(parts) < 11:
            continue
        # --- np tag ---
        npv = -1
        rq = 0.0
        for tg in parts[11:]:
            if tg[:5] == b"np:i:":
                npv = int(tg[5:])
            elif tg[:5] == b"rq:f:":
                rq = float(tg[5:])
        if npv < a.min_np:
            continue
        # --- qual (field 11, index 10) ---
        qual = parts[10]
        if not qual:
            continue
        arr = np.frombuffer(qual, dtype=np.uint8).astype(np.int16).copy()
        if arr.size == 0:
            continue
        lenv = min(arr.size, a.maxlen)
        code = arr[:lenv] - 33
        valid = (code >= 0) & (code <= 99)
        nreads += 1
        b = int(valid.sum())
        if b == 0:
            continue
        q30 = int(np.count_nonzero(valid & (code >= 30)))
        q20 = int(np.count_nonzero(valid & (code >= 20)))
        bases_total += b; q30sum += q30; q20sum += q20
        d = per_np.setdefault(npv, [0, 0, 0, 0]); d[0] += 1; d[1] += b; d[2] += q30; d[3] += q20
        rd = rq_hist.setdefault(int(rq * 100 // 5), [0, 0]); rd[0] += b; rd[1] += q30
        ld = len_hist.setdefault(arr.size, [0, 0]); ld[0] += b; ld[1] += q30
        pos = np.arange(lenv, dtype=np.int16)
        posbase += np.bincount(pos[valid], minlength=a.maxlen)[:a.maxlen]
        pos30   += np.bincount(pos[valid & (code >= 30)], minlength=a.maxlen)[:a.maxlen]

    p.stdout.close(); p.wait()
    if nreads == 0:
        print("No reads matched np>=%d (scanned %d)" % (a.min_np, total)); return

    print("\n==== SUMMARY (np >= %d) ====" % a.min_np)
    print("reads: %d   total valid bases: %d" % (nreads, bases_total))
    print("OVERALL  Q20 = %.2f%%   Q30 = %.2f%%" % (100 * q20sum / bases_total, 100 * q30sum / bases_total))

    print("\n==== Q30 vs np ====")
    print("np     reads      bases      Q20%     Q30%")
    for npv in sorted(per_np):
        r, b, q30, q20 = per_np[npv]
        print("%-4d %9d %11d %7.2f %7.2f" % (npv, r, b, 100 * q20 / b, 100 * q30 / b))

    print("\n==== Q30 vs read length (top by bases) ====")
    print("len        bases     Q30%")
    for L in sorted(len_hist, key=lambda x: -len_hist[x][0])[:15]:
        b, q30 = len_hist[L]; print("%6d %10d %7.2f" % (L, b, 100 * q30 / b))

    print("\n==== Q30 vs read-level rq (5% bins) ====")
    print("rq_lo      bases     Q30%")
    for rb in sorted(rq_hist):
        b, q30 = rq_hist[rb]; print("%5.2f %10d %7.2f" % (rb * 0.05, b, 100 * q30 / b))

    print("\n==== Q30 by position (60 buckets over %d bp) ====" % a.maxlen)
    nb = 60
    edge = np.linspace(0, a.maxlen, nb + 1).astype(int)
    bp = posbase[1:]; qp = pos30[1:]
    for i in range(nb):
        s, e = edge[i], max(edge[i] + 1, edge[i + 1])
        bb = int(bp[s:e].sum()); qq = int(qp[s:e].sum())
        if bb > 0:
            print("pos %5d-%5d  Q30=%6.2f%%" % (s, e, 100 * qq / bb))

if __name__ == "__main__":
    main()
