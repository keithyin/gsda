import re, glob, os
import xml.etree.ElementTree as ET
from collections import Counter

# ---- reference parse (as validated) ----
std = open("/data1/ccs_data/20260813-corn-ssr-analysis/corn-ssr-standard.md", encoding='utf-8').read()
locus_names = dict(re.findall(r'\|\s*(PM\d+)\s*\|\s*([\w.]+)\s*\|', std))

def docx_to_lines(path):
    tree = ET.parse(path)
    ns = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
    return [''.join(t.text or '' for t in p.iter(f'{{{ns}}}t')) for p in tree.getroot().iter(f'{{{ns}}}p')]

refs = {}
for path in glob.glob('/tmp/corn_ssr_extract/*/word/document.xml'):
    lines = docx_to_lines(path)
    cur, parts = None, []
    for line in lines:
        line = line.strip()
        m = re.match(r'PM(\d+)[：:]\s*', line)
        if m:
            num = int(m.group(1)); key = f"PM{num:02d}"
            if cur: refs[cur] = ''.join(parts)
            cur = key
            rest = line[m.end():]
            name = locus_names.get(key)
            if name and rest.startswith(name): rest = rest[len(name):]
            parts = [rest]
        else:
            if cur: parts.append(line)
    if cur: refs[cur] = ''.join(parts)
clean = lambda s: re.sub(r'[^ACGT]', '', s.upper())
refs = {k: clean(v) for k, v in refs.items() if 700 <= len(clean(v)) <= 2500}

def find_ssr(seq):
    best = None
    for u in range(2, 7):
        for start in range(len(seq)-u+1):
            unit = seq[start:start+u]
            j, c = start, 0
            while j+u <= len(seq) and seq[j:j+u] == unit:
                c += 1; j += u
            if c >= 4 and (best is None or c > best[1]):
                best = (unit, c, start)
    return best

ssr_info = {}
for k, seq in refs.items():
    b = find_ssr(seq)
    if b:
        unit, cnt, st = b
        end = st + len(unit)*cnt
        L = 22
        left = seq[st-L:st] if st >= L else None
        right = seq[end:end+L] if end+L <= len(seq) else None
        if left and right:
            ssr_info[k] = (unit, left, right, cnt)

print(f"可锚定位点数: {len(ssr_info)}")
rc = lambda s: s.translate(str.maketrans('ACGT','TGCA'))[::-1]

def count_repeats(seg, motif):
    """count longest consecutive motif repeats in seg"""
    best = 0
    m = motif
    for rot in range(len(m)):
        mm = m[rot:] + m[:rot]
        for start in range(len(seg)):
            j, c = start, 0
            while j+len(mm) <= len(seg) and seg[j:j+len(mm)] == mm:
                c += 1; j += len(mm)
            if c > best: best = c
    return best

def assign_and_count(read):
    """return (locus, cn) or None"""
    strands = [read, rc(read)]
    for strand in strands:
        for locus, (motif, left, right, _refcn) in ssr_info.items():
            li = strand.find(left)
            if li < 0: continue
            # right flank must be after SSR region
            seg = strand[li+len(left):]
            ri = seg.find(right)
            if ri < 0: continue
            if ri < 2 or ri > 200: continue   # SSR region sanity
            ssr_seg = seg[:ri]
            n = count_repeats(ssr_seg, motif)
            return locus, n
    return None

# ---- run on pilot samples ----
fq_dir = "/data1/ccs_data/20260813-corn-ssr-analysis/20260702_251201Y0001_Run0003/Group_0/barcodes_reads_fastq_amplicon"
pilot = ["Group_0_Adaptor-barcode201-0_26070201422997-MP19548-1-1",
         "Group_0_Adaptor-barcode202-0_26070201422997-MP19548-1-2",
         "Group_0_Adaptor-barcode203-0_26070201422997-MP19548-1-3",
         "Group_0_Adaptor-barcode206-0_26070201422997-MP19548-1-6"]

all_results = {}
for name in pilot:
    seqs = []
    with open(os.path.join(fq_dir, name + ".fastq")) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1: seqs.append(line.strip())
    per_locus = {}
    unmapped = 0
    for s in seqs:
        res = assign_and_count(s)
        if res is None:
            unmapped += 1
            continue
        loc, cn = res
        per_locus.setdefault(loc, []).append(cn)
    all_results[name.split('-1-')[1]] = (len(seqs), unmapped, per_locus)
    print(f"\n=== {name} : {len(seqs)} reads, 未归属 {unmapped} ===")

# ---- summarize per sample × locus ----
print("\n\n========== 汇总：每样品×位点 copy number 一致性 ==========")
print(f"{'样品':<5}{'位点':<7}{'基序':<10}{'深度':<6}{'主CN':<6}{'主支持%':<8}{'±1%':<7}{'≥±2%':<7}{'主峰清晰'}")
for smp, (nreads, unmapped, per_locus) in sorted(all_results.items()):
    for loc, cns in sorted(per_locus.items(), key=lambda x: int(x[0][2:])):
        c = Counter(cns); total = len(cns)
        if total < 5: continue
        mode, mode_n = c.most_common(1)[0]
        s1 = c.get(mode+1,0)+c.get(mode-1,0)
        s2 = sum(v for k,v in c.items() if k not in (mode, mode+1, mode-1))
        motif = ssr_info[loc][0]
        clear = '✅' if mode_n/total >= 0.6 else ('⚠️' if mode_n/total >= 0.4 else '❌')
        print(f"{smp:<5}{loc:<7}{motif:<10}{total:<6}{mode:<6}{mode_n/total*100:<8.1f}{s1/total*100:<7.1f}{s2/total*100:<7.1f}{clear}")
