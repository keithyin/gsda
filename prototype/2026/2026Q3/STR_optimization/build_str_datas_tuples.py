#!/usr/bin/env python3
"""为 str-optimization 第二/三/四批构造 (adapter, smc, ref, bed, vcf) 5 元组。

布局约定（见 .claude/skills/smc_barcode_split 产物）：
  <batch>/<RUN>/barcode-partitioned/
      <RUN>_called-demuxed-v4-.smc_all_reads.BarcodeNN.bam   # SMC 共识
      <RUN>_called-demuxed.BarcodeNN.bam                     # subread / adapter
参考取自各批 `*一代测序/.../merged_output/`：
  第二批  STR<plasmid>.fa
  第三批  <date>STR<plasmid>.intersect.fa   (缺 32-1 → 回退 barcode_ref_align_report/refs/)
  第四批  260827STR-<plasmid>.intersect.fa  (缺 93-1/95-3 → 回退 report/refs/)
barcode 对应关系：
  第二/三批  plasmid 列 + barcode 列 `24标签-N` -> Barcode(N)
  第四批     plasmid  barcode 列已是裸整数
"""

import os
import re
import sys

ROOT = "/data1/ccs_data/str-optimization"

# barcode 少于该读段数的是噪声孔位（本三批实测都是 1–3 条），不出元组
MIN_READS = 100

BATCHES = [
    # (批目录, run 目录, mapping 文件, merged_output 目录, refs 回退目录, barcode 命名宽度)
    ("second-batch-of-data", "20260805_250804Y0004_Run0001",
     "plasmid_2_barcode.tsv",
     "STR第二批一代测序/STR第二批一代测序/merged_output",
     None, 2),
    ("second-batch-of-data", "20260805_250804Y0004_Run0002",
     "plasmid_2_barcode.tsv",
     "STR第二批一代测序/STR第二批一代测序/merged_output",
     None, 2),
    ("third-batch-of-data", "20260815_250804Y0004_Run0002",
     "plasmid_name_2_barcode.tsv",
     "STR第三批一代测序/STR第三批一代测序/merged_output",
     "barcode_ref_align_report/refs", 2),
    ("fourth-batch-of-data", "20260831_250302Y0001_Run0001",
     "plasmid_name_2_barcode.tsv",
     "STR第四批一代测序/merged_output",
     "report/refs", 2),
]


def norm_barcode(raw):
    """`24标签-7` / `7` / `Barcode07` -> int 7"""
    m = re.search(r'(\d+)\s*$', raw.strip())
    if not m:
        raise ValueError(f"无法解析 barcode: {raw!r}")
    return int(m.group(1))


def read_mapping(path, run=None):
    """返回 [(plasmid, barcode_int)]；第二批的 tsv 带 RUN号 列，按 run 过滤。"""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        run_col = header.index("RUN号") if "RUN号" in header else None
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if run_col is not None and cols[run_col] != run:
                continue
            rows.append((cols[header.index("plasmid")], norm_barcode(cols[header.index("barcode")])))
    return rows


def find_ref(mo_dir, batch_dir, run, fallback_dir, plasmid):
    """在 merged_output 里找该 plasmid 的参考；再退到 <run>/report*/refs/<plasmid>.fa。"""
    pat = re.compile(r'(?:^|STR[-_?])?' + re.escape(plasmid) + r'\.intersect\.fa$|^STR' +
                     re.escape(plasmid) + r'\.fa$')
    for name in sorted(os.listdir(os.path.join(ROOT, batch_dir, mo_dir))):
        if not pat.search(name):
            continue
        # 防止 3N-1 之类的名字被 -1 / 2-1 之类误匹配：STR 与 plasmid 之间不能有别的数字
        stem = re.sub(r'(\.intersect)?\.fa$', '', name)
        stem = re.sub(r'^\d{6}', '', stem)           # 去日期前缀 260805 / 260731
        stem = re.sub(r'^STR-?', '', stem)            # 去 STR / STR-
        if stem == plasmid:
            return os.path.join(ROOT, batch_dir, mo_dir, name), None
    if fallback_dir:
        cand = os.path.join(ROOT, batch_dir, run, fallback_dir, plasmid + ".fa")
        if os.path.exists(cand):
            return cand, "report/refs 回退"
    return None, None


def bam_count(path):
    import subprocess
    return int(subprocess.check_output(["samtools", "view", "-c", path]).decode().strip())


def unmapped_runs():
    """mapping 文件里出现、但数据不在 /data1 上的 run（第二批 Run0007 在 72 服务器）。"""
    notes = []
    seen = {(b, r) for b, r, *_ in BATCHES}
    for batch, mapfile in {("second-batch-of-data", "plasmid_2_barcode.tsv"),
                           ("third-batch-of-data", "plasmid_name_2_barcode.tsv"),
                           ("fourth-batch-of-data", "plasmid_name_2_barcode.tsv")}:
        path = os.path.join(ROOT, batch, mapfile)
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if "RUN号" not in header:
                continue
            col = header.index("RUN号")
            runs = {line.rstrip("\n").split("\t")[col] for line in fh if line.strip()}
            for run in sorted(runs - {r for b, r in seen if b == batch}):
                why = "run 目录不在 /data1" if not os.path.exists(os.path.join(ROOT, batch, run)) else "run 目录存在但未纳入"
                notes.append((batch, run, "-", "-", why))
    return notes


def main():
    tuples, skipped = [], unmapped_runs()
    for batch, run, mapfile, mo_dir, fallback_dir, width in BATCHES:
        run_dir = os.path.join(ROOT, batch, run)
        part = os.path.join(run_dir, "barcode-partitioned")
        for plasmid, bc in read_mapping(os.path.join(ROOT, batch, mapfile), run=run):
            bcn = f"Barcode{bc:0{width}d}"
            smc = os.path.join(part, f"{run}_called-demuxed-v4-.smc_all_reads.{bcn}.bam")
            adapter = os.path.join(part, f"{run}_called-demuxed.{bcn}.bam")
            missing = [p for p in (smc, adapter) if not os.path.exists(p)]
            if missing:
                skipped.append((batch, run, bcn, plasmid, "BAM 缺失: " + ", ".join(missing)))
                continue
            ref, note = find_ref(mo_dir, batch, run, fallback_dir, plasmid)
            if ref is None:
                skipped.append((batch, run, bcn, plasmid, "merged_output 无一代参考"))
                continue
            reads = bam_count(smc)
            if reads < MIN_READS:
                skipped.append((batch, run, bcn, plasmid, f"噪声孔位 ({reads} reads)"))
                continue
            tuples.append((batch, run, bcn, plasmid, adapter, smc, ref, reads, note))

    print("datas = [")
    cur = None
    for batch, run, bcn, plasmid, adapter, smc, ref, reads, note in tuples:
        if (batch, run) != cur:
            cur = (batch, run)
            print(f"    # ---- {batch} / {run} / {bcn}={plasmid} ({reads} reads)")
        tail = "  # ref 取自 report/refs" if note else ""
        print(f'        ("{adapter}",')
        print(f'         "{smc}", "{ref}", "None", "None"),{tail}')
    print("]\n")
    print("# 未生成（共 %d 条）:" % len(skipped))
    for batch, run, bcn, plasmid, why in skipped:
        print(f"#   {batch} {run} {bcn}={plasmid}: {why}")
    print(f"\n# 元组总数: {len(tuples)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
