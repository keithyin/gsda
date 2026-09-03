"""按 demuxed FASTQ 的 barcode 分组，把合并的 SMC-all-reads BAM 拆成每个 barcode 一个 BAM。

逻辑：demuxed 目录下每个 `Adaptor-barcode<NNN>-<k>.fastq` 对应一个 barcode，其 FASTQ
header qname 有 3 段（`Run/ID/pos-pos`），而合并 BAM 的 qname 有 2 段（`Run/ID`）。
两者的对应关系为 `"/".join(fastq_qname.split("/")[:2]) == bam_qname`（去掉第 3 段）。
据此把 BAM 里每条 read 按 barcode 路由到对应的输出 BAM（文件名 = 对应 FASTQ 换成 .bam）。
"""

import argparse
import os
import pathlib
from multiprocessing import cpu_count
from typing import Dict, List

import pysam
from tqdm import tqdm


def build_barcode_index(demuxed_dir: str) -> Dict[str, List[str]]:
    """遍历 demuxed 目录，建立 bam_qname(prefix) -> [输出 bam 名] 的映射。

    输出 bam 名 = 对应 FASTQ 文件名把 .fastq 换成 .bam。
    同一 prefix 可能出现在多个 barcode（防御性处理），映射值用 list 收集。
    """
    demuxed_dir = pathlib.Path(demuxed_dir)
    prefix_to_outputs: Dict[str, set] = {}
    for fastq_path in sorted(demuxed_dir.glob("*.fastq")):
        out_name = fastq_path.name[: -len(".fastq")] + ".bam"
        with open(fastq_path) as f:
            for line in f:
                if line[0] != "@":
                    continue
                qname = line[1:].strip()
                prefix = "/".join(qname.split("/")[:2])
                prefix_to_outputs.setdefault(prefix, set()).add(out_name)

    return {p: sorted(s) for p, s in prefix_to_outputs.items()}


def split_bam_by_barcode(demuxed_dir: str, bam: str, output_dir: str, threads: int) -> None:
    index = build_barcode_index(demuxed_dir)
    print(f"📂 demuxed 目录: {demuxed_dir}")
    print(f"📂 待拆分的 BAM : {bam}")
    print(f"📂 输出目录     : {output_dir}")
    print(f"🔢 建立的 prefix→barcode 映射条目数: {len(index):,}")

    os.makedirs(output_dir, exist_ok=True)

    writers: Dict[str, pysam.AlignmentFile] = {}
    counts: Dict[str, int] = {}
    n_total = 0
    n_routed = 0
    n_unrouted = 0

    with pysam.AlignmentFile(
        filename=bam, mode="rb", check_sq=False, threads=threads
    ) as reader:
        for record in tqdm(reader.fetch(until_eof=True), desc="routing reads"):
            n_total += 1
            out_names = index.get(record.query_name)
            if not out_names:
                n_unrouted += 1
                continue
            for out_name in out_names:
                writer = writers.get(out_name)
                if writer is None:
                    writer = pysam.AlignmentFile(
                        os.path.join(output_dir, out_name),
                        "wb",
                        template=reader,
                    )
                    writers[out_name] = writer
                    counts[out_name] = 0
                writer.write(record)
                counts[out_name] += 1
                n_routed += 1

    for writer in writers.values():
        writer.close()
    writers.clear()

    print(f"\n📊 总计 reads      : {n_total:,}")
    print(f"📊 已路由 reads    : {n_routed:,}")
    print(f"📊 未路由 reads    : {n_unrouted:,}（qname 不在任何 barcode 中）")
    print(f"📊 输出 BAM 数量   : {len(counts)}")
    for out_name in sorted(counts):
        print(f"   {out_name}: {counts[out_name]:,}")
    print(f"✅ 拆分完成，输出目录: {output_dir}")


def main_cli():
    parser = argparse.ArgumentParser(
        prog="split-bam-by-barcode",
        description="按 demuxed FASTQ 的 barcode 分组拆分合并的 SMC-all-reads BAM",
    )
    parser.add_argument("demuxed_dir", type=str, help="demuxed FASTQ 所在目录")
    parser.add_argument("bam", type=str, help="待拆分的合并 BAM 路径")
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="每个 barcode 的 BAM 输出目录（不存在会自动创建）",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=cpu_count() // 2,
        help=f"pysam 读取线程数。Default: {cpu_count() // 2}",
    )
    args = parser.parse_args()
    split_bam_by_barcode(
        demuxed_dir=args.demuxed_dir,
        bam=args.bam,
        output_dir=args.output_dir,
        threads=args.threads,
    )


if __name__ == "__main__":
    main_cli()
