import subprocess
import pathlib
import pysam
from tqdm import tqdm
import argparse
from glob import glob
import multiprocessing as mp


class RecordInfo:
    def __init__(self, name, length, rq):
        self.name = name
        self.length = length
        self.rq = rq


def read_bam_file(bam_file: str, rq_thr: float):
    infos = {}
    # 使用 pysam 读取 BAM 及其 rq 标签
    with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False, threads=mp.cpu_count()) as in_bam:
        for record in tqdm(in_bam.fetch(until_eof=True), desc=f"reading BAM: {bam_file}"):
            rq = 1.0
            if record.has_tag("rq"):
                rq = record.get_tag("rq")

            if rq < rq_thr:
                continue

            infos[record.query_name] = RecordInfo(
                record.query_name, length=len(record.query_sequence), rq=rq)
    return infos


def process_fastx_to_fasta(input_file: str, output_fasta: str):
    """
    使用 pysam.FastxFile 读取 FASTQ 或 FASTA，
    并提取长度信息，同时写出一个标准的 FASTA 文件。
    """
    infos = {}
    with pysam.FastxFile(input_file) as fx:
        out_fa = None
        if input_file.endswith(".fq") or input_file.endswith(".fastq"):
            out_fa = open(output_fasta, "w")

        for entry in tqdm(fx, desc=f"processing {input_file}"):
            # 记录元数据
            infos[entry.name] = RecordInfo(
                entry.name, length=len(entry.sequence), rq=1.0)
            # 写入 FASTA 格式
            if out_fa is not None:
                out_fa.write(f">{entry.name}\n{entry.sequence}\n")
        if out_fa is not None:
            out_fa.flush()
            out_fa.close()

    return infos


def gff_reader(fname):
    infos = {}
    if not pathlib.Path(fname).exists():
        return infos
    with open(fname, mode="r", encoding="utf8") as in_file:
        for line in tqdm(in_file, desc=f"reading GFF: {fname}"):
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            items = line.split("\t")
            if len(items) < 5:
                continue
            key = items[0]
            start = int(items[3])
            end = int(items[4])
            length = end - start + 1
            infos.setdefault(key, 0)
            infos[key] += length
    return infos


def main(args):
    input_files = []
    for f_pat in args.files:
        input_files.extend(glob(f_pat))

    if not input_files:
        print("No files found. Please check your input paths.")
        return

    print(f"Processing files: {input_files}")
    report_inner = ""

    for f_path in input_files:
        p = pathlib.Path(f_path)
        suffix = p.suffix.lower()
        fasta_file = p.parent.joinpath(f"{p.stem}.fasta")
        gff_file = p.parent.joinpath(f"{p.stem}.gff")

        # 1. 根据不同类型读取信息并统一生成/准备 FASTA
        if suffix == ".bam":
            record_infos = read_bam_file(bam_file=f_path, rq_thr=args.rq_thr)
            cmd_conv = f"samtools fasta {f_path} > {fasta_file}"
            print(f"running {cmd_conv}")
            subprocess.check_call(cmd_conv, shell=True)
        elif suffix in [".fastq", ".fq", ".fasta", ".fa"]:
            # 使用 pysam.FastxFile 处理，不再调用 shell 命令转换
            record_infos = process_fastx_to_fasta(f_path, str(fasta_file))
        else:
            print(f"Skipping unsupported file type: {suffix}")
            continue

        # 2. 索引 FASTA (tr-finder 通常需要 faidx)
        cmd_idx = f"samtools faidx {fasta_file}"
        print(f"running {cmd_idx}")
        subprocess.check_call(cmd_idx, shell=True)

        # 3. 运行 tr-finder
        cmd_tr = f"tr-finder called {fasta_file} --unitAndRepeats 1-4,2-3,3-3 -o {gff_file}"
        print(f"running {cmd_tr}")
        subprocess.check_call(cmd_tr, shell=True)

        # 4. 统计结果
        gff_infos = gff_reader(gff_file)

        tot_len = 0
        tr_len = 0
        for record_key, record_info in record_infos.items():
            tot_len += record_info.length
            tr_len += gff_infos.get(record_key, 0)

        ratio = (tr_len / tot_len * 100) if tot_len > 0 else 0
        res_line = f"File: {p.name} | Ratio: {ratio:.2f} %"
        report_inner += f"- {res_line}\n"

    report_str = f"""
================= HomoAndStr Report =================
{report_inner}
=====================================================
"""
    print(report_str)
    return report_str


def main_cli():
    parser = argparse.ArgumentParser(
        description="HomoAndStr Ratio Tool (Supports BAM/FASTQ/FASTA)")
    parser.add_argument("files", nargs="+",
                        help="Input files (.bam, .fq, .fastq, .fa, .fasta)")
    parser.add_argument('--rq-thr', type=float, default=0.95,
                        help="Read Quality threshold (BAM only)")
    args = parser.parse_args()
    main(args=args)


if __name__ == "__main__":
    main_cli()
