import subprocess
import pathlib

import pysam
from tqdm import tqdm
import argparse
from glob import glob


class RecordInfo:
    def __init__(self, name, length, rq):
        self.name = name
        self.length = length
        self.rq = rq


def read_bam_file(bam_file: str, rq_thr: float):
    infos = {}
    with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False, threads=40) as in_bam:
        for record in tqdm(in_bam.fetch(until_eof=True), desc=f"reading {bam_file}"):
            rq = 1.0
            if record.has_tag("rq"):
                rq = record.get_tag("rq")

            if rq < rq_thr:
                continue

            infos[record.query_name] = RecordInfo(
                record.query_name, length=len(record.query_sequence), rq=rq)

    return infos


def gff_reader(fname):
    infos = {}
    with open(fname, mode="r", encoding="utf8") as in_file:
        for line in tqdm(in_file, desc=f"reading {fname}"):
            line = line.strip()
            if line.startswith("#"):
                continue
            items = line.split("\t")
            key = items[0]
            start = int(items[3])
            end = int(items[4])
            length = end - start + 1
            infos.setdefault(key, 0)
            infos[key] += length
    return infos


def main():

    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("bam_file", nargs="+")
    parser.add_argument('--rq-thr', type=float, default=0.95, dest="rq_thr")
    args = parser.parse_args()
    
    bam_files = []
    for f_pat in args.bam_file:
        bam_files.extend(glob(f_pat))
    print(f"processing bam files: {bam_files}")
    
    for bam_file in bam_files:

        record_infos = read_bam_file(bam_file=bam_file, rq_thr=args.rq_thr)

        p = pathlib.Path(bam_file)
        fasta_file = p.parent.joinpath(f"{p.stem}.fasta")

        cmd = f"samtools fasta {bam_file} > {fasta_file}"
        print(f"running {cmd}")

        subprocess.check_call(cmd, shell=True)

        cmd = f"samtools faidx {fasta_file}"
        print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        gff_file = p.parent.joinpath(f"{p.stem}.gff")
        cmd = f"tr-finder called {fasta_file} --unitAndRepeats 1-4,2-3,3-3 -o {gff_file}"
        print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        gff_infos = gff_reader(gff_file)

        tot_len = 0
        tr_len = 0
        for record_key, record_info in tqdm(record_infos.items(), desc=f"counting ..."):
            tot_len += record_info.length
            tr_len += gff_infos.get(record_key, 0)

        print(f"tr_ratio:{(tr_len / tot_len) * 100:.2f} %")


if __name__ == "__main__":
    main()
