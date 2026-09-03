
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import argparse
import pysam
import multiprocessing as mp


def read_pos_from_vcf_file(fname, interested_names):
    positions = []
    with open(fname, mode="r") as f:
        for line in tqdm(f, desc=f"reading {fname}"):
            line = line.strip()
            if line.startswith("#"):
                continue

            items = line.split("\t")
            if items[0].strip() in interested_names:
                pos = int(items[1])
                positions.append(pos)

    return np.array(positions)


def get_interested_channels(bam_path, length_low, length_high):
    names = set()
    tot = 0
    filtered = 0
    with pysam.AlignmentFile(filename=bam_path, check_sq=False, mode="rb", threads=mp.cpu_count()) as bam:
        for record in tqdm(bam.fetch(until_eof=True), desc=f"reading {bam_path}"):
            tot += 1
            if record.query_length >= length_low and record.query_length < length_high:
                names.add(record.query_name)
                filtered += 1
    print(f"TOT:{tot}, filtered:{filtered}")

    return names


def main():

    parser = argparse.ArgumentParser("")
    parser.add_argument("--smc-bam", required=True, dest="smc_bam")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--length-low", default=1300, dest="length_low")
    parser.add_argument("--length-high", default=1700, dest="length_high")
    args = parser.parse_args()
    
    smc_bam = args.smc_bam
    vcf_file = args.vcf
    names = get_interested_channels(smc_bam, length_low=args.length_low, length_high=args.length_high)
    positions = read_pos_from_vcf_file(vcf_file, names)

    # ================================
    #  绘制位点分布直方图
    # ================================
    plt.figure(figsize=(12, 5))
    plt.hist(positions, bins=2000)  # bins 可根据点数调整
    plt.xlabel("Genomic Position (POS)")
    plt.ylabel("Variant Count")
    plt.title("Variant Site Position Distribution")
    plt.tight_layout()
    plt.show()
    plt.savefig("vcf_position.jpg")


if __name__ == "__main__":
    main()
    pass
