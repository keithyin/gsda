import pysam
from tqdm import tqdm
import numpy as np

def q30_base_ratio_in_q30_channel(bam_file: str):
    q30_base_cnt = 0
    tot_base_cnt = 0
    with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False) as bam_in:
        for record in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {bam_file}"):
            rq = float(record.get_tag("rq"))
            if rq < 0.999:
                continue
            qual = record.query_qualities
            qual = np.array(qual)
            q30_base_cnt += (qual >= 30).sum()
            tot_base_cnt += len(qual)

    print(f"q30_base_ratio_in_q30_channel:{q30_base_cnt / tot_base_cnt * 100. :3f}%")

def main():
    fname = "/data/ccs_data/20250804-ludaopei/20250721_240601Y0005_Run0001.polish2.smc_all_reads.bam"
    fname = "/data/ccs_data/20250804-ludaopei/20250721_240601Y0005_Run0001.smc_all_reads.bam"
    q30_base_ratio_in_q30_channel(fname)


if __name__ == "__main__":
    main()
    pass
             