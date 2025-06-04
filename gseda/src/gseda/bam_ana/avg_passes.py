
import pysam
from tqdm import tqdm
def compute_avg_np_for_high_rq(bam_path, rq_threshold=0.99):
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=40)
    np_values = []
    for read in tqdm(bam, desc=f"reading {bam_path}"):
        if read.has_tag("rq") and read.get_tag("rq") >= rq_threshold:
            if read.has_tag("np"):
                np_values.append(read.get_tag("np"))
    bam.close()
    
    if not np_values:
        print("No records found with rq >= threshold and np tag.")
        return None
    
    avg_np = sum(np_values) / len(np_values)
    return avg_np

def main():
    # 使用示例
    bam_file = "/data/ccs_data/little-mouse/mouse-smc.q20.bam"
    average_np = compute_avg_np_for_high_rq(bam_file)
    if average_np is not None:
        print(f"Average np for rq >= 0.99: {average_np:.2f}")

if __name__ == "__main__":
    main()