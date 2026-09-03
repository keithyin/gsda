# 基于混样的结果进行进行错误分析
import pysam
import multiprocessing as mp
import pathlib
from tqdm import tqdm
import subprocess

def extract_sample_channels(bam_file: str, cq_thr=0):
    channels = set()
    with pysam.AlignmentFile(bam_file, mode="rb", threads=mp.cpu_count(), check_sq=False) as bam_in:
        for record in bam_in.fetch(until_eof=True):
            query_name = record.query_name
            if record.has_tag("cq"):
                cq = record.get_tag("cq")
                if cq < cq_thr:
                    continue
            assert isinstance(query_name, str)
            channel_id = int(query_name.split("/")[1].split("_")[0])
            channels.add(channel_id)
    return channels

def filter_bam(bam_path: str, channel_ids, infix="subset"):
    p = pathlib.Path(bam_path)
    stem = p.stem
    parent = p.parent
    tmp_bam_path = parent.joinpath(f"{stem}.{infix}.bam")
    
    tot = 0
    output_cnt = 0

    with pysam.AlignmentFile(bam_path, mode="rb", check_sq=False, threads=mp.cpu_count() // 2) as bam_file:
        with pysam.AlignmentFile(str(tmp_bam_path), mode="wb", check_sq=False, threads=mp.cpu_count() // 2, header=bam_file.header) as bam_out:
            for record in tqdm(bam_file.fetch(until_eof=True), desc=f"{bam_path} -> {str(tmp_bam_path)}"):
                ch = int(record.get_tag("ch"))
                tot += 1
                if ch in channel_ids:
                    bam_out.write(record)
                    output_cnt += 1
    print(f"tot:{tot}, output_cnt:{output_cnt}. output_ratio:{output_cnt / tot * 100.:.2f}%")
    return str(tmp_bam_path)


def main():
    # barcode7 -> barcode227-0
    barcode = "227-0"
    # smc2consensus_alignment_bam = f"/data1/ccs_data/20251229-saisuofei/20251229_241201Y0002_Run0001/Group_0/barcodes_reads_cons_gen_amplicon/Consensus/Bam/Group_0_Adaptor-barcode{barcode}.sort.bam"
    bam_file = f"/data1/ccs_data/20260108-saisuofei-fuce/20260107_240601Y0088_Run0005/barcode_remover/Barcode7.bam"
    sbr_bam = "/data1/ccs_data/20260108-saisuofei-fuce/20260107_240601Y0088_Run0005/20260107_240601Y0088_Run0005_adapter.bam"
    smc_bam = "/data1/ccs_data/20260108-saisuofei-fuce/20260107_240601Y0088_Run0005/20260107_240601Y0088_Run0005.smc_all_reads.bam"
    ref_fa = f"/data1/ccs_data/20251229-saisuofei/barcode{barcode}.ref.fasta"
    
    channel_ids = extract_sample_channels(bam_file, cq_thr=20.)
    sbr_subset_bam = filter_bam(sbr_bam, channel_ids)
    smc_subset_bam = filter_bam(smc_bam, channel_ids)
    prefix = f"saisuofei-barcode-{barcode}-sbr2smc"
    cmd = f"asrtc --ref-fa {ref_fa} -q {sbr_subset_bam} -t {smc_subset_bam} -p {prefix}-new-param --ref-range 80:100 --rq-range 0.99:1.1 --np-range 7:10000 -m 4 -M 10 -o 4,48 -e 3,1"
    print(f"running cmd:{cmd}")
    subprocess.check_call(cmd, shell=True)
    

if __name__ == "__main__":
    main()
    pass