import pysam
import multiprocessing as mp
import pathlib
from typing import Set
from tqdm import tqdm
import subprocess


def extract_query_name_from_bam(bam_path):
    smc_names = set()
    with pysam.AlignmentFile(filename=bam_path, mode="rb", threads=mp.cpu_count()) as bam_in:
        for record in bam_in.fetch(until_eof=True):
            smc_names.add(record.query_name)
    return smc_names


def dump_subset_smc_bam(smc_bam_path: str, smc_names: Set[str]) -> str:
    p = pathlib.Path(smc_bam_path)

    subset_bam_path = p.parent.joinpath(f"{p.stem}.subset.bam")
    with pysam.AlignmentFile(filename=smc_bam_path, mode="rb", threads=mp.cpu_count() // 2, check_sq=False) as bam_in:
        with pysam.AlignmentFile(filename=str(subset_bam_path), mode="wb", threads=mp.cpu_count() // 2, check_sq=False, header=bam_in.header) as bam_out:
            for record in tqdm(bam_in.fetch(until_eof=True), desc=f"dumping {smc_bam_path}"):
                if record.query_name in smc_names:
                    bam_out.write(record)
    return str(subset_bam_path)


def run_asrtc(ref_fa, subreads_bam, smc_bam, prefix):
    cmd = f"asrtc --ref-fa {ref_fa} -q {subreads_bam} -t {smc_bam} --np-range 0:100 --rq-range 0.99:1.1 -p {prefix} -m 2 -M 4 -o 4,24 -e 2,1"
    print(f"cmd:{cmd}")
    subprocess.check_call(cmd, shell=True)


def main():

    bam_path = "/data1/ccs_data/saisuofei/20251204_241201Y0002_Run0001/Group_0/barcodes_reads_cons_gen_amplicon/Consensus/Bam/Group_0_Adaptor-barcode203-0.sort.bam"
    smc_bam_path = "/data1/ccs_data/saisuofei/20251204_241201Y0002_Run0001.smc_all_reads.bam"
    subreads_bam_path = "/data1/ccs_data/saisuofei/20251204_241201Y0002_Run0001_adapter.bam"
    ref_fasta = "/data1/ccs_data/saisuofei/Group_0_Adaptor-barcode203-0/Group_0_Adaptor-barcode203-0/ref.fasta"
    smc_names = extract_query_name_from_bam(bam_path=bam_path)
    print(f"num smc reads:{len(smc_names)}")
    smc_subset_bam_path = dump_subset_smc_bam(smc_bam_path, smc_names)
    run_asrtc(ref_fa=ref_fasta, subreads_bam=subreads_bam_path,
              smc_bam=smc_subset_bam_path, prefix="saisuofei-barcode-203")

    pass


if __name__ == "__main__":
    main()
    pass
