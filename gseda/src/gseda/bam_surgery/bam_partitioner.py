
import pysam
from multiprocessing import cpu_count
import os
import pathlib
from tqdm import tqdm


def partition(bam_file: str, num_par: int):

    thread_per_bam = cpu_count() // num_par
    thread_per_bam = max(1, thread_per_bam)

    out_dir = os.path.dirname(bam_file)
    stem = pathlib.Path(bam_file).stem
    out_filenames = []
    for idx in range(0, num_par):
        out_filenames.append(os.path.join(out_dir, f"{stem}.{idx}.bam"))

    with pysam.AlignmentFile(bam_file, mode="rb", check_sq=False, threads=thread_per_bam) as bam_in:

        out_bams = [pysam.AlignmentFile(
            fname, mode="wb", check_sq=False, threads=thread_per_bam, header=bam_in.header) for fname in out_filenames]

        for (cnt, record) in tqdm(enumerate(bam_in.fetch(until_eof=True)), desc=f"dumping {bam_file}"):
            idx = cnt % num_par
            out_bams[idx].write(record)

    for out_bam in out_bams:
        out_bam.close()


def main():
    bam_file = "/data/512-data-for-metric/20250825_240601Y0012_Run0001.called.valid.bam"
    par = 4
    partition(bam_file=bam_file, num_par=par)
    


if __name__ == "__main__":
    main()
