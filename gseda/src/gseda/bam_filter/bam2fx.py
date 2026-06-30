import argparse
import multiprocessing as mp
from pathlib import Path
import math
import pysam
from tqdm import tqdm


def bam_to_fastx(bam_path: str, out_path: str, fmt: str, min_rq: float):
    """Dump reads with rq >= min_rq to a single FASTA/FASTQ file."""
    tot = 0
    dumped = 0

    ext = "fastq" if fmt == "fq" else "fasta"
    with pysam.AlignmentFile(
        bam_path, "rb", threads=mp.cpu_count(), check_sq=False
    ) as bam_file, open(out_path, "w") as fout:
        for read in tqdm(
            bam_file.fetch(until_eof=True), desc=f"dumping {bam_path} to {out_path}"
        ):
            tot += 1
            try:
                rq = read.get_tag("rq")
                if rq < min_rq:
                    continue
            except KeyError:
                pass

            name = read.query_name
            seq = read.query_sequence
            qual = read.qual
            if seq is None or qual is None:
                continue

            dumped += 1
            if ext == "fastq":
                fout.write(f"@{name}\n{seq}\n+\n{qual}\n")
            else:
                fout.write(f">{name}\n{seq}\n")

    print(f"Tot:{tot}, dumped:{dumped}, ratio:{dumped / tot: .4f}")
    print(f"转换完成，输出文件: {out_path}")


def main_cli():
    parser = argparse.ArgumentParser(
        description="Convert BAM to FASTA/FASTQ files grouped by rq threshold.")
    parser.add_argument("bam", help="Input BAM file path")
    parser.add_argument(
        "fastx", choices=["fa", "fq"],
        help="Output format: 'fa' for FASTA, 'fq' for FASTQ",
    )
    parser.add_argument(
        "--phreqs", type=int, default=[20, 30], nargs="+",
        help="Phred-scaled rq thresholds (default: [20, 30]). "
             "Produces one output file per threshold.",
    )

    args = parser.parse_args()
    bam_path_abs = str(Path(args.bam).resolve())
    parent = Path(bam_path_abs).parent
    stem = Path(bam_path_abs).stem
    for thresh in args.phreqs:
        out_path = str(parent / f"{stem}.q{thresh}.{args.fastx}")
        bam_to_fastx(args.bam, out_path, args.fastx,
                     1.0 - math.pow(10.0, -thresh / 10.0))


if __name__ == "__main__":

    main_cli()
