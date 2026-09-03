"""
Run minimap2 alignment (fastq query -> fasta target), then sort and index the BAM.

Supports FASTQ or BAM as query input. When BAM is given, it is converted to
FASTQ on the fly via ``samtools fastq`` and the temporary file is cleaned up
on exit.

Pipeline:
    minimap2  ->  unsorted BAM
    samtools sort  ->  sorted BAM
    samtools index  ->  .bai
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile
import atexit
logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)

# --- temp file management ---
_temp_files: list[str] = []


def _clean_temps() -> None:
    for path in _temp_files:
        try:
            os.remove(path)
            log.info("Cleaned temp file: %s", path)
        except OSError:
            pass


atexit.register(_clean_temps)


def _try_bam2fq(bam_path: str) -> str:
    """Convert a BAM file to a temporary FASTQ file and return its path."""
    tmpfd, tmppath = tempfile.mkstemp(suffix=".fq", prefix="bam2fq_")
    os.close(tmpfd)
    _temp_files.append(tmppath)
    log.info("Converting BAM -> FASTQ  %s -> %s", bam_path, tmppath)
    cmd = ["samtools", "fastq", "-F", "2304", "-0", tmppath, bam_path]
    print(" ".join(cmd))
    # -F 2304 (0x900): skip unmapped + supplementary reads
    result = subprocess.run(
        cmd,
        capture_output=True,
    )

    if result.returncode != 0:
        log.error("samtools fastq failed (rc=%d)\n%s",
                  result.returncode, result.stderr.decode(errors="replace"))
        sys.exit(result.returncode)
    return tmppath


def _detect_format(path: str, explicit: str | None) -> str:
    """Return 'fq' or 'bam' based on explicit arg or file extension."""
    if explicit and explicit != "auto":
        return explicit
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    if ext == "bam":
        return "bam"
    return "fq"


def main():
    parser = argparse.ArgumentParser(
        prog="minimap2_align",
        description="Minimap2 alignment (fastq -> fasta), sort, and index BAM.",
    )
    parser.add_argument("query", help="Input query file (FASTQ or BAM)")
    parser.add_argument("target", help="Target reference FASTA file")
    parser.add_argument("-o", "--out", default="aligned",
                        help="Output BAM base name (default: aligned)")
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count()
                        or 1, help="Number of threads (default: auto)")
    parser.add_argument("--preset", default="map-ont", choices=["splice", "map-pb", "map-ont", "sr", "hisat2"],
                        help="minimap2 preset (default: map-ont)")
    parser.add_argument("--format", choices=["auto", "fq", "bam"], default="auto",
                        help="Query format (default: auto, inferred from extension)")
    args = parser.parse_args()

    query_format = _detect_format(args.query, args.format)

    # --- Handle BAM query: convert to FASTQ ---
    query_file = args.query
    if query_format == "bam":
        query_file = _try_bam2fq(args.query)

    # --- Step 1: minimap2 alignment ---
    log.info("Step 1/3: minimap2 alignment  %s -> %s", query_file, args.target)
    minimap_cmd = [
        "minimap2",
        "-ax", args.preset,
        "-t", str(args.threads),
        "--eqx",
        "--secondary=no",
        "-a",          # SAM output (text)
        args.target,
        query_file,
    ]
    minimap_proc = subprocess.Popen(
        minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    sort_cmd = [
        "samtools", "sort", "-@", str(args.threads),
        "-o", f"{args.out}.sorted.bam",
        "-",
    ]
    sort_proc = subprocess.Popen(
        sort_cmd, stdin=minimap_proc.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
    )
    minimap_proc.stdout.close()
    _, minimap_err = minimap_proc.communicate()
    _, sort_err = sort_proc.communicate()

    if minimap_proc.returncode != 0:
        log.error("minimap2 failed (rc=%d):\n%s",
                  minimap_proc.returncode, minimap_err.decode(errors="replace"))
        sys.exit(minimap_proc.returncode)

    if sort_proc.returncode != 0:
        log.error("samtools sort failed (rc=%d):\n%s",
                  sort_proc.returncode, sort_err.decode(errors="replace"))
        sys.exit(sort_proc.returncode)

    if sort_err:
        log.info("samtools sort: %s", sort_err.decode(
            errors="replace").strip())

    # --- Step 3: index ---
    sorted_bam = f"{args.out}.sorted.bam"
    bai_path = sorted_bam + ".bai"
    if os.path.exists(bai_path):
        log.info("Removing existing index: %s", bai_path)
        os.remove(bai_path)
    log.info("Step 3/3: samtools index  %s", sorted_bam)
    idx_result = subprocess.run(
        ["samtools", "index", sorted_bam],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if idx_result.returncode != 0:
        log.error("samtools index failed (rc=%d):\n%s",
                  idx_result.returncode, idx_result.stderr.decode(errors="replace"))
        sys.exit(idx_result.returncode)

    log.info("Done. Output: %s  %s.bai", sorted_bam, bai_path)


if __name__ == "__main__":
    main()
