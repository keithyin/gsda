#!/bin/bash

# Usage: ./process_all_runs.sh [run1 run2 ...]
# If no arguments, process all subdirectories that contain *_adapter.bam

BASE_DIR=`pwd`
echo "BASE_DIR:${BASE_DIR}"

# Determine which runs to process
if [ $# -gt 0 ]; then
    runs=("$@")   # use command line arguments
else
    # Find all directories directly under BASE_DIR that contain a *_adapter.bam file
    mapfile -t runs < <(find "$BASE_DIR" -maxdepth 2 -mindepth 1 -type f -name "*_adapter.bam" -exec dirname {} \; | sort -u)
    if [ ${#runs[@]} -eq 0 ]; then
        echo "Error: No run directories with *_adapter.bam found in $BASE_DIR" >&2
        exit 1
    fi
fi

echo "process runs:${runs[@]}"


# Create the conversion script once in BASE_DIR (reused by all runs)
cat > "$BASE_DIR/bam2fx.py" << 'EOF'
import pysam
import argparse
from tqdm import tqdm
import multiprocessing as mp
import pathlib

def bam_to_fastx(bam_path, fastx_path: str, rq_threshold):
    tot = 0
    dumped = 0
    output_fq = fastx_path.endswith("fastq")
    from_name = pathlib.Path(bam_path).name
    to_name = pathlib.Path(fastx_path).name
    with pysam.AlignmentFile(
        bam_path, "rb", threads=mp.cpu_count(), check_sq=False
    ) as bam_file, open(fastx_path, "w") as fastx_out:
        for read in tqdm(
            bam_file.fetch(until_eof=True), desc=f"dumping {from_name} to {to_name}"
        ):
            tot += 1

            try:
                rq = read.get_tag("rq")
                if rq < rq_threshold:
                    continue
            except KeyError:
                pass

            name = read.query_name
            seq = read.query_sequence
            qual = read.qual

            dumped += 1

            if seq is None or qual is None:
                continue
            if output_fq:
                fastx_out.write(f"@{name}\n{seq}\n+\n{qual}\n")
            else:
                fastx_out.write(f">{name}\n{seq}\n")

    print(f"Tot:{tot}, dumped:{dumped}, ratio:{dumped / tot: .4f}")
    print(f"转换完成，输出文件: {fastx_path}")


def main_cli():
    parser = argparse.ArgumentParser(description="Convert BAM to FASTQ with rq filter.")
    parser.add_argument("bam", help="Input BAM file path")
    parser.add_argument("fastx", help="Output FASTX file path. fasta/fastq")
    parser.add_argument(
        "--rq", type=float, default=0.0, help="Minimum rq threshold (default: 0.0)"
    )

    args = parser.parse_args()
    assert args.fastx.endswith("fastq") or args.fastx.endswith(
        "fasta"), "only fastq/fasta are supported"
    bam_to_fastx(args.bam, args.fastx, args.rq)


if __name__ == "__main__":
    main_cli()
EOF

# Process each run
for run_dir in "${runs[@]}"; do
    # Ensure run_dir is an absolute path
    if [[ "$run_dir" != /* ]]; then
        run_dir="$BASE_DIR/$run_dir"
    fi

    if [ ! -d "$run_dir" ]; then
        echo "Warning: Directory $run_dir does not exist, skipping." >&2
        continue
    fi

    echo "=========================================="
    echo "Processing run: $run_dir"
    echo "=========================================="

    cd "$run_dir" || { echo "Cannot cd into $run_dir, skipping"; continue; }

    # Find exactly one _adapter.bam
    adapter_files=( *_adapter.bam )
    if [ ${#adapter_files[@]} -eq 1 ]; then
        adapter_file="${adapter_files[0]}"
    else
        echo "Error: found ${#adapter_files[@]} adapter files in $run_dir, expected exactly one. Skipping." >&2
        continue
    fi

    output_prefix="${adapter_file%_adapter.bam}.bystrand"
    echo "Processing $adapter_file -> $output_prefix"

    IMAGE_BASE="192.168.3.38:5000/algo/adacus"
    # Get latest image
    LATEST_IMAGE=$(docker images --format "table {{.Repository}}:{{.Tag}}\t{{.CreatedAt}}" | \
                   grep "^${IMAGE_BASE}:" | \
                   sort -k2 -r | \
                   head -n1 | \
                   awk '{print $1}')

    # Run the smc command (mount current run directory to /data)
    docker run -i -v "$(pwd):/data" --gpus all \
        "$LATEST_IMAGE" \
        smc \
        "/data/${adapter_file}" \
        "/data/${output_prefix}" \
        --maxLength 10000000 \
        --minLength 10 \
        --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4 \
        --byStrand 1

    bystrand_filename="${output_prefix}.smc_all_reads.bam"

    # Run conversions using the shared bam2fx.py script (absolute path)
    PY_CMD="/home/user/anaconda3/envs/conda_env38/bin/python"
    for thresh in "0.84151" "0.9" "0.99" "0.999"; do
        case $thresh in
            0.84151) suffix="q8" ;;
            0.9)     suffix="q10" ;;
            0.99)    suffix="q20" ;;
            0.999)   suffix="q30" ;;
        esac
        docker exec -i gseq_product_dna \
            "$PY_CMD" "$BASE_DIR/bam2fx.py" \
            "$(pwd)/${bystrand_filename}" \
            "$(pwd)/${output_prefix}.${suffix}.fastq" \
            --rq "$thresh"
    done

    echo "Finished run: $run_dir"
done

echo "All runs processed."