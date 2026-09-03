#!/usr/bin/env bash
# SNP / Insertion / Deletion 位点 Q 值分布分析管线驱动脚本。
# 对 5 个 barcode 依次执行: 合并 -> Q30 过滤 -> 变体 calling (SNP/INS/DEL)
# -> 各类型 Q 值分布统计 + VAF>10% 位点 Q 基本统计。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"          # snp_q_analysis/
RESULTS_DIR="$ROOT_DIR/results"

DATA_DIR="/data1/ccs_data/20260805-xunming"
RUN1_DIR="$DATA_DIR/run1"
RUN2_DIR="$DATA_DIR/run2"

BARCODES="Barcode01 Barcode02 Barcode03 Barcode04 Barcode05"
VAR_TYPES="SNP INS DEL"
MIN_Q=30
MIN_IDENTITY=0.9
VAF_THRESHOLD=0.10
THREADS=8

mkdir -p "$RESULTS_DIR"/{merged,filtered,variant,qdist,vaf_q_stats}

for barcode in $BARCODES; do
    echo ""
    echo "==================== $barcode ===================="

    # 步骤 1: 合并 run1/run2 相同 barcode
    python3 "$SCRIPT_DIR/01_merge_fastq.py" "$barcode" \
        --run1-dir "$RUN1_DIR" \
        --run2-dir "$RUN2_DIR" \
        -o "$RESULTS_DIR/merged/$barcode.fastq"

    # 步骤 2a: Q30 过滤 (readsQ>=30)
    python3 "$SCRIPT_DIR/02_q30_filter.py" "$RESULTS_DIR/merged/$barcode.fastq" \
        -o "$RESULTS_DIR/filtered/$barcode.q30.fastq" \
        --min-q "$MIN_Q"

    # 步骤 2b: 比对参考 -> 变体位点 (SNP/INS/DEL)
    python3 "$SCRIPT_DIR/03_variant_calling.py" \
        --fastq "$RESULTS_DIR/filtered/$barcode.q30.fastq" \
        --ref "$DATA_DIR/${barcode}_300_draft.fasta" \
        --barcode "$barcode" \
        -o "$RESULTS_DIR/variant/$barcode.variants.tsv" \
        --min-identity "$MIN_IDENTITY" \
        -t "$THREADS"
done

echo ""
echo "==================== 各类型 Q 值分布 ===================="
for barcode in $BARCODES; do
    for vtype in $VAR_TYPES; do
        python3 "$SCRIPT_DIR/04_snp_q_dist.py" "$RESULTS_DIR/variant/$barcode.variants.tsv" \
            --type "$vtype" \
            --out-prefix "$RESULTS_DIR/qdist/${barcode}_${vtype}_qdist"
    done
done

echo ""
echo "全部完成。结果目录: $RESULTS_DIR"
echo "说明: VAF>10% 位点 Q 统计 (05) 需按各 barcode 的 full-length read 数传 --n-reads，例如:"
echo "  python3 $SCRIPT_DIR/05_snp_vaf_q_stats.py $RESULTS_DIR/variant/Barcode01.variants.tsv \\"
echo "      --type INS --n-reads 94516 --out-prefix $RESULTS_DIR/vaf_q_stats/Barcode01"
