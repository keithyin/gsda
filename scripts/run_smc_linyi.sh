#!/bin/bash

HOST_WORK_DIR="/data1/EurusResV3"
CONTAINER_WORK_DIR="/data"

# 支持包含子目录的文件列表
FILE_LIST=(
    "20260202_250301Y0003_Run0002/20260202_250301Y0003_Run0002_adapter.bam"
)

cd "$HOST_WORK_DIR"

for input_file in "${FILE_LIST[@]}"; do
    input_file=$(echo "$input_file" | xargs)
    [[ -z "$input_file" || ! -f "$input_file" ]] && echo "跳过: $input_file" && continue
    
    base_name="${input_file%.*}"
    output_file="${base_name}.byStrand"
    output_dir=$(dirname "$output_file")
    
    [[ "$output_dir" != "." ]] && mkdir -p "$output_dir"
    
    echo "处理: $input_file -> $output_file"
    
    docker run -it -v "$HOST_WORK_DIR:$CONTAINER_WORK_DIR" --gpus all \
        192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.5_barcode_remover1.0.3_smicing0.5.2_bmi_0.1.5 \
        smc \
            "$CONTAINER_WORK_DIR/$input_file" \
            "$CONTAINER_WORK_DIR/$output_file" \
            --maxLength 10000000 \
            --minLength 10 \
            --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4 \
            --byStrand 1 
    
    [ $? -eq 0 ] && echo "✓ 成功: $input_file" || echo "✗ 失败: $input_file"
    echo "----------------------------------------"
done

echo "完成！"