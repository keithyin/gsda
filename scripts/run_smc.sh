#!/bin/bash

# 设置工作目录
WORK_DIR="/data1/ccs_data/20251216-hla-case-study/zj_HLA"

# 切换到工作目录
cd "$WORK_DIR"

# 查找所有以 called.bam 结尾的文件
for input_file in *_adapter.bam; do
    if [[ -f "$input_file" ]]; then
        # 提取基础文件名（去掉 _called.bam 后缀）
        base_name="${input_file%_adapter.bam}"
        
        # 设置输出文件名
        output_file="${base_name}.hmm.byStrand"
        
        echo "处理文件: $input_file -> $output_file"
        
        # 执行 docker 命令
        docker run -it -v /data1:/data1 --gpus all \
            192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.4_barcode_remover1.0.3_smicing0.4.4_bmi_0.1.4 \
             smc \
                ${WORK_DIR}/${input_file} \
                ${WORK_DIR}/${output_file} \
                --maxLength 10000000 \
                --minLength 10 \
                --hmmModelSpec Kit__500_Chem__1_BC__1_PW3_v4 \
                --byStrand 1 
        
        # 检查命令是否执行成功
        if [ $? -eq 0 ]; then
            echo "✓ 成功处理: $input_file"
        else
            echo "✗ 处理失败: $input_file"
        fi
        echo "----------------------------------------"
    fi
done

echo "所有文件处理完成！"
