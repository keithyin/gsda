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
        output_file="${base_name}"
        
        echo "处理文件: $input_file -> $output_file"
        
        # 执行 docker 命令
        docker run -it -v /data1:/data1 --gpus all \
            192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.4_barcode_remover1.0.3_smicing0.4.4_bmi_0.1.4 \
            smicing consensus \
                ${WORK_DIR}/${input_file} \
                ${WORK_DIR}/${output_file} \
                --maxLength 10000000 \
                --minLength 10 \
                --modelSpecs /root/models/2025Q1-stage1-polish-cnn-dw-smicing-epo50.tar.gz \
                --modelSpecs /root/models/hg-2025Q1-icing-selfattn-dw-qhead-ls0.01-smicing-epo32-fp16-2o.tar.gz \
                --modelSpecs /root/models/hg-2025Q1-icing-selfattn-dw-qhead-ls0.01-smicing-epo32-fp16-2o.tar.gz \
                --enableCalibration 1 \
                --poaEdUnifyStrand 1 
        
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


#!/bin/bash

# 获取所有以subreads.bam结尾的文件列表
files=(*adapter.bam)

# 检查是否找到文件
if [ ${#files[@]} -eq 0 ]; then
    echo "未找到任何以subreads.bam结尾的文件"
    exit 1
fi

echo "找到以下待处理文件:"
printf '%s\n' "${files[@]}"

# 遍历处理每个文件
for file in "${files[@]}"; do
    echo "正在处理文件: $file"
    fstem="${file%.bam}"
    echo "$file -> $fstem"

    sync
    sudo echo 3 > /proc/sys/vm/drop_caches

    docker run --rm -t --gpus all \
        -v `pwd`:/data \
        192.168.3.38:5000/algo/adacus:smc5.5.0_adapter_demux0.0.5_barcode_remover1.0.3_smicing0.6.1_bmi_0.1.5_gseda_1.19.2 \
        /bin/bash -c " \
            LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2  smicing consensus /data/$file /data/$fstem.new \
                --smcMinPasses 1 \
                --smcMaxPasses 25 \
                --maxLength 10000000 --minLength 10 \
                --modelSpecs /root/models/2025Q1-stage1-polish-cnn-dw-smicing-epo50.tar.gz \
                --modelSpecs /root/models/hg-2025Q1-icing-selfattn-dw-qhead-ls0.01-smicing-epo32-fp16-2o.tar.gz \
                --poaEdUnifyStrand 1 \
        "
done

echo "所有文件处理完成"
