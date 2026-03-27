
#!/bin/bash

# 指定目标目录
TARGET_DIR="/data1/ccs_data/chicken-data"

# 检查目录是否存在
if [ ! -d "$TARGET_DIR" ]; then
    echo "错误: 目录 $TARGET_DIR 不存在。"
    exit 1
fi

# 进入目录
cd "$TARGET_DIR" || exit

# 遍历所有符合条件的文件
for bam_file in *new.smc_all_reads.bam; do
    # 检查是否有匹配的文件，防止通配符未匹配时报错
    [ -e "$bam_file" ] || continue

    # 生成对应的 fastq 文件名：将 .bam 替换为 .fastq
    # ${bam_file%.bam} 表示去掉文件名末尾的 .bam
    fastq_file="${bam_file%.bam}.fastq"

    echo "正在处理: $bam_file -> $fastq_file"

    # 执行转换命令
    bam2fx-cli "$bam_file" "$fastq_file"
done

echo "任务完成！"