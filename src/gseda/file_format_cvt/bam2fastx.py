import pysam
import os
import multiprocessing
from tqdm import tqdm
import argparse

def bam_to_fastx(input_bam_path, output_file_path, rq_threshold=None):
    """
    将 BAM 文件转换为 FASTA 或 FASTQ 格式，并可选地根据 RQ 标签进行过滤。

    Args:
        input_bam_path (str): 输入 BAM 文件的路径。
        output_file_path (str): 输出 FASTA (.fa, .fasta) 或 FASTQ (.fq, .fastq) 文件的路径。
        rq_threshold (float, optional): RQ 标签的最低阈值。
                                        只有 RQ 值大于或等于此阈值的读取才会被写入。
                                        如果为 None，则不进行 RQ 过滤。
    """
    # 检查输入文件是否存在
    if not os.path.exists(input_bam_path):
        print(f"错误：输入 BAM 文件不存在于 '{input_bam_path}'")
        return

    # 根据输出文件扩展名确定输出格式
    output_extension = os.path.splitext(output_file_path)[1].lower()
    if output_extension in ['.fa', '.fasta']:
        output_format = 'fasta'
    elif output_extension in ['.fq', '.fastq']:
        output_format = 'fastq'
    else:
        print(
            f"错误：不支持的输出文件扩展名 '{output_extension}'。请使用 .fa/.fasta 或 .fq/.fastq。")
        return

    try:
        # 打开 BAM 文件进行读取
        # 'rb' 表示以二进制读取模式打开 BAM 文件
        with pysam.AlignmentFile(input_bam_path, "rb", check_sq=False, threads=multiprocessing.cpu_count()) as bam_file:
            # 打开输出文件进行写入
            # 'w' 表示写入模式，会覆盖现有文件
            with open(output_file_path, "w") as out_file:
                read_count = 0
                filtered_count = 0

                for read in tqdm(bam_file, desc=f"reading {input_bam_path}"):

                    # 检查 RQ 标签并进行过滤
                    if rq_threshold is not None:
                        try:
                            # 尝试获取 RQ 标签的值
                            # RQ 标签通常存储为浮点数
                            rq_value = read.get_tag('rq')
                            if rq_value < rq_threshold:
                                filtered_count += 1
                                continue  # 跳过不符合阈值的读取
                        except KeyError:
                            # 如果读取没有 RQ 标签，可以选择跳过或包含它
                            # 这里选择跳过，你可以根据需要修改
                            # print(f"警告：读取 '{read.query_name}' 没有 'RQ' 标签，将被跳过。")
                            filtered_count += 1
                            continue

                    read_count += 1

                    if output_format == 'fasta':
                        # FASTA 格式：>read_name\nsequence
                        out_file.write(
                            f">{read.query_name}\n{read.query_sequence}\n")
                    elif output_format == 'fastq':
                        # FASTQ 格式：@read_name\nsequence\n+\nquality_scores (Phred+33)
                        # pysam.AlignedSegment.query_qualities 返回 Phred 质量分数（整数列表）
                        # 需要转换为 ASCII 字符表示 (Phred+33)
                        quality_string = "".join(
                            [chr(q + 33) for q in read.query_qualities]
                        )
                        out_file.write(
                            f"@{read.query_name}\n{read.query_sequence}\n+\n{quality_string}\n"
                        )

                print(f"转换完成：'{input_bam_path}' -> '{output_file_path}'")
                print(f"总共处理了 {read_count + filtered_count} 条读取。")
                print(
                    f"根据 RQ 阈值 ({rq_threshold if rq_threshold is not None else '无'}) 过滤了 {filtered_count} 条读取。")
                print(f"成功写入了 {read_count} 条读取。")

    except Exception as e:
        print(f"处理文件时发生错误：{e}")


def main_cli():
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("inp")
    parser.add_argument("oup")
    parser.add_argument("--rq-thr", default=None, type=float, dest="rq_thr")
    
    args = parser.parse_args()
    bam_to_fastx(args.inp, args.oup, rq_threshold=args.rq_thr)
    
    pass

# --- 示例用法 ---
if __name__ == "__main__":
    main_cli()
