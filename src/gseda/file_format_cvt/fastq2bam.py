import pysam
import os
import math
import multiprocessing
from tqdm import tqdm
import argparse


def calc_rq(qualities):
    """
    根据 Phred 质量分数计算 reads 平均正确率 (RQ)。

    Phred score Q = -10 * log10(P_error)，其中 P_error 是碱基调用错误概率。
    单碱基正确率 = 1 - 10^(-Q/10)
    RQ = 所有碱基正确率的均值（整个 read 的平均准确度）

    Args:
        qualities (list[int]): Phred 质量分数列表（整数）

    Returns:
        float: RQ 值 (0.0 ~ 1.0)，例如 0.99 = 99% 平均正确率
    """
    if not qualities:
        return 0.0
    rq = sum(1 - math.pow(10, -q / 10) for q in qualities) / len(qualities)
    return round(rq, 6)


def extract_ch(qname):
    """
    从 qname 路径中抽取 channel 编号。

    PacBio qname 格式: run_id/ch_num/seq_num
    例如: 20260529_250302Y0004_Run0001/1163584/2 -> ch=1163584

    Args:
        qname (str): read name，可能包含 '/' 分隔的 channel 路径

    Returns:
        str or None: channel 编号（字符串），未找到返回 None
    """
    parts = qname.split("/")
    if len(parts) >= 2:
        # 取第二个 /-separated segment 作为 channel
        return parts[1]
    return None


def fastx_to_bam(input_fastx_path, output_bam_path, ref_fasta=None, read_group=None, threads=0):
    """
    将 FASTQ 或 FASTA 文件转换为 BAM 格式，所有 reads 以未映射状态写入。

    - FASTQ: 自动计算 RQ（Phred 质量平均正确率），提取 channel 编号
    - FASTA: 无质量分数，不计算 RQ

    Args:
        input_fastx_path (str): 输入 FASTQ (.fq/.fastq) 或 FASTA (.fa/.fasta) 文件路径。
        output_bam_path (str): 输出 BAM 文件路径。
        ref_fasta (str, optional): 参考基因组 FASTA 路径，用于生成 BAM header 中的 @SQ 行。
        read_group (str, optional): Read Group ID（写入 RG tag）。
        threads (int): 线程数 (默认 0 = 全部 CPU)。
    """
    # 检查输入文件是否存在
    if not os.path.exists(input_fastx_path):
        print(f"错误：输入文件不存在于 '{input_fastx_path}'")
        return

    # 根据扩展名判断输入格式（FASTQ 有质量分数，FASTA 无）
    input_extension = os.path.splitext(input_fastx_path)[1].lower()
    is_fastq = input_extension in ['.fq', '.fastq']

    # 计算线程数
    num_threads = threads if threads > 0 else multiprocessing.cpu_count()

    try:
        # 构建 BAM header
        header_dict = {'HD': {'VN': '1.5', 'SO': 'unknown'}}
        if read_group is not None:
            header_dict['RG'] = [{'ID': read_group}]
        if ref_fasta is not None:
            if not os.path.exists(ref_fasta):
                print(f"错误：参考基因组文件不存在于 '{ref_fasta}'")
                return
            # 通过 pysam.FastaFile 读取参考序列名称和长度，构建 @SQ 行
            with pysam.FastaFile(ref_fasta) as fa:
                ref_names = fa.references
                ref_lengths = fa.lengths
                header_dict['SQ'] = [
                    {'SN': name, 'LN': length}
                    for name, length in zip(ref_names, ref_lengths)
                ]
        header = pysam.AlignmentHeader.from_dict(header_dict)

        # 打开输入文件（pysam.FastxFile 自动检测 FASTQ/FASTA 格式）
        read_count = 0

        with pysam.FastxFile(input_fastx_path) as fx:
            # 打开 BAM 文件进行写入
            # 'wb' 表示二进制写入模式
            with pysam.AlignmentFile(
                output_bam_path, "wb", header=header,
                threads=num_threads, check_sq=False,
            ) as bam_out:

                for entry in tqdm(fx, desc=f"writing {output_bam_path}"):
                    # 使用 pysam.AlignedSegment 构造 read
                    read = pysam.AlignedSegment()
                    read.query_name = entry.name
                    read.query_sequence = entry.sequence if entry.sequence else ""
                    read.is_unmapped = True

                    # 从 qname 中抽取 channel 编号
                    ch = extract_ch(entry.name)
                    if ch is not None:
                        try:
                            read.set_tag('ch', int(ch), value_type='i')
                        except ValueError:
                            pass  # 非数字 channel，跳过

                    # FASTQ 包含质量分数，FASTA 不包含
                    if is_fastq and entry.quality is not None:
                        # entry.quality 是 Phred+33 编码的字符串，需要转换为整数列表
                        read.query_qualities = [ord(c) - 33 for c in entry.quality]

                        # 计算 RQ（reads 平均正确率）并写入标签
                        rq_value = calc_rq(read.query_qualities)
                        read.set_tag('rq', rq_value, value_type='d')

                    # 可选写入 RG tag
                    if read_group is not None:
                        read.set_tag('RG', read_group)

                    bam_out.write(read)
                    read_count += 1

        print(f"转换完成：'{input_fastx_path}' -> '{output_bam_path}'")
        print(f"成功写入了 {read_count} 条读取。")

    except Exception as e:
        print(f"处理文件时发生错误：{e}")


def main_cli():
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("inp", help="输入 FASTQ/FASTA 文件")
    parser.add_argument("oup", help="输出 BAM 文件")
    parser.add_argument("--ref", default=None, help="参考基因组 FASTA（可选，用于生成 @SQ header）")
    parser.add_argument("--rg-id", default=None, dest="rg_id", help="Read Group ID")
    parser.add_argument("--threads", default=0, type=int, help="线程数 (0 = 全部 CPU)")

    args = parser.parse_args()
    fastx_to_bam(args.inp, args.oup, ref_fasta=args.ref, read_group=args.rg_id, threads=args.threads)


# --- 示例用法 ---
if __name__ == "__main__":
    main_cli()
