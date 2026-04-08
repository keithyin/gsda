#!/usr/bin/env python3
"""
minimap2_pipeline.py
使用 minimap2 进行比对，通过管道传递给 samtools sort 生成 BAM 文件，并建立索引。
支持直接输入 FASTQ/FASTA 或 BAM 文件作为 query。
"""

import argparse
import subprocess
import sys
import logging
from shutil import which
import multiprocessing as mp

logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def check_tool(tool_name):
    """检查工具是否在 PATH 中可用"""
    if which(tool_name) is None:
        logger.error(f"未找到 {tool_name}，请确保已安装并加入 PATH")
        sys.exit(1)


def is_bam_file(filepath):
    """通过文件扩展名判断是否为 BAM 文件"""
    return filepath.lower().endswith('.bam')


def run_pipeline(reference, query, output_bam, preset, minimap2_threads, sort_threads, fastq_threads, extra_minimap2_opts):
    # 检查必要工具
    check_tool("minimap2")
    check_tool("samtools")

    # 判断是否需要先将 BAM 转 FASTQ
    if is_bam_file(query):
        logger.info(f"检测到 BAM 输入文件: {query}，将使用 samtools fastq 转换为 FASTQ 流")
        check_tool("samtools")  # 已经检查过，但确保
        # 构建 samtools fastq 命令，输出到 stdout
        fastq_cmd = ["samtools", "fastq", "-@", str(fastq_threads), query]
        # minimap2 从 stdin 读取 query
        minimap2_input = "-"
        # 需要启动额外的 fastq 进程
        use_bam = True
    else:
        fastq_cmd = None
        minimap2_input = query
        use_bam = False

    # 构建 minimap2 命令
    minimap2_cmd = [
        "minimap2",
        "-ax", preset,
        "-t", str(minimap2_threads),
        reference, minimap2_input
    ]
    if extra_minimap2_opts:
        minimap2_cmd.extend(extra_minimap2_opts)

    # 构建 samtools sort 命令
    sort_cmd = [
        "samtools", "sort",
        "-@", str(sort_threads),
        "-o", output_bam,
        "-"
    ]

    logger.info(f"minimap2 命令: {' '.join(minimap2_cmd)}")
    logger.info(f"samtools sort 命令: {' '.join(sort_cmd)}")
    if use_bam:
        logger.info(f"samtools fastq 命令: {' '.join(fastq_cmd)}")

    # 构建进程管道
    if use_bam:
        # 管道: samtools fastq | minimap2 | samtools sort
        proc_fastq = subprocess.Popen(
            fastq_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        proc_minimap2 = subprocess.Popen(
            minimap2_cmd, stdin=proc_fastq.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        proc_fastq.stdout.close()  # 允许 fastq 进程在 sort 结束后收到 SIGPIPE
        proc_sort = subprocess.Popen(sort_cmd, stdin=proc_minimap2.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        proc_minimap2.stdout.close()
        # 等待所有进程
        fastq_ret = proc_fastq.wait()
        minimap2_ret = proc_minimap2.wait()
        sort_ret = proc_sort.wait()
        # 收集错误输出（可选）
        _, fastq_err = proc_fastq.communicate()
        _, minimap2_err = proc_minimap2.communicate()
        _, sort_err = proc_sort.communicate()
    else:
        # 管道: minimap2 | samtools sort
        proc_minimap2 = subprocess.Popen(
            minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        proc_sort = subprocess.Popen(sort_cmd, stdin=proc_minimap2.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        proc_minimap2.stdout.close()
        minimap2_ret = proc_minimap2.wait()
        sort_ret = proc_sort.wait()
        _, minimap2_err = proc_minimap2.communicate()
        _, sort_err = proc_sort.communicate()
        fastq_ret = 0
        fastq_err = ""

    # 错误处理
    if use_bam and fastq_ret != 0:
        logger.error(f"samtools fastq 失败，返回码 {fastq_ret}")
        if fastq_err:
            logger.error(f"samtools fastq 错误输出:\n{fastq_err}")
        sys.exit(1)

    if minimap2_ret != 0:
        logger.error(f"minimap2 失败，返回码 {minimap2_ret}")
        if minimap2_err:
            logger.error(f"minimap2 错误输出:\n{minimap2_err}")
        sys.exit(1)

    if sort_ret != 0:
        logger.error(f"samtools sort 失败，返回码 {sort_ret}")
        if sort_err:
            logger.error(f"samtools sort 错误输出:\n{sort_err}")
        sys.exit(1)

    logger.info(f"比对和排序完成，BAM 文件: {output_bam}")

    # 构建索引
    logger.info(f"为 {output_bam} 建立索引")
    index_cmd = ["samtools", "index", output_bam]
    try:
        subprocess.run(index_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"samtools index 失败: {e.stderr}")
        sys.exit(1)

    logger.info("索引建立完成")
    logger.info(f"最终输出: {output_bam} 和 {output_bam}.bai")


def main():
    parser = argparse.ArgumentParser(
        description="使用 minimap2 比对并生成排序后的 BAM 及索引，支持 BAM 格式的 query 自动转换")
    parser.add_argument("-r", "--reference", required=True,
                        help="参考基因组文件 (FASTA)")
    parser.add_argument("-q", "--query", required=True,
                        help="查询序列文件 (FASTA/FASTQ 或 BAM)")
    parser.add_argument("-o", "--output", required=True,
                        help="输出 BAM 文件路径 (例如: alignments.bam)")
    parser.add_argument("-p", "--preset", default="map-ont",
                        choices=["map-ont", "map-pb", "map-hifi", "sr"],
                        help="minimap2 预设 (默认: map-ont)")
    parser.add_argument("-t", "--minimap2-threads", type=int, default=mp.cpu_count(),
                        help="minimap2 线程数 (默认: 4)")
    parser.add_argument("-s", "--sort-threads", type=int, default=2,
                        help="samtools sort 线程数 (默认: 2)")
    parser.add_argument("-f", "--fastq-threads", type=int, default=1,
                        help="当 query 为 BAM 时，samtools fastq 使用的线程数 (默认: 1)")
    parser.add_argument("--extra", nargs=argparse.REMAINDER,
                        help="传递给 minimap2 的额外参数 (例如 --MD --secondary=no)")

    args = parser.parse_args()

    run_pipeline(
        reference=args.reference,
        query=args.query,
        output_bam=args.output,
        preset=args.preset,
        minimap2_threads=args.minimap2_threads,
        sort_threads=args.sort_threads,
        fastq_threads=args.fastq_threads,
        extra_minimap2_opts=args.extra
    )


if __name__ == "__main__":
    main()
