import os
import pysam
import logging
import argparse
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def filter_called_bam(smc_bam, called_bam, output_bam):
    """
    对bam文件进行过滤
    :param smc_bam: smc之后的bam文件
    :param called_bam: basecaller之后的bam文件
    :param output_bam: 过滤之后的output bam文件
    :param temp_dir: 临时目录
    :return:
    """
    # 调用日志系统

    # 检查输入文件路径是否存在
    if not os.path.exists(smc_bam):
        logging.error(f"Input SMC BAM file {smc_bam} does not exist.")
        raise FileNotFoundError(
            f"Input SMC BAM file {smc_bam} does not exist.")  # 抛出异常

    if not os.path.exists(called_bam):
        logging.error(f"Input called BAM file {called_bam} does not exist.")
        raise FileNotFoundError(
            f"Input called BAM file {called_bam} does not exist.")  # 抛出异常

    smc_read_names = set()
    try:
        with pysam.AlignmentFile(smc_bam, "rb", check_sq=False, threads=os.cpu_count()) as a_bam:
            for read in tqdm(a_bam.fetch(until_eof=True), desc=f"reading {smc_bam}"):
                base_name = read.query_name.split()[0]
                # Adjust the parsing according to your read name format
                channel_number = base_name.split('/')[1]
                smc_read_names.add("read_" + channel_number)
        logging.info(
            f"Extracted {len(smc_read_names)} unique read names from {smc_bam}.")
    except Exception as e:
        logging.error(
            f"Failed to extract read names from {smc_bam}. Error: {e}")
        raise  # 重新抛出异常

    with pysam.AlignmentFile(filename=called_bam, mode="rb", check_sq=False, threads=os.cpu_count() // 2) as in_bam:
        with pysam.AlignmentFile(filename=output_bam, mode="wb", check_sq=False, threads=os.cpu_count() // 2, header=in_bam.header) as out_bam:
            cnt = 0
            for read in tqdm(in_bam.fetch(until_eof=True), desc=f"dumping {called_bam} to {output_bam}"):
                if read.query_name in smc_read_names:
                    cnt += 1 
                    out_bam.write(read)
            logging.info(f"dump {cnt} records")

def main(smc_bam: str, called_bam: str, output_bam: str):
    filter_called_bam(smc_bam=smc_bam, called_bam=called_bam,
                      output_bam=output_bam)



def main_cli():

    parser = argparse.ArgumentParser(prog="parser")
    parser.add_argument("--called-bam", type=str,
                        required=True, dest="called_bam")
    parser.add_argument("--smc-bam", type=str, required=True, dest="smc_bam")
    parser.add_argument("--out-bam", type=str, default=None, dest="out_bam")

    args = parser.parse_args()
    called_bam = args.called_bam
    smc_bam = args.smc_bam
    assert isinstance(called_bam, str)
    out_bam = args.out_bam
    if out_bam is None:
        out_bam = "{}.output.bam".format(called_bam.rsplit(".")[0])
    main(smc_bam, called_bam, out_bam)        

if __name__ == "__main__":
    main_cli()
