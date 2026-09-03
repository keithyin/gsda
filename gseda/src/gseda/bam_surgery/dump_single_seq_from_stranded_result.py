import pysam
from typing import Dict
from tqdm import tqdm

def extract_interested_channels(fwd_rev_input_bam: str):
    bam_in = pysam.AlignmentFile(fwd_rev_input_bam, "rb", check_sq=False, threads=40)
    # 存储每个 channel 的前向和反向 reads
    channel_reads: Dict[str, Dict[str, pysam.AlignedSegment]] = {}
    
    interested_channels = set()

    for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {fwd_rev_input_bam}"):
        if not read.has_tag('rq') or not read.has_tag('ch'):  # 确保必要的 tags 存在
            continue

        channel = read.get_tag('ch')
        rq = read.get_tag('rq')
        direction = 'fwd' if read.query_name.endswith('fwd') else 'rev'

        if channel not in channel_reads:
            channel_reads[channel] = {'fwd': 0.0, 'rev': 0.0}

        channel_reads[channel][direction] = rq

    # 遍历收集的 reads 并进行过滤
    for channel, reads in tqdm(channel_reads.items(), desc='Filtering reads'):
        fwd_rq = reads['fwd']
        rev_rq = reads['rev']
        if fwd_rq >= 0.999 or rev_rq >= 0.999:
            # bam_out.write(fwd_read)
            interested_channels.add(channel)

    bam_in.close()
    return interested_channels

def filter_bam(input_bam: str, fwd_rev_bam: str, output_bam: str):
    """
        通过 fwd_rev_bam 将 Q30及以上 的channels 信息抽取出来。
        然后 将 input_bam 中 hit 上述 channels 的 reads 写入 output_bam

    Args:
        input_bam (str): smc --byStrand=false 生成出来的 bam 文件
        fwd_rev_bam (str): smc --byStrand=true 生成出来的 bam 文件
        output_bam (str): 过滤后的结果文件
    """
    bam_in = pysam.AlignmentFile(input_bam, "rb", check_sq=False, threads=40)
    bam_out = pysam.AlignmentFile(output_bam, "wb", check_sq=False, threads=40, header=bam_in.header)

    interested_channels = extract_interested_channels(fwd_rev_bam)

    for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {input_bam}"):
        if not read.has_tag('rq') or not read.has_tag('ch'):  # 确保必要的 tags 存在
            continue
        channel = read.get_tag('ch')
        if channel in interested_channels:
            bam_out.write(read)

    bam_in.close()
    bam_out.close()


def main():
    """
    
    """
    input_bam = '/data/ccs_data/case-study/20250310-lowQ30/Output.smc_all_reads.bam'
    fwd_rev_bam = '/data/ccs_data/case-study/20250310-lowQ30/output-bystrand.smc_all_reads.bam'
    output_bam = '/data/ccs_data/case-study/20250310-lowQ30/output.smc_all_reads.Q30.bam'
    filter_bam(input_bam, fwd_rev_bam, output_bam)


if __name__ == "__main__":
    main()