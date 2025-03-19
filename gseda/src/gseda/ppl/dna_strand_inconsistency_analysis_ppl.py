
import pysam
from tqdm import tqdm
from typing import Dict, Set, Tuple, List
import subprocess
import numpy as np
import pathlib
import os


class InterestLocus:
    def __init__(self, pos, pos_base):
        self.pos = pos
        self.pos_base = pos_base
        self.fwd_pos = None
        self.fwd_base = None
        
        self.fwd_sbr_top1_base = None
        
        
        # reverse consensus 对应的 base
        self.rev_pos = None
        self.rev_base = None
        self.rev_sbr_top1_base = None
        
    def set_fwd(self, pos, base):
        self.fwd_pos = pos
        self.fwd_base = base
    
    def set_rev(self, pos, base):
        self.rev_pos = pos
        self.rev_base = base
        
    def valid(self):
        return self.fwd_pos is not None and self.rev_pos is not None
    
    def valid2(self):
        return self.fwd_sbr_top1_base is not None and self.rev_sbr_top1_base is not None
    
    def inconsistent(self):
        assert self.valid() and self.valid2()
        if self.fwd_base == self.fwd_sbr_top1_base and self.rev_base == self.rev_sbr_top1_base:
            if self.fwd_base != self.rev_base:
                return True
        return False


class Qual:
    def __init__(self, qual):
        self.qual = qual
    def get_baseq(self, pos, rev):
        pos = pos if not rev else len(self.qual) - pos - 1
        return self.qual[pos]

class Read:
    def __init__(self, record: pysam.AlignedSegment):
        self.name = record.query_name
        self.seq = record.query_sequence
        self.qual = np.array(record.query_qualities.tolist())
    

def asts_alignment(query_bam: str, target_bam: str):
    o_dir = os.path.dirname(target_bam)
    o_name = "{}-TO-{}.align".format(pathlib.Path(query_bam).stem, pathlib.Path(target_bam).stem)
    o_prefix = f"{o_dir}/{o_name}"
    
    o_name = f"{o_prefix}.bam"
    
    if os.path.exists(o_name):
        return o_name
    
    cmd = f"asts -q {query_bam} -t {target_bam} -p {o_prefix}"
    subprocess.check_call(cmd, shell=True)
    return o_name


def extract_interested_channels(fwd_rev_input_bam: str, threshold: float = 0.999):
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
        if fwd_rq >= threshold or rev_rq >= threshold:
            # bam_out.write(fwd_read)
            interested_channels.add(channel)

    bam_in.close()
    return interested_channels


def dump_smc_bam(input_bam: str, infix: str, channels: Set[int]) -> str:    
    output_bam = "{}.{}.bam".format(input_bam.rsplit(".", maxsplit=1)[0], infix)
    bam_in = pysam.AlignmentFile(input_bam, "rb", check_sq=False, threads=40)
    bam_out = pysam.AlignmentFile(output_bam, "wb", check_sq=False, threads=40, header=bam_in.header)

    for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {input_bam}"):
        if not read.has_tag('rq') or not read.has_tag('ch'):  # 确保必要的 tags 存在
            continue
        channel = read.get_tag('ch')
        if channel in channels:
            bam_out.write(read)

    bam_in.close()
    bam_out.close()
    return output_bam


def seperate_fwd_rev_bam(input_bam: str) -> Tuple[str, str]:
    fwd_bam = "{}.fwd.bam".format(input_bam.rsplit(".", maxsplit=1)[0])
    rev_bam = "{}.rev.bam".format(input_bam.rsplit(".", maxsplit=1)[0])
    bam_in = pysam.AlignmentFile(input_bam, "rb", check_sq=False, threads=40)
    fwd_bam_out = pysam.AlignmentFile(fwd_bam, "wb", check_sq=False, threads=40, header=bam_in.header)
    rev_bam_out = pysam.AlignmentFile(rev_bam, "wb", check_sq=False, threads=40, header=bam_in.header)

    for read in tqdm(bam_in.fetch(until_eof=True), desc=f"reading {input_bam}"):
        if read.query_name.endswith('fwd'):
            fwd_bam_out.write(read)
        elif read.query_name.endswith('rev'):
            rev_bam_out.write(read)

    bam_in.close()
    fwd_bam_out.close()
    rev_bam_out.close()
    return fwd_bam, rev_bam


def process_read_worker(
    batch: List[Tuple[str, 'Read']], 
    fwd_rev_name_to_qual: Dict[str, 'Qual'], 
    fwd_rev2ori_bam: str
):
    interested_loci_batch = {}

    with pysam.AlignmentFile(fwd_rev2ori_bam, "rb", threads=1) as fwd_rev2ori:
        for name, read in batch:
            indices, = np.nonzero((read.qual < 15).astype(np.int32))
            indices_set = set(indices.tolist())
            if len(indices_set) == 0:
                continue

            interested_loci = []
            for plp_col in fwd_rev2ori.pileup(
                name, 
                start=max(indices.min() - 10, 0), 
                end=min(indices.max() + 10, len(read.seq)), 
                min_base_quality=1
            ):
                if plp_col.reference_pos not in indices_set:
                    continue

                interested_locus = InterestLocus(plp_col.reference_pos, read.seq[plp_col.reference_pos])

                for pileupread in plp_col.pileups:
                    if pileupread.alignment.is_reverse or pileupread.is_del:
                        continue

                    query_pos = pileupread.query_position
                    query_name = pileupread.alignment.query_name
                    baseq = fwd_rev_name_to_qual[query_name].get_baseq(query_pos, pileupread.alignment.is_reverse)

                    if baseq < 25:
                        continue

                    if query_name.endswith("fwd"):
                        interested_locus.set_fwd(query_pos, pileupread.alignment.query_sequence[query_pos])
                    else:
                        interested_locus.set_rev(query_pos, pileupread.alignment.query_sequence[query_pos])

                if interested_locus.valid():
                    interested_loci.append(interested_locus)

            if interested_loci:
                interested_loci_batch[name] = interested_loci

    return interested_loci_batch

def main():
    """
    
    1. extract Q30 channels from fwd_rev.bam(--byStrand=True) & dump them into a file
    2. dump Q30 channels from origin.bam(--byStrand=False)
    3. seperate Q30 channels in fwd_rev.bam into fwd.bam & rev.bam
    4. alignment
        1. align fwd_rev.q30.bam to origin.q30.bam
        2. align sbr.bam to fwd.q30.bam and and rev.q30.bam
    5. using pileup to analysis result    
    """
    origin_bam = "/data/ccs_data/case-study/20250310-lowQ30/Output.smc_all_reads.bam"
    fwd_rev_bam = "/data/ccs_data/case-study/20250310-lowQ30/output-bystrand.smc_all_reads.bam"
    sbr_bam = "/data/ccs_data/case-study/20250310-lowQ30/Output.subreads.bam"
    
    interested_channels = extract_interested_channels(fwd_rev_bam)
    origin_bam = dump_smc_bam(origin_bam, "q30", interested_channels)
    fwd_rev_bam = dump_smc_bam(fwd_rev_bam, "q30", interested_channels)
    (fwd_bam, rev_bam) = seperate_fwd_rev_bam(fwd_rev_bam)
    
    fwd_rev2ori_bam = asts_alignment(fwd_rev_bam, origin_bam)
    sbr2fwd_bam = asts_alignment(sbr_bam, fwd_bam)
    sbr2rev_bam = asts_alignment(sbr_bam, rev_bam)
    
    pass

if __name__ == "__main__":
    main()