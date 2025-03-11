import mappy as mp
import pysam
from collections import defaultdict
from tqdm import tqdm
from typing import List, Tuple

def calculate_identity(cigar_tuple):
    """
    Calculate the identity of a read based on CIGAR tuple with 'eqx' operations.

    Args:
        cigar_tuple (list): A list of (operation, length) tuples from a CIGAR string.

    Returns:
        float: The calculated identity (matches / alignment length).
    """
    matches = 0
    alignment_length = 0
    
    for length, op in cigar_tuple:
        if op == 7:  # Exact match (EQ)
            matches += length
            alignment_length += length
        elif op == 8:  # Mismatch (X)
            alignment_length += length
        elif op in {0, 1, 2}:  # Match (M), Insertion (I), Deletion (D)
            alignment_length += length

    if alignment_length == 0:
        return 0.0
    return matches / alignment_length

def pin_start_end(cigar_tuple: List[Tuple[int, int]]):
    if len(cigar_tuple)>0:
        if cigar_tuple[0][1] != 7:
            cigar_tuple.pop(0)
    if len(cigar_tuple) >0:
        if cigar_tuple[-1][1] != 7:
            cigar_tuple.pop(-1)
    return cigar_tuple

def analyze_bam(bam_file: str):
    """
    Analyze a BAM file for strand inconsistency based on 'channel' tag.

    Args:
        bam_file (str): Path to the input BAM file.
    """
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40)

    # Group reads by 'channel' tag
    reads_by_channel = defaultdict(lambda: {'fwd': None, 'rev': None})

    channels = set()
    
    for read in tqdm(bam.fetch(until_eof=True), desc=f"reading {bam_file}"):
        channel = read.get_tag('channel')
        name = read.query_name
        channels.add(channel)
        if name.endswith('fwd'):
            reads_by_channel[channel]['fwd'] = read.query_sequence
        elif name.endswith('rev'):
            reads_by_channel[channel]['rev'] = read.query_sequence

    bam.close()

    # Align and check inconsistencies
    # inconsistent_reads = []
    fwd_rev_channels = 0
    inconsistent_channels = 0
    for channel, reads in tqdm(reads_by_channel.items(), desc='Analyzing channels'):
        if reads['fwd'] is not None and reads['rev'] is not None:
            fwd_rev_channels += 1
            fwd_seq = reads['fwd']
            rev_seq = reads['rev']

            # Use mappy to align
            aligner = mp.Aligner(seq=fwd_seq, extra_flags = 67108864)  # Default to nucleotide alignment
            
            for hit in aligner.map(rev_seq):
                # print(hit.cigar)
                # print(hit.NM)
                cigar = pin_start_end(hit.cigar)
                identity = calculate_identity(cigar)
                if 0.99 < identity < 0.999:
                    inconsistent_channels += 1
    print("tot:{}, fwd_rev_channels:{}, inconsistent_channels:{}, inconsistent_ratio:{}".format(len(channels), fwd_rev_channels, inconsistent_channels, inconsistent_channels/fwd_rev_channels))

    # Print inconsistent reads
    # if inconsistent_reads:
    #     print('Inconsistent reads found:')
    #     for channel, fwd_name, rev_name in inconsistent_reads:
    #         print(f'Channel: {channel}, Fwd: {fwd_name}, Rev: {rev_name}')
    # else:
    #     print('No inconsistencies found.')
        
    
def main():
    analyze_bam("/data/ccs_data/case-study/20250310-lowQ30/output-bystrand.smc_all_reads.bam")        

if __name__ == "__main__":
    main()