import pysam
from typing import List, Set


class LocusBaseCounter:
    def __init__(self):
        self.ref_base = "-"
        self.inner = {"A": 0, "C": 0, "G":0, "T": 0, "-": 0}
        self.deletion = 0
        self.a = 0
        self.c = 0
        self.g = 0
        self.t = 0

    def update(self, base):
        self.inner[base] += 1
        
    def preety_string(self):
        tot = sum([cnt for _, cnt in self.inner.items()])
        
        items = sorted(list(self.inner.items()), key=lambda x: x[1], reverse=True)
        items = [(item[0], item[1] / tot) for item in items]
        
        info = ""
        for k, v in items:
            info += f"{k}:{v*100:.2f}%\t"
        return f"ref_base:{self.ref_base}\ttot:{tot}\t{info}"


def variant_ratio(bam_file: str, list_locus: Set[int], contig_name: str, ref_seq: str):

    start = min(list_locus)
    stop = max(list_locus) + 1

    all_locus_counters = {locus: LocusBaseCounter() for locus in list_locus}

    with pysam.AlignmentFile(bam_file, mode="rb") as bam_in:

        for plp_col in bam_in.pileup(contig=contig_name, start=start, stop=stop, min_base_quality=1):
            if (
                plp_col.reference_pos < start
                or plp_col.reference_pos >= stop
            ):
                continue
            if plp_col.reference_pos not in list_locus:
                continue
            all_locus_counters[plp_col.reference_pos].ref_base = ref_seq[plp_col.reference_pos]

            for query in plp_col.pileups:
                query_seq = query.alignment.query_sequence
                if query.is_refskip:
                    continue
                if query.is_del:
                    all_locus_counters[plp_col.reference_pos].update(base="-")

                    continue

                qpos = query.query_position
                base = query_seq[qpos]

                all_locus_counters[plp_col.reference_pos].update(base=base)

    all_locus_counters = list(all_locus_counters.items())
    all_locus_counters = sorted(all_locus_counters, key=lambda x: x[0])
    for (locus, counter) in all_locus_counters:
        print(f"locus:{locus}\t{counter.preety_string()}")
    
    
def main():
    list_locus = set([74, 76, 148, 193])
    contig_name = "refN"
    
    ref_seq = "CACAGAGCATCCGAGGACAAAACTTCCATTTGCTGCTGGTGGACGAGGCGCACTTCATCAAGAAAGAGGCCTTCAACACCATTCTGGGGTTCCTGGCCCAGAATACCACCAAGATCATCTTCATATCGTCCACCAACACTACCAGTGACGCCACGTGTTTCTTGACGCGCCTCAACAATGCGCCCTTTGACATGCTCAACGTGGTCTCTTACGTGTGCGAAGAGCACCTGCACAGCTTCACGGAGAAGGGCGACGCCACGGCGTGTCCTTGCTACCGACTGCACAAGCCCACCTTCATCAGCCTCAACTCGCAGGTGCGCAAGACGGCCAACATGTTTATGCCGGGCGCTTTCATGGACGAGATCATCGGCGGTACCAATAAAATCTCGCAGAACACCGTGCTCATCACGGACCAGAGCCGCGAAGAGTTCGATATTTTGCGTTACAGCACGCTCAACACCAACGCCTACGATTATTTCGGCAAGACGCTTTACGTGTATCTGGACCCGGCCTTCACCACCAACCGCAAGGCCTCGGGCACGGGCGTGGCGGCCGTAGGCGCCTACCGACACCAGTTTCTCATTTACGGCCTAGAGCATTTCTTTTTGCGCGACCTCTCCGAGAGTTCTGAGGTAGCCATCGCCGAGTGCGCGGCGCACATGATCATCTCGGTGCTGAGCCTGCACCCTTACCTGGACGAACTGCGTATCGCCGTGGAGGGCAACACCAACCAGGCGGCGGCCGTGCGCATCGCCTGCCTCATCCGACAGAGCGTGCAGAGCAGCACGCTCATCCGCGTGCTCTTCTACCACACGCCCGACCAGAACCACATCGAACAGCCCTTCTACCTCATGGGCCGCGACAAGGCGCTGGCCGTGGAACAGTTCATCTCGCGTTTCAACTCGGGCTACATCAAAGCCTCGCAAGAG"
    bam_files = [
        "/data/ccs_data/20250804-ludaopei/Barcodes_20250804_240601Y0005_Run0002/cmv-b3-hmm.bam",
        "/data/ccs_data/20250804-ludaopei/icing-out/Barcodes_icing_20250804_240601Y0005_Run0002/cmv-b3-icing.bam",
        "/data/ccs_data/20250804-ludaopei/Barcodes_20250804_240601Y0005_Run0002/cmv-b4-hmm.bam",
        "/data/ccs_data/20250804-ludaopei/icing-out/Barcodes_icing_20250804_240601Y0005_Run0002/cmv-b4-icing.bam",
        "/data/ccs_data/20250804-ludaopei/Barcodes_20250804_240601Y0005_Run0002/cmv-b5-hmm.bam",
        "/data/ccs_data/20250804-ludaopei/icing-out/Barcodes_icing_20250804_240601Y0005_Run0002/cmv-b5-icing.bam",
        "/data/ccs_data/20250804-ludaopei/Barcodes_20250804_240601Y0005_Run0002/cmv-b6-hmm.bam",
        "/data/ccs_data/20250804-ludaopei/icing-out/Barcodes_icing_20250804_240601Y0005_Run0002/cmv-b6-icing.bam",
        "/data/ccs_data/20250804-ludaopei/Barcodes_20250804_240601Y0005_Run0002/cmv-b7-hmm.bam",
        "/data/ccs_data/20250804-ludaopei/icing-out/Barcodes_icing_20250804_240601Y0005_Run0002/cmv-b7-icing.bam",
        "/data/ccs_data/20250804-ludaopei/Barcodes_20250804_240601Y0005_Run0002/cmv-b8-hmm.bam",
        "/data/ccs_data/20250804-ludaopei/icing-out/Barcodes_icing_20250804_240601Y0005_Run0002/cmv-b8-icing.bam",
        
        
    ]
    for bam_file in bam_files:
        print(bam_file)
        variant_ratio(bam_file, list_locus, contig_name, ref_seq)

if __name__ == "__main__":
    main()
