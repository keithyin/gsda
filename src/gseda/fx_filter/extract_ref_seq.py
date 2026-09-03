from pathlib import Path


def extract_chr9(input_fasta: str, output_fasta: str):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        write = False
        for line in infile:
            if line.startswith(">"):
                if "chr9" in line:
                    write = True
                    outfile.write(line)  # 写入 fasta header
                else:
                    write = False
            elif write:
                outfile.write(line)


# 示例调用
extract_chr9(
    input_fasta="/data/ccs_data/HG002/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1-chr22.fasta",
    output_fasta="/data/REF_GENOMES/chr9.fasta",
)
