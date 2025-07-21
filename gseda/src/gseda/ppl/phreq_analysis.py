import os
import sys
from pathlib import Path

cur_file = os.path.abspath(__file__)
cur_file = Path(cur_file)
sys.path.insert(0, cur_file.parent.as_posix())
sys.path.insert(0, cur_file.parent.parent.as_posix())

import argparse
import subprocess
from fact_table_ana import pred_baseq_and_emp_q, rq_iy_analysis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", type=str, required=True, help="query file")
    parser.add_argument("-t", type=str, required=True, help="target file")
    parser.add_argument("-f", action="store_true", help="force to overwrite the output bam file")
    parser.add_argument("--qname-suffix", type=str, default=None, help="only the qname endswith $qname_suffix will be considered", dest="qname_suffix")

    args = parser.parse_args()

    """
    gsmm2 align -q $query_file -t $ref_file -p outputbam_prefix --noMar
    gsetl --outdir $outdir aligned-bam --bam $aligned_bam --ref-file $ref_file
    """

    file_stem = Path(args.q).stem

    oup_bam_prefix = f"{args.q}.gsmm2-aligned"
    oup_bam = f"{oup_bam_prefix}.bam"
    if args.f and os.path.exists(oup_bam):
        os.remove(oup_bam)
    if not os.path.exists(oup_bam):
        ## 1. do alignment
        cmd_str = f"gsmm2 align -q {args.q} -t {args.t} -p {oup_bam_prefix} --noMar"
        if args.qname_suffix:
            cmd_str += f" --qname-suffix {args.qname_suffix}"
        subprocess.check_call(cmd_str, shell=True)

    ## 2. gsetl extract fact table
    gsetl_outdir = f"{oup_bam}-gsetl"
    cmd_str = f"gsetl -f --outdir {gsetl_outdir} aligned-bam --bam {oup_bam} --ref-file {args.t} --factPolyInfo 0 "
    subprocess.check_call(cmd_str, shell=True)

    ## 3. do analysis

    baseq_ana_args = {
        "fact_table": f"{gsetl_outdir}/fact_baseq_stat.csv",
        "o_prefix": file_stem,
    }
    pred_baseq_and_emp_q.main(argparse.Namespace(**baseq_ana_args))

    readsq_ana_args = {
        "fact_table": f"{gsetl_outdir}/fact_aligned_bam_bam_basic.csv",
        "o_prefix": file_stem,
    }
    rq_iy_analysis.main(argparse.Namespace(**readsq_ana_args))


if __name__ == "__main__":
    main()
