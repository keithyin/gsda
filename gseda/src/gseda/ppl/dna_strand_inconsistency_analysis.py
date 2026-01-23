
import subprocess
import pathlib
import os
import argparse
import json

def asts_alignment(query_bam: str, target_bam: str, np_thr: int, rq_thr: float):
    o_dir = os.path.dirname(target_bam)
    o_name = "{}-TO-{}.align".format(pathlib.Path(query_bam).stem,
                                     pathlib.Path(target_bam).stem)
    o_prefix = f"{o_dir}/{o_name}"

    o_name = f"{o_prefix}.bam"

    if os.path.exists(o_name):
        return o_name

    cmd = f"asts -q {query_bam} -t {target_bam} -p {o_prefix} --np-range {np_thr}:10000000 --rq-range {rq_thr}:1.1"
    subprocess.check_call(cmd, shell=True)
    return o_name


def main(args):

    smc_bam = args.smc_bam
    sbr_bam = args.sbr_bam

    sbr2smc_alignment = asts_alignment(
        sbr_bam, smc_bam, np_thr=args.np_thr, rq_thr=args.rq_thr)

    cmd = f"base_mismatch_identification v2 --smc-bam {smc_bam} --sbr2smc-bam {sbr2smc_alignment} --origin-bam-rq-thr {args.rq_thr} --origin-bam-np-thr {args.np_thr}"
    subprocess.check_call(cmd, shell=True)
    
    smc_bam_path = pathlib.Path(smc_bam)
    stem = smc_bam_path.stem
    root = smc_bam_path.parent
    result_json_path = root.joinpath(f"{stem}.strand-consistency.json")
    with open(result_json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    inconsistence_ratio = data["inconsistence_ratio"]
    
    if os.path.exists(sbr2smc_alignment):
        os.remove(sbr2smc_alignment)

    report_str = f"""
================= DNA Strand Inconsistency Analysis =================
- InconsistenceRatio: {inconsistence_ratio * 100:.2f}%
=====================================================================
"""

    return report_str
    pass


def main_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smc-bam", type=str, required=True,
                        dest="smc_bam", help="byStrand=False bam file")
    parser.add_argument("--sbr-bam", type=str, required=True,
                        dest="sbr_bam", help="subreads bam file")
    parser.add_argument("--np-thr", type=int, dest="np_thr", default=14,
                        help="only the channel that np≥np_thr will be processed")
    parser.add_argument("--rq-thr", type=float, dest="rq_thr",
                        default=0.95,
                        help="only the channel that rq≥rq_thr will be processed")

    args = parser.parse_args()
    main(args=args)


if __name__ == "__main__":
    main_cli()
