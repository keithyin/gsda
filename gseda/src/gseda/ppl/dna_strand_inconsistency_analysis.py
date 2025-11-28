
import subprocess
import pathlib
import os
import argparse


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
    smc_bam = args.smc_bam
    sbr_bam = args.sbr_bam

    sbr2smc_alignment = asts_alignment(
        sbr_bam, smc_bam, np_thr=args.np_thr, rq_thr=args.rq_thr)

    cmd = f"base_mismatch_identification v2 --smc-bam {smc_bam} --sbr2smc-bam {sbr2smc_alignment} --origin-bam-rq-thr {args.rq_thr} --origin-bam-np-thr {args.np_thr}"
    subprocess.check_call(cmd, shell=True)
    if os.path.exists(sbr2smc_alignment):
        os.remove(sbr2smc_alignment)

    pass


if __name__ == "__main__":
    main()
