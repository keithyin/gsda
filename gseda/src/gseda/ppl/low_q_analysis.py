
import os
import argparse
import pathlib
import multiprocessing as mp
import sys
cur_path = pathlib.Path(os.path.abspath(__file__))
cur_dir = cur_path.parent
prev_dir = cur_path.parent.parent
prev_prev_dir = cur_dir.parent.parent.parent
sys.path.append(str(cur_dir))
sys.path.append(str(prev_dir))
sys.path.append(str(prev_prev_dir))

print(f"cur_dir:{cur_dir}")
print(f"prev_dir:{prev_dir}")
print(f"prev_prev_dir:{prev_prev_dir}")
print(f"sys.path={sys.path}")

import homo_and_str_region_coverage  # noqa: E402
import dna_strand_inconsistency_analysis  # noqa: E402
import asts_analysis  # noqa: E402
import inverted_repeat_ratio  # noqa: E402


def main(args):
    """
    低Q20 Q30 yield 的问题分析脚本。
    当前进行以下分析：
        1. 是否存在大量正反链不一致
        2. subreads 2 smc 的比对结果是不是很差
        3. 是否大量的 homo+str

    Args:
        args (_type_): 参数
    """
    sbr_bam = args.sbr_bam
    smc_bam = args.smc_bam

    report = "\n\n=================Summary=================\n\n"

    asts_param = {
        "sbr_bam": sbr_bam,
        "smc_bam": smc_bam,
        "rq_range": "0.95:1.0",
        "np_range": "5:1000000"
    }
    asts_param = argparse.Namespace(**asts_param)
    report += asts_analysis.main(asts_param)

    dna_strand_in_consisteny_param = {
        "smc_bam": smc_bam,
        "sbr_bam": sbr_bam,
        "np_thr": 14,
        "rq_thr": 0.95
    }
    dna_strand_in_consisteny_param = argparse.Namespace(
        **dna_strand_in_consisteny_param)
    report += dna_strand_inconsistency_analysis.main(
        dna_strand_in_consisteny_param)

    homo_str_regin_cov_param = {
        "bam_file": [smc_bam],
        "rq_thr": 0.95
    }
    homo_str_regin_cov_param = argparse.Namespace(**homo_str_regin_cov_param)
    report += homo_and_str_region_coverage.main(homo_str_regin_cov_param)

    # ir_ratio_param = {
    #     ""
    # }

    report += inverted_repeat_ratio.main(bam_path=smc_bam,
                                         threads=mp.cpu_count())

    print(report)
    pass


def main_cli():
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("--sbr-bam", dest="sbr_bam")
    parser.add_argument("--smc-bam", dest="smc_bam")

    args = parser.parse_args()
    main(args=args)
    pass


if __name__ == "__main__":
    main_cli()
    pass
