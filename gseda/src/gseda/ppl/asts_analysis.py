import subprocess
import argparse
import pathlib
import pysam
from multiprocessing import cpu_count
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np


def asts_bam_identity(asts_bam_file: str):
    identity_count = {}
    coverage_count = {}
    identity_m_coverage_count = {}

    with pysam.AlignmentFile(asts_bam_file, mode="rb", threads=cpu_count()) as in_bam:
        for record in tqdm(in_bam.fetch(), desc=f"reading {asts_bam_file}"):
            identity = record.get_tag("iy")
            coverage = record.get_tag("ec")
            id_cv = identity * coverage

            identity = int(identity * 100)
            coverage = int(coverage * 100)
            id_cv = int(id_cv * 100)
            # q = int(round((-10 * math.log10(1 - identity))))
            identity_count.setdefault(identity, 0)
            identity_count[identity] += 1

            coverage_count.setdefault(coverage, 0)
            coverage_count[coverage] += 1

            identity_m_coverage_count.setdefault(id_cv, 0)
            identity_m_coverage_count[id_cv] += 1

    return identity_count, coverage_count, identity_m_coverage_count


def _probabilities(counts):
    # 将 {value: frequency} 换算成 0-100 全量概率向量
    all_values = list(range(0, 101))
    total_count = sum([v for _, v in counts.items()])
    frequencies = [counts.get(value, 0) for value in all_values]
    return all_values, [freq/total_count for freq in frequencies]


def top5_report(counts, tag):
    # 统计：返回频次最高的前 5 个值的报告字符串
    all_values, probabilities = _probabilities(counts)
    top5_indices = sorted(range(len(probabilities)),
                          key=lambda i: probabilities[i],
                          reverse=True)[:5]

    report_str_inner = ""
    # 打印结果
    for i, idx in enumerate(top5_indices):
        prob = probabilities[idx]
        value = all_values[idx]  # 获取对应的原始值
        inner = f"- {tag} Rank: {i+1}, Value={value}, Freq={prob*100:.2f}% \n"
        report_str_inner += inner
    return report_str_inner


def draw_axis(counts, ax, tag):
    # 绘图：在 ax 上画出频数直方图
    all_values, probabilities = _probabilities(counts)
    # 频数直方图
    bars1 = ax.bar(all_values, probabilities, color='skyblue',
                   edgecolor='navy', alpha=0.7, width=0.8)
    ax.set_xlabel('(0-100)', fontsize=12)
    ax.set_ylabel('Freq', fontsize=12)
    ax.set_title(tag, fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    ax.set_xticks(range(0, 100, 5))  # 每5个刻度显示一个
    ax.set_xlim(0, 100)

    # 在频数柱子上方显示数值（只显示非零值）
    for bar, freq in zip(bars1, probabilities):
        if freq > 0:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(probabilities)*0.01,
                    f'{freq:.2f}', ha='center', va='bottom', fontsize=8)
    return ax


def threshold_report(counts, tag, thresholds):
    # 统计：value 分别大于各阈值的 read 占比（累计比例）
    total_count = sum([v for _, v in counts.items()])
    lines = []
    for t in thresholds:
        pass_count = sum(v for k, v in counts.items() if k > t)
        ratio = pass_count / total_count if total_count else 0.0
        lines.append(f"- {tag} > {t}: {pass_count}/{total_count} ({ratio*100:.2f}%) \n")
    return "".join(lines)


def build_report(identity_count, coverage_count, id_cv_count):
    # 组装文本报告
    total_count = sum([v for _, v in identity_count.items()])
    report_str_inner = f"- Subreads Processed: {total_count} \n"

    report_str_inner += top5_report(identity_count, "Identity")
    report_str_inner += threshold_report(
        identity_count, "Identity", [91, 90, 89, 88, 87, 86, 85, 75, 50])
    report_str_inner += "- \n"

    report_str_inner += top5_report(coverage_count, "Coverage")
    report_str_inner += threshold_report(
        coverage_count, "Coverage",
        [99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 75, 50])
    report_str_inner += "- \n"

    report_str_inner += top5_report(id_cv_count, "Identity * Coverage")

    return f"""
================= Asts Report =================
{report_str_inner}
===============================================
"""


def plot(identity_count, coverage_count, id_cv_count, dest_filepath):
    # 创建子图
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 10))

    draw_axis(identity_count, ax1, tag="Identity")
    draw_axis(coverage_count, ax2, tag="Coverage")
    draw_axis(id_cv_count, ax3, tag="Identity * Coverage")

    plt.tight_layout()
    plt.savefig(dest_filepath)
    return


def main_cli():
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("sbr_bam")
    parser.add_argument("smc_bam")
    parser.add_argument("--np-range", type=str,
                        default="5:100000000", dest="np_range", help="start:end")
    parser.add_argument("--rq-range", type=str,
                        default="0.95:1.0", dest="rq_range", help="start:end")

    args = parser.parse_args()
    main(args=args)
    pass


def smc_has_stranded_reads(smc_bam: str) -> bool:
    with pysam.AlignmentFile(smc_bam, mode="rb") as bam:
        for record in bam:
            name = record.query_name
            if name.endswith("fwd") or name.endswith("rev"):
                return True
    return False


def concat_bams(bam_files, out_bam: str):
    cmd = "samtools cat " + " ".join(bam_files) + f" -o {out_bam}"
    subprocess.check_call(cmd, shell=True)
    subprocess.check_call(f"samtools index {out_bam}", shell=True)


def main(args):

    subreads_bam = args.sbr_bam
    smc_bam = args.smc_bam

    smc_path = pathlib.Path(smc_bam)
    prefix = smc_path.parent.joinpath(f"{smc_path.stem}.asts")

    shared_flags = (f"--np-range {args.np_range} "
                    f"--rq-range {args.rq_range} "
                    f"--ptTags dw")

    if smc_has_stranded_reads(smc_bam):
        for strand, flag in (("fwd", "--fwd-only"), ("rev", "--rev-only")):
            prefix_strand = f"{prefix}.{strand}"
            asts_cmd = (f"asts -q {subreads_bam} -t {smc_bam} "
                        f"-p {prefix_strand} {flag} {shared_flags}")
            print(f"running: {asts_cmd}")
            subprocess.check_call(asts_cmd, shell=True)
            print("")
        concat_bams([f"{prefix}.fwd.bam", f"{prefix}.rev.bam"],
                    f"{prefix}.bam")
    else:
        asts_cmd = (f"asts -q {subreads_bam} -t {smc_bam} -p {prefix} "
                    f"{shared_flags}")
        print(f"running: {asts_cmd}")
        subprocess.check_call(asts_cmd, shell=True)
        print("")

    identity_count, coverage_count, id_cv_count = asts_bam_identity(
        f"{prefix}.bam")
    dest_img_path = f"{prefix}.q_dist.jpg"

    plot(identity_count=identity_count, coverage_count=coverage_count,
         id_cv_count=id_cv_count, dest_filepath=dest_img_path)
    print(f"dest_image_path: {dest_img_path}")

    report_str = build_report(identity_count=identity_count,
                              coverage_count=coverage_count,
                              id_cv_count=id_cv_count)
    print(report_str)
    return report_str


if __name__ == "__main__":
    main_cli()
