import bam_concordance
from typing import List, Tuple
import polars as pl
import subprocess
import logging
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import pathlib

filepath = os.path.abspath(__file__)
cur_dir = filepath.rsplit("/", maxsplit=1)[0]
sys.path.append(cur_dir)


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def polars_env_init():
    os.environ["POLARS_FMT_TABLE_ROUNDED_CORNERS"] = "1"
    os.environ["POLARS_FMT_MAX_COLS"] = "100"
    os.environ["POLARS_FMT_MAX_ROWS"] = "300"
    os.environ["POLARS_FMT_STR_LEN"] = "100"


def set_tick_size(ax):
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    ax.xaxis.label.set_size(14)
    ax.yaxis.label.set_size(14)


def q2phreq_expr(inp_name, oup_name=None):
    oup_name = oup_name if oup_name is not None else inp_name
    return (
        -10.0
        * (
            1
            - pl.when(pl.col(inp_name) > (1 - 1e-10))
            .then(1 - 1e-10)
            .otherwise(pl.col(inp_name))
        ).log10()
    ).alias(oup_name)


def extract_filename(filepath: str) -> str:
    p = pathlib.Path(filepath)
    return p.stem


def do_alignment(
    bam_file: str, ref_fasta: str, outdir: str, force: bool = False, hifi=False
) -> str:
    rq_t = 0
    if hifi:
        rq_t = 0.99

    res_bam_prefix = "{}/{}.aligned".format(outdir, extract_filename(bam_file))
    result_bam = f"{res_bam_prefix}.bam"

    if force and os.path.exists(result_bam):
        os.remove(result_bam)

    preset = "map-hifi" if hifi else "map-ont"
    if not os.path.exists(result_bam):
        cmd = f"""gsmm2 --preset {preset} align -q {bam_file} -t {ref_fasta} -p {res_bam_prefix} --noMar --rq-range {rq_t}:1.1"""

        logging.info("cmd: %s", cmd)
        subprocess.check_call(cmd, shell=True)
    else:
        logging.info(f"{result_bam} exists, skip the alignment procedure")

    return result_bam


def align_and_stat(
    bam_file: str,
    ref_fasta: str,
    outdir: str,
    force: bool = False,
    rq_thr=0,
    hc_regions=None,
    hc_variants=None,
) -> str:
    res_csv_path = "{}/{}.aligned.metric.csv".format(
        outdir, extract_filename(bam_file))
    if force and os.path.exists(res_csv_path):
        os.remove(res_csv_path)

    aligned_bam = do_alignment(
        bam_file=bam_file,
        ref_fasta=ref_fasta,
        force=force,
        hifi=rq_thr > 0.98,
        outdir=outdir,
    )

    if not os.path.exists(res_csv_path):
        cmd = f"bam_stat bam-concordance {ref_fasta} {aligned_bam} "
        if hc_regions is not None:
            cmd += f" --hcregions {hc_regions}"
        if hc_variants is not None:
            cmd += f" --hcvariants {hc_variants}"
        logging.info(f"runing {cmd}")
        subprocess.check_call(cmd, shell=True)

    else:
        logging.info(
            f"{res_csv_path} exists, skip the bam_concordance procedure")

    return res_csv_path


def channel_level_aggr(
    csv_file: str, idx: int, np_thr: int = 25, sbr_qv_thr=0, rq_thr=None
) -> pl.DataFrame:
    logging.info(f"reading {csv_file}")
    df = pl.read_csv(csv_file, separator="\t")
    df = (
        df.select(
            pl.col("channel_id").cast(pl.Int64),
            pl.col("readLengthBp").cast(pl.Int64),
            pl.col("subreadPasses").cast(pl.Int64),
            pl.col("predictedConcordance").cast(pl.Float32),
            pl.col("concordance").cast(pl.Float32),
            pl.col("concordanceQv").cast(pl.Float32),
            pl.col("matchBp").cast(pl.Int64),
            pl.col("mismatchBp").cast(pl.Int64),
            pl.col("nonHpInsertionBp").cast(pl.Int64),
            pl.col("nonHpDeletionBp").cast(pl.Int64),
            pl.col("hpInsertionBp").cast(pl.Int64),
            pl.col("hpDeletionBp").cast(pl.Int64),
        )
        .with_columns(
            [
                (pl.col("nonHpInsertionBp") +
                 pl.col("hpInsertionBp")).alias("insBp"),
                (pl.col("nonHpDeletionBp") + pl.col("hpDeletionBp")).alias("delBp"),
            ]
        )
        .with_columns(
            [
                (
                    pl.col("matchBp")
                    + pl.col("mismatchBp")
                    + pl.col("insBp")
                    + pl.col("delBp")
                ).alias("alignedSpanLen")
            ]
        )
    )

    if rq_thr is not None:
        df = df.filter(pl.col("predictedConcordance") >= rq_thr)

    if idx == 0:  # surbeads
        df = (
            df.group_by(["channel_id"])
            .agg(
                [
                    pl.col("readLengthBp").sum(),
                    pl.col("matchBp").sum(),
                    pl.len().alias("subreadPasses"),
                    pl.col("predictedConcordance").mean(),
                    # pl.col("concordance").mean(),
                    # pl.col("concordanceQv").mean(),
                    pl.col("mismatchBp").sum(),
                    pl.col("nonHpInsertionBp").sum(),
                    pl.col("nonHpDeletionBp").sum(),
                    pl.col("hpInsertionBp").sum(),
                    pl.col("hpDeletionBp").sum(),
                    pl.col("insBp").sum(),
                    pl.col("delBp").sum(),
                    pl.col("alignedSpanLen").sum(),
                ]
            )
            .with_columns(
                [
                    (
                        pl.col("matchBp")
                        / (
                            pl.col("matchBp")
                            + pl.col("mismatchBp")
                            + pl.col("insBp")
                            + pl.col("delBp")
                        )
                    ).alias("concordance")
                ]
            )
            .with_columns(
                [q2phreq_expr(inp_name="concordance",
                              oup_name="concordanceQv")]
            )
            .filter(pl.col("concordanceQv") >= sbr_qv_thr)
        )

    df = df.with_columns(
        [
            pl.when(pl.col("subreadPasses") > np_thr)
            .then(pl.lit(25))
            .otherwise(pl.col("subreadPasses"))
            .alias("subreadPasses")
        ]
    )
    return df


def np_level_aggr(channel_level_df: pl.DataFrame) -> pl.DataFrame:

    df = (
        channel_level_df.group_by(["subreadPasses"])
        .agg(
            [
                pl.len().alias("numChannels"),
                pl.col("concordanceQv").min().alias("Qv_min"),
                pl.quantile("concordanceQv", quantile=0.25).alias("Qv_25"),
                pl.quantile("concordanceQv", quantile=0.5).alias("Qv_50"),
                pl.quantile("concordanceQv", quantile=0.75).alias("Qv_75"),
                pl.col("concordanceQv").max().alias("Qv_max"),
            ]
        )
        .sort(["subreadPasses"], descending=[False])
    )

    return df


def joined_np_level_aggr(
    channel_level_dfs: List[pl.DataFrame], primary_np_idx: int
) -> pl.DataFrame:
    joined_df, not_none_indices = _join_dfs(
        channel_level_dfs=channel_level_dfs)

    stats_exprs = [pl.len().alias("numChannels")]
    for idx in not_none_indices:
        qv_name = f"concordanceQv_{idx}"
        expr = pl.concat_str(
            pl.col(qv_name).min().round(2),
            pl.quantile(qv_name, quantile=0.01).round(2),
            pl.quantile(qv_name, quantile=0.05).round(2),
            pl.quantile(qv_name, quantile=0.25).round(2),
            pl.quantile(qv_name, quantile=0.50).round(2),
            pl.col(qv_name).max().alias("Qv_max").round(2),
            separator=", ",
        ).alias(f"Qv{idx}_0_1_5_25_50_100")

        stats_exprs.append(expr)
        stats_exprs.append(
            pl.col(f"concordance_{idx}").median().alias(f"M_r_Median_{idx}")
        )

    df = (
        joined_df.rename({f"subreadPasses_{primary_np_idx}": "subreadPasses"})
        .group_by(["subreadPasses"])
        .agg(stats_exprs)
        .sort(["subreadPasses"], descending=[False])
    )

    return df


def global_aggr(channel_level_df: pl.DataFrame) -> pl.DataFrame:
    return (
        channel_level_df.select(
            pl.len().alias("numChannels"),
            pl.col("readLengthBp").sum().alias("Tp"),
            pl.col("matchBp").sum().alias("ETp"),
            pl.col("alignedSpanLen").sum().alias("TotAlignedSpan"),
            pl.col("nonHpInsertionBp").sum().alias("NHIns"),
            pl.col("hpInsertionBp").sum().alias("HIns"),
            pl.col("nonHpDeletionBp").sum().alias("NHDel"),
            pl.col("hpDeletionBp").sum().alias("HDel"),
            pl.col("mismatchBp").sum().alias("MM"),
            (pl.col("concordanceQv") >= 30).sum().alias("≥Q30"),
            (pl.col("concordanceQv") >= 30).mean().round(4).alias("≥Q30_R"),
            (pl.col("matchBp") / pl.col("alignedSpanLen")
             ).median().alias("M_r_Median"),
        )
        .with_columns(
            [
                (pl.col("NHIns") / pl.col("TotAlignedSpan")).alias("NHIns_R"),
                (pl.col("HIns") / pl.col("TotAlignedSpan")).alias("HIns_R"),
                (pl.col("NHDel") / pl.col("TotAlignedSpan")).alias("NHDel_R"),
                (pl.col("HDel") / pl.col("TotAlignedSpan")).alias("HDel_R"),
                (pl.col("MM") / pl.col("TotAlignedSpan")).alias("MM_R"),
                (pl.col("ETp") / pl.col("TotAlignedSpan")).alias("M_R"),
                (pl.col("ETp") / pl.col("Tp")).alias("Valid_R"),
                (
                    (pl.col("ETp") + pl.col("MM") +
                     pl.col("NHIns") + pl.col("HIns"))
                    / pl.col("Tp")
                ).alias("Cov_R"),
            ]
        )
        .drop(["NHIns", "HIns", "NHDel", "HDel", "MM"])
        .with_columns([(-10.0 * (1.0 - pl.col("M_R")).log10()).alias("phreq")])
    )


def joined_global_aggr(channel_level_dfs: List[pl.DataFrame]) -> pl.DataFrame:
    joined_df, not_none_indices = _join_dfs(
        channel_level_dfs=channel_level_dfs)
    stat_exprs = [pl.len().alias("numChannels")]

    qv_median_name_fmt = "QvMedian_{}"
    be_q30_name_fmt = "≥Q30_{}"
    tp_etp_name_fmt = "Tp,ETp_{}"

    for idx in not_none_indices:
        tp_name = f"readLengthBp_{idx}"
        e_tp_name = f"matchBp_{idx}"
        aligned_span_len_name = f"alignedSpanLen_{idx}"
        qv_name = f"concordanceQv_{idx}"
        exprs = [
            pl.sum(tp_name),
            pl.sum(e_tp_name),
            pl.sum(aligned_span_len_name),
            pl.col(qv_name).median().round(2).alias(
                qv_median_name_fmt.format(idx)),
            (pl.col(qv_name).ge(pl.lit(30)).sum() / pl.len())
            .round(2)
            .alias(be_q30_name_fmt.format(idx)),
        ]
        stat_exprs.extend(exprs)

    aggr_df = joined_df.select(*stat_exprs)

    stat_exprs = [pl.col("numChannels")]
    for idx in not_none_indices:
        tp_name = f"readLengthBp_{idx}"
        e_tp_name = f"matchBp_{idx}"
        aligned_span_len_name = f"alignedSpanLen_{idx}"
        exprs = [
            pl.concat_str(pl.col(e_tp_name), pl.col(tp_name), separator=", ").alias(
                tp_etp_name_fmt.format(idx)
            ),
            pl.col(qv_median_name_fmt.format(idx)),
            pl.col(be_q30_name_fmt.format(idx)),
            (
                pl.lit(-10)
                * (
                    pl.lit(1)
                    - pl.col(e_tp_name) / pl.col(aligned_span_len_name)
                    + pl.lit(1e-8)
                ).log10()
            ).alias(f"Qv_{idx}"),
        ]
        stat_exprs.extend(exprs)

    return aggr_df.select(*stat_exprs)


def joined_global_aggr_error_info(
    channel_level_dfs: List[pl.DataFrame],
) -> pl.DataFrame:
    joined_df, not_none_indices = _join_dfs(
        channel_level_dfs=channel_level_dfs)
    stat_exprs = [pl.len().alias("numChannels")]

    qv_median_name_fmt = "QvMedian_{}"
    be_q30_name_fmt = "≥Q30_{}"

    """
    pl.col("mismatchBp").sum(),
    pl.col("nonHpInsertionBp").sum(),
    pl.col("nonHpDeletionBp").sum(),
    pl.col("hpInsertionBp").sum(),
    pl.col("hpDeletionBp").sum(),
    """

    for idx in not_none_indices:
        tp_name = f"readLengthBp_{idx}"
        e_tp_name = f"matchBp_{idx}"
        mm_bp_name = f"mismatchBp_{idx}"
        non_hp_ins_bp_name = f"nonHpInsertionBp_{idx}"
        non_hp_del_bp_name = f"nonHpDeletionBp_{idx}"
        hp_ins_bp_name = f"hpInsertionBp_{idx}"
        hp_del_bp_name = f"hpDeletionBp_{idx}"

        aligned_span_len_name = f"alignedSpanLen_{idx}"
        qv_name = f"concordanceQv_{idx}"
        exprs = [
            pl.sum(tp_name),
            pl.sum(e_tp_name),
            pl.sum(aligned_span_len_name),
            pl.sum(mm_bp_name),
            pl.sum(non_hp_ins_bp_name),
            pl.sum(non_hp_del_bp_name),
            pl.sum(hp_ins_bp_name),
            pl.sum(hp_del_bp_name),
            pl.col(qv_name).median().round(2).alias(
                qv_median_name_fmt.format(idx)),
            (pl.col(qv_name).ge(pl.lit(30)).sum() / pl.len())
            .round(2)
            .alias(be_q30_name_fmt.format(idx)),
        ]
        stat_exprs.extend(exprs)

    aggr_df = joined_df.select(*stat_exprs)

    stat_exprs = [pl.col("numChannels")]
    for idx in not_none_indices:
        tp_name = f"readLengthBp_{idx}"
        e_tp_name = f"matchBp_{idx}"
        aligned_span_len_name = f"alignedSpanLen_{idx}"
        exprs = [
            (pl.col(non_hp_ins_bp_name) / pl.col(aligned_span_len_name)).alias(
                f"NHIns_{idx}"
            ),
            (pl.col(hp_ins_bp_name) / pl.col(aligned_span_len_name)).alias(
                f"HIns_{idx}"
            ),
            (pl.col(non_hp_del_bp_name) / pl.col(aligned_span_len_name)).alias(
                f"NHDel_{idx}"
            ),
            (pl.col(hp_del_bp_name) / pl.col(aligned_span_len_name)).alias(
                f"HDel_{idx}"
            ),
            (pl.col(mm_bp_name) /
             pl.col(aligned_span_len_name)).alias(f"MM_{idx}"),
            (pl.col(e_tp_name) /
             pl.col(aligned_span_len_name)).alias(f"M_{idx}"),
        ]
        stat_exprs.extend(exprs)

    return aggr_df.select(*stat_exprs)


def np_level_plotter(df: pl.DataFrame, filename: str):
    res_png_name = filename.rsplit(".", maxsplit=1)[0] + ".np_qv_box_plot.png"
    df = (
        df.with_columns(
            [
                pl.when(pl.col("subreadPasses").ge(20))
                .then(20)
                .otherwise(pl.col("subreadPasses"))
                .alias("newSubreadPasses")
            ]
        )
        .drop(["subreadPasses"])
        .rename({"newSubreadPasses": "subreadPasses"})
        .sort(["subreadPasses"], descending=[False])
    )

    fig = plt.figure(figsize=(20, 10))
    axes = fig.add_subplot(1, 1, 1)
    plt.sca(axes)
    sns.boxplot(df, x="subreadPasses", y="concordanceQv")
    set_tick_size(axes)
    plt.grid(visible=True, color="r", linestyle="--")
    logging.info(f"saving fig to {res_png_name}")
    fig.savefig(fname=res_png_name)


def joined_np_level_plotter(
    channel_level_dfs: List[pl.DataFrame],
    primary_np_idx: int,
    skip_subreads_df: bool = False,
):
    res_png_name = "joined_np_qv_box_plot.png"
    joined_df, not_none_indices = _join_dfs(
        channel_level_dfs=channel_level_dfs)
    joined_df = (
        joined_df.rename({f"subreadPasses_{primary_np_idx}": "subreadPasses"})
        .with_columns(
            [
                pl.when(pl.col("subreadPasses").ge(20))
                .then(20)
                .otherwise(pl.col("subreadPasses"))
                .alias("newSubreadPasses")
            ]
        )
        .drop(["subreadPasses"])
        .rename({"newSubreadPasses": "subreadPasses"})
        .sort(["subreadPasses"], descending=[False])
    )

    col2rows = []
    qv_name_fmt = "concordanceQv_{}"
    for idx in not_none_indices:
        if skip_subreads_df and idx == 0:
            continue
        col2rows.append(
            joined_df.select(
                "subreadPasses",
                pl.col(qv_name_fmt.format(idx)).alias("concordanceQv"),
                pl.lit(idx).alias("tag"),
            )
        )
    df = pl.concat(col2rows).sort(["subreadPasses"], descending=[False])
    fig = plt.figure(figsize=(20, 10))
    axes = fig.add_subplot(1, 1, 1)
    plt.sca(axes)
    sns.boxplot(df, x="subreadPasses", y="concordanceQv", hue="tag")
    set_tick_size(axes)

    plt.grid(visible=True, color="r", linestyle="--")
    logging.info(f"saving fig to {res_png_name}")
    fig.savefig(fname=res_png_name)
    pass


def _join_dfs(channel_level_dfs: List[pl.DataFrame]) -> Tuple[pl.DataFrame, List[int]]:
    """inner join multiple dfs according to their channel_id"""
    joined_df = None
    not_none_indices = []
    for idx, df in enumerate(channel_level_dfs):
        if df is None:
            continue
        not_none_indices.append(idx)
        df = df.rename(
            {
                ori_name: f"{ori_name}_{idx}"
                for ori_name in df.columns
                if ori_name != "channel_id"
            }
        )
        if joined_df is None:
            joined_df = df
            continue
        joined_df = joined_df.join(df, on="channel_id", how="inner", suffix="")
    return joined_df, not_none_indices


def main(args):

    outdir = (
        args.outdir
        if args.outdir is not None
        else "{}_metric".format(args.s_bam.rsplit(".", maxsplit=1)[0])
    )

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    bams = []
    bams.append(args.s_bam)
    bams.extend(args.smc_bams)

    basic_stat_paths = []

    for idx, bam_file in enumerate(bams):
        basic_stat_path = None
        if bam_file is not None:
            basic_stat_path = align_and_stat(
                bam_file=bam_file,
                ref_fasta=args.reffasta,
                outdir=outdir,
                rq_thr=0 if idx == 0 else args.smc_rq_thr,
                hc_regions=args.hc_regions,
                hc_variants=args.hc_variants,
            )
        basic_stat_paths.append(basic_stat_path)

    channel_level_dfs = []
    for idx, csv_file in enumerate(basic_stat_paths):
        channel_level = None
        if csv_file is not None:
            channel_level = channel_level_aggr(
                csv_file=csv_file,
                idx=idx,
                np_thr=args.np_thr,
                sbr_qv_thr=args.sbr_qv_thr,
                rq_thr=None if idx == 0 else args.smc_rq_thr,
            )
        channel_level_dfs.append(channel_level)

    logging.info("JOINED_NP_LEVEL_AGGR")
    if len(channel_level_dfs) > 1:
        joined_np_level_df = joined_np_level_aggr(
            channel_level_dfs=channel_level_dfs, primary_np_idx=args.primary_np_idx
        )
        print(joined_np_level_df)

    logging.info("JOINED_GLOBAL_AGGR")
    if len(channel_level_dfs) > 1:
        joined_global_df = joined_global_aggr(
            channel_level_dfs=channel_level_dfs)
        print(joined_global_df)
        joined_global_df = joined_global_aggr_error_info(
            channel_level_dfs=channel_level_dfs
        )
        print(joined_global_df)

    logging.info("GLOBAL_AGGR")
    for idx, csv_file in enumerate(basic_stat_paths):
        if csv_file is not None:
            df = channel_level_aggr(csv_file=csv_file, idx=idx)
            print(idx, bams[idx], "\n", global_aggr(df))


if __name__ == "__main__":
    polars_env_init()

    parser = argparse.ArgumentParser(prog="subreads_smc_funnel_stat")
    parser.add_argument("reffasta", metavar="ref.fasta", type=str)
    parser.add_argument("--s_bam", metavar="subreads.bam",
                        type=str, default=None)
    parser.add_argument(
        "--smc_bams", metavar="smc1.bam smc2.bam ...", nargs="+", type=str, default=[]
    )
    parser.add_argument("--primary_np_idx", metavar=1, type=int, default=1)
    parser.add_argument(
        "--np_thr", metavar=25, type=int, default=25, help="np = min(np, np_thr)"
    )
    parser.add_argument(
        "--sbr_qv_thr",
        metavar=0,
        type=float,
        default=0,
        help="only the channel that subreads_qv ≥ sbr_qv_thr will be reserved",
    )

    parser.add_argument(
        "--smc-rq-thr",
        metavar=0.99,
        type=float,
        default=0.99,
        help="only the channel that smc.rq ≥ smc_rq_thr will be reserved",
        dest="smc_rq_thr",
    )

    parser.add_argument("--skip_sbr_df_plot", action="store_true")
    parser.add_argument(
        "--hc-regions",
        dest="hc_regions",
        default=None,
        type=str,
        help="hc regions bed file",
    )
    parser.add_argument(
        "--hc-variants",
        dest="hc_variants",
        default=None,
        type=str,
        help="hc variants vcf file",
    )

    parser.add_argument("--outdir", default=None,
                        type=str, help="ana output dir")

    args_ = parser.parse_args()
    main(args_)
