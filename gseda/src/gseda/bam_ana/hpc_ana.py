import pathlib
import os
import subprocess
import sys
import pysam
from tqdm import tqdm
import polars as pl
import argparse

pre_dir = os.path.abspath(__file__).rsplit("/", maxsplit=2)[0]
print(pre_dir)
sys.path.append(pre_dir)

from fact_table_ana.utils import read_fastx_file
from fact_table_ana.polars_init import polars_env_init


def do_alignment(hpc_bam: str, hpc_ref: str):
    out_root = os.path.dirname(hpc_bam)

    o_path = "{}/{}-{}.aligned".format(
        out_root, pathlib.Path(hpc_bam).stem, pathlib.Path(hpc_ref).stem
    )

    cmd = f"gsmm2 align -q {hpc_bam} -t {hpc_ref} --noMar -p {o_path}"
    subprocess.check_call(cmd, shell=True)
    return f"{o_path}.bam"


def read_bam_data(aligned_bam: str, hpc_ref_cnt: str, hpc_ref_seq: str):
    ref_hpc_seq_data = read_fastx_file(hpc_ref_seq)

    ref_hpc_cnt_data = read_fastx_file(hpc_ref_cnt)
    ref_hpc_cnt_data = {
        qname: list(map(int, seq[0].split(",")))
        for qname, seq in ref_hpc_cnt_data.items()
    }

    counter = {}

    with pysam.AlignmentFile(aligned_bam, mode="rb", threads=40) as bam_h:
        for refname in ref_hpc_cnt_data.keys():
            ref_seq = ref_hpc_seq_data[refname][0]
            ref_hpc_cnt = ref_hpc_cnt_data[refname]

            for i, record in enumerate(
                tqdm(bam_h.fetch(refname), desc=f"reading {refname}")
            ):
                if i > 50000:
                    break
                aligned_pairs = record.get_aligned_pairs(matches_only=True)
                query_hpc_cnt = record.get_tag("dw")
                query_seq = record.query_sequence
                for qpos, rpos in aligned_pairs:
                    if query_seq[qpos] != ref_seq[rpos]:
                        continue

                    if ref_hpc_cnt[rpos] <= 1:
                        continue
                    counter.setdefault(
                        f"{ref_seq[rpos]}{ref_hpc_cnt[rpos]}", []
                    ).append(query_hpc_cnt[qpos])
    keys = []
    values = []
    for k, counts in counter.items():
        for c in counts:
            keys.append(k)
            values.append(c)

    return pl.DataFrame({"key": keys, "value": values})


def ana(df: pl.DataFrame):
    df = (
        df.group_by(["key", "value"])
        .agg([pl.len().alias("cnt")])
        .with_columns(
            [
                (pl.col("cnt") / pl.col("cnt").sum().over(pl.col("key"))).alias(
                    "ratio"
                ),
                pl.col("key").str.slice(1).cast(pl.Int32).alias("ref_hp_cnt"),
            ]
        )
        .filter((pl.col("ref_hp_cnt") - pl.col("value")).abs() <= 1)
    ).sort(by=["key", "value"])

    df = df.with_columns(
        [
            (pl.col("ref_hp_cnt") + 1).alias("ref_hp_cnt_ADD_1"),
            (pl.col("ref_hp_cnt") - 1).alias("ref_hp_cnt_MINUS_1"),
        ]
    )

    df_eq = df.filter(pl.col("ref_hp_cnt") == pl.col("value"))
    df_minus_one = df.filter(pl.col("ref_hp_cnt_MINUS_1") == pl.col("value"))
    df_add_one = df.filter(pl.col("ref_hp_cnt_ADD_1") == pl.col("value"))

    all = df_eq.join(df_minus_one, on="key", suffix="_minus_one").join(
        df_add_one, on="key", suffix="_add_one"
    )

    print(all.head(1000))


def main(args):

    # align_res = do_alignment(args.hpc_bam, args.hpc_ref)
    align_res = "/data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.hpc-ref_Saureus_ATCC25923.m.new.corrected.hpc.aligned.bam"
    df = read_bam_data(align_res, args.hpc_ref_cnt, args.hpc_ref)
    ana(df=df)

    pass


if __name__ == "__main__":
    polars_env_init()

    cli_params = {
        "hpc_bam": "/data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.subreads.hpc.bam",
        "hpc_ref": "/data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.hpc.fa",
        "hpc_ref_cnt": "/data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.hpc.txt",
    }

    main(argparse.Namespace(**cli_params))
