import os
import sys

cur_dir = os.path.abspath(__file__).rsplit("/", maxsplit=1)[0]
sys.path.insert(0, cur_dir)

import pysam
import utils
import polars as pl
import argparse
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt


class BaseQStat:
    def __init__(self, pred_q):
        self.pred_q = pred_q
        self.eq = 0
        self.diff = 0
        self.insertion = 0
        self.deletion = 0
        self.depth = 0

    def add_eq(self, num=1):
        self.eq += num

    def add_diff(self, num=1):
        self.diff += num

    def add_insertion(self, num=1):
        self.insertion += num

    def add_deletion(self, num=1):
        self.deletion += num


def stat_one_record(aligned_pairs, qual, baseq2baseq_stat, query_seq, refseq, ref_end):
    ref_pos_cursor = None
    query_pos_cursor = None

    for qpos, rpos in aligned_pairs:

        if qpos is not None:
            query_pos_cursor = qpos
        if rpos is not None:
            ref_pos_cursor = rpos

        if ref_pos_cursor is None:
            continue

        if query_pos_cursor is None:
            continue

        baseq_stat = baseq2baseq_stat.setdefault(
            qual[query_pos_cursor], BaseQStat(qual[query_pos_cursor])
        )
        assert isinstance(baseq_stat, BaseQStat)
        if rpos is None:
            baseq_stat.add_insertion()
            if ref_pos_cursor == (ref_end - 1):
                break
            continue

        if qpos is None:
            baseq_stat.add_deletion()
        else:
            if refseq[rpos] == query_seq[qpos]:
                baseq_stat.add_eq()
            else:
                baseq_stat.add_diff()

        if ref_pos_cursor == (ref_end - 1):
            break


def stat(aligned_bam_file: str, ref_file: str):
    ref_data = utils.read_fastx_file(ref_file)

    baseq2baseq_stat = {}

    with pysam.AlignmentFile(aligned_bam_file, mode="rb") as bam_h:
        for refname, refseq in ref_data.items():

            for record in tqdm(
                bam_h.fetch(contig=refname), desc=f"processing {refname}"
            ):
                if record.is_secondary or record.is_supplementary or record.is_unmapped:
                    continue

                ref_pos_cursor = None
                query_pos_cursor = None

                ref_start = record.reference_start
                ref_end = record.reference_end
                query_start = record.query_alignment_start
                query_end = record.query_alignment_end

                qual = record.query_qualities
                query_seq = record.query_sequence

                stat_one_record(
                    record.get_aligned_pairs(),
                    qual,
                    baseq2baseq_stat,
                    query_seq,
                    refseq,
                    ref_end,
                )
                # for qpos, rpos in record.get_aligned_pairs():

                #     if qpos is not None:
                #         query_pos_cursor = qpos
                #     if rpos is not None:
                #         ref_pos_cursor = rpos

                #     if ref_pos_cursor is None:
                #         continue

                #     if query_pos_cursor is None:
                #         continue

                #     baseq_stat = baseq2baseq_stat.setdefault(
                #         qual[query_pos_cursor], BaseQStat(qual[query_pos_cursor])
                #     )
                #     assert isinstance(baseq_stat, BaseQStat)
                #     if rpos is None:
                #         baseq_stat.add_insertion()
                #         if ref_pos_cursor == (ref_end - 1):
                #             break
                #         continue

                #     if qpos is None:
                #         baseq_stat.add_deletion()
                #     else:
                #         if refseq[rpos] == query_seq[qpos]:
                #             baseq_stat.add_eq()
                #         else:
                #             baseq_stat.add_diff()

                #     if ref_pos_cursor == (ref_end - 1):
                #         break
    baseqs = []
    eqs = []
    diffs = []
    insertions = []
    deletions = []

    for bq, stat in baseq2baseq_stat.items():
        assert isinstance(stat, BaseQStat)
        baseqs.append(bq)
        eqs.append(stat.eq)
        diffs.append(stat.diff)
        insertions.append(stat.insertion)
        deletions.append(stat.deletion)

    return pl.DataFrame(
        {"baseq": baseqs, "eq": eqs, "diff": diffs, "ins": insertions, "del": deletions}
    )


def main(args):
    """Deprected. use pred-baseq-and-emp-q.py instead"""
    df = stat(args.aln_bam, args.ref_file)
    df = df.with_columns(
        [
            (pl.col("eq") / (pl.col("eq") + pl.col("diff") + pl.col("ins"))).alias(
                "emp_rq"
            )
        ]
    ).with_columns([utils.q2phreq_expr("emp_rq", "emp_phreq")])
    figure = plt.figure(figsize=(20, 10))
    axs = figure.add_subplot(1, 1, 1)
    sns.scatterplot(df.to_pandas(), x="baseq", y="emp_phreq", ax=axs)

    print(df.head(10))
    figure.savefig(fname="baseq2empq")


if __name__ == "__main__":

    params = {
        "aln_bam": "/data/ccs_data/ccs_eval2024q4/output-all/smc2ref.bam",
        "ref_file": "/data/ccs_data/MG1655.fa",
    }

    main(argparse.Namespace(**params))
