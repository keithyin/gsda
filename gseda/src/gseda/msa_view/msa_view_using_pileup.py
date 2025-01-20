"""
this script is used to visualize the msa.
the result can be used in two ways:
    * copy the fasta info and past to the jalview, jalview is used for the plot. OR
    * generate the msa picture directly

the whole pipeline is:
1) do alignment (pairwise alignment)  mini_align -P -m -r ref_filepath.fa -i inp.fa -t 1 -p calls2draft
2) use this script to generate msa pic

"""

import argparse
import pysam
import numpy as np
import tempfile
import os
from typing import Mapping
from tqdm import tqdm


def read_fastx_file(fname: str):
    fh = pysam.FastxFile(fname)
    res = {}
    for entry in fh:
        res[entry.name] = entry.sequence
    return res


def read_bam_file(bam_file: str) -> Mapping[str, str]:
    res = {}
    with pysam.AlignmentFile(bam_file, mode="rb", threads=40, check_sq=False) as bam_h:
        for record in tqdm(
            bam_h.fetch(until_eof=True), desc=f"read_bam_file:>> reading {bam_file}"
        ):
            res[record.query_name] = record.query_sequence

    return res


class Name2Seq:
    def __init__(self, fname: str):
        if (
            fname.endswith("fa")
            or fname.endswith("fasta")
            or fname.endswith("fna")
            or fname.endswith("fq")
            or fname.endswith("fastq")
        ):
            self.qname2seq = read_fastx_file(fname)
        elif fname.endswith("bam"):
            self.qname2seq = read_bam_file(fname)
        else:
            raise ValueError(f"invalid file format, {fname}")

    def fetch(self, name):
        return self.qname2seq[name]


def build_query_name(align_seg: pysam.AlignedSegment):
    return f"{align_seg.query_name}_SE_{align_seg.query_alignment_start}_{align_seg.query_alignment_end}"


class ResultMatrix:
    """msa alignment matrix"""

    def __init__(self, ref_start, ref_end, query_names, ref_name=None) -> None:
        ref_len = ref_end - ref_start
        self.num_records = len(query_names) + 1  # +ref
        self.query_names = sorted(query_names)
        self.query2idx = {
            query_name: idx for idx, query_name in enumerate(self.query_names, start=1)
        }
        self.ref_name = "REF" if ref_name is None else ref_name
        self.ref_name = f"{self.ref_name}_{ref_start}_{ref_end}"
        self.matrix = ResultMatrix.init_matrix(ref_len * 2, self.num_records)

        self.ref_end = ref_end  # exclusive
        self.ref_start = ref_start

        self.cur_row = 0
        self.cur_ref_pos = self.ref_start
        self.first_query_is_rev = None

    def update(self, pileup_col: pysam.PileupColumn, ref: str = None):
        """
        query1 query2 query3 , ....

        A        A
        .        A

        """
        self.cur_ref_pos = pileup_col.reference_pos

        max_ins = 0
        position_names = []
        for query in pileup_col.pileups:
            position_names.append(query.alignment.query_name)
            if query.indel > 0:
                max_ins = max([max_ins, query.indel])
        # print(self.cur_ref_pos, "\n", "\n".join(sorted(position_names)))

        self.extend_matrix(max_offset=max_ins + 1)
        if ref is not None:
            self.matrix[self.cur_row, 0] = ref[pileup_col.reference_pos]

        # print(self.cur_ref_pos, len(pileup_col.pileups))
        for query in pileup_col.pileups:
            if query.is_refskip:
                continue
            cur_query_max_offset = max([0, query.indel])
            q_name = build_query_name(query.alignment)
            query_pos = self.query2idx[q_name]
            if query_pos == 1:
                self.first_query_is_rev = query.alignment.is_reverse
                # print(f"{query.alignment.query_name} --> {query.alignment.is_reverse}")

            qpos = query.query_position_or_next

            if query.is_del:
                qpos -= 1

            # print(pileup_col.reference_pos)
            # if q_name == "read_99766/99766/subread/0_SE_97_1160":
            #     print(
            #         pileup_col.reference_pos,
            #         qpos,
            #         query.alignment.query_sequence[
            #             qpos : (qpos + cur_query_max_offset + 1)
            #         ],
            #         query.alignment.is_reverse,
            #     )

            for offset in range(cur_query_max_offset + 1):
                if query.is_del and offset == 0:
                    continue
                self.matrix[self.cur_row + offset, query_pos] = (
                    query.alignment.query_sequence[qpos + offset]
                )

        self.cur_row += max_ins + 1

    def extend_matrix(self, max_offset):
        """if the matrix row have run out, extend the matrix"""
        if (self.cur_row + max_offset) >= self.matrix.shape[0]:
            extended_row = 2 * (self.ref_end - self.cur_ref_pos) + max_offset
            self.matrix = np.concatenate(
                [self.matrix, ResultMatrix.init_matrix(extended_row, self.num_records)],
                axis=0,
            )

    def get_raw_result(self):
        return self.matrix[: self.cur_row]

    def get_query_names(self):
        names = [self.ref_name]
        names.extend(self.query_names)
        return names

    def get_result(self):
        """may the matrix has invalid rows, trim it and return"""
        return self.get_raw_result().transpose()

    def get_result_str(self):
        """valid matrix to string"""
        res = self.get_result()
        names = [self.ref_name]
        names.extend(self.query_names)
        result_strs = []
        for row_idx in range(res.shape[0]):

            q_name = names[row_idx]
            result_strs.append(f">{q_name}")
            result_strs.append("".join(res[row_idx].tolist()))

        return "\n".join(result_strs)

    @staticmethod
    def init_matrix(num_rows, num_cols):
        """init matrix using "." """
        matrix = np.empty(shape=[num_rows, num_cols], dtype=np.str_)
        matrix.fill(".")
        return matrix


def extract_reference(ref_filename):
    """extract reference from fasta file"""
    with open(ref_filename, mode="r", encoding="utf8") as file:
        lines = file.readlines()
        ref_name = lines[0].split(" ")[0][1:]
        return (ref_name.strip(), "".join(lines[1:]))


def bam2fa4jalview(
    aligned_bam_filename,
    ref_filename,
    ref_name,
    interested_ref_start=None,
    interested_ref_end=None,
):
    """generate the fasta info that can be used in jalview from bam file and [ref file]
    Params:
        aligned_bam_filename
        contig
        ref_filename
    """
    samfile = pysam.AlignmentFile(aligned_bam_filename, mode="rb", threads=40)

    fastx_data = Name2Seq(ref_filename)

    ref_seq = fastx_data.fetch(ref_name)

    ref_start = 2**32
    ref_end = 0
    query_names = []
    for query in samfile.fetch(contig=ref_name):
        if (
            interested_ref_start is not None
            and query.reference_end <= interested_ref_start
        ):
            continue
        if (
            interested_ref_end is not None
            and query.reference_start >= interested_ref_end
        ):
            continue
        ref_start = min([query.reference_start, ref_start])
        ref_end = max([query.reference_end, ref_end])

        query_names.append(build_query_name(query))

    interested_ref_start = (
        ref_start if interested_ref_start is None else interested_ref_start
    )
    interested_ref_end = ref_end if interested_ref_end is None else interested_ref_end

    result_matrix = ResultMatrix(
        ref_start=interested_ref_start,
        ref_end=interested_ref_end,
        query_names=query_names,
        ref_name=ref_name,
    )

    for pileup_col in samfile.pileup(contig=ref_name, min_base_quality=1):
        if (
            pileup_col.reference_pos < interested_ref_start
            or pileup_col.reference_pos >= interested_ref_end
        ):
            continue
        result_matrix.update(pileup_col=pileup_col, ref=ref_seq)

    return result_matrix


def plot_msa_align(inp_filename, oup_filename=None):
    """plot the msa align according to the fasta file"""

    from pymsaviz import MsaViz

    mv = MsaViz(
        inp_filename,
        wrap_length=150,
        show_count=True,
        show_grid=True,
        color_scheme="Identity",
    )
    if oup_filename is None:
        oup_filename = f"{inp_filename}.png"
    mv.savefig(oup_filename)


def main(args):
    res = bam2fa4jalview(
        args.bam,
        ref_filename=args.ref_fasta,
        ref_name=args.ref_name,
        interested_ref_start=args.start,
        interested_ref_end=args.end,
    )

    res_str = res.get_result_str()
    if args.o_fasta is not None:
        with open(args.o_fasta, "w", encoding="utf8") as file:
            file.write(res_str)
    else:
        print(res_str)

    if args.o_pic is not None:
        if args.o_fasta is not None:
            plot_msa_align(inp_filename=args.o_fasta, oup_filename=args.o_pic)

        else:
            with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
                tmp.write(res_str)
                tmp.close()
                print(f"temp file name: {tmp.name}")
                plot_msa_align(inp_filename=tmp.name, oup_filename=args.o_pic)
            os.remove(tmp.name)
    return res.first_query_is_rev


def main_cli():
    p = argparse.ArgumentParser(
        "mas_view",
        description="""
        the whole pipeline is:
            1) do alignment (pairwise alignment), minimap2 is an option;
            2) use this script to generate msa pic
    """,
    )
    p.add_argument("--bam", type=str, help="subreads2smc alignment", required=True)
    p.add_argument(
        "--ref-fastx-or-bam",
        type=str,
        help="smc.fa/smc.fq/.bam,",
        required=True,
        dest="ref_fasta",
    )

    p.add_argument("--ref-name", required=True, type=str, dest="ref_name")
    p.add_argument(
        "--ref-start", type=int, default=None, help="contig end", dest="start"
    )
    p.add_argument("--ref-end", type=int, default=None, help="contig start", dest="end")

    p.add_argument(
        "--o-fasta",
        default=None,
        help="output fasta file, if not provided, this content will be output to the stdout",
        dest="o_fasta",
    )
    p.add_argument(
        "--o-pic",
        default=None,
        help="visualization picture file path, if not provided, the plot procedure will be skipped",
        dest="o_pic",
    )
    args = p.parse_args()
    main(args)


if __name__ == "__main__":
    main_cli()
