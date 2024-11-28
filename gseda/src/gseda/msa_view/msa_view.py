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
from typing import Dict


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


def init_matrix(num_rows, num_cols):
    """init matrix using "." """
    matrix = np.empty(shape=[num_rows, num_cols], dtype=np.str_)
    matrix.fill(".")
    return matrix


class ResultMatrix:
    """msa alignment matrix"""

    def __init__(
        self,
        ref_start,
        ref_end,
        query_names,
        refpos2length: Dict[int, int],
        ref_name=None,
    ) -> None:
        self.num_records = len(query_names) + 1  # +ref
        self.query_names = sorted(query_names)
        self.query2idx = {
            query_name: idx for idx, query_name in enumerate(self.query_names, start=1)
        }
        self.ref_name = "REF" if ref_name is None else ref_name
        self.ref_name = f"{self.ref_name}_{ref_start}_{ref_end}"
        self.matrix = init_matrix(self.num_records, sum(refpos2length.values()))

        refpos2length_list = sorted(list(refpos2length.items()), key=lambda x: x[0])
        rpos2matrix_col = [[refpos2length_list[0][0], 0]]
        for i in range(1, len(refpos2length_list)):
            cur_item = refpos2length_list[i]
            rpos2matrix_col.append(
                [cur_item[0], refpos2length_list[i - 1][1] + rpos2matrix_col[-1][1]]
            )

        print(refpos2length_list)
        print(rpos2matrix_col)
        self.rpos2matrix_col = {
            rpos: matrix_col for rpos, matrix_col in rpos2matrix_col
        }

        self.ref_end = ref_end  # exclusive
        self.ref_start = ref_start

        print(f"self.end={self.ref_end}")

    def update(self, record: pysam.AlignedSegment, ref: str = None):

        idx = self.query2idx[build_query_name(record)]

        rpos_cursor = None
        # qpos_cursor = None
        offset = 0
        query_seq = record.query_sequence

        ref_aligned = []
        query_aligned = []

        qpos_start = None
        qpos_end = None

        for qpos, rpos in record.get_aligned_pairs():
            if rpos is not None:
                rpos_cursor = rpos
            if rpos_cursor is None:
                continue
            if rpos_cursor < self.ref_start or rpos_cursor >= self.ref_end:
                continue

            if qpos_start is None and qpos is not None:
                qpos_start = qpos

            if qpos is not None:
                qpos_end = qpos

            if rpos_cursor not in self.rpos2matrix_col:
                print(
                    rpos_cursor, " not in ", sorted(list(self.rpos2matrix_col.keys()))
                )
                raise ValueError()

            ref_aligned.append("-" if rpos is None else ref[rpos])
            query_aligned.append("-" if qpos is None else query_seq[qpos])

            matrix_init_col = self.rpos2matrix_col[rpos_cursor]
            if rpos is None:
                offset += 1
            else:
                offset = 0
                self.matrix[0, matrix_init_col] = ref[rpos]

            matrix_col = matrix_init_col + offset
            if qpos is not None:
                self.matrix[idx, matrix_col] = query_seq[qpos]

        # ref_aligned = "".join(ref_aligned)
        # query_aligned = "".join(query_aligned)
        # info = f"qname:{record.query_name}\n{ref_aligned}\n{query_aligned}"
        # print(info)

        seq_len = record.query_length
        if record.is_reverse:
            qpos_start, qpos_end = seq_len - qpos_end, seq_len - qpos_start

        called_start = None
        called_end = None
        if record.has_tag("be"):
            shift = record.get_tag("be")[0]
            called_start = qpos_start + shift
            called_end = qpos_end + shift

        print(
            f"{record.query_name}: sbr:{qpos_start}-{qpos_end}, called:{called_start}-{called_end}"
        )

    def get_raw_result(self):
        return self.matrix

    def get_query_names(self):
        names = [self.ref_name]
        names.extend(self.query_names)
        return names

    def get_result(self):
        """may the matrix has invalid rows, trim it and return"""
        return self.get_raw_result()

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


def build_ref_pos_maxins(
    sam_file: pysam.AlignmentFile, contig: str, ref_start: int, ref_end: int
):
    rpos2max_ins = {rpos: 0 for rpos in range(ref_start, ref_end)}

    for query in sam_file.fetch(contig=contig):
        rpos_cursor = None
        cur_query_ins = 0
        for _, rpos in query.get_aligned_pairs():

            if rpos is not None:
                rpos_cursor = rpos

            if rpos_cursor is None:
                continue
            if rpos_cursor < ref_start:
                continue

            if rpos_cursor >= ref_end:
                break

            if rpos is None:
                cur_query_ins += 1
            else:
                if rpos_cursor > ref_start:
                    rpos2max_ins[rpos_cursor - 1] = max(
                        cur_query_ins, rpos2max_ins[rpos_cursor - 1]
                    )
                cur_query_ins = 0

        if rpos_cursor > ref_start and rpos_cursor <= ref_end:
            rpos2max_ins[rpos_cursor - 1] = max(
                cur_query_ins, rpos2max_ins[rpos_cursor - 1]
            )

    return rpos2max_ins


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

    print(
        f"interested_ref_start={interested_ref_start}, interested_ref_end={interested_ref_end}"
    )

    rpos2maxins = build_ref_pos_maxins(
        samfile,
        contig=ref_name,
        ref_start=interested_ref_start,
        ref_end=interested_ref_end,
    )

    rpos2length = {pos: ins + 1 for pos, ins in rpos2maxins.items()}

    result_matrix = ResultMatrix(
        ref_start=interested_ref_start,
        ref_end=interested_ref_end,
        query_names=query_names,
        refpos2length=rpos2length,
        ref_name=ref_name,
    )

    for query in samfile.fetch(contig=ref_name):
        if query.reference_end <= interested_ref_start:
            continue
        if query.reference_start >= interested_ref_end:
            continue
        result_matrix.update(query, ref=ref_seq)

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


if __name__ == "__main__":
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
