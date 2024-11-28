import pysam
from typing import Mapping, Tuple
from tqdm import tqdm
import polars as pl


def read_bam_file(bam_file: str) -> Mapping[str, Tuple[str, int]]:
    res = {}
    with pysam.AlignmentFile(bam_file, mode="rb", threads=40, check_sq=False) as bam_h:
        for record in tqdm(
            bam_h.fetch(until_eof=True), desc=f"read_bam_file:>> reading {bam_file}"
        ):
            res[record.query_name] = (record.query_sequence, f"00_{record.query_name}")

    return res


def read_fastx_file(fname: str) -> Mapping[str, Tuple[str, int]]:
    fh = pysam.FastxFile(fname)
    res = {}
    for entry in fh:
        res[entry.name] = (entry.sequence, f"00_{entry.name}")
    return res


def q2phreq_expr(inp_name, oup_name=None):
    oup_name = oup_name if oup_name is not None else inp_name
    return (
        -10.0
        * (
            1
            - pl.when(pl.col(inp_name) > (1 - 1e-6))
            .then(1 - 1e-6)
            .otherwise(pl.col(inp_name))
        ).log10()
    ).alias(oup_name)
