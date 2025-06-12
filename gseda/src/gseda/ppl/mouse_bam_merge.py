"""
将多个 adapter.bam 和 smc_all_reads.bam 文件合并成一个. dw 和 ar 特征需要从 called bam 中进行透传
"""

from typing import List, Mapping
import pysam
from tqdm import tqdm
import pathlib


def merge_smc_bams(bams: List[str], o_bam_name: str) -> Mapping[str, str]:
    old2new = {}
    out_bam = pysam.AlignmentFile(
        o_bam_name,
        mode="wb",
        check_sq=False,
        threads=40,
        header={
            "HD": {"VN": "1.5", "SO": "unknown"},
        },
    )
    for bam_path in bams:
        with pysam.AlignmentFile(
            bam_path, mode="rb", check_sq=False, threads=40
        ) as in_bam:
            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"dumping {bam_path}"
            ):
                rq = record.get_tag("rq")
                if rq < 0.968:
                    continue

                qname = record.query_name.rsplit("/", maxsplit=1)[0]
                new_qname = len(old2new)
                old2new[qname] = f"{new_qname}"
                record.query_name = old2new[qname]
                record.set_tag("ch", new_qname)
                out_bam.write(record)

    if out_bam is not None:
        out_bam.close()
    return old2new


class Channel:
    def __init__(self, dw, ar):
        self.dw = dw
        self.ar = ar


def reset_dw_ar_of_adapter_bam(called_bam: str, adapter_bam: str, new_adapter_bam: str):
    ch2dw_ar = {}
    with pysam.AlignmentFile(
        called_bam, mode="rb", check_sq=False, threads=40
    ) as in_bam:
        for record in tqdm(in_bam.fetch(until_eof=True), desc=f"reading {called_bam}"):
            ch = int(record.query_name.rsplit(sep="_", maxsplit=1)[1])
            dw = record.get_tag("dw")
            ar = record.get_tag("ar")
            ch2dw_ar[ch] = Channel(dw=dw, ar=ar)

    with pysam.AlignmentFile(
        adapter_bam, mode="rb", check_sq=False, threads=40
    ) as in_bam:
        with pysam.AlignmentFile(
            new_adapter_bam, mode="wb", check_sq=False, threads=40, header=in_bam.header
        ) as out_bam:
            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"dumping {adapter_bam}"
            ):
                ch = int(record.get_tag("ch"))
                dw_ar = ch2dw_ar[ch]
                start, end = record.get_tag("be")
                record.set_tag("dw", dw_ar.dw[start:end])
                record.set_tag("ar", dw_ar.ar[start:end])
                out_bam.write(record)


def dump_merged_adapter_bam(
    adapter_bam_paths: List[str], old2new: Mapping[str, str], out_bam_path: str
):
    out_bam = pysam.AlignmentFile(
        out_bam_path,
        mode="wb",
        check_sq=False,
        threads=40,
        header={
            "HD": {"VN": "1.5", "SO": "unknown"},
        },
    )

    for adapter_bam_path in adapter_bam_paths:
        with pysam.AlignmentFile(
            adapter_bam_path, mode="rb", check_sq=False, threads=40
        ) as in_bam:
            run_name = in_bam.header["RG"][0]["rn"]
            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"dumping {adapter_bam_path}"
            ):
                sbr_idx = record.query_name.rsplit("/", maxsplit=1)[-1]
                ch = record.get_tag("ch")
                key = f"{run_name}/{ch}"
                if key not in old2new:
                    continue
                new_ch = old2new[key]
                record.query_name = f"{new_ch}/{sbr_idx}"
                record.set_tag("ch", int(new_ch))
                out_bam.write(record)

    out_bam.close()


def dump_merged_called_bam(
    bam_paths: List[str], old2new: Mapping[str, str], out_bam_path: str
):
    out_bam = pysam.AlignmentFile(
        out_bam_path,
        mode="wb",
        check_sq=False,
        threads=40,
        header={
            "HD": {"VN": "1.5", "SO": "unknown"},
        },
    )

    for bam_path in bam_paths:
        with pysam.AlignmentFile(
            bam_path, mode="rb", check_sq=False, threads=40
        ) as in_bam:
            run_name = in_bam.header["RG"][0]["rn"]
            for record in tqdm(
                in_bam.fetch(until_eof=True), desc=f"dumping {bam_path}"
            ):
                ch = record.query_name.split("_", maxsplit=1)[1]
                key = f"{run_name}/{ch}"
                if key not in old2new:
                    continue
                new_ch = old2new[key]
                record.query_name = f"read_{new_ch}"
                record.set_tag("ch", int(new_ch))
                out_bam.write(record)

    out_bam.close()
    pass


def main():
    called_bams = [
        "/data/ccs_data/wga/20250605_250302Y0004_Run0001_called.bam",
        "/data/ccs_data/wga/20250606_250302Y0002_Run0004_called.bam",
    ]

    sbr_bams = [
        "/data/ccs_data/wga/20250605_250302Y0004_Run0001_adapter.bam",
        "/data/ccs_data/wga/20250606_250302Y0002_Run0004_adapter.bam",
    ]

    smc_bams = [
        "/data/ccs_data/wga/20250605_250302Y0004_Run0001.smc_all_reads.bam",
        "/data/ccs_data/wga/20250606_250302Y0002_Run0004.smc_all_reads.bam",
    ]
    old2new = merge_smc_bams(smc_bams, "/data/ccs_data/wga/wga-smc.bam")

    print(list(old2new)[:10])

    sbr_new_bams = []
    # for called_bam_path, sbr_bam_path in zip(called_bams, sbr_bams):
    #     sbr_p = pathlib.Path(sbr_bam_path)
    #     new_sbr_bam_path = "{}/{}.new.bam".format(sbr_p.parent, sbr_p.stem)
    #     reset_dw_ar_of_adapter_bam(called_bam_path, sbr_bam_path, new_sbr_bam_path)
    #     sbr_new_bams.append(new_sbr_bam_path)

    # dump_merged_adapter_bam(
    #     sbr_new_bams, old2new, "/data/ccs_data/wga//wga-adapter.bam"
    # )

    dump_merged_called_bam(
        called_bams, old2new, "/data/ccs_data/wga/wga-called.bam"
    )


if __name__ == "__main__":
    main()
