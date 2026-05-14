from tqdm import tqdm
import pysam
import os
import argparse
import pathlib
import logging
import sys
from itertools import islice
import multiprocessing as mp
from multiprocessing import Pool
sys.path.append(os.path.abspath(__file__).rsplit("/", maxsplit=1)[0])
import reads_quality_stats_hp  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def _parse_range(range_str: str):
    parts = range_str.split(":")
    return float(parts[0]), float(parts[1])


def _extract_partition(args) -> dict:
    """Extract a partition of channels from subreads.bam. Returns {ch: count, ...}."""
    partition_channels, subreads_bam, outdir = args
    valid_channels = set(partition_channels)
    valid_bams = []
    with pysam.AlignmentFile(subreads_bam, "rb", check_sq=False, threads=1) as in_bam:
        header = in_bam.header
        out_handle = None
        current_ch = None

        for record in in_bam.fetch(until_eof=True):
            ch = int(record.get_tag("ch"))

            if ch not in valid_channels:
                continue

            if ch != current_ch:
                # Flush previous channel
                if out_handle is not None:
                    out_handle.close()
                out_path = os.path.join(
                    outdir, f"{pathlib.Path(subreads_bam).stem}.ch{ch}.filtered.bam")
                valid_bams.append(out_path)
                out_handle = pysam.AlignmentFile(
                    out_path, "wb", header=header, check_sq=False, threads=1)

                current_ch = ch

            out_handle.write(record)

        if out_handle is not None:
            out_handle.close()
    return valid_bams


def _process_partition(partition_bams: list, ref_fa: str) -> list:
    """Process a partition of channel BAMs with reads_quality_stats_hp."""
    results = []
    for ch_bam in partition_bams:
        logging.info(f"processing {ch_bam} with ref_fa={ref_fa}")
        aggr_csv, fact_csv = reads_quality_stats_hp.main(
            bam_file=ch_bam, ref_fa=ref_fa, force=True, threads=10)
        results.append((ch_bam, (aggr_csv, fact_csv)))
        logging.info(f"  aggr={aggr_csv}, fact={fact_csv}")
    return results


def chunked(iterable, batch_size, *args):
    it = iter(iterable)

    while True:
        batch = list(islice(it, batch_size))

        if not batch:
            break

        yield (batch, *args)


def dump_channel_bams(
    subreads_bam: str,
    smc_bam: str,
    np_range: str,
    rq_range: str,
    outdir: str,
    ref_fa: str,
) -> list:
    np_min, np_max = _parse_range(np_range)
    rq_min, rq_max = _parse_range(rq_range)

    # Step 1: scan smc.bam, collect channel stats (ch -> np, rq)
    # SMC BAM: reads from same channel are contiguous, each read has ch/np/rq tags
    tot = 0
    valid_channels = set()
    with pysam.AlignmentFile(smc_bam, "rb", check_sq=False, threads=os.cpu_count()) as bam:
        for record in tqdm(bam.fetch(until_eof=True), desc="scanning smc.bam"):
            ch = int(record.get_tag("ch"))
            np = int(record.get_tag("np"))
            rq = float(record.get_tag("rq"))
            tot += 1
            if np_min <= np <= np_max and rq_min <= rq <= rq_max:
                valid_channels.add(ch)

    logging.info(
        f"{len(valid_channels)} channels pass filters: np_range={np_range}, rq_range={rq_range}")

    if not valid_channels:
        logging.warning("No valid channels found, nothing to do")
        return []

    # Create output dir
    os.makedirs(outdir, exist_ok=True)

    # Step 3: streaming extraction from subreads.bam
    # Reads from the same channel are contiguous, so we open/close output
    # handles as channel boundaries are encountered instead of keeping all open.
    valid_bams = []

    with Pool(processes=os.cpu_count() // 2) as pool:
        valid_bams = []
        for result in tqdm(pool.imap(_extract_partition, chunked(valid_channels, 100, subreads_bam, outdir), chunksize=1), desc=f"extrating channel subreads"):
            valid_bams.extend(result)

    print(valid_bams)

    # Process channel BAMs with multiprocessing, each process handles a partition
    n_workers = os.cpu_count() // 10
    n_channels = len(valid_bams)
    partition_size = max(1, n_channels // n_workers)
    partitions = [
        valid_bams[i:i + partition_size]
        for i in range(0, len(valid_bams), partition_size)
    ]
    logging.info(
        f"processing {n_channels} channels with {len(partitions)} partitions")
    with Pool(processes=n_workers) as pool:
        all_results = pool.starmap(
            _process_partition, [(p, ref_fa) for p in partitions]
        )
    results = [r for sublist in all_results for r in sublist]

    return results


def main_cli():
    parser = argparse.ArgumentParser(prog="channel_subreads_ref_analysis")
    parser.add_argument("subreads_bam", help="subreads BAM file")
    parser.add_argument("smc_bam", help="SMC BAM file (has ch/np/rq tags)")
    parser.add_argument(
        "ref_fa", help="reference FASTA (reserved for future use)")
    parser.add_argument("--np-range", default="5:100000000", dest="np_range",
                        help="np range start:end")
    parser.add_argument("--rq-range", default="0.95:1.1", dest="rq_range",
                        help="rq range start:end")
    parser.add_argument("--outdir", default=None, dest="out_dir",
                        help="output directory")
    args = parser.parse_args()

    if args.out_dir is None:
        out_dir = os.path.join(os.path.dirname(
            args.subreads_bam) or ".", "channel_subreads")
    else:
        out_dir = args.out_dir

    results = dump_channel_bams(
        subreads_bam=args.subreads_bam,
        smc_bam=args.smc_bam,
        np_range=args.np_range,
        rq_range=args.rq_range,
        outdir=out_dir,
        ref_fa=args.ref_fa,
    )
    logging.info(f"done: {len(results)} channels processed")


if __name__ == "__main__":
    """
        /data1/ccs_data/202603-rna-data/rna_data/RNA0.5K

        /data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/20260327_240601Y0004_Run0001_demuxed.bam /data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/20260327_240601Y0004_Run0001_demuxed.smc_all_reads.bam /data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/ref_0.5K-nopoly.fa --np-range 3:100 --rq-range 0.99:1.1
        
        
        
        
        python channel_subreads_ref_analysis.py /data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/20260327_240601Y0004_Run0001_demuxed.bam  /data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/20260327_240601Y0004_Run0001_demuxed.smc_all_reads.bam /data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/ref_0.5K-nopoly.fa --np-range 5:100 --rq-range 0.99:1.1
    """
    main_cli()
    # for a in chunked([1, 2, 3, 4, 5, 6, 7], 2, "g", "c"):
    # print(a)
