import os
import argparse
import logging
import polars as pl
from tqdm import tqdm
from typing import List

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def analyze_motif_prob(
    metric_dir: str,
    motifs: List[str],
    called: int,
    tag: str = "pure+mixed",
    ratio_lo: float = 0.0,
    ratio_hi: float = 1.0,
) -> dict:
    """Scan all aggr.csv files under metric_dir, filter by motif/called/tag/ratio,
    and return motif counts among matching rows."""
    motif_counts: dict[str, int] = {m: 0 for m in motifs}
    total = 0
    files_scanned = 0

    pbar = tqdm(desc="processing")

    for dirpath, _, filenames in os.walk(metric_dir):
        for fname in filenames:
            if not fname.endswith("-hp-aggr.csv"):
                continue
            pbar.update(1)
            fpath = os.path.join(dirpath, fname)

            df = (
                pl.scan_csv(fpath, separator="\t")
                .filter(
                    pl.col("motif").is_in(motifs)
                    & (pl.col("tag") == tag)
                    & (pl.col("called") == called)
                    & (pl.col("ratio_within_motif_tag") < ratio_hi)
                    & (pl.col("ratio_within_motif_tag") > ratio_lo)
                )
                .collect()
            )
            total += len(df)
            for motif in motifs:
                motif_counts[motif] += int(df.filter(pl.col("motif")
                                           == motif).height)
            files_scanned += 1
    pbar.close()

    print(f"files_scanned: {files_scanned}")
    print(f"rows matching: {total}")
    for m in motifs:
        prob = motif_counts[m] / total if total > 0 else 0.0
        print(f"P(motif={m} | ...) = {motif_counts[m]} / {total} = {prob:.6f}")

    return {"motif_counts": motif_counts, "total": total, "files_scanned": files_scanned}


def analyze_c4_c5_called5_prob(metric_dir: str) -> dict:
    """Convenience wrapper: P(motif=C(4) or C(5) | called==5, ratio_within_motif_tag<0.3)."""
    result = analyze_motif_prob(
        metric_dir,
        motifs=["(C)5", "(C)4"],
        called=5,
        tag="pure+mixed",
        ratio_lo=0.2,
        ratio_hi=0.3,
    )
    c4 = result["motif_counts"].get("(C)4", 0)
    c5 = result["motif_counts"].get("(C)5", 0)
    total = result["total"]
    c4_prob = c4 / total if total > 0 else 0.0
    c5_prob = c5 / total if total > 0 else 0.0

    print(f"P(motif=C(4) | ...) = {c4} / {total} = {c4_prob:.6f}")
    print(f"P(motif=C(5) | ...) = {c5} / {total} = {c5_prob:.6f}")

    return {"c4_count": c4, "c5_count": c5, "total": total,
            "c4_prob": c4_prob, "c5_prob": c5_prob}


def main_cli():
    """
    接着 channel_subreads_ref_analysis.py 的结果继续进行分析
    channel_subreads_ref_analysis.py 会生成一个 channel_subreads 文件夹，这个工具的输入就是这个文件夹
    """
    parser = argparse.ArgumentParser(prog="channel_subreads_ref_analysis2")
    parser.add_argument(
        "metric_dir", nargs="?", default=None,
        help="channel_subreads directory (default: auto-detect)")
    parser.add_argument("--motifs", nargs="+", default=None,
                        help="motif values to analyze (default: (C)5 (C)4)")
    parser.add_argument("--called", type=int, default=None,
                        help="called value to filter on (default: 5)")
    parser.add_argument("--tag", default="pure+mixed",
                        help="tag filter (default: pure+mixed)")
    parser.add_argument("--ratio-lo", type=float, default=0.2,
                        help="lower bound for ratio_within_motif_tag (default: 0.2)")
    parser.add_argument("--ratio-hi", type=float, default=0.3,
                        help="upper bound for ratio_within_motif_tag (default: 0.3)")
    parser.add_argument("--general", action="store_true",
                        help="use the general analyze_motif_prob function")
    args = parser.parse_args()

    metric_dir = args.metric_dir or "/data1/ccs_data/202603-rna-data/rna_data/RNA0.5K/channel_subreads"
    print(f"motifs:{args.motifs}")

    if args.general:
        motifs = args.motifs or ["(C)5", "(C)4"]
        analyze_motif_prob(
            metric_dir,
            motifs=motifs,
            called=args.called or 5,
            tag=args.tag,
            ratio_lo=args.ratio_lo,
            ratio_hi=args.ratio_hi,
        )
    else:
        analyze_c4_c5_called5_prob(metric_dir)


if __name__ == "__main__":
    main_cli()
