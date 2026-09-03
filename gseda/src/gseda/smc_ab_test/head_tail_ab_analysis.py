#!/usr/bin/env python3
"""Compare SMC strategy changes by measuring edit distance on head/tail windows of aligned reads.

Two input modes:
  file mode:  --exp-fq <file> --ctrl-fq <file> — match reads by qname between the two files
  folder mode: --exp-dirs <d1,d2,...> --ctrl-dirs <d1,d2,...> — glob *.fastq/*.fq in each dir,
               correlate by stem name, then match by qname

For each matching read pair (exp_seq, ctrl_seq):
  - head_win = seq[:window]
  - tail_win = seq[-window:]
  - edit_dist(exp_head, ctrl_head), edit_dist(exp_tail, ctrl_tail)

Output:
  - JSON stats file with per-pair and aggregate statistics
  - PNG histograms of edit distance distributions
"""

import argparse
import json
import logging
import multiprocessing as mp
import os
import pathlib
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Tuple

import edlib
import mappy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from scipy import stats
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def _sanitize(obj):
    """Recursively replace NaN/Inf floats with None for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    elif isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    return obj


def _compute_strand_stats(ed_results: list) -> dict:
    """Count strand orientations from ed_results."""
    strands = [r.get("strand", "fwd") for r in ed_results]
    c = Counter(strands)
    n = len(strands) or 1
    return {s: {"count": cnt, "pct": round(cnt / n * 100)} for s, cnt in c.items()}


# ---------------------------------------------------------------------------
# FASTQ parsing
# ---------------------------------------------------------------------------

def read_fastq_qnames(path: str) -> dict:
    """Return {qname: (seq, qual)} from a fastq file."""
    records = {}
    with pysam.FastxFile(filename=path) as fh:
        for rec in fh:
            records[rec.name] = (rec.sequence, rec.quality)
    return records


def glob_fastq_files(directory: str) -> list:
    """Glob *.fastq and *.fq files from a directory."""
    p = pathlib.Path(directory)
    files = []
    for ext in ("*.fastq", "*.fq", "*.FastQ", "*.FASTQ"):
        files.extend(p.glob(ext))
    for ext in ("*.fastq.gz", "*.fq.gz", "*.FastQ.gz"):
        files.extend(p.glob(ext))
    return sorted(files)


def _revcomp(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def _detect_strand_orientation(exp_seq: str, ctrl_seq: str) -> Tuple[str, str]:
    """Return (aligned_exp, aligned_ctrl) where both are in the same orientation.

    mappy.Aligner(seq=ctrl_seq) loads ctrl into an in-memory index.
    aligner.map(exp_seq) automatically probes both strands of ctrl —
    the primary hit's hit.strand tells us the relationship:
      1 → exp aligns to forward strand (same direction)
     -1 → exp aligns to reverse complement of ctrl (opposite)

    Falls back to edlib when mappy finds no decent hit.
    """
    aligner = mappy.Aligner(
        seq=ctrl_seq, extra_flags=67108864, k=15, w=3, best_n=20, n_threads=1,
    )

    for hit in aligner.map(exp_seq):
        if not hit.is_primary:
            continue
        qlen = hit.q_en - hit.q_st
        identity, total_b = 0.0, 0
        for length, op in hit.cigar:
            if op == 7:   # '='
                identity += length; total_b += length
            elif op == 8:  # 'X'
                total_b += length
        identity /= total_b if total_b else 1.0

        # Need substantial coverage and reasonable identity to trust mappy's strand call
        if identity >= 0.85 and qlen / max(len(exp_seq), 1) > 0.9:
            if hit.strand == -1:
                return _revcomp(exp_seq), ctrl_seq
            return exp_seq, ctrl_seq

    # No decent mappy hit — fall back to edlib on first 200bp
    rc_exp = _revcomp(exp_seq)
    f_ed = edlib.align(exp_seq[:200], ctrl_seq[:200], mode="NW", task="distance")
    r_ed = edlib.align(rc_exp[:200], ctrl_seq[:200], mode="NW", task="distance")
    if (r_ed["editDistance"] if r_ed["editDistance"] != -1 else 999) < (f_ed["editDistance"] if f_ed["editDistance"] != -1 else 999):
        return rc_exp, ctrl_seq
    return exp_seq, ctrl_seq


# ---------------------------------------------------------------------------
# Read matching and feature extraction
# ---------------------------------------------------------------------------

def match_read_pairs(
    exp_records: dict, ctrl_records: dict
) -> list:
    """Match reads by qname between two fastq record dicts.

    Returns list of (qname, aligned_exp_seq, aligned_ctrl_seq, strand) for
    shared qnames only.  strand is 'fwd' or 'rc'.
    """
    common = set(exp_records.keys()) & set(ctrl_records.keys())
    pairs = []
    for qname in sorted(common):
        exp_seq, _ = exp_records[qname]
        ctrl_seq, _ = ctrl_records[qname]

        if len(exp_seq) >= 10 and len(ctrl_seq) >= 10:
            aligned_exp, aligned_ctrl = _detect_strand_orientation(exp_seq, ctrl_seq)
            # Determine strand label: was exp RC'd to match ctrl?
            rc_of_exp = _revcomp(exp_seq)
            strand = "rc" if (aligned_exp == rc_of_exp) else "fwd"
            pairs.append((qname, aligned_exp, aligned_ctrl, strand))
    return pairs


def extract_windows(
    pairs: list, window: int
) -> list:
    """Extract head/tail windows from matched read pairs.

    Returns list of (qname, exp_head, exp_tail, ctrl_head, ctrl_tail, strand).
    """
    result = []
    for entry in pairs:
        if len(entry) == 4:
            qname, exp_seq, ctrl_seq, strand = entry
        else:
            qname, exp_seq, ctrl_seq = entry
            strand = "fwd"
        result.append((
            qname,
            exp_seq[:window],
            exp_seq[-window:],
            ctrl_seq[:window],
            ctrl_seq[-window:],
            strand,
        ))
    return result


# ---------------------------------------------------------------------------
# Edit distance computation (parallel)
# ---------------------------------------------------------------------------

def _compute_ed_single(item: Tuple[str, str, str, str, str, str]) -> dict:
    """Compute edit distances for one read pair's head and tail windows."""
    qname, exp_head, exp_tail, ctrl_head, ctrl_tail, strand = item

    # Head edit distance between aligned exp_seq and ctrl_seq
    head_res = edlib.align(exp_head, ctrl_head, mode="NW", task="distance")
    head_ed = head_res["editDistance"] if head_res["editDistance"] != -1 else 0

    # Tail edit distance between aligned exp_seq and ctrl_seq
    tail_res = edlib.align(exp_tail, ctrl_tail, mode="NW", task="distance")
    tail_ed = tail_res["editDistance"] if tail_res["editDistance"] != -1 else 0

    return {
        "qname": qname,
        "head_edit_distance": head_ed,
        "tail_edit_distance": tail_ed,
        "strand": strand,
        "head_len": len(exp_head),
        "tail_len": len(exp_tail),
    }


def compute_all_ed(pairs: list, window: int, n_workers: int = 4) -> list:
    """Compute edit distances for all pairs using multiprocessing."""
    if not pairs:
        return []

    tasks = [
        (q, eh, et, ch, ct, st)
        for entry in pairs
        for (q, eh, et, ch, ct, st) in ([entry + ("fwd",) if len(entry) == 5 else entry])
    ]

    results = []
    with mp.Pool(processes=n_workers) as pool:
        it = tqdm(
            pool.imap(_compute_ed_single, tasks, chunksize=5),
            total=len(tasks),
            desc="Computing edit distances",
        )
        for r in it:
            results.append(r)

    return results


def _compute_all_ed_sequential(pairs: list, window: int) -> list:
    """Compute edit distances for all pairs sequentially (single-threaded)."""
    if not pairs:
        return []
    return [_compute_ed_single(
        entry + ("fwd",) if len(entry) == 5 else entry
    ) for entry in pairs]


# ---------------------------------------------------------------------------
# Statistics computation
# ---------------------------------------------------------------------------

def compute_percentiles(ed_list: list) -> dict:
    """Compute key percentile values."""
    arr = np.array(ed_list)
    return {
        "p10": float(np.percentile(arr, 10)),
        "p25": float(np.percentile(arr, 25)),
        "p50": float(np.median(arr)),
        "p75": float(np.percentile(arr, 75)),
        "p90": float(np.percentile(arr, 90)),
        "p95": float(np.percentile(arr, 95)),
        "max": float(np.max(arr)),
    }


def compute_stats_for_array(ed_list: list) -> dict:
    """Compute mean/median/std for a single array of edit distances."""
    arr = np.array(ed_list)
    n = len(arr)
    return {
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr, ddof=1)) if n > 1 else 0.0,
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "n": n,
    }


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def _plot_ed_ratio(
    ed_values: list,
    yaxis: plt.Axes,
    title_text: str,
) -> None:
    """Draw a histogram of ED values with y-axis as ratio; annotate zero/one-ED ratio."""
    n = len(ed_values) or 1
    # Compute ED→ratio explicitly (not via plt.hist counts) so bars show ratios directly.
    ed_counts = Counter(ed_values)
    xs = sorted(ed_counts.keys())
    ys = [ed_counts[x] / n * 100 for x in xs]

    yaxis.bar(xs, ys, width=0.8, alpha=0.6, color="#e74c3c")

    # Annotate ED=0 and ED=1 if present
    ed_labels = {ed: f"{ed_counts[ed] / n * 100:.1f}%" for ed in (0, 1) if ed in ed_counts}
    for i, (ed, pct) in enumerate(ed_labels.items()):
        y_val = ys[xs.index(ed)]
        yaxis.text(
            ed, y_val + 2,
            f"ED={ed}: {pct}",
            ha="center", va="bottom", fontsize=10, color="#e74c3c",
        )

    # Set sensible y-limits so annotations aren't clipped
    max_y = max(ys) if ys else 0
    yaxis.set_ylim(bottom=0, top=max(max_y * 1.15, max_y + 5))

    yaxis.set_title(title_text)
    yaxis.set_xlabel("Edit Distance")
    yaxis.set_ylabel("Ratio (%)")
    yaxis.grid(True, alpha=0.3, axis="y")


def plot_histograms(
    ed_results: list,
    output_prefix: str,
    title: str = "",
) -> None:
    """Plot edit distance histograms for head and tail windows."""
    head_ed = [r["head_edit_distance"] for r in ed_results]
    tail_ed = [r["tail_edit_distance"] for r in ed_results]

    win_size = ed_results[0]["head_len"] if ed_results else 50
    n = len(head_ed) or 1

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    title_head = f"{title} — Head (window={win_size}, n={n})" if title else f"Head window edit distance (window={win_size}, n={n})"
    _plot_ed_ratio(head_ed, axes[0], title_head)

    title_tail = f"{title} — Tail (window={win_size}, n={n})" if title else f"Tail window edit distance (window={win_size}, n={n})"
    _plot_ed_ratio(tail_ed, axes[1], title_tail)

    plt.tight_layout()
    out_path = f"{output_prefix}_ed_hist.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    logging.info("Histogram saved to %s", out_path)


def plot_aggregate_ed_dist(
    all_ed_results: list,
    output_prefix: str,
    title: str = "",
) -> None:
    """Aggregate edit distances across all file pairs and plot overall head/tail distributions."""
    head_ed = []
    tail_ed = []
    for r in all_ed_results:
        if isinstance(r, dict) and r.get("error"):
            logging.warning("Skipping error result: %s", r.get("error"))
            continue
        if isinstance(r, dict) and "per_read" in r:
            pr = r["per_read"]
        else:
            pr = [r]
        head_ed.extend(x["head_edit_distance"] for x in pr)
        tail_ed.extend(x["tail_edit_distance"] for x in pr)

    if not all(v for v in (head_ed, tail_ed)):
        logging.warning("No aggregated data to plot")
        return

    n = len(head_ed) or 1
    win_size = 50

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    title_h = f"{title} — Head (aggregate, n={n})" if title else f"Head window edit distance (aggregated, n={n})"
    _plot_ed_ratio(head_ed, axes[0], title_h)

    title_t = f"{title} — Tail (aggregate, n={n})" if title else f"Tail window edit distance (aggregated, n={n})"
    _plot_ed_ratio(tail_ed, axes[1], title_t)

    plt.tight_layout()
    out_path = f"{output_prefix}_ed_dist.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    logging.info("Aggregated ed distribution saved to %s", out_path)


# ---------------------------------------------------------------------------
# Folder-mode matching helper
# ---------------------------------------------------------------------------

def correlate_folder_pairs(exp_dirs: list, ctrl_dirs: list) -> list:
    """Correlate exp_dir/ctrl_dir pairs by filename stem.

    Returns list of (exp_file, ctrl_file) tuples for matching files.
    """
    def _collect(dirs):
        stem_map = defaultdict(list)
        for d in dirs:
            for f in glob_fastq_files(d):
                stem = pathlib.Path(f).stem
                if stem.endswith(".gz"):
                    stem = stem[:-3]
                stem_map[stem].append(str(f))
        return stem_map

    exp_map = _collect(exp_dirs)
    ctrl_map = _collect(ctrl_dirs)

    pairs = []
    common_stems = set(exp_map.keys()) & set(ctrl_map.keys())
    for stem in sorted(common_stems):
        for ef in exp_map[stem]:
            for cf in ctrl_map[stem]:
                pairs.append((ef, cf))

    return pairs


# ---------------------------------------------------------------------------
# Main analysis pipeline
# ---------------------------------------------------------------------------

def _aggregate_all(all_results: list, file_pairs: list) -> dict:
    """Aggregate stats across all file pairs into a single summary."""
    total_exp = total_ctrl = total_matched = 0
    all_head_ed = []
    all_tail_ed = []
    all_strands = Counter()
    per_pair = []

    for r in all_results:
        head_ed = [x["head_edit_distance"] for x in r.get("per_read", [])]
        tail_ed = [x["tail_edit_distance"] for x in r.get("per_read", [])]
        total_exp += r.get("summary", {}).get("n_exp_reads", 0)
        total_ctrl += r.get("summary", {}).get("n_ctrl_reads", 0)
        total_matched += r.get("summary", {}).get("n_matched_pairs", 0)
        all_head_ed.extend(head_ed)
        all_tail_ed.extend(tail_ed)
        sd = r.get("summary", {}) or {}
        for k, v in (sd.get("strand_distribution") or {}).items():
            if isinstance(v, int):
                all_strands[k] += v

    n = len(all_head_ed) or 1
    agg = {
        "n_file_pairs": len(file_pairs),
        "total_exp_reads": total_exp,
        "total_ctrl_reads": total_ctrl,
        "total_matched_pairs": total_matched,
        "head_stats": compute_stats_for_array(all_head_ed),
        "tail_stats": compute_stats_for_array(all_tail_ed),
        "strand_distribution": dict(all_strands),
    }

    # Paired comparison from pooled data
    head_arr = np.array(all_head_ed)
    tail_arr = np.array(all_tail_ed)
    mean_diff = float(np.mean(head_arr - tail_arr))
    paired_t_result = stats.ttest_rel(head_arr, tail_arr)
    diff_arr = head_arr - tail_arr
    nonzero_mask = diff_arr != 0
    if np.sum(nonzero_mask) >= 2:
        wilcoxon_result = stats.wilcoxon(diff_arr[nonzero_mask])
        wilcoxon_pvalue = float(wilcoxon_result.pvalue)
    else:
        wilcoxon_result = None
        wilcoxon_pvalue = np.nan

    agg["head_vs_tail_comparison"] = {
        "mean_diff_head_minus_tail": mean_diff,
        "paired_t_statistic": float(paired_t_result.statistic),
        "paired_t_pvalue": float(paired_t_result.pvalue),
        "wilcoxon_statistic": float(wilcoxon_result.statistic) if wilcoxon_result else None,
        "wilcoxon_pvalue": wilcoxon_pvalue,
        "percent_reads_ed_head_gt_tail": float(np.mean(head_arr > tail_arr) * 100),
    }

    exp_mean_ed = float(np.mean(head_arr))
    significant_head_more = bool(paired_t_result.pvalue < 0.05 and mean_diff > 0)
    agg["divergence_summary"] = {
        "overall_mean_edit_distance": exp_mean_ed,
        "strategy_affects_head_more": significant_head_more,
        "significant_at_005": bool(paired_t_result.pvalue < 0.05),
    }

    # Keep per-pair summaries (not full per_read) for reference
    agg["per_pair"] = [
        {k: v for k, v in r.items() if k != "per_read"}
        for r in all_results
    ]
    return agg


def run_analysis(
    exp_fq_path: str,
    ctrl_fq_path: str,
    window: int = 50,
    n_workers: int = 4,
    use_pool: bool = True,
    output_prefix: str = "head_tail_ab_results",
) -> dict:
    """Run the full head/tail AB analysis on a single pair of fastq files.

    Returns the stats dict (JSON/PNG writing is done by the caller).

    For each matched read pair (exp_seq, ctrl_seq):
      - head_ed = edit_distance(exp[:window], ctrl[:window])
      - tail_ed = edit_distance(exp[-window:], ctrl[-window:])

    use_pool: if True, use mp.Pool for query-level parallelism (file mode);
              if False, process queries sequentially (folder mode already has
              its own ProcessPoolExecutor at the pair level).
    """
    logging.info("Reading %s ...", exp_fq_path)
    exp_records = read_fastq_qnames(exp_fq_path)
    logging.info("  -> %d reads", len(exp_records))

    logging.info("Reading %s ...", ctrl_fq_path)
    ctrl_records = read_fastq_qnames(ctrl_fq_path)
    logging.info("  -> %d reads", len(ctrl_records))

    pairs = match_read_pairs(exp_records, ctrl_records)
    logging.info("Matched %d read pairs", len(pairs))

    if not pairs:
        logging.error("No matched reads — check file contents or --window size")
        return {"error": "no_matched_pairs"}

    windowed = extract_windows(pairs, window)
    if use_pool:
        ed_results = compute_all_ed(windowed, window, n_workers)
    else:
        ed_results = _compute_all_ed_sequential(windowed, window)

    if not ed_results:
        logging.error("No valid edit distance results")
        return {"error": "no_ed_results"}

    # Extract per-pair values
    head_ed = [r["head_edit_distance"] for r in ed_results]
    tail_ed = [r["tail_edit_distance"] for r in ed_results]
    head_arr = np.array(head_ed)
    tail_arr = np.array(tail_ed)

    # Paired comparison: does head diverge more than tail? (ed_head > ed_tail?)
    paired_t_result = stats.ttest_rel(head_arr, tail_arr)

    diff_arr = head_arr - tail_arr
    nonzero_mask = diff_arr != 0
    if np.sum(nonzero_mask) >= 2:
        wilcoxon_result = stats.wilcoxon(diff_arr[nonzero_mask])
        wilcoxon_pvalue = float(wilcoxon_result.pvalue)
    else:
        wilcoxon_result = None
        wilcoxon_pvalue = np.nan

    # Strategy-change judgment
    mean_diff = float(np.mean(head_arr - tail_arr))
    exp_mean_ed = float(np.mean(head_arr))  # overall avg divergence (exp vs ctrl)
    significant_head_more = bool(paired_t_result.pvalue < 0.05 and mean_diff > 0)

    stats_output = {
        "input": {
            "exp_fq": exp_fq_path,
            "ctrl_fq": ctrl_fq_path,
            "window_size": window,
        },
        "summary": {
            "n_exp_reads": len(exp_records),
            "n_ctrl_reads": len(ctrl_records),
            "n_matched_pairs": len(ed_results),
            "strand_distribution": _sanitize(_compute_strand_stats(ed_results)),
        },
        # Edit distance stats for each region (head/tail)
        "head_stats": {
            **compute_stats_for_array(head_ed),
            **compute_percentiles(head_ed),
        },
        "tail_stats": {
            **compute_stats_for_array(tail_ed),
            **compute_percentiles(tail_ed),
        },
        # Paired comparison: head vs tail (key for judging strategy change)
        "head_vs_tail_comparison": {
            "mean_diff_head_minus_tail": mean_diff,
            "paired_t_statistic": float(paired_t_result.statistic),
            "paired_t_pvalue": float(paired_t_result.pvalue),
            "wilcoxon_statistic": float(wilcoxon_result.statistic) if wilcoxon_result else None,
            "wilcoxon_pvalue": wilcoxon_pvalue,
            "percent_reads_ed_head_gt_tail": float(np.mean(head_arr > tail_arr) * 100),
        },
        # Overall divergence between exp and ctrl (for decision-making)
        "divergence_summary": {
            "overall_mean_edit_distance": exp_mean_ed,
            "strategy_affects_head_more": significant_head_more,
            "significant_at_005": bool(paired_t_result.pvalue < 0.05),
            # Simple judgment: if overall divergence is low AND head/tail diff is not significant
            # the strategy change is safe to deploy; otherwise investigate
        },
    }

    # Sanitize NaN/Inf to None so JSON output is valid
    stats_output = _sanitize(stats_output)

    # Add per-read details to JSON (for debugging / downstream)
    stats_output["per_read"] = ed_results

    return stats_output


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main_cli():
    parser = argparse.ArgumentParser(
        description="Compare SMC strategy changes via head/tail edit distance analysis.",
    )

    # Either file mode or folder mode
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--exp-fq", help="Path to experimental fastq file")
    grp.add_argument("--exp-dirs", help="Comma-separated list of experiment group directories")
    grp.add_argument("--exp-dir", help="Single experiment directory")

    grp2 = parser.add_mutually_exclusive_group(required=True)
    grp2.add_argument("--ctrl-fq", help="Path to control fastq file")
    grp2.add_argument("--ctrl-dirs", help="Comma-separated list of control group directories")
    grp2.add_argument("--ctrl-dir", help="Single control directory")

    parser.add_argument("-o", "--output-prefix", default="head_tail_ab_results", help="Output file prefix")
    parser.add_argument("--window", type=int, default=50, help="Head/tail window size (default: 50)")
    parser.add_argument("--workers", type=int, default=min(mp.cpu_count(), 8), help="Number of parallel workers")

    args = parser.parse_args()

    if args.exp_dir:
        exp_dirs = [args.exp_dir]
    elif args.exp_dirs:
        exp_dirs = [d.strip() for d in args.exp_dirs.split(",")]
    else:
        exp_dirs = [args.exp_fq]

    if args.ctrl_dir:
        ctrl_dirs = [args.ctrl_dir]
    elif args.ctrl_dirs:
        ctrl_dirs = [d.strip() for d in args.ctrl_dirs.split(",")]
    else:
        ctrl_dirs = [args.ctrl_fq]

    # File mode vs folder mode
    if len(exp_dirs) == 1 and os.path.isfile(exp_dirs[0]) and len(ctrl_dirs) == 1 and os.path.isfile(ctrl_dirs[0]):
        logging.info("File mode: %s vs %s", exp_dirs[0], ctrl_dirs[0])
        result = run_analysis(
            exp_fq_path=exp_dirs[0],
            ctrl_fq_path=ctrl_dirs[0],
            window=args.window,
            n_workers=args.workers,
            output_prefix=args.output_prefix,
        )
        if isinstance(result, dict) and result.get("error"):
            sys.exit(1)

        # Write outputs (file mode only)
        json_path = f"{args.output_prefix}_stats.json"
        with open(json_path, "w") as f:
            json.dump(result, f, indent=2)
        logging.info("Stats saved to %s", json_path)

        ed_results = result.get("per_read", [])
        plot_histograms(ed_results, args.output_prefix, title=f"{os.path.basename(exp_dirs[0])} vs {os.path.basename(ctrl_dirs[0])}")
    else:
        # Folder mode — correlate by filename stem
        file_pairs = correlate_folder_pairs(exp_dirs, ctrl_dirs)
        if not file_pairs:
            logging.error("No matching files found across directories")
            sys.exit(1)

        n_workers_folder = min(args.workers, len(file_pairs))
        logging.info("Folder mode: %d file pairs, using %d workers", len(file_pairs), n_workers_folder)

        all_results = []
        with ProcessPoolExecutor(max_workers=n_workers_folder) as executor:
            future_to_pair = {
                executor.submit(
                    run_analysis,
                    exp_fq_path=ef,
                    ctrl_fq_path=cf,
                    window=args.window,
                    n_workers=args.workers,  # query-level parallelism within each worker
                    use_pool=False,           # file mode: still parallelize queries internally
                    output_prefix=args.output_prefix,
                ): (ef, cf)
                for ef, cf in file_pairs
            }
            for future in tqdm(as_completed(future_to_pair), total=len(file_pairs), desc="Processing pairs"):
                ef, cf = future_to_pair[future]
                try:
                    result = future.result()
                    if not isinstance(result, dict) or result.get("error"):
                        logging.warning("Skipping pair with error: %s vs %s", ef, cf)
                        continue
                    all_results.append(result)
                except Exception as exc:
                    logging.error("Pair %s vs %s generated an exception: %s", ef, cf, exc)

        # Aggregate stats across all pairs
        agg = _aggregate_all(all_results, file_pairs)
        agg_path = f"{args.output_prefix}_aggregate.json"
        with open(agg_path, "w") as f:
            json.dump(agg, f, indent=2)
        logging.info("Aggregate stats saved to %s", agg_path)

        # Aggregate edit distance histogram across all pairs
        plot_aggregate_ed_dist(all_results, args.output_prefix, title="Aggregated")


if __name__ == "__main__":
    main_cli()
