#!/usr/bin/env python3
"""
将比对到 ref 某个区间的 所有的 query 的序列抽取出来。

输入是比对的 bam 文件，以及 ref name & start/end 区间。
输出为 FASTA 格式（仅包含 query 中对应 ref start-end 区间的子序列）。

使用 pysam
"""

import argparse
import sys
from collections import Counter
import pysam
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')  # headless backend for servers without display
    import matplotlib.pyplot as plt
    _HAS_MATPLOTLIB = True
except ImportError:
    _HAS_MATPLOTLIB = False


def contig_cnt(seq, motif):
    """
    Count non-overlapping occurrences of `motif` in `seq`.

    Greedy left-to-right: after a match at position i, the next search
    resumes at i + len(motif), so overlapping matches are not counted.

    Examples
        >>> contig_cnt("GGCGGGGGGC", "GGC")
        2
        >>> contig_cnt("GGCGGCGGGC", "GGC")
        3
    """
    if not motif:
        return 0

    n = 0
    start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1:
            break
        n += 1
        start = pos + len(motif)  # advance past this match to avoid overlaps
    return n


def extract_subseq(read, ref_start, ref_end):
    """
    Extract the query subsequence corresponding to reference interval
    [ref_start, ref_end).

    Rules
    -----
    * Read must fully span the reference interval.
    * Only aligned query bases are emitted.
    * Deletions (D/N) contribute nothing.
    * Insertions are ignored because they have no reference coordinate.
    """

    # Read must fully cover the requested interval
    if read.reference_start > ref_start:
        return ""

    if read.reference_end < ref_end:
        return ""

    query = read.query_sequence
    chars = []
    q_cursor = -1
    r_cursor = -1
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if qpos is not None:
            q_cursor = qpos

        if rpos is not None:
            r_cursor = rpos

        if q_cursor is None and r_cursor is None:
            continue

        if r_cursor < (ref_start - 1):
            continue

        if r_cursor == (ref_start - 1) and rpos is not None:
            continue

        if r_cursor < ref_start:
            continue

        if r_cursor >= ref_end:
            break

        # deletion
        if qpos is None:
            continue

        chars.append(query[qpos])

    return "".join(chars)


def print_ref_view(ref_name, ref_seq, start, end):
    """Print the reference sequence with 10bp flanking and target boxed in []."""
    flank = 10
    view_start = max(0, start - flank)
    view_end = min(len(ref_seq), end + flank)
    prefix = ref_seq[view_start:start]
    target = ref_seq[start:end]
    suffix = ref_seq[end:view_end]
    boxed = f"{prefix}[{target}]{suffix}"
    print(f">{ref_name}_{start}_{end}", file=sys.stderr)
    # Split into 60bp-wide chunks for readability
    for i in range(0, len(boxed), 60):
        print(boxed[i: i + 60], file=sys.stderr)


def main_cli():
    parser = argparse.ArgumentParser(
        prog="extract_align_ref_region",
        description="Extract query subsequences from BAM that overlap a reference region.",
    )
    parser.add_argument("bam_file", help="aligned BAM file")
    parser.add_argument("--ref-name", required=True,
                        help="reference name (contig/chromosome)")
    parser.add_argument("--start", type=int, required=True,
                        help="region start (0-based, inclusive)")
    parser.add_argument("--end", type=int, required=True,
                        help="region end (exclusive)")
    parser.add_argument("--ref-file", default=None,
                        help="fasta file containing the reference sequence")
    parser.add_argument("-o", "--output", default=None,
                        help="output FASTA file (default: stdout)")
    parser.add_argument("--output-plot", default="ggc_cnts_distribution.png",
                        help="output path for the plot image (default: ggc_cnts_distribution.png)")

    args = parser.parse_args()

    end_val = args.end  # default, may be clamped by ref file length

    if args.ref_file is not None:
        with pysam.FastxFile(args.ref_file) as fx:
            found = False
            for entry in fx:
                if entry.name == args.ref_name:
                    ref_seq = entry.sequence
                    found = True
                    break
        if not found:
            print(
                f"error: reference '{args.ref_name}' not found in {args.ref_file}", file=sys.stderr)
            sys.exit(1)
        end_val = min(args.end, len(ref_seq))
        if end_val < args.end:
            print(
                f"warning: target end ({args.end}) exceeds ref length ({len(ref_seq)}), using {end_val}", file=sys.stderr)
        print_ref_view(args.ref_name, ref_seq, args.start, end_val)

    # sys.exit(1)
    fout = open(args.output, "w",
                encoding="utf8") if args.output is not None else sys.stdout

    motif = "GC"
    with pysam.AlignmentFile(args.bam_file, "rb", check_sq=False) as bam:
        n = 0
        total_len = 0
        reads_lengths = []
        ggc_cnts = []
        for read in bam.fetch(contig=args.ref_name, start=args.start, end=end_val):
            subseq = extract_subseq(read, args.start, end_val)
            if len(subseq) > 0:
                total_len += len(subseq)
                reads_lengths.append(len(subseq))
                # if read.query_name.startswith("20260714_250701Y0006_Run0001/423882"):
                if "T" in subseq:
                    fout.write(f">{read.query_name}: {len(subseq)} \n{subseq}\n")
                n += 1

                ggc_cnt = contig_cnt(seq=subseq, motif=motif)
                ggc_cnts.append(ggc_cnt)

        print(f"""
              length_mean:     {total_len  / n} 
              length_median:   {np.median(reads_lengths)}
              length_mean/3:   {total_len / 3 / n}
              length_median/3: {np.median(reads_lengths) / 3 }
              
              """)

        reads_lengths = np.array(reads_lengths)
        ggc_cnts = np.array(ggc_cnts)
        print(f"{motif}_cnts, median: ", np.median(ggc_cnts))
        print(f"{motif}_cnts, mean: ", np.mean(ggc_cnts))
        print(f"Num Records: {len(ggc_cnts)}")

        # Draw GGC count distribution plot (opt-in via --plot)
        if _HAS_MATPLOTLIB and len(ggc_cnts) > 0:
            cnt = Counter(ggc_cnts)
            keys = sorted(cnt.keys())
            values = [cnt[k] for k in keys]
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.bar(keys, values, width=0.6, edgecolor='black')

            ggc_arr = np.array(ggc_cnts)
            mean_val = np.mean(ggc_arr)
            median_val = np.median(ggc_arr)
            xmax = max(keys) if keys else 0
            ymax = max(values) if values else 0
            ax.axvline(mean_val, color='red', linestyle='--', linewidth=1.2,
                       label=f'mean={mean_val:.2f}')
            ax.axvline(median_val, color='green', linestyle='-.', linewidth=1.2,
                       label=f'median={median_val:.2f}')
            ax.legend(loc='upper left')

            # Add vertical line annotations at the top of the bars
            for x, lbl in [(mean_val, 'mean'), (median_val, 'median')]:
                if xmax > 0:
                    ax.annotate(lbl, xy=(x, ymax), xytext=(x, ymax + ymax * 0.1),
                                ha='center', fontsize=9, color={'mean': 'red', 'median': 'green'}[lbl])

            ax.set_xlabel(f'{motif} count')
            ax.set_ylabel('Frequency (number of reads)')
            ax.set_title(f'{motif} Motif Count Distribution')
            fig.tight_layout()
            fig.savefig(args.output_plot, dpi=150)
            plt.close(fig)
            print(f"plot saved to {args.output_plot}", file=sys.stderr)
        elif not _HAS_MATPLOTLIB:
            print("warning: matplotlib not installed; install with 'pip install matplotlib' to enable plotting", file=sys.stderr)

    if args.output is not None:
        fout.close()

    if n == 0:
        print("warning: no reads found in region", file=sys.stderr)


if __name__ == "__main__":
    main_cli()
