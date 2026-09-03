#!/usr/bin/env python3
"""Search one or more adapter BAMs for reads harboring a query sequence.

For every read, an alignment of the full query against a contiguous substring
of the read is computed with edlib (mode ``HW`` + task ``path``: the entire
query must align somewhere inside the read, allowing mismatches and
insertions/deletions). A read is reported when the edit distance is
<= ``(1 - sim) * len(query)``, i.e. the substring similarity is >= ``--sim``.

Search scope: every ``*.adapter.bam`` file under the search root (including
subdirectories). Results are written as a TSV, sorted by descending
similarity. ``read_aln_start``/``read_aln_end`` are the 1-based inclusive
coordinates of the aligned window on the read.
"""

import argparse
import multiprocessing as mp
from pathlib import Path

import edlib
import pysam
from tqdm import tqdm


def _collect_adapter_bams(root: Path) -> list:
    """Return every *.adapter.bam under `root`, sorted for determinism."""
    bams = sorted(str(p) for p in root.rglob("*.adapter.bam"))
    if not bams:
        raise SystemExit(f"no *.adapter.bam found under {root}")
    return bams


def search_read(query: str, max_edit: int, seq: str, name: str):
    """Return a result tuple if `seq` contains the query with <= `max_edit` edits."""
    if not seq:
        return None
    # k= gives a banded search: distance -1 means no alignment <= max_edit.
    result = edlib.align(query=query, target=seq, mode="HW", task="path", k=max_edit)
    if result["editDistance"] is None or result["editDistance"] < 0:
        return None
    aln_start, aln_end = result["locations"][0]  # 0-based on the read
    return (
        name,
        len(seq),
        result["editDistance"],
        (1.0 - result["editDistance"] / len(query)) * 100.0,
        aln_start + 1,  # 1-based, inclusive
        aln_end,
    )


def search_bam(
    bam_path: str,
    query: str,
    max_edit: int,
    threads: int,
    show_progress: bool,
):
    """Scan one BAM for reads containing the query; return per-read results."""
    results = []
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=threads) as aln:
        iterator = aln.fetch(until_eof=True)
        if show_progress:
            iterator = tqdm(iterator, desc=str(Path(bam_path).name))
        for read in iterator:
            res = search_read(query, max_edit, read.query_sequence, read.query_name)
            if res is not None:
                results.append(res)
    return bam_path, results


def search_adapter_bams(
    root: str,
    query: str,
    sim: float,
    threads: int,
    show_progress: bool = True,
):
    """Search all *.adapter.bam under `root` and return a sorted result list.

    Each result is (bam_path, read_name, read_len, edits, similarity_pct,
    read_aln_start, read_aln_end).
    """
    query = query.upper()
    max_edit = int((1.0 - sim) * len(query))
    if max_edit < 0:
        raise ValueError(f"similarity {sim} impossible for a {len(query)}-bp query")

    bams = _collect_adapter_bams(Path(root))
    pool = mp.Pool(min(threads, len(bams)))
    rows = []
    try:
        for bam_path, results in pool.starmap(
            search_bam,
            [(bam, query, max_edit, 1, show_progress) for bam in bams],
        ):
            for r in results:
                rows.append((bam_path,) + r)
    finally:
        pool.close()
        pool.join()

    rows.sort(key=lambda r: r[4], reverse=True)
    return rows


def write_tsv(results, out_path: str) -> None:
    """Write results as a sorted TSV."""
    header = (
        "bam\tread_name\tread_len\tedits\tsimilarity_pct\t"
        "read_aln_start\tread_aln_end"
    )
    lines = [header]
    lines += [
        "\t".join([
            bam,
            name,
            str(read_len),
            str(edits),
            f"{sim_pct:.2f}",
            str(aln_start),
            str(aln_end),
        ])
        for (bam, name, read_len, edits, sim_pct, aln_start, aln_end) in results
    ]
    Path(out_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main_cli():
    parser = argparse.ArgumentParser(
        description="Search *.adapter.bam files under a directory for reads "
                    "containing a query sequence with similarity >= --sim.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "root",
        help="Directory to search recursively for *.adapter.bam files",
    )
    parser.add_argument(
        "query",
        help="Query sequence (DNA/RNA, case-insensitive)",
    )
    parser.add_argument(
        "--sim", type=float, default=0.9,
        help="Minimum similarity of the query against a read substring",
    )
    parser.add_argument(
        "-o", "--output",
        help="Output TSV path (default: ./seq_search_results.tsv)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=mp.cpu_count(),
        help="Number of parallel BAM-processing workers",
    )
    parser.add_argument(
        "--no-progress", action="store_true",
        help="Disable per-BAM progress bars",
    )
    args = parser.parse_args()

    out_path = args.output or "./seq_search_results.tsv"
    show_progress = not args.no_progress

    print(f"query length: {len(args.query)} bp")
    print(f"threshold: similarity >= {args.sim} "
          f"(edits <= {(1.0 - args.sim) * len(args.query):.0f})", flush=True)

    results = search_adapter_bams(
        args.root, args.query, args.sim, args.threads, show_progress
    )

    write_tsv(results, out_path)
    print(f"\n{len(results)} reads matched; results written to {out_path}")


if __name__ == "__main__":
    main_cli()
