"""Dump query names of reads with deletions at specified reference positions."""

import argparse
from collections import defaultdict


def dump_deletion_reads(bam_file: str, positions_str: str, num: int = 10) -> None:
    # Parse comma-separated positions
    req_positions: set[int] = set(
        int(p.strip()) for p in positions_str.split(",") if p.strip()
    )

    import pysam

    with pysam.AlignmentFile(bam_file, "rb", threads=40) as bam:
        # Collect all contig names to iterate over
        ref_names = bam.references
        ref_lengths = bam.lengths

        pos2names: dict[int, list[str]] = defaultdict(list)
        pending: set[int] = set(req_positions)

        for i, ref_name in enumerate(ref_names):
            ref_end = ref_lengths[i]
            if not any(0 <= p < ref_end for p in pending):
                continue

            min_pos = min(p for p in pending if 0 <= p < ref_end)
            max_pos = max(p for p in pending if 0 <= p < ref_end)

            for pileup_col in bam.pileup(
                contig=ref_name,
                start=min_pos,
                end=max_pos + 1,
                min_base_quality=0,
                truncate_ends=False,
            ):
                pos = pileup_col.reference_pos
                if pos not in pending:
                    continue

                for query in pileup_col.pileups:
                    if query.is_del:
                        qname = query.alignment.query_name
                        pos2names[pos].append(qname)

                        if len(pos2names[pos]) >= num:
                            pending.discard(pos)
                            break

                if not pending:
                    break

        # Print results sorted by position
        for pos in sorted(pos2names):
            print(f"REF_POS {pos}:")
            for qname in pos2names[pos][:num]:
                print(qname)
            if len(pos2names[pos]) < num:
                print(f"  (only {len(pos2names[pos])} found)")


def main_cli():
    parser = argparse.ArgumentParser(
        description="Dump query names of reads with deletions at specified reference positions."
    )
    parser.add_argument("--bam", required=True, help="Input aligned BAM file.")
    parser.add_argument(
        "--positions", required=True, help="Comma-separated reference positions (e.g. 100,200,300)."
    )
    parser.add_argument(
        "--num",
        type=int,
        default=10,
        help="Max query names per position (default: 10).",
    )
    args = parser.parse_args()
    dump_deletion_reads(args.bam, args.positions, args.num)


if __name__ == "__main__":
    main_cli()
