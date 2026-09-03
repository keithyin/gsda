import mappy as mp
import argparse
import re


def calculate_identity_eqx(cigar_str):
    """
    Calculate the identity based on the CIGAR string with `=` and `X` used for match and mismatch respectively.

    Parameters:
    cigar_str (str): The CIGAR string representing the alignment (e.g. '8=2X4=1D5=')

    Returns:
    float: The percentage identity based on the CIGAR string.
    """
    match_count = 0
    mismatch_count = 0
    total_bases = 0

    # Regular expression to match segments in the CIGAR string (e.g. '8=', '2X', '4=', etc.)
    cigar_regex = re.compile(r"(\d+)([=XDI])")

    for match in cigar_regex.finditer(cigar_str):
        length = int(match.group(1))  # Length of this segment
        operation = match.group(2)  # Operation type (=, X, D, etc.)

        if operation == "=":  # Exact match
            match_count += length
            total_bases += length
        elif operation == "X":  # Mismatch
            mismatch_count += length
            total_bases += length
        elif operation == "I":  # Insertion
            total_bases += length
        elif operation == "D":  # Deletion
            total_bases += length
        # For other operations like 'S' (soft clip), 'H' (hard clip), etc., don't count towards total bases.

    if total_bases == 0:
        return 0.0

    # Identity is the fraction of matching bases out of total bases considered in the alignment.
    identity = (match_count / total_bases) * 100
    return identity


def merge_intervals_and_calculate_length(intervals):
    if len(intervals) == 0:
        return [], 0, 0
    # 排序区间，按起始位置排序
    intervals.sort(key=lambda x: x[0])

    merged = []
    total_overlap_length = 0
    prev_start, prev_end = intervals[0]

    for current_start, current_end in intervals[1:]:
        # 如果当前区间与上一个区间重叠或相邻
        if current_start <= prev_end:
            # 计算重叠部分的长度
            overlap_length = min(prev_end, current_end) - current_start
            total_overlap_length += overlap_length

            # 合并区间，更新结束位置
            prev_end = max(prev_end, current_end)
        else:
            # 如果没有重叠，直接添加上一个区间到合并区间
            merged.append([prev_start, prev_end])
            prev_start, prev_end = current_start, current_end

    # 添加最后一个区间
    merged.append([prev_start, prev_end])

    # 计算每个合并后区间的长度
    merged_length = sum([end - start for start, end in merged])

    return merged, merged_length, total_overlap_length


def find_uncovered_part(big_interval, merged_intervals):
    uncovered_parts = []
    big_start, big_end = big_interval

    # 如果没有合并区间或合并区间为空
    if not merged_intervals:
        return [big_interval]

    # 遍历合并区间，计算未被覆盖的部分
    # 1. 检查大区间开始到第一个合并区间之前的部分
    if big_start < merged_intervals[0][0]:
        uncovered_parts.append([big_start, merged_intervals[0][0]])

    # 2. 检查两个合并区间之间的部分
    for i in range(1, len(merged_intervals)):
        prev_end = merged_intervals[i - 1][1]
        current_start = merged_intervals[i][0]

        if prev_end < current_start:
            uncovered_parts.append([prev_end, current_start])

    # 3. 检查最后一个合并区间到大区间结束的部分
    if merged_intervals[-1][1] < big_end:
        uncovered_parts.append([merged_intervals[-1][1], big_end])

    return uncovered_parts


class SingleQueryHits:
    def __init__(
        self, query_name: str, ref_name: str, query_length: int, ref_length: int, hits
    ):
        self.query_name = query_name
        self.ref_name = ref_name
        self.query_length = query_length
        self.ref_length = ref_length
        self.hits = hits

        self.identity = None
        self.query_coverage = None
        self.query_uncovered = None
        self.query_ovlp = None

        self.ref_coverage = None
        self.ref_uncovered = None
        self.ref_ovlp = None

    def fill_field(self):
        cigar_strs = [hit.cigar_str for hit in self.hits]
        self.identity = calculate_identity_eqx(cigar_str="".join(cigar_strs))

        query_start_ends = [[hit.q_st, hit.q_en] for hit in self.hits]
        merged, merged_length, total_overlap_length = (
            merge_intervals_and_calculate_length(query_start_ends)
        )
        self.query_ovlp = total_overlap_length / max(merged_length, 1)
        self.query_uncovered = find_uncovered_part(
            [0, self.query_length], merged_intervals=merged
        )
        self.query_coverage = merged_length / self.query_length

        ref_start_ends = [[hit.r_st, hit.r_en] for hit in self.hits]
        merged, merged_length, total_overlap_length = (
            merge_intervals_and_calculate_length(ref_start_ends)
        )
        self.ref_ovlp = total_overlap_length / max(merged_length, 1)
        self.ref_uncovered = find_uncovered_part(
            [0, self.ref_length], merged_intervals=merged
        )
        self.ref_coverage = merged_length / self.ref_length

    def to_str(self) -> str:
        res = f"""
        -------------------------------------------------------------------
        - query_name: {self.query_name}
        - ref_name: {self.ref_name}
        - identity: {self.identity}
        - 
        - query_coverage: {self.query_coverage}
        - query_uncovered: {self.query_uncovered}
        - query_ovlp: {self.query_ovlp}
        -
        - ref_coverage: {self.ref_coverage}
        - ref_uncovered: {self.ref_uncovered}
        - ref_ovlp: {self.ref_ovlp}
        -------------------------------------------------------------------
        """
        return res


def main(ref_fastx: str, query_fastx: str, preset="map-ont"):
    # 67108864 eqx
    aligner = mp.Aligner(ref_fastx, preset=preset, extra_flags=67108864)
    assert aligner.n_seq == 1
    ref_seq_name = aligner.seq_names[0]
    ref_seq = aligner.seq(ref_seq_name)

    for name, seq, _ in mp.fastx_read(query_fastx):
        hits = list(aligner.map(seq))
        single_query_hits = SingleQueryHits(
            query_name=name,
            ref_name=ref_seq_name,
            query_length=len(seq),
            ref_length=len(ref_seq),
            hits=hits,
        )
        single_query_hits.fill_field()

        print(single_query_hits.to_str())


def main_cli():
    parser = argparse.ArgumentParser(prog="WGA")
    parser.add_argument("ref", type=str, help="ref fastx")
    parser.add_argument("query", type=str, help="query fastx")
    parser.add_argument("--preset", type=str, default="map-ont")
    args = parser.parse_args()

    main(args.ref, args.query, args.preset)


if __name__ == "__main__":
    main_cli()
