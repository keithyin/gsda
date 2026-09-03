#!/usr/bin/env python3

from collections import defaultdict

def print_summary_table(hist, min_len=2):
    """
    输出：
    Repeat Count    A    C    G    T    合计
    2              26   53   35   25   139
    3               5    9    9    5    28
    ...
    """

    # 所有出现过的 repeat length
    lengths = sorted({
        length
        for base_hist in hist.values()
        for length in base_hist.keys()
    })

    print("\n" + "=" * 60)
    print("Homopolymer Summary")
    print("=" * 60)

    print(
        f"{'Repeat Count':>12}"
        f"{'A':>8}"
        f"{'C':>8}"
        f"{'G':>8}"
        f"{'T':>8}"
        f"{'Total':>10}"
    )

    for length in lengths:
        if length < min_len:
            continue

        a = hist.get("A", {}).get(length, 0)
        c = hist.get("C", {}).get(length, 0)
        g = hist.get("G", {}).get(length, 0)
        t = hist.get("T", {}).get(length, 0)

        total = a + c + g + t

        print(
            f"{length:>12}"
            f"{a:>8}"
            f"{c:>8}"
            f"{g:>8}"
            f"{t:>8}"
            f"{total:>10}"
        )


def find_homopolymers(seq: str):
    """
    返回所有 homopolymer run

    每个元素:
    {
        "base": 'A',
        "start": 123,
        "end": 126,
        "length": 4
    }

    start/end 为 0-based, end 为 inclusive
    """

    seq = seq.upper()

    if not seq:
        return []

    runs = []

    start = 0
    base = seq[0]

    for i in range(1, len(seq)):
        if seq[i] != base:
            runs.append({
                "base": base,
                "start": start,
                "end": i - 1,
                "length": i - start,
            })

            start = i
            base = seq[i]

    runs.append({
        "base": base,
        "start": start,
        "end": len(seq) - 1,
        "length": len(seq) - start,
    })

    return runs


def build_histogram(runs):
    """
    per-base histogram

    hist['A'][3] = AAA 出现次数
    """

    hist = defaultdict(lambda: defaultdict(int))

    for r in runs:
        hist[r["base"]][r["length"]] += 1

    return hist


def build_total_histogram(runs):
    total = defaultdict(int)

    for r in runs:
        total[r["length"]] += 1

    return total


def print_histogram(hist):
    print("=" * 60)
    print("Per-base homopolymer histogram")
    print("=" * 60)

    for base in "ACGT":
        print(f"\n{base}")

        if base not in hist:
            print("  None")
            continue

        for length in sorted(hist[base]):
            print(f"  {length:2d} : {hist[base][length]}")


def print_total_histogram(total):
    print("\n" + "=" * 60)
    print("Overall histogram")
    print("=" * 60)

    for length in sorted(total):
        print(f"{length:2d} : {total[length]}")


def print_long_runs(runs, min_len=2):
    print("\n" + "=" * 60)
    print(f"Homopolymers (length >= {min_len})")
    print("=" * 60)

    print(f"{'Base':<5} {'Start':>8} {'End':>8} {'Length':>8}")

    for r in runs:
        if r["length"] >= min_len:
            print(
                f"{r['base']:<5}"
                f"{r['start']:>8}"
                f"{r['end']:>8}"
                f"{r['length']:>8}"
            )


def print_longest(runs):
    longest = max(runs, key=lambda x: x["length"])

    print("\n" + "=" * 60)
    print("Longest homopolymer")
    print("=" * 60)

    print(longest)


def main():

    seq = """
CCAAGATGGAGTGAAAGTGTCTCAGGGGAATCGACGCGATCGAGATCGGTCGATATTCCCCCGGTCCACCTGGTCTCTCCACCTTAGGTACAAAGACGGTTCGGCACTGCCGTAGAATTTCGGGTATTTCGCCTCGTGCCATCCATGCATTGAACATTTCCGCCTTCAAGTGCACAGGAACCGCACGCCACTGACCCGAACGTATACCGTCCGGGCCCGGCGAAGTTCGCCAGTCAAAGCGGGAGGCCTTGATCTCTTCCACCGAGATCGGCTTCCACAGCTGGGTGTAGTCGCGATTAGACCCAATCGCTCAGTCAGGTAGTGATGCTTCTTGCGGTGGCCATCCGGTATCATAGACAGAACCGCGCTTTTCCTAGGTTTCAGGCGTAGGCCCATCTGCCTACCGACACAGTCCACAGCAGAGATGGACTCCTGCATCCCTACCTTCGACCCCGCAAGCAGGACTAGGTCGTCAGCATAGGCCAGAGCGGACACGAGTTCCATCTCCAACCTATACCCGACCCTCTCCGGCAGGGAAGCCAGGATGAGGTCCATCACCACGTTGAAGAGTATCGGCGACAGAGGGTCCCCTTGACGAACCCCTCGTCCCACTTTTACAGGGCTGCTCATTTCATTGTTCACGGCTAAGGTGGTGGACGCCGTATCGTATAGGTGAGCAATGTAGCCGCAGAACTGTTCGGGCATGCCCCTCAACCTCAGCAATTCGACAAGTGCCTCGTGAGACACTGTGTCAAATGCCTTGGCGAAGTCTAGCACCGCCACGTGACATTCCCGCAGCTTCTTCCTGCTATCCCCAAGCACCGCGTCCAGTACTGCGGAATTCTCCAGCGTACCGTCGGCGCAGATAAATCCGCGCTGTCGTGCATCAGGGGGGCAGCAAGCCAACAGCCTCCGGG
    """

    seq = "".join(seq.split())

    runs = find_homopolymers(seq)

    hist = build_histogram(runs)

    total = build_total_histogram(runs)

    print_histogram(hist)

    print_total_histogram(total)

    print_longest(runs)

    print_long_runs(runs, min_len=2)
    print_summary_table(hist)


if __name__ == "__main__":
    main()