

def plp_counter(seqs):
    """Count nucleotide/gap frequency per position across MSA sequences.

    Returns a 2D array of shape [n_positions, 5] — each row is one alignment
    column, the five columns are (A, C, G, T, gap) counts.

    Args:
        seqs: list of MSA sequence strings (same length), e.g.
              ["ACGT-", "AGG--", ...]

    Returns:
        tuple (data, labels):
            data  — list[list[int]], shape [n_positions, 5]
            labels — ['A', 'C', 'G', 'T', 'gap'] for the column mapping
    """
    if not seqs:
        return [], ['A', 'C', 'G', 'T', 'gap']

    n_cols = len(seqs[0])
    assert all(len(s) == n_cols for s in seqs), "All sequences must have the same length"

    labels = ['A', 'C', 'G', 'T', 'gap']
    data: list[list[int]] = []
    for col in range(n_cols):
        counts = [0, 0, 0, 0, 0]
        for s in seqs:
            ch = s[col].upper()
            if ch == '-':
                counts[4] += 1
            elif ch == 'A':
                counts[0] += 1
            elif ch == 'C':
                counts[1] += 1
            elif ch == 'G':
                counts[2] += 1
            elif ch == 'T':
                counts[3] += 1
        data.append(counts)

    return data, labels

def plp_ratio(seqs):
    """Per-position nucleotide/gap ratios — same shape as plp_counter but as floats.

    Each row sums to 1.0 (or close via floating-point).  gap column is the PLP
    deletion rate at that position.

    Args:
        seqs: list of MSA sequence strings (same length)

    Returns:
        tuple (data, labels):
            data  — list[list[float]], shape [n_positions, 5]
            labels — ['A', 'C', 'G', 'T', 'gap']
    """
    if not seqs:
        return [], ['A', 'C', 'G', 'T', 'gap']

    n_cols = len(seqs[0])
    assert all(len(s) == n_cols for s in seqs), "All sequences must have the same length"

    n_seqs = len(seqs)
    labels = ['A', 'C', 'G', 'T', 'gap']
    data: list[list[float]] = []
    for col in range(n_cols):
        counts = [0.0, 0.0, 0.0, 0.0, 0.0]
        for s in seqs:
            ch = s[col].upper()
            if ch == '-':
                counts[4] += 1.0
            elif ch == 'A':
                counts[0] += 1.0
            elif ch == 'C':
                counts[1] += 1.0
            elif ch == 'G':
                counts[2] += 1.0
            elif ch == 'T':
                counts[3] += 1.0
        data.append([c / n_seqs for c in counts])

    return data, labels
    


def main():

    seqs = [
        "GGGAGGTATT-CTAT-GAAG-----A-CGCATAG",
        "GGGAGGTA-T-CTAT-GAAG-----A-CGCATAG",
        "GGGAGGTATT-CTAT-GAAG-----A-CGCATAG",
        "GGGAGG-ATT-CT-T-G-AG-----A-CGCATA-",
        "GGGAGGTATT-CTAT-GAAGGGGGT--CGCATAG",
        "GGGAGGTATT-CTAT-GAAG-----A-GGCAT--",
        "-GGAGGTATT-CTAT-GAAG-----A-CGC-TAG",
        "GGGA-GT-TTCCTATGGAAGGGGGTA-CGCATAG",
        "GGGAGGTATTCCTAT-GAAG-----A-CGCATAG",
        "GGG-GGTATT-CTAT-GAAG-----A-CGCATA-",
        "GGGAGGTATT-CTAT-GAAG-----ACCGCATA-",
        "GGTAGGT-TT-CTAT-GAAG-----A-CGCATAG",
        "GGGAGGTATT-CTAT-GAAG-----ACCGCATAG",
        "GGGAGGT-TT-CTAT-GAAG-----A-CGCATA-",
        "GGGAGGTATT-CTAT-GAAG-----A-CGCATAG",
        "GGGAGGTATT-CTAT-GAAG-----A-CGCATAG",
        "GGGAGGTATT-C-AT-GAAG-------CGCATAG",
        "GGGAGGTATT-CTAT-GAAG-----A-CGCATAG",
        "GGGAGGTATT-CTAT-GAAG-----A-CGCATA-"
    ]

    data, labels = plp_counter(seqs)
    n_rows = len(data)
    n_cols = len(data[0]) if data else 0
    header = "  ".join(labels) + "  total"
    print(f"[{n_rows}, {n_cols}]   {header}")
    for i, row in enumerate(data):
        total = sum(row)
        print(f"{i:>4d}   {'  '.join(str(v) for v in row):<16s}   {total}")
        
    print(plp_counter(seqs=seqs[1:]))
    print(plp_ratio(seqs=seqs[1:])[0])


if __name__ == "__main__":
    main()
