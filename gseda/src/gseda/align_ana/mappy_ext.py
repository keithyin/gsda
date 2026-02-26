import re


def revcomp(seq):
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def parse_cigar(cigar):
    """把 CIGAR 拆成 [(len, op), ...]"""
    num = ""
    for c in cigar:
        if c.isdigit():
            num += c
        else:
            yield int(num), c
            num = ""


def blast_like_alignment(query, ref, aln, line_width=110):
    q_pos = aln.q_st
    if aln.strand == -1:
        q_pos = len(query) - aln.q_en
        query = revcomp(query)
    r_pos = aln.r_st

    q_line = []
    r_line = []
    mid_line = []

    for length, op in parse_cigar(aln.cigar_str):
        if op == "M" or op == "=" or op == "X":
            for i in range(length):
                qb = query[q_pos]
                rb = ref[r_pos]
                q_line.append(qb)
                r_line.append(rb)
                mid_line.append("|" if qb == rb else " ")
                q_pos += 1
                r_pos += 1

        elif op == "I":  # 插入到 reference
            for i in range(length):
                qb = query[q_pos]
                q_line.append(qb)
                r_line.append("-")
                mid_line.append(" ")
                q_pos += 1

        elif op == "D":  # reference 多出来
            for i in range(length):
                rb = ref[r_pos]
                q_line.append("-")
                r_line.append(rb)
                mid_line.append(" ")
                r_pos += 1

        elif op == "S":  # soft clip，通常忽略
            q_pos += length

    # 分段打印（像 blastn）
    q_str = "".join(q_line)
    r_str = "".join(r_line)
    m_str = "".join(mid_line)

    for i in range(0, len(q_str), line_width):
        q_seg = q_str[i:i+line_width]
        m_seg = m_str[i:i+line_width]
        r_seg = r_str[i:i+line_width]

        q_start = aln.q_st + i + 1
        q_end = q_start + len(q_seg) - 1
        r_start = aln.r_st + i + 1
        r_end = r_start + len(r_seg.replace("-", "")) - 1

        print(f"Sbjct  {r_start:<6} {r_seg} {r_end}")
        print(f"              {m_seg}")
        print(f"Query  {q_start:<6} {q_seg} {q_end}")

        print()


def calculate_identity_from_cigar(cigar_string):
    pattern = r'(\d+)([=IDX])'
    matches = re.findall(pattern, cigar_string)

    match_count = 0
    mismatch_count = 0
    total_aligned = 0

    for length_str, operation in matches:
        length = int(length_str)
        if operation == '=':
            match_count += length
            total_aligned += length
        elif operation == 'X':
            mismatch_count += length
            total_aligned += length
        elif operation in ('I', 'D'):
            total_aligned += length

    if total_aligned == 0:
        return 0.0
    return match_count / total_aligned
