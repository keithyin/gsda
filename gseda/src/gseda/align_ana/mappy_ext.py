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
    """
    生成类似 BLAST 的双序列比对文本输出。
    参数：
        query, ref: 字符串，0‑based 索引。
        aln: 对象，包含 q_st, q_en, r_st, r_en, strand, cigar_str 等属性。
        line_width: 每行显示的碱基数。
    """
    # 复制坐标，避免修改原始对象
    q_st = aln.q_st
    q_en = aln.q_en
    r_st = aln.r_st
    strand = aln.strand
    cigar_str = aln.cigar_str

    # 负链处理：将 query 替换为反向互补序列，并调整坐标
    if strand == -1:
        query = revcomp(query)                     # 反向互补
        # 原始坐标 q_st, q_en 映射到新序列上
        new_q_st = len(query) - q_en               # 注意：len(query) 现在是原长度
        new_q_en = len(query) - q_st
        q_st, q_en = new_q_st, new_q_en
        # 将比对视为正链
        strand = 1

    # 现在可以统一按正链处理
    q_pos = q_st
    r_pos = r_st

    q_line = []       # 显示的 query 碱基
    r_line = []       # 显示的参考碱基
    mid_line = []     # 匹配/错配符号
    q_coords = []     # 每个位置对应的原始 query 坐标（0‑based）
    r_coords = []     # 每个位置对应的参考坐标（0‑based）

    for length, op in parse_cigar(cigar_str):
        if op in ("M", "=", "X"):          # 匹配或错配
            for _ in range(length):
                qb = query[q_pos]
                rb = ref[r_pos]
                q_line.append(qb)
                r_line.append(rb)
                mid_line.append("|" if qb == rb else " ")
                q_coords.append(q_pos)
                r_coords.append(r_pos)
                q_pos += 1
                r_pos += 1
        elif op == "I":                    # 插入（query 有，参考无）
            for _ in range(length):
                qb = query[q_pos]
                q_line.append(qb)
                r_line.append("-")
                mid_line.append(" ")
                q_coords.append(q_pos)
                r_coords.append(None)
                q_pos += 1
        elif op == "D":                    # 缺失（query 无，参考有）
            for _ in range(length):
                rb = ref[r_pos]
                q_line.append("-")
                r_line.append(rb)
                mid_line.append(" ")
                q_coords.append(None)
                r_coords.append(r_pos)
                r_pos += 1
        elif op == "S":                    # soft clip，仅移动 query 坐标
            q_pos += length
        # 其他操作（H、P 等）可忽略

    # 分块打印
    total_len = len(q_line)
    for i in range(0, total_len, line_width):
        q_seg = q_line[i:i+line_width]
        m_seg = mid_line[i:i+line_width]
        r_seg = r_line[i:i+line_width]
        qc_seg = q_coords[i:i+line_width]
        rc_seg = r_coords[i:i+line_width]

        # 查询段的起止坐标（转换为 1‑based）
        q_start = q_end = None
        for coord in qc_seg:
            if coord is not None:
                if q_start is None:
                    q_start = coord
                q_end = coord
        if q_start is None:
            continue
        q_start += 1
        q_end += 1

        # 参考段的起止坐标
        r_start = r_end = None
        for coord in rc_seg:
            if coord is not None:
                if r_start is None:
                    r_start = coord
                r_end = coord
        if r_start is None:
            continue
        r_start += 1
        r_end += 1
        if aln.strand == -1:
            tmp = len(query) - q_end
            q_end = len(query) - q_start
            q_start = tmp
            q_start, q_end = q_end + 1, q_start + 1
            
        print(f"Sbjct  {r_start:<6} {''.join(r_seg)} {r_end}")
        print(f"              {''.join(m_seg)}")
        print(f"Query  {q_start:<6} {''.join(q_seg)} {q_end}")
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
