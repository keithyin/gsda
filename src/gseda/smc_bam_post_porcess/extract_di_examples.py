#!/usr/bin/env python3
"""Extract example reads demonstrating high-phred `di` bases in homopolymer
context and print them in a clear, scannable format.

Each example shows:
  1. The consensus sequence around the `di` position, with the homopolymer
     underlined
  2. What subreads support (the `di` base)
  3. What `insert_di` would change
"""

import pysam

IN = (
    "/data1/ccs_data/ecoli-data-for-marketing/"
    "20250623_250302Y0002_Run0004_adapter.withdel.smc_all_reads.bam"
)


def parse_di(di_string):
    """Tolerant 5-field parser."""
    if not di_string:
        return []
    out = []
    for seg in str(di_string).split(";"):
        seg = seg.strip()
        if not seg:
            continue
        parts = seg.split(",")
        if len(parts) != 5:
            continue
        try:
            phreq = float(parts[4])
        except (ValueError, IndexError):
            continue
        out.append((int(parts[0]), parts[1].upper(),
                    float(parts[2]), int(parts[3]), phreq))
    return out


def homopolymer_around(seq, pos, base):
    """Find homopolymer run that spans or touches `pos` for `base`.

    Returns (base_type, run_start, run_end, run_length, is_full_span).
    is_full_span means the insertion point sits *within* the homopolymer.
    """
    b = base.upper()
    left_ok = (pos > 0 and seq[pos - 1].upper() == b)
    if left_ok:
        start = pos - 1
        while start > 0 and seq[start - 1].upper() == b:
            start -= 1
        if pos < len(seq) and seq[pos].upper() == b:
            end = pos
            while end + 1 < len(seq) and seq[end + 1].upper() == b:
                end += 1
            return b, start, end, end - start + 1, True
        return b, start, pos, pos - start + 1, False
    right_ok = (pos < len(seq) and seq[pos].upper() == b)
    if right_ok:
        end = pos + 1
        while end < len(seq) and seq[end].upper() == b:
            end += 1
        return b, pos, end, end - pos + 1, False
    return None, -1, -1, 0, False


def display_example(ex):
    """Print a single example in a clear, scannable format."""
    seq = ex["seq"]
    dbase = ex["dbase"]
    hp_base = ex["hp_base"]
    hp_start = ex["hp_start"]
    hp_end = ex["hp_end"]
    hp_len = ex["hp_len"]
    pos = ex["dpos"]
    n_di = ex["n_di"]
    phred = int(ex["phred"])

    left_start = max(0, hp_start - 3)
    right_end = min(len(seq), hp_end + 4)
    left_seq = seq[left_start:hp_start]
    right_seq = seq[hp_end + 1:right_end]
    n_support = int(round(ex["frac"] * ex["depth"]))
    n_oppose = ex["depth"] - n_support

    # Build visual markers
    w = 70

    print()
    print(tilde_bar(w))
    print("  Read    : " + ex["name"][:60])
    print("  rq      : {:.4f}  np   : {}  n_di      : {}".format(ex["rq"], ex["np"], n_di))
    print("  frac    : {:.3f}  depth  : {}  phred     : Q{}".format(ex["frac"], ex["depth"], phred))
    print()

    hp_region = hp_base * hp_len
    full_before = left_seq + hp_region
    marker_start = hp_start - left_start
    print("  Homopolymer '{}' x {} at position {:>4d}:{:>4d}".format(
        hp_base, hp_len, hp_start, hp_end))
    print()
    # Show sequence with underline
    print("  Consensus: [...{}{}{}]...".format(left_seq, hp_region, right_seq))
    print("             {:>{}}{}{}".format(
        "", len(left_seq), "~" * hp_len, " " * (len(right_seq))))

    # Show where di position lands
    di_off = pos - hp_start
    print()
    print("  >>> Subreads want to insert '{}' at offset +{:>2d}".format(dbase, di_off))
    print("      {:>{}}{}{}".format(
        "", di_off + 14, "^", " " * max(0, 14 - hp_len)))
    # The └ symbol with appropriate width
    print("      {:>{}}+-- within homopolymer: {} '{}'".format(
        "", di_off + 13, hp_len, hp_base))

    if ex["is_full"]:
        print()
        print("  * INSIDE THE HOMOPOLYMER:")
        print("    Subreads support {} copies of '{}', consensus only has {}.")
        print("    → SMN underestimated the homopolymer length by 1.")
        print("    (count goes from {} → {} for '{}')".format(
            hp_len, hp_len + 1, hp_base))
    else:
        print()
        print("  * NEIGHBOR-EXTENDING:")
        msg1 = "    The inserted '{}' matches a neighboring consensus base,"
        msg2 = "    meaning it should extend the nearby homopolymer."
        print(msg1.format(dbase))
        print(msg2)

    # Show before / after
    print()
    before_str = left_seq + hp_region + right_seq
    new_hp = hp_base * (hp_len + 1)
    after_str = left_seq + new_hp + right_seq
    print("  Before:  [...{}{}]...".format(left_seq, hp_region))
    print("  After:   [...{}{}]...".format(left_seq, new_hp))
    underline = " " * len(left_seq) + "~" * (hp_len + 1)
    # Fix spacing if right_seq is too long
    if len(after_str) > 27:
        underline += " " * (len(right_seq) - len("..."))

    print()
    print("  Subreads: {}/{} agree, {}/{} disagree".format(
        n_support, ex["depth"], n_oppose, ex["depth"]))
    # phred interpretation
    error_rate = 10 ** (-phred / 10.0)
    print("  → Confidence: Q{} means ~{:.2e} error rate for calling a *different* "
          "base".format(phred, error_rate))


def tilde_bar(width):
    return "~" * width


def main():
    examples = []
    seen = set()

    with pysam.AlignmentFile(IN, "rb", check_sq=False, threads=8) as f:
        for record in f.fetch(until_eof=True):
            if not record.has_tag("np"):
                continue
            if int(record.get_tag("np")) < 3:
                continue
            if not record.has_tag("rq"):
                continue
            if float(record.get_tag("rq")) < 0.9:
                continue
            if not record.has_tag("di"):
                continue

            units = parse_di(record.get_tag("di"))
            if not units:
                continue

            seq = record.query_sequence or ""
            if not seq:
                continue

            for (pos, base, frac, depth, phred) in units:
                if phred <= 20:
                    continue
                hp = homopolymer_around(seq, pos, base)
                hp_base, hp_start, hp_end, hp_len, is_full = hp
                if hp_base is None:
                    continue

                key = record.query_name + "|" + str(pos)
                if key in seen:
                    continue
                seen.add(key)

                examples.append({
                    "name": record.query_name,
                    "rq": float(record.get_tag("rq")),
                    "np": int(record.get_tag("np")),
                    "dpos": pos,
                    "n_di": len(units),
                    "dbase": base,
                    "frac": frac,
                    "depth": depth,
                    "phred": phred,
                    "seq": seq,
                    "hp_base": hp_base,
                    "hp_start": hp_start,
                    "hp_end": hp_end,
                    "hp_len": hp_len,
                    "is_full": is_full,
                })
                if len(examples) >= 6:
                    break
            if len(examples) >= 6:
                break

    return examples


examples = main()

print("=" * 80)
print("HIGH-PHRED (>Q20) `di` BASES -- CONCRETE EXAMPLES")
print("Each shows a homopolymer where SMN consensus missed a copy")
print("that subreads confidently support.")
print("=" * 80)

for ex in examples:
    display_example(ex)

# Write TSV
with open("/tmp/di_examples.tsv", "w") as f:
    header = (
        "example\tread_name\trq\tnp\tn_di\tdi_pos\tdi_base\t"
        "di_phred\tdi_frac\tdi_depth\thomopolymer\thp_length"
    )
    f.write(header + "\n")
    for i, ex in enumerate(examples, 1):
        row = (
            "{i}\t{name}\t{rq:.4f}\t{np}\t{n_di}\t{dpos}\t{dbase}\t"
            "Q{phred}\t{frac:.3f}\t{depth}\t{hp_base}\t{hp_len}".format(
            i=i, name=ex["name"][:40], rq=ex["rq"], np=ex["np"],
            n_di=ex["n_di"], dpos=ex["dpos"], dbase=ex["dbase"],
            phred=int(ex["phred"]), frac=ex["frac"],
            depth=ex["depth"], hp_base=ex["hp_base"], hp_len=ex["hp_len"]))
        f.write(row + "\n")

print()
print("~" * 80)
print("Table -> /tmp/di_examples.tsv ({} examples)".format(len(examples)))
print("~" * 80)
