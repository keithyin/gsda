"""将正反向 Sanger 测序读段通过 mappy 比对取公共重叠区合并为单一序列。
仅有单端时直接输出该读段；比对失败或重叠不足时回退为较长读段。"""

import argparse
import logging
import os
import re

import mappy

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)


def read_seq(path):
    """Read a .seq file and return the raw sequence string with whitespace removed."""
    with open(path, "r") as fh:
        return re.sub(r"\s+", "", fh.read())


def primer_strand(primer):
    """Map a primer name to 'F' or 'R'; return None if unrecognised."""
    if primer == "T7":
        return "F"
    if primer == "T7T":
        return "R"
    if primer.startswith("M13F"):
        return "F"
    if primer.startswith("M13R"):
        return "R"
    return None


def discover_samples(input_dir):
    """Scan input_dir for .seq files; return {sample: {'F': seq, 'R': seq}}."""
    samples = {}
    for fname in sorted(os.listdir(input_dir)):
        if not fname.endswith(".seq"):
            continue
        stem = fname[:-4]
        parts = stem.split("_")
        if len(parts) < 4:
            log.warning("Skipping file with unexpected filename: %s", fname)
            continue
        sample = parts[0]
        primer = parts[1]
        strand = primer_strand(primer)
        if strand is None:
            log.warning("Unrecognised primer '%s' in %s – skipping.", primer, fname)
            continue
        seq = read_seq(os.path.join(input_dir, fname))
        samples.setdefault(sample, {})[strand] = seq
    return samples


def parse_cigar(cigar_str):
    """Parse a CIGAR string into a list of (length, op) tuples."""
    return [(int(m.group(1)), m.group(2))
            for m in re.finditer(r"(\d+)([=XIDSMH])", cigar_str)]


def revcomp(seq):
    """Return the reverse complement of a DNA sequence string."""
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def align(f_seq, r_seq):
    """Align r_seq (query) against f_seq (reference) with mappy.

    Returns the first primary hit or None.
    """
    aligner = mappy.Aligner(
        seq=f_seq,
        extra_flags=67108864,
        k=7,
        w=5,
        best_n=10,
        n_threads=1,
        min_cnt=2,
        min_chain_score=4,
    )
    hits = aligner.map(r_seq)
    for hit in hits:
        if hit.is_primary:
            return hit
    return None


def build_columns(f_seq, r_seq, hit, min_overlap):
    """Walk the hit CIGAR and build per-column (F-base, Rrc-base) pairs.

    Returns (columns, overlap_len) on success, or (None, 0) on failure.
    """
    rrc = revcomp(r_seq)
    columns = []
    f_pos = hit.r_st
    r_pos = len(r_seq) - hit.q_en

    for length, op in parse_cigar(hit.cigar_str):
        for _ in range(length):
            if op == "=":
                columns.append((f_seq[f_pos], rrc[r_pos]))
                f_pos += 1
                r_pos += 1
            elif op == "X":
                columns.append((f_seq[f_pos], rrc[r_pos]))
                f_pos += 1
                r_pos += 1
            elif op == "I":
                columns.append(("-", rrc[r_pos]))
                r_pos += 1
            elif op == "D":
                columns.append((f_seq[f_pos], "-"))
                f_pos += 1
            elif op == "S":
                r_pos += 1
            # H (hold) is a no-op

    overlap_len = len(columns)
    if overlap_len < min_overlap:
        return None, 0
    return columns, overlap_len


def is_perfect(col):
    """A column is perfect iff both bases are non-gap and equal."""
    return col[0] != "-" and col[1] != "-" and col[0] == col[1]


def trim_ends(columns):
    """Trim non-perfect-match columns from both ends.

    Returns the kept slice or None if nothing remains.
    """
    left = 0
    right = len(columns)
    while left < right and not is_perfect(columns[left]):
        left += 1
    while right > left and not is_perfect(columns[right - 1]):
        right -= 1
    kept = columns[left:right]
    return kept if kept else None


def hit_identity(hit):
    """Compute identity as #(=) / ( #(=) + #(X) ) from the CIGAR. None if no hit."""
    if hit is None:
        return None
    eq = 0
    x = 0
    for length, op in parse_cigar(hit.cigar_str):
        if op == "=":
            eq += length
        elif op == "X":
            x += length
    total = eq + x
    if total == 0:
        return None
    return eq / total


def merge_pair(f_seq, r_seq, min_overlap):
    """Merge forward and reverse reads into a single sequence.

    Pure function: strings in, dict out. No filesystem access.
    Returns {'sequence': str, 'mode': str, 'identity': float | None}.
    """
    # SINGLE
    if f_seq is None:
        return {"sequence": r_seq, "mode": "SINGLE", "identity": None}
    if r_seq is None:
        return {"sequence": f_seq, "mode": "SINGLE", "identity": None}

    hit = align(f_seq, r_seq)

    # FALLBACK helper
    def _fallback():
        seq = f_seq if len(f_seq) >= len(r_seq) else r_seq
        return {"sequence": seq, "mode": "FALLBACK", "identity": None}

    if hit is None:
        return _fallback()

    columns, _ = build_columns(f_seq, r_seq, hit, min_overlap)
    if columns is None:
        return _fallback()

    kept = trim_ends(columns)
    if kept is None:
        return _fallback()

    final_seq = "".join(c[0] if c[0] != "-" else c[1] for c in kept)
    return {"sequence": final_seq, "mode": "INTERSECT", "identity": hit_identity(hit)}


def wrap(seq, width=60):
    """Wrap a sequence string into a list of lines of at most *width* chars."""
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def write_outputs(outdir, results):
    """Write per-sample .fa files and a combined merged.fa into outdir."""
    os.makedirs(outdir, exist_ok=True)
    merged_parts = []
    for sample in sorted(results):
        r = results[sample]
        lines = [">" + sample] + wrap(r["sequence"])
        content = "\n".join(lines) + "\n"
        with open(os.path.join(outdir, sample + ".fa"), "w") as fh:
            fh.write(content)
        merged_parts.append(content)
    with open(os.path.join(outdir, "merged.fa"), "w") as fh:
        fh.write("".join(merged_parts))


def main():
    parser = argparse.ArgumentParser(description="Merge paired-end Sanger sequencing reads.")
    parser.add_argument("input_dir", help="Directory containing .seq files")
    parser.add_argument("--outdir", default=None, help="Output directory (default: <input_dir>/merged_output)")
    parser.add_argument("--min-overlap", type=int, default=50, help="Minimum overlap length (default 50)")
    args = parser.parse_args()

    outdir = args.outdir or os.path.join(args.input_dir, "merged_output")
    min_overlap = args.min_overlap

    samples = discover_samples(args.input_dir)
    if not samples:
        log.warning("No valid samples found in %s", args.input_dir)
        return

    results = {}
    mode_counts = {}

    for sample in sorted(samples):
        f_seq = samples[sample].get("F")
        r_seq = samples[sample].get("R")
        result = merge_pair(f_seq, r_seq, min_overlap)
        results[sample] = result

        f_len = len(f_seq) if f_seq else 0
        r_len = len(r_seq) if r_seq else 0
        ident_str = "{:.4f}".format(result["identity"]) if result["identity"] is not None else "-"
        log.info("%s | F:%d | R:%d | mode:%s | out:%d | identity:%s",
                 sample, f_len, r_len, result["mode"], len(result["sequence"]), ident_str)

        mode_counts[result["mode"]] = mode_counts.get(result["mode"], 0) + 1

    write_outputs(outdir, results)

    log.info("=== Summary ===")
    log.info("Total samples: %d", len(results))
    for mode in sorted(mode_counts):
        log.info("  %s: %d", mode, mode_counts[mode])
    log.info("Output directory: %s", outdir)


if __name__ == "__main__":
    main()
