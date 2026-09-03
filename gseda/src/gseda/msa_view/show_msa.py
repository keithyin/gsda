"""MSA chromatogram -- generate AB1-like peak plots from MSA alignments.

Each column of the MSA is modelled as four Gaussian peaks (A/C/G/T),
whose heights encode how common each base is at that position across all reads.
The result looks like a real Sanger sequencing trace with mixed signals at variant sites.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

_AB1_COLORS = {"A": "#008000", "C": "#0000FF", "G": "#A52A2A", "T": "#FF0000"}


_DNA_REVERSE = str.maketrans("ACGTacgtNn-.", "TGCAtgcaNn-.")
def _revcomp(seq: str) -> str:
    return seq.translate(_DNA_REVERSE)


def plot_msa_chromatogram(
    msa: Sequence[str],
    *,
    consensus_seq: Optional[str] = None,
    read_labels: Optional[Sequence[str]] = None,
    reverse: bool = False,
    min_coverage: int = 2,
    peak_fwhm: float = 1.2,
    output_path: str | None = None,
    dpi: int = 150,
    fig_width_inches: float | None = None,
    include_read_labels: bool = True,
) -> plt.Figure:
    """Plot an MSA as an AB1-like chromatogram.

    Parameters
    ----------
    msa :
        Iterable of alignment rows (gaps / dots allowed).
    consensus_seq :
        Optional reference shown as text tracker at the bottom.
        Defaults to majority-rule consensus from *msa*.
    read_labels :
        Labels for each row.  Defaults to Read-0, Read-1 ...
    reverse :
        Reverse-complement before plotting (e.g. opposite-strand reads).
    min_coverage :
        Only show peaks where >= this many non-gap bases exist.
    peak_fwhm :
        Full-width at half-max of each Gaussian peak (in MSA-columns).
        Smaller = sharper / more Sanger-like.
    output_path :
        If given, save a PNG copy to this file.
    dpi :
        Output resolution.
    fig_width_inches :
        Override figure width; auto-calculated when None.
    include_read_labels :
        Whether to show read labels on the left side of each row.

    Returns
    -------
    The matplotlib Figure object.
    """

    if not msa:
        raise ValueError("MSA must contain at least one row")

    seqs = [s.upper().replace(".", "-") for s in msa]
    width = len(seqs[0])
    if any(len(s) != width for s in seqs):
        raise ValueError("All MSA rows must have the same length")

    if reverse:
        seqs = [_revcomp(s) for s in seqs]

    n_reads = len(seqs)
    read_labels = list(read_labels or [f"Read-{i}" for i in range(n_reads)])

    # majority-rule consensus
    if consensus_seq is None:
        cs_list = []
        for col in range(width):
            counts: Dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
            for s in seqs:
                b = s[col]
                if b in counts:
                    counts[b] += 1
            cs_list.append(max(counts, key=counts.get))
        consensus_seq = "".join(cs_list)

    n_cols = width

    # coverage per column (for dynamic scaling)
    cov = np.zeros(n_cols, dtype=float)
    for s in seqs:
        for j in range(n_cols):
            if s[j] in _AB1_COLORS:
                cov[j] += 1

    filt_pos = list(np.where(cov >= min_coverage)[0])
    max_amp = float(max(cov[filt_pos])) if filt_pos else 1.0
    peak_sigma = peak_fwhm / np.sqrt(2 * np.log(2))

    # -- figure layout --------------------------------------------------------
    base_h = max(0.7, 3.5 / (n_reads + 1))
    gap_h = max(0.04, 0.25 / max(n_reads + 1, 1))
    n_plot_rows = n_reads + (1 if consensus_seq else 0)
    total_h = n_plot_rows * base_h + n_plot_rows * gap_h

    fig_w = fig_width_inches or max(8, n_cols / 6)
    fig, axes_list = plt.subplots(
        n_plot_rows, 1, figsize=(fig_w, total_h),
        gridspec_kw={"height_ratios": [base_h] * n_plot_rows},
    )
    if n_plot_rows == 1:
        axes_list = [axes_list]

    base_map = {"A": 0, "C": 1, "G": 2, "T": 3}

    # --- read rows ---
    for row_idx in range(n_reads):
        ax = axes_list[row_idx]
        s = seqs[row_idx]
        label = read_labels[row_idx]

        # draw a Gaussian peak for every base at every coverage column
        y_max_local = 0.0
        for col_idx in filt_pos:
            b = s[col_idx]
            if b not in base_map:
                continue
            amp = cov[col_idx] / max_amp
            # evaluate Gaussian along x centered on this column
            vals = amp * np.exp(-0.5 * ((np.arange(n_cols, dtype=float) - col_idx) / peak_sigma) ** 2)
            color = _AB1_COLORS[b]
            ax.plot(np.arange(n_cols), vals, color=color, linewidth=0.6)
            y_max_local = max(y_max_local, amp)

        ax.set_xlim(-0.5, n_cols - 0.5)
        ax.set_ylim(0, max(y_max_local * 1.2, 0.1))
        if include_read_labels:
            ax.set_ylabel(label, fontsize=8, ha="left", va="center", color="#333")
        ax.tick_params(axis="y", length=0)
        ax.set_yticks([])
        ax.grid(True, alpha=0.10, axis="x")
        ax.spines[["top", "right", "bottom"]].set_visible(False)

    # --- consensus tracker row at the bottom ---
    cons_row = n_reads  # index into axes_list
    if consensus_seq:
        cax = axes_list[cons_row]
        y_max_c = 0.0

        for col_idx in filt_pos:
            counts: Dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
            for s in seqs:
                b = s[col_idx]
                if b in counts:
                    counts[b] += 1
            maj = max(counts, key=counts.get)
            amp = counts[maj] / max_amp
            vals = amp * np.exp(-0.5 * ((np.arange(n_cols, dtype=float) - col_idx) / peak_sigma) ** 2)
            cax.plot(np.arange(n_cols), vals, color=_AB1_COLORS[maj], linewidth=0.6)
            y_max_c = max(y_max_c, amp)

        cax.set_xlim(-0.5, n_cols - 0.5)
        cax.set_ylim(0, max(y_max_c * 1.2, 0.1))
        mixed_count = 0
        for col_idx in filt_pos:
            counts: Dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
            for s in seqs:
                b = s[col_idx]
                if b in counts:
                    counts[b] += 1
            maj = max(counts, key=counts.get)
            minority_sum = sum(v for k, v in counts.items() if k != maj) / cov[col_idx]
            is_gap = seqs[0][col_idx] == "-"
            is_mixed = minority_sum >= 0.15 and cov[col_idx] >= min_coverage * 2
            txt_color = _AB1_COLORS[maj] if not is_mixed else f"{_AB1_COLORS[maj]}AA"
            if is_gap:
                txt_color = "#888"
            cax.text(
                col_idx + 0.5, 0, consensus_seq[col_idx],
                ha="center", va="top", fontsize=9, fontweight="bold",
                color=txt_color, family="monospace",
            )
        cax.set_ylabel("consensus overview  (bars on mixed)", fontsize=7.5, color="#666")
        cax.tick_params(axis="y", length=0)
        cax.set_yticks([])
        cax.tick_params(axis="x", labelsize=7)
        cax.grid(True, alpha=0.10, axis="x")
        cax.spines[["top", "right", "bottom"]].set_visible(False)

    # --- shared x-axis on last row only ---
    axes_list[-1].tick_params(axis="x", which="both", direction="inout", length=4, labelsize=7)

    fig.text(0.5, 0.02, "Consensus sequence (majority-rule)", ha="center", fontsize=8, color="#555")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    return fig


# --- CLI demo ----------------------------------------------------------

def main():
    import sys

    msa_data = [
        "-GGGAGGTATT-CTAT-GAAG-----A-CGCATAG",
        "-GGGAGGTA-T-CTATGAAG-----A-CGCATA-",
        "-GGGAGGTATT-CTAT-GAAG-----ACGCATA-",
    ]
    max_len = max(len(s) for s in msa_data)
    msa_data = [s + "-" * (max_len - len(s)) if len(s) < max_len else s for s in msa_data]

    label = " ".join(sys.argv[1:]) if len(sys.argv) > 1 else "demo"
    out = f"/tmp/msa_chromatogram_{label}.png"
    fig = plot_msa_chromatogram(msa_data, output_path=out, dpi=200, fig_width_inches=16)
    plt.close(fig)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
