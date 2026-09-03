#!/usr/bin/env python3
"""Analyze `di` (missing-base) units: compare high-phred (>Q20) vs the rest.

Builds one row per `di` unit directly from an smc consensus BAM (or loads a
precomputed `.npz` with the same schema), then prints per-cohort metrics and the
full Q (phred) distribution as a console histogram, and renders a 2x3 figure
comparing the >Q20 cohort against <=Q20 on phred distribution, base composition,
support fraction, homopolymer context, and the local consensus quality at the
insertion position.

`di` tag format (authoritative, per insert_di.md):
    "pos,base,frac,depth,phreq;..."   (pos may repeat)

The input BAM is NOT an alignment BAM (unmapped smc consensus), so it is
opened with `check_sq=False` and iterated via `fetch(until_eof=True)`.

Row schema (one row per `di` unit):
    phreq  float   the unit's phreq
    base   int8    "ACGT".index(base)  -> A=0,C=1,G=2,T=3
    frac   float   support fraction
    depth  int     support depth
    lq     float   mean consensus qual over seq[pos-2 .. pos+1] (the 4
                   existing bases bracketing the insertion; clamped to ends)
    homo   int8    1 if seq[pos-1]==base or seq[pos]==base (homopolymer ctx)
    ndi    int     number of `di` units in this read (repeated per unit)
    rq     float   whole-read `rq` tag
"""

from __future__ import annotations

import argparse
import array
import pathlib
import sys

import numpy as np
import pysam
from tqdm import tqdm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# --- reuse the tolerant `di` parser from the sibling module --------------------
try:
    from .insert_di import parse_di
except ImportError:  # pragma: no cover - only if run as a standalone script
    from insert_di import parse_di

_BN = ["A", "C", "G", "T"]
_BIDX = {b: i for i, b in enumerate(_BN)}

# --- project palette (validated, see references/palette.md) -------------------
COLOR = "#2a78d6"
CO = "#eb6834"
SUR = "#fcfcfb"
GRID = "#e1e0d9"
INK = "#0b0b0b"


# ---------------------------------------------------------------------------
# Collection
# ---------------------------------------------------------------------------

def collect_di_units(
    input_bam: str,
    min_np: int | None = None,
    min_rq: float | None = None,
    threads: int = 4,
    show_progress: bool = True,
) -> dict[str, np.ndarray]:
    """Scan the smc consensus BAM; return one row (dict of arrays) per `di` unit.

    A record is kept only if its `np` tag is >= min_np and its `rq` tag is
    >= min_rq; an absent tag is treated as failing that filter. min_* = None
    disables the respective filter.
    """
    ph = array.array("f"); ba = array.array("b"); fr = array.array("f")
    de = array.array("i"); lq = array.array("f"); ho = array.array("b")
    nd = array.array("i"); rq = array.array("f")
    n_read = 0

    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, threads=threads) as f:
        iterator = f.fetch(until_eof=True)
        if show_progress:
            iterator = tqdm(iterator, desc=f"collect {pathlib.PurePosixPath(input_bam).name}")

        for r in iterator:
            n_read += 1
            if min_np is not None and (not r.has_tag("np") or int(r.get_tag("np")) < min_np):
                continue
            if min_rq is not None and (not r.has_tag("rq") or float(r.get_tag("rq")) < min_rq):
                continue
            if not r.has_tag("di"):
                continue
            units = parse_di(r.get_tag("di"))
            if not units:
                continue

            seq = r.query_sequence or ""
            rlen = len(seq)
            quals = r.query_qualities
            Qa = np.asarray(quals) if quals is not None else None
            ndi = len(units)
            r_q = float(r.get_tag("rq")) if r.has_tag("rq") else 0.0

            for (pos, base, frac, depth, phreq) in units:
                # local consensus qual: the 4 existing bases bracketing the insertion
                if Qa is not None and rlen:
                    lo = max(0, pos - 2)
                    hi = min(rlen, pos + 2)
                    lv = float(Qa[lo:hi].mean()) if hi > lo else 0.0
                else:
                    lv = 0.0
                hv = 1 if (pos > 0 and seq[pos - 1] == base) \
                          or (pos < rlen and seq[pos] == base) else 0
                ph.append(phreq); ba.append(_BIDX.get(base, 0)); fr.append(frac)
                de.append(depth); lq.append(lv); ho.append(hv)
                nd.append(ndi); rq.append(r_q)

    print(f"records scanned: {n_read:,}; di units collected: {len(ph):,}", file=sys.stderr)
    return {
        "phreq": np.array(ph, "f"), "base": np.array(ba, "i1"),
        "frac": np.array(fr, "f"), "depth": np.array(de, "i"),
        "lq": np.array(lq, "f"), "homo": np.array(ho, "i1"),
        "ndi": np.array(nd, "i"), "rq": np.array(rq, "f"),
    }


def load_di_units(npz_path: str) -> dict[str, np.ndarray]:
    """Load a precomputed `.npz` (same schema) into a dict of arrays."""
    d = np.load(npz_path)
    return {k: d[k] for k in ("phreq", "base", "frac", "depth", "lq", "homo", "ndi", "rq")}


# ---------------------------------------------------------------------------
# Reporting (stdout)
# ---------------------------------------------------------------------------

def print_phred_distribution(phreq: np.ndarray, q_cut: int = 20, bar_w: int = 48) -> None:
    """Print the `di` phred values as an integer-binned console histogram.

    Bins are floor(Q), matching the 1-unit bins of the figure. Bars are linear
    against the peak bin, i.e. proportional to the raw per-bin probability. The
    distribution spans several orders of magnitude (peak ~1e6 at the modal Q,
    the Q45 spike ~1e4), so everything outside the peak region floors to a
    single `#`. The Q20 cohort split is additionally reported on the raw values
    (a floored bin straddling the cut — e.g. an exact 20.0 landing in bin 20 —
    is not a clean divider).
    """
    if phreq.size == 0:
        return
    vals, cnts = np.unique(np.floor(phreq).astype(np.int64), return_counts=True)
    total = int(phreq.size)
    cmax = int(cnts.max())

    print()
    print("di phred (Q) distribution  [bin = floor(Q), bar = count, linear to peak]:")
    print("  {:>6s} {:>10s} {:>9s} {:>9s}".format("Q", "count", "share%", "cum%"))
    cum = 0
    for v, cnt in zip(vals, cnts):
        n = int(cnt)
        cum += n
        bar = "#" * max(1, int(round(bar_w * n / cmax)))
        print("  {:>6d} {:>10,d} {:>9.4f} {:>9.4f}  {}".format(
            v, n, 100 * n / total, 100 * cum / total, bar))
        if v == q_cut - 1:
            print("  {:>6s} {:>10s} {:>9s} {:>9s}  {}".format(
                "-", "", "", "",
                f"Q{q_cut} cut (bin {q_cut} straddles it; exact split below)"))

    n_hi = int((phreq > q_cut).sum())
    n_lo = total - n_hi
    print("  " + "-" * 70)
    print("  exact split on raw Q:  <={cut} {n_lo:,} ({p_lo:.3f}%)   "
          ">{cut} {n_hi:,} ({p_hi:.3f}%)".format(
              cut=q_cut, n_lo=n_lo, n_hi=n_hi,
              p_lo=100 * n_lo / total, p_hi=100 * n_hi / total))
    top = int(vals.max())
    print("  occupied bins: {}/{} (Q{a}..Q{b}); peak bin Q{p} with {n:,}".format(
        len(vals), top + 1, a=int(vals.min()), b=top, p=int(vals[cnts.argmax()]),
        n=cmax))


def report(d: dict[str, np.ndarray]) -> None:
    """Print cohort counts and the high-vs-low metric table."""
    phreq = d["phreq"]; base = d["base"]; frac = d["frac"]; depth = d["depth"]
    lq = d["lq"]; homo = d["homo"]; ndi = d["ndi"]; rq = d["rq"]

    hi = phreq > 20
    lo = ~hi

    print("=" * 72)
    print("total di  units          : {:,}".format(len(phreq)))
    print("  >Q20  (high)           : {:,}  ({:.3f}% of all)".format(
        int(hi.sum()), 100 * hi.mean()))
    print("  <=Q20 (low)            : {:,}  ({:.3f}%)".format(
        int(lo.sum()), 100 * lo.mean()))
    print("  phred  min / median / max: {:.1f} / {:.1f} / {:.1f}".format(
        phreq.min(), np.median(phreq), phreq.max()))

    n_hi = int(hi.sum())
    n_hi_lqhi = int((hi & (lq > 20)).sum())
    n_hi_lqlo = int((hi & (lq <= 20)).sum())
    print("  of hi-phred, local-qual also >Q20: {:,} ({:.1f}%)".format(
        n_hi_lqhi, 100 * n_hi_lqhi / n_hi if n_hi else 0))
    print("  of hi-phred, local-qual   <=Q20   : {:,} ({:.1f}%)".format(
        n_hi_lqlo, 100 * n_hi_lqlo / n_hi if n_hi else 0))
    print()
    print("  {:>26s} {:>12s} {:>12s} {:>8s}".format("metric", ">Q20", "<=Q20", "ratio"))

    def row(name, f):
        a, b = float(f(hi)), float(f(lo))
        r = (a / b) if b else float("nan")
        print("  {:>26s} {:>12.3f} {:>12.3f} {:>8.2f}".format(name, a, b, r))

    row("median depth",        lambda m: np.median(depth[m]))
    row("mean frac",           lambda m: np.mean(frac[m]))
    row("median local-qual",   lambda m: np.median(lq[m]))
    row("share homopol.",      lambda m: np.mean(homo[m]))
    row("mean read rq",        lambda m: np.mean(rq[m]))
    row("avg di per read",     lambda m: np.mean(ndi[m]))

    print()
    print("base comp (share):")
    print("  >Q20 ", {b: 100 * int(round(100 * np.mean(base[hi] == i))) / 100
                     for i, b in enumerate(_BN)})
    print("  <=Q20", {b: 100 * int(round(100 * np.mean(base[lo] == i))) / 100
                     for i, b in enumerate(_BN)})
    print_phred_distribution(phreq)
    print("=" * 72)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def figure(d: dict[str, np.ndarray], out_png: str, dpi: int = 140) -> str:
    """Render the 2x3 high-vs-low comparison figure and save it as PNG."""
    phreq = d["phreq"]; base = d["base"]; frac = d["frac"]; depth = d["depth"]
    lq = d["lq"]; homo = d["homo"]; ndi = d["ndi"]; rq = d["rq"]

    hi = phreq > 20
    lo = ~hi
    n_hi = int(hi.sum())

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.set_facecolor(SUR)

    def sa(ax):
        ax.set_facecolor(SUR)
        ax.set_axisbelow(True)
        ax.grid(True, axis="y", color=GRID, lw=0.8, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(GRID)
        ax.tick_params(colors=INK, labelsize=8)

    # --- (1,1) phred distribution (split color at Q20) -------------------------
    edges = np.arange(0, phreq.max() + 2, 1)
    c = edges[:-1] + 0.5
    cnt, _ = np.histogram(phreq, bins=edges)
    for i in range(len(cnt)):
        if cnt[i] == 0:
            continue
        ccol = CO if c[i] > 20.5 else COLOR
        axes[0, 0].bar(c[i], cnt[i], width=1, color=ccol, alpha=0.9)
    axes[0, 0].axvline(20.5, color=INK, lw=1, ls="--", alpha=0.4)
    axes[0, 0].text(21.5, int(cnt.max() * 0.88), "Q20", color=INK, fontsize=8)
    sa(axes[0, 0])
    axes[0, 0].set_xlabel("di phred score")
    axes[0, 0].set_ylabel("# di units")
    axes[0, 0].set_title("di phred  (blue = >Q20, orange = <=Q20)",
                         fontsize=11, fontweight="bold")

    # --- (1,2) local qual, share-normalized overplot ----------------------------
    edges_lq = np.arange(0, 62, 1)
    c_lq = edges_lq[:-1] + 0.5
    hi_cnt = np.histogram(lq[hi], bins=edges_lq)[0].astype(float)
    lo_cnt = np.histogram(lq[lo], bins=edges_lq)[0].astype(float)
    hi_sh = hi_cnt / (hi_cnt.sum() or 1)
    lo_sh = lo_cnt / (lo_cnt.sum() or 1)
    ax01 = axes[0, 1]
    ax01.plot(c_lq, hi_sh, color=CO, lw=1.5, alpha=0.9, label=">Q20")
    ax01.plot(c_lq, lo_sh, color=COLOR, lw=1.5, alpha=0.8, label="<=Q20")
    ax01.fill_between(c_lq, hi_sh, alpha=0.15, color=CO)
    ax01.fill_between(c_lq, lo_sh, alpha=0.15, color=COLOR)
    ax01.set_xlabel("local consensus qual (avg +/-2 bp)")
    ax01.set_ylabel("normalized share")
    sa(ax01); ax01.legend(fontsize=9, frameon=False)
    ax01.set_title("local quality at `di` position (share norm.)",
                   fontsize=11, fontweight="bold", color=INK)

    # --- (1,3) frac (support fraction) ------------------------------------------
    fe = np.linspace(0, 1, 41)
    hx, _ = np.histogram(frac[hi], bins=fe); hx = hx.astype(float); hx /= (hx.sum() or 1)
    lx, _ = np.histogram(frac[lo], bins=fe); lx = lx.astype(float); lx /= (lx.sum() or 1)
    ax02 = axes[0, 2]
    ax02.plot(fe[:-1] + 0.0125, hx, color=CO, lw=1.5, label=">Q20")
    ax02.plot(fe[:-1] + 0.0125, lx, color=COLOR, lw=1.5, label="<=Q20")
    ax02.set_xlabel("frac (support fraction)")
    ax02.set_ylabel("normalized share")
    sa(ax02); ax02.legend(fontsize=9, frameon=False)
    ax02.set_title("support fraction", fontsize=11, fontweight="bold", color=INK)

    # --- (2,1) key metrics text table ------------------------------------------
    ax = axes[1, 0]
    ax.axis("off")
    sa(ax)
    lines = [
        "Characteristics of high-phred (>Q20)",
        "di bases",
        "",
        "share of all di       {:,} {:.3f}%".format(n_hi, 100 * hi.mean()),
        "median depth          {:.1f}".format(np.median(depth[hi])),
        "frac (support)        {:.4f}  vs  {:.4f}  (<=Q20)".format(
            np.mean(frac[hi]), np.mean(frac[lo])),
        "homopolymer context   {:.1f}%  vs  {:.1f}%  (<=Q20)".format(
            100 * np.mean(homo[hi]), 100 * np.mean(homo[lo])),
        "read rq               {:.4f}  vs  {:.4f}".format(
            np.mean(rq[hi]), np.mean(rq[lo])),
        "avg di per read       {:.1f}  vs  {:.1f}".format(
            np.mean(ndi[hi]), np.mean(ndi[lo])),
        "base A C G T         {} {} {} {}".format(*(
            str(int(round(100 * np.mean(base[hi] == i)))) + "%" for i in range(4))),
    ]
    ax.text(0.02, 0.5, "\n".join(lines), transform=ax.transAxes, fontsize=9,
            color=INK, family="monospace", va="center")
    ax.set_title("Key metrics (>Q20 vs <=Q20)", fontsize=11, fontweight="bold", color=INK)

    # --- (2,2) local qual count histogram (log) ---------------------------------
    hi_cnt2 = np.histogram(lq[hi], bins=edges_lq)[0]
    lo_cnt2 = np.histogram(lq[lo], bins=edges_lq)[0]
    mask = (hi_cnt2 > 0) | (lo_cnt2 > 0)
    c2 = c_lq[mask]
    ax21 = axes[1, 1]
    sa(ax21)
    ax21.bar(c2, hi_cnt2[mask], width=0.5, color=CO, alpha=0.85)
    ax21.bar(c2, lo_cnt2[mask], width=0.5, color=COLOR, alpha=0.6)
    ax21.set_yscale("log")
    ax21.set_xlabel("local consensus qual (avg +/-2 bp)")
    ax21.set_ylabel("# di units")
    ax21.set_title("local qual count (cohort-specific)", fontsize=11, fontweight="bold",
                   color=INK)
    ax21.text(0.02, 0.95,
              "median  >Q20 = {:.0f}\n     <=Q20 = {:.0f}".format(
                  np.median(lq[hi]), np.median(lq[lo])),
              transform=ax21.transAxes, fontsize=9, color=INK, va="top",
              bbox=dict(boxstyle="round", fc=SUR, ec=GRID))

    # --- (2,3) empty ------------------------------------------------------------
    axes[1, 2].axis("off"); sa(axes[1, 2])

    fig.suptitle("High-phred (>Q20) `di` bases -- what's special about them?",
                 fontsize=13, fontweight="bold", color=INK)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor=SUR)
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _default_output(input_path: str) -> str:
    """`<input dir>/<input stem>.di_q_dist.png` — next to the input file.

    Same convention as the sibling tools (`di_q45_dist.py`, `di_vs_seq_qual.py`):
    the input name minus its extension is the prefix, `di_q_dist` the infix.
    Takes the BAM or the `.npz` path, whichever was used as input.
    """
    p = pathlib.PurePosixPath(input_path)
    stem = p.name
    for ext in (".bam", ".npz"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    return str(p.with_name(f"{stem}.di_q_dist.png"))


def main_cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare >Q20 vs <=Q20 `di` units from an smc consensus BAM "
                    "(or a precomputed .npz): print cohort metrics + save a 2x3 figure."
    )
    parser.add_argument(
        "input_bam", nargs="?", default=None,
        help="Input smc consensus BAM (unmapped). Provide this OR --npz.",
    )
    parser.add_argument(
        "--npz", type=str, default=None,
        help="Load a precomputed di_units .npz instead of building from a BAM.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output PNG path (default: alongside the input as "
             "<input stem>.di_q_dist.png).",
    )
    parser.add_argument(
        "--min-np", type=int, default=3,
        help="Keep records with np >= --min-np (default: 3).",
    )
    parser.add_argument(
        "--min-rq", type=float, default=None,
        help="Keep records with rq >= --min-rq (e.g. 0.9). Default: no rq filter.",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for BAM I/O (default: 4).",
    )
    parser.add_argument(
        "--dpi", type=int, default=140,
        help="Figure DPI (default: 140).",
    )
    parser.add_argument(
        "--no-progress", action="store_true",
        help="Disable the tqdm progress bar.",
    )

    args = parser.parse_args(argv)

    if args.npz:
        d = load_di_units(args.npz)
    elif args.input_bam:
        d = collect_di_units(
            args.input_bam,
            min_np=args.min_np,
            min_rq=args.min_rq,
            threads=args.threads,
            show_progress=not args.no_progress,
        )
    else:
        parser.error("provide input_bam or --npz")
        return 1

    if len(d["phreq"]) == 0:
        print("No `di` units found — nothing to analyze.", file=sys.stderr)
        return 1

    # exactly one of the two is set here (the else-branch above errors otherwise)
    out = args.output or _default_output(args.npz or args.input_bam)

    report(d)
    saved = figure(d, out, dpi=args.dpi)
    print("\nfigure ->", saved)
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
