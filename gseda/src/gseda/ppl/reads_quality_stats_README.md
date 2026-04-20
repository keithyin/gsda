# Reads Quality Statistics Tools

This directory contains a suite of tools for evaluating the quality of BAM files, specifically focused on read lengths, alignment metrics, and motif-level precision. These tools typically act as wrappers around the `gsmm2` toolset.

## Overview of Tools

| File | Backend Tool | Primary Purpose | Key Focus |
| :--- | :--- | :--- | :--- |
| `reads_quality_stats_v2.py` | `gsmm2-aligned-metric` | General Quality Reporting | Standard alignment and quality metrics. |
| `reads_quality_stats_v3.py` | `gsmm2-aligned-metric` | General Quality Reporting | Updated version of v2 with refined filtering. |
| `reads_quality_stats_hp.py` | `gsmm2-metric` (`--mode hp-v2`) | High-Precision Motif Analysis | Base-calling precision for specific motifs. |
| `reads_quality_stats_hp_tr.py`| `gsmm2-metric` (`--mode hp-tr-v2`) | High-Precision TR Analysis | Error rates (Match/Mismatch/Indel) for motifs. |

---

## Detailed Functional Analysis

### 1. General Quality Statistics (`v2` & `v3`)
These files provide a comprehensive overview of the sequencing run quality.

- **Core Metrics**:
    - **Read Statistics**: Total reads, total bases, N50, and read length distributions (Avg, Min, P25, P50, P75, Max).
    - **Quality Ratios**: Percentage of bases passing various thresholds ($\ge Q8, \ge Q10, \ge Q15, \ge Q20, \ge Q30$) and the `4xQ20` metric.
    - **Alignment Metrics**: Global and per-read `queryCoverage`, `identity`, `mmRate` (mismatch), `insRate` (insertion), and `delRate` (deletion).
- **Advanced Analysis**:
    - **Segment Analysis**: Analyzes how reads are split into segments and the overlap ratios (`qOvlpRatio`, `rOvlpRatio`).
    - **Gap Analysis**: Identifies and distributes gaps in alignments.
- **Difference (v2 vs v3)**: `v3` is an iterative improvement over `v2`, primarily refining how it handles certain arguments and internal processing.

### 2. High-Precision Motif Analysis (`hp`)
Designed for deep-dives into the accuracy of base calling for specific genomic motifs.

- **Core Logic**:
    - Uses `gsmm2-metric` in `hp-v2` mode to generate detailed motif-level facts.
    - **True vs. Called**: The business logic focuses on comparing the "True Base" (expected) vs the "Called Base" (observed) for each motif.
    - **Ratios**: Computes the distribution of called bases relative to the true base for each motif, allowing for the identification of systematic base-calling biases.

### 3. High-Precision TR Analysis (`hp_tr`)
Specialized for Tandem Repeats (TR) or similar repetitive structures.

- **Core Logic**:
    - Uses `gsmm2-metric` in `hp-tr-v2` mode.
    - **Error Rate Focus**: Unlike `hp.py` which looks at base distributions, `hp_tr.py` focuses on the rates of different event types per motif.
- **Key Metrics**:
    - `eq_rate`: Match rate.
    - `diff_rate`: Mismatch rate.
    - `ins_rate`: Insertion rate.
    - `del_rate`: Deletion rate.

---

## Comparison Summary

| Feature | Standard (`v2`/`v3`) | High-Precision (`hp`) | High-Precision TR (`hp_tr`) |
| :--- | :--- | :--- | :--- |
| **Analysis Level** | Read / Global | Motif | Motif |
| **Primary Goal** | General Quality Assessment | Base-calling Accuracy | Indel/Substitution Rates |
| **Key Output** | Q-stats, N50, Coverage | True vs Called Dist. | Match/Mismatch/Indel Rates |
| **Backend Mode** | `gsmm2-aligned-metric` | `gsmm2-metric --mode hp-v2` | `gsmm2-metric --mode hp-tr-v2` |
