"""Plot mismatch / insertion / deletion error distribution per reference position."""

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

# data & plotting
import polars as pl  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot per-position mismatch / insertion / deletion error distribution."
    )
    parser.add_argument(
        "csv_path", help="Path to the tab-separated locus info CSV")
    parser.add_argument("-o", "--output", default="error_distribution.png",
                        help="Output PNG path (default: error_distribution.png)")
    args = parser.parse_args()

    # read data
    df = pl.read_csv(args.csv_path, separator="\t")
    print(f"Loaded {len(df)} rows, columns: {df.columns}")

    # keep only rows with at least one error
    mask = (pl.col("diff") + pl.col("ins") + pl.col("del")) > 0
    df_err = df.filter(mask).sort("pos").select([
        "pos", "diff", "ins", "del", "depth", "aroundBases",
        "curBase", "curIsHomo",
    ])

    n_rows = len(df_err)
    print(f"Error rows: {n_rows}")
    print(
        f"Total diff: {int(df_err['diff'].sum())}, ins: {int(df_err['ins'].sum())}, del: {int(df_err['del'].sum())}")

    if n_rows == 0:
        print("No error positions found — nothing to plot.")
        return
    
    TOP_N = 50

    # ---- print top50 aroundBases for diff/ins/del ----
    for err_type in ("diff", "ins", "del"):
        label = {"diff": "Mismatch", "ins": "Insertion",
                 "del": "Deletion"}[err_type]
        top_n_dicts = df_err.sort(
            err_type, descending=True).head(TOP_N).to_dicts()
        print(f"\n--- Top {TOP_N} {label} (count > 0) ---")
        for i, row in enumerate(top_n_dicts):
            cnt = int(row[err_type])
            if cnt == 0:
                break
            pos_val = row["pos"]
            ab = row.get("aroundBases", "") or ""
            print(
                f"  #{i+1:2d}  pos={pos_val:>8d}  {err_type:>3s}={cnt:>5d}  aroundBases={ab}")

    # ---- aggregate errors by curBase + curIsHomo, compare homo vs non-homo ----
    # count total positions for each (base, curIsHomo) from the FULL dataframe (incl. zero-error rows)
    base_homo_aggr = df.group_by(["curBase", "curIsHomo"]).len().sort([
        "curBase", "curIsHomo"])
    base_aggr = df.group_by(["curBase"]).len()
    tot_pos = len(df)

    print(base_aggr)
    print(tot_pos)

    base_homo_aggr = base_homo_aggr.to_dicts()

    print(base_homo_aggr)

    agg_df = df_err.group_by(["curBase", "curIsHomo"]).agg([
        pl.col("diff").sum().alias("diff"),
        pl.col("ins").sum().alias("ins"),
        pl.col("del").sum().alias("del"),
        pl.len().alias("n_pos_err"),  # positions with at least one error
    ]).sort(["curBase", "curIsHomo"])

    # build shared lookup table from agg_df + n_pos_full
    print(agg_df)
    comp_rows = []
    for row in agg_df.to_dicts():
        base = row["curBase"]
        homo_or_not = "homo" if row["curIsHomo"] == 1 else "non-homo"
        n_tot = None
        for nr in base_homo_aggr:
            if nr["curBase"] == base and nr["curIsHomo"] == row["curIsHomo"]:
                n_tot = int(nr["len"])
                break
        comp_rows.append({
            "base": base,
            "region": homo_or_not,
            "diff": int(row["diff"]),
            "ins": int(row["ins"]),
            "del": int(row["del"]),
            "errors": int(row["diff"]) + int(row["ins"]) + int(row["del"]),
            "n_pos": n_tot,
            "n_pos_err": int(row["n_pos_err"]),
        })

    # ---- Table 1: error aggregation by curBase + curIsHomo ----
    print("\n=== Error aggregation by curBase + curIsHomo ===")
    t1_rows = [{"Base": r["base"], "Region": r["region"], "diff": r["diff"], "ins": r["ins"],
                "del": r["del"], "total": r["errors"], "n_pos_err": r["n_pos_err"], "n_pos_tot": r["n_pos"]}
               for r in comp_rows]
    print(pl.DataFrame(t1_rows).to_pandas().to_string(index=False))

    # ---- Table 2: homopolymer vs non-homopolymer corrected ratio per base ----
    print("\n=== Homopolymer vs Non-Homopolymer comparison (per base) ===")
    corr_rows = []
    for base in ["A", "C", "G", "T"]:
        homo_row = next(
            (r for r in comp_rows if r["base"] == base and r["region"] == "homo"), None)
        nonhomo_row = next(
            (r for r in comp_rows if r["base"] == base and r["region"] == "non-homo"), None)
        if homo_row is None or nonhomo_row is None:
            continue

        n_homo_tot, n_nonhomo_tot = homo_row["n_pos"], nonhomo_row["n_pos"]
        pos_ratio = (n_homo_tot / (n_nonhomo_tot + n_homo_tot)
                     ) if n_nonhomo_tot > 0 else float("inf")

        # total line
        homo_error_cnt, non_homo_error_cnt = homo_row["errors"], nonhomo_row["errors"]
        # rate_ht = ht / n_homo_tot if n_homo_tot > 0 else 0
        # rate_nt = nt / n_nonhomo_tot if n_nonhomo_tot > 0 else 0
        homo_error_cnt_corrected = int(homo_error_cnt / pos_ratio)
        non_homo_error_cnt_corrected = int(non_homo_error_cnt / (1-pos_ratio))

        corr_rows.append({"Base": base, "Type": "[total]", "h_pos": n_homo_tot, "n_pos": n_nonhomo_tot,
                          "pos_ratio": pos_ratio, "homo_error_cnt": homo_error_cnt, "non_homo_error_cnt": non_homo_error_cnt,
                          "homo_error_cnt_corrected": homo_error_cnt_corrected, "non_homo_error_cnt_corrected": non_homo_error_cnt_corrected,
                          "homoVSnonhomo": homo_error_cnt_corrected / non_homo_error_cnt_corrected
                          })

        # per-type lines
        for err_name in ("diff", "ins", "del"):
            homo_error_cnt = homo_row[err_name]
            non_homo_error_cnt = nonhomo_row[err_name]
            homo_error_cnt_corrected = int(homo_error_cnt / pos_ratio)
            non_homo_error_cnt_corrected = int(non_homo_error_cnt / (1-pos_ratio))

            corr_rows.append({"Base": base, "Type": err_name, "h_pos": None, "n_pos": None,
                              "pos_ratio": None,
                              "homo_error_cnt": homo_error_cnt, "non_homo_error_cnt": non_homo_error_cnt,
                              "homo_error_cnt_corrected": homo_error_cnt_corrected,
                              "non_homo_error_cnt_corrected": non_homo_error_cnt_corrected,
                              "homoVSnonhomo": homo_error_cnt_corrected / non_homo_error_cnt_corrected
                              })

    if corr_rows:
        print(pl.DataFrame(corr_rows).to_pandas().to_string(index=False))

    # ---- Table 3: cross-base homopolymer error comparison (baseline-normalized) ----
    print("\n=== Cross-base homopolymer error comparison ===")

    corr_rows = pl.DataFrame(corr_rows)
    corr_rows = corr_rows.select([pl.col("Base"), pl.col("Type"), (pl.col(
        "homo_error_cnt_corrected")+pl.col("non_homo_error_cnt_corrected")).alias("error_cnt")])

    print(corr_rows.to_pandas().to_string(index=False))

    print(base_aggr)

    #     cross_rows = []
    # # total line per base (sum of diff+ins+del)

    # cross-base comparison: rows = Type, columns = base (raw + corrected per base)
    propensity_scores = {}
    for base in ["A", "C", "G", "T"]:
        ps = base_aggr.filter(pl.col("curBase").eq(pl.lit(base))).select(
            (pl.col("len") / pl.lit(tot_pos)).alias("ps")).item(row=0, column="ps")
        propensity_scores[base] = ps

    # aggregate into one row per Type, with all 4 bases in columns
    by_type: dict[str, dict] = {}
    for r in corr_rows.to_dicts():
        base, err_type = r["Base"], r["Type"]
        raw = int(r["error_cnt"])
        corrected = int(raw / propensity_scores[base])
        if err_type not in by_type:
            by_type[err_type] = {"Type": err_type}
        by_type[err_type][f"raw.{base}"] = raw
        by_type[err_type][f"corr.{base}"] = corrected

    table = pl.DataFrame(list(by_type.values()))
    # sort by custom Type order
    type_order = {"[total]": 0, "diff": 1, "ins": 2, "del": 3}
    table = table.sort(pl.col("Type").map_elements(
        lambda t: type_order.get(t, 9), return_dtype=pl.Int32))
    col_order = [pl.col(col_name) for col_name in ["Type",
                                                   "raw.A", "raw.C", "raw.G", "raw.T", "corr.A", "corr.C", "corr.G", "corr.T"]]
    table = table.select(col_order)

    print(table.to_pandas().to_string(index=False))

    #     ratio = /tot_pos

    #     homo_raw = h_row["errors"]
    #     n_pos_nonhomo = n_row["n_pos"]
    #     n_errors_nonhomo = n_row["errors"]
    #     tot_error = homo_raw + n_errors_nonhomo

    #     corr_cnt = int(non_homo_rate * min_homo_n)

    #     cross_rows.append({
    #         "Base": base, "Type": "total",
    #         "n_pos_homo": h_row["n_pos"], "n_pos_nonhomo": n_row["n_pos"],
    #         "homo_error_cnt": homo_raw,
    #         "non_homo_rate": non_homo_rate,
    #         "corr_error_cnt": corr_cnt,
    #     })

    # ---- plot (3 subplots) ----
    fig, axes = plt.subplots(3, 1, figsize=(16, 9), sharex=True)

    x = df_err["pos"].to_list()
    bar_width = 0.6

    colors = ["#b71c1c", "#2e7d32", "#1565c0"]
    labels = ("Mismatch (diff)", "Insertion (ins)", "Deletion (del)")
    y_cols = ["diff", "ins", "del"]

    for ax, color, label, y_col in zip(axes, colors, labels, y_cols):
        ax.bar(x, df_err[y_col].to_list(),
               width=bar_width, color=color, alpha=1.0)
        ax.set_ylabel(label, fontsize=11)

    axes[-1].set_xlabel("Reference Position", fontsize=11)
    fig.suptitle("Error Distribution per Reference Position",
                 fontsize=13, fontweight="bold", y=0.98)

    # summary box in upper-right corner of first subplot
    total_errors = int(df_err["diff"].sum() +
                       df_err["ins"].sum() + df_err["del"].sum())
    info_text = (
        f"Total errors: {total_errors}\n"
        f"Mismatches / Ins / Del: "
        f"{int(df_err['diff'].sum())} / {int(df_err['ins'].sum())} / {int(df_err['del'].sum())}\n"
        f"Error positions: {n_rows}"
    )
    axes[0].text(0.98, 0.95, info_text, transform=axes[0].transAxes, ha="right", va="top",
                 fontsize=9, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.6))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
