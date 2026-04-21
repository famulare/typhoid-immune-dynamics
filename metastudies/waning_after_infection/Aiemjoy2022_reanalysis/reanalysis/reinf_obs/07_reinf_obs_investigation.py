"""
07: Investigate the reinf_obs column by applying the Aiemjoy reinfection
    detection algorithm ourselves and comparing to the CSV flag.

Aiemjoy definition (paper p. e580, appendix 4 p. 2):
  Suspected reinfection = ≥3-fold increase in ≥2 antigen-isotype combinations
  at visits ≥3 months from fever onset, UNLESS the absolute difference < 1 EU.

We apply this to the raw wide-format CSV (all 7 antigen-isotypes), then:
  1. Compare our classification to reinf_obs
  2. Facet Vi IgG trajectories by our reinfection classification × n_obs
  3. Facet Vi IgG trajectories by reinf_obs × n_obs (moved from 04)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/reinf_obs")

ANTIGENS = ["HlyE_IgA", "HlyE_IgG", "LPS_IgA", "LPS_IgG", "MP_IgA", "MP_IgG", "Vi_IgG"]
MIN_FOLD = 3.0
MIN_MONTHS_DAYS = 90  # ≥3 months
MIN_ABS_DIFF = 1.0


def apply_reinfection_algo(raw):
    """Apply Aiemjoy reinfection detection to each subject.

    For each subject, check all consecutive visit pairs where the later visit
    is ≥90 days from fever onset. Count how many antigen-isotypes show a
    ≥3-fold increase with absolute difference ≥1 EU. If ≥2 antigens qualify
    at any visit pair, flag as suspected reinfection.

    Returns a DataFrame with index_id, reinf_algo (bool), n_ag_triggered,
    and trigger_visit (first visit number where reinfection detected).
    """
    results = []
    for idx, row in raw.iterrows():
        sid = row["index_id"]

        # Gather (visit_number, time, {antigen: value}) for all visits with time data
        visits = []
        for v in range(1, 8):
            t = row.get(f"TimeInDays_visit{v}", np.nan)
            if pd.isna(t):
                continue
            vals = {}
            for ag in ANTIGENS:
                val = row.get(f"{ag}_visit{v}", np.nan)
                if pd.notna(val):
                    vals[ag] = val
            visits.append((v, t, vals))

        # Check all pairs: earlier visit vs later visit ≥90 days from fever onset
        reinf_detected = False
        best_n_ag = 0
        trigger_visit = None

        for i in range(len(visits)):
            for j in range(i + 1, len(visits)):
                v_early, t_early, vals_early = visits[i]
                v_late, t_late, vals_late = visits[j]

                # Later visit must be ≥3 months from fever onset
                if t_late < MIN_MONTHS_DAYS:
                    continue

                # Count antigen-isotypes with ≥3-fold increase AND abs diff ≥1
                n_ag_triggered = 0
                for ag in ANTIGENS:
                    if ag not in vals_early or ag not in vals_late:
                        continue
                    val_early = vals_early[ag]
                    val_late = vals_late[ag]
                    if val_early <= 0:
                        continue
                    fold = val_late / val_early
                    abs_diff = val_late - val_early
                    if fold >= MIN_FOLD and abs_diff >= MIN_ABS_DIFF:
                        n_ag_triggered += 1

                if n_ag_triggered >= 2:
                    reinf_detected = True
                    if n_ag_triggered > best_n_ag:
                        best_n_ag = n_ag_triggered
                        trigger_visit = v_late

        results.append({
            "index_id": sid,
            "reinf_algo": reinf_detected,
            "n_ag_triggered": best_n_ag,
            "trigger_visit": trigger_visit,
        })

    return pd.DataFrame(results)


def compute_vi_slopes(long):
    """Compute first-to-last log10 slope for each longitudinal Vi subject."""
    slopes = []
    for sid, grp in long.groupby("index_id"):
        if not grp["has_longitudinal"].iloc[0]:
            continue
        grp = grp.sort_values("days_since_fever_onset")
        t0, eu0 = grp.iloc[0]["days_since_fever_onset"], grp.iloc[0]["vi_igg_eu"]
        t1, eu1 = grp.iloc[-1]["days_since_fever_onset"], grp.iloc[-1]["vi_igg_eu"]
        dt = t1 - t0
        if dt <= 0 or eu0 <= 0 or eu1 <= 0:
            continue
        slopes.append({
            "index_id": sid,
            "reinf_obs": grp.iloc[0]["reinf_obs"],
            "n_obs": len(grp),
            "duration_days": dt,
            "eu_start": eu0,
            "eu_end": eu1,
            "fold_change": eu1 / eu0,
            "direction": "declining" if eu1 < eu0 else "rising",
        })
    return pd.DataFrame(slopes)


def plot_faceted_trajectories(long, slopes, row_col, row_labels, row_filter_fn,
                              npts_vals, dir_colors, title, filename):
    """Generic faceted trajectory plot: rows × n_obs columns, colored by direction."""
    n_cols = len(npts_vals)
    n_rows = len(row_labels)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows),
                             sharex=True, sharey=True)
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    for col, npts in enumerate(npts_vals):
        for row, rlabel in enumerate(row_labels):
            ax = axes[row, col]
            sub_ids = row_filter_fn(slopes, npts, rlabel)
            sub_slopes = slopes[slopes["index_id"].isin(sub_ids)]
            sub_data = long[long["index_id"].isin(sub_ids)]

            for sid, grp in sub_data.groupby("index_id"):
                grp = grp.sort_values("days_since_fever_onset")
                direction = sub_slopes[sub_slopes["index_id"] == sid]["direction"].iloc[0]
                ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                        color=dir_colors[direction], linewidth=1.0, alpha=0.5,
                        marker="o", markersize=3)

            ax.set_yscale("log")
            ax.set_xlim(-10, 1200)
            ax.set_ylim(10, 5000)
            ax.grid(True, alpha=0.3)

            if row == 0:
                total = len(slopes[slopes["n_obs"] == npts])
                ax.set_title(f"{npts} observations (n={total})", fontsize=11)
            if col == 0:
                ax.set_ylabel(f"{rlabel} (n={len(sub_ids)})\nVi IgG (EU)")
            else:
                ax.text(0.02, 0.95, f"n={len(sub_ids)}", transform=ax.transAxes,
                        fontsize=9, va="top", fontweight="bold")
            if row == n_rows - 1:
                ax.set_xlabel("Days since fever onset")

    fig.suptitle(title, fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {filename}")


def main():
    raw = pd.read_csv(DATA_DIR / "TypoidCaseData_github_09.30.21.csv")
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    slopes = compute_vi_slopes(long)

    print("=" * 60)
    print("REINFECTION ALGORITHM — APPLYING AIEMJOY DEFINITION")
    print("=" * 60)

    # Apply algorithm to ALL subjects in CSV
    reinf = apply_reinfection_algo(raw)
    print(f"Total subjects: {len(reinf)}")
    print(f"Suspected reinfection (our algo): {reinf['reinf_algo'].sum()}")
    print(f"Paper reports: 37")
    print()

    # Cross-tab with reinf_obs
    merged = reinf.merge(raw[["index_id", "reinf_obs"]], on="index_id")
    ct = pd.crosstab(
        merged["reinf_obs"].fillna("NaN"),
        merged["reinf_algo"].map({True: "algo=reinfected", False: "algo=clean"}),
        margins=True,
    )
    print("Cross-tab: reinf_obs vs our algorithm")
    print(ct)
    print()

    # Focus on Vi IgG subjects
    vi_ids = set(long["index_id"].unique())
    reinf_vi = reinf[reinf["index_id"].isin(vi_ids)]
    print(f"Vi IgG subjects: {len(reinf_vi)}")
    print(f"  algo reinfected: {reinf_vi['reinf_algo'].sum()}")
    print()

    # Merge algo result into slopes
    slopes = slopes.merge(reinf[["index_id", "reinf_algo", "n_ag_triggered"]], on="index_id")

    # Print summary by algo status
    for status in [True, False]:
        sub = slopes[slopes["reinf_algo"] == status]
        label = "Algo: reinfected" if status else "Algo: clean"
        n_dec = (sub["direction"] == "declining").sum()
        n_ris = (sub["direction"] == "rising").sum()
        print(f"{label} (n={len(sub)}): {n_dec} declining, {n_ris} rising, "
              f"median FC={sub['fold_change'].median():.3f}")
    print()

    # Cross-tab reinf_obs vs algo for longitudinal Vi subjects only
    ct2 = pd.crosstab(
        slopes["reinf_obs"].fillna("NaN"),
        slopes["reinf_algo"].map({True: "algo=reinfected", False: "algo=clean"}),
        margins=True,
    )
    print("Cross-tab (longitudinal Vi only): reinf_obs vs algo")
    print(ct2)
    print()

    # =====================================================================
    # Setup
    # =====================================================================
    npts_vals = sorted(slopes["n_obs"].unique())
    dir_colors = {"declining": "#d62728", "rising": "#2ca02c"}

    # =====================================================================
    # Figure 1: Faceted by reinf_obs (moved from 04)
    # =====================================================================
    def reinf_obs_filter(sl, npts, rlabel):
        if rlabel == "reinf_obs = 1":
            return sl[(sl["n_obs"] == npts) & (sl["reinf_obs"] == 1)]["index_id"]
        else:
            return sl[(sl["n_obs"] == npts) & (sl["reinf_obs"].isna())]["index_id"]

    plot_faceted_trajectories(
        long, slopes,
        row_col="reinf_obs",
        row_labels=["reinf_obs = 1", "reinf_obs = NaN"],
        row_filter_fn=reinf_obs_filter,
        npts_vals=npts_vals,
        dir_colors=dir_colors,
        title="Vi IgG trajectories: n_obs × reinf_obs\n"
              "(reinf_obs=1 ≈ 'evaluable for reinfection screening', not 'reinfected')",
        filename="07_faceted_by_reinf_obs.png",
    )

    # =====================================================================
    # Figure 2: Faceted by our reinfection algorithm
    # =====================================================================
    def algo_filter(sl, npts, rlabel):
        if rlabel == "Algo: suspected reinfection":
            return sl[(sl["n_obs"] == npts) & (sl["reinf_algo"])]["index_id"]
        else:
            return sl[(sl["n_obs"] == npts) & (~sl["reinf_algo"])]["index_id"]

    plot_faceted_trajectories(
        long, slopes,
        row_col="reinf_algo",
        row_labels=["Algo: suspected reinfection", "Algo: no reinfection"],
        row_filter_fn=algo_filter,
        npts_vals=npts_vals,
        dir_colors=dir_colors,
        title="Vi IgG trajectories: n_obs × Aiemjoy reinfection algorithm\n"
              "(≥3-fold rise in ≥2 antigen-isotypes at ≥3 months, abs diff ≥1 EU)",
        filename="07_faceted_by_reinf_algo.png",
    )

    # =====================================================================
    # Figure 3: Combined 2×2 — reinf_obs × algo, ignoring n_obs
    # =====================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
    combos = [
        (0, 0, "reinf_obs=1, algo=clean"),
        (0, 1, "reinf_obs=1, algo=reinfected"),
        (1, 0, "reinf_obs=NaN, algo=clean"),
        (1, 1, "reinf_obs=NaN, algo=reinfected"),
    ]
    for row, col, label in combos:
        ax = axes[row, col]
        if "obs=1" in label:
            mask = slopes["reinf_obs"] == 1
        else:
            mask = slopes["reinf_obs"].isna()
        if "algo=reinfected" in label:
            mask = mask & slopes["reinf_algo"]
        else:
            mask = mask & ~slopes["reinf_algo"]

        sub_ids = slopes[mask]["index_id"]
        sub_slopes = slopes[slopes["index_id"].isin(sub_ids)]
        sub_data = long[long["index_id"].isin(sub_ids)]

        for sid, grp in sub_data.groupby("index_id"):
            grp = grp.sort_values("days_since_fever_onset")
            direction = sub_slopes[sub_slopes["index_id"] == sid]["direction"].iloc[0]
            ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                    color=dir_colors[direction], linewidth=1.0, alpha=0.5,
                    marker="o", markersize=3)

        ax.set_yscale("log")
        ax.set_xlim(-10, 1200)
        ax.set_ylim(10, 5000)
        ax.grid(True, alpha=0.3)
        n_dec = (sub_slopes["direction"] == "declining").sum()
        n_ris = (sub_slopes["direction"] == "rising").sum()
        ax.set_title(f"{label}\n(n={len(sub_ids)}: {n_dec} dec, {n_ris} ris)", fontsize=10)
        if col == 0:
            ax.set_ylabel("Vi IgG (EU)")
        if row == 1:
            ax.set_xlabel("Days since fever onset")

    fig.suptitle("Vi IgG trajectories: reinf_obs × Aiemjoy reinfection algorithm",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "07_reinf_obs_vs_algo.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 07_reinf_obs_vs_algo.png")


if __name__ == "__main__":
    main()
