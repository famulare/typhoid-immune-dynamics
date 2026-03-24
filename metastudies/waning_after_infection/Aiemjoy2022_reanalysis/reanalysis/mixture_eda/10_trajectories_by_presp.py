"""
10: Vi IgG trajectories stratified by Step 4 P(responder) quintile.

Lines colored by fold-change direction (red=declining, green=rising).
Dot shape: circle = not flagged by Aiemjoy reinfection algo,
           star = flagged as suspected reinfection (≥3× in ≥2 Ag at ≥3mo).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module
m08 = import_module("08_fold_change_mixture")

sys.path.insert(0, str(Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/reinf_obs")))
m07 = import_module("07_reinf_obs_investigation")

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")


def main():
    # Load and fit Step 4
    df = m08.load_fold_change_data()
    x = df["log2_fc"].values
    t0 = df["t_start"].values
    t1 = df["t_end"].values

    p1 = m08.fit_two_gaussian(x)
    p4 = m08.fit_teunis_untrunc(x, t0, t1, p1)  # Working model (Step 4a, untruncated)

    df["p_resp"] = 1 - p4["p_noise"]
    df["fc"] = 2 ** df["log2_fc"]
    df["direction"] = np.where(df["log2_fc"] >= 0, "rising", "declining")

    # Quintiles of P(responder)
    df["p_resp_quintile"] = pd.qcut(df["p_resp"], 5, labels=False) + 1
    quintile_bounds = df.groupby("p_resp_quintile")["p_resp"].agg(["min", "max"])

    # Aiemjoy reinfection algo
    raw = pd.read_csv(DATA_DIR / "TypoidCaseData_github_09.30.21.csv")
    reinf = m07.apply_reinfection_algo(raw)
    df = df.merge(reinf[["index_id", "reinf_algo"]], on="index_id", how="left")

    # Load long-format for trajectories
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")

    # Classify: ≥2 samples within 200 days of fever onset (acute follow-up)
    acute_ids = set()
    for sid, grp in long[long["has_longitudinal"]].groupby("index_id"):
        n_within_200 = (grp["days_since_fever_onset"] <= 200).sum()
        if n_within_200 >= 2:
            acute_ids.add(sid)
    df["acute_followup"] = df["index_id"].isin(acute_ids)

    print(f"Subjects: {len(df)}")
    print(f"  ≥2 samples within 200d: {df['acute_followup'].sum()}")
    print(f"  Only long-gap samples:  {(~df['acute_followup']).sum()}")
    print(f"\nP(resp) quintile bounds:")
    for q, row in quintile_bounds.iterrows():
        n = (df["p_resp_quintile"] == q).sum()
        n_reinf = df[df["p_resp_quintile"] == q]["reinf_algo"].sum()
        print(f"  Q{q}: P ∈ [{row['min']:.3f}, {row['max']:.3f}], n={n}, reinf_algo={n_reinf}")

    from matplotlib.lines import Line2D
    dir_colors = {"declining": "#d62728", "rising": "#2ca02c"}

    # =====================================================================
    # Figure 1: 1×5 quintile plot (original)
    # =====================================================================
    fig, axes = plt.subplots(1, 5, figsize=(25, 5), sharey=True, sharex=True)

    for qi, ax in enumerate(axes, 1):
        sub_ids = df[df["p_resp_quintile"] == qi]
        bounds = quintile_bounds.loc[qi]
        n_dec = (sub_ids["direction"] == "declining").sum()
        n_ris = (sub_ids["direction"] == "rising").sum()
        n_reinf = sub_ids["reinf_algo"].sum()

        for _, row in sub_ids.iterrows():
            sid = row["index_id"]
            grp = long[long["index_id"] == sid].sort_values("days_since_fever_onset")
            color = dir_colors[row["direction"]]
            is_reinf = row["reinf_algo"]
            ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                    color=color, linewidth=0.8, alpha=0.5)
            marker = "*" if is_reinf else "o"
            size = 60 if is_reinf else 15
            ax.scatter(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                       color=color, marker=marker, s=size, zorder=5 if is_reinf else 4,
                       alpha=0.8 if is_reinf else 0.5,
                       edgecolors="black", linewidths=0.3 if is_reinf else 0.2)

        ax.set_yscale("log"); ax.set_ylim(10, 5000); ax.set_xlim(-10, 1200)
        ax.grid(True, alpha=0.3); ax.set_xlabel("Days since fever onset")
        ax.set_title(f"Q{qi}: P ∈ [{bounds['min']:.2f}, {bounds['max']:.2f}]\n"
                     f"n={len(sub_ids)} ({n_dec}↓ {n_ris}↑, {n_reinf}★)", fontsize=9)
        if qi == 1:
            ax.set_ylabel("Vi IgG (EU)")

    legend_elements = [
        Line2D([0], [0], color="#d62728", lw=1.5, label="Declining (FC<1)"),
        Line2D([0], [0], color="#2ca02c", lw=1.5, label="Rising (FC≥1)"),
        Line2D([0], [0], marker="o", color="gray", lw=0, markersize=5, label="No reinf flag"),
        Line2D([0], [0], marker="*", color="gray", lw=0, markersize=10, label="Aiemjoy reinf algo"),
    ]
    axes[-1].legend(handles=legend_elements, fontsize=7, loc="upper right")
    fig.suptitle("Vi IgG trajectories by P(responder) quintile\n"
                 "(red=declining, green=rising, ★=Aiemjoy reinfection algo)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "10_trajectories_by_presp.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 10_trajectories_by_presp.png")

    # =====================================================================
    # Figure 2: 2×5 — rows = acute follow-up vs not, columns = quintile
    # =====================================================================
    row_labels = [("≥2 samples ≤200d", True), ("Long-gap only", False)]
    fig, axes = plt.subplots(2, 5, figsize=(25, 10), sharey=True, sharex=True)

    for col, qi in enumerate(range(1, 6)):
        bounds = quintile_bounds.loc[qi]
        for row_idx, (row_label, is_acute) in enumerate(row_labels):
            ax = axes[row_idx, col]
            sub = df[(df["p_resp_quintile"] == qi) & (df["acute_followup"] == is_acute)]
            n_dec = (sub["direction"] == "declining").sum()
            n_ris = (sub["direction"] == "rising").sum()
            n_reinf = sub["reinf_algo"].sum()

            for _, r in sub.iterrows():
                sid = r["index_id"]
                grp = long[long["index_id"] == sid].sort_values("days_since_fever_onset")
                color = dir_colors[r["direction"]]
                is_reinf = r["reinf_algo"]
                ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                        color=color, linewidth=0.8, alpha=0.5)
                marker = "*" if is_reinf else "o"
                size = 60 if is_reinf else 15
                ax.scatter(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                           color=color, marker=marker, s=size,
                           zorder=5 if is_reinf else 4,
                           alpha=0.8 if is_reinf else 0.5,
                           edgecolors="black", linewidths=0.3 if is_reinf else 0.2)

            ax.set_yscale("log"); ax.set_ylim(10, 5000); ax.set_xlim(-10, 1200)
            ax.grid(True, alpha=0.3)

            if row_idx == 0:
                ax.set_title(f"Q{qi}: P ∈ [{bounds['min']:.2f}, {bounds['max']:.2f}]",
                             fontsize=9)
            if col == 0:
                ax.set_ylabel(f"{row_label}\n(n={len(sub)})\nVi IgG (EU)")
            else:
                ax.text(0.02, 0.95, f"n={len(sub)}\n{n_dec}↓ {n_ris}↑ {n_reinf}★",
                        transform=ax.transAxes, fontsize=8, va="top",
                        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
            if row_idx == 1:
                ax.set_xlabel("Days since fever onset")

    axes[-1, -1].legend(handles=legend_elements, fontsize=7, loc="upper right")
    fig.suptitle("Vi IgG trajectories: acute follow-up × P(responder) quintile\n"
                 "(top: ≥2 samples within 200d, bottom: long-gap only)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "10_trajectories_acute_vs_long.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 10_trajectories_acute_vs_long.png")


    # =====================================================================
    # Figure 3: 3-4 point subjects only, stratified by first-two-point gap
    # Rows: first two points within 180d vs not
    # Color: P(resp) quintile
    # Shape: circle vs star (Aiemjoy reinf algo)
    # =====================================================================

    # Get observation count per subject from the full long data
    nobs = long[long["has_longitudinal"]].groupby("index_id").size()
    multipoint_ids = set(nobs[nobs >= 3].index)

    # First-two-point gap for these subjects
    gap_map = {}
    for sid in multipoint_ids:
        grp = long[long["index_id"] == sid].sort_values("days_since_fever_onset")
        gap_map[sid] = grp.iloc[1]["days_since_fever_onset"] - grp.iloc[0]["days_since_fever_onset"]

    df_mp = df[df["index_id"].isin(multipoint_ids)].copy()
    df_mp["first_gap"] = df_mp["index_id"].map(gap_map)
    df_mp["gap_le_180"] = df_mp["first_gap"] <= 180
    df_mp["n_obs"] = df_mp["index_id"].map(nobs)

    print(f"\n3-4 point subjects in fold-change data: {len(df_mp)}")
    print(f"  First-two gap ≤180d: {df_mp['gap_le_180'].sum()}")
    print(f"  First-two gap >180d: {(~df_mp['gap_le_180']).sum()}")

    # Quintile colors
    quintile_cmap = plt.cm.viridis
    q_colors = {q: quintile_cmap(q / 5.0) for q in range(1, 6)}

    row_defs = [("First two points ≤180d apart", True),
                ("First two points >180d apart", False)]
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True, sharey=True)

    for row_idx, (row_label, is_short_gap) in enumerate(row_defs):
        ax = axes[row_idx]
        sub = df_mp[df_mp["gap_le_180"] == is_short_gap]

        for _, r in sub.iterrows():
            sid = r["index_id"]
            grp = long[long["index_id"] == sid].sort_values("days_since_fever_onset")
            qi = r["p_resp_quintile"]
            color = q_colors[qi]
            is_reinf = r["reinf_algo"]

            ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                    color=color, linewidth=1.0, alpha=0.6)
            marker = "*" if is_reinf else "o"
            size = 80 if is_reinf else 20
            ax.scatter(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                       color=color, marker=marker, s=size,
                       zorder=5 if is_reinf else 4,
                       alpha=0.9 if is_reinf else 0.6,
                       edgecolors="black", linewidths=0.4 if is_reinf else 0.2)

        ax.set_yscale("log")
        ax.set_ylim(10, 5000)
        ax.set_xlim(-10, 1200)
        ax.grid(True, alpha=0.3)
        ax.set_ylabel("Vi IgG (EU)")

        n_by_q = sub["p_resp_quintile"].value_counts().sort_index()
        q_str = ", ".join([f"Q{q}:{n_by_q.get(q, 0)}" for q in range(1, 6)])
        ax.set_title(f"{row_label} (n={len(sub)}, {q_str})", fontsize=10)

    axes[1].set_xlabel("Days since fever onset")

    # Legend: quintile colors + marker shapes
    legend_elements = []
    for q in range(1, 6):
        bounds = quintile_bounds.loc[q]
        legend_elements.append(
            Line2D([0], [0], color=q_colors[q], lw=2,
                   label=f"Q{q}: P∈[{bounds['min']:.2f},{bounds['max']:.2f}]"))
    legend_elements.append(Line2D([0], [0], marker="o", color="gray", lw=0,
                                   markersize=6, label="No reinf flag"))
    legend_elements.append(Line2D([0], [0], marker="*", color="gray", lw=0,
                                   markersize=12, label="Aiemjoy reinf algo"))
    axes[0].legend(handles=legend_elements, fontsize=7, loc="upper right", ncol=2)

    fig.suptitle("3-4 point Vi IgG trajectories: first-two-point gap × P(resp) quintile\n"
                 "(color=quintile, ★=Aiemjoy reinfection algo)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "10_trajectories_multipoint.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 10_trajectories_multipoint.png")


    # =====================================================================
    # Figure 4: 3-4 point subjects with ≥2 points by day 180
    # Columns: P(resp) quintile
    # Color: bldculres (typhi vs paratyphi)
    # Shape: circle vs star (Aiemjoy reinf algo)
    # =====================================================================

    # Add bldculres from raw CSV
    raw_serovar = raw.set_index("index_id")["bldculres"]

    # Select: ≥3 obs AND ≥2 points by day 180
    good_ids = set()
    for sid in multipoint_ids:
        grp = long[long["index_id"] == sid].sort_values("days_since_fever_onset")
        if (grp["days_since_fever_onset"] <= 180).sum() >= 2:
            good_ids.add(sid)

    df_good = df[df["index_id"].isin(good_ids)].copy()
    df_good["bldculres"] = df_good["index_id"].map(raw_serovar)

    # Recompute quintiles within this subset using the FULL-data quintile assignments
    # (keep the same quintile boundaries from the full data for comparability)

    serovar_colors = {"typhi": "#9133be", "paratyphi": "#2ca02c"}

    print(f"\nFigure 4: 3-4 pt subjects with ≥2 by day 180: {len(df_good)}")
    print(f"  typhi: {(df_good['bldculres']=='typhi').sum()}, "
          f"paratyphi: {(df_good['bldculres']=='paratyphi').sum()}")

    # How many quintiles are populated?
    populated_quintiles = sorted(df_good["p_resp_quintile"].unique())
    n_cols = len(populated_quintiles)

    fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 6), sharey=True, sharex=True)
    if n_cols == 1:
        axes = [axes]

    for col, qi in enumerate(populated_quintiles):
        ax = axes[col]
        sub = df_good[df_good["p_resp_quintile"] == qi]
        bounds = quintile_bounds.loc[qi]

        n_typhi = (sub["bldculres"] == "typhi").sum()
        n_para = (sub["bldculres"] == "paratyphi").sum()
        n_reinf = sub["reinf_algo"].sum()

        for _, r in sub.iterrows():
            sid = r["index_id"]
            grp = long[long["index_id"] == sid].sort_values("days_since_fever_onset")
            color = serovar_colors[r["bldculres"]]
            is_reinf = r["reinf_algo"]

            ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                    color=color, linewidth=1.0, alpha=0.6)
            marker = "*" if is_reinf else "o"
            size = 80 if is_reinf else 20
            ax.scatter(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                       color=color, marker=marker, s=size,
                       zorder=5 if is_reinf else 4,
                       alpha=0.9 if is_reinf else 0.6,
                       edgecolors="black", linewidths=0.4 if is_reinf else 0.2)

        ax.set_yscale("log"); ax.set_ylim(10, 5000); ax.set_xlim(-10, 1200)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Days since fever onset")
        ax.set_title(f"Q{qi}: P ∈ [{bounds['min']:.2f}, {bounds['max']:.2f}]\n"
                     f"n={len(sub)} ({n_typhi} T, {n_para} P, {n_reinf}★)",
                     fontsize=10)
        if col == 0:
            ax.set_ylabel("Vi IgG (EU)")

    legend_elements = [
        Line2D([0], [0], color="#9133be", lw=2, label="S. Typhi"),
        Line2D([0], [0], color="#2ca02c", lw=2, label="S. Paratyphi A"),
        Line2D([0], [0], marker="o", color="gray", lw=0, markersize=6, label="No reinf flag"),
        Line2D([0], [0], marker="*", color="gray", lw=0, markersize=12, label="Aiemjoy reinf algo"),
    ]
    axes[-1].legend(handles=legend_elements, fontsize=8, loc="upper right")

    fig.suptitle(f"Vi IgG trajectories: ≥3 points, ≥2 by day 180 (n={len(df_good)})\n"
                 f"Colored by blood culture result, faceted by P(resp) quintile",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "10_trajectories_multipoint_serovar.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 10_trajectories_multipoint_serovar.png")


if __name__ == "__main__":
    main()
