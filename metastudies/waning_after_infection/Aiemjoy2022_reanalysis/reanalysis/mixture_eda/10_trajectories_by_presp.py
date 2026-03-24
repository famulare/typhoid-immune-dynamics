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

    print(f"Subjects: {len(df)}")
    print(f"P(resp) quintile bounds:")
    for q, row in quintile_bounds.iterrows():
        n = (df["p_resp_quintile"] == q).sum()
        n_reinf = df[df["p_resp_quintile"] == q]["reinf_algo"].sum()
        print(f"  Q{q}: P ∈ [{row['min']:.3f}, {row['max']:.3f}], n={n}, reinf_algo={n_reinf}")

    # Plot: 1 row × 5 columns (one per quintile)
    fig, axes = plt.subplots(1, 5, figsize=(25, 5), sharey=True, sharex=True)
    dir_colors = {"declining": "#d62728", "rising": "#2ca02c"}

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

            # Line
            ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                    color=color, linewidth=0.8, alpha=0.5)

            # Points: stars for reinf_algo, circles otherwise
            if is_reinf:
                ax.scatter(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                           color=color, marker="*", s=60, zorder=5, alpha=0.8,
                           edgecolors="black", linewidths=0.3)
            else:
                ax.scatter(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                           color=color, marker="o", s=15, zorder=4, alpha=0.5,
                           edgecolors="black", linewidths=0.2)

        ax.set_yscale("log")
        ax.set_ylim(10, 5000)
        ax.set_xlim(-10, 1200)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Days since fever onset")
        ax.set_title(f"Q{qi}: P(resp) ∈ [{bounds['min']:.2f}, {bounds['max']:.2f}]\n"
                     f"n={len(sub_ids)} ({n_dec}↓ {n_ris}↑, {n_reinf}★)",
                     fontsize=9)
        if qi == 1:
            ax.set_ylabel("Vi IgG (EU)")

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="#d62728", lw=1.5, label="Declining (FC<1)"),
        Line2D([0], [0], color="#2ca02c", lw=1.5, label="Rising (FC≥1)"),
        Line2D([0], [0], marker="o", color="gray", lw=0, markersize=5, label="No reinf flag"),
        Line2D([0], [0], marker="*", color="gray", lw=0, markersize=10, label="Aiemjoy reinf algo"),
    ]
    axes[-1].legend(handles=legend_elements, fontsize=7, loc="upper right")

    fig.suptitle("Vi IgG trajectories by Step 4 P(responder) quintile\n"
                 "(red=declining, green=rising, ★=Aiemjoy reinfection algo flagged)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "10_trajectories_by_presp.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 10_trajectories_by_presp.png")


if __name__ == "__main__":
    main()
