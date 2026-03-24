"""
04: Slope analysis of Vi IgG trajectories using full CSV dataset.

Expanded version of the PDF-based analysis with more subjects,
exact values, stratification by serovar/age/reinfection status.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis")


def compute_slopes(long):
    """Compute first-to-last log10 slope for each longitudinal subject."""
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
            "serovar": grp.iloc[0]["serovar"],
            "age": grp.iloc[0]["age"],
            "age_group": "<5" if grp.iloc[0]["age"] < 5 else ("5-15" if grp.iloc[0]["age"] <= 15 else "16+"),
            "reinf_obs": grp.iloc[0]["reinf_obs"],
            "n_obs": len(grp),
            "t_start": t0,
            "t_end": t1,
            "duration_days": dt,
            "eu_start": eu0,
            "eu_end": eu1,
            "log10_slope_per_day": (np.log10(eu1) - np.log10(eu0)) / dt,
            "fold_change": eu1 / eu0,
            "direction": "declining" if eu1 < eu0 else "rising",
        })
    return pd.DataFrame(slopes)


def print_summary(slopes, label="All"):
    """Print slope summary statistics."""
    n = len(slopes)
    if n == 0:
        print(f"  {label}: no data")
        return
    n_dec = (slopes["direction"] == "declining").sum()
    n_ris = (slopes["direction"] == "rising").sum()
    ls = slopes["log10_slope_per_day"]
    fc = slopes["fold_change"]
    print(f"  {label} (n={n}): {n_dec} declining, {n_ris} rising")
    print(f"    Median log10 slope/day: {ls.median():.6f}, mean: {ls.mean():.6f}")
    print(f"    IQR: [{ls.quantile(0.25):.6f}, {ls.quantile(0.75):.6f}]")
    print(f"    Median fold change: {fc.median():.3f}, IQR: [{fc.quantile(0.25):.3f}, {fc.quantile(0.75):.3f}]")
    print(f"    Median duration: {slopes['duration_days'].median():.0f} days")


def main():
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    slopes = compute_slopes(long)

    print("=" * 60)
    print("SLOPE ANALYSIS — FULL CSV DATASET")
    print("=" * 60)
    print_summary(slopes, "All longitudinal")
    print()
    for sv in ["typhi", "paratyphi"]:
        print_summary(slopes[slopes["serovar"] == sv], f"  {sv}")
    print()
    for ag in ["<5", "5-15", "16+"]:
        print_summary(slopes[slopes["age_group"] == ag], f"  age {ag}")
    print()
    long_dur = slopes[slopes["duration_days"] > 200]
    print_summary(long_dur, "Duration > 200 days")

    # =====================================================================
    # Figure 1: Slope distributions (2x2)
    # =====================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    for sv, color in [("typhi", "#9133be"), ("paratyphi", "#2ca02c")]:
        sub = slopes[slopes["serovar"] == sv]
        ax.hist(sub["log10_slope_per_day"], bins=25, alpha=0.6, color=color,
                label=f"{sv} (n={len(sub)})", edgecolor="black", linewidth=0.5)
    ax.axvline(0, color="red", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Log10(EU) slope per day")
    ax.set_ylabel("Count")
    ax.set_title("Slope distribution by serovar")
    ax.legend(fontsize=9)

    ax = axes[0, 1]
    for sv, color in [("typhi", "#9133be"), ("paratyphi", "#2ca02c")]:
        sub = slopes[slopes["serovar"] == sv]
        ax.hist(sub["fold_change"], bins=25, alpha=0.6, color=color,
                label=f"{sv}", edgecolor="black", linewidth=0.5)
    ax.axvline(1.0, color="red", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Fold change (end / start)")
    ax.set_ylabel("Count")
    ax.set_title("Fold change distribution")
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    colors = ["#d62728" if d == "declining" else "#2ca02c" for d in slopes["direction"]]
    ax.scatter(slopes["duration_days"], slopes["log10_slope_per_day"],
               c=colors, s=20, alpha=0.5)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5)
    ax.set_xlabel("Observation duration (days)")
    ax.set_ylabel("Log10(EU) slope per day")
    ax.set_title("Slope vs. duration (red=declining, green=rising)")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.scatter(slopes["eu_start"], slopes["log10_slope_per_day"],
               c=colors, s=20, alpha=0.5)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("Starting ELISA units")
    ax.set_ylabel("Log10(EU) slope per day")
    ax.set_title("Slope vs. starting antibody level")
    ax.grid(True, alpha=0.3)

    n_dec = (slopes["direction"] == "declining").sum()
    n_ris = (slopes["direction"] == "rising").sum()
    fig.suptitle(f"Vi IgG waning signal: {len(slopes)} trajectories from CSV\n"
                 f"({n_dec} declining, {n_ris} rising — "
                 f"median fold change {slopes['fold_change'].median():.3f})",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "04_slope_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 04_slope_analysis.png")

    # =====================================================================
    # Figure 2: Faceted trajectories (npoints × direction)
    # =====================================================================
    npts_vals = sorted(slopes["n_obs"].unique())
    directions = ["declining", "rising"]
    dir_colors = {"declining": "#d62728", "rising": "#2ca02c"}

    n_cols = len(npts_vals)
    fig, axes = plt.subplots(2, n_cols, figsize=(5 * n_cols, 8), sharex=True, sharey=True)

    for col, npts in enumerate(npts_vals):
        for row, direction in enumerate(directions):
            ax = axes[row, col]
            sub_ids = slopes[(slopes["n_obs"] == npts) & (slopes["direction"] == direction)]["index_id"]
            sub_data = long[long["index_id"].isin(sub_ids)]
            for sid, grp in sub_data.groupby("index_id"):
                grp = grp.sort_values("days_since_fever_onset")
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
                dir_label = "Declining" if direction == "declining" else "Rising / flat"
                ax.set_ylabel(f"{dir_label} (n={len(sub_ids)})\nVi IgG (EU)")
            else:
                ax.text(0.02, 0.95, f"n={len(sub_ids)}", transform=ax.transAxes,
                        fontsize=9, va="top", fontweight="bold", color=dir_colors[direction])
            if row == 1:
                ax.set_xlabel("Days since fever onset")

    fig.suptitle("Vi IgG trajectories from CSV: time points × direction",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "04_faceted_trajectories.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 04_faceted_trajectories.png")

    # =====================================================================
    # Figure 3: Spaghetti colored by serovar
    # =====================================================================
    fig, ax = plt.subplots(figsize=(10, 6))
    for sid, grp in long[long["has_longitudinal"]].groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        color = "#9133be" if grp["serovar"].iloc[0] == "typhi" else "#2ca02c"
        ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                color=color, linewidth=0.8, alpha=0.4, marker=".", markersize=2)
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("Anti-Vi IgG (ELISA units)")
    ax.set_yscale("log")
    ax.set_ylim(10, 5000)
    ax.set_title(f"All longitudinal Vi IgG trajectories (n={len(slopes)})\n"
                 f"Purple = S. Typhi, Green = S. Paratyphi A")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "04_spaghetti_by_serovar.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 04_spaghetti_by_serovar.png")


if __name__ == "__main__":
    main()
