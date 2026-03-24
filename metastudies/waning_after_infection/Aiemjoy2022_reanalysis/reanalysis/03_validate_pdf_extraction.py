"""
03: Validate PDF extraction accuracy against CSV ground truth.

Match PDF-extracted trajectories to CSV subjects and quantify extraction error.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
EXTRACTION_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/fig_s11_extraction")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis")


def main():
    # Load CSV long-format
    csv_long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    csv_typhi = csv_long[(csv_long["serovar"] == "typhi") & csv_long["has_longitudinal"]].copy()

    # Load PDF extraction
    pdf = pd.read_csv(EXTRACTION_DIR / "panel_A_trajectories.csv")
    pdf_indiv = pdf[pdf["trajectory_type"] == "individual"].copy()

    print(f"CSV: {csv_typhi['index_id'].nunique()} S. Typhi longitudinal subjects, "
          f"{len(csv_typhi)} obs")
    print(f"PDF: {pdf_indiv['trajectory_id'].nunique()} extracted trajectories, "
          f"{len(pdf_indiv)} points")

    # Build per-subject summaries for matching
    csv_subj = []
    for sid, grp in csv_typhi.groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        csv_subj.append({
            "index_id": sid,
            "n_obs": len(grp),
            "t_first": grp["days_since_fever_onset"].iloc[0],
            "t_last": grp["days_since_fever_onset"].iloc[-1],
            "eu_first": grp["vi_igg_eu"].iloc[0],
            "eu_last": grp["vi_igg_eu"].iloc[-1],
        })
    csv_subj = pd.DataFrame(csv_subj)

    pdf_subj = []
    for tid, grp in pdf_indiv.groupby("trajectory_id"):
        grp = grp.sort_values("days_since_fever_onset")
        pdf_subj.append({
            "trajectory_id": tid,
            "n_points": len(grp),
            "t_first": grp["days_since_fever_onset"].iloc[0],
            "t_last": grp["days_since_fever_onset"].iloc[-1],
            "eu_first": grp["elisa_units"].iloc[0],
            "eu_last": grp["elisa_units"].iloc[-1],
        })
    pdf_subj = pd.DataFrame(pdf_subj)

    # Match by nearest (t_first, eu_first) — greedy matching
    matches = []
    csv_used = set()
    for _, p in pdf_subj.iterrows():
        best_dist = float("inf")
        best_csv_id = None
        for _, c in csv_subj.iterrows():
            if c["index_id"] in csv_used:
                continue
            if c["n_obs"] != p["n_points"]:
                continue
            # Distance in (log_time, log_eu) space
            dt = abs(p["t_first"] - c["t_first"]) / max(c["t_first"], 1)
            deu = abs(np.log10(max(p["eu_first"], 1)) - np.log10(max(c["eu_first"], 1)))
            dist = dt + deu
            if dist < best_dist:
                best_dist = dist
                best_csv_id = c["index_id"]

        if best_dist < 0.5:  # reasonable tolerance
            matches.append({
                "trajectory_id": p["trajectory_id"],
                "index_id": best_csv_id,
                "n_points": p["n_points"],
                "pdf_t_first": p["t_first"],
                "csv_t_first": csv_subj.loc[csv_subj["index_id"] == best_csv_id, "t_first"].iloc[0],
                "pdf_eu_first": p["eu_first"],
                "csv_eu_first": csv_subj.loc[csv_subj["index_id"] == best_csv_id, "eu_first"].iloc[0],
                "match_dist": best_dist,
            })
            csv_used.add(best_csv_id)

    matches = pd.DataFrame(matches)
    print(f"\nMatched: {len(matches)}/{len(pdf_subj)} PDF trajectories "
          f"({100*len(matches)/len(pdf_subj):.0f}%)")

    if len(matches) == 0:
        print("No matches found — skipping detailed validation")
        return

    # Now do point-level matching for matched subjects
    point_matches = []
    for _, m in matches.iterrows():
        pdf_pts = pdf_indiv[pdf_indiv["trajectory_id"] == m["trajectory_id"]].sort_values("days_since_fever_onset")
        csv_pts = csv_typhi[csv_typhi["index_id"] == m["index_id"]].sort_values("days_since_fever_onset")

        for (_, pp), (_, cp) in zip(pdf_pts.iterrows(), csv_pts.iterrows()):
            point_matches.append({
                "pdf_days": pp["days_since_fever_onset"],
                "csv_days": cp["days_since_fever_onset"],
                "pdf_eu": pp["elisa_units"],
                "csv_eu": cp["vi_igg_eu"],
            })

    pm = pd.DataFrame(point_matches)

    # Error statistics
    pm["day_error"] = pm["pdf_days"] - pm["csv_days"]
    pm["eu_error_pct"] = (pm["pdf_eu"] - pm["csv_eu"]) / pm["csv_eu"] * 100
    pm["eu_ratio"] = pm["pdf_eu"] / pm["csv_eu"]

    print(f"\nPoint-level validation (n={len(pm)} matched points):")
    print(f"  Day error: median={pm['day_error'].median():.1f}, "
          f"IQR=[{pm['day_error'].quantile(0.25):.1f}, {pm['day_error'].quantile(0.75):.1f}]")
    print(f"  EU error %: median={pm['eu_error_pct'].median():.1f}%, "
          f"IQR=[{pm['eu_error_pct'].quantile(0.25):.1f}%, {pm['eu_error_pct'].quantile(0.75):.1f}%]")
    print(f"  EU ratio: median={pm['eu_ratio'].median():.3f}, "
          f"IQR=[{pm['eu_ratio'].quantile(0.25):.3f}, {pm['eu_ratio'].quantile(0.75):.3f}]")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    ax.scatter(pm["csv_days"], pm["pdf_days"], s=15, alpha=0.5, color="#9133be")
    lims = [0, max(pm["csv_days"].max(), pm["pdf_days"].max()) * 1.05]
    ax.plot(lims, lims, "k--", linewidth=0.5)
    ax.set_xlabel("CSV days since fever onset")
    ax.set_ylabel("PDF-extracted days")
    ax.set_title(f"Time extraction accuracy\nMedian error: {pm['day_error'].median():.1f} days")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.scatter(pm["csv_eu"], pm["pdf_eu"], s=15, alpha=0.5, color="#9133be")
    lims = [10, max(pm["csv_eu"].max(), pm["pdf_eu"].max()) * 1.1]
    ax.plot(lims, lims, "k--", linewidth=0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("CSV ELISA units")
    ax.set_ylabel("PDF-extracted ELISA units")
    ax.set_title(f"Value extraction accuracy\nMedian error: {pm['eu_error_pct'].median():.1f}%")
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    finite_err = pm["eu_error_pct"].replace([np.inf, -np.inf], np.nan).dropna()
    ax.hist(finite_err, bins=30, color="#9133be", alpha=0.7, edgecolor="black")
    ax.axvline(0, color="red", linestyle="--")
    ax.set_xlabel("ELISA unit error (%)")
    ax.set_ylabel("Count")
    n_inf = len(pm) - len(finite_err)
    ax.set_title(f"Distribution of extraction errors\n({n_inf} infinite values from CSV zeros excluded)")

    fig.suptitle(f"PDF extraction validation: {len(matches)} trajectories matched, "
                 f"{len(pm)} points compared",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "03_pdf_extraction_validation.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved 03_pdf_extraction_validation.png")


if __name__ == "__main__":
    main()
