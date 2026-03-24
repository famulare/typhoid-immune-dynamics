"""
06: Descriptive summary of the Vi IgG dataset.

Characterize the full dataset: observation patterns, time coverage,
value distributions, missing data.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis")


def main():
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    print(f"Loaded {len(long)} observations from {long['index_id'].nunique()} subjects")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # 1. Observations per subject
    ax = axes[0, 0]
    obs_counts = long.groupby("index_id").size()
    ax.bar(obs_counts.value_counts().sort_index().index,
           obs_counts.value_counts().sort_index().values,
           color="#9133be", alpha=0.7, edgecolor="black")
    ax.set_xlabel("Number of Vi IgG observations")
    ax.set_ylabel("Number of subjects")
    ax.set_title(f"Observations per subject (n={len(obs_counts)})")
    for x, y in obs_counts.value_counts().sort_index().items():
        ax.text(x, y + 2, str(y), ha="center", fontsize=9)

    # 2. Time coverage: first and last observation
    ax = axes[0, 1]
    longitudinal = long[long["has_longitudinal"]]
    first_t = longitudinal.groupby("index_id")["days_since_fever_onset"].min()
    last_t = longitudinal.groupby("index_id")["days_since_fever_onset"].max()
    duration = last_t - first_t
    ax.hist(duration, bins=30, color="#9133be", alpha=0.7, edgecolor="black")
    ax.set_xlabel("Follow-up duration (days)")
    ax.set_ylabel("Count")
    ax.set_title(f"Follow-up duration (longitudinal, n={len(duration)})\n"
                 f"Median {duration.median():.0f}d, range {duration.min():.0f}-{duration.max():.0f}d")
    ax.axvline(duration.median(), color="red", linestyle="--")

    # 3. First observation time
    ax = axes[0, 2]
    ax.hist(first_t, bins=30, color="#2ca02c", alpha=0.7, edgecolor="black")
    ax.set_xlabel("First observation (days since fever onset)")
    ax.set_ylabel("Count")
    ax.set_title(f"Timing of first Vi IgG sample\n"
                 f"Median {first_t.median():.0f}d")
    ax.axvline(first_t.median(), color="red", linestyle="--")

    # 4. Vi IgG distribution by visit number
    ax = axes[1, 0]
    visit_data = []
    for v in sorted(long["visit_number"].unique()):
        vals = long[long["visit_number"] == v]["vi_igg_eu"]
        vals = vals[vals > 0]
        if len(vals) > 0:
            visit_data.append(vals.values)
    bp = ax.boxplot(visit_data, labels=[f"V{v}" for v in sorted(long["visit_number"].unique())],
                    patch_artist=True, showfliers=True,
                    flierprops=dict(markersize=3, alpha=0.3))
    for patch in bp["boxes"]:
        patch.set_facecolor("#9133be")
        patch.set_alpha(0.5)
    ax.set_yscale("log")
    ax.set_xlabel("Visit number")
    ax.set_ylabel("Vi IgG (ELISA units)")
    ax.set_title("Vi IgG by visit number")
    ax.grid(True, alpha=0.3, axis="y")

    # 5. By age group
    ax = axes[1, 1]
    subj = long.drop_duplicates("index_id").copy()
    subj["age_group"] = pd.cut(subj["age"], bins=[0, 5, 15, 100], labels=["<5", "5-15", "16+"])
    age_data = [long[long["index_id"].isin(subj[subj["age_group"] == ag]["index_id"])]["vi_igg_eu"].values
                for ag in ["<5", "5-15", "16+"]]
    age_data = [d[d > 0] for d in age_data]
    bp = ax.boxplot(age_data, labels=["<5\n(n={})".format(len(age_data[0])),
                                       "5-15\n(n={})".format(len(age_data[1])),
                                       "16+\n(n={})".format(len(age_data[2]))],
                    patch_artist=True, showfliers=True,
                    flierprops=dict(markersize=3, alpha=0.3))
    for patch in bp["boxes"]:
        patch.set_facecolor("#9133be")
        patch.set_alpha(0.5)
    ax.set_yscale("log")
    ax.set_ylabel("Vi IgG (ELISA units)")
    ax.set_title("Vi IgG by age group (all visits)")
    ax.grid(True, alpha=0.3, axis="y")

    # 6. By serovar
    ax = axes[1, 2]
    for sv, color, offset in [("typhi", "#9133be", -0.15), ("paratyphi", "#2ca02c", 0.15)]:
        sub = long[(long["serovar"] == sv) & (long["vi_igg_eu"] > 0)]
        jitter = np.random.RandomState(42).uniform(-0.1, 0.1, len(sub))
        ax.scatter(sub["days_since_fever_onset"] + jitter * 10, sub["vi_igg_eu"],
                   color=color, s=8, alpha=0.3, label=f"{sv} (n={len(sub)})")
    ax.set_yscale("log")
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("Vi IgG (ELISA units)")
    ax.set_title("All Vi IgG observations by serovar")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Aiemjoy 2022 Vi IgG Dataset Summary (from CSV)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "06_descriptive_summary.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 06_descriptive_summary.png")


if __name__ == "__main__":
    main()
