"""
09: Compare mixture-model P(responder) to fold-change threshold rules.

The mixture model (Step 4) gives a continuous probability of recent
Vi-boosting infection. Threshold rules (e.g., FC ≥ 3×) are the standard
seroincidence approach. They should imperfectly correlate, illustrating
the censoring limitations of hard thresholds.

Comparisons:
  1. Vi IgG fold change thresholds (2×, 3×, 4×)
  2. Aiemjoy multi-antigen reinfection rule (≥3× in ≥2 antigens at ≥3mo)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys

# Import Step 4 fitting from script 08
sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module
m08 = import_module("08_fold_change_mixture")

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")

# Import reinfection algo from script 07
sys.path.insert(0, str(Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/reinf_obs")))
m07 = import_module("07_reinf_obs_investigation")


def main():
    # Load data and fit models
    df = m08.load_fold_change_data()
    x = df["log2_fc"].values
    t0 = df["t_start"].values
    t1 = df["t_end"].values

    p1 = m08.fit_two_gaussian(x)
    z = m08.tau_ratio_covariate(t0, t1)
    p4 = m08.fit_teunis_mixture(x, t0, t1, p1)

    df["p_resp"] = 1 - p4["p_noise"]
    df["fc"] = 2 ** df["log2_fc"]

    # Add serovar
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    serovar_map = long.drop_duplicates("index_id").set_index("index_id")["serovar"]
    df["serovar"] = df["index_id"].map(serovar_map)

    # Apply Aiemjoy reinfection algo
    raw = pd.read_csv(DATA_DIR / "TypoidCaseData_github_09.30.21.csv")
    reinf = m07.apply_reinfection_algo(raw)
    df = df.merge(reinf[["index_id", "reinf_algo"]], on="index_id", how="left")

    # Vi IgG threshold classifications
    for thresh in [2, 3, 4]:
        df[f"vi_fc_ge_{thresh}x"] = df["fc"] >= thresh

    print(f"Subjects: {len(df)} ({(df['serovar']=='typhi').sum()} typhi, "
          f"{(df['serovar']=='paratyphi').sum()} paratyphi)")
    print()

    # Cross-tabs
    print("=" * 60)
    print("CROSS-TABS: P(responder) > 0.5 vs threshold rules")
    print("=" * 60)
    df["mixture_resp"] = df["p_resp"] > 0.5

    for rule_col, rule_name in [
        ("vi_fc_ge_2x", "Vi FC ≥ 2×"),
        ("vi_fc_ge_3x", "Vi FC ≥ 3×"),
        ("vi_fc_ge_4x", "Vi FC ≥ 4×"),
        ("reinf_algo", "Multi-antigen ≥3× in ≥2 Ag"),
    ]:
        ct = pd.crosstab(df["mixture_resp"].map({True: "Mixture+", False: "Mixture-"}),
                         df[rule_col].map({True: f"{rule_name}+", False: f"{rule_name}-"}),
                         margins=True)
        print(f"\n{rule_name}:")
        print(ct)
        # Agreement rate
        agree = ((df["mixture_resp"] == df[rule_col]).sum()) / len(df)
        print(f"  Agreement: {agree:.1%}")

    # By serovar
    print("\n" + "=" * 60)
    print("BY SEROVAR")
    print("=" * 60)
    for sv in ["typhi", "paratyphi"]:
        sub = df[df["serovar"] == sv]
        print(f"\n  {sv} (n={len(sub)}):")
        print(f"    Mixture P>0.5: {(sub['mixture_resp']).sum()} ({100*(sub['mixture_resp']).mean():.0f}%)")
        for thresh in [2, 3, 4]:
            n = sub[f"vi_fc_ge_{thresh}x"].sum()
            print(f"    Vi FC ≥ {thresh}×:   {n} ({100*n/len(sub):.0f}%)")
        n_reinf = sub["reinf_algo"].sum()
        print(f"    Multi-Ag rule:  {n_reinf} ({100*n_reinf/len(sub):.0f}%)")

    # =====================================================================
    # Figure
    # =====================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Panel A: P(resp) vs Vi FC, colored by serovar
    ax = axes[0, 0]
    for sv, color, marker in [("typhi", "#9133be", "o"), ("paratyphi", "#2ca02c", "s")]:
        sub = df[df["serovar"] == sv]
        ax.scatter(sub["fc"], sub["p_resp"], c=color, s=25, alpha=0.6,
                   marker=marker, edgecolors="black", linewidths=0.3, label=sv)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, label="Mixture P=0.5")
    for thresh, ls in [(2, ":"), (3, "--"), (4, "-.")]:
        ax.axvline(thresh, color="red", linestyle=ls, linewidth=0.8, label=f"FC={thresh}×")
    ax.set_xscale("log")
    ax.set_xlabel("Vi IgG fold change (last/first)")
    ax.set_ylabel("P(responder) — mixture model")
    ax.set_title("Mixture probability vs fold-change thresholds")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(True, alpha=0.3)

    # Panel B: P(resp) distribution by threshold classification (3×)
    ax = axes[0, 1]
    above = df[df["vi_fc_ge_3x"]]["p_resp"]
    below = df[~df["vi_fc_ge_3x"]]["p_resp"]
    ax.hist(below, bins=20, alpha=0.6, color="blue", label=f"FC < 3× (n={len(below)})",
            edgecolor="black", linewidth=0.5)
    ax.hist(above, bins=20, alpha=0.6, color="red", label=f"FC ≥ 3× (n={len(above)})",
            edgecolor="black", linewidth=0.5)
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("P(responder)")
    ax.set_ylabel("Count")
    ax.set_title("Mixture P(resp) by 3× threshold\n"
                 f"({len(below[below>0.5])}/{len(below)} sub-threshold assigned as responders)")
    ax.legend(fontsize=9)

    # Panel C: P(resp) vs starting EU, colored by serovar
    ax = axes[1, 0]
    for sv, color, marker in [("typhi", "#9133be", "o"), ("paratyphi", "#2ca02c", "s")]:
        sub = df[df["serovar"] == sv]
        ax.scatter(sub["eu_start"], sub["p_resp"], c=color, s=25, alpha=0.6,
                   marker=marker, edgecolors="black", linewidths=0.3, label=sv)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Starting Vi IgG (EU)")
    ax.set_ylabel("P(responder)")
    ax.set_title("Mixture assignment vs starting antibody\n(lower starting EU → more likely responder)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel D: Sensitivity analysis — what fraction detected at each threshold?
    ax = axes[1, 1]
    thresholds = np.arange(1.0, 6.1, 0.1)
    for sv, color, ls in [("typhi", "#9133be", "-"), ("paratyphi", "#2ca02c", "--"), (None, "black", "-")]:
        if sv is None:
            sub = df
            label = f"All (n={len(sub)})"
        else:
            sub = df[df["serovar"] == sv]
            label = f"{sv} (n={len(sub)})"
        # Fraction with P(resp) > 0.5 that would be caught by each threshold
        resp = sub[sub["p_resp"] > 0.5]
        if len(resp) > 0:
            fracs = [100 * (resp["fc"] >= t).mean() for t in thresholds]
            ax.plot(thresholds, fracs, color=color, linestyle=ls, linewidth=2, label=label)
    ax.set_xlabel("Fold-change threshold")
    ax.set_ylabel("% of mixture-responders detected")
    ax.set_title("Sensitivity of threshold rule\n(% of P(resp)>0.5 subjects caught)")
    ax.axvline(3, color="red", linestyle="--", linewidth=0.8, label="3× threshold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1, 6)
    ax.set_ylim(0, 100)

    fig.suptitle("Mixture model vs threshold rules for Vi IgG seroconversion",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "09_threshold_vs_mixture.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 09_threshold_vs_mixture.png")


if __name__ == "__main__":
    main()
