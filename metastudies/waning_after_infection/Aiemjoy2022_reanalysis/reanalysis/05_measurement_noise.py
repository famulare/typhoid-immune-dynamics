"""
05: Measurement noise analysis using full CSV dataset.

Same framework as the PDF-based analysis but with more subjects
and exact values.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis")

CV_VALUES = [0.15, 0.20, 0.25, 0.30]
N_SIM = 10000


def main():
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    longitudinal = long[long["has_longitudinal"]].copy()

    # Compute first-to-last fold change per subject
    fcs = []
    for sid, grp in longitudinal.groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        eu0, eu1 = grp.iloc[0]["vi_igg_eu"], grp.iloc[-1]["vi_igg_eu"]
        dt = grp.iloc[-1]["days_since_fever_onset"] - grp.iloc[0]["days_since_fever_onset"]
        if dt > 0 and eu0 > 0 and eu1 > 0:
            fcs.append({"dt": dt, "log2_fc": np.log2(eu1 / eu0)})
    fcs = pd.DataFrame(fcs)

    obs_log2 = fcs["log2_fc"].values
    durations = fcs["dt"].values

    print("=" * 60)
    print("MEASUREMENT NOISE ANALYSIS — FULL CSV")
    print("=" * 60)
    print(f"Trajectories: {len(fcs)}")
    print(f"Duration: median {np.median(durations):.0f} days, "
          f"range {durations.min():.0f}-{durations.max():.0f}")
    print(f"Observed log2(FC): median={np.median(obs_log2):.3f}, SD={np.std(obs_log2):.3f}")

    rng = np.random.default_rng(42)

    # Analyze for each CV
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes_flat = axes.flatten()

    results = []
    for i, cv in enumerate(CV_VALUES):
        sigma = np.sqrt(np.log(1 + cv ** 2))
        sigma_log2 = sigma / np.log(2)
        sigma_fc = np.sqrt(2) * sigma_log2

        # Simulate null fold-changes
        null_log2 = rng.normal(0, sigma_fc, size=len(fcs) * N_SIM)

        # KS test
        ks_stat, ks_p = stats.ks_2samp(obs_log2, null_log2[:len(obs_log2) * 10])

        # Per-trajectory p-values
        p_values = []
        for obs in obs_log2:
            p = 2 * (1 - stats.norm.cdf(abs(obs), 0, sigma_fc))
            p_values.append(p)
        n_sig = sum(1 for p in p_values if p < 0.05)

        # Population-level minimum detectable half-life
        z_a = stats.norm.ppf(0.975)
        z_p = stats.norm.ppf(0.80)
        median_dt = np.median(durations)
        n = len(fcs)
        min_rate_log2 = (z_a + z_p) * sigma_fc / (median_dt * np.sqrt(n))
        min_rate_log10 = min_rate_log2 * np.log10(2)
        min_hl = np.log10(2) / min_rate_log10

        print(f"\n--- CV = {cv} ---")
        print(f"  Expected SD of log2(FC): {sigma_fc:.4f}")
        print(f"  Observed SD: {np.std(obs_log2):.4f} (ratio: {np.std(obs_log2)/sigma_fc:.2f})")
        print(f"  KS test: stat={ks_stat:.3f}, p={ks_p:.3f}")
        print(f"  Trajectories with p<0.05: {n_sig}/{len(fcs)} ({100*n_sig/len(fcs):.0f}%)")
        print(f"  Min detectable half-life (population, 80% power): {min_hl:.0f} days ({min_hl/365:.1f} years)")

        results.append({"cv": cv, "sigma_fc": sigma_fc, "ks_p": ks_p,
                         "n_sig": n_sig, "min_hl_days": min_hl})

        # Plot
        ax = axes_flat[i]
        ax.hist(null_log2[:5000], bins=80, density=True, alpha=0.4, color="gray",
                label=f"Null (CV={cv})")
        ax.hist(obs_log2, bins=30, density=True, alpha=0.7, color="#9133be",
                label="Observed", histtype="step", linewidth=2)
        ax.axvline(0, color="black", linestyle="--", linewidth=0.5)
        n_consistent = len(fcs) - n_sig
        ax.set_title(f"CV = {cv}\n"
                     f"{n_consistent}/{len(fcs)} consistent with null\n"
                     f"KS p={ks_p:.3f} | Min HL={min_hl:.0f}d ({min_hl/365:.1f}y)",
                     fontsize=9)
        ax.set_xlabel("Log2(fold change)")
        ax.set_ylabel("Density")
        ax.legend(fontsize=7)
        ax.set_xlim(-4, 4)

    fig.suptitle(f"Vi IgG fold-changes vs. measurement noise null (n={len(fcs)} from CSV)\n"
                 f"Observed SD={np.std(obs_log2):.3f}",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "05_measurement_noise.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 05_measurement_noise.png")

    # Summary table
    results = pd.DataFrame(results)
    print(f"\n{'='*60}")
    print("SUMMARY: Minimum detectable Vi IgG half-life")
    print(f"{'='*60}")
    print(f"n={len(fcs)} trajectories, median duration={np.median(durations):.0f} days")
    for _, r in results.iterrows():
        print(f"  CV={r['cv']:.2f}: min detectable HL = {r['min_hl_days']:.0f} days "
              f"({r['min_hl_days']/365:.1f} years)")


if __name__ == "__main__":
    main()
