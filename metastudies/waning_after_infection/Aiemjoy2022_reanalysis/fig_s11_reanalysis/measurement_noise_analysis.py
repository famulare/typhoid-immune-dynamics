"""
Measurement noise analysis for Aiemjoy 2022 Vi IgG trajectories.

Given the assay CV, what fraction of observed trajectory slopes are
consistent with zero true change? What is the minimum detectable
waning rate for this assay and study design?

The core idea: each ELISA measurement is the true value multiplied by
lognormal noise. A fold-change computed from two noisy measurements
compounds both errors. We simulate the distribution of observed
fold-changes under the null hypothesis of no true waning, then compare
to the actual observed distribution.
"""

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

EXTRACTION_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/fig_s11_extraction")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/fig_s11_reanalysis")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Plausible CV range for Vi IgG ELISA
# HlyE/LPS CVs from Table S3 range 0.11-0.31; Vi not reported
# Vi measured at 1:100 (lower dilution than HlyE 1:500), potentially noisier
CV_VALUES = [0.15, 0.20, 0.25, 0.30]
N_SIM = 10000  # simulations per trajectory


def load_trajectories():
    """Load individual trajectories from extracted CSV."""
    with open(EXTRACTION_DIR / "panel_A_trajectories.csv") as f:
        rows = list(csv.DictReader(f))

    trajs = {}
    for r in rows:
        if r["trajectory_type"] != "individual":
            continue
        tid = int(r["trajectory_id"])
        trajs.setdefault(tid, []).append((
            float(r["days_since_fever_onset"]),
            float(r["elisa_units"]),
        ))
    for pts in trajs.values():
        pts.sort()
    return trajs


def observed_fold_changes(trajs):
    """Compute observed first-to-last fold change for each trajectory."""
    fcs = []
    for tid, pts in trajs.items():
        if len(pts) < 2:
            continue
        t0, eu0 = pts[0]
        t1, eu1 = pts[-1]
        dt = t1 - t0
        if dt > 0 and eu0 > 0 and eu1 > 0:
            fcs.append({
                "tid": tid,
                "dt": dt,
                "eu_start": eu0,
                "eu_end": eu1,
                "fold_change": eu1 / eu0,
                "log2_fc": np.log2(eu1 / eu0),
            })
    return fcs


def simulate_null_fold_changes(fcs, cv, rng):
    """Simulate fold-changes under H0: no true waning.

    Model: observed = true * exp(N(0, sigma^2))
    where sigma = sqrt(log(1 + CV^2)) for lognormal noise.

    Under H0, both measurements sample the same true value,
    so observed FC = exp(e2) / exp(e1) = exp(e2 - e1)
    where e1, e2 ~ N(0, sigma^2) iid.

    The log2(FC) under null is N(0, 2*sigma^2 / ln(2)^2).
    """
    sigma = np.sqrt(np.log(1 + cv ** 2))

    # For each real trajectory, simulate N_SIM null fold-changes
    # preserving the actual number of time points
    null_log2_fcs = []
    for fc in fcs:
        # Two independent measurement errors
        e1 = rng.normal(0, sigma, N_SIM)
        e2 = rng.normal(0, sigma, N_SIM)
        sim_log2_fc = (e2 - e1) / np.log(2)
        null_log2_fcs.append(sim_log2_fc)

    return null_log2_fcs


def compute_consistency_with_null(fcs, null_log2_fcs):
    """For each observed trajectory, compute the probability of seeing
    a fold-change at least as extreme under the null."""
    p_values = []
    for i, fc in enumerate(fcs):
        obs = fc["log2_fc"]
        null = null_log2_fcs[i]
        # Two-sided: fraction of null sims with |log2_fc| >= |obs|
        p = np.mean(np.abs(null) >= np.abs(obs))
        p_values.append(p)
    return p_values


def minimum_detectable_waning_rate(cv, durations, power=0.80, alpha=0.05):
    """What true waning rate (log10 EU/day) would be detectable with
    given power, for a single trajectory of given duration?

    Under H1: true log2(FC) = rate * dt (negative = waning)
    Under H0: observed log2(FC) ~ N(0, 2*sigma^2/ln2^2)

    For a one-sample z-test on the mean of n trajectories:
    detectable effect = (z_alpha + z_power) * SD / sqrt(n)

    But we can also ask per-trajectory: what rate makes the
    signal-to-noise ratio = 1 for a given duration?
    SNR = |rate * dt| / sqrt(2) / sigma_log2
    """
    sigma = np.sqrt(np.log(1 + cv ** 2))
    sigma_log2 = sigma / np.log(2)  # SD of a single measurement in log2 space
    sigma_fc = np.sqrt(2) * sigma_log2  # SD of log2(fold-change) from 2 measurements

    # Per-trajectory: minimum |true log2 FC| detectable at SNR=1
    # |true_log2_fc| > sigma_fc  =>  |rate_log2_per_day| > sigma_fc / dt
    rates = []
    for dt in durations:
        # Convert to log10 per day (more interpretable)
        min_rate_log2 = sigma_fc / dt
        min_rate_log10 = min_rate_log2 * np.log10(2)
        rates.append({"dt": dt, "min_rate_log10_per_day": min_rate_log10,
                       "half_life_days": np.log10(2) / min_rate_log10 if min_rate_log10 > 0 else np.inf})

    # Population-level: detectable mean rate across n trajectories
    n = len(durations)
    median_dt = np.median(durations)
    z_alpha = stats.norm.ppf(1 - alpha / 2)
    z_power = stats.norm.ppf(power)
    # Minimum detectable mean log2 rate = (z_a + z_p) * sigma_fc / (median_dt * sqrt(n))
    min_mean_log2_rate = (z_alpha + z_power) * sigma_fc / (median_dt * np.sqrt(n))
    min_mean_log10_rate = min_mean_log2_rate * np.log10(2)
    pop_half_life = np.log10(2) / min_mean_log10_rate if min_mean_log10_rate > 0 else np.inf

    return rates, {
        "n_trajectories": n,
        "median_duration_days": median_dt,
        "min_mean_log10_rate_per_day": min_mean_log10_rate,
        "implied_half_life_days": pop_half_life,
        "cv": cv,
        "power": power,
        "alpha": alpha,
    }


# =============================================================================
# Plots
# =============================================================================

def plot_null_comparison(fcs, cv_results):
    """Compare observed fold-change distribution to null simulations."""
    n_cvs = len(cv_results)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    obs_log2_fcs = [fc["log2_fc"] for fc in fcs]

    for i, (cv, (null_log2_fcs, p_values)) in enumerate(cv_results.items()):
        ax = axes[i]

        # Pool all null simulations for histogram
        all_null = np.concatenate(null_log2_fcs)

        ax.hist(all_null, bins=80, density=True, alpha=0.4, color="gray",
                label=f"Null (CV={cv})")
        ax.hist(obs_log2_fcs, bins=30, density=True, alpha=0.7, color="#9133be",
                label="Observed", histtype="step", linewidth=2)
        ax.axvline(0, color="black", linestyle="--", linewidth=0.5)

        # KS test
        ks_stat, ks_p = stats.ks_2samp(obs_log2_fcs, all_null[:len(obs_log2_fcs) * 10])

        n_consistent = sum(1 for p in p_values if p > 0.05)
        ax.set_title(f"CV = {cv}\n"
                     f"{n_consistent}/{len(p_values)} trajectories consistent with null (p>0.05)\n"
                     f"KS test vs null: p={ks_p:.3f}",
                     fontsize=10)
        ax.set_xlabel("Log2(fold change)")
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)
        ax.set_xlim(-4, 4)

    fig.suptitle("Observed Vi IgG fold-changes vs. measurement noise null\n"
                 "(null = no true waning, all variation from assay CV)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "null_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved null_comparison.png")


def plot_detectable_rates(cv_rates):
    """Plot minimum detectable waning rate as function of CV and duration."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: per-trajectory minimum detectable half-life vs duration
    ax = axes[0]
    for cv, (rates, pop) in cv_rates.items():
        dts = [r["dt"] for r in rates]
        hls = [r["half_life_days"] for r in rates]
        ax.scatter(dts, hls, s=15, alpha=0.5, label=f"CV={cv}")

    ax.set_xlabel("Observation duration (days)")
    ax.set_ylabel("Minimum detectable half-life (days)\n(per-trajectory, SNR=1)")
    ax.set_title("Per-trajectory detection limit")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 3000)
    ax.grid(True, alpha=0.3)
    ax.axhline(365, color="red", linestyle="--", alpha=0.5, label="1 year")
    ax.axhline(730, color="orange", linestyle="--", alpha=0.5, label="2 years")

    # Right: population-level minimum detectable half-life vs CV
    ax = axes[1]
    cvs = sorted(cv_rates.keys())
    pop_hls = [cv_rates[cv][1]["implied_half_life_days"] for cv in cvs]
    ax.bar(range(len(cvs)), pop_hls, color="#9133be", alpha=0.7)
    ax.set_xticks(range(len(cvs)))
    ax.set_xticklabels([f"CV={cv}" for cv in cvs])
    ax.set_ylabel("Minimum detectable half-life (days)\n(population mean, 80% power)")

    for j, hl in enumerate(pop_hls):
        ax.text(j, hl + 5, f"{hl:.0f}d\n({hl/365:.1f}y)", ha="center", fontsize=9)

    n = cv_rates[cvs[0]][1]["n_trajectories"]
    med_dur = cv_rates[cvs[0]][1]["median_duration_days"]
    ax.set_title(f"Population-level detection limit\n"
                 f"(n={n} trajectories, median duration={med_dur:.0f}d, "
                 f"α=0.05, power=0.80)")
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle("Minimum detectable Vi IgG waning rate given assay noise",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "detectable_waning_rates.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved detectable_waning_rates.png")


def plot_p_value_distribution(cv_results):
    """Plot distribution of per-trajectory p-values for each CV."""
    fig, ax = plt.subplots(figsize=(8, 5))

    for cv, (_, p_values) in sorted(cv_results.items()):
        p_sorted = np.sort(p_values)
        ax.plot(np.linspace(0, 1, len(p_sorted)), p_sorted,
                label=f"CV={cv}", linewidth=2)

    ax.plot([0, 1], [0, 1], "k--", linewidth=0.5, label="Uniform (expected under null)")
    ax.set_xlabel("Expected quantile")
    ax.set_ylabel("Observed p-value")
    ax.set_title("P-P plot: per-trajectory p-values vs. null\n"
                 "(diagonal = data indistinguishable from pure noise)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "pvalue_pp_plot.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved pvalue_pp_plot.png")


# =============================================================================
# Main
# =============================================================================

def main():
    rng = np.random.default_rng(42)
    trajs = load_trajectories()
    fcs = observed_fold_changes(trajs)
    durations = [fc["dt"] for fc in fcs]

    print("=" * 60)
    print("MEASUREMENT NOISE ANALYSIS")
    print("=" * 60)
    print(f"Trajectories with computable fold-change: {len(fcs)}")
    print(f"Duration range: {min(durations):.0f}-{max(durations):.0f} days "
          f"(median {np.median(durations):.0f})")

    obs_log2 = [fc["log2_fc"] for fc in fcs]
    print(f"Observed log2(FC): median={np.median(obs_log2):.3f}, "
          f"SD={np.std(obs_log2):.3f}")

    # Run analysis for each CV
    cv_results = {}  # cv -> (null_log2_fcs, p_values)
    cv_rates = {}    # cv -> (per_traj_rates, pop_stats)

    for cv in CV_VALUES:
        sigma = np.sqrt(np.log(1 + cv ** 2))
        sigma_log2 = sigma / np.log(2)
        sigma_fc = np.sqrt(2) * sigma_log2

        print(f"\n--- CV = {cv} ---")
        print(f"  Lognormal sigma: {sigma:.4f}")
        print(f"  SD of single measurement (log2): {sigma_log2:.4f}")
        print(f"  SD of fold-change (log2): {sigma_fc:.4f}")
        print(f"  Observed SD of log2(FC): {np.std(obs_log2):.4f}")
        print(f"  Ratio observed/expected SD: {np.std(obs_log2) / sigma_fc:.2f}")

        null_log2_fcs = simulate_null_fold_changes(fcs, cv, rng)
        p_values = compute_consistency_with_null(fcs, null_log2_fcs)

        n_consistent = sum(1 for p in p_values if p > 0.05)
        print(f"  Trajectories consistent with null (p>0.05): "
              f"{n_consistent}/{len(p_values)} ({100*n_consistent/len(p_values):.0f}%)")

        # KS test: is the observed FC distribution distinguishable from null?
        all_null = np.concatenate(null_log2_fcs)
        ks_stat, ks_p = stats.ks_2samp(obs_log2, all_null[:len(obs_log2) * 10])
        print(f"  KS test (observed vs null): stat={ks_stat:.3f}, p={ks_p:.3f}")

        cv_results[cv] = (null_log2_fcs, p_values)

        # Minimum detectable rates
        rates, pop = minimum_detectable_waning_rate(cv, durations)
        cv_rates[cv] = (rates, pop)
        print(f"  Population-level min detectable half-life: "
              f"{pop['implied_half_life_days']:.0f} days "
              f"({pop['implied_half_life_days']/365:.1f} years)")

    # Plots
    print("\nGenerating plots...")
    plot_null_comparison(fcs, cv_results)
    plot_detectable_rates(cv_rates)
    plot_p_value_distribution(cv_results)

    print(f"\nAll outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
