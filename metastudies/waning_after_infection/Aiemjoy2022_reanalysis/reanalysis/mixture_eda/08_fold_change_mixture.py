"""
08: Fold-change mixture model analysis for Vi IgG.

Three models of increasing complexity fitted to log2(fold change) data:
  Step 1: Two-Gaussian mixture (EM algorithm)
  Step 2: Gaussian + Skew-Normal mixture (MLE)
  Step 3: Gaussian + time-dependent Normal mixture regression
          (waning linear in log2(duration) — power-law decay)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, optimize
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")


def load_fold_change_data():
    """Load longitudinal Vi IgG data and compute fold changes."""
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    records = []
    for sid, grp in long[long["has_longitudinal"]].groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        t0, eu0 = grp.iloc[0]["days_since_fever_onset"], grp.iloc[0]["vi_igg_eu"]
        t1, eu1 = grp.iloc[-1]["days_since_fever_onset"], grp.iloc[-1]["vi_igg_eu"]
        dt = t1 - t0
        if dt <= 0 or eu0 <= 0 or eu1 <= 0:
            continue
        records.append({
            "index_id": sid,
            "log2_fc": np.log2(eu1 / eu0),
            "duration_days": dt,
            "log2_duration": np.log2(dt),
            "eu_start": eu0,
        })
    return pd.DataFrame(records)


def sigma_log2_to_cv(sigma_log2):
    """Convert SD on log2 scale (for fold-change of two measurements) to natural-scale CV."""
    sigma_per_meas = sigma_log2 / np.sqrt(2)
    sigma_ln = sigma_per_meas * np.log(2)
    return np.sqrt(np.exp(sigma_ln**2) - 1)


# =========================================================================
# Step 1: Two-Gaussian mixture via EM
# =========================================================================

def fit_two_gaussian(x, n_init=10, max_iter=200, tol=1e-8):
    """Fit 2-component Gaussian mixture via EM algorithm."""
    n = len(x)
    best_ll = -np.inf
    best_params = None

    rng = np.random.default_rng(42)

    for init in range(n_init):
        # Initialize: random split
        if init == 0:
            # First init: split at median
            mask = x < np.median(x)
        else:
            mask = rng.random(n) < rng.uniform(0.3, 0.7)

        mu = np.array([x[mask].mean(), x[~mask].mean()])
        sigma = np.array([x[mask].std(), x[~mask].std()])
        sigma = np.maximum(sigma, 0.1)
        pi = np.array([mask.mean(), 1 - mask.mean()])

        for _ in range(max_iter):
            # E-step
            log_resp = np.zeros((n, 2))
            for k in range(2):
                log_resp[:, k] = np.log(pi[k] + 1e-300) + stats.norm.logpdf(x, mu[k], sigma[k])
            log_resp -= np.max(log_resp, axis=1, keepdims=True)
            resp = np.exp(log_resp)
            resp /= resp.sum(axis=1, keepdims=True)

            # M-step
            nk = resp.sum(axis=0)
            pi_new = nk / n
            mu_new = (resp * x[:, None]).sum(axis=0) / nk
            sigma_new = np.sqrt((resp * (x[:, None] - mu_new)**2).sum(axis=0) / nk)
            sigma_new = np.maximum(sigma_new, 1e-6)

            if np.max(np.abs(mu_new - mu)) < tol and np.max(np.abs(sigma_new - sigma)) < tol:
                mu, sigma, pi = mu_new, sigma_new, pi_new
                break
            mu, sigma, pi = mu_new, sigma_new, pi_new

        # Log-likelihood
        ll = np.sum(np.log(pi[0] * stats.norm.pdf(x, mu[0], sigma[0]) +
                           pi[1] * stats.norm.pdf(x, mu[1], sigma[1]) + 1e-300))
        if ll > best_ll:
            best_ll = ll
            best_params = (mu.copy(), sigma.copy(), pi.copy())

    mu, sigma, pi = best_params
    # Order: component 0 = smaller mean (noise)
    order = np.argsort(mu)
    mu, sigma, pi = mu[order], sigma[order], pi[order]

    n_params = 5
    params = {
        "mu1": mu[0], "sigma1": sigma[0], "mu2": mu[1], "sigma2": sigma[1],
        "pi_noise": pi[0],
        "ll": best_ll,
        "aic": 2 * n_params - 2 * best_ll,
        "bic": n_params * np.log(n) - 2 * best_ll,
    }

    # Posterior assignments
    pdf0 = pi[0] * stats.norm.pdf(x, mu[0], sigma[0])
    pdf1 = pi[1] * stats.norm.pdf(x, mu[1], sigma[1])
    params["p_noise"] = pdf0 / (pdf0 + pdf1)

    return params


# =========================================================================
# Step 2: Gaussian + Skew-Normal mixture
# =========================================================================

def neg_ll_gauss_skewnorm(theta, x):
    """Negative log-likelihood for Gaussian + Skew-Normal mixture."""
    mu1, log_sigma1, xi, log_omega, alpha, logit_pi = theta
    sigma1 = np.exp(log_sigma1)
    omega = np.exp(log_omega)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    pdf_noise = stats.norm.pdf(x, mu1, sigma1)
    pdf_signal = stats.skewnorm.pdf(x, alpha, xi, omega)

    mixture = pi_noise * pdf_noise + (1 - pi_noise) * pdf_signal
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def fit_gauss_skewnorm(x, gauss_params):
    """Fit Gaussian + Skew-Normal mixture via MLE, initialized from Gaussian fit."""
    x0 = [
        gauss_params["mu1"],
        np.log(gauss_params["sigma1"]),
        gauss_params["mu2"],
        np.log(gauss_params["sigma2"]),
        0.0,  # alpha=0 (no skew initially)
        np.log(gauss_params["pi_noise"] / (1 - gauss_params["pi_noise"])),
    ]

    result = optimize.minimize(
        neg_ll_gauss_skewnorm, x0, args=(x,),
        method="Nelder-Mead", options={"maxiter": 50000, "xatol": 1e-8, "fatol": 1e-8},
    )

    mu1, log_sigma1, xi, log_omega, alpha, logit_pi = result.x
    sigma1 = np.exp(log_sigma1)
    omega = np.exp(log_omega)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    ll = -result.fun
    n_params = 6
    params = {
        "mu1": mu1, "sigma1": sigma1,
        "xi": xi, "omega": omega, "alpha": alpha,
        "pi_noise": pi_noise,
        "ll": ll,
        "aic": 2 * n_params - 2 * ll,
        "bic": n_params * np.log(len(x)) - 2 * ll,
        "converged": result.success,
    }

    pdf_noise = stats.norm.pdf(x, mu1, sigma1)
    pdf_signal = stats.skewnorm.pdf(x, alpha, xi, omega)
    mixture = pi_noise * pdf_noise + (1 - pi_noise) * pdf_signal
    params["p_noise"] = (pi_noise * pdf_noise) / mixture

    return params


# =========================================================================
# Step 3: Gaussian + time-dependent Normal (mixture regression)
# =========================================================================

def neg_ll_gauss_timedep(theta, x, log2_dt):
    """Negative log-likelihood for Gaussian + time-dependent Normal mixture.

    Component 1 (noise): N(mu1, sigma1) — no time dependence
    Component 2 (signal): N(mu2 - beta * log2(dt), sigma2) — power-law waning
    """
    mu1, log_sigma1, mu2, log_sigma2, beta, logit_pi = theta
    sigma1 = np.exp(log_sigma1)
    sigma2 = np.exp(log_sigma2)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    pdf_noise = stats.norm.pdf(x, mu1, sigma1)
    mu2_effective = mu2 - beta * log2_dt
    pdf_signal = stats.norm.pdf(x, mu2_effective, sigma2)

    mixture = pi_noise * pdf_noise + (1 - pi_noise) * pdf_signal
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def fit_gauss_timedep(x, log2_dt, gauss_params):
    """Fit Gaussian + time-dependent Normal, initialized from Gaussian fit."""
    x0 = [
        gauss_params["mu1"],
        np.log(gauss_params["sigma1"]),
        gauss_params["mu2"],
        np.log(gauss_params["sigma2"]),
        0.0,  # beta=0 (no waning initially)
        np.log(gauss_params["pi_noise"] / (1 - gauss_params["pi_noise"])),
    ]

    result = optimize.minimize(
        neg_ll_gauss_timedep, x0, args=(x, log2_dt),
        method="Nelder-Mead", options={"maxiter": 50000, "xatol": 1e-8, "fatol": 1e-8},
    )

    mu1, log_sigma1, mu2, log_sigma2, beta, logit_pi = result.x
    sigma1 = np.exp(log_sigma1)
    sigma2 = np.exp(log_sigma2)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    ll = -result.fun
    n_params = 6
    params = {
        "mu1": mu1, "sigma1": sigma1,
        "mu2": mu2, "sigma2": sigma2,
        "beta": beta, "pi_noise": pi_noise,
        "ll": ll,
        "aic": 2 * n_params - 2 * ll,
        "bic": n_params * np.log(len(x)) - 2 * ll,
        "converged": result.success,
    }

    pdf_noise = stats.norm.pdf(x, mu1, sigma1)
    mu2_eff = mu2 - beta * log2_dt
    pdf_signal = stats.norm.pdf(x, mu2_eff, sigma2)
    mixture = pi_noise * pdf_noise + (1 - pi_noise) * pdf_signal
    params["p_noise"] = (pi_noise * pdf_noise) / mixture

    return params


# =========================================================================
# Plotting — one figure per step + summary
# =========================================================================

def plot_step1(x, p):
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: histogram + fit
    ax = axes[0]
    ax.hist(x, bins=30, density=True, alpha=0.4, color="gray", edgecolor="black", linewidth=0.5)
    pdf_n = p["pi_noise"] * stats.norm.pdf(xgrid, p["mu1"], p["sigma1"])
    pdf_s = (1 - p["pi_noise"]) * stats.norm.pdf(xgrid, p["mu2"], p["sigma2"])
    ax.plot(xgrid, pdf_n, "b-", linewidth=2, label=f"Noise (π={p['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_s, "r-", linewidth=2, label=f"Signal (π={1-p['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_n + pdf_s, "k--", linewidth=1.5, label="Mixture")
    ax.axvline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("log2(fold change)")
    ax.set_ylabel("Density")
    cv = sigma_log2_to_cv(p["sigma1"])
    ax.set_title(f"Two-Gaussian mixture\n"
                 f"Noise: μ={p['mu1']:.3f}, σ={p['sigma1']:.3f} (CV≈{cv:.2f})\n"
                 f"Signal: μ={p['mu2']:.3f}, σ={p['sigma2']:.3f}")
    ax.legend(fontsize=9)
    ax.set_xlim(xgrid[0], xgrid[-1])

    # Panel B: posterior assignment strip
    ax = axes[1]
    p_resp = 1 - p["p_noise"]
    order = np.argsort(x)
    ax.scatter(x[order], p_resp[order], c=p_resp[order], cmap="RdYlBu_r",
               s=20, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2(fold change)")
    ax.set_ylabel("P(responder)")
    n_resp = (p_resp > 0.5).sum()
    ax.set_title(f"Posterior component assignment\n"
                 f"{n_resp}/{len(x)} classified as responders (P>0.5)")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Step 1: Two-Gaussian Mixture", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step1_two_gaussian.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 08_step1_two_gaussian.png")


def plot_step2(x, p1, p2):
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: histogram + both fits
    ax = axes[0]
    ax.hist(x, bins=30, density=True, alpha=0.4, color="gray", edgecolor="black", linewidth=0.5)

    # Step 1 fit (dashed)
    mix1 = (p1["pi_noise"] * stats.norm.pdf(xgrid, p1["mu1"], p1["sigma1"]) +
            (1 - p1["pi_noise"]) * stats.norm.pdf(xgrid, p1["mu2"], p1["sigma2"]))
    ax.plot(xgrid, mix1, "k:", linewidth=1.5, alpha=0.5, label="Step 1 (Gauss+Gauss)")

    # Step 2 fit
    pdf_n = p2["pi_noise"] * stats.norm.pdf(xgrid, p2["mu1"], p2["sigma1"])
    pdf_s = (1 - p2["pi_noise"]) * stats.skewnorm.pdf(xgrid, p2["alpha"], p2["xi"], p2["omega"])
    ax.plot(xgrid, pdf_n, "b-", linewidth=2, label=f"Noise (π={p2['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_s, "r-", linewidth=2, label=f"Signal (π={1-p2['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_n + pdf_s, "k--", linewidth=1.5, label="Mixture")
    ax.axvline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("log2(fold change)")
    ax.set_ylabel("Density")
    cv = sigma_log2_to_cv(p2["sigma1"])
    ax.set_title(f"Gaussian + Skew-Normal\n"
                 f"Noise: μ={p2['mu1']:.3f}, σ={p2['sigma1']:.3f} (CV≈{cv:.2f})\n"
                 f"Skew-Normal: α={p2['alpha']:.2f}")
    ax.legend(fontsize=8)
    ax.set_xlim(xgrid[0], xgrid[-1])

    # Panel B: posterior assignment comparison
    ax = axes[1]
    p_resp_1 = 1 - p1["p_noise"]
    p_resp_2 = 1 - p2["p_noise"]
    ax.scatter(p_resp_1, p_resp_2, c=x, cmap="coolwarm", s=20, alpha=0.7,
               edgecolors="black", linewidths=0.3, vmin=-2, vmax=2)
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.8)
    ax.set_xlabel("P(responder) — Step 1")
    ax.set_ylabel("P(responder) — Step 2")
    ax.set_title("Assignment comparison\n(color = log2 FC)")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    delta_aic = p1["aic"] - p2["aic"]
    ax.text(0.05, 0.95, f"ΔAIC = {delta_aic:.1f}\n({'Skew-Normal better' if delta_aic > 0 else 'Gaussian better'})",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle("Step 2: Gaussian + Skew-Normal Mixture", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step2_skew_normal.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 08_step2_skew_normal.png")


def plot_step3(df, p3):
    x = df["log2_fc"].values
    log2_dt = df["log2_duration"].values
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Panel A: scatter colored by P(responder)
    ax = axes[0, 0]
    p_resp = 1 - p3["p_noise"]
    sc = ax.scatter(df["duration_days"], x, c=p_resp, cmap="RdYlBu_r",
                    s=25, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    dt_grid = np.array([30, 60, 90, 180, 365, 730, 1100])
    mu2_grid = p3["mu2"] - p3["beta"] * np.log2(dt_grid)
    ax.plot(dt_grid, mu2_grid, "r-", linewidth=2, label="Signal mean")
    ax.fill_between(dt_grid, mu2_grid - p3["sigma2"], mu2_grid + p3["sigma2"],
                     color="red", alpha=0.15, label="±1σ signal")
    ax.axhline(p3["mu1"], color="blue", linestyle="--", linewidth=1, label=f"Noise mean")
    ax.set_xlabel("Duration (days)")
    ax.set_ylabel("log2(fold change)")
    ax.set_xscale("log")
    ax.set_title(f"Mixture regression: log2(FC) vs log2(duration)\n"
                 f"β={p3['beta']:.4f} (power-law exponent)")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, alpha=0.3)

    # Panel B: fitted densities at representative durations
    ax = axes[0, 1]
    ax.hist(x, bins=30, density=True, alpha=0.3, color="gray", edgecolor="black", linewidth=0.5)
    dt_examples = [60, 180, 365]
    colors_dt = ["#e41a1c", "#ff7f00", "#984ea3"]
    for dt_val, color in zip(dt_examples, colors_dt):
        mu2_eff = p3["mu2"] - p3["beta"] * np.log2(dt_val)
        pdf_n = p3["pi_noise"] * stats.norm.pdf(xgrid, p3["mu1"], p3["sigma1"])
        pdf_s = (1 - p3["pi_noise"]) * stats.norm.pdf(xgrid, mu2_eff, p3["sigma2"])
        ax.plot(xgrid, pdf_n + pdf_s, color=color, linewidth=1.5, label=f"Δt={dt_val}d")
    ax.set_xlabel("log2(fold change)")
    ax.set_ylabel("Density")
    ax.set_title("Conditional mixture density at representative durations")
    ax.legend(fontsize=9)
    ax.set_xlim(xgrid[0], xgrid[-1])

    # Panel C: scatter colored by P(responder), x=starting EU
    ax = axes[1, 0]
    sc = ax.scatter(df["eu_start"], x, c=p_resp, cmap="RdYlBu_r",
                    s=25, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    ax.set_xscale("log")
    ax.set_xlabel("Starting Vi IgG (EU)")
    ax.set_ylabel("log2(fold change)")
    ax.set_title("Component assignment vs starting antibody level")
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.grid(True, alpha=0.3)

    # Panel D: posterior assignment histogram
    ax = axes[1, 1]
    ax.hist(p_resp, bins=30, color="#9133be", alpha=0.7, edgecolor="black")
    ax.axvline(0.5, color="red", linestyle="--", linewidth=1.5)
    ax.set_xlabel("P(responder)")
    ax.set_ylabel("Count")
    n_resp = (p_resp > 0.5).sum()
    ax.set_title(f"Posterior assignment distribution\n"
                 f"{n_resp}/{len(x)} responders (P>0.5)")

    cv = sigma_log2_to_cv(p3["sigma1"])
    fig.suptitle(f"Step 3: Mixture Regression (n={len(x)})\n"
                 f"Noise CV≈{cv:.2f}, Signal β={p3['beta']:.4f}, π_noise={p3['pi_noise']:.2f}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step3_mixture_regression.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 08_step3_mixture_regression.png")


def plot_summary(x, p1, p2, p3):
    """Compact summary comparison of all three models."""
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for ax, p, label in [
        (axes[0], p1, "Step 1: Two Gaussians"),
        (axes[1], p2, "Step 2: Gauss + Skew-Normal"),
    ]:
        ax.hist(x, bins=30, density=True, alpha=0.4, color="gray", edgecolor="black", linewidth=0.5)
        pdf_n = p["pi_noise"] * stats.norm.pdf(xgrid, p["mu1"], p["sigma1"])
        if "alpha" in p:
            pdf_s = (1 - p["pi_noise"]) * stats.skewnorm.pdf(xgrid, p["alpha"], p["xi"], p["omega"])
        else:
            pdf_s = (1 - p["pi_noise"]) * stats.norm.pdf(xgrid, p["mu2"], p["sigma2"])
        ax.plot(xgrid, pdf_n, "b-", linewidth=2)
        ax.plot(xgrid, pdf_s, "r-", linewidth=2)
        ax.plot(xgrid, pdf_n + pdf_s, "k--", linewidth=1.5)
        cv = sigma_log2_to_cv(p["sigma1"])
        ax.set_title(f"{label}\nAIC={p['aic']:.1f}, CV≈{cv:.2f}", fontsize=10)
        ax.set_xlabel("log2(fold change)")
        ax.set_xlim(xgrid[0], xgrid[-1])

    # Step 3: show at median duration
    ax = axes[2]
    ax.hist(x, bins=30, density=True, alpha=0.4, color="gray", edgecolor="black", linewidth=0.5)
    for dt_val, color, ls in [(60, "#e41a1c", "-"), (180, "#ff7f00", "-"), (365, "#984ea3", "-")]:
        mu2_eff = p3["mu2"] - p3["beta"] * np.log2(dt_val)
        pdf_n = p3["pi_noise"] * stats.norm.pdf(xgrid, p3["mu1"], p3["sigma1"])
        pdf_s = (1 - p3["pi_noise"]) * stats.norm.pdf(xgrid, mu2_eff, p3["sigma2"])
        ax.plot(xgrid, pdf_n + pdf_s, color=color, linewidth=1.5, label=f"Δt={dt_val}d")
    cv3 = sigma_log2_to_cv(p3["sigma1"])
    ax.set_title(f"Step 3: Mixture Regression\nAIC={p3['aic']:.1f}, CV≈{cv3:.2f}, β={p3['beta']:.4f}",
                 fontsize=10)
    ax.set_xlabel("log2(fold change)")
    ax.legend(fontsize=8)
    ax.set_xlim(xgrid[0], xgrid[-1])

    fig.suptitle(f"Model comparison (n={len(x)})", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_model_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 08_model_comparison.png")


# =========================================================================
# Main
# =========================================================================

def main():
    df = load_fold_change_data()
    x = df["log2_fc"].values
    log2_dt = df["log2_duration"].values

    print(f"Loaded {len(df)} fold-change observations")
    print(f"  log2(FC): median={np.median(x):.3f}, SD={np.std(x):.3f}")
    print(f"  Duration: median={np.median(df['duration_days']):.0f} days, "
          f"range {df['duration_days'].min():.0f}-{df['duration_days'].max():.0f}")

    # =================================================================
    # Step 1: Two-Gaussian mixture
    # =================================================================
    print("\n" + "=" * 60)
    print("STEP 1: Two-Gaussian mixture")
    print("=" * 60)
    p1 = fit_two_gaussian(x)
    cv1 = sigma_log2_to_cv(p1["sigma1"])
    print(f"  Noise:  μ={p1['mu1']:.4f}, σ={p1['sigma1']:.4f}, π={p1['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv1:.3f}")
    print(f"  Signal: μ={p1['mu2']:.4f}, σ={p1['sigma2']:.4f}, π={1-p1['pi_noise']:.3f}")
    print(f"  LL={p1['ll']:.2f}, AIC={p1['aic']:.2f}, BIC={p1['bic']:.2f}")
    plot_step1(x, p1)

    # =================================================================
    # Step 2: Gaussian + Skew-Normal
    # =================================================================
    print("\n" + "=" * 60)
    print("STEP 2: Gaussian + Skew-Normal")
    print("=" * 60)
    p2 = fit_gauss_skewnorm(x, p1)
    cv2 = sigma_log2_to_cv(p2["sigma1"])
    print(f"  Noise:  μ={p2['mu1']:.4f}, σ={p2['sigma1']:.4f}, π={p2['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv2:.3f}")
    print(f"  Signal: ξ={p2['xi']:.4f}, ω={p2['omega']:.4f}, α={p2['alpha']:.4f}")
    sn_mean = p2["xi"] + p2["omega"] * np.sqrt(2 / np.pi) * p2["alpha"] / np.sqrt(1 + p2["alpha"]**2)
    print(f"          Skew-normal mean = {sn_mean:.4f}")
    print(f"  LL={p2['ll']:.2f}, AIC={p2['aic']:.2f}, BIC={p2['bic']:.2f}")
    print(f"  Converged: {p2['converged']}")
    delta_aic = p1["aic"] - p2["aic"]
    print(f"  ΔAIC vs Step 1: {delta_aic:.2f} ({'Step 2 better' if delta_aic > 0 else 'Step 1 better'})")
    plot_step2(x, p1, p2)

    # =================================================================
    # Step 3: Gaussian + time-dependent Normal
    # =================================================================
    print("\n" + "=" * 60)
    print("STEP 3: Gaussian + time-dependent Normal (power-law waning)")
    print("=" * 60)
    p3 = fit_gauss_timedep(x, log2_dt, p1)
    cv3 = sigma_log2_to_cv(p3["sigma1"])
    print(f"  Noise:  μ={p3['mu1']:.4f}, σ={p3['sigma1']:.4f}, π={p3['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv3:.3f}")
    print(f"  Signal: μ₂={p3['mu2']:.4f}, σ₂={p3['sigma2']:.4f}, β={p3['beta']:.4f}")
    print(f"          β is slope of log2(FC) vs log2(days): power-law exponent")
    for dt_val in [60, 180, 365]:
        mu2_eff = p3["mu2"] - p3["beta"] * np.log2(dt_val)
        print(f"          Signal mean at Δt={dt_val}d: {mu2_eff:.4f} log2(FC) = {2**mu2_eff:.3f}× FC")
    print(f"  LL={p3['ll']:.2f}, AIC={p3['aic']:.2f}, BIC={p3['bic']:.2f}")
    print(f"  Converged: {p3['converged']}")
    delta_aic_31 = p1["aic"] - p3["aic"]
    delta_aic_32 = p2["aic"] - p3["aic"]
    print(f"  ΔAIC vs Step 1: {delta_aic_31:.2f}")
    print(f"  ΔAIC vs Step 2: {delta_aic_32:.2f}")
    plot_step3(df, p3)

    # =================================================================
    # Model comparison
    # =================================================================
    print("\n" + "=" * 60)
    print("MODEL COMPARISON")
    print("=" * 60)
    print(f"{'Model':<35} {'LL':>8} {'AIC':>8} {'BIC':>8} {'CV_est':>7}")
    print("-" * 70)
    for name, p, cv in [
        ("1. Two Gaussians", p1, cv1),
        ("2. Gaussian + Skew-Normal", p2, cv2),
        ("3. Gaussian + time-dep Normal", p3, cv3),
    ]:
        print(f"{name:<35} {p['ll']:>8.2f} {p['aic']:>8.2f} {p['bic']:>8.2f} {cv:>7.3f}")

    # =================================================================
    # Component assignments summary (best model)
    # =================================================================
    best = min([(p1, "Step 1"), (p2, "Step 2"), (p3, "Step 3")], key=lambda x: x[0]["aic"])
    best_p, best_name = best
    print(f"\n{'='*60}")
    print(f"COMPONENT ASSIGNMENTS ({best_name}, P(responder) > 0.5)")
    print(f"{'='*60}")
    df["p_responder"] = 1 - best_p["p_noise"]
    responders = df[df["p_responder"] > 0.5]
    noise = df[df["p_responder"] <= 0.5]
    print(f"  Responders: {len(responders)}/{len(df)} ({100*len(responders)/len(df):.0f}%)")
    if len(responders) > 0:
        print(f"    Median FC: {2**responders['log2_fc'].median():.3f}")
        print(f"    Median duration: {responders['duration_days'].median():.0f} days")
        print(f"    Median starting EU: {responders['eu_start'].median():.1f}")
    print(f"  Noise: {len(noise)}/{len(df)} ({100*len(noise)/len(df):.0f}%)")
    if len(noise) > 0:
        print(f"    Median FC: {2**noise['log2_fc'].median():.3f}")
        print(f"    Median duration: {noise['duration_days'].median():.0f} days")
        print(f"    Median starting EU: {noise['eu_start'].median():.1f}")

    plot_summary(x, p1, p2, p3)


if __name__ == "__main__":
    main()
