"""
08: Fold-change mixture model analysis for Vi IgG.

Four models of increasing complexity fitted to log2(fold change) data:
  Step 1: Two-Gaussian mixture (MLE)
  Step 2: Gaussian + Skew-Normal mixture (MLE)
  Step 3: Gaussian + time-dependent Normal — log2(Δt) covariate, signal only
  Step 4: Teunis power-law mixture — log2(t1/t0) covariate, both components

All models enforce: sigma2 = sqrt(sigma1^2 + sigma_bio^2)
  — the signal component inherits measurement noise and adds biological variance.

Step 4 physics:
  For Teunis power-law decay y(τ) = y_peak·(1+β·τ)^(-α), with τ = t - t_peak.
  We tested the full form with β estimated, but β diverges to ∞ — the data
  can't separate β from α. In the large-β limit, (1+β·τ)^(-α) → (β·τ)^(-α)
  ∝ τ^(-α), so the fold change reduces to:
      FC = (τ₁/τ₀)^(-α)
      log2(FC) = -α · log2(τ₁/τ₀)
  where τ = max(t - t_peak, 1 day).

  We use this power-law form directly with:
    - t_peak = 15 days (from HlyE IgG, see comment in code)
    - τ floored at 1 day for pre-peak observations
    - FC = 1 at τ₁ = τ₀ automatic (no intercept needed for noise component)
    - α is the Teunis power-law waning exponent
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, optimize
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")


def load_fold_change_data():
    """Load longitudinal Vi IgG data and compute fold changes.

    Uses first two observations (consecutive visits) rather than first-to-last.
    This captures the acute-phase response more cleanly: the first two visits
    are closest to infection and most informative about boosting, rather than
    conflating initial boost with long-term waning in 3-4 point subjects.
    """
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    records = []
    for sid, grp in long[long["has_longitudinal"]].groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        t0, eu0 = grp.iloc[0]["days_since_fever_onset"], grp.iloc[0]["vi_igg_eu"]
        t1, eu1 = grp.iloc[1]["days_since_fever_onset"], grp.iloc[1]["vi_igg_eu"]
        dt = t1 - t0
        if dt <= 0 or eu0 <= 0 or eu1 <= 0:
            continue
        log2_fc = np.log2(eu1 / eu0)
        # Filter indx387: 291→0.6 EU in 95 days (500× decline), almost
        # certainly a lab error. No Vi IgG biological process does this.
        if abs(log2_fc) > 6:
            continue
        records.append({
            "index_id": sid,
            "log2_fc": np.log2(eu1 / eu0),
            "duration_days": dt,
            "log2_duration": np.log2(dt),
            "t_start": t0,
            "t_end": t1,
            "log2_time_ratio": np.log2(t1 / t0),  # Teunis power-law covariate
            "eu_start": eu0,
        })
    return pd.DataFrame(records)


def sigma_log2_to_cv(sigma_log2):
    """Convert SD on log2 scale (for fold-change of two measurements) to natural-scale CV."""
    sigma_per_meas = sigma_log2 / np.sqrt(2)
    sigma_ln = sigma_per_meas * np.log(2)
    return np.sqrt(np.exp(sigma_ln**2) - 1)


def compute_sigma2(sigma1, sigma_bio):
    """Signal SD = sqrt(measurement^2 + biological^2)."""
    return np.sqrt(sigma1**2 + sigma_bio**2)


def compute_posteriors(x, pi_noise, pdf_noise, pdf_signal):
    """Compute posterior P(noise) for each observation."""
    mixture = pi_noise * pdf_noise + (1 - pi_noise) * pdf_signal
    mixture = np.maximum(mixture, 1e-300)
    return (pi_noise * pdf_noise) / mixture, mixture


def compute_ic(ll, n_params, n):
    """AIC and BIC."""
    return 2 * n_params - 2 * ll, n_params * np.log(n) - 2 * ll


# =========================================================================
# Step 1: Two-Gaussian mixture (MLE)
# =========================================================================

def neg_ll_two_gauss(theta, x):
    # Reparametrize mu2 = mu1 + exp(log_delta) to enforce mu2 > mu1 during
    # optimization. This resolves the label-switching ambiguity at the
    # model level rather than by post-hoc relabeling, which would break
    # the asymmetric sigma2 = sqrt(sigma1^2 + sigma_bio^2) constraint.
    mu1, log_sigma1, log_delta, log_sigma_bio, logit_pi = theta
    mu2 = mu1 + np.exp(log_delta)
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    pdf_n = stats.norm.pdf(x, mu1, sigma1)
    pdf_s = stats.norm.pdf(x, mu2, sigma2)
    mixture = pi_noise * pdf_n + (1 - pi_noise) * pdf_s
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def fit_two_gaussian(x):
    """Fit 2-component Gaussian mixture via MLE with sigma2 = sqrt(sigma1^2 + sigma_bio^2)."""
    # Multi-start optimization
    best_ll = -np.inf
    best_result = None
    rng = np.random.default_rng(42)

    for _ in range(20):
        mu1_init = rng.normal(np.median(x) - 0.3, 0.3)
        delta_init = max(rng.normal(0.8, 0.3), 0.1)
        x0 = [mu1_init, np.log(0.4), np.log(delta_init), np.log(0.8), 0.0]

        result = optimize.minimize(
            neg_ll_two_gauss, x0, args=(x,),
            method="Nelder-Mead", options={"maxiter": 50000, "xatol": 1e-8, "fatol": 1e-8},
        )
        if -result.fun > best_ll:
            best_ll = -result.fun
            best_result = result

    mu1, log_sigma1, log_delta, log_sigma_bio, logit_pi = best_result.x
    mu2 = mu1 + np.exp(log_delta)
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    n_params = 5
    aic, bic = compute_ic(best_ll, n_params, len(x))

    pdf_n = stats.norm.pdf(x, mu1, sigma1)
    pdf_s = stats.norm.pdf(x, mu2, sigma2)
    p_noise, _ = compute_posteriors(x, pi_noise, pdf_n, pdf_s)

    return {
        "mu1": mu1, "sigma1": sigma1, "mu2": mu2,
        "sigma_bio": sigma_bio, "sigma2": sigma2,
        "pi_noise": pi_noise,
        "ll": best_ll, "aic": aic, "bic": bic,
        "p_noise": p_noise,
    }


# =========================================================================
# Step 2: Gaussian + Skew-Normal mixture
# =========================================================================

def neg_ll_gauss_skewnorm(theta, x):
    mu1, log_sigma1, xi, log_sigma_bio, alpha, logit_pi = theta
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    omega = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    pdf_n = stats.norm.pdf(x, mu1, sigma1)
    pdf_s = stats.skewnorm.pdf(x, alpha, xi, omega)
    mixture = pi_noise * pdf_n + (1 - pi_noise) * pdf_s
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def fit_gauss_skewnorm(x, p1):
    x0 = [
        p1["mu1"], np.log(p1["sigma1"]),
        p1["mu2"], np.log(p1["sigma_bio"]),
        0.0,  # alpha=0
        np.log(p1["pi_noise"] / (1 - p1["pi_noise"])),
    ]
    result = optimize.minimize(
        neg_ll_gauss_skewnorm, x0, args=(x,),
        method="Nelder-Mead", options={"maxiter": 50000, "xatol": 1e-8, "fatol": 1e-8},
    )

    mu1, log_sigma1, xi, log_sigma_bio, alpha, logit_pi = result.x
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    omega = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))
    ll = -result.fun
    aic, bic = compute_ic(ll, 6, len(x))

    pdf_n = stats.norm.pdf(x, mu1, sigma1)
    pdf_s = stats.skewnorm.pdf(x, alpha, xi, omega)
    p_noise, _ = compute_posteriors(x, pi_noise, pdf_n, pdf_s)

    return {
        "mu1": mu1, "sigma1": sigma1,
        "xi": xi, "omega": omega, "alpha": alpha,
        "sigma_bio": sigma_bio,
        "pi_noise": pi_noise,
        "ll": ll, "aic": aic, "bic": bic,
        "converged": result.success, "p_noise": p_noise,
    }


# =========================================================================
# Step 3: Gaussian + time-dependent Normal (signal only)
# =========================================================================

def neg_ll_gauss_timedep(theta, x, log2_dt):
    mu1, log_sigma1, mu2, log_sigma_bio, beta, logit_pi = theta
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    pdf_n = stats.norm.pdf(x, mu1, sigma1)
    pdf_s = stats.norm.pdf(x, mu2 - beta * log2_dt, sigma2)
    mixture = pi_noise * pdf_n + (1 - pi_noise) * pdf_s
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def fit_gauss_timedep(x, log2_dt, p1):
    x0 = [
        p1["mu1"], np.log(p1["sigma1"]),
        p1["mu2"], np.log(p1["sigma_bio"]),
        0.0,  # beta=0
        np.log(p1["pi_noise"] / (1 - p1["pi_noise"])),
    ]
    result = optimize.minimize(
        neg_ll_gauss_timedep, x0, args=(x, log2_dt),
        method="Nelder-Mead", options={"maxiter": 50000, "xatol": 1e-8, "fatol": 1e-8},
    )

    mu1, log_sigma1, mu2, log_sigma_bio, beta, logit_pi = result.x
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))
    ll = -result.fun
    aic, bic = compute_ic(ll, 6, len(x))

    pdf_n = stats.norm.pdf(x, mu1, sigma1)
    pdf_s = stats.norm.pdf(x, mu2 - beta * log2_dt, sigma2)
    p_noise, _ = compute_posteriors(x, pi_noise, pdf_n, pdf_s)

    return {
        "mu1": mu1, "sigma1": sigma1,
        "mu2": mu2, "sigma_bio": sigma_bio, "sigma2": sigma2,
        "beta": beta, "pi_noise": pi_noise,
        "ll": ll, "aic": aic, "bic": bic,
        "converged": result.success, "p_noise": p_noise,
    }


# =========================================================================
# Step 4: Teunis power-law mixture (both components, τ-ratio covariate)
#
# Power-law decay in the large-βτ limit: y(τ) ∝ τ^(-α)
# Fold change: FC = (τ₁/τ₀)^(-α), log2(FC) = -α * log2(τ₁/τ₀)
# where τ = max(t - t_peak, 1 day).
#
# We tested the full Teunis form (1+β·τ)^(-α) with β estimated, but
# β diverges to ∞ — the data can't distinguish β from α in the fold-
# change. The full form collapses to the power-law, so we use it directly.
#
# t_peak = 15 days post fever onset (fixed).
#   Rationale: Aiemjoy Table S1 reports t1 (time to peak, days since fever
#   onset) for the well-identified antigens: HlyE IgA=19.6d, HlyE IgG=15.6d,
#   LPS IgA=11.1d, LPS IgG=5.5d. Vi IgG t1=2.9d is unreliable (model can't
#   identify the Vi peak). HlyE IgG (15.6d) has the clearest waning signal
#   and the best-identified parameters in the same subjects/platform. We use
#   15d as a round number. This is days from fever onset, which is roughly
#   coincident with or a few days after bacteremia onset.
#
# τ is floored at 1 day. Observations before t_peak (e.g., t=0 on fever
# onset day) are treated as "1 day post-peak" rather than τ=0, which
# would cause log2(τ₁/τ₀) to diverge.
#
# Noise component: N(-α₁ * log2(τ₁/τ₀), σ₁)
#   - No intercept: FC = 1 when τ₁ = τ₀ (automatic)
#   - α₁ is the power-law waning exponent for the non-boosted population
#
# Signal component: B₀ ~ TruncNorm(μ_bio, σ_bio, lower=0), then waned + noise
#   - B₀ is the initial boost at t_peak (always ≥ 0 on log2 scale)
#   - Observed: X = B₀ - α₂·z + ε, where ε ~ N(0, σ₁)
#   - So (X + α₂·z) ~ TruncNorm(μ_bio, σ_bio, ≥0) ⊛ N(0, σ₁)
#   - The truncation is on the INITIAL BOOST, not the observed FC.
#     Subjects who boosted then waned past baseline (FC < 1) are allowed.
#   - Closed-form: evaluate truncnorm_conv_pdf at (x + α₂·z) with mean μ_bio
#   - μ_bio is the mean of the un-truncated initial boost distribution
#   - α₂ captures waning in the boosted population
#
# 6 parameters: α₁, σ₁, μ_bio, α₂, σ_bio, π
# (t_peak fixed at 15 days)
# =========================================================================

T_PEAK = 15.0  # days since fever onset; see rationale above


def tau_ratio_covariate(t0, t1, tp=T_PEAK):
    """Compute z = log2(τ₁/τ₀) where τ = max(t - tp, 1).

    The 1-day floor avoids log2(0) for subjects observed at or before
    t_peak and prevents the degenerate β→∞ solution seen with the full
    Teunis form.
    """
    tau0 = np.maximum(t0 - tp, 1.0)
    tau1 = np.maximum(t1 - tp, 1.0)
    return np.log2(tau1 / tau0)


def truncnorm_conv_pdf(x, mu, sigma1, sigma_bio):
    """PDF of TruncNorm(mu, sigma_bio, lower=0) convolved with N(0, sigma1).

    The biological fold change B ~ TruncNorm(mu, sigma_bio, lower=0) represents
    a boost that is always ≥ 0 on the log2 scale (FC ≥ 1). Observation noise
    ε ~ N(0, sigma1) is added on top: X = B + ε.

    Closed form (derivation: complete the square in the convolution integral):
      f(x) = N(x; mu, sigma2) · Φ(m/v) / (1 - Φ(-mu/sigma_bio))
    where:
      sigma2 = sqrt(sigma1^2 + sigma_bio^2)
      m = (x·sigma_bio^2 + mu·sigma1^2) / sigma2^2
      v = sigma1·sigma_bio / sigma2
      Φ(-mu/sigma_bio) normalizes the truncated normal

    Limits:
      - When mu >> 0, truncation is irrelevant → recovers N(x; mu, sigma2)
      - When x → -∞, Φ(m/v) → 0 → f(x) → 0 (left tail dies, gives monotonic P(resp))
    """
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    m = (x * sigma_bio**2 + mu * sigma1**2) / sigma2**2
    v = sigma1 * sigma_bio / sigma2

    trunc_norm = 1 - stats.norm.cdf(-mu / sigma_bio)
    # Floor at 0.01: don't allow >99% of bio mass below 0.
    # When mu << 0, nearly all bio mass is truncated away, making the
    # PDF blow up. This floor prevents degenerate solutions where the
    # optimizer pushes mu_bio → -∞ to exploit the singularity.
    trunc_norm = np.maximum(trunc_norm, 0.01)

    return stats.norm.pdf(x, mu, sigma2) * stats.norm.cdf(m / v) / trunc_norm


def neg_ll_teunis_untrunc(theta, x, z):
    """Untruncated Gaussian signal — the working model (Step 4a)."""
    alpha1, log_sigma1, mu_bio, alpha2, log_sigma_bio, logit_pi = theta
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    pdf_n = stats.norm.pdf(x, -alpha1 * z, sigma1)
    pdf_s = stats.norm.pdf(x, mu_bio - alpha2 * z, sigma2)

    mixture = pi_noise * pdf_n + (1 - pi_noise) * pdf_s
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def neg_ll_teunis_mixture(theta, x, z):
    """Truncated signal — experimental (Step 4b)."""
    alpha1, log_sigma1, mu_bio, alpha2, log_sigma_bio, logit_pi = theta
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    pdf_n = stats.norm.pdf(x, -alpha1 * z, sigma1)
    pdf_s = truncnorm_conv_pdf(x + alpha2 * z, mu_bio, sigma1, sigma_bio)

    mixture = pi_noise * pdf_n + (1 - pi_noise) * pdf_s
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def fit_teunis_untrunc(x, t0, t1, p1):
    """Fit untruncated Teunis mixture — the WORKING MODEL (Step 4a).

    Both components Gaussian, τ-ratio covariate, σ₂ = √(σ₁² + σ_bio²).
    """
    z = tau_ratio_covariate(t0, t1)
    x0 = [
        0.025, np.log(p1["sigma1"]),
        p1["mu2"], -0.1,
        np.log(p1["sigma_bio"]),
        np.log(p1["pi_noise"] / (1 - p1["pi_noise"])),
    ]
    result = optimize.minimize(
        neg_ll_teunis_untrunc, x0, args=(x, z),
        method="Nelder-Mead", options={"maxiter": 100000, "xatol": 1e-8, "fatol": 1e-8},
    )

    alpha1, log_sigma1, mu_bio, alpha2, log_sigma_bio, logit_pi = result.x
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))
    ll = -result.fun
    aic, bic = compute_ic(ll, 6, len(x))

    pdf_n = stats.norm.pdf(x, -alpha1 * z, sigma1)
    pdf_s = stats.norm.pdf(x, mu_bio - alpha2 * z, sigma2)
    p_noise, _ = compute_posteriors(x, pi_noise, pdf_n, pdf_s)

    return {
        "alpha1": alpha1, "sigma1": sigma1,
        "mu_bio": mu_bio, "alpha2": alpha2,
        "sigma_bio": sigma_bio, "sigma2": sigma2,
        "t_peak": T_PEAK, "pi_noise": pi_noise,
        "ll": ll, "aic": aic, "bic": bic,
        "converged": result.success, "p_noise": p_noise,
    }


def fit_teunis_mixture(x, t0, t1, p1):
    """Fit truncated Teunis mixture — experimental (Step 4b).

    Multi-start optimization to avoid local minima — the truncated signal
    component creates a more complex likelihood surface than the Gaussian version.
    """
    z = tau_ratio_covariate(t0, t1)
    rng = np.random.default_rng(42)
    best_ll = -np.inf
    best_result = None

    # Start from several initializations
    inits = [
        # From Step 1 estimates
        [0.1, np.log(p1["sigma1"]), p1["mu2"], 0.1, np.log(p1["sigma_bio"]),
         np.log(p1["pi_noise"] / (1 - p1["pi_noise"]))],
        # Near untruncated Step 4 solution (the target neighborhood)
        [0.025, np.log(0.43), -0.09, -0.11, np.log(0.86), np.log(0.38 / 0.62)],
        # Variations around untruncated solution
        [0.05, np.log(0.40), 0.0, -0.05, np.log(0.90), -0.3],
        [0.02, np.log(0.45), -0.2, -0.15, np.log(0.80), -0.5],
        [0.03, np.log(0.40), 0.2, 0.0, np.log(0.70), 0.5],
        # Small σ_bio (truncation matters more)
        [0.03, np.log(0.40), 0.5, -0.10, np.log(0.50), -0.3],
        [0.03, np.log(0.40), 0.3, -0.05, np.log(0.60), 0.0],
    ]
    # Add random perturbations
    for _ in range(16):
        base = inits[rng.integers(len(inits))].copy()
        perturbed = [b + rng.normal(0, 0.3) for b in base]
        inits.append(perturbed)

    for x0 in inits:
        result = optimize.minimize(
            neg_ll_teunis_mixture, x0, args=(x, z),
            method="Nelder-Mead", options={"maxiter": 100000, "xatol": 1e-8, "fatol": 1e-8},
        )
        if -result.fun > best_ll:
            best_ll = -result.fun
            best_result = result

    result = best_result

    alpha1, log_sigma1, mu_bio, alpha2, log_sigma_bio, logit_pi = result.x
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))
    ll = -result.fun
    aic, bic = compute_ic(ll, 6, len(x))

    pdf_n = stats.norm.pdf(x, -alpha1 * z, sigma1)
    pdf_s = truncnorm_conv_pdf(x + alpha2 * z, mu_bio, sigma1, sigma_bio)
    p_noise, _ = compute_posteriors(x, pi_noise, pdf_n, pdf_s)

    return {
        "alpha1": alpha1, "sigma1": sigma1,
        "mu_bio": mu_bio, "alpha2": alpha2,
        "sigma_bio": sigma_bio, "sigma2": sigma2,
        "t_peak": T_PEAK, "pi_noise": pi_noise,
        "ll": ll, "aic": aic, "bic": bic,
        "converged": result.success, "p_noise": p_noise,
    }


# =========================================================================
# Step 5: Fold-rise-informed mixture (Stage A)
#
# Replaces the free μ_bio with the fold-rise model from
# scratch/cohort_incidence_model_proof_of_concept.R:
#
#   fold_rise = 10^(μ₀ * (1 - (log10(CoP_pre) - log10(CoP_min))
#                              / (log10(CoP_max) - log10(CoP_min))))
#
# On the log2 scale, the expected initial boost for a responder with
# starting EU = eu_start is:
#   B₀_mean(eu_start) = log2(fold_rise(eu_start))
#                      = μ₀/log10(2) * (1 - log10(eu_start)/log10(CoP_max))
#   (using CoP_min = 1, so log10(CoP_min) = 0)
#
# The signal component becomes:
#   B₀ ~ TruncNorm(B₀_mean(eu_start), σ_bio, lower=0), then waned + noise
#   Observed: X = B₀ - α₂·z + ε
#
# Noise component unchanged: N(-α₁·z, σ₁)
#
# Parameters: μ₀ (boost intensity), log10(CoP_max), α₁, α₂, σ₁, σ_bio, π
# 7 parameters. The fold-rise structure predicts the eu_start-FC correlation
# that was ρ=-0.62 in the data.
# =========================================================================

COP_MIN = 1.0  # EU/ml, baseline for unexposed (below LOD)


def fold_rise_log2(eu_start, mu_0, log10_cop_max):
    """Expected log2(fold rise) from the fold-rise model.

    fold_rise = 10^(μ₀ * (1 - log10(eu_start) / log10(CoP_max)))
    log2(fold_rise) = μ₀ / log10(2) * (1 - log10(eu_start) / log10(CoP_max))

    Clamped to ≥ 0: subjects at or above CoP_max get fold_rise = 1 (no boost).
    """
    frac = np.log10(np.maximum(eu_start, COP_MIN)) / log10_cop_max
    return np.maximum(mu_0 / np.log10(2) * (1 - frac), 0.0)


def neg_ll_foldrise_mixture(theta, x, z, eu_start):
    alpha1, log_sigma1, mu_0, log10_cop_max, alpha2, log_sigma_bio, logit_pi = theta
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))

    # Noise: waning only
    pdf_n = stats.norm.pdf(x, -alpha1 * z, sigma1)

    # Signal: fold-rise-predicted boost, truncated at 0, then waned + noise
    b0_mean = fold_rise_log2(eu_start, mu_0, log10_cop_max)
    pdf_s = truncnorm_conv_pdf(x + alpha2 * z, b0_mean, sigma1, sigma_bio)

    mixture = pi_noise * pdf_n + (1 - pi_noise) * pdf_s
    mixture = np.maximum(mixture, 1e-300)
    return -np.sum(np.log(mixture))


def fit_foldrise_mixture(x, t0, t1, eu_start, p1):
    """Fit fold-rise-informed mixture. Uses τ-ratio covariate + eu_start."""
    z = tau_ratio_covariate(t0, t1)
    rng = np.random.default_rng(42)
    best_ll = -np.inf
    best_result = None

    inits = [
        # From cohort model defaults: mu_0=2.5, CoP_max=10^3.5
        [0.025, np.log(0.43), 2.5, 3.5, -0.11, np.log(0.5), np.log(0.38 / 0.62)],
        # Variations
        [0.03, np.log(0.40), 2.0, 3.0, -0.05, np.log(0.6), -0.3],
        [0.02, np.log(0.45), 3.0, 3.5, -0.15, np.log(0.4), -0.5],
        [0.01, np.log(0.40), 1.5, 3.0, 0.0, np.log(0.8), 0.0],
        [0.05, np.log(0.40), 2.5, 4.0, -0.10, np.log(0.5), 0.3],
    ]
    for _ in range(20):
        base = inits[rng.integers(len(inits))].copy()
        perturbed = [b + rng.normal(0, 0.2) for b in base]
        inits.append(perturbed)

    for x0 in inits:
        result = optimize.minimize(
            neg_ll_foldrise_mixture, x0, args=(x, z, eu_start),
            method="Nelder-Mead", options={"maxiter": 100000, "xatol": 1e-8, "fatol": 1e-8},
        )
        if -result.fun > best_ll:
            best_ll = -result.fun
            best_result = result

    result = best_result
    alpha1, log_sigma1, mu_0, log10_cop_max, alpha2, log_sigma_bio, logit_pi = result.x
    sigma1 = np.exp(log_sigma1)
    sigma_bio = np.exp(log_sigma_bio)
    sigma2 = compute_sigma2(sigma1, sigma_bio)
    pi_noise = 1 / (1 + np.exp(-logit_pi))
    ll = -result.fun
    aic, bic = compute_ic(ll, 7, len(x))

    b0_mean = fold_rise_log2(eu_start, mu_0, log10_cop_max)
    pdf_n = stats.norm.pdf(x, -alpha1 * z, sigma1)
    pdf_s = truncnorm_conv_pdf(x + alpha2 * z, b0_mean, sigma1, sigma_bio)
    p_noise, _ = compute_posteriors(x, pi_noise, pdf_n, pdf_s)

    return {
        "alpha1": alpha1, "sigma1": sigma1,
        "mu_0": mu_0, "log10_cop_max": log10_cop_max,
        "cop_max": 10**log10_cop_max,
        "alpha2": alpha2, "sigma_bio": sigma_bio, "sigma2": sigma2,
        "t_peak": T_PEAK, "pi_noise": pi_noise,
        "ll": ll, "aic": aic, "bic": bic,
        "converged": result.success, "p_noise": p_noise,
        "b0_mean": b0_mean,  # per-subject predicted boost
    }


# =========================================================================
# Plotting — one figure per step + summary
# =========================================================================

def plot_step1(x, p):
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    ax.hist(x, bins=30, density=True, alpha=0.4, color="gray", edgecolor="black", linewidth=0.5)
    pdf_n = p["pi_noise"] * stats.norm.pdf(xgrid, p["mu1"], p["sigma1"])
    pdf_s = (1 - p["pi_noise"]) * stats.norm.pdf(xgrid, p["mu2"], p["sigma2"])
    ax.plot(xgrid, pdf_n, "b-", lw=2, label=f"Noise (π={p['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_s, "r-", lw=2, label=f"Signal (π={1-p['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_n + pdf_s, "k--", lw=1.5, label="Mixture")
    ax.axvline(0, color="gray", linestyle=":", lw=0.8)
    ax.set_xlabel("log2(fold change)"); ax.set_ylabel("Density")
    cv = sigma_log2_to_cv(p["sigma1"])
    ax.set_title(f"Two-Gaussian mixture\nNoise: μ₁={p['mu1']:.3f}, σ₁={p['sigma1']:.3f} (CV≈{cv:.2f})\n"
                 f"Signal: μ₂={p['mu2']:.3f}, σ₂={p['sigma2']:.3f} (σ_bio={p['sigma_bio']:.3f})")
    ax.legend(fontsize=9); ax.set_xlim(xgrid[0], xgrid[-1])
    ax = axes[1]
    p_resp = 1 - p["p_noise"]
    ax.scatter(x, p_resp, c=p_resp, cmap="RdYlBu_r", s=20, alpha=0.7,
               edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    ax.axhline(0.5, color="gray", linestyle="--", lw=0.8)
    ax.set_xlabel("log2(fold change)"); ax.set_ylabel("P(responder)")
    ax.set_title(f"Posterior assignment\n{(p_resp > 0.5).sum()}/{len(x)} responders (P>0.5)")
    ax.grid(True, alpha=0.3)
    fig.suptitle("Step 1: Two-Gaussian Mixture (σ₂ = √(σ₁² + σ_bio²))", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step1_two_gaussian.png", dpi=150, bbox_inches="tight"); plt.close()
    print("Saved 08_step1_two_gaussian.png")


def plot_step2(x, p1, p2):
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    ax.hist(x, bins=30, density=True, alpha=0.4, color="gray", edgecolor="black", linewidth=0.5)
    mix1 = (p1["pi_noise"] * stats.norm.pdf(xgrid, p1["mu1"], p1["sigma1"]) +
            (1 - p1["pi_noise"]) * stats.norm.pdf(xgrid, p1["mu2"], p1["sigma2"]))
    ax.plot(xgrid, mix1, "k:", lw=1.5, alpha=0.5, label="Step 1")
    pdf_n = p2["pi_noise"] * stats.norm.pdf(xgrid, p2["mu1"], p2["sigma1"])
    pdf_s = (1 - p2["pi_noise"]) * stats.skewnorm.pdf(xgrid, p2["alpha"], p2["xi"], p2["omega"])
    ax.plot(xgrid, pdf_n, "b-", lw=2, label=f"Noise (π={p2['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_s, "r-", lw=2, label=f"Signal (π={1-p2['pi_noise']:.2f})")
    ax.plot(xgrid, pdf_n + pdf_s, "k--", lw=1.5, label="Mixture")
    ax.set_xlabel("log2(fold change)"); ax.set_ylabel("Density")
    cv = sigma_log2_to_cv(p2["sigma1"])
    ax.set_title(f"Gaussian + Skew-Normal\nσ₁={p2['sigma1']:.3f} (CV≈{cv:.2f}), "
                 f"ω={p2['omega']:.3f}, α={p2['alpha']:.2f}")
    ax.legend(fontsize=8); ax.set_xlim(xgrid[0], xgrid[-1])
    ax = axes[1]
    ax.scatter(1 - p1["p_noise"], 1 - p2["p_noise"], c=x, cmap="coolwarm", s=20, alpha=0.7,
               edgecolors="black", linewidths=0.3, vmin=-2, vmax=2)
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("P(responder) — Step 1"); ax.set_ylabel("P(responder) — Step 2")
    ax.set_title("Assignment comparison (color = log2 FC)")
    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05); ax.grid(True, alpha=0.3)
    delta = p1["aic"] - p2["aic"]
    ax.text(0.05, 0.95, f"ΔAIC = {delta:.1f}", transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    fig.suptitle("Step 2: Gaussian + Skew-Normal (ω = √(σ₁² + σ_bio²))", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step2_skew_normal.png", dpi=150, bbox_inches="tight"); plt.close()
    print("Saved 08_step2_skew_normal.png")


def plot_step3(df, p3):
    x = df["log2_fc"].values
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    p_resp = 1 - p3["p_noise"]
    ax = axes[0, 0]
    sc = ax.scatter(df["duration_days"], x, c=p_resp, cmap="RdYlBu_r",
                    s=25, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    dt_grid = np.array([30, 60, 90, 180, 365, 730, 1100])
    mu2_grid = p3["mu2"] - p3["beta"] * np.log2(dt_grid)
    ax.plot(dt_grid, mu2_grid, "r-", lw=2, label="Signal mean")
    ax.fill_between(dt_grid, mu2_grid - p3["sigma2"], mu2_grid + p3["sigma2"], color="red", alpha=0.15)
    ax.axhline(p3["mu1"], color="blue", linestyle="--", lw=1, label="Noise mean")
    ax.set_xlabel("Duration (days)"); ax.set_ylabel("log2(fold change)"); ax.set_xscale("log")
    ax.set_title(f"Signal β={p3['beta']:.4f}"); ax.legend(fontsize=7, loc="upper right"); ax.grid(True, alpha=0.3)
    ax = axes[0, 1]
    ax.hist(x, bins=30, density=True, alpha=0.3, color="gray", edgecolor="black", linewidth=0.5)
    for dt_val, color in [(60, "#e41a1c"), (180, "#ff7f00"), (365, "#984ea3")]:
        mu2_eff = p3["mu2"] - p3["beta"] * np.log2(dt_val)
        pdf_n = p3["pi_noise"] * stats.norm.pdf(xgrid, p3["mu1"], p3["sigma1"])
        pdf_s = (1 - p3["pi_noise"]) * stats.norm.pdf(xgrid, mu2_eff, p3["sigma2"])
        ax.plot(xgrid, pdf_n + pdf_s, color=color, lw=1.5, label=f"Δt={dt_val}d")
    ax.set_xlabel("log2(fold change)"); ax.set_ylabel("Density")
    ax.set_title("Conditional density"); ax.legend(fontsize=9); ax.set_xlim(xgrid[0], xgrid[-1])
    ax = axes[1, 0]
    sc = ax.scatter(df["eu_start"], x, c=p_resp, cmap="RdYlBu_r",
                    s=25, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)"); ax.set_xscale("log")
    ax.set_xlabel("Starting Vi IgG (EU)"); ax.set_ylabel("log2(fold change)")
    ax.set_title("Assignment vs starting antibody"); ax.axhline(0, color="gray", linestyle=":", lw=0.8); ax.grid(True, alpha=0.3)
    ax = axes[1, 1]
    ax.hist(p_resp, bins=30, color="#9133be", alpha=0.7, edgecolor="black")
    ax.axvline(0.5, color="red", linestyle="--", lw=1.5)
    ax.set_xlabel("P(responder)"); ax.set_ylabel("Count")
    ax.set_title(f"{(p_resp > 0.5).sum()}/{len(x)} responders (P>0.5)")
    cv = sigma_log2_to_cv(p3["sigma1"])
    fig.suptitle(f"Step 3: Signal Time-Dependent (n={len(x)})\nCV≈{cv:.2f}, σ_bio={p3['sigma_bio']:.3f}, β={p3['beta']:.4f}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step3_mixture_regression.png", dpi=150, bbox_inches="tight"); plt.close()
    print("Saved 08_step3_mixture_regression.png")


def plot_step4(df, p4):
    x = df["log2_fc"].values
    t0 = df["t_start"].values
    t1 = df["t_end"].values
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    z = tau_ratio_covariate(t0, t1)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    p_resp = 1 - p4["p_noise"]

    # Panel A: scatter vs Teunis covariate z
    ax = axes[0, 0]
    sc = ax.scatter(z, x, c=p_resp, cmap="RdYlBu_r",
                    s=25, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    z_grid = np.linspace(z.min(), z.max(), 100)
    ax.plot(z_grid, -p4["alpha1"] * z_grid, "b-", lw=2, label=f"Noise (α₁={p4['alpha1']:.3f})")
    ax.fill_between(z_grid, -p4["alpha1"] * z_grid - p4["sigma1"],
                     -p4["alpha1"] * z_grid + p4["sigma1"], color="blue", alpha=0.1)
    ax.plot(z_grid, p4["mu_bio"] - p4["alpha2"] * z_grid, "r-", lw=2,
            label=f"Signal (α₂={p4['alpha2']:.3f})")
    ax.fill_between(z_grid, p4["mu_bio"] - p4["alpha2"] * z_grid - p4["sigma2"],
                     p4["mu_bio"] - p4["alpha2"] * z_grid + p4["sigma2"], color="red", alpha=0.1)
    ax.set_xlabel(f"Teunis covariate z = log2(τ₁/τ₀)  [tp={p4['t_peak']:.0f}d, τ=max(t-tp, 1)]")
    ax.set_ylabel("log2(fold change)")
    ax.set_title("Teunis power-law mixture"); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel B: scatter vs duration (for comparison)
    ax = axes[0, 1]
    sc = ax.scatter(df["duration_days"], x, c=p_resp, cmap="RdYlBu_r",
                    s=25, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    ax.set_xlabel("Duration Δt (days)"); ax.set_ylabel("log2(fold change)"); ax.set_xscale("log")
    ax.set_title("Same assignments vs duration"); ax.axhline(0, color="gray", linestyle=":", lw=0.8)
    ax.grid(True, alpha=0.3)

    # Panel C: component densities at a representative (t0, t1)
    ax = axes[1, 0]
    ax.hist(x, bins=30, density=True, alpha=0.3, color="gray", edgecolor="black", linewidth=0.5)
    examples = [(14, 100), (14, 365), (30, 365)]
    colors_ex = ["#e41a1c", "#ff7f00", "#984ea3"]
    for (t0_ex, t1_ex), color in zip(examples, colors_ex):
        z_ex = tau_ratio_covariate(t0_ex, t1_ex)
        mu_n = -p4["alpha1"] * z_ex
        pdf_n = p4["pi_noise"] * stats.norm.pdf(xgrid, mu_n, p4["sigma1"])
        # Signal: truncation on initial boost, then waned by alpha2*z
        pdf_s = (1 - p4["pi_noise"]) * truncnorm_conv_pdf(
            xgrid + p4["alpha2"] * z_ex, p4["mu_bio"], p4["sigma1"], p4["sigma_bio"])
        ax.plot(xgrid, pdf_n + pdf_s, color=color, lw=1.5, label=f"t₀={t0_ex}→t₁={t1_ex}d")
    ax.set_xlabel("log2(fold change)"); ax.set_ylabel("Density")
    ax.set_title("Conditional density at representative (t₀, t₁)")
    ax.legend(fontsize=8); ax.set_xlim(xgrid[0], xgrid[-1])

    # Panel D: assignment vs starting EU
    ax = axes[1, 1]
    sc = ax.scatter(df["eu_start"], x, c=p_resp, cmap="RdYlBu_r",
                    s=25, alpha=0.7, edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)"); ax.set_xscale("log")
    ax.set_xlabel("Starting Vi IgG (EU)"); ax.set_ylabel("log2(fold change)")
    ax.set_title("Assignment vs starting antibody"); ax.axhline(0, color="gray", linestyle=":", lw=0.8)
    ax.grid(True, alpha=0.3)

    cv = sigma_log2_to_cv(p4["sigma1"])
    fig.suptitle(f"Step 4: Teunis Power-Law Mixture (n={len(x)}, tp={p4['t_peak']:.0f}d)\n"
                 f"CV≈{cv:.2f}, α₁={p4['alpha1']:.3f} (waning), "
                 f"α₂={p4['alpha2']:.3f}, μ_bio={p4['mu_bio']:.3f}",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step4_teunis_mixture.png", dpi=150, bbox_inches="tight"); plt.close()
    print("Saved 08_step4_teunis_mixture.png")


def plot_step5(df, p5):
    x = df["log2_fc"].values
    eu = df["eu_start"].values
    t0 = df["t_start"].values
    t1 = df["t_end"].values
    z = tau_ratio_covariate(t0, t1)
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    p_resp = 1 - p5["p_noise"]

    # Panel A: P(resp) vs starting EU, colored by fold change
    ax = axes[0, 0]
    sc = ax.scatter(eu, p_resp, c=x, cmap="coolwarm", s=25, alpha=0.7,
                    edgecolors="black", linewidths=0.3, vmin=-2, vmax=2)
    plt.colorbar(sc, ax=ax, label="log2(FC)")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Starting Vi IgG (EU)")
    ax.set_ylabel("P(responder)")
    ax.set_title("Fold-rise model: P(resp) vs starting titer\n(color = observed log2 FC)")
    ax.grid(True, alpha=0.3)

    # Panel B: predicted B₀ vs observed FC, colored by P(resp)
    ax = axes[0, 1]
    sc = ax.scatter(p5["b0_mean"], x, c=p_resp, cmap="RdYlBu_r", s=25, alpha=0.7,
                    edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    ax.plot([0, 10], [0, 10], "k:", linewidth=0.8, label="1:1")
    ax.set_xlabel("Predicted initial boost B₀ = log2(fold_rise(eu_start))")
    ax.set_ylabel("Observed log2(FC)")
    ax.set_title("Fold-rise prediction vs observed FC")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel C: fold-rise model curve + data
    ax = axes[1, 0]
    eu_grid = np.logspace(0, 4, 200)
    fr_grid = fold_rise_log2(eu_grid, p5["mu_0"], p5["log10_cop_max"])
    ax.plot(eu_grid, 2**fr_grid, "r-", linewidth=2, label=f"Fold-rise (μ₀={p5['mu_0']:.2f})")
    # Show observed FC colored by P(resp)
    sc = ax.scatter(eu, 2**x, c=p_resp, cmap="RdYlBu_r", s=25, alpha=0.6,
                    edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Starting Vi IgG (EU)")
    ax.set_ylabel("Fold change (natural scale)")
    ax.axhline(1, color="gray", linestyle=":", linewidth=0.8)
    ax.set_title(f"Fold-rise model fit\nμ₀={p5['mu_0']:.2f}, CoP_max={p5['cop_max']:.0f} EU")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel D: P(resp) vs Teunis covariate z
    ax = axes[1, 1]
    sc = ax.scatter(z, x, c=p_resp, cmap="RdYlBu_r", s=25, alpha=0.7,
                    edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="P(responder)")
    ax.set_xlabel("Teunis covariate z = log2(τ₁/τ₀)")
    ax.set_ylabel("log2(fold change)")
    ax.set_title("Component assignments vs time ratio")
    ax.grid(True, alpha=0.3)

    cv = sigma_log2_to_cv(p5["sigma1"])
    fig.suptitle(f"Step 5: Fold-Rise-Informed Mixture (n={len(x)})\n"
                 f"CV≈{cv:.2f}, μ₀={p5['mu_0']:.2f}, CoP_max={p5['cop_max']:.0f}, "
                 f"α₁={p5['alpha1']:.3f}, α₂={p5['alpha2']:.3f}",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_step5_foldrise_mixture.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 08_step5_foldrise_mixture.png")


def plot_summary(x, p1, p2, p3, p4, p5, df):
    xgrid = np.linspace(x.min() - 1, x.max() + 1, 500)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    bins = 30

    def hist_bg(ax):
        ax.hist(x, bins=bins, density=True, alpha=0.4, color="gray",
                edgecolor="black", linewidth=0.5)
        ax.set_xlabel("log2(FC)"); ax.set_xlim(xgrid[0], xgrid[-1])

    # Step 1
    ax = axes[0, 0]; hist_bg(ax)
    pdf_n = p1["pi_noise"] * stats.norm.pdf(xgrid, p1["mu1"], p1["sigma1"])
    pdf_s = (1 - p1["pi_noise"]) * stats.norm.pdf(xgrid, p1["mu2"], p1["sigma2"])
    ax.plot(xgrid, pdf_n, "b-", lw=2); ax.plot(xgrid, pdf_s, "r-", lw=2)
    ax.plot(xgrid, pdf_n + pdf_s, "k--", lw=1.5)
    ax.set_title(f"Step 1: Two Gaussians\nAIC={p1['aic']:.1f}, CV≈{sigma_log2_to_cv(p1['sigma1']):.2f}")

    # Step 2
    ax = axes[0, 1]; hist_bg(ax)
    pdf_n = p2["pi_noise"] * stats.norm.pdf(xgrid, p2["mu1"], p2["sigma1"])
    pdf_s = (1 - p2["pi_noise"]) * stats.skewnorm.pdf(xgrid, p2["alpha"], p2["xi"], p2["omega"])
    ax.plot(xgrid, pdf_n, "b-", lw=2); ax.plot(xgrid, pdf_s, "r-", lw=2)
    ax.plot(xgrid, pdf_n + pdf_s, "k--", lw=1.5)
    ax.set_title(f"Step 2: Gauss + Skew-Normal\nAIC={p2['aic']:.1f}, α={p2['alpha']:.2f}")

    # Step 3 at Δt=180d
    ax = axes[0, 2]; hist_bg(ax)
    mu2_eff = p3["mu2"] - p3["beta"] * np.log2(180)
    pdf_n = p3["pi_noise"] * stats.norm.pdf(xgrid, p3["mu1"], p3["sigma1"])
    pdf_s = (1 - p3["pi_noise"]) * stats.norm.pdf(xgrid, mu2_eff, p3["sigma2"])
    ax.plot(xgrid, pdf_n, "b-", lw=2); ax.plot(xgrid, pdf_s, "r-", lw=2)
    ax.plot(xgrid, pdf_n + pdf_s, "k--", lw=1.5)
    ax.set_title(f"Step 3: Signal time-dep\nAIC={p3['aic']:.1f}, β={p3['beta']:.4f}")

    # Step 4 at representative (t0=14, t1=180)
    ax = axes[1, 0]; hist_bg(ax)
    z_repr = tau_ratio_covariate(14, 180)
    mu_n = -p4["alpha1"] * z_repr
    pdf_n = p4["pi_noise"] * stats.norm.pdf(xgrid, mu_n, p4["sigma1"])
    pdf_s = (1 - p4["pi_noise"]) * truncnorm_conv_pdf(
        xgrid + p4["alpha2"] * z_repr, p4["mu_bio"], p4["sigma1"], p4["sigma_bio"])
    ax.plot(xgrid, pdf_n, "b-", lw=2); ax.plot(xgrid, pdf_s, "r-", lw=2)
    ax.plot(xgrid, pdf_n + pdf_s, "k--", lw=1.5)
    ax.set_title(f"Step 4: Teunis τ-ratio\nAIC={p4['aic']:.1f}, CV≈{sigma_log2_to_cv(p4['sigma1']):.2f}")

    # Step 5: MC-marginalized over empirical (eu_start, z) distribution
    ax = axes[1, 1]; hist_bg(ax)
    eu = df["eu_start"].values
    z_all = tau_ratio_covariate(df["t_start"].values, df["t_end"].values)
    # For each xgrid point, average the per-subject PDF over all subjects
    pdf_n5 = np.zeros_like(xgrid)
    pdf_s5 = np.zeros_like(xgrid)
    for i in range(len(df)):
        pdf_n5 += p5["pi_noise"] * stats.norm.pdf(xgrid, -p5["alpha1"] * z_all[i], p5["sigma1"])
        b0_mean_i = fold_rise_log2(eu[i], p5["mu_0"], p5["log10_cop_max"])
        pdf_s5 += (1 - p5["pi_noise"]) * truncnorm_conv_pdf(
            xgrid + p5["alpha2"] * z_all[i], b0_mean_i, p5["sigma1"], p5["sigma_bio"])
    pdf_n5 /= len(df)
    pdf_s5 /= len(df)
    ax.plot(xgrid, pdf_n5, "b-", lw=2); ax.plot(xgrid, pdf_s5, "r-", lw=2)
    ax.plot(xgrid, pdf_n5 + pdf_s5, "k--", lw=1.5)
    cv5 = sigma_log2_to_cv(p5["sigma1"])
    ax.set_title(f"Step 5: Fold-rise (marginalized)\nAIC={p5['aic']:.1f}, μ₀={p5['mu_0']:.2f}, "
                 f"CoP_max={p5['cop_max']:.0f}")

    # Step 5 detail: fold-rise curve
    ax = axes[1, 2]
    eu_grid = np.logspace(0, 4, 200)
    fr_grid = fold_rise_log2(eu_grid, p5["mu_0"], p5["log10_cop_max"])
    ax.plot(eu_grid, 2**fr_grid, "r-", linewidth=2, label=f"μ₀={p5['mu_0']:.2f}")
    p_resp5 = 1 - p5["p_noise"]
    sc = ax.scatter(eu, 2**x, c=p_resp5, cmap="RdYlBu_r", s=20, alpha=0.6,
                    edgecolors="black", linewidths=0.3, vmin=0, vmax=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.axhline(1, color="gray", linestyle=":", lw=0.8)
    ax.set_xlabel("Starting EU"); ax.set_ylabel("Fold change")
    ax.set_title(f"Step 5: Fold-rise model\nCoP_max={p5['cop_max']:.0f} EU")
    ax.legend(fontsize=8)

    fig.suptitle(f"Model comparison (n={len(x)}, all with σ₂ = √(σ₁² + σ_bio²))", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_model_comparison.png", dpi=150, bbox_inches="tight"); plt.close()
    print("Saved 08_model_comparison.png")


# =========================================================================
# Main
# =========================================================================

def main():
    df = load_fold_change_data()
    x = df["log2_fc"].values
    log2_dt = df["log2_duration"].values
    t0 = df["t_start"].values
    t1 = df["t_end"].values

    print(f"Loaded {len(df)} fold-change observations")
    print(f"  log2(FC): median={np.median(x):.3f}, SD={np.std(x):.3f}")
    print(f"  Duration: median={np.median(df['duration_days']):.0f} days, "
          f"range {df['duration_days'].min():.0f}-{df['duration_days'].max():.0f}")
    print(f"  t_start: median={np.median(t0):.0f}, range {t0.min():.0f}-{t0.max():.0f}")

    # Step 1
    print("\n" + "=" * 60)
    print("STEP 1: Two-Gaussian mixture")
    print("=" * 60)
    p1 = fit_two_gaussian(x)
    cv1 = sigma_log2_to_cv(p1["sigma1"])
    print(f"  Noise:  μ₁={p1['mu1']:.4f}, σ₁={p1['sigma1']:.4f}, π={p1['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv1:.3f}")
    print(f"  Signal: μ₂={p1['mu2']:.4f}, σ₂={p1['sigma2']:.4f}, σ_bio={p1['sigma_bio']:.4f}")
    print(f"  LL={p1['ll']:.2f}, AIC={p1['aic']:.2f}, BIC={p1['bic']:.2f}")
    plot_step1(x, p1)

    # Step 2
    print("\n" + "=" * 60)
    print("STEP 2: Gaussian + Skew-Normal")
    print("=" * 60)
    p2 = fit_gauss_skewnorm(x, p1)
    cv2 = sigma_log2_to_cv(p2["sigma1"])
    print(f"  Noise:  μ₁={p2['mu1']:.4f}, σ₁={p2['sigma1']:.4f}, π={p2['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv2:.3f}")
    print(f"  Signal: ξ={p2['xi']:.4f}, ω={p2['omega']:.4f}, α={p2['alpha']:.4f}, σ_bio={p2['sigma_bio']:.4f}")
    print(f"  LL={p2['ll']:.2f}, AIC={p2['aic']:.2f}, BIC={p2['bic']:.2f}")
    print(f"  ΔAIC vs Step 1: {p1['aic'] - p2['aic']:.2f}")
    plot_step2(x, p1, p2)

    # Step 3
    print("\n" + "=" * 60)
    print("STEP 3: Signal time-dependent (log2 Δt)")
    print("=" * 60)
    p3 = fit_gauss_timedep(x, log2_dt, p1)
    cv3 = sigma_log2_to_cv(p3["sigma1"])
    print(f"  Noise:  μ₁={p3['mu1']:.4f}, σ₁={p3['sigma1']:.4f}, π={p3['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv3:.3f}")
    print(f"  Signal: μ₂={p3['mu2']:.4f}, σ₂={p3['sigma2']:.4f}, β={p3['beta']:.4f}, σ_bio={p3['sigma_bio']:.4f}")
    print(f"  LL={p3['ll']:.2f}, AIC={p3['aic']:.2f}, BIC={p3['bic']:.2f}")
    print(f"  ΔAIC vs Step 1: {p1['aic'] - p3['aic']:.2f}")
    plot_step3(df, p3)

    # Step 4
    print("\n" + "=" * 60)
    print("STEP 4: Teunis power-law mixture (τ-ratio covariate)")
    print(f"  t_peak = {T_PEAK:.0f}d (fixed, from HlyE IgG in Aiemjoy Table S1)")
    print(f"  Noise:  N(-α₁·z, σ₁)        where z = log2(τ₁/τ₀)")
    print(f"  Signal: N(μ_bio - α₂·z, σ₂)  and τ = max(t - t_peak, 1)")
    print("=" * 60)
    p4 = fit_teunis_mixture(x, t0, t1, p1)
    cv4 = sigma_log2_to_cv(p4["sigma1"])
    print(f"  Noise:  α₁={p4['alpha1']:.4f}, σ₁={p4['sigma1']:.4f}, π={p4['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv4:.3f}")
    print(f"          α₁ is the Teunis power-law waning exponent")
    for t0_ex, t1_ex in [(14, 100), (14, 365), (30, 365)]:
        z_ex = tau_ratio_covariate(t0_ex, t1_ex)
        mu_n = -p4["alpha1"] * z_ex
        print(f"          Noise mean at t₀={t0_ex}→t₁={t1_ex}d (z={z_ex:.2f}): {mu_n:.4f} = {2**mu_n:.3f}× FC")
    print(f"  Signal: μ_bio={p4['mu_bio']:.4f}, α₂={p4['alpha2']:.4f}, σ₂={p4['sigma2']:.4f}, σ_bio={p4['sigma_bio']:.4f}")
    for t0_ex, t1_ex in [(14, 100), (14, 365), (30, 365)]:
        z_ex = tau_ratio_covariate(t0_ex, t1_ex)
        mu_s = p4["mu_bio"] - p4["alpha2"] * z_ex
        print(f"          Signal mean at t₀={t0_ex}→t₁={t1_ex}d: {mu_s:.4f} = {2**mu_s:.3f}× FC")
    print(f"  LL={p4['ll']:.2f}, AIC={p4['aic']:.2f}, BIC={p4['bic']:.2f}")
    print(f"  Converged: {p4['converged']}")
    print(f"  ΔAIC vs Step 1: {p1['aic'] - p4['aic']:.2f} ({'Step 4 better' if p1['aic'] > p4['aic'] else 'Step 1 better'})")
    plot_step4(df, p4)

    # Step 5
    eu_start = df["eu_start"].values
    print("\n" + "=" * 60)
    print("STEP 5: Fold-rise-informed mixture")
    print(f"  Signal boost = fold_rise(eu_start; μ₀, CoP_max)")
    print(f"  B₀ ~ TruncNorm(log2(fold_rise), σ_bio, ≥0), then waned + noise")
    print(f"  Noise unchanged: N(-α₁·z, σ₁)")
    print("=" * 60)
    p5 = fit_foldrise_mixture(x, t0, t1, eu_start, p1)
    cv5 = sigma_log2_to_cv(p5["sigma1"])
    print(f"  Noise:  α₁={p5['alpha1']:.4f}, σ₁={p5['sigma1']:.4f}, π={p5['pi_noise']:.3f}")
    print(f"          Implied assay CV = {cv5:.3f}")
    print(f"  Signal: μ₀={p5['mu_0']:.4f} (fold-rise intensity)")
    print(f"          CoP_max={p5['cop_max']:.1f} EU (log10={p5['log10_cop_max']:.2f})")
    print(f"          α₂={p5['alpha2']:.4f} (signal waning exponent)")
    print(f"          σ_bio={p5['sigma_bio']:.4f}, σ₂={p5['sigma2']:.4f}")
    # Show predicted boost at representative starting titers
    for eu_ex in [50, 200, 500, 1000]:
        b0 = fold_rise_log2(eu_ex, p5["mu_0"], p5["log10_cop_max"])
        print(f"          eu_start={eu_ex}: predicted log2(boost)={b0:.2f} = {2**b0:.1f}× fold rise")
    print(f"  LL={p5['ll']:.2f}, AIC={p5['aic']:.2f}, BIC={p5['bic']:.2f}")
    print(f"  Converged: {p5['converged']}")
    print(f"  ΔAIC vs Step 1: {p1['aic'] - p5['aic']:.2f} ({'Step 5 better' if p1['aic'] > p5['aic'] else 'Step 1 better'})")
    print(f"  ΔAIC vs Step 4: {p4['aic'] - p5['aic']:.2f}")
    plot_step5(df, p5)

    # Model comparison
    print("\n" + "=" * 60)
    print("MODEL COMPARISON (all with σ₂ = √(σ₁² + σ_bio²))")
    print("=" * 60)
    print(f"{'Model':<42} {'k':>3} {'LL':>8} {'AIC':>8} {'BIC':>8} {'CV':>6} {'σ_bio':>6}")
    print("-" * 85)
    for name, p, cv, k in [
        ("1. Two Gaussians", p1, cv1, 5),
        ("2. Gauss + Skew-Normal", p2, cv2, 6),
        ("3. Signal time-dep (log2 Δt)", p3, cv3, 6),
        ("4. Teunis power-law (τ-ratio)", p4, cv4, 6),
        ("5. Fold-rise-informed", p5, cv5, 7),
    ]:
        print(f"{name:<42} {k:>3} {p['ll']:>8.2f} {p['aic']:>8.2f} {p['bic']:>8.2f} "
              f"{cv:>6.3f} {p['sigma_bio']:>6.3f}")

    # Component assignments (best AIC)
    all_models = [(p1, "Step 1"), (p2, "Step 2"), (p3, "Step 3"),
                  (p4, "Step 4"), (p5, "Step 5")]
    best_p, best_name = min(all_models, key=lambda x: x[0]["aic"])
    print(f"\n{'='*60}")
    print(f"COMPONENT ASSIGNMENTS ({best_name}, P(responder) > 0.5)")
    print(f"{'='*60}")
    df["p_responder"] = 1 - best_p["p_noise"]
    responders = df[df["p_responder"] > 0.5]
    noise = df[df["p_responder"] <= 0.5]
    print(f"  Responders: {len(responders)}/{len(df)} ({100*len(responders)/len(df):.0f}%)")
    if len(responders) > 0:
        print(f"    Median FC: {2**responders['log2_fc'].median():.3f}, "
              f"duration: {responders['duration_days'].median():.0f}d, "
              f"starting EU: {responders['eu_start'].median():.1f}")
    print(f"  Noise: {len(noise)}/{len(df)} ({100*len(noise)/len(df):.0f}%)")
    if len(noise) > 0:
        print(f"    Median FC: {2**noise['log2_fc'].median():.3f}, "
              f"duration: {noise['duration_days'].median():.0f}d, "
              f"starting EU: {noise['eu_start'].median():.1f}")

    plot_summary(x, p1, p2, p3, p4, p5, df)


if __name__ == "__main__":
    main()
