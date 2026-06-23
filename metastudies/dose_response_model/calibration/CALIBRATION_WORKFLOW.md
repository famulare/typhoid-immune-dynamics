# Calibration workflow — iterative ladder

How the dose-response Stan calibration is meant to be built up, in deliberate
steps. Each step must sample cleanly and pass diagnostics before the next is
added. Terminology reconciles with `joint_inference_plan.md` (which defines
**Tier 1** = 25 obs / 11 params incl. study RE, and **Tier 2** = +6 Oxford
shedding / +η, 13 params).

| Step | = plan concept | obs | active params | status |
|---|---|---:|---|---|
| **1. Minimal Tier 1** | sub-milestone *below* repo Tier 1 | 25 | 10 (6 bio + 4 nuisance); `sigma_study` inert | **clean: 0/4000 divergences, R-hat<1.01, ESS>1700** |
| **2. Tier 1 complete** | repo **Tier 1** | 25 | 11 (+ study RE) | not started |
| **3. Tier 2** | repo **Tier 2** | 31 | 13 (+ η_lo, κ) | not started |

Run with [fit_dose_response.R](fit_dose_response.R) (cmdstanr + CmdStan 2.39).

## Workflow tooling (added 2026-06-23 — Buffalo-style, R-native)
- [priors.yaml](priors.yaml) — **single source** for priors; the `.stan` reads
  hyperparameters from its data block (change a value + refit, no recompile).
  [priors.R](priors.R) maps yaml → Stan data + generic prior samplers/densities.
- [data_prep.R](data_prep.R) — `build_stan_data()` builds the flat per-observation
  layout (one `obs_prob()` in the `.stan` handles all 5 likelihood groups).
- [diagnostics.R](diagnostics.R) — `diagnose_fit()` emits the standard battery
  (trace/pairs/rank/energy/density/rhat/neff via bayesplot + posterior),
  prior-vs-posterior overlay, PPC fed by Stan `p_pred` (no R likelihood mirror),
  priorsense power-scaling, `summary.md`/`summary.csv`/`results.json`/`fit.rds`.
- [simulate_recovery.R](simulate_recovery.R) — known-truth recovery harness
  (`generate_quantities` simulator → refit), SBC-lite coverage, and the controlled
  cliff-vs-reparam divergence attribution.
- [run_scenarios.R](run_scenarios.R) — `run_scenario()`/`summarize_scenarios()`
  cross-run comparison (row-filter + prior-override sensitivities, loo on grouped
  units). Deps: bayesplot, loo, priorsense, yaml, gridExtra (+ cmdstanr/posterior).

## Step 1 — Minimal Tier 1 (current)
Goal: a clean-sampling 25-observation fit (compiles, 0 divergences, R-hat < 1.01,
ESS > 400) with prior + posterior predictive checks. Oxford shedding/η and the
study RE are deliberately excluded; `delta`, the Maryland mixture, and `φ` are in.

**Divergence pathology resolved (2026-06-23).** Root cause was the ordering-constraint
cliff (`(log10_N50_fevginf - log10_N50_inf) ~ normal(0,1) T[0,]`), which produced
~99% divergences. Fixed by the prior-preserving reparam
`log10_N50_fevginf = log10_N50_inf + d_fev`, `d_fev<lower=0>` (the three original
N50 prior terms are retained; change-of-variables Jacobian = 1). The controlled
attribution in [simulate_recovery.R](simulate_recovery.R) confirms it:

| N50 parameterization (same data/priors/likelihood) | divergences |
|---|---|
| cliff `T[0,]` on the difference | 3959 / 4000 (99.0%) |
| reparam `d_fev <lower=0>` offset | **0 / 4000 (0.0%)** |

The unified `obs_prob()` refactor is verified against the original per-group
likelihood by [test_obs_prob_parity.R](test_obs_prob_parity.R) (the only guard the
recovery harness cannot provide — run it before trusting the harness).

**Remaining (science decisions, out of scope of the workflow upgrade):**
1. Constant-φ cap — the PPC now visibly shows high-dose Hornick fever (H-F-8/H-F-9)
   underfit because φ=0.25 caps fitted fever (reviewer-2 MC1). Needs dose-dependent φ.
2. δ prior-data tension — `log10_delta` posterior ≈1.6 vs prior mean 3.5; the
   `delta_prior_*` scenarios show residual prior pull. priorsense flags it as
   prior-data conflict.

Expected limitation even when clean: γ_inf is only weakly identified without the
Oxford vaccine shedding contrast (plan §2.6); `alpha_fevginf` shows the highest
priorsense prior-sensitivity (~0.56), consistent with weak fever-heterogeneity ID.

## Step 2 — Tier 1 complete (add study random effect)
Wire the design's study-level random effect (`sigma_study`, currently inert): add
a `study` index per observation to the data, use a non-centered
`z_study ~ normal(0,1)` scaled by `sigma_study`, applied on the logit of the
Maryland trial probabilities. This is repo-canonical Tier 1.

## Step 3 — Tier 2 (Oxford shedding + η)
Restore the 6 Oxford shedding rows (`tier2_active==1`, `N_ox_inf>0`). **Decide
first:** η Option A (parametric `eta_lo`, `κ` — what the Stan code currently
implements) vs Option C (fixed per-obs `eta_fixed_optC` from the CSV).
Recommendation: Option A. Tier 2 is preferred for γ_inf identification.

## Cross-cutting / deferred (track across steps)
- **CoP titer model:** non-naive Oxford `CoP` (1.15, 1.12, 5.0, 2.0) are
  placeholders; γ posteriors are conditional on them until the external
  anti-Vi → CoP map exists.
- **Hornick 10³ (H-F-3 = 0/14):** reviewer-2 MC4 flags 0/14 vs 9/14; run as a
  sensitivity once the base fit is clean.
- **`Gil-F-rest`** derived by subtraction — fit with/without as a sensitivity.
