# Calibration workflow — iterative ladder

How the dose-response Stan calibration is meant to be built up, in deliberate
steps. Each step must sample cleanly and pass diagnostics before the next is
added. Terminology reconciles with `joint_inference_plan.md` (which defines
**Tier 1** = 25 obs / 11 params incl. study RE, and **Tier 2** = +6 Oxford
shedding / +η, 13 params).

| Step | = plan concept | obs | active params | status |
|---|---|---:|---|---|
| **1. Minimal Tier 1** | sub-milestone *below* repo Tier 1 | 25 | 10 (6 bio + 4 nuisance); `sigma_study` inert | runs; **blocked** by divergence pathology |
| **2. Tier 1 complete** | repo **Tier 1** | 25 | 11 (+ study RE) | not started |
| **3. Tier 2** | repo **Tier 2** | 31 | 13 (+ η_lo, κ) | not started |

Run with [fit_dose_response.R](fit_dose_response.R) (cmdstanr + CmdStan 2.39).

## Step 1 — Minimal Tier 1 (current)
Goal: a clean-sampling 25-observation fit (compiles, 0 divergences, R-hat < 1.01,
ESS > 400) with prior + posterior predictive checks. Oxford shedding/η and the
study RE are deliberately excluded; `delta`, the Maryland mixture, and `φ` are in.

**Blocked.** The model compiles and samples but produces ~99% divergent
transitions. Root cause isolated to the ordering-constraint cliff at
`typhoid_dose_response.stan` L132; see
[tier1_pathology_diagnosis.md](tier1_pathology_diagnosis.md). Two things must be
resolved before this step yields a trustworthy posterior:
1. Reparameterize the N50 ordering constraint (cliff → `<lower=0>` offset).
2. Resolve the constant-φ cap that makes high-dose Hornick fever unfittable.

Expected limitation even when clean: γ_inf is only weakly identified without the
Oxford vaccine shedding contrast (plan §2.6, "Oxford alone is underconstrained").

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
