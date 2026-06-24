# Tier 1.5 plan — VaccZyme EU/mL immunity axis + individual-level Darton

**Status:** in progress (2026-06-23). Tracked as GitHub issue (see lab notebook D2).
Builds on Tier 1 (`CALIBRATION_WORKFLOW.md`, `tier1_lab_notebook.md`).

## Why

Tier 1 samples cleanly, but **immunity is the soft spot**: γ (immunity scaling) is
identified almost entirely through two *placeholder* Jin CoPs (5.0/2.0); φ is a free
scalar pressing the [0,1] boundary; and the richest evidence in the repo —
**Darton's 91-subject individual-level anti-Vi titres + 5-threshold fever data, all
on the VaccZyme EU/mL scale** (`analysis_data/darton_individual_endpoints.csv`) — is
extracted but consumed by nothing in the fit. Only one Darton row (placebo fever
20/30, placeholder CoP=1.15) is in the Tier-1 likelihood.

Verified this session: both Jin 2017 (p.2474) and Darton 2016 (p.4) use the
commercial **VaccZyme anti-Vi IgG ELISA (EU/mL, LLD 7.4)** — one comparable scale.
Decision: adopt **VaccZyme EU/mL as the definitional CoP unit** and pull Darton's
individual data into the likelihood. Goal of the first increment: **see whether
identifiability improves** (γ_fev, φ, the δ↔N50 ridge).

## Full Tier 1.5 vision (the upgrade)

Individual-level Darton sub-likelihood over all 91 per-protocol subjects:
- per-subject anti-Vi titre (EU/mL) → fever/infection outcomes (Bernoulli);
- 5-threshold fever (37–39 °C) → identify φ(T) **as data** (the Oxford ladder φ was
  originally eyeballed from);
- bacteremia + stool per subject → separate **γ_inf vs γ_fevginf** (the cascade);
- vaccine-group terms for the non-Vi arms (M01ZH09, Ty21a protect via non-anti-Vi
  mechanisms → would confound the anti-Vi axis if pooled naively);
- stretch: digitize Jin Fig S3 for individual high-titre points.

**Increment ladder:** minimal → +φ → +cascade → +vaccine-terms → +Jin-digitize.

## Increment built first: minimal + φ

### C0. VaccZyme EU/mL immunity axis (value-swap; foundation)
- CoP = anti-Vi IgG (EU/mL) / naive-ref, **naive-ref = 3.7** (the extract's <LLD
  imputation). CoP=1 at naive. Enters unchanged as `CoP^gamma` (per-log₁₀-titre
  power law); the reference is fit-invariant under this form (absorbed into N50/γ).
- Oxford group titres → real GMTs: Jin Vi-TT **563** (CoP 152), Vi-PS **141** (38),
  control **8.0** (2.2); Waddington/Gibani naive **3.7** (1.0).
- γ re-scales (CoP range 1–152, not 1–5) → γ posterior much smaller; expected.
  γ_fev prior stays **vague** — the Darton individual data identifies it.
- Maryland latent `CoP_susc`/`CoP_imm` now anti-Vi-equivalent EU/mL-relative →
  priors rescaled: `CoP_imm ~ lognormal(log 5, 0.7)`, `CoP_susc ~ lognormal(0, 0.3)`
  (weakly-informative; `pi_susc` unchanged).

### C1. Darton placebo individual titre→fever (minimal)
- Replace the D-F-plac group binomial with the **30 placebo subjects as n=1 rows**
  (`CoP_i = vi_igg_prechallenge_i / 3.7`, `y = fever_td_i`). Fits the existing flat
  `obs_prob()` unchanged. Built programmatically from the S1 extract in `data_prep.R`
  (single source of truth; not hand-transcribed). Observed gradient: undetectable
  72% fever vs detectable 58% (n=30, low-titre).

### C2. Jin (group, EU/mL)
- Jin arms stay group binomials at real GMTs (563/141/8) — the high-titre anchor.
  Vi-PS Jensen-bias caveat noted (Fig S3 digitization deferred).

→ **First fit = C0+C1+C2**, then read identifiability before C3.

### C3. Darton multi-threshold → φ(T) (the φ piece)
- Add the placebo subjects' 5 threshold-exceedance indicators.
- Model φ(T) = P(peak temp ≥ T | fever at Oxford reference) as a monotone-decreasing
  1-parameter curve (exponential decay; open to logistic) anchored at the reference,
  fit to the Darton ladder. Maryland reads φ at each study threshold (Hornick 39.4,
  Levine/Gilman 38.3) — **replacing the free `phi_md` scalar.**
- Caveat: ladder is at one (Oxford) dose; 39.4 is a modest extrapolation beyond the
  39.0 max, now with propagated uncertainty. (Dose-dependent φ(T,D) = future.)

## Files
- `dose_response_data.csv` — Oxford group CoP → EU/mL-relative; D-F-plac group row
  removed (replaced by individual rows injected in `data_prep.R`).
- `data_prep.R` — inject Darton placebo individual rows (`NAIVE_VI_REF = 3.7`); (C3)
  per-subject thresholds.
- `priors.yaml` — Maryland CoP priors rescaled; γ_fev vague; (C3) φ(T) hyperparams.
- `typhoid_dose_response.stan` — CoP already documented as EU/mL/naive; (C3) φ(T)
  sub-model + threshold likelihood, retire `phi_md`.
- `dose_response_curves.R` / `diagnostics.R` — figure extended for individual Darton
  points + Jin EU/mL + (C3) the φ(T) curve.

## Verification — SEE identifiability
- Compile + 4-chain sample; target 0 divergences, R-hat ≤ 1.01.
- vs the current Beta(1,1) Tier-1 fit: does **γ_fevginf** move prior-dominated →
  data-identified (priorsense)? does **φ** stop edge-pressing (C3)? does **δ↔N50**
  loosen? `γ_inf` likely still wide pre-cascade (expected).
- Bespoke dose-response figure with 30 individual Darton points + Jin EU/mL spanning
  the titre axis; post-fit γ_fev sanity vs the Darton per-log₁₀ slope (validation).
- Preserve the prior Tier-1 fit for comparison (don't overwrite `results/tier1`).
