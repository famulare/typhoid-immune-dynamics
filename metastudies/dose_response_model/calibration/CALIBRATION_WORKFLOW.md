# Calibration workflow — iterative ladder

How the dose-response Stan calibration is meant to be built up, in deliberate
steps. Each step must sample cleanly and pass diagnostics before the next is
added. Terminology reconciles with `joint_inference_plan.md` (which defines
**Tier 1** = 25 obs / 11 params incl. study RE, and **Tier 2** = +6 Oxford
shedding / +η, 13 params).

| Step | = plan concept | obs | active params | status |
|---|---|---:|---|---|
| **1. Tier 1 (floated φ)** | sub-milestone *below* repo Tier 1 | 25 | 11 (6 bio + 5 nuisance incl. `phi_md`); `sigma_study`/`eta_lo`/`κ` inert | **clean: 0/4000 div, R-hat≤1.002, ESS>1700; φ̂≈0.97 (uniform prior, edge-pressing), δ̂≈280×** |
| **1.5. EU/mL axis + individual Darton** | immunity upgrade ([issue #15](https://github.com/famulare/typhoid-immune-dynamics/issues/15), `tier1.5_plan.md`) | 25 grp + 30 indiv | CoP in VaccZyme EU/mL; Darton placebo individual; (+φ) φ(T) from threshold ladder | building minimal+φ |
| **2. Tier 1 complete** | repo **Tier 1** | 25 | 12 (+ study RE) | not started |
| **3. Tier 2** | repo **Tier 2** | 31 | 14 (+ η_lo, κ) | not started |

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
- [reference_points.csv](reference_points.csv) — **single source** for the
  hand-tuned cohort-model point estimates (from
  `scratch/cohort_incidence_model_proof_of_concept.R` and
  `..._high_dose_explore_scratch.R`), mapped to Stan parameters **by equation**
  (both use the identical beta-Poisson `P = 1 − (1 + D·(2^(1/α)−1)/N50)^(−α/CoP^γ)`).
  `diagnose_fit()` loads it by default (`load_reference_points()` in diagnostics.R)
  and draws them as labeled orange triangles on the x-axis of the matching
  `prior_posterior.png` facet — so they **always appear when the plot is
  regenerated**. They are reference markers only and do **not** inform the prior.
  Off-panel values are pinned at the facet edge with the value flagged "off-scale"
  so a far-tail point can't squash the densities. **Two read caveats baked into
  the table's `note` column:** (1) **dose-frame** — the script N50s are in raw
  ingested-bacilli (one curve blending Hornick milk + Oxford bicarb, no `delta`),
  whereas calibration `N50_inf`/`N50_fevginf` are bicarb-equivalent CFU with a
  separate `log10_delta` offset, so the N50 gap is frame-confounded, not pure
  miscalibration; (2) **fever object** — the script's `*_fever_given_dose` params
  describe the P(fever|dose) **marginal**, while `*_fevginf` describe P(fever|
  infection) **conditional** (mapped by functional-form role only). The infection
  triple (`N50_inf`, `alpha_inf`, `gamma_inf`) is an exact same-object mapping.
  To drop/add a marker, edit the CSV — no code change. (The fold-rise `mu_0`,
  2.5 vs 1.25, also differs between scripts but is an immunity/boosting parameter
  outside this dose-response calibration's parameter set, so it is not plotted.)
- [simulate_recovery.R](simulate_recovery.R) — known-truth recovery harness
  (`generate_quantities` simulator → refit), SBC-lite coverage, and the controlled
  cliff-vs-reparam divergence attribution. `recover_from_prior()` is the standard
  single-prior-draw check (draw one truth from the prior → simulate synthetic data
  matching the real covariate table → refit → recovery report); run it via
  [recover_from_prior_example.R](recover_from_prior_example.R), output to
  `results/recovery/tier1/`. Per-tier reuse = a `*_REPORT_PARS`/`*_INERT_PARS` pair
  + the `tier_col`; the simulator/diagnostics are unchanged (likelihood lives only
  in the `.stan`).
- [run_scenarios.R](run_scenarios.R) — `run_scenario()`/`summarize_scenarios()`
  cross-run comparison (row-filter + prior-override sensitivities, loo on grouped
  units). Deps: bayesplot, loo, priorsense, yaml, gridExtra (+ cmdstanr/posterior).

## CoP units (definitional)

The correlate of protection (CoP) is **denominated in anti-Vi IgG titre on the
commercial VaccZyme ELISA scale** (The Binding Site; EU/mL, LLD 7.4) — the modern
standard for typhoid anti-Vi. **Verified** that both modern Oxford inputs use this
exact assay: Jin 2017 (p.2474) and Darton 2016 (p.4), same LLD. CoP is expressed
**relative to the naive reference** so CoP = 1 at naive (anti-Vi < LLD, imputed
~3.7 EU/mL). It enters the dose-response as `CoP^gamma` — a per-log₁₀-titre power
law — so **γ is the protection slope, anchorable to Darton's HR = 0.29 per log₁₀
anti-Vi** (same scale, no cross-assay calibration). Maryland is latent
anti-Vi-*equivalent* EU/mL (anti-Vi unmeasured pre-VaccZyme).

Reference titres (source of truth): Jin Vi-TT 562.9, Vi-PS 140.5, control 8.0
EU/mL (measured GMTs); Darton placebo per-subject in S1; naive 3.7 (LLD-imputed).

**Status:** the *units definition* is wired into the contracts and the `.stan`.
The per-group **value-swap + refit** (replacing the interim CoP placeholders
1.15/1.12/5.0/2.0 with titre-relative values, rescaling the Maryland latent CoP
priors, and re-anchoring γ to the Darton HR) is the immediate next step, pending
three choices logged in `tier1_lab_notebook.md` D1 (naive reference; Maryland
latent prior scale; Darton per-subject vs GMT).

## Step 1 — Tier 1 with floated φ and δ (current)

Goal: a clean-sampling 25-observation fit (compiles, 0 divergences, R-hat < 1.01,
ESS > 400) with prior + posterior predictive checks. Oxford shedding/η and the
study RE are deliberately excluded. **Both Maryland nuisance scales now float:**
`delta` (always was estimated) and — as of 2026-06-23 — `phi_md`, a single
**estimated scalar** Maryland fever definition-sensitivity (`Beta(1,1)` uniform
prior; a deliberate loosening from plan §7's `Beta(5,5)`, see below), replacing the
fixed per-obs `φ` data (0.25/0.65). The Maryland mixture is also in.

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

**φ resolution (2026-06-23).** The constant-φ cap is removed: `φ=0.25` was the
*low-dose asymptote* of the definition-sensitivity function applied dose-wide, so
it capped Maryland fitted fever at 0.25 while Hornick's high-dose data are 0.89–0.95
(structurally unfittable, not just biased; reviewer-2 MC1). Rather than dose-dependent
φ (+1–2 params on already-thin ID), Tier 1 now **estimates a single scalar `phi_md`**
identified by the Hornick dose-range plateau (at saturation `P_fev→1` so `fitted→φ`,
which directly measures φ near the plateau). Plan §7 specified `φ ~ Beta(5,5)`, but
that is too strong here, so Tier 1 uses **`Beta(1,1)` (uniform)** — a deliberate
departure (see two-stage history below). Decision log:

- **Single global `phi_md`, not per-study.** Gilman/Levine are single-dose (all 10⁵)
  so they cannot identify their own φ (confounded with the 10⁵ attack rate); only
  Hornick's dose range can. The old 0.25 (Hornick) vs 0.65 (Gilman/Levine) split was
  itself unaudited LLM-derived; cross-study level differences move to the study RE
  (Step 2) / residual. Per-study or ratio-preserving φ is a documented alternative.
- **Prior choice.** `Beta(5,5)` (90% mass [0.25,0.75], ~10 pseudo-obs at 0.5) pulled
  φ̂ to 0.886 — the prior's 99.8th pctile — and still left H-F-9 underfit, so it was
  loosened to `Beta(1,1)`. Fallback if φ destabilizes at the [0,1] edge: a
  mildly-informed elicited prior (e.g. `Beta(2,2)`).

**Re-fit result (2026-06-23, `Beta(1,1)`).** 0/4000 divergences, R-hat ≤ 1.002,
ESS > 1700 (φ_md ess_tail ~1250). See `results/tier1/dose_response_fit.png` for the
bespoke dose-response PPC.

- **`phi_md` = 0.966** (median 0.975, 90% CI 0.91–0.998) — floats up to the plateau
  and **presses the [0,1] boundary**: the data wants ≈no definitional suppression, so
  the whole Maryland↔Oxford gap is carried by δ, not φ.
- **`log10_delta` rose to 2.46** (δ ≈ 280×; was 1.6 / 42× under fixed φ) — continuing
  toward the prior mean 3.5. The φ-cap and δ tension were the same artifact, confirmed.
- **`alpha_fevginf` re-identified** (priorsense prior-sensitivity 0.56 → ~0.06); the
  fixed cap had been corrupting it.
- **Key residual — the binding constraint SHIFTED off φ.** High-dose Hornick is *still*
  underfit ~0.12 (H-F-8 0.889 → fit 0.77; H-F-9 0.952 → fit 0.83) **even with φ≈0.97**.
  Cause is no longer φ: the Maryland mixture's immune component (~39% at CoP_imm≈2.2)
  plus the slow `P_fev|inf` saturation (small `alpha_fevginf`) cap the *predicted*
  plateau at ~0.86. A tighter φ prior would lower this further — so the elicited φ
  prior addresses the boundary/ID, but the high-dose misfit is a separate structural
  question (is immunity fully overwhelmed at 10⁹? does fever|infection saturate faster?).
- δ↔N50 ridge now visible (r ≈ −0.76/−0.79) — dose-scale confound, a Tier-2 concern.
- Gilman/Levine 10⁵ fever scatter (±0.2 residual) and the H-FgI-7 conditional overfit
  (0.57 obs vs 0.75 fit) are study-RE / fever-process items for Steps 2+.

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
