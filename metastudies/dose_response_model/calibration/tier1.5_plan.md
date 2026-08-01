# Tier 1.5 plan — VaccZyme EU/mL immunity axis + individual-level Darton

**Status:** C0–C3, `+cascade`, `+Jin`-as-CoP-anchor, and `+vaccine-terms` are ALL DONE. The
resulting configurations are `t1-indiv` (no vaccine arms) and `t1-indiv-vax` (+ M01ZH09/Ty21a)
(see `TIER_LADDER.md`). Last updated 2026-07-31.
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
- Bespoke dose-response figure with the individual Darton points (the +cascade increment
  made this 30 infection + 26 fever|infection = 56 rows) + Jin EU/mL spanning
  the titre axis; post-fit γ_fev sanity vs the Darton per-log₁₀ slope (validation).
- Preserve prior fits for comparison. The driver writes `results/<spec>__<stage>/`, so a
  new model stage never overwrites the previous one; keep a deliberate snapshot as
  `results/<spec>__<stage>__<tag>/`. Pre-2026-07-31 run dirs were deleted as non-reproducible.

---

## Later increments — status + designs

- **+cascade — DONE** (commit `0088b9e`): Darton placebo per-subject infection
  (`bact_or_stool`, group 6 `ox_inf_indiv`) + fever|infection (`fever_td` among infected,
  group 7 `ox_fevginf_indiv`). Result: the γ_inf/γ_fevginf split is titre-range-limited
  (Darton clusters at low anti-Vi), slopes near-identical (~0.15–0.20).
- **+Jin — resolved as CoP-anchor** (commit `9c5126e`): Jin's published logistic OR
  (0.37/log₁₀ anti-Vi, verified) informs the γ prior in lieu of digitizing Fig S3.
  Higher-fidelity +Jin-digitize (WebPlotDigitizer/human or deposited trial data) stays open.
- **+vaccine-terms — DONE** (2026-07-31; see results after the design below).

## +vaccine-terms increment (planned — design + Darton VE, 2026-07-31)

Add Darton's two non-anti-Vi vaccine arms (currently extracted but unused). Motivation and
the exact published efficacy, so the design rests on numbers not memory.

### The unused data (Darton 2016, individual endpoints CSV)
| arm | n | infection (bact_or_stool) | fever (TD) | anti-Vi median (max) | mechanism |
|---|---|---|---|---|---|
| Placebo (in fit) | 30 | 0.87 | 0.67 | 3.7 (62) | — |
| **M01ZH09** | 31 | 0.68 | 0.58 | 3.7 (204) | live oral, Ty2 ΔaroC ΔssaV |
| **Ty21a** | 30 | 0.53 | 0.43 | 3.7 (155), 29/30 measured | live oral, **Vi-negative** |

**CORRECTION 2026-07-31** (caught by Mike): the table above originally said Ty21a anti-Vi
is "NA (not assayed)" — wrong. That claim came from the STALE grouped CSV row's note
("Cannot map to anti-Vi CoP"), which describes the *vaccine's mechanism* (Ty21a doesn't
raise anti-Vi), over-read as "no titre data exists." The individual endpoints file has a
real, measured `vi_igg_prechallenge` for 29/30 Ty21a subjects (mostly at the naive floor
3.7, but 8/29 detectable/elevated up to 155 — presumably pre-existing unrelated exposure,
same as Placebo). Only 1 subject is genuinely missing (dropped, not imputed). This changes
the design below: Ty21a is individualized the SAME way as Placebo (own titre -> CoP^gamma),
NOT forced to CoP=1.

### Darton published vaccine efficacy (Table 2, p.10) [from extract]
| endpoint | VE M01ZH09 [95% CI] | VE Ty21a [95% CI] |
|---|---|---|
| Primary TD (unadjusted) | 13% [−29, 41] | 35% [−5, 60] |
| Primary TD (**adj. for baseline anti-Vi**) | 19% [−17, 43] | 31% [−8, 55] |
| Fever ≥38.0 °C (adj.) | 19% [−27, 48] | **48% [4, 72]** |
| Any bacteraemia (adj.) | 28% [−7, 52] | **41% [2, 64]** |
| Bacteraemia or stool positive (unadj.) | 22% [−3, 41] | **38% [12, 57]** |

**Key findings** [observed]: (1) **M01ZH09 is weak / non-significant** on every endpoint
(all CIs cross 0) — matches Darton's headline. (2) **Ty21a is moderate and SIGNIFICANT** on
the infection/bacteraemia endpoints (bact-or-stool 38% [12,57]; adj. bacteraemia 41% [2,64];
adj. fever≥38 48% [4,72]). (3) **Adjusting for baseline anti-Vi barely changes the VE**
(Ty21a TD 35%→31%, M01ZH09 13%→19%) → the protection is **demonstrably NOT anti-Vi-mediated**.
This is the clean justification for a *separate* protection channel.

### The confounding problem (why they can't just be added)
Ty21a is protected (infection 0.53 vs placebo 0.87) while sitting at **naive anti-Vi**. On the
CoP = anti-Vi/3.7 axis they land at CoP≈1 with low attack rates; pooled into the cascade the
model would flatten **γ** (mis-attribute cell-mediated protection to the anti-Vi slope) or
inflate the baseline. So the arms need a per-vaccine term that absorbs their protection.

### Proposed parameterization — SUPERSEDED, see IMPLEMENTED below
*(original text, kept for the record; the "Ty21a: CoP=1" line was wrong per the correction
above and was NOT what got built)*
- Add M01ZH09 (31) + Ty21a (30) as individual cascade rows, **same groups as placebo**
  (6 = infection `bact_or_stool`; 7 = fever|inf `fever_td` among infected). +61 subjects.
- **Per-vaccine protection factor `V_v`, a separate channel from CoP^γ:**
  `exponent = −alpha / (CoP^gamma · V_v)`, with `V_placebo = 1`, `V_v ≥ 1` protective,
  `log V_v ~ Normal(·)` weakly-informative. Mirrors the CoP^γ form (protection divides the
  exponent → lowers the attack rate) but as an **independent** per-vaccine magnitude, so it
  **stays out of γ** (the whole point). ~~Ty21a: `CoP = 1` (anti-Vi NA) → `V_Ty21a` carries
  all its protection~~ WRONG, see correction above -- Ty21a uses its own measured titre too.
- **Params:** +2 (`V_M01ZH09`, `V_Ty21a`); optionally per-endpoint (+2) if inf vs fever VE
  diverge (Ty21a bacteraemia 41% vs fever≥38 48% — thin at n≈30; start single, split only if
  a residual demands it).
- **Priors:** weakly-informative — each arm's own attack-rate contrast vs placebo identifies
  `V_v` cleanly (61 subjects, no cross-arm confound). Optionally anchor on the Darton VE the
  way +Jin anchored γ.
- **Harness:** groups 6/7 are already in the parity gate + `model_math.R`. Adding `V_v` needs
  a `vaccine_id` covariate + the V lookup in `obs_prob`, plus `PARAM_NAMES`/`vecs`/`TRUTH`
  and `priors.yaml`; the `assert_fitted_params_match` guard flags the stale lists.

### Wrinkles (as originally written; see correction above)
- ~~Ty21a anti-Vi = NA~~ → WRONG. Both arms use their own titre; only 1 Ty21a subject
  (genuinely missing) is dropped.
- Treatment-truncation (η) ignored at Tier 1, same as placebo.
- Distinct from the model's placeholder `VE_fev_ViTT/ViPS` GQ (those are anti-Vi *Vi*-vaccines;
  this is the orthogonal non-anti-Vi case).

## +vaccine-terms IMPLEMENTED (2026-07-31) — new tier `t1-indiv-vax`

**What actually got built**, correcting the design above: M01ZH09 (31 subj, 21 infected)
and Ty21a (29 subj after dropping 1 missing titre, 16 infected) individualize the SAME way
as Placebo — own anti-Vi titre feeds the shared `CoP^gamma` channel via cascade groups 6/7
(`darton_arm_individual_rows()` in `data_prep.R`, generalized from the Placebo-only
`darton_placebo_individual_rows()`). A new `vaccine_id` covariate (0/1/2) routes groups 6/7
through `beta_poisson_vax()`, a new Stan/model_math kernel with an extra multiplicative
factor: `exponent = -alpha/(CoP^gamma * V)`, `V=1` recovering the plain kernel exactly.
`V_v = exp(log_V_v)`, `log_V_v ~ Normal(0,1)` (V=1/log_V=0 = "no additional effect"; the
data can move it either direction). New tier `t1-indiv-vax` (stage `phi-rho-vax`, N_obs=177
= 80 + 52 M01ZH09 + 45 Ty21a). Full parity gate green (5 tiers + 1164-row synthetic grid).

One regression caught mid-build: the generalization initially derived Placebo's obs_id tag
from the arm NAME ("Placebo"->"placebo"), which would have silently renamed the pre-existing
`D-I-plac-*` rows and broken `curve_specs.R`'s regex + `cohort_id` "OX-DAR-2013-PLAC". Fixed
by special-casing Placebo to keep its historical "plac" tag; M01ZH09/Ty21a get their
naturally-derived tags ("m01zh09"/"ty21a", NOT "m01" as first guessed -- `gsub` strips
non-alphanumerics, of which "M01ZH09" has none). Also fixed: `.ppc_rows()` in
`diagnostics.R` grouped Darton individual rows by `(study, likelihood_group, group,
dose_cfu)` only, which would have POOLED all three arms' PPC points into one misleading
marker (same trial, same 18200 CFU dose) -- added `cohort_id` to the grouping (a no-op for
every pre-existing single-cohort case).

### Fit results (t1-indiv-vax, 4x1000, adapt_delta 0.9)
Clean: **0/4000 divergences**, R-hat <= 1.003 on all 22 params, min E-BFMI 0.94.
Curve-parity check (figure math vs Stan p_pred) 4.4e-8, both prior and posterior.

**V_M01ZH09 = 1.40** (median, 90% CI [1.00, 1.87]) — CI lower bound sits almost exactly at
1 (no additional effect). **V_Ty21a = 1.84** (90% CI [1.30, 2.58]) — clearly excludes 1.
**This reproduces Darton's own published asymmetry** (M01ZH09 weak/CIs-cross-null on every
endpoint; Ty21a significant on infection/bacteraemia) [Table 2 above] -- but now as a joint
posterior estimate from the shared dose-response model, not Darton's separate regression.

**priorsense: both vaccine params strongly DATA-DRIVEN.** `log_V_M01ZH09` prior-sensitivity
0.011 vs likelihood 0.099 (~9x); `log_V_Ty21a` 0.018 vs 0.085 (~5x) -- no prior-data-conflict
or prior-dominance flag on either. Each arm's own contrast vs Placebo (same challenge,
individualized identically) identifies `V_v` cleanly, as hoped.

**Honest finding -- the shared-biology shift.** Adding 92 more individual cascade subjects
(M01ZH09+Ty21a) moved `alpha_inf` (0.278->0.203), `alpha_fevginf` (0.250->0.298), and
`gamma_inf` (0.172->0.118) by roughly **0.75-0.92 OLD-posterior-SD** -- a real, non-trivial
shift, and precisely the risk flagged in "The confounding problem" above (if `V_v` doesn't
fully absorb each arm's non-anti-Vi protection, the residual can leak into the shared
biology). Checked directly rather than waved away: **fitted medians on 10 key non-Darton
rows** (Hornick H-F-5/8/9, Gilman Hlo/Hhi, Levine 1/2, Jin ctrl/ViTT/ViPS) **moved by at
most ~0.03** between the pre- and post-vaccine-terms fits -- no degradation anywhere. Read:
`alpha_inf`/`alpha_fevginf`/`gamma_inf` sit on a weakly-identified ridge (posterior SDs
shrank ~30-40% with the added data), and this looks like more data resolving a previously
uncertain point on that ridge rather than the vaccine effect corrupting the shared curve --
but this is a judgment call on ambiguous evidence, not a clean gate pass, and is flagged as
such for review.

**PPC per-arm separation verified**: Ty21a's pooled PPC points (q1 naive-titre stratum:
11/21 infected obs 0.524 vs fit 0.636; 8/11 fever|inf obs 0.727 vs fit 0.689) land in the
right neighborhood, sampling-noise-consistent at n=11-21, and are no longer conflated with
Placebo/M01ZH09 (the `.ppc_rows()` fix).

**Not done / open:** per-endpoint split of `V_v` (infection vs fever|inf) if a residual
demands it; anchoring `log_V_v` priors on the Darton VE numbers (currently uninformed,
Normal(0,1)); the anti-Vi-vs-cell-mediated decomposition as its own analysis question
(V_Ty21a IS that decomposition's headline number now, not just a design proposal).

### Value / priority [honest]
- **Does NOT extend the anti-Vi axis** (both arms ~naive) → **zero help to the γ-split or the
  immunity slope** — the reason it's ranked below +Jin.
- **Does buy:** (a) +61 backbone subjects at the Oxford dose; (b) non-anti-Vi VE as a model
  quantity to validate against Darton Table 2; (c) enables the **anti-Vi-vs-cell-mediated
  protection decomposition** — a *distinct scientific aim* (Ty21a's significant,
  anti-Vi-independent protection is a clean natural experiment), arguably its own analysis
  rather than part of the dose-response/CoP calibration.
- **Verdict:** park unless the decomposition question is the goal. If pursued it's low-risk
  (`V_v` well-identified by its own arm) and self-contained — no threat to the existing fit.

## Beyond Tier 1.5: see tier2_plan.md

This document stops at `t1-indiv-vax` (stage `phi-rho-vax`). The next increment — Oxford
`ox_inf` shedding restored (η, already implemented at C3/pre-Tier-1.5) + ψ
infection-definition correction (`joint_inference_plan.md` Sec 2.8, adopted 2026-07-31,
implemented 2026-07-31) — is `tier2_plan.md`, LOCKED 2026-07-31. It produces `t2-indiv`
and `t2-indiv-vax` (stages `phi-rho-eta-psi` / `phi-rho-eta-psi-vax`); no new grouped Tier 2
configuration is being built (Tier 2 is individualized-Darton only going forward).
