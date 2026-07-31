# Tier 1.5 plan — VaccZyme EU/mL immunity axis + individual-level Darton

**Status:** C0–C3, `+cascade` and `+Jin`-as-CoP-anchor are DONE; `+vaccine-terms` parked. The resulting configuration is `t1-indiv` (see `TIER_LADDER.md`). Last updated 2026-07-31.
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
- **+vaccine-terms — PLANNED** (design below).

## +vaccine-terms increment (planned — design + Darton VE, 2026-07-31)

Add Darton's two non-anti-Vi vaccine arms (currently extracted but unused). Motivation and
the exact published efficacy, so the design rests on numbers not memory.

### The unused data (Darton 2016, individual endpoints CSV)
| arm | n | infection (bact_or_stool) | fever (TD) | anti-Vi median (max) | mechanism |
|---|---|---|---|---|---|
| Placebo (in fit) | 30 | 0.87 | 0.67 | 3.7 (62) | — |
| **M01ZH09** | 31 | 0.68 | 0.58 | 3.7 (204) | live oral, Ty2 ΔaroC ΔssaV |
| **Ty21a** | 30 | 0.53 | 0.43 | **NA (not assayed)** | live oral, **Vi-negative** |

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

### Proposed parameterization
- Add M01ZH09 (31) + Ty21a (30) as individual cascade rows, **same groups as placebo**
  (6 = infection `bact_or_stool`; 7 = fever|inf `fever_td` among infected). +61 subjects.
- **Per-vaccine protection factor `V_v`, a separate channel from CoP^γ:**
  `exponent = −alpha / (CoP^gamma · V_v)`, with `V_placebo = 1`, `V_v ≥ 1` protective,
  `log V_v ~ Normal(·)` weakly-informative. Mirrors the CoP^γ form (protection divides the
  exponent → lowers the attack rate) but as an **independent** per-vaccine magnitude, so it
  **stays out of γ** (the whole point). Ty21a: `CoP = 1` (anti-Vi NA) → `V_Ty21a` carries all
  its protection; M01ZH09: `CoP = measured × V_M01ZH09`.
- **Params:** +2 (`V_M01ZH09`, `V_Ty21a`); optionally per-endpoint (+2) if inf vs fever VE
  diverge (Ty21a bacteraemia 41% vs fever≥38 48% — thin at n≈30; start single, split only if
  a residual demands it).
- **Priors:** weakly-informative — each arm's own attack-rate contrast vs placebo identifies
  `V_v` cleanly (61 subjects, no cross-arm confound). Optionally anchor on the Darton VE the
  way +Jin anchored γ.
- **Harness:** groups 6/7 are already in the parity gate + `model_math.R`. Adding `V_v` needs
  a `vaccine_id` covariate + the V lookup in `obs_prob`, plus `PARAM_NAMES`/`vecs`/`TRUTH`
  and `priors.yaml`; the `assert_fitted_params_match` guard flags the stale lists.

### Wrinkles
- **Ty21a anti-Vi = NA** → it *must* be a pure vaccine effect (CoP=1); asymmetric with
  M01ZH09 (which has a measured, mostly-naive titre). This is a data fact, not a choice.
- Treatment-truncation (η) ignored at Tier 1, same as placebo.
- Distinct from the model's placeholder `VE_fev_ViTT/ViPS` GQ (those are anti-Vi *Vi*-vaccines;
  this is the orthogonal non-anti-Vi case).

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
