# Tier 1 lab notebook — dose-response calibration

Running scientific log of the **Tier 1** iterations. Picks up *after* the
mechanical bug fixes + Buffalo-style workflow refactor (those are in
`tier1_pathology_diagnosis.md` / `CALIBRATION_WORKFLOW.md` and not re-litigated
here). Starting point: the model compiles and samples cleanly (cliff fixed via the
`d_fev<lower=0>` reparam, 0/4000 divergences), the likelihood is refactored
(`obs_prob()`, flat data, priors-as-data), **φ is still fixed at 0.25/0.65**, and
δ shows prior-data tension. Everything below is one continuous session, 2026-06-23.

Provenance tags: [observed] measured/in-fit, [inferred], [Mike], [from paper].

---

## E0 — Starting state (post-refactor, fixed φ)
- Sampling clean: 0/4000 div, R-hat ≤ 1.002 (after the N50 cliff fix). [observed]
- φ fixed as per-obs data: 0.25 (Hornick), 0.65 (Gilman/Levine). [observed]
- δ (`log10_delta`) ≈ 1.6 (δ ≈ 42×) vs prior `normal(3.5,0.7)` — prior-data
  conflict flagged by priorsense. [observed]
- `alpha_fevginf` prior-dominated (priorsense prior-sens 0.56). [observed]
- **Open live bug:** constant φ caps Maryland fitted fever at φ; PPC shows
  high-dose Hornick (H-F-8/9) underfit.

## E1 — Why φ=0.25 is there, and why it's wrong (provenance)
- φ=0.25 derivation [from `joint_inference_plan.md` §2.5]: extrapolated temperature
  sensitivity at Hornick's ≥103°F/39.4°C threshold (~0.30, *extrapolated past the
  end of the Oxford ladder*) × a duration penalty (~0.83) ≈ 0.25.
- What φ *is*: definition-sensitivity = fraction of the inclusive Oxford-composite
  TD cases that also meet Hornick's strict sustained-high-fever bar.
- The bug [inferred]: 0.25 is the **low-dose asymptote** of φ(T), applied as a
  dose-**constant**. Likelihood is `y ~ Binomial(n, φ·p_fev)`, so fitted ≤ φ. But
  at saturating dose everyone develops florid disease → φ must → 1. Hornick raw
  (strict definition): 10³ 0/14, 10⁵ 32/116 (0.28), **10⁸ 8/9 (0.89), 10⁹ 40/42
  (0.95)** [from Hornick 1966]. So 0.25 is structurally impossible at high dose,
  not merely biased.
- Why it survived: only Hornick spans the dose range to saturation; Gilman/Levine
  are all single-dose (10⁵, φ=0.65) so the cap never bit there.
- Plan §7 had actually specified `φ ~ Beta(5,5)` (estimated); the §2.5 fixed
  0.25/0.65 "initial exploration" shortcut overrode it. The 0.25/0.65 split was
  unaudited LLM-derived. [Mike: "this is LLM drift that never got audited"]

## E2 — Reframe: the bug is "φ LOCKED," not "φ constant" [Mike]
- At saturation `p_fev → 1` so `fitted → φ`; the **Hornick high-dose plateau is a
  direct measurement of φ ≈ 0.9**. So a single *estimated* scalar suffices — no
  need for dose-dependent φ (which would cost +1–2 params on already-thin ID).
- δ is in Tier 1 and **covaries with φ** [Mike]: both suppress Maryland fitted
  fever, so high-φ solutions pair with high-δ solutions. Pinning φ low forced δ
  low → the φ-cap and the δ prior-data tension are plausibly the same artifact.
- Oxford anchors vs Hornick raw (the cross-era gap that δ + φ must explain):
  Oxford 10³ bicarb naive = 11/20 (55%), 10⁴ = 13/20 (65%); Hornick 10³ milk =
  0/14 (0%). Same nominal dose, ~55-pt gap, split between δ (milk delivers fewer
  viable organisms) and φ (definition). [observed]
- Decision: estimate a **single global scalar `phi_md`** (not per-study —
  Gilman/Levine single-dose can't identify their own φ; only Hornick's dose range
  can). Cross-study level differences → study RE (Step 2) / residual.

## E3 — Iteration 1: float `phi_md ~ Beta(5,5)` (plan §7)
Implementation: new `real<lower=0,upper=1> phi_md`; routed through `obs_prob()`,
`lprior`, priors.yaml/data_prep; φ dropped from data. Refit [observed]:
- 0/4000 div, R-hat ≤ 1.002, ESS > 1700.
- **φ̂ = 0.886** (90% CI 0.81–0.95). High-dose Hornick now fittable (H-F-8→0.74,
  H-F-9→0.79 vs old ≤0.25 ceiling).
- **δ rose 1.6 → 2.30** (δ ≈ 200×) toward its prior; priorsense flipped δ from
  "prior-data conflict" to "strong prior / weak likelihood" → **φ-cap and δ-tension
  confirmed to be the same artifact** (E2 prediction).
- `alpha_fevginf` re-identified: prior-sens 0.56 → 0.06 (post 0.91 → 0.24). The cap
  had been corrupting the fever-heterogeneity parameter.
- Problem: φ̂=0.886 sits at the **prior's 99.8th percentile** (Beta(5,5): mean 0.5,
  sd 0.15, 90% CI [0.25,0.75], ~10 pseudo-obs at 0.5). Prior too strong; priorsense
  flags φ conflict; H-F-9 still underfit because the prior holds φ down. [observed]

## E4 — Iteration 2: `phi_md ~ Beta(1,1)` (uniform) [Mike: 5,5 too strong]
Span check confirmed Beta(5,5) too informative (above). Loosened to uniform; "let
the plateau speak; informed prior held in reserve if convergence suffers." Refit
[observed]:
- 0/4000 div, R-hat ≤ 1.002 (φ_md ess_tail ~1250 — fine). Convergence held.
- **φ̂ = 0.966** (median 0.975, 90% CI 0.91–0.998) — floats to the plateau and
  **presses the [0,1] boundary**: the data wants ≈no definitional suppression.
- **δ rose further to 2.46** (δ ≈ 280×), still climbing toward the 3.5 prior.
- **KEY FINDING — the binding constraint shifted OFF φ.** High-dose Hornick is
  *still* ~0.12 underfit (H-F-8 0.89→0.77, H-F-9 0.95→0.83) **even at φ≈0.97**.
  Cause is no longer φ: the predicted plateau is capped at ~0.86 by the Maryland
  **mixture's immune component** (~39% at CoP_imm≈2.2) + the slow `P_fev|inf`
  saturation (small `alpha_fevginf`). A tighter φ prior would *lower* the plateau,
  not fix this. [inferred from the per-obs decomposition]
- **Hornick small-N caveat** [Mike, confirmed]: H-F-8 n=9 → Wilson 95% [0.57,0.98]
  *contains* the fit; H-F-9 n=42 → Wilson [0.84,0.99] *overlaps* the fit [0.76,0.89].
  So the "underfit" is consistent with sampling noise — **downgraded** from
  residual-to-chase to watch.

## E5 — Bespoke diagnostic figure
`dose_response_curves.R` → `results/tier1/dose_response_fit.png`, wired into the
driver. 3 panels (Oxford fever / Maryland fever / Maryland infection): posterior
median + 90% ribbon over a dose grid, observed points with Wilson 95% CIs, Stan
`p_pred` (×), and the φ̂ ceiling line. Makes the plateau/residual legible; curves
are an R visualization mirror (inference PPC still from Stan `p_pred`).

## E6 — Prior over-regularization: which params to watch [Mike Q: "flatter than ideal?"]
priorsense (Beta(1,1) fit), by *prior* sensitivity [observed]:
- **`alpha_inf` 0.29 — top.** Prior `lognormal(-1.5,0.8)` centers α at 0.22
  (shallow beta-Poisson); data wants steeper (post 0.39). **This is the "too flat"
  lever** — shallow α ⇒ slow rise toward the plateau.
- `N50_inf`/`N50_fevginf` 0.26/0.25, `delta` 0.21 — dose-scale block, prior-pulled.
- `gamma_inf` 0.16, "strong prior / weak likelihood" — barely data-informed (by
  design; needs Tier 2 Oxford shedding).
- `alpha_fevginf` now fine (0.08, data-identified).

## E7 — CoP taxonomy [Mike Q: are we fitting Hornick CoP? real vs guess vs fit?]
- **No Hornick-specific CoP.** Hornick + Levine + unstratified Gilman share ONE
  latent Maryland mixture, fit [observed]: `pi_susc` 0.61 (CI 0.38–0.82),
  `CoP_susc` 0.84 (0.61–1.11), `CoP_imm` 2.24 (1.21–3.82).
- `CoP_imm`/`CoP_susc` are **natural-scale** params (lognormal *prior*); they enter
  as **`CoP^gamma`** so protection is sub-linear (2.24^0.52 ≈ 1.5×, not 2.24×).
- Three tiers of CoP handling:
  - **real (definitional):** naive Oxford = CoP 1.0 (Waddington, Gibani; naive by design).
  - **guess (placeholder, treated as data):** non-naive Oxford — Darton placebo 1.15,
    Jin control 1.12, **Jin Vi-TT 5.0, Jin Vi-PS 2.0**. Hand-mapped, unfitted.
  - **fit (latent):** Maryland mixture (3 params); Gilman strata use *real* H-antibody
    data for component *assignment* but the CoP magnitude is fit.
- **Load-bearing caveat:** γ (immunity scaling) is identified almost entirely
  through the two placeholder Jin vaccine CoPs (5.0/2.0) → γ is conditional on
  guesses. This is the model's soft spot, and why γ_inf reads prior-dominated.

## E8 — Where did the Jin CoPs come from? [Mike Q]
- Jin **data is real**: attack rates (Vi-TT 13/37, Vi-PS 13/35) and post-vax
  anti-Vi IgG GMTs **Vi-TT 562.9, Vi-PS 140.5 EU/mL** [from Jin 2017].
- `CoP = 5.0 / 2.0` are **not in Jin** — eyeballed stand-ins for the output of an
  unbuilt `CoP = g(anti-Vi IgG)` titer model. CSV says verbatim "PLACEHOLDER —
  needs titer model mapping." The intended anchor is **Darton's HR = 0.29 per log₁₀
  anti-Vi** (specified in the plan as a γ constraint, not implemented). 5.0/2.0 are
  order-plausible and correctly ordered, but unfitted.

## E9 — Assay verification: both studies use VaccZyme [Mike directive]
- **Jin 2017** (Procedures, p.2474) [from paper]: "Anti-Vi IgG titres were measured
  using a commercial ELISA kit (**VaccZyme, The Binding Site, Birmingham, UK**)
  according to the manufacturer's guidelines." LLD 7.4 EU/mL (<LLD → 3.7 before
  log). The "adapted from VaccZyme" caveat applies *only* to the IgG **subclass**
  assay (IgG1/2/3, LLD 1.56), NOT the total anti-Vi IgG used for the GMTs.
- **Darton 2016** (Methods, p.4) [from paper]: "IgG responses to Vi were measured
  using a commercial ELISA kit (**VaccZyme, The Binding Site Ltd, Birmingham, UK**)
  according to the manufacturer's instructions." Same LLD 7.4 EU/mL; the HR 0.29 is
  "per 1 log₁₀ increase in anti-Vi IgG titre" on this scale.
- **Conclusion:** Jin and Darton are on the *identical* commercial VaccZyme anti-Vi
  IgG EU/mL scale (LLD 7.4). The GMTs are directly comparable; the Darton HR applies
  to Jin's titres with **no cross-assay calibration**. So the 5.0/2.0 placeholders
  are not broken by a scale mismatch — they're just an un-derived eyeball of a
  well-posed mapping.

## D1 — Decision (this moment): VaccZyme anti-Vi IgG EU/mL is the definitional CoP unit [Mike]
Adopt **VaccZyme anti-Vi IgG titre (EU/mL)** as the definitional immunity unit
throughout — the modern standard for typhoid anti-Vi. Rationale: both modern Oxford
inputs are on it; it replaces the dimensionless CoP + the 5.0/2.0 placeholders with
real, comparable, measured numbers.

Key realization for the link [inferred]: the existing `CoP^gamma` form, with **CoP
redefined as the titre ratio `anti-Vi / naive-ref`**, *is* a per-log₁₀-titre power
law — i.e., exactly the functional form Darton's "HR per log₁₀ anti-Vi" implies, with
**γ = the Darton protection slope**. So this is a re-denomination, not a new link
function. Reference titres (VaccZyme EU/mL):

| group | anti-Vi IgG (EU/mL) | basis |
|---|---|---|
| Jin Vi-TT | 562.9 | measured GMT |
| Jin Vi-PS | 140.5 | measured GMT |
| Jin control (day 28) | 8.0 | measured GMT |
| Darton placebo | ~baseline (per-subject in S1) | measured (S1) |
| Oxford naive (Waddington, Gibani) | <7.4 → 3.7 | LLD imputation |
| Maryland (Hornick/Levine/Gilman) | not measured (pre-VaccZyme) | **latent** anti-Vi-equiv EU/mL |

**Open choices reserved for Mike before the value-swap refit:**
1. naive reference level (3.7 vs 7.4 EU/mL) — affects interpretation of CoP=1, not
   the fit (absorbed into N50/γ under the power law).
2. Maryland latent prior on the *new* EU/mL-relative scale (currently lognormal
   centered ~1–2.7 on the old dimensionless scale — needs rescaling).
3. Darton placebo: per-subject S1 titres (avoids the Vi-PS Jensen bias the plan
   flagged) vs group GMT.

## Current model snapshot (end of session)
- Tier 1, 25 obs, 11 active params; clean (0 div). φ_md ≈ 0.97 (Beta(1,1),
  edge-pressing). δ ≈ 280×. Mixture: π_susc 0.61, CoP_susc 0.84, CoP_imm 2.24.
- Open threads: (a) elicit an informed φ prior? (φ edge-pressing); (b) high-dose
  Hornick residual = immunity/fever-saturation, *not* φ (and small-N anyway);
  (c) `alpha_inf` prior pulling curves flat; (d) δ↔N50 ridge; (e) γ weak pending
  Tier-2 Oxford shedding; (f) study RE = Step 2; (g) Gilman own-panel offered.
- **Immediate next:** wire VaccZyme EU/mL as the definitional CoP unit (D1), then
  the value-swap refit once the three reserved choices are set.

## D2 — Tier 1.5 decided: pull Darton individual data into the likelihood [Mike]

The "as much as possible" question (is Darton in the likelihood?) surfaced that the
91-subject individual-level Darton data — per-subject anti-Vi titres (EU/mL), 5
fever thresholds, bacteremia/stool — is extracted but used by **nothing**; only the
placebo group binomial (20/30, placeholder CoP) is in the fit. Decision [Mike]: this
is a **big shift**, filed for the record and built incrementally.

- **GitHub issue #15** (`famulare/typhoid-immune-dynamics`) — full individual-level
  upgrade scope + increment ladder.
- **Plan stage:** `calibration/tier1.5_plan.md`; ladder row added to `CALIBRATION_WORKFLOW.md`.
- **Today (identifiability goal):** build **minimal + φ** — EU/mL value-swap (Jin
  563/141/8, naive 3.7; Maryland latent CoP priors rescaled), Darton placebo as 30
  n=1 individual rows (`CoP_i = titre_i/3.7`, from the S1 extract), Jin group at real
  GMTs, then Darton multi-threshold → φ(T). The individual data makes identifiability
  *easier* (titre→fever slope as data, not the placeholder-CoP fudge), which is why
  we proceed now rather than waiting for the full upgrade.
- Reserved choices set as plan defaults: naive-ref 3.7 (fit-invariant); Maryland
  `CoP_imm ~ lognormal(log 5, 0.7)`, `CoP_susc ~ lognormal(0, 0.3)`; φ(T) exponential
  decay; stage C0+C1+C2 → fit → C3.
- Darton placebo gradient (sanity) [observed]: undetectable anti-Vi 72% fever (18) vs
  detectable 58% (12) — modest low-titre protection signal, n=30.

## Identifiability & prior-dependence status (Tier 1 → Tier 1.5)

Pathfinding snapshot for detailed review. All [observed] from priorsense
(prior-sensitivity / likelihood-sensitivity) on the clean fits. Tier 1 = Beta(1,1),
dimensionless CoP; Tier 1.5 = VaccZyme EU/mL axis + individual Darton + Jin.

**Strong / data-identified**
- **φ_md** — strongly data-driven (lik 0.26→0.38); but *edge-pressing* at the [0,1]
  boundary under the uniform prior (φ̂≈0.97). Being replaced by φ(T) in C3.
- **titre→fever MARGINAL slope** — newly data-informed in Tier 1.5 via individual
  Darton + Jin on the EU/mL axis: `gamma_fevginf` lik (0.085) now > prior (0.069),
  CI [0.045,0.29]. Combined γ_inf+γ_fevginf ≈ 0.36 ⇒ ~0.44 protection/decade,
  in-sample-consistent with Darton OR 0.37 / HR 0.29 (consistency, not out-of-sample
  validation — Darton is in the likelihood).
- **N50_inf / N50_fevginf** — pinned by Oxford naive fever (but entangled in the ridge).
- **δ (log10_delta)** — data-leaning in T1.5 (lik 0.22 > prior 0.13); δ̂≈315×.
- **pi_susc** — Gilman-informed (Beta(7,4)).

**Weak / prior-dependent (the open items)**
- **γ_inf** — *prior-dominated both tiers* (T1.5 prior 0.21 ≫ lik 0.065). The
  infection-vs-fever split of the immunity slope is unidentified without individual
  *infection* endpoints → the +cascade increment (issue #15). The marginal slope is
  pinned; the split is the prior.
- **CoP_imm** (Maryland latent immune anti-Vi-equiv) — *flipped to prior-dominated*
  in T1.5 (sens 0.14→0.43). **Not data-identifiable** (no 1960s serology). Highest-
  leverage **open prior elicitation**.
- **alpha_inf** — prior-leaning (sens 0.41). Infection-curve *steepness*; needs dose
  spread (Maryland via δ). The "flatter than ideal" lever Mike flagged.
- alpha_fevginf — was prior-dominated pre-φ-float; now data-leaning (lik 0.35).

**Open prior-elicitation / prior-choice items (for review)**
1. **CoP_imm / CoP_susc** — Maryland latent immune/susc anti-Vi-equiv (EU/mL).
   Unidentifiable; pure elicitation. Current (my weakly-informative picks):
   `CoP_imm ~ lognormal(log 5, 0.7)` (~5× naive ≈ 18 EU/mL), `CoP_susc ~ lognormal(0, 0.3)`.
   **Highest-leverage unelicited choice.**
2. **γ priors** recentered to median ~0.2 on the titre scale (my pick); γ_inf rides it.
3. **φ_md prior** Beta(1,1) uniform — edge-pressing; informed prior in reserve
   (C3 replaces the scalar with φ(T) from the Darton ladder).
4. **naive-ref 3.7 EU/mL** — fit-invariant under CoP^γ; low stakes.

**Structural identifiability limits (not prior)**
- **δ↔N50 ridge** (r ≈ −0.82) — dose-scale confound; needs dose spread; persists T1→T1.5.
- **γ_inf vs γ_fevginf split** — needs individual infection endpoints (+cascade).
- **φ dose-independence** — φ(T) calibrated at one (Oxford) dose; Hornick 39.4
  extrapolates beyond the 39.0 ladder max.

## C3 finding — the Darton ladder validates φ₀, and proves φ must be dose-dependent

Worked the Darton placebo temperature ladder (the +φ piece). Among the 20 TD+
placebo subjects: ≥38 = 16, ≥38.5 = 10, ≥39 = 8 [observed]. Two results:

1. **The original eyeballed φ values are right — as the LOW-DOSE asymptote.** The
   data-derived definition-sensitivities are φ(38.3) ≈ 0.66 and φ(39.4) ≈ 0.30
   (interp/extrap), essentially the hand-picked 0.65 / 0.25. So φ₀(T) was never the
   bug; *constancy across dose* was. The arc closes.
2. **A dose-INDEPENDENT φ(T) re-introduces the cap.** Plugging φ(39.4) ≈ 0.30 into
   the Maryland likelihood caps high-dose Hornick: in the current fit p_fev_mix(10⁹)
   ≈ 0.87, so fitted fever → 0.30·0.87 ≈ **0.26 vs observed 0.95** — the original
   catastrophe. The free `phi_md ≈ 0.97` only fits the saturating points *because*
   it's near 1; the ladder says φ should be 0.30 at that threshold. **That
   contradiction (φ≈0.30 at Oxford dose vs φ→1 needed at Hornick's saturating dose)
   is the data-grounded proof that φ is dose-dependent** [inferred, arithmetic].

**Decision [Mike]: bank dose-dependent φ(T,D) for next session** (fork option b).
The proper fix is the plan §2.5 form `φ(T,D) = φ₀(T) + (1−φ₀(T))·P_fev_naive(D)^β_φ`
— φ₀(T) from the ladder, → 1 at saturating dose. Retires `phi_md`, adds φ₀-level/λ/β_φ;
**β_φ rests on Hornick's 4 multi-dose fever points (thin-data ID — watch it).** A
dose-independent φ(T) is a dead-end (re-caps), so we go straight to dose-dependent.
Task-listed in `tier1.5_harness_handoff.md`; the harness parity/recovery guards must
re-sync to the new φ params after C3 (parity gate is GREEN for the current phi_md state).

## NEXT SESSION — pick up here (2026-06-24)

State: branch `dose-response-tier1-resurrection` (not merged to main). Tier 1.5
minimal+φ (C0+C1+C2) is committed and clean (0 div); the prior Tier-1 (Beta(1,1),
dimensionless) fit is preserved at `results/tier1_pre_eumL/` for comparison.
Full individual-Darton upgrade tracked as **issue #15**.

**The main task — C3: dose-dependent φ(T,D)** (banked tonight; see C3 finding above):
- Implement plan §2.5 form: `φ(T,D) = φ₀(T) + (1−φ₀(T))·P_fev_naive(D_eff)^β_φ`.
  φ₀(T) = Darton-ladder asymptote (≈0.30 at 39.4 / 0.66 at 38.3); → 1 at saturating dose.
- Retire `phi_md`; add φ₀-level + λ (ladder) + β_φ (dose-rise). **Watch β_φ — it rests
  on Hornick's 4 multi-dose fever points only (thin ID).**
- Per-obs temperature threshold needed in the data (Hornick 39.4, Levine/Gilman 38.3);
  map by study in `data_prep.R`. Darton ladder counts (n_TD=20, ≥38=16, ≥38.5=10, ≥39=8)
  feed the φ₀ sub-likelihood.
- **Then re-sync the harness guards** (they'll go red — that's the guard working):
  parity `PARAM_NAMES`/`vecs`/`obs_prob_R` (replicate φ(T,D)) and recovery
  `PARAM_NAMES`/`TRUTH_REALISTIC`. See `tier1.5_harness_handoff.md`.

**Other dangling bits:**
- **Plots**: add a titre→protection (CoP-axis) panel to `dose_response_curves.R` — the
  current dose-axis panels stack the 30 Darton individuals at one dose (y∈{0,1}), which
  doesn't show the immunity slope we gained.
- **Recovery `TRUTH_REALISTIC`** still on the pre-EU/mL scale → re-sync (delta~2.5,
  CoP_imm~5–7, gamma~0.2, alpha_inf~0.4) so point recovery tests a representative point.
- **Open prior elicitation (your call)**: `CoP_imm` (Maryland latent immune anti-Vi-equiv)
  is prior-dominated + unidentifiable — highest leverage. Also confirm the reserved
  picks: naive-ref 3.7 (fit-invariant); Maryland CoP priors (`CoP_imm ~ lognormal(log5,0.7)`,
  `CoP_susc ~ lognormal(0,0.3)`); γ recenter to median ~0.2 (titre scale).
- **Detail review** you flagged: the EU/mL value-swap (Jin 563/141/8→CoP, Darton
  individual rows), the γ_inf still-prior-dominated split, the δ↔N50 ridge.

**Don't re-derive:** VaccZyme verified for both Jin & Darton (EU/mL, LLD 7.4); the
φ-cap↔δ-tension-same-artifact result; φ₀ validates the original 0.25/0.65.

---

## C3 IMPLEMENTED — dose-dependent φ(T,D) (2026-07-15) [observed]

Implemented the banked `φ(T,D) = φ₀(T) + (1−φ₀(T))·P_fev_naive(D_eff)^β_φ`. Design
decisions this session [Mike]: **β_φ floats** (not fixed=1); **CoP_imm ~ Exponential,
mean 50 EU/mL absolute** (= 13.5 on the CoP ratio scale, rate 0.074) replacing the
lognormal. Wiring: φ₀(T)=inv_logit(phi0_a−phi0_b·(T−T_ref)), T_ref=38.0, pinned by a
**decoupled** Darton-ladder binomial sub-likelihood (16/10/8 of 20 at ≥38/38.5/39);
β_φ + Hornick multi-dose supply the dose-lift. Retired `phi_md`. `exponential` family
added to priors.R. Per-obs threshold mapped by study in data_prep.R. Parity gate GREEN
(max |Δp|=5e-9 over 3 vecs × 54 rows). Harnesses + fit driver + curves re-synced.

**Fit (4×1000, adapt_delta 0.9): CLEAN — 0/4000 div, R-hat ≤ 1.004, E-BFMI 0.74.**

- **φ(T,D) works.** φ(39.4, D) = 0.31 → 0.70 → 0.97 across 10³/10⁵/10⁹; φ₀(39.4)=0.22.
  So φ→~1 at saturation (resolves the constant-φ cap) while φ₀(38.3)=0.69 / φ₀(39.4)=0.22
  reproduce the hand-derived 0.66/0.30. The C3 arc closes as designed.
- **β_φ = 0.76 (median 0.67, 90% CI 0.27–1.59).** Floated, CI spans 1, **no pathology**
  (the thin-ID risk did not bite the geometry). BUT priorsense: **β_φ prior-dominated**
  (prior 0.20 ≫ lik 0.055) — the dose-lift *shape* is essentially the prior; the data
  barely constrains it. Honest limitation, as flagged. β_φ<1 lean = φ rises slightly
  faster than the naive fever curve.
- **φ₀ data-identified** (phi0_a lik 0.098>prior 0.040; phi0_b lik 0.085>prior 0.054) —
  the ladder does its job. ess_tail ~440–635 (only 3 counts; fine).
- **δ now likelihood-leaning** (lik 0.31 > prior 0.25), δ̂≈143–209× (log10_delta 2.17,
  DOWN from the minimal+φ ~2.46 that was climbing toward the 3.5 prior). Dose-dependent
  φ further unwound the δ↔φ entanglement. [improvement]
- **CoP_imm prior-dominated** (prior 0.37 ≫ lik 0.050); posterior mean 14.2 ≈ prior mean
  13.5 → confirmed **unidentified**, now cleanly carried by the Exponential prior.
- **alpha_inf** still prior-leaning (prior 0.38, lik 0.12; the "too flat" lever), γ_inf
  prior-dominated (needs +cascade) — both expected/open.

**Residuals to watch (NOT C3 regressions):**
1. High-dose Hornick marginal still ~0.09 underfit: H-F-8 obs 0.889 / fit 0.786; H-F-9
   obs 0.952 / fit 0.860. Even at φ≈0.97 the fit is capped by the mixture immune
   component + fever saturation (the E4 finding), within Hornick small-N Wilson. Slight
   improvement vs minimal+φ (H-F-9 0.83→0.86).
2. **NEW mild tension: Hornick conditional H-FgI-7** (P(fever|inf) at 10⁷) obs 0.571 /
   fit 0.761 — model over-predicts the conditional. φ(39.4,10⁷)≈0.85 pushes fever|inf
   high vs the observed 0.57. One point (n=28); flag for the +cascade increment (which
   adds individual infection endpoints and may re-balance the inf/fever split at 10⁷).
3. Gilman susceptible stratum Gil-F-Hlo obs 0.611 / fit 0.420 (underfit); Levine trial
   scatter (0.25–0.55 around a pooled 0.358) → study RE (Step 2). Pre-existing.

Backups: minimal+φ fit preserved at `results/tier1_minimal_phi/`; pre-EU/mL at
`results/tier1_pre_eumL/`. C3 fit at `results/tier1/`.

**Next (issue #15 ladder): +cascade** — Darton per-subject bacteremia/stool → split
γ_inf vs γ_fevginf (γ_inf still prior-dominated); may also relieve residual #2.

---

## β_φ PINNED=1 + +cascade (2026-07-15) [observed]

Two changes batched into one refit [Mike: "c" (pin β_φ) + continue +cascade]:

**(1) β_φ pinned to 1.** The C3 fit showed β_φ prior-dominated (data can't identify the
dose-lift *shape* from 4 Hornick points). Simplified to `φ(T,D)=φ₀(T)+(1−φ₀(T))·P_fev_naive`
— retains the correct limits (φ₀→1 at saturation) + inherited dose scale, drops the
un-identified shape param. One fewer parameter.

**(2) +cascade.** Darton placebo per-subject endpoints as a PROPER cascade: group 6
(ox_inf_indiv, P_inf, y=bact_or_stool) for all 30; group 7 (ox_fevginf_indiv, P_fev|inf,
y=fever_td) for the 26 infected. Replaces the C1 composite-fever rows (which conflated
the layers). Infection=bact_or_stool (broadest; no eta at Tier 1 — flagged). No new
params. 80 obs total (was 54). Parity GREEN (max |Δp|=5e-9 × 80 rows × 3 vecs).

**Fit: 1/4000 div (0.03%, noise), R-hat ≤ 1.004, ESS fine.** vs C3 (isolating cascade):
- **γ_inf prior-dependence HALVED** (priorsense prior sens 0.168→0.075; now ≈ lik 0.056,
  roughly balanced). γ_inf 0.182 [0.057,0.316]. The individual infection endpoints pulled
  it off the prior — the +cascade mechanism working.
- **γ_fevginf now likelihood-leaning** (prior 0.038 < lik 0.054, diagnosis "-"). Data-identified.
- **alpha_inf prior-dependence halved** (0.378→0.204), CI tightened 0.32→0.28 [0.18,0.46] —
  addresses the "too flat α" open item; individual infection data pinned the steepness.
- **δ likelihood-leaning** (lik 0.28 > prior 0.23), dropped to ~114× (log10_delta 2.06).
- Individual endpoints fit well: ox_inf lo/hi-titre obs 0.89/0.83 → fit 0.89/0.80;
  ox_fevginf 0.81/0.70 → 0.85/0.76 (captures the weak immunity gradient).

**HONEST LIMITATION: the γ_inf vs γ_fevginf SPLIT only weakly resolves.** Both sit ~0.16–0.18
(near the prior median 0.2); γ_inf is not cleanly < γ_fevginf despite the raw data hint
(inf 0.89→0.83 weaker than fev|inf 0.81→0.70). Cause: Darton's titre range is thin
(median 3.7, only 12/30 detectable, max ~62) — not enough immunity spread to separate the
two slopes. +cascade reduced prior-dependence + sharpened α/δ, but SEPARATING the γ's needs
a longer titre axis → **+vaccine-terms** (M01ZH09/Ty21a arms) and/or **+Jin-digitize**
(high-titre Vi points). CoP_imm still prior-dominated (unidentified, Exp prior carries it).

Backups: C3 (float β_φ) results superseded (committed as code d1880a8); minimal+φ at
`results/tier1_minimal_phi/`; current +cascade fit at `results/tier1/`.

## Consolidation + NEXT SESSION (2026-07-15 close) [Mike: "c for now, b tomorrow"]

Consolidated (no new modeling increment):
- **Titre→protection figure** `results/tier1/titre_protection.png` (new `plot_titre_protection()`
  in `dose_response_curves.R`, wired into the driver). CoP-axis, 3 facets (P_inf /
  P_fev|inf / composite fever) with the CoP^γ curve + Darton individuals (grey, jittered
  0/1) + Jin vaccine groups (blue, Wilson). **Shows the story visually:** Darton clusters
  at low titre (3.7–60), Jin anchors the high end (141/563) on the composite curve and
  sits on the ribbon; the P_inf and P_fev|inf slopes are shallow and near-identical →
  the γ-split is titre-range-limited. Motivates +Jin-digitize.
- **model_structure.md** mermaid diagram refreshed (groups 1..7, phi(T,D), ladder, cascade).
- Harness handoff (`tier1.5_harness_handoff.md`) closed out.

**NEXT SESSION — +Jin-digitize (issue #15 stretch) [b, tomorrow]:**
- Digitize Jin 2015 **Fig S3** (per-subject anti-Vi titre → outcome) to add Jin vaccine-arm
  subjects as INDIVIDUAL rows on the anti-Vi axis (CoP up to ~152). Currently Jin is only
  3 GROUP fever points (ox_fev, CoP 2.2/38/152). Individuals would populate the high-titre
  end of the P_inf and P_fev|inf facets → the spread needed to SEPARATE γ_inf vs γ_fevginf
  (the +cascade split that stalled on Darton's thin low-titre range).
- Route Jin individuals through the same cascade groups (6 = infection, 7 = fever|inf) if
  Fig S3 gives both endpoints per subject; else fever-composite (group 1) individuals.
  Jensen-bias caveat (group GMT vs individual) is exactly what this retires.
- No new parameters expected; re-run parity + fit; check whether γ_inf and γ_fevginf
  separate and whether the Hornick conditional over-prediction (residual #2) eases.
- **Don't re-derive:** β_φ pinned=1 (was prior-dominated); CoP_imm Exp(mean 50 EU/mL abs)
  unidentified→prior-carried; the cascade decomposition (g6/g7) + parity are GREEN.

Commits this session: d1880a8 (C3 dose-dependent φ), 0088b9e (+cascade + β_φ=1 pin),
+ consolidation (plots + docs). Branch `dose-response-tier1-resurrection`, not pushed.

## +Jin-CoP-anchor — published OR informs the γ prior (in lieu of digitization) [Mike]

**Decision:** Jin Fig S3 (individual anti-Vi titre → diagnosis probability) digitization
is **feasible but NOT done** — it needs WebPlotDigitizer + human axis calibration (the
individual points exist only in the scatter; Jin tabulates no per-subject titres). Two
reasons it's low-value anyway: (a) the figure's quantitative content is *already a
published number* (the logistic OR); (b) it is titre→**diagnosis** (composite TD), so it
would tighten the combined γ, NOT separate γ_inf vs γ_fevginf. Instead we **use Jin's
published logistic result to inform the γ prior.**

- **Jin 2017 adjusted OR = 0.37 (95% CI 0.15–0.88) per log₁₀ anti-Vi IgG (EU/mL)** for
  typhoid diagnosis [from paper, p.2477]. Increment **VERIFIED**: Fig S3 x-axis is
  "Log10 anti-Vi IgG (ELISA units per mL)" (pypdf text extract of the main paper).
  Unadjusted OR 0.35 (0.21–0.59). Same VaccZyme EU/mL scale as the model.
- **Map OR → γ:** single-layer beta-Poisson at the control op point (P≈0.77, u=−ln(1−P)=1.47):
  `d logit(P)/d log₁₀CoP = −γ·ln(10)·u/(1−e^−u) = −γ·4.40`. OR 0.37 ⇒ γ_eff = ln(0.37)/−4.40
  = **0.23**; CI ⇒ γ_eff ∈ ~[0.03, 0.43]. Both Jin (0.37) and Darton (HR 0.29/log₁₀)
  land γ_eff ≈ 0.2–0.3 for the COMPOSITE slope.
- **Prior update** (`priors.yaml`, no recompile — hyperparams are data): γ_inf & γ_fevginf
  `sdlog 0.9 → 0.7`, median unchanged at 0.20 (two concordant published anchors justify the
  mild tightening; 95% ~[0.05,0.8] ≈ OR-implied range). Anchor is on the COMBINED slope,
  expressed as concordant priors on both layers — the C3 fit already reproduces Jin's
  composite with γ_inf~0.20/γ_fevginf~0.15, so this doesn't fight the data.

**If we later want the actual scatter:** WebPlotDigitizer pass (human) or locate deposited
TyVAC/VAST trial data; wire extracted points as individual rows (I integrate, don't invent).
Tracked under issue #15 (+Jin-digitize remains open as the higher-fidelity option).

## +vaccine-terms — Darton VE pulled + design documented (2026-07-31) [Mike]

Pulled Darton 2016 Table 2 vaccine efficacies (M01ZH09, Ty21a); full table + increment
design in `tier1.5_plan.md` ("+vaccine-terms increment (planned)"). Headline [from extract]:
M01ZH09 weak/non-significant (VE 13–28%, CIs cross 0); **Ty21a moderate + significant** on
infection/bacteraemia (bact-or-stool 38% [12,57]); and **adjusting for baseline anti-Vi
barely changes VE** → protection is non-anti-Vi-mediated. Proposed: add the two arms as
cascade rows (groups 6/7) with a separate per-vaccine factor `V_v` in the exponent
(`−alpha/(CoP^gamma·V_v)`), keeping vaccine protection OUT of γ; Ty21a CoP=1 (anti-Vi NA).
Verdict: parked — it does NOT extend the anti-Vi axis (no γ-split help); real value is the
anti-Vi-vs-cell-mediated protection *decomposition*, a distinct aim.

## PPC representation — individual rows pooled to arm-level rates (2026-07-31) [Mike]

**Problem [Mike]:** `ppc.png` mixed two representations on one axis. The 56 Darton
individual rows (`ox_inf_indiv` n=30, `ox_fevginf_indiv` n=26) are n=1, so their "observed
attack rate" is exactly 0 or 1 — they piled up on the x=0 and x=1 edges as label spaghetti
while every other marker was a grouped-binomial arm rate. Not interpretable side by side.

**Fix** (`diagnostics.R`, `.ppc_rows()` / `.ppc_quantile_bins()`): `*_indiv` groups are
pooled to arm-level dot-and-whisker like everything else. Observed = `sum(y)/sum(n)`;
fitted = the **n-weighted mean of the per-subject `p_pred` within each draw**, then 5/50/95
quantiles. That is exactly the quantity `p_pred` already reports for a grouped row, so
pooled and grouped markers are on the same footing. Per-subject rows are retained in
`ppc.csv` under `level == "individual"` (nothing lost); only `level == "arm"` is plotted.

**CoP stratification, and why 3 terciles are not achievable [observed]:** whole-arm pooling
averages away the titre spread that is the whole reason the individual rows exist, so each
arm is split into `indiv_bins = 3` ascending CoP strata. It realises as **2**, because
**18/30 (inf) and 16/26 (fev|inf) subjects sit exactly at the `<LLD` imputation floor
CoP = 1** (`vi_igg_prechallenge` → `NAIVE_VI_REF`). Subjects sharing a CoP get identical
model probabilities, so splitting that tie block would manufacture two markers differing
only by an arbitrary partition of one stratum. Binning is therefore tie-safe (ties never
split; bins under `indiv_bin_min = 4` merge into their smaller neighbour), and the realised
cut is **naive-floor vs detectable anti-Vi** — the informative contrast anyway:

| marker | n | obs | CoP range | fitted (median, 90%) |
|---|---|---|---|---|
| `D-I-plac q1` | 18 | 0.889 | 1.0 – 1.0 | 0.886 [0.825, 0.943] |
| `D-I-plac q2` | 12 | 0.833 | 2.14 – 16.9 | 0.796 [0.724, 0.863] |
| `D-FgI-plac q1` | 16 | 0.812 | 1.0 – 1.0 | 0.850 [0.790, 0.905] |
| `D-FgI-plac q2` | 10 | 0.700 | 2.40 – 16.9 | 0.755 [0.683, 0.825] |

Gradient sign is right in both endpoints and all four bins straddle the diagonal.
**Limit:** this cannot test a gradient *within* the floor — it is censored, not resolved.
The continuous titre check stays in `titre_protection.png`.

**Also fixed [observed]:** `simulate_recovery.R:165` passed `obs = obs` (the REAL y) to
`diagnose_fit()` while the fit was to `y_sim`, so the recovery `ppc.png` compared real
observed rates against synthetic-truth predictions. `obs_syn` was already built two lines
above for exactly this and used for `extra_plots` only → now `obs = obs_syn`. Existing
`results/recovery*/ppc.png` are stale until refit.

**Regeneration:** `ppc_from_dir(run_dir)` (mirrors `figures_from_dir()`) rebuilds
`ppc.png` + `ppc.csv` from a run dir's `fit.rds` + `stan_data.rds`, no refit; `plot_ppc()`
now hard-errors if `nrow(obs) != ncol(p_pred)`. Done for `results/tier1` and
`results/tier1_prior`. `results/tier1_minimal_phi`, `results/tier1_pre_eumL` and the five
`results/scenarios/*` have `fit.rds` but **no `stan_data.rds`**, so their row selection is
not recoverable — their `ppc.png` stay old-style until refit.

**Provenance note:** the `diagnostics.R` half of this change landed inside commit `68c799b`
("Add cohort_id…") and the `simulate_recovery.R` fix inside `7c0c48b` ("Model figure
suite…") — swept in by a concurrent session committing the shared worktree. `git log`
on those files will mis-attribute; this entry is the record.
