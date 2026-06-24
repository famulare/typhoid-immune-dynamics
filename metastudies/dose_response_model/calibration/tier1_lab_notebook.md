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
