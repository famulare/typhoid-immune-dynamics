# Tier 1 Stan calibration — resurrection + pathology diagnosis

**Date:** 2026-06-23
**Author:** Mike Famulare + Claude (Opus 4.8)
**Status:** Model now COMPILES and SAMPLES; first Tier-1 fit is NOT yet valid
(divergence pathology diagnosed, root cause isolated, fix deferred to iteration).

This documents what it took to get `typhoid_dose_response.stan` to run at all, and
the diagnosis of why the first run does not yet produce a trustworthy posterior.
It is a checkpoint, not a finished calibration.

---

## 1. What was fixed to make it run (mechanical resurrection)

The committed Stan file was an explicit "UNTESTED DRAFT" and had never compiled.
Three changes were required:

1. **Compile blocker — `_lp` suffix misuse.** `beta_poisson_lp` returned a plain
   probability but used Stan's *reserved* `_lp` suffix, which is illegal to call
   from a non-`_lp` function (`maryland_mixture`) and from `generated quantities`.
   Renamed `beta_poisson_lp` → `beta_poisson` everywhere. `[observed: stanc error]`
2. **Likelihood double-count — Hornick conditional.** The old line
   `n_hornick_cond ~ binomial(30, p_inf)` re-added the 28/30 infection datum that
   is already carried by the `md_inf` row `H-I-7` (reviewer-2 concern MC2).
   Removed; the conditional `16 ~ binomial(28, p_cond)` is kept. `n_hornick_cond`
   (=28) is now used only as the conditional denominator.
3. **Added a `prior_only` data flag** wrapping the likelihood, for prior
   predictive checks.

Toolchain installed: `cmdstanr` 0.9.0 + CmdStan 2.39.0 (R 4.5.3 via pixi).
The driver is [fit_dose_response.R](fit_dose_response.R); it assembles the 25
`tier1_active` observations from [dose_response_data.csv](dose_response_data.csv)
into the Stan groups (ox_fev=7, md_fev=11, md_inf=6, hornick_cond=1) with Oxford
shedding disabled (`N_ox_inf=0`).

The data is clean and maps mechanically; no data work was needed for this step.

---

## 2. The symptom

First Tier-1 fit (4 chains, 1000+1000, `adapt_delta=0.9`):
**3959 / 4000 (99%) divergent transitions.** Sampler finishes in ~1 s, R-hat
1.02–1.09, ESS 35–175. So it "runs" but does not explore the posterior — the fit
is not usable.

---

## 3. Diagnosis: elimination ladder

Each row is an independent 4-chain refit; `div` = divergent transitions / 4000.
`[observed]` — all numbers are measured, not inferred.

| Experiment | div / 4000 | Reading |
|---|---:|---|
| baseline (`adapt_delta=0.9`) | 3959 | — |
| `adapt_delta=0.99` | 3981 | **not** a step-size issue (worse) |
| φ = 1 for all md_fev (remove φ cap) | 3715 | φ cap is **not** the driver |
| drop H-F-8, H-F-9 (saturating points) | 3918 | saturating points **not** the driver |
| remove `2^(1/α)` term from scale | 3728 | beta-Poisson exponent **not** the driver |
| **Oxford fever ONLY** (no mixture, no δ, no φ) | **3422** | pathology is in the **core**, not Maryland |
| + Maryland fever (mixture) | 3444 | mixture adds little |
| + Maryland infection | 3990 | second aggravator (see §5) |
| **remove ordering constraint, Oxford-only** | **0** | **root cause** |
| **remove ordering constraint, full Tier 1** | **0** | **root cause confirmed** |

The decisive finding: Oxford fever *alone* — 7 binomial points, naive-ish, no
mixture, no cross-era machinery — already diverges 85%. So the cause is in the
core parameterization shared by every sub-model. Removing one line takes full
Tier 1 from 99% → **0%** divergent.

---

## 4. ROOT CAUSE — ordering-constraint density cliff

[typhoid_dose_response.stan L132](typhoid_dose_response.stan):

```stan
(log10_N50_fevginf - log10_N50_inf) ~ normal(0, 1) T[0, ];
```

This is a truncated sampling statement (`T[0,]`) on the **difference of two
unconstrained parameters**. The truncation makes the log-density **−∞** whenever
the difference goes below 0 — a hard **cliff**, not a smooth penalty. HMC/NUTS
integrates Hamiltonian trajectories; hitting a −∞ wall makes the leapfrog
integrator diverge. Because both `log10_N50_*` are free reals whose typical set
straddles the boundary, nearly every trajectory hits the cliff. `[inferred from
the 0-divergence result; mechanism is standard Stan behavior]`

**Recommended fix (for iteration, not yet applied):** impose the ordering
*smoothly and hard* via a non-negative offset parameter:

```stan
parameters {
  real log10_N50_inf;
  real<lower=0> d_fev;            // log10 gap, fever threshold above infection
}
transformed parameters {
  real log10_N50_fevginf = log10_N50_inf + d_fev;
}
model {
  d_fev ~ normal(0, 1);          // half-normal via the lower=0 bound
}
```

The `<lower=0>` bound is handled by Stan's unconstraining transform (log scale),
so there is no cliff — the constraint is built into the geometry. This preserves
the design intent (N50_fev|inf ≥ N50_inf) with zero divergences.

**Left active in the committed model on purpose** (with a flag comment at L132)
so the reparameterization is a deliberate joint decision, not a silent rewrite of
the scientific model.

---

## 5. SECONDARY findings (independent of the divergence fix)

These do **not** cause the divergences but will bias or shape the posterior once
sampling works. They need attention before any Tier-1 result is trusted.

### 5a. Constant φ = 0.25 caps Hornick fever below the data `[observed]`
The Maryland fever likelihood is `y ~ binomial(n, φ · p_fev)`, so the maximum
achievable fitted rate is φ. With φ = 0.25 (Hornick), three points exceed the cap:

| obs | dose (milk CFU) | observed rate | φ cap |
|---|---:|---:|---:|
| H-F-5 | 1e5 | 0.276 | 0.25 |
| H-F-8 | 1e8 | **0.889** | 0.25 |
| H-F-9 | 1e9 | **0.952** | 0.25 |

At saturating dose essentially everyone develops fever under *any* definition, so
φ must → 1 at high dose; a **constant** φ across Hornick's 6-log dose range is the
problem. This is exactly the design's own flagged risk (reviewer-2 MC1:
"constant φ approximation dangerous across the Maryland dose range; use
dose-dependent φ(D_eff)"). The `eta_fixed_optC` column hints the same machinery
exists for shedding. **Decision needed:** dose-dependent φ(D_eff), or revisit the
φ values / definition mapping.

**RESOLVED 2026-06-23 — floated φ as an estimated scalar (not dose-dependent).**
Key reframing (Mike): the bug is not "φ constant" but "φ *locked* at its low-dose
asymptote." At saturating dose `p_fev → 1` so `fitted → φ`, i.e. Hornick's high-dose
plateau (0.89, 0.95) is a *direct measurement* of φ ≈ 0.9 — so a single **estimated**
scalar suffices and is identified by the dose range; dose-dependent φ (+1–2 params)
isn't required and would strain the thin ID budget. φ is estimated as a single
*global* scalar (not per-study: Gilman/Levine are single-dose at 10⁵ and cannot
identify their own φ; only Hornick's dose range can — plan §7 had specified an
estimated φ all along, which the §2.5 fixed 0.25/0.65 overrode).

Prior went through two stages. **`Beta(5,5)`** (plan §7) gave φ̂ = 0.886 but pulled
it to the prior's 99.8th percentile and still left H-F-9 underfit, so it was loosened
to **`Beta(1,1)` (uniform)**. Final re-fit (Beta(1,1)): 0/4000 divergences,
**φ̂ = 0.966 (median 0.975, 90% CI 0.91–0.998)** — φ floats to the plateau and
*presses the [0,1] boundary* (the data wants ≈no definitional suppression). Bonus, as
predicted: freeing φ pulled `log10_delta` 1.6 → 2.46 (δ ≈ 280×, toward its prior) and
re-identified `alpha_fevginf` (priorsense prior-sensitivity 0.56 → ~0.06) — the cap
had been corrupting both.

**Residual — the binding constraint shifted OFF φ.** High-dose Hornick is *still*
~0.12 underfit (H-F-8 0.889 → fit 0.77; H-F-9 0.952 → fit 0.83) **even at φ̂ ≈ 0.97**.
The model's predicted plateau is capped at ~0.86 by the Maryland mixture's immune
component (~39% at CoP_imm ≈ 2.2) plus the slow fever-given-infection saturation
(small `alpha_fevginf`), not by φ. So: (a) a tighter/elicited φ prior addresses the
boundary-pressing/identifiability but would *lower* the plateau, not fix this residual;
(b) the high-dose misfit is a separate structural question — is immunity fully
overwhelmed at 10⁹ (Hornick 1970 reports ~100%)? does fever|infection saturate faster?
Visual: `results/tier1/dose_response_fit.png`. See `CALIBRATION_WORKFLOW.md` Step 1.

### 5b. `delta`–`N50` ridge `[observed]`
Posterior correlations: `log10_delta` vs `log10_N50_inf` ≈ −0.74, vs
`log10_N50_fevginf` ≈ −0.80. The Maryland likelihood depends on
dose/(δ·N50), so δ and N50 are separated only by the Oxford data (δ=1 there).
With ~7 Oxford fever points this is weak — matches the design note "Oxford alone
is underconstrained" (git `f6b7c91`). Not a blocker once the cliff is fixed, but
it governs δ/N50 identifiability and the prior–data tension below.

### 5c. Prior–data tension on `delta` `[observed, flagged]`
Even in the (divergent, so provisional) fit, `log10_delta` posterior sits near
~1.7–1.8 while its prior is `normal(3.5, 0.7)` — the data pulls δ to ~50–120×,
the prior wants ~3000×. Re-check after the cliff fix; if it persists, either the
δ prior or the Maryland/Oxford dose-scale bridge needs reconsideration.

**Update 2026-06-23 — substantially relaxed by floating φ (see §5a).** With the
cliff fixed *and* φ floated, `log10_delta` rose to **2.30** (δ ≈ 200×) and
priorsense flipped it from "prior-data conflict" to "strong prior / weak
likelihood." The tension was largely an artifact of the locked-low φ: φ and δ
covary (both suppress Maryland fitted fever), so pinning φ at 0.25 forced δ low.
Not fully resolved (still below the 3.5 prior mean), and now entangled with the
δ↔N50 ridge (§5b), but no longer a standalone conflict.

### 5d. Inert parameters in this Step-1 config
`sigma_study`, `eta_lo`, `kappa` are declared and given priors but enter no
likelihood term in Tier 1 (study RE unimplemented; η is Tier 2). They sample
their priors — do not interpret them. `sigma_study` belongs in repo-canonical
Tier 1 (Step 2 of the workflow).

---

## 6. State left for the next iteration

- **Fix the cliff** (§4) — should be first; turns the model from unusable to
  clean-sampling with a known-good reparameterization.
- **Resolve φ** (§5a) — required before any Tier-1 number is credible.
- Then re-examine δ identifiability/prior (§5b–c), add the study RE (Step 2), and
  proceed to Tier 2 (Oxford shedding + η). See [CALIBRATION_WORKFLOW.md](CALIBRATION_WORKFLOW.md).

Ruled out as divergence causes: `adapt_delta`, the φ cap, the saturating Hornick
points, the `2^(1/α)` term. The single cause is the ordering-constraint cliff.
