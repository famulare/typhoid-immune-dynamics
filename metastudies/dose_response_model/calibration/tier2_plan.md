# Tier 2 plan — η-correction + ψ-correction + Oxford shedding, on top of `t1-indiv-vax`

**Status:** LOCKED 2026-07-31 [Mike]. Builds on `tier1.5_plan.md` (`t1-indiv-vax`,
stage `phi-rho-vax`, N_obs=177). Target: new tier `t2-indiv-vax`, stage
`phi-rho-eta-psi-vax`, N_obs=182 (177 + 5 restored Oxford `ox_inf` grouped rows).
Also updates `t2-indiv` (individualized, non-vax) the same way; `t2-grouped` is
**retired, not built out** — Tier 2 is indiv-only going forward, no new grouped
configuration (decision below).

**Decisions locked 2026-07-31:**
1. Design decision A: **multiplicative**, `psi_late = psi_stool * frac_late`,
   `frac_late` gets its own weakly-informative prior. Ordering (`psi_late <=
   psi_stool`) is structural, not an estimated-then-checked constraint.
2. Design decision B: explicit `psi_active` (0/1) Stan data flag, threaded through
   `tier_specs.R` like `include_vaccine_arms`, as the third gate kind in
   `derive_stage_token()` alongside `needs_group2`/`needs_vaccine_rows`.
3. **No new `t2-grouped`.** Tier 2 is individualized-Darton only: `t2-indiv` and
   `t2-indiv-vax`. `t2-grouped` stays in the registry as a **retired** entry (history
   preserved, `status_declared` changed to reflect "will not be built," not deleted
   outright — matches the repo's own principle that a blocked/retired rung keeps its
   reason rather than disappearing).

Two open scientific tensions from `TIER_LADDER.md`'s `t2-grouped`/`t2-indiv` block
are being resolved by **explicit decision, not by code fix**:
1. η Option A (parametric, monotone-decreasing-in-dose) vs Option C (fixed,
   non-monotone across studies) — **decision: proceed with Option A**, which is
   already coded (`eta_detection()`, `typhoid_dose_response.stan:98`). See "The η
   confound" below for the caveat this carries forward, undissolved.
2. `joint_inference_plan.md` §2.6 excluding Oxford shedding on treatment-truncation
   grounds vs the Tier 2 design restoring it via η — **resolved by re-reading**: η
   *is* the treatment-truncation correction, so restoring the data with η attached
   is not in tension with §2.6's concern, it's the direct implementation of it. §2.6
   needs a one-line edit noting it's superseded/operationalized by §2.8's η, not
   contradicted. (Doc drift — see "Doc drift to fix" below.)

ψ (`joint_inference_plan.md` §2.8) was adopted 2026-07-31 but never implemented —
that's the other half of Tier 2, and the harder engineering piece.

## The η confound [Mike, flagged 2026-07-31]

`eta_fixed_optC` values by study/dose: Waddington 1e3→1.00, Waddington 1e4→0.62,
Darton 18200→0.94, Jin 20000→0.92. **Within** Waddington (the only paired-dose
study) the two points are monotone decreasing, consistent with Option A. The
apparent non-monotonicity only appears **across** studies at the higher doses, where
dose is confounded with per-study treatment/diagnosis protocol (τ_treat differs by
trial, not just physiology). So Option A isn't falsified by the raw curve shape —
but a single global η(D) fit across 3 studies at effectively 3 dose levels (one
study contributes 2) risks absorbing between-study protocol variation into what
looks like a dose effect. **This is accepted and documented, not resolved** — same
treatment as `CoP_imm`/ψ_late being prior-carried. Falsifiable read: if `kappa`'s
posterior is mostly prior-carried (priorsense) or trades off hard against per-study
structure we don't have parameters for, that's this confound showing up in the fit,
not a bug.

## 1. Restore Oxford `ox_inf` (η already implemented, mechanics already work)

`eta_lo`/`kappa` and `obs_prob()` group 2 (`ox_inf`) are already in the .stan
(inert at Tier 1 because no `tier2_active` row is in the data — `needs_group2` gate
in `tier_specs.R`). Building `t2-indiv-vax` = `t1-indiv-vax`'s tier_col/individualize
/vax switches, but reading `tier2_active` rows too:

- Restored: `W-I-3` (1e3), `W-I-4` (1e4), `J-I-ctrl`/`J-I-ViTT`/`J-I-ViPS` (2e4) — 5 rows.
- **Not restored**: `D-I-plac` (Darton grouped shedding) — same reasoning as
  `t2-indiv`: it's a double-count of the 30 `ox_inf_indiv` rows for the same 30
  volunteers (already in `t1-indiv-vax` via `+cascade`). `data_prep.R:134` already
  documents this: "Under ψ ... the stool-vs-broad contrast is carried by a decoupled
  sub-likelihood on the S1 cross-tab, so neither grouped row is needed."
- N_obs: 177 + 5 = **182**.
- No new Stan code for this piece — it's a registry entry (`TIER_SPECS` in
  `tier_specs.R`) plus confirming `build_tier_data("t2-indiv-vax", ...)` produces
  the right `expect` counts. Genuinely mechanical.

## 2. ψ-correction — new Stan code, and one gating design decision

ψ multiplies `P_inf` for `md_inf` (group 4) rows by a definition-sensitivity factor,
relative to the broadest definition (`bact_or_stool`/blood-or-stool, ψ=1):

| row(s) | definition | ψ |
|---|---|---|
| `H-I-7` (Hornick) | stool-or-blood culture | 1 (reference, no correction) |
| `Lev-I-*` (Levine) | stool-only, any time | `psi_stool` |
| `Gil-I-ctrl` (Gilman) | stool-only, late (4–30d) | `psi_late` |

Plus a **decoupled sub-likelihood** exactly like the existing Darton φ0 ladder
(`typhoid_dose_response.stan:223-229, 346-349, 389-391`) — same pattern, new term:

```
psi_crosstab_y = 15   # stool-positive AMONG Darton placebo bact_or_stool+ (26) -- see correction below
psi_crosstab_n = 26
psi_crosstab_dose = 18200

target += binomial_lpmf(psi_crosstab_y | psi_crosstab_n,
                         psi_stool * eta_detection(psi_crosstab_dose, N50_inf, eta_lo, kappa));
```

This is the §2.8-documented `psi_stool · eta(18200)` confound made explicit in code
(not hidden) — the sub-likelihood pins the *product*, and only η's dose-dependence
(from the other Oxford doses, item 1 above) separates the two factors.

**Correction to §2.8's stated numbers [Mike, caught mid-implementation 2026-07-31].**
`joint_inference_plan.md` §2.8 states "19 ~ binomial(26, psi_stool)" using the
**marginal totals** (19 total stool_positive, 26 total bact_or_stool+ among the 30
Darton Placebo subjects). The subject-level cross-tab
(`analysis_data/darton_individual_endpoints.csv`) shows these **do not nest**: only
**15** of the 26 bact_or_stool+ subjects are also stool_positive — 4 subjects have
`stool_positive=1` with `bact_or_stool=0` *and* `bacteremia=0`, so `bact_or_stool`
is not simply `bacteremia | stool_positive` in this extract (source of the
discrepancy not yet traced to `darton_s1_extract.R`). **Decision: use the true
nested intersection, 15/26 (anchor ≈0.58), not the doc's marginal 19/26 (≈0.73).**
`priors.yaml`'s `psi_stool` prior updated accordingly (weakly-informative, not
tightly anchored on this same count — would double-count the crosstab likelihood
term's own information). §2.8 needs its own correction note (doc drift, item 4
below) — **flagging as an open task, do not silently treat 19/26 as still valid
anywhere else in the repo.**

### Design decision A — how ψ_late stays ≤ ψ_stool

§2.8: "a narrower window cannot detect more" — Gilman late-shedding sensitivity
should be bounded above by Levine's any-time-stool sensitivity, with no cross-tab
of its own to pin it (fully prior-carried per §2.8's caveat 2).

**Proposed:** parameterize `psi_late = psi_stool * frac_late`, `frac_late ~
Beta(a,b)` (its own weakly-informative prior, e.g. mean ~0.8), which makes the
ordering **structural** rather than a soft constraint to check post-hoc. Reported
parameter becomes `frac_late` instead of a free `psi_late`; `tier_report_pars()`
needs the substitution. Alternative: two independent `psi_stool`/`psi_late` free
params with an `ordered` or inline `<= ` constraint and let the data + prior sort it
out — simpler Stan, but the ordering is then only enforced by the constraint, not
implied by the mechanism. **Flagging for you to pick** — no strong reason for one
over the other from the docs; I lean toward the multiplicative form because it
matches "narrower cannot exceed broader" as a structural fact rather than an
estimated-then-checked one, but it's a modeling choice, not a derivation.

### Design decision B — how ψ avoids leaking into Tier 1 fits

Unlike η (gated for free: `ox_inf` rows only exist in `tier2_active` data) and vax
(gated for free: `vaccine_id != 0` rows only exist in vax tiers), **ψ's target rows
(`md_inf`: Hornick/Levine/Gilman) are already present in every active tier**,
including `t1-indiv-vax`. If `psi_stool`/`psi_late` are simply declared as Stan
parameters and group 4 is unconditionally multiplied by ψ, **every existing t1-\*
fit changes** the moment the new .stan is used — silently, because
`STAGE_INCREMENTS`' gates (`needs_group2`, `needs_vaccine_rows`) are both
row-presence gates, and there is no row whose presence signals "ψ should apply
here."

**Proposed:** a genuinely new gate type — a Stan `data` scalar
`psi_active` (0/1), threaded through like `include_vaccine_arms`, defaulting to 0
for every `t1-*` tier (ψ inert, group 4 = plain `maryland_mixture()`, bit-identical
to today) and 1 for `t2-*` tiers. `derive_stage_token()` in `tier_specs.R` gets a
third gate kind (`needs_psi_flag`, driven by the data list's `psi_active` rather
than a row-presence test) alongside `needs_group2`/`needs_vaccine_rows`. This
preserves the registry's invariant (a stage is derived + asserted, not just
labelled) without requiring a signaling row that doesn't exist. **Flagging this
because it's a small extension to `tier_specs.R`'s design, worth you seeing before
I touch that file** — the alternative (add a dummy `psi_def` per-row covariate that
tiers before t2 set to a value that maps to identity) is more convoluted for the
same result.

## 3. Registry changes (`tier_specs.R`)

- `STAGE_INCREMENTS`: add `psi` token, `pars = c("psi_stool", "frac_late")` (or
  `psi_late`, per decision A), gated per decision B.
- New tier `t2-indiv-vax`: `tier_col = "tier2_active"`, `individualize_darton =
  TRUE`, `include_vaccine_arms = TRUE`, `psi_active = TRUE`, stage
  `phi-rho-eta-psi-vax`, `expect$N_obs = 182L` + group counts (177's groups + `ox_inf
  = 5`).
- `t2-grouped`/`t2-indiv` (non-vax): once η+ψ are implemented generically (not
  vax-specific), these become mechanically fittable too. Proposal: **unblock all
  three t2-\* tiers together**, rewriting their `blocked_reason` → `note` as
  "Option A + ψ adopted despite documented confounds (η cross-study; ψ_late
  prior-carried) — same treatment as CoP_imm," rather than deleting the history.
  Confirm you want all three live, not just the vax one.

## 4. Doc drift to fix (per your go-ahead)

- `joint_inference_plan.md` §2.6: add a forward-reference noting the "exclude
  Oxford shedding" recommendation is operationalized, not overridden, by §2.8's η
  once Tier 2 restores the rows — the two sections currently read as contradicting
  each other.
- `tier1.5_plan.md`: currently silent on Tier 2 entirely (its "Later increments"
  section stops at `+vaccine-terms`); add a short forward-pointer to this file.
- `TIER_LADDER.md`: regenerate (`write_tier_ladder_md()`) once the registry changes
  land — it's generated, not hand-edited.
- `CALIBRATION_WORKFLOW.md:17`: the `t2-grouped`/`t2-indiv` row ("declared and
  BLOCKED") needs updating once/if unblocked per §3 above.

## 5. Verification plan (mirrors `+vaccine-terms`'s own verification)

- Parity gate (figure math `model_math.R` vs Stan `p_pred`) extended to groups
  4 (ψ) and the new ψ-crosstab log_lik tail entry, both prior and posterior.
- `validate_all_tiers()` picks up `t2-indiv-vax` with `stanc=TRUE`, asserting
  N_obs=182 and the group table.
- 4-chain sample, target 0 divergences, R-hat ≤ 1.01.
- priorsense on `psi_stool`, `frac_late`/`psi_late`, `eta_lo`, `kappa` — expect
  `frac_late` (or `psi_late`) to show as prior-dominated per §2.8's own prediction;
  flag if `psi_stool` or `eta_lo`/`kappa` show prior-data-conflict (would sharpen
  the confound discussion above from "accepted" to "needs a second look").
- Same shared-biology-shift check as `+vaccine-terms`: fitted medians on the 10
  non-Darton reference rows, before vs after — flag if movement exceeds the
  ~0.03 (absolute) precedent set by that increment.
- Implied-ψ corroboration table (already in §2.8, computed from t1 residuals:
  H-I-7 ≈1.06, Levine ≈0.87, Gilman ≈0.89) as a prior-predictive sanity check
  against the fitted posterior.

## Open questions for you before I touch `typhoid_dose_response.stan`

1. Decision A: multiplicative `psi_late = psi_stool * frac_late`, or independent
   params with an ordering constraint?
2. Decision B: the `psi_active` data-flag gate — good, or do you want a different
   mechanism?
3. Unblock `t2-grouped`/`t2-indiv` (non-vax) alongside `t2-indiv-vax`, or leave
   those two as-is and only build the vax config?
