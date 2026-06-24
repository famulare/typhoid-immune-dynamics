# Handoff: sync the calibration harnesses with the model parameter block

**From:** Claude Code session, 2026-06-23 (the rhat.png-bug thread)
**To:** whoever is doing the Tier 1.5 model work
**Why you're reading this:** the model gained an estimated parameter (`phi_md`) but
two harness scripts still assume the old parameter set. I fixed the mechanical and
robustness parts; the remaining items need a *modeling decision* (a `phi_md` truth
value, and re-syncing the parity R reference to the model) — those are yours.

## The drift, in one line

The model `parameters{}` block now has **14** params; the two
`generate_quantities(fitted_params=)` harnesses hardcoded a **13**-param list
missing `phi_md`, so they died with the cryptic
`"Mismatch between model and fitted_parameters csv file"`.

- Current model params (source of truth — `names(mod$variables()$parameters)`):
  `log10_N50_inf, d_fev, alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
  log10_delta, pi_susc, CoP_imm, CoP_susc, phi_md, eta_lo, kappa, sigma_study`
- `phi_md` was promoted from a fixed per-obs data value (0.25/0.65) to an estimated
  scalar parameter — `typhoid_dose_response.stan` line 169; comment at line 129:
  "phi is no longer data." [observed]

## ALREADY DONE — please don't redo

1. **rhat figure bug** — `diagnostics.R` (`plot_battery`) previously fed a
   multi-var `draws_array` to `posterior::rhat()`, which collapses to one spurious
   scalar (the old `rhat.png` showed a single dark-blue R̂≈1.36 "non-convergence"
   point that was an artifact). Now uses `summarise_draws(draws, "rhat")` →
   per-parameter named vector. The plot derives its variable set from the draws, so
   it's robust to the parameter block changing. All 8 `rhat.png` regenerated.

2. **New fail-loud guard** `assert_fitted_params_match(mod, provided)` in
   `data_prep.R` (sourced by both harnesses). It diffs your supplied names against
   `mod$variables()$parameters` and `stop()`s with the missing/extra names *before*
   `generate_quantities`. Wired into `simulate_y` (`simulate_recovery.R`) and before
   the gq call in `test_obs_prob_parity.R`.
   **This is your friend:** when you add any new Tier 1.5 parameter, run a harness
   and it tells you exactly which list is stale, e.g.:
   `missing (in model, not provided): phi_md → update the PARAM_NAMES / truth vector`.

3. **SBC made param-set-agnostic** — `recover_repeated` now builds its truth from
   `names(mod$variables()$parameters)` and prior-samples each (every model param has
   a prior in `priors.yaml`, verified). Any new param with a prior flows through
   automatically; no hardcoded list to maintain. Verified it runs (k=2, 0 div).

## YOUR TODOs — these need a modeling call

### A. Point recovery — `simulate_recovery.R`
`recover_once` now fails the guard with `missing: phi_md` (intentional: loud, not
cryptic). To fix:
- Add `"phi_md"` to `PARAM_NAMES` **and** `phi_md = <value>` to `TRUTH_REALISTIC`
  (`make_truth` requires the two to be `setequal`).
- The value is your judgment — tier1 posterior mean was ~0.97; pick a realistic
  recovery truth.
- Note: `PARAM_NAMES` here is now used *only* by the `recover_once`/`make_truth`
  path (`recover_repeated` no longer uses it). Optional cleanup: derive it from the
  model too and delete the constant.

### B. Parity gate (M4, LOAD-BEARING) — `test_obs_prob_parity.R` — currently RED
**Two edits, both required.** Edit 1 alone will not make it pass.

1. Add `"phi_md"` to `PARAM_NAMES` (~line 64) and a value at the **same index** in
   each of the three `vecs` (~lines 68–70). phi_md ∈ (0,1); any valid values
   exercise parity.
2. **Update the R reference `obs_prob_R`** (the independent transcription it checks
   Stan against). It currently reads the old fixed per-row `phi <- row$phi` in the
   `md_fev` (g==3) and `hornick_cond` (g==5) branches. The model now uses the scalar
   parameter `phi_md`. Change `phi <- row$phi` → `phi <- p$phi_md` in both branches
   (and drop the unused `row$phi`). Without this the gate will FAIL on a *real*
   numeric disagreement (or error, if `build_stan_data` no longer emits `row$phi`) —
   not a name error.
3. Verify: `Rscript test_obs_prob_parity.R` → expect `PARITY PASS` (exit 0).

## Heads up

- My edits stack on your uncommitted WIP in `data_prep.R`, `simulate_recovery.R`,
  and `test_obs_prob_parity.R` — you'll see the new `assert_fitted_params_match`
  calls and the model-derived SBC truth.
- `results/` is gitignored, so the regenerated figures won't appear in
  `git status`; the only tracked code change from my side is in these `.R` files.

## Quick reference

| Thing | Where |
|---|---|
| Guard | `data_prep.R` → `assert_fitted_params_match` |
| Point recovery | `simulate_recovery.R` → `PARAM_NAMES`, `TRUTH_REALISTIC`, `make_truth` |
| SBC (already fixed) | `simulate_recovery.R` → `recover_repeated` |
| Parity gate | `test_obs_prob_parity.R` → `PARAM_NAMES`, `vecs`, `obs_prob_R` |
| Model params (truth) | `names(cmdstan_model("typhoid_dose_response.stan", compile=FALSE)$variables()$parameters)` |

---

## UPDATE 2026-06-23 (Tier 1.5 minimal+φ session) — A & B resolved

- **B (parity gate) — DONE, GREEN.** Added `phi_md` to `PARAM_NAMES`/`vecs` and
  changed `obs_prob_R` `phi <- row$phi` → `phi <- p$phi_md` in the g==3/g==5
  branches. `Rscript test_obs_prob_parity.R` → **PARITY PASS** (max |Δp|=5e-9 over
  3 vectors × **54 rows**, incl. the 30 Darton individual rows).
- **A (point recovery) — `phi_md` added** to `PARAM_NAMES`/`TRUTH_REALISTIC`
  (=0.97) so the guard passes. *Caveat:* the other truth values are still on the
  pre-EU/mL scale (see below).

### Tier 1.5 remaining (the model task list — "get through Tier 1.5")

1. **C3 — dose-dependent φ(T,D)** [BANKED tonight; do next session]. The Darton
   threshold ladder gives φ₀(T) ≈ 0.30 (39.4) / 0.66 (38.3) — but a *dose-independent*
   φ(T) re-introduces the high-dose Hornick cap (would fit H-F-9 at ~0.27 vs 0.95).
   Implement the plan §2.5 form `φ(T,D)=φ₀(T)+(1−φ₀(T))·P_fev_naive(D)^β_φ`: φ₀(T) from
   the ladder, → 1 at saturating dose. Adds φ₀-level, λ, β_φ (retires `phi_md`);
   **β_φ rests on Hornick's 4 multi-dose fever points — thin-data ID, watch it.** See
   `tier1_lab_notebook.md` (C3 finding) + `tier1.5_plan.md`.
2. **Harness re-sync after C3** — C3 changes the parameter block again, so BOTH
   `PARAM_NAMES`/`vecs`/`obs_prob_R` (parity) and `PARAM_NAMES`/`TRUTH_REALISTIC`
   (recovery) need the new φ params; `obs_prob_R` must replicate φ(T,D). The guard
   will name what's stale.
3. **Recovery truth re-sync to EU/mL scale** — update `TRUTH_REALISTIC` (delta~2.5,
   CoP_imm~5–7, gamma~0.2, alpha_inf~0.4) so point recovery tests a representative point.
4. **Tier 1.5 plots** — `dose_response_curves.R` dose-axis panels stack the 30 Darton
   individuals at one dose (y∈{0,1}); add a **titre→protection (CoP-axis) panel** —
   anti-Vi EU/mL on x, fever prob on y, individuals + Jin GMT points + the CoP^γ curve —
   the right view for the immunity slope.
5. **Open prior elicitation** (review, not code): `CoP_imm` (Maryland latent immune
   anti-Vi-equiv) is prior-dominated and unidentifiable — highest-leverage unelicited
   choice. See `tier1_lab_notebook.md` "Identifiability & prior-dependence status."
