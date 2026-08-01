# Dose-response project progress checklist

**Contract document**: `dose_response_extraction_contract.md`
**Last updated**: 2026-08-01

## Current status: complete for the locked working deliverable

The dose-response metastudy has reached a reproducible project-completion
checkpoint. The default downstream artifact is the **Tier 2 individualized-vaccine
fit**, `t2-indiv-vax`, stage `phi-rho-eta-psi-vax`. All planned extraction,
specification, implementation, diagnostic, and handoff work needed for this
working result is complete.

This completion statement is deliberately bounded. It means the repository has a
canonical model, data, fit, diagnostics, figures, and documented assumptions. It
does not mean that the model is fully identified, that every residual is explained,
or that the remaining sensitivity analyses are unnecessary.

## Locked Tier 2 result

- [x] Tier 2 individualized vaccine is the default and canonical working fit.
- [x] Oxford shedding is restored with the `eta(D)` treatment-truncation correction.
- [x] Infection-definition sensitivity `psi_d` is active; the Darton nested
  cross-tab uses 15/26, not the marginal 19/26.
- [x] The shared beta-binomial overdispersion parameter
  `grand_overdispersion_rho` is implemented; the study/cohort random effect was
  withdrawn and `sigma_study` deleted.
- [x] Darton placebo, M01ZH09, and Ty21a arms are individualized where supported.
- [x] The preferred fit presents 182 Stan observations: 29 grouped and 153
  individual.
- [x] Prior-predictive and posterior-predictive artifacts are retained.

### Retained fit checks

| Check | Retained result |
|---|---:|
| Posterior draws | 4 chains x 1,000 warmup + 1,000 sampling |
| Divergences | 1 / 4,000 (0.025%) |
| Max-treedepth hits | 0 |
| Minimum E-BFMI | 0.947 |
| Maximum reported R-hat | 1.003 |
| Figure/Stan curve parity | max absolute difference 3.6e-08 |

The one posterior divergence is an explicit residual review item. The result is
locked as the working baseline, not described as diagnostically perfect.

## Phase completion

- [x] **Phase 1 - Setup and infrastructure.** Folder structure, extraction
  conventions, outcome notes, and calibration scaffolding are in place.
- [x] **Phase 2 - First-pass extraction.** All 20 listed papers have extracts and
  triage entries; source provenance is retained in `input_papers/` and `extracts/`.
- [x] **Phase 3 - Reference model specification.** Biology, latent states, DAG,
  equations, and reference-to-practical simplifications are documented in
  `dose_response_model_specification.md`.
- [x] **Phase 4 - Joint review.** Review, verification, and decision work was
  consolidated into `joint_inference_plan.md` and the linked meeting/verification
  notes rather than kept as a separate per-paper table.
- [x] **Phase 5 - Normalization and schema design.** The canonical analysis CSV
  and data dictionary are complete. YAML schemas and per-extract YAML conversion
  were intentionally skipped in favor of the directly validated CSV.
- [x] **Phase 6 - Calibration problem design.** Likelihood, outcome-definition
  maps, latent Maryland immunity, no-double-counting rules, and identifiability
  assumptions are documented and locked.
- [x] **Phase 7 - Priors.** `calibration/priors.yaml` is the single source of
  prior hyperparameters, consumed as Stan data and covered by prior-predictive
  diagnostics. Prior-carried parameters remain explicitly labeled.
- [x] **Phase 8 - Fit, validate, and document.** The Stan model, tier registry,
  fit driver, posterior/prior runs, diagnostic battery, PPCs, parity check, figure
  suite, and handoff documents are retained.

## Residual limitations (not hidden by completion)

- `N50` and the medium bridge `delta` trade off (r = -0.61 / -0.69 on the log10
  scale at `t2-indiv-vax`) but are NOT jointly unidentified: each marginal
  contracts about 2x from its prior, the milk-frame product `N50_inf*delta`
  contracts 3.3x, and the split between them still contracts 1.7x. The caution is
  interpretive -- `delta` absorbs every era difference the vehicle stands in for --
  not a claim that only the product is estimable. (Corrected 2026-08-01; the
  earlier "structurally confounded" wording, and the hardcoded -0.72/-0.78 in the
  `delta_bridge` caption, overstated it and matched no fit in the repo.)
- The Maryland immunity mixture is latent and partly prior-carried.
- `eta(D)`, `psi_stool`, and the late-shedding fraction are only partly separated
  by the available Oxford data. `kappa` IS prior-dominated as expected (priorsense
  prior 0.63 / likelihood 0.05); `frac_late` is NOT (prior 0.08 / likelihood 0.20)
  -- Gilman's `psi_def=2` row constrains it more than `joint_inference_plan.md`
  Sec 2.8 assumed. `eta_lo` meets `tier2_plan.md` Sec 5's stated trigger for a
  second look at the eta confound, though priorsense flags 16 of 26 parameters at
  similar magnitudes, so the signal is not specific. OPEN.
- `gamma_inf` falls monotonically as individualized Darton rows enter the
  likelihood (t1-indiv 0.166, t1-indiv-vax 0.111, t2-indiv-vax 0.081) while
  `gamma_fevginf` rises. The fitted titre slope is shallower than every direct
  contrast in the corpus, including both Jin Vi-vaccine arms. Individualization
  puts 153 thin-titre-range rows against 3 arm-level rows that carry the high-titre
  information. A down-weighted / grouped-Darton sensitivity fit is the check. OPEN.
- A single posterior divergence remains in the preferred fit.
- The Gibani rechallenge/susceptibility paradox, incubation, seroconversion, and
  other unused outcomes remain outside the current likelihood.

These are scientific interpretation and future-sensitivity boundaries, not missing
implementation steps for the locked working result.

## Canonical handoff paths

1. [README](README.md) - folder map and fresh-clone rebuild ladder.
2. [Onboarding one-pager](onboarding_one_pager.md) - concise current orientation.
3. [Model specification](dose_response_model_specification.md) - biology and data
   mapping.
4. [Joint inference plan](joint_inference_plan.md) - locked likelihood and
   identifiability assumptions.
5. [Tier ladder](calibration/TIER_LADDER.md) - generated tier counts and stages.
6. [Posterior summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md)
   - canonical parameters, diagnostics, and figures.
7. [Run manifest](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/run_manifest.json)
   - input hashes, sampler settings, code/environment provenance, and output list.

## Provenance note

The retained run manifest records the fit generation time, input hashes, sampler
settings, and code SHA. The artifacts were generated from the working tree before
the result-retention and default-selection commits; the manifest records that
working-tree state explicitly. The result is therefore a committed, auditable
working artifact, while byte-identical regeneration should use the manifest rather
than assume the current `HEAD` alone reproduces every output byte.

## Session log

| Date | Milestone | Outcome |
|---|---|---|
| 2026-02-03 to 2026-02-05 | Phases 1-7 consolidated | Source extraction, verification, reference model, joint plan, and canonical CSV established. |
| 2026-06-23 | Tier 1 resurrection and pathology diagnosis | N50 ordering cliff reparameterized; constant-phi misspecification diagnosed; workflow tooling added. |
| 2026-07-31 | Tier 1 overdispersion and Tier 2 implementation | `grand_overdispersion_rho`, Oxford `eta`, infection `psi`, vaccine terms, diagnostics, and canonical figures completed. |
| 2026-08-01 | Default and documentation lock | `t2-indiv-vax` made the discoverable default; result, onboarding, and project-completion documents refreshed. |
