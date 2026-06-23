# Dose-Response Extraction Progress Checklist

**Contract document**: `dose_response_extraction_contract.md`
**Last updated**: 2026-06-23

---

## Current status (2026-06-23)

Phases 1–3 complete. Phases 4–7 were carried out in **consolidated form** inside
`calibration/joint_inference_plan.md` (+ `reviewer2_response.md`, the data
`notes/verification/` set, and `notes/meeting_famulare_buffalo_dose_response_identifiability*.md`)
rather than as the separately-named files the original checklist anticipated — so
some boxes below are unchecked even though the work exists; pointers are annotated
inline.

Phase 8 (fit) is underway and currently **blocked on two root-caused defects**:

1. **Sampling pathology** — the model compiles and samples but produces ~99%
   divergences; root cause isolated to the ordering-constraint *density cliff* at
   `typhoid_dose_response.stan` L132 (removing it → 0/4000). Fix = `<lower=0>`
   offset reparameterization.
2. **Constant-φ misspecification** — `φ=0.25` is the *low-dose asymptote* of the
   definition-sensitivity function applied as a dose-constant; it caps Maryland
   fever at 0.25 while Hornick's high-dose data are 0.89–0.95 (structurally
   unfittable, not just biased). The plan's own default (and Reviewer 2) call for
   *dose-dependent* φ; the constant is the "initial exploration" shortcut promoted
   into the implementation.

Both are documented in `calibration/tier1_pathology_diagnosis.md`. A
workflow-improvement step is in progress to make this kind of debugging
repeatable. Next: reparameterize the constraint, implement dose-dependent φ, then
re-fit Step 1 of `calibration/CALIBRATION_WORKFLOW.md`.

---

## Phase 1: Setup & Infrastructure

- [x] Create folder structure
  - [x] `extracts/`
  - [x] `schemas/`
  - [x] `analysis_data/`
  - [x] `notes/`
  - [x] `calibration/`
- [x] Create `notes/outcome_mapping.md` (initial draft)
- [x] Create `notes/paper_triage.md` (empty template)
- [x] Create `notes/cross_cutting_observations.md` (empty template)

---

## Phase 2: Batch First-Pass Extraction

### Paper Extractions (20 papers)

| # | Paper | Extract Created | Triage Updated | Status |
|---|-------|-----------------|----------------|--------|
| 1 | Dahora et al. 2019 | [x] | [x] | complete |
| 2 | Darton et al. 2012 | [x] | [x] | complete |
| 3 | Darton et al. 2016 | [x] | [x] | complete |
| 4 | Darton et al. 2017 | [x] | [x] | complete |
| 5 | Dupont et al. 1971 | [x] | [x] | complete |
| 6 | Gibani et al. 2019 | [x] | [x] | complete |
| 7 | Gibani et al. 2020 | [x] | [x] | complete |
| 8 | Gilman et al. 1977 | [x] | [x] | complete |
| 9 | Glynn & Bradley 1992 | [x] | [x] | complete |
| 10 | Glynn et al. 1995 | [x] | [x] | complete |
| 11 | Hornick et al. 1967 | [x] | [x] | complete |
| 12 | Hornick (Appraisal) | [x] | [x] | complete |
| 13 | Hornick 1966 | [x] | [x] | complete |
| 14 | Hornick & Snyder 1970 | [x] | [x] | complete |
| 15 | Jin et al. 2017 | [x] | [x] | complete |
| 16 | Juel et al. 2018 | [x] | [x] | complete |
| 17 | Levine et al. 1976 | [x] | [x] | complete |
| 18 | Waddington et al. 2014 (review) | [x] | [x] | complete |
| 19 | Waddington et al. 2014 (outpatient) | [x] | [x] | complete |
| 20 | Woodward 1980 | [x] | [x] | complete |

### Cross-Cutting Documents

- [x] `notes/cross_cutting_observations.md` needs completion in Phase 3
- [x] `notes/paper_triage.md` completed with all papers
- [ ] `notes/outcome_mapping.md` to be updated with all observed outcome types in Phase 3

---

## Phase 3: Reference Model Specification

**Purpose**: Define the most elaborate model the literature supports—the benchmark for reasoning about simplifications.

### Deliverables

- [x] `notes/dose_response_model_specification.md` written with full prose + equations
- [x] DAG created showing latent processes → observables
- [ ] Explicit list of simplifications needed to reach practical model (Section 6.5, after Phase 4)
- [ ] Mapping from reference model components to available data (Section 7, after Phase 4)

### Components to Specify

- [x] Latent biological processes (infection cascade, outcomes hierarchy)
- [x] Severity gradations
- [x] Joint/conditional relationships between outcomes
- [x] Immunity representation (full mechanistic → CoP simplification pathway)
- [x] Dose-response framework
- [x] Observational model (P(obs | latent) for each measurement type)
- [x] Heterogeneity sources (study-level and individual-level)

---

## Phase 4: Paper-by-Paper Joint Review

### Review Sessions

| Paper | Reviewed | [ASSISTANT-PROPOSED] Resolved | [OPEN] Resolved | Finalized |
|-------|----------|------------------------------|-----------------|-----------|
| (to be filled during Phase 4 review) | | | | |

### Review Order (based on triage):
1. **Core multi-dose papers**: Hornick 1966, Hornick & Snyder 1970, Waddington 2014 (outpatient)
2. **Support papers with immunity**: Darton 2016, Jin 2017, Gibani 2020, Gilman 1977, Dupont 1971, Levine 1976
3. **Remaining papers**: Dahora 2019, Darton 2017, Gibani 2019, Juel 2018, Glynn 1995, Glynn & Bradley 1992, Woodward 1980
4. **Exclusions to confirm**: Darton 2012, Waddington 2014 (review)

---

## Phase 5: Normalization & Schema Design

> **Done in consolidated form**, not as separate files: the working model lives in
> `calibration/joint_inference_plan.md` (§2) and `notes/yolo_working_model_notes.md`;
> the YAML-schema step was skipped in favor of building the analysis CSV directly.

- [~] Reference model collapsed to working model (in `joint_inference_plan.md` §2, not `notes/working_model.md`)
  - [x] Estimable vs simplified components documented (`joint_inference_plan.md` §2, Table at §11)
  - [x] Simplifications justified (data limitation vs parsimony)
  - [x] Working model specification written (consolidated into `joint_inference_plan.md`)
  - [ ] `notes/dose_response_model_specification.md` Sections 7-8 updated
- [ ] YAML schema designed — **skipped** (CSV built directly)
- [ ] All extracts converted to YAML — **skipped**
- [x] Analysis-ready CSV compiled (`calibration/dose_response_data.csv`, 37 obs)
- [~] Data dictionary created (column semantics in CSV `notes` field + `joint_inference_plan.md` §4)

---

## Phase 6: Calibration Problem Design

> **Done in consolidated form** inside `calibration/joint_inference_plan.md` and the
> identifiability meeting notes.

- [~] `notes/outcome_mapping.md` finalized (definition mapping in `joint_inference_plan.md` §2.5)
- [x] Likelihood design written (`joint_inference_plan.md` §5–6, not `calibration/likelihood_design.md`)
- [x] Identifiability memo (`notes/meeting_famulare_buffalo_dose_response_identifiability*.md` + `joint_inference_plan.md` §2.6, MC3)
- [x] Heterogeneity structure decided (study RE + Maryland mixture; `joint_inference_plan.md` §2.3–2.4)
- [x] Latent immunity model specified (Maryland mixture / latent CoP; `joint_inference_plan.md` §2.4)

---

## Phase 7: Prior Specification

> Priors now in `calibration/priors.yaml` as the **single source of truth**; the
> Stan model reads hyperparameters from its data block (change + refit, no recompile).
> Values from `joint_inference_plan.md` §7, still flagged **"TO BE EXAMINED"** but now
> screened by a prior predictive run + priorsense power-scaling sensitivity.

- [x] `calibration/priors.yaml` created (consumed by `priors.R` → Stan data + overlay/sampler)
- [~] All prior choices documented with rationale (`priors.yaml` comments + `joint_inference_plan.md` §7; priorsense flags δ + alpha_fevginf sensitivity)

---

## Phase 8: Fit, Validate, Document

Workflow ladder + status: see `calibration/CALIBRATION_WORKFLOW.md`.

- [x] Model compiles & samples (cmdstanr + CmdStan 2.39); driver `calibration/fit_dose_response.R`
- [x] MCMC diagnostics passed — **0/4000 divergent, R-hat<1.01, ESS>1700**. Fixed via prior-preserving `<lower=0>` offset reparam (`d_fev`); controlled attribution in `simulate_recovery.R` confirms cliff 99.0% → reparam 0.0%
- [x] Likelihood refactor verified — unified `obs_prob()` matches the original per-group likelihood (`test_obs_prob_parity.R`, the guard recovery cannot provide)
- [ ] Implement dose-dependent φ — constant φ=0.25 caps Hornick high-dose fever below the data; now **visible in the PPC** (H-F-8/H-F-9). Open (`tier1_pathology_diagnosis.md` §5a; Reviewer 2 MC1; plan §2.5)
- [x] Posterior predictive checks completed — rate-space PPC fed by Stan `p_pred` (no R mirror); `diagnostics.R`
- [~] Sensitivity analyses — harness in place (`run_scenarios.R`: row-filter + prior-override scenarios, loo on grouped units); structural-Stan variants (share α/γ, single CoP_md, drop φ) deferred
- [ ] Final documentation written
- [x] **Workflow upgrade (Buffalo-style, R-native)**: `diagnostics.R`, `simulate_recovery.R`, `run_scenarios.R`, `priors.yaml`/`priors.R`, `data_prep.R` (plan `please-plan-a-b-smooth-ritchie.md`)

---

## Session Log

| Date | Session Summary | Stopping Point | Next Steps |
|------|-----------------|----------------|------------|
| 2026-02-03 | Created contract and checklist | Ready to begin Phase 1 | Create folder structure, begin batch extraction |
| 2026-02-03 | Completed Phase 1 setup | Phase 1 complete | Begin Phase 2 batch extraction |
| 2026-02-03 | Completed Phase 2 batch extraction (all 20 papers) | Phase 2 complete | Begin Phase 3 reference model |
| 2026-02-03 | Fixed mislabeled Hornick 1970 Part 1 PDF filename; merged duplicate extracts | Cleanup complete | Phase 3 reference model |
| 2026-02-03 | Inserted Phase 3 (Reference Model Specification); renumbered phases 3-7 → 4-8 | Contract/checklist updated | Begin Phase 3 reference model |
| 2026-02-03 | Added Phase 5.1 (collapse reference model to working model); created `notes/reference_model.md` draft | Contract/checklist updated | Begin Phase 3 reference model |
| 2026-02-04 | Completed Phase 3 model specification; renamed to `dose_response_model_specification.md`; cleaned contract/spec separation; aligned `outcome_mapping.md`; added README | Phase 3 complete | Phase 4 joint review |
| 2026-02–05 (bridge, from git history) | Phases 4–7 consolidated into `calibration/joint_inference_plan.md`: data verification (`notes/verification/`), Reviewer 2 response, Darton S1 individual-level extraction, removal of fabricated Waddington 10⁵ arm, Oxford-shedding exclusion, identifiability meeting; canonical `dose_response_data.csv` (37 obs) committed; untested Tier-2 Stan skeleton committed | Plan + data + skeleton in place | Run the Stan inference |
| 2026-06-23 | **Resurrection + pathology diagnosis.** Installed cmdstanr/CmdStan 2.39; fixed 3 defects so the model compiles & samples (`_lp` rename, Hornick double-count, `prior_only` flag); wrote driver `fit_dose_response.R`. First Tier-1 fit: ~99% divergent. Root-caused to the ordering-constraint density cliff (L132) by elimination (removing it → 0/4000). Separately root-caused the constant-φ misspecification (φ=0.25 is the low-dose asymptote applied as a dose-constant; caps Hornick high-dose fever below the data). Docs: `tier1_pathology_diagnosis.md`, `CALIBRATION_WORKFLOW.md`. Branch `dose-response-tier1-resurrection`. | Phase 8 blocked on 2 root-caused defects; workflow-improvement step in progress | Reparameterize constraint (`<lower=0>` offset); implement dose-dependent φ; re-fit Step 1 |
| 2026-06-23 | **Buffalo-style workflow upgrade + divergence fix.** Adopted Vince Buffalo's Stan-workflow patterns R-native (plan `please-plan-a-b-smooth-ritchie.md`). Stan refactor: unified `obs_prob()` (one source for all 5 likelihood groups), flat per-obs data, `lprior` accumulator, per-obs `p_pred`/`y_rep`/`log_lik` in GQ, prior-preserving N50 reparam (`d_fev<lower=0>`). New R-native infra: `diagnostics.R` (bayesplot/posterior battery + PPC from Stan p_pred + priorsense), `simulate_recovery.R` (known-truth recovery + SBC-lite + cliff-vs-reparam attribution), `run_scenarios.R` (cross-run comparison + loo on grouped units), `priors.yaml`/`priors.R` (single-source priors as Stan data), `data_prep.R`. Parity-gated (`test_obs_prob_parity.R`). **Result: Tier-1 now 0/4000 divergences (was 99%); attribution confirms it's the N50 geometry alone.** Independent codex review incorporated. | Step 1 unblocked & clean; φ-cap + δ-tension remain as science decisions | Implement dose-dependent φ (MC1); wire study RE (Step 2); decide η option (Step 3); run full SBC (k≥40) |

