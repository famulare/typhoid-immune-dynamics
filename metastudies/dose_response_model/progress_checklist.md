# Dose-Response Extraction Progress Checklist

**Contract document**: `dose_response_extraction_contract.md`
**Last updated**: 2026-02-03

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

### Paper Extractions (21 papers)

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
| 14 | Hornick & Snyder 2008 | [x] | [x] | complete (duplicate of 1970 Part 1) |
| 15 | Hornick & Snyder 1970 | [x] | [x] | complete |
| 16 | Jin et al. 2017 | [x] | [x] | complete |
| 17 | Juel et al. 2018 | [x] | [x] | complete |
| 18 | Levine et al. 1976 | [x] | [x] | complete |
| 19 | Waddington et al. 2014 (review) | [x] | [x] | complete |
| 20 | Waddington et al. 2014 (outpatient) | [x] | [x] | complete |
| 21 | Woodward 1980 | [x] | [x] | complete |

### Cross-Cutting Documents

- [x] `notes/cross_cutting_observations.md` needs completion in Phase 3
- [x] `notes/paper_triage.md` completed with all papers
- [ ] `notes/outcome_mapping.md` to be updated with all observed outcome types in Phase 3

---

## Phase 3: Paper-by-Paper Joint Review

### Review Sessions

| Paper | Reviewed | [ASSISTANT-PROPOSED] Resolved | [OPEN] Resolved | Finalized |
|-------|----------|------------------------------|-----------------|-----------|
| (to be filled during Phase 3 review) | | | | |

### Review Order (based on triage):
1. **Core multi-dose papers**: Hornick 1966, Waddington 2014 (outpatient)
2. **Support papers with immunity**: Darton 2016, Jin 2017, Gibani 2020, Gilman 1977, Dupont 1971, Levine 1976
3. **Remaining papers**: Dahora 2019, Darton 2017, Gibani 2019, Juel 2018, Glynn 1995, Glynn & Bradley 1992, Woodward 1980
4. **Exclusions to confirm**: Darton 2012, Hornick 2008 (duplicate), Waddington 2014 (review)

---

## Phase 4: Normalization & Schema Design

- [ ] YAML schema designed
- [ ] All extracts converted to YAML
- [ ] Analysis-ready CSV compiled
- [ ] Data dictionary created

---

## Phase 5: Calibration Problem Design

- [ ] `notes/outcome_mapping.md` finalized
- [ ] `calibration/likelihood_design.md` written
- [ ] `notes/identifiability_memo.md` written
- [ ] Heterogeneity structure decided
- [ ] Latent immunity model specified

---

## Phase 6: Prior Specification

- [ ] `calibration/priors.yaml` created
- [ ] All prior choices documented with rationale

---

## Phase 7: Fit, Validate, Document

- [ ] Model implemented in Stan/brms
- [ ] MCMC diagnostics passed
- [ ] Posterior predictive checks completed
- [ ] Sensitivity analyses completed
- [ ] Final documentation written

---

## Session Log

| Date | Session Summary | Stopping Point | Next Steps |
|------|-----------------|----------------|------------|
| 2026-02-03 | Created contract and checklist | Ready to begin Phase 1 | Create folder structure, begin batch extraction |
| 2026-02-03 | Completed Phase 1 setup | Phase 1 complete | Begin Phase 2 batch extraction |
| 2026-02-03 | Completed Phase 2 batch extraction (all 21 papers) | Phase 2 complete | Begin Phase 3 joint review |

