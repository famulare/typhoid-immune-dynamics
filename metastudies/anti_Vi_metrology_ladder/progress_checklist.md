# Anti-Vi Metrology Bridge Progress Checklist

## Phase 1: Setup & Infrastructure
- [x] Folder structure verified/created
- [x] Extraction template defined (`extracts/_TEMPLATE.md`)
- [x] `assay_mapping.md` initialized
- [x] Decision bookkeeping convention documented (in contract)

## Phase 2: Batch First-Pass Extraction
- [ ] Barrett 1983 (ELISA development) extracted
- [ ] Barrett 1983 (carrier identification) extracted
- [ ] Lee 2020 (PLOS NTD) extracted
- [ ] NIBSC 16/138 IFU extracted
- [ ] Rijpkema 2018 (Vaccine) extracted (if accessible)
- [ ] WHO BS.2017.2307 extracted (if useful)
- [x] `cross_cutting_observations.md` created (placeholder structure)
- [x] `paper_triage.md` created (placeholder structure, to be filled during extraction)

## Phase 3: Reference Bridge Specification
- [x] `metrology_bridge_specification.md` drafted (skeleton with Sections 5-7 pending data)
- [x] Bridge graph created (ASCII in specification)
- [x] Latent immunity construct defined
- [ ] Simplification principles documented (after Phase 4)

## Phase 4: Paper-by-Paper Joint Review
- [ ] Barrett 1983 reviewed with user
- [ ] Lee 2020 reviewed with user
- [ ] NIBSC/WHO documentation reviewed
- [ ] All [ASSISTANT-PROPOSED] resolved
- [ ] All [OPEN] items dispositioned
- [ ] Extracts finalized

## Phase 5: Normalization & Schema Design
- [ ] Working model defined
- [ ] YAML schema designed
- [ ] Extracts converted to YAML
- [ ] Paired datasets compiled

## Phase 6: Bridge Model Design
- [ ] Assay mapping finalized
- [ ] Likelihood structure documented
- [ ] Heterogeneity structure specified
- [ ] Identifiability memo completed

## Phase 7: Prior Specification
- [ ] Priors specified with rationale
- [ ] Prior predictive checks run

## Phase 8: Fit, Validate, Document
- [ ] Models implemented
- [ ] Convergence diagnostics passed
- [ ] Posterior predictive checks
- [ ] Sensitivity analyses
- [ ] Stop conditions evaluated
- [ ] Final documentation complete
