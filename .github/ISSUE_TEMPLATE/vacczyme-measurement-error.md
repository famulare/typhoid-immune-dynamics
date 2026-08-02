---
name: VaccZyme Measurement Error Investigation
about: Track investigation of anti-IgG assay measurement error modeling
title: "Investigate VaccZyme Anti-IgG Assay Measurement Error Modeling (CV = 0.37)"
labels: ["investigation", "model"]
---

## Summary

Investigate and model the measurement error in the VaccZyme anti-IgG assay used to observe the Correlate of Protection (CoP) in the typhoid dose-response model. The coefficient of variation (CV) is approximately 0.37.

## Current State

Per the dose-response model specification document (§5.7.5), the VaccZyme assay observes anti-Vi IgG where:
$$\text{CoP} = \text{anti-Vi IgG} / 3.7 \text{ EU/mL}$$

where 3.7 EU/mL is the imputed limit of detection (LLD). The assay results currently enter the model without explicit characterization of measurement error, which contributes to model uncertainty, particularly given the left-censored distribution:
- 64 of 90 Darton per-subject rows sit at exactly CoP = 1.0 (the LLD)
- The distribution is left-censored with an atom at the LLD

## Why This Matters

1. **CoP inference**: Assay measurement error directly affects inference of the Correlate of Protection, which conditions all latent biological processes in the model
2. **Substantial variability**: A CV of 0.37 represents ~37% relative measurement variability that should be explicitly propagated through the likelihood
3. **High-titre end**: Critical for the high-titre region of the CoP axis, which is carried by:
   - Three grouped Jin arms (§7)
   - 26 titre-informative Darton subjects
4. **Titre slope identification**: Affects identification of the titre slope parameter (discussed in §8.5)

## Investigation Scope

- [ ] **Review assay documentation**
   - Validate the CV = 0.37 figure from primary sources
   - Document assay replicability and precision
   
- [ ] **Characterize measurement error structure**
   - Determine error model (e.g., log-normal, normal-on-log-scale)
   - Evaluate replicate variance patterns across titre ranges
   - Document censoring mechanism at LLD
   
- [ ] **Determine parameterization**
   - Single shared error term vs. assay-specific parameters
   - Interaction with left-censoring at LLD
   - Potential titre-dependent variance (heteroscedasticity)
   
- [ ] **Impact assessment**
   - Effect on CoP inference (posterior SD, coverage)
   - Effect on dose-response parameter posteriors
   - Effect on LOO model comparison (when evidence is available)
   - Effect on titre slope identification
   
- [ ] **Documentation & implementation**
   - Update `dose_response_model_specification.md` §5.7.5
   - Propose prior specification in `calibration/priors.yaml`
   - Plan implementation for next Tier 2 or Tier 3 update

## Related Documentation

- **Primary**: `metastudies/dose_response_model/dose_response_model_specification.md`
  - §5.7.5 (CoP formalization and VaccZyme assay)
  - §8.5 (Identifiability concerns, titre slope)
  
- **Supporting**:
  - `calibration/priors.yaml` — immunity parameter priors
  - `calibration/jin_within_arm_titre_plan.md` — draft quadrature for within-arm titre variation
  - Data: `analysis_data/darton_individual_endpoints.csv` — Darton per-subject titres

## Expected Outcomes

1. Validated CV value and measurement error model specification
2. Impact quantification on posterior uncertainty
3. Recommended parameterization for model implementation
4. Updated specification documents with error model incorporated
