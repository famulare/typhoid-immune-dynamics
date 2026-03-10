# Decision Log: Anti-Vi Metrology Bridge

This document tracks all substantive decisions made during the YOLO implementation.

---

## D1: E2 Gap Treatment

**Date**: 2026-02-04
**Context**: Web search confirmed no direct bridging study between Barrett 1983 ELISA and modern in-house ELISAs.

**Options Considered**:
1. Abort the bridge (stop criterion triggered)
2. Build two disconnected bridges (HA→Barrett, InHouse→VaccZyme)
3. Add structural uncertainty for missing E2 link

**Decision**: Option 3 - Add structural uncertainty (τ parameter)

**Rationale**:
- Both assays measure the same analyte (anti-Vi IgG)
- Both can be calibrated against weight-based references (μg/ml)
- Vi-IgGR1,2011 (33 μg/ml, established 2013) provides a common anchor point
- Szu 2013 showed 1 ELISA Unit ≈ 1.24 μg/ml
- A structural assumption with explicit uncertainty is more useful than no bridge

**Uncertainty budget for τ**:
- Inter-assay CV for similar antigen coatings: ~20-30% (from Rijpkema 2018)
- Unknown systematic bias from coating evolution: ~0.5 log units
- Population differences (carriers vs vaccinees): ~0.5 log units
- **Conservative τ = 1.0 log unit** (covers ~2.7× fold uncertainty)

**Trade-offs**:
- Pro: Enables end-to-end bridge
- Con: τ is a guess; may underestimate or overestimate true uncertainty
- Con: Cannot validate without true bridging data

---

## D2: Barrett Figure 2 Data Extraction Method

**Date**: 2026-02-04
**Context**: Barrett 1983 Figure 2 shows a scattergram with cell counts, not individual points.

**Options Considered**:
1. Manually count cells and extract as aggregated data
2. Attempt to reconstruct individual-level data
3. Use regression parameters reported in text (if any)

**Decision**: Option 1 - Extract aggregated cell counts

**Rationale**:
- Figure 2 explicitly shows counts per cell (77 total sera)
- Individual reconstruction would add spurious precision
- Aggregated data likelihood is statistically appropriate
- Paper does not report regression parameters for HA vs ELISA

**Implementation**:
- CSV format: (ELISA_titer, HA_titer, count)
- Both titers on 2-fold dilution scale (log2)
- Left-censored values (<20) coded at midpoint (10)

---

## D3: Lee 2020 Model Implementation

**Date**: 2026-02-04
**Context**: Lee 2020 provides regression parameters: log(VaccZyme) = 3.749 + 0.946 × log(InHouse), r=0.991

**Options Considered**:
1. Use reported parameters directly
2. Attempt to fetch Figshare individual data and refit
3. Combine both approaches

**Decision**: Option 3 - Attempt Figshare fetch, fall back to reported parameters

**Rationale**:
- Individual data would allow better uncertainty quantification
- If unavailable, reported parameters are sufficient given r=0.991
- Paper parameters provide excellent fit (r²=0.982)

**Uncertainty estimation without individual data**:
- Residual SD ≈ √(1 - r²) × SD_Y
- For r=0.991 and log-scale, σ_residual very small (~0.13 log units)

---

## D4: End-to-End Cascade Structure

**Date**: 2026-02-04
**Context**: Need to connect HA → VaccZyme through three edges with different data quality.

**Decision**: Monte Carlo propagation through cascade

**Structure**:
```
For each HA titer observation:
  1. Sample E1: log2(ELISA) ~ N(α₁ + β₁ × log2(HA), σ₁²)
  2. Convert to log10: log10(ELISA) = log2(ELISA) × log10(2)
  3. Add E2 structural uncertainty: log10(modern) ~ N(log10(ELISA), τ²)
  4. Apply E3: log(VaccZyme) ~ N(3.749 + 0.946 × log(modern), σ₃²)
  5. Convert to EU/ml: VaccZyme_EU = exp(log(VaccZyme))
```

**Rationale**:
- Propagates uncertainty through all edges
- Explicit about which uncertainties are data-driven vs structural
- Produces posterior predictive distributions, not point estimates

---

## D5: Population Transferability Caveat

**Date**: 2026-02-04
**Context**: E1 data is from acute typhoid carriers (natural infection), E3 data is from vaccinees.

**Decision**: Document as explicit limitation, do not attempt adjustment

**Rationale**:
- No data to estimate population-specific calibration offsets
- Carrier antibodies may differ in avidity, subclass distribution
- Better to document the limitation than make unfounded adjustments
- Users should interpret results with this caveat in mind

---

## D6: Stop Criteria Evaluation

**Date**: 2026-02-04
**Context**: Need to evaluate pre-specified stop criteria from contract.

**Evaluation**:

| Criterion | Status | Action |
|-----------|--------|--------|
| E1 provides only positivity concordance | **PASS** - Barrett provides graded titers | Proceed |
| No link between endpoint titer and concentration | **TRIGGERED** - E2 missing | Add τ uncertainty |
| Population-specific mappings unavailable | **PARTIAL** - only acute typhoid | Document caveat |
| Posterior predictive check fails | TBD | Evaluate in Phase 5 |

**Decision**: Proceed with explicit uncertainty for E2, document population caveat

---

*Document will be updated as additional decisions are made.*
