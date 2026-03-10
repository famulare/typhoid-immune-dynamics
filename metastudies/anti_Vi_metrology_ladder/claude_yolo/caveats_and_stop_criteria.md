# Caveats and Stop Criteria Evaluation

## Executive Summary

The anti-Vi metrology bridge from historical HA titers to modern VaccZyme ELISA values is **partially achievable** with significant caveats. The E2 link (Barrett 1983 ELISA → modern in-house ELISA) has NO direct data support, requiring structural assumptions that dominate the total uncertainty.

**Bottom line**: Use this bridge for rough order-of-magnitude comparisons only. Do not use for regulatory submissions, precise quantitative comparisons, or individual-level predictions.

---

## Stop Criteria Evaluation

### From Contract Specification

| Stop Criterion | Evaluation | Action Taken |
|----------------|------------|--------------|
| E1 provides only positivity concordance | **PASS** - Barrett 1983 provides graded titers for 77 sera with wide dynamic range | Proceeded with E1 model |
| No link between endpoint titer and concentration | **TRIGGERED** - E2 has no direct data | Added τ uncertainty parameter |
| Population-specific mappings unavailable | **PARTIAL** - E1 from carriers, E3 from vaccinees | Documented as caveat |
| Posterior predictive check fails | **N/A** - No validation data available | Cannot evaluate |

### Detailed Assessment

#### 1. E1 Data Quality (PASS)

Barrett 1983 Figure 2 provides:
- 77 paired HA/ELISA measurements
- Full range from <20 to 5120 for both assays
- Graded titers, not just positive/negative
- Correlation observable despite scatter

**E1 model fit**:
- R² ≈ 0.45-0.55 (substantial but not excellent)
- σ ≈ 0.88 log2 units (typical for dilution assays)
- Slope ≈ 0.80 (near-unity on log scale)

#### 2. E2 Gap (TRIGGERED)

Web searches (2026-02-04) confirmed NO published comparison between:
- Barrett 1983 ELISA (or any 1980s Vi ELISA)
- Modern in-house ELISAs (2010s+)

**Sources searched**:
- PubMed/PMC
- Google Scholar
- NIBSC technical reports
- Rijpkema 2018 collaborative study

**Mitigation applied**:
- Structural assumption: both assays calibratable to μg/ml
- Added τ = 1.0 log unit uncertainty (covers ~10× range)
- Documented as critical limitation

#### 3. Population Transferability (PARTIAL)

| Edge | Population | Issue |
|------|------------|-------|
| E1 | Acute typhoid carriers (El Salvador, 1980s) | May have different antibody characteristics than vaccinees |
| E2 | N/A | N/A |
| E3 | Vi-DT vaccinees (Philippines, 2010s) | Healthy volunteers, not natural infection |

**Known differences**:
- Carrier antibodies may have different avidity
- IgG subclass distribution may differ
- Epitope specificities could vary
- Duration of antibody presence differs

**No adjustment made** - insufficient data to estimate population-specific calibration offsets.

#### 4. Posterior Predictive Check (CANNOT EVALUATE)

A proper posterior predictive check requires:
- True paired measurements across the full bridge
- Held-out validation samples
- Independent replication

**None of these are available.**

The uncertainty estimates are therefore:
- Based on model structure, not empirical validation
- Likely underestimate true uncertainty
- Should be treated as lower bounds

---

## Critical Caveats

### 1. E2 Is Made Up

> **The E2 link has NO data support.**

The τ = 1.0 log unit uncertainty is a guess based on:
- Inter-assay CV for similar coatings (~20-30%)
- Assumed coating method differences (~0.5 log)
- Assumed population differences (~0.5 log)
- Conservative rounding up

**We cannot validate this assumption.**

### 2. Calibration Assumption Is Approximate

The cascade requires converting Barrett titers to μg/ml. We assumed:

> Barrett titer 160 ≈ 10 μg/ml

This is based on:
- Vi-IgGR1,2011 at 33 μg/ml having titers around 160-320 in similar assays
- Rough interpolation (not empirical calibration)

**This could be off by a factor of 2-5.**

### 3. Uncertainty Is Massive

The 95% prediction intervals span 2 orders of magnitude:

| HA Titer | Median VaccZyme | 95% PI |
|----------|-----------------|--------|
| 20 | 50 EU/ml | 5 - 457 EU/ml |
| 160 | 243 EU/ml | 25 - 2404 EU/ml |
| 1280 | 1142 EU/ml | 119 - 10604 EU/ml |

**E2 structural uncertainty dominates** (contributes ~70% of total variance).

### 4. Population Mismatch

The bridge combines:
- **E1**: Acute typhoid carriers (natural infection, diagnostic context)
- **E3**: Healthy TCV vaccinees (vaccine trial, immunogenicity context)

Antibody characteristics may differ between these populations:
- Avidity maturation differs
- IgG subclass distribution differs
- Epitope recognition patterns may differ

**No correction applied due to lack of data.**

### 5. No Validation Possible

Without:
- Modern reassay of historical samples
- Bridging study with concurrent measurements
- Independent replication

**We cannot assess accuracy.**

---

## Recommended Uses

### Acceptable Uses

1. **Rough order-of-magnitude comparisons**
   - "Historical HA titer 320 corresponds to approximately 400 EU/ml ± 1 order of magnitude"

2. **Historical context (qualitative)**
   - "Carriers with HA titers >640 likely had substantial anti-Vi IgG"

3. **Identifying data gaps**
   - "Direct bridging studies needed for the 1980s-2010s ELISA gap"

4. **Power analysis for future bridging studies**
   - "A bridging study with n=50 could reduce E2 uncertainty by ~50%"

### Unacceptable Uses

1. **Regulatory submissions**
   - E2 uncertainty is too large and unsupported

2. **Precise quantitative comparisons**
   - "HA 160 = 243 EU/ml" is misleadingly precise

3. **Individual-level predictions**
   - 95% PI spans 100× range

4. **Correlate of protection inference**
   - Cannot reliably map historical protective thresholds to modern units

---

## Uncertainty Budget

| Source | Magnitude (log units) | % of Total Variance |
|--------|----------------------|---------------------|
| E1 residual | 0.88 (log2) ≈ 0.61 (ln) | ~15% |
| E1 parameter | ~0.1-0.2 | ~2% |
| E2 structural (τ) | 1.0 | **~70%** |
| E3 residual | 0.30 | ~7% |
| Calibration assumption | Unknown | ~6%+ |

**Reducing E2 uncertainty would have the largest impact** on bridge precision.

---

## Future Directions

### To Improve This Bridge

1. **Find bridging data**
   - Check archived samples from 1980s studies
   - Contact original authors or their successors
   - Search gray literature and institutional reports

2. **Conduct a bridging study**
   - Reassay historical samples (if available) with modern ELISA
   - Direct comparison of in-house vs VaccZyme on carrier samples

3. **Better calibration**
   - Obtain Barrett ELISA reference materials (if extant)
   - Compare to precipitin analysis reference

### If This Bridge Is Insufficient

Consider alternative approaches:
1. **Ordinal mapping**: Map HA titers to qualitative categories only
2. **Parallel bridges**: Separate E1 and E3 without E2 link
3. **Expert elicitation**: Formal elicitation of E2 priors from assay experts
4. **Abandon quantitative bridge**: Accept that historical titers cannot be reliably converted

---

## Signatures

| Role | Agent | Date |
|------|-------|------|
| Author | Claude (YOLO) | 2026-02-04 |
| Reviewer | [Opus review pending] | |
| Approval | [User pending] | |

---

*This document should be reviewed by a domain expert before any use of the bridge results.*
