# Lee 2020 Data Extraction

## Source

Lee EY, Park JY, Kim DR, Song M, Sahastrabuddhe S, Kim H, Chon Y, Yang JS. Comparison of anti-Vi IgG responses between two clinical studies of typhoid Vi conjugate vaccines (Vi-DT vs Vi-TT). **PLoS Negl Trop Dis.** 2020;14(3):e0008171.

DOI: https://doi.org/10.1371/journal.pntd.0008171

## Key Finding: E3 Bridge Parameters

### Regression: In-house ELISA → VaccZyme ELISA

From **Figure 2A** (page 7):

```
log(VaccZyme EU/ml) = 3.749 + 0.946 × log(In-house μg/ml)
```

| Parameter | Value | 95% CI |
|-----------|-------|--------|
| Intercept | 3.749 | Not reported |
| Slope | 0.946 | Not reported |
| Pearson r | 0.991 | P < 0.001 |
| R² | 0.982 | (from r²) |
| n | 48 | |

**Note**: The regression uses natural log (base e) on both axes.

### Alternative Regression (Vi IgGR1,2011 as reference)

From **Figure 1A** (page 6):

```
log(VaccZyme μg/ml) = 0.337 + 0.971 × log(In-house μg/ml)
```

| Parameter | Value |
|-----------|-------|
| Intercept | 0.337 |
| Slope | 0.971 |
| Pearson r | 0.991 |

This alternative uses the same units (μg/ml) on both axes, showing near-unity relationship when both are calibrated to Vi-IgGR1,2011.

## Sample Characteristics

### 48 Bridging Samples (Table 1)

| Age Group | Low | Medium | High | Total |
|-----------|-----|--------|------|-------|
| Adults (18-45 yrs) | 5 | 5 | 5 | 15 |
| Adolescents (6-17 yrs) | 4 | 6 | 6 | 16 |
| Children (2-5 yrs) | 5 | 7 | 5 | 17 |
| **Total** | 14 | 18 | 16 | **48** |

**Stratification criteria**:
- Low: Anti-Vi IgG < 20.279 μg/ml
- Medium: 20.279 ≤ Anti-Vi IgG < 51.143 μg/ml
- High: Anti-Vi IgG ≥ 51.143 μg/ml

### Summary Statistics (Table 2)

| Assay | GMT ± SD | Median | Min | Max | CV |
|-------|----------|--------|-----|-----|-----|
| In-house (μg/ml) | 10.293 ± 14.574 | 31.107 | 0.033 | 187.588 | 1.042 |
| VaccZyme (μg/ml) | 13.469 ± 13.807 | 39.430 | 0.020 | 246.350 | 0.971 |
| VaccZyme (EU/ml) | 385.492 ± 12.927 | 1191.135 | 0.992 | 5694.640 | 0.881 |

**Note**: SD here is geometric SD (multiplicative factor).

## Uncertainty Estimation

### Residual Standard Deviation

From r = 0.991:
```
r² = 0.982
1 - r² = 0.018
```

If we assume SD_Y (log VaccZyme) ≈ 2.5 log units (estimated from range):
```
σ_residual ≈ √(1 - r²) × SD_Y ≈ √0.018 × 2.5 ≈ 0.34 log units
```

More conservatively, from the data range:
```
log(5694.64) - log(0.992) ≈ 8.65 - 0 = 8.65 log units
SD_Y ≈ 8.65 / 4 ≈ 2.16  (assuming range ≈ 4 SD)
σ_residual ≈ √0.018 × 2.16 ≈ 0.29 log units
```

**Conservative estimate**: σ_residual ≈ 0.3 log units (natural log)

In log10 scale: σ_residual ≈ 0.3 / ln(10) ≈ 0.13 log10 units

## Assay Methods

### In-house ELISA (IVI)
- Plates precoated with 10 μg/ml poly-L-lysine
- Vi antigen (SK Bioscience): 2 μg/ml absorbed overnight at 37°C
- AP-conjugated mouse anti-human IgG detection
- Calibrated against Vi-IgGR1,2011 (33 μg/ml)
- Designated "in-house ELISA 1" in Rijpkema 2018 collaborative study

### VaccZyme ELISA
- Commercial kit (Binding Site, UK)
- Kit calibrators or Vi-IgGR1,2011 as reference
- Calibrator range: 7.4 - 600 EU/ml
- Sensitivity: 7.4 EU/ml

## Population

- Vi-DT Phase 1 trial participants (Philippines)
- Healthy Filipino adults and children
- Pre- and post-vaccination samples (Day 0, Day 28)
- Excluded: elevated liver function, other vaccine history, extreme values

## Key Limitations for Bridging

1. **Vaccinee population**: All samples from TCV trial participants, not natural infection
2. **Selection bias**: Samples stratified by titer level, excluded extremes
3. **Age effects**: Model includes age-group interaction terms (not used in simple regression)
4. **Individual data not available**: Only summary statistics and regression parameters in paper

## Extracted Parameters for E3 Model

```r
# E3 Bridge: In-house ELISA → VaccZyme ELISA
e3_intercept <- 3.749  # log scale
e3_slope <- 0.946
e3_r <- 0.991
e3_n <- 48
e3_sigma_residual <- 0.30  # log units (natural log), conservative estimate

# Conversion function
# log(VaccZyme_EU) = e3_intercept + e3_slope * log(InHouse_ug)
# VaccZyme_EU = exp(e3_intercept) * InHouse_ug^e3_slope
```

## Reference Standard Information

### Vi-IgGR1,2011 (CBER Reference)
- Concentration: 33 μg/ml anti-Vi IgG
- Established by precipitin analysis
- Conversion: 1 EU ≈ 1.24 μg/ml (from Szu 2013)

### NIBSC 16/138 (International Standard)
- Assigned: 100 IU/ampoule
- Vi-IgGR1,2011: 163 IU relative to 16/138

---
*Extracted 2026-02-04*
