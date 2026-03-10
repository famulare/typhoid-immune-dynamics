# Anti-Vi Metrology Bridge: Detailed Graph

## Visual Representation

```
┌─────────────────────────────────────────────────────────────────────────────────────────────┐
│                          ANTI-Vi METROLOGY BRIDGE (YOLO IMPLEMENTATION)                      │
├─────────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                              │
│   ┌──────────────┐         E1          ┌──────────────┐        E2         ┌──────────────┐ │
│   │   Vi HA      │═══════════════════▶│ Barrett ELISA │══════════════════▶│ Modern       │ │
│   │   Titer      │  Barrett 1983       │   Titer      │   STRUCTURAL      │ In-house     │ │
│   │ (reciprocal) │  n=77 pairs         │ (reciprocal) │   ASSUMPTION      │ ELISA        │ │
│   │              │  acute carriers     │              │   τ = 1.0 log     │ (μg/ml)      │ │
│   └──────────────┘                     └──────────────┘                   └──────┬───────┘ │
│                                                                                   │         │
│                                                                                   │ E3      │
│                                                                          Lee 2020 │         │
│                                                                          n=48     │         │
│                                                                          r=0.991  │         │
│                                                                                   ▼         │
│                                                                           ┌──────────────┐ │
│                                                                           │  VaccZyme    │ │
│                                                                           │  ELISA       │ │
│                                                                           │  (EU/ml)     │ │
│                                                                           └──────┬───────┘ │
│                                                                                   │         │
│                                                                                   │ E4      │
│                                                                       NIBSC 16/138│         │
│                                                                          IFU      │         │
│                                                                                   ▼         │
│                                                                           ┌──────────────┐ │
│                                                                           │    IU/ml     │ │
│                                                                           │ (IS-relative)│ │
│                                                                           └──────────────┘ │
│                                                                                              │
├─────────────────────────────────────────────────────────────────────────────────────────────┤
│ Legend:  ═══▶ Data-supported edge    ═══▶ Structural assumption (no direct data)            │
└─────────────────────────────────────────────────────────────────────────────────────────────┘
```

## Edge Details

### E1: Vi HA → Barrett ELISA

| Property | Value |
|----------|-------|
| **Source Node** | Vi Hemagglutination Titer |
| **Target Node** | Barrett ELISA Titer |
| **Data Source** | Barrett et al. 1983, J Clin Microbiol |
| **Figure/Table** | Figure 2 (Part 1, p. 626) |
| **Sample Size** | 77 paired sera |
| **Population** | Acute typhoid patients and carriers (Washington state outbreak) |
| **Scale** | Both log2 (2-fold dilution series) |
| **Lower Limit** | <20 (coded as 10 for midpoint) |
| **Model** | log2(ELISA) ~ Normal(α + β × log2(HA), σ²) |

**Citation**: Barrett TJ, Grubb B, Gruschkau H, et al. Enzyme-linked immunosorbent assay for detection of human antibodies to Salmonella typhi Vi antigen. J Clin Microbiol. 1983;18(3):625-632.

### E2: Barrett ELISA → Modern In-house ELISA

| Property | Value |
|----------|-------|
| **Source Node** | Barrett ELISA Titer (1983) |
| **Target Node** | Modern In-house ELISA (2010s) |
| **Data Source** | **NONE** |
| **Sample Size** | N/A |
| **Model** | log10(modern) ~ Normal(log10(Barrett) + μ_offset, τ²) |
| **Structural Assumption** | Common calibration to weight-based units |
| **Uncertainty (τ)** | 1.0 log unit (default) |

**Justification for τ = 1.0**:
- Inter-assay CV for similar coatings: ~20-30%
- Coating method differences (glutaraldehyde vs poly-L-lysine): ~0.5 log
- Population differences (carriers vs vaccinees): ~0.5 log
- Conservative total: 1.0 log unit covers ~10-fold uncertainty range

**Common anchor**: Both assays can theoretically be calibrated to μg/ml via:
- Precipitin analysis with purified Vi (Barrett era)
- Vi-IgGR1,2011 reference (modern era, 33 μg/ml = 1.0 ELISA Unit × 33)

### E3: Modern In-house ELISA → VaccZyme

| Property | Value |
|----------|-------|
| **Source Node** | Modern In-house ELISA (μg/ml) |
| **Target Node** | VaccZyme ELISA (EU/ml) |
| **Data Source** | Lee et al. 2020, PLoS Negl Trop Dis |
| **Figure/Table** | Figure 1, Table 2 |
| **Sample Size** | 48 paired sera |
| **Population** | Vi-DT vaccinees (Phase 1 trial) |
| **Correlation** | r = 0.991 (Pearson, log scale) |
| **Model** | log(VaccZyme) = 3.749 + 0.946 × log(InHouse) + ε |
| **Residual SD** | σ ≈ 0.13 log units (estimated from r) |

**Citation**: Lee C, Yang JS, Park HK, et al. Comparison of anti-Vi IgG responses between two clinical studies of typhoid Vi conjugate vaccines (Vi-DT vs Vi-TT). PLoS Negl Trop Dis. 2020;14(4):e0008171.

### E4: VaccZyme → IS-relative IU

| Property | Value |
|----------|-------|
| **Source Node** | VaccZyme EU/ml |
| **Target Node** | IU/ml (relative to NIBSC 16/138) |
| **Data Source** | NIBSC 16/138 Instructions for Use |
| **Conversion** | Structural (calibrator-defined) |
| **Reference** | NIBSC 16/138 = 100 IU/ampoule |

**Note**: The VaccZyme kit calibrators are already traceable to international standards, so this edge is primarily for documentation completeness.

**Citation**: NIBSC. Instructions for Use: 1st International Standard for Anti-Typhoid Capsular Vi Polysaccharide IgG (Human). NIBSC code 16/138.

## Unit Conversions

| Source Unit | Target Unit | Conversion |
|-------------|-------------|------------|
| HA titer (reciprocal) | log2(HA) | log2 transformation |
| ELISA titer (reciprocal) | log2(ELISA) | log2 transformation |
| log2(ELISA) | log10(ELISA) | multiply by log10(2) ≈ 0.301 |
| In-house μg/ml | log(μg/ml) | natural log |
| log(VaccZyme EU) | VaccZyme EU/ml | exp transformation |
| VaccZyme EU/ml | IU/ml | 1:1 when calibrated to 16/138 |

## Uncertainty Budget

| Edge | Source | Magnitude | Notes |
|------|--------|-----------|-------|
| E1 | Regression residual | σ₁ ~ 0.5-1.0 log2 | To be estimated from data |
| E2 | Structural assumption | τ = 1.0 log10 | No data; conservative guess |
| E3 | Regression residual | σ₃ ~ 0.13 log | From r = 0.991 |
| E4 | Calibrator precision | Negligible | Kit-defined |

**Total uncertainty** (end-to-end) dominated by:
1. E2 structural uncertainty (τ)
2. E1 measurement uncertainty (σ₁)
3. E3 measurement uncertainty (σ₃)

## Populations and Transferability

| Edge | Population | Transferability Concern |
|------|------------|------------------------|
| E1 | Acute typhoid carriers | May not apply to vaccinees |
| E2 | N/A | N/A |
| E3 | Vi-DT vaccinees | May not apply to carriers |
| E4 | Mixed (collaborative study) | Broadly applicable |

**Critical caveat**: The bridge mixes carrier data (E1) with vaccinee data (E3). Antibody characteristics may differ between these populations (avidity, subclass, epitope specificity).
