# Anti-Vi Metrology Bridge: YOLO Implementation

## Executive Summary

This directory contains the autonomous ("YOLO") implementation of a metrology bridge from historical Vi hemagglutination (HA) titers to modern VaccZyme ELISA values.

**Bottom line**: A continuous bridge is partially achievable, but requires explicit structural uncertainty for a missing 40-year gap (E2) between Barrett's 1983 ELISA and modern in-house ELISAs.

## Bridge Graph

```
                    E1 (Barrett 1983)           E2 (MISSING)              E3 (Lee 2020)
Vi HA Titer ─────────────────────> Barrett ELISA ═══════════════════> Modern In-house ─────────────────────> VaccZyme EU
 (reciprocal)                       (reciprocal)    [structural        ELISA (μg/ml)                          (EU/ml)
                                                     assumption]
                                                         │
                                                         ▼
                                              E4 (NIBSC 16/138 IFU)
                                              VaccZyme → IS-relative IU
```

**Legend:**
- Solid arrows (───>): Data-supported edges with quantitative bridging
- Double lines (═══>): Structural assumption only (no direct data)

## Edge Summary

| Edge | Source | Target | Data | Quality | Status |
|------|--------|--------|------|---------|--------|
| E1 | Vi HA | Barrett ELISA 1983 | Barrett 1983 Fig 2 (77 pairs) | Digitizable | **AVAILABLE** |
| E2 | Barrett ELISA | Modern in-house | None | N/A | **MISSING** |
| E3 | In-house ELISA | VaccZyme | Lee 2020 (48 pairs, r=0.991) | Excellent | **AVAILABLE** |
| E4 | VaccZyme | IS-relative IU | NIBSC 16/138 IFU | Structural | **AVAILABLE** |

## Critical Gap: E2

The ~40-year gap between Barrett's 1983 ELISA and modern in-house ELISAs has **no direct bridging study**. Web searches (2026-02-04) found:
- No published comparison between early 1980s ELISAs and modern platforms
- No WHO/NIBSC technical reports with historical calibration chains
- The Rijpkema 2018 collaborative study (PMC6238147) only covers modern assays

### Mitigation Strategy

Both Barrett 1983 and Lee 2020 in-house ELISAs can theoretically be calibrated against weight-based units (μg/ml). From Szu 2013:
- Vi-IgGR1,2011: 33 μg/ml anti-Vi IgG by precipitin analysis
- Conversion: 1 ELISA Unit ≈ 1.24 μg/ml

This allows treating E2 as a calibration offset with added uncertainty (τ) representing:
1. Unknown systematic bias between assay generations
2. Coating antigen differences (glutaraldehyde-fixed Vi vs poly-L-lysine)
3. Antibody detection method evolution
4. Population-level differences between carriers (Barrett) and vaccinees (Lee)

**Default assumption**: τ = 1 log unit (i.e., ~2.7× uncertainty factor)

## Directory Contents

```
claude_yolo/
├── README.md                 # This file
├── decision_log.md           # All decisions with rationale
├── bridge_graph.md           # Detailed bridge graph with citations
├── caveats_and_stop_criteria.md  # Limitations and stop conditions
├── extraction/               # Paper data extraction
│   ├── barrett_1983_extraction.md
│   ├── barrett_1983_fig2_data.csv
│   ├── lee_2020_extraction.md
│   └── nibsc_16138_extraction.md
├── analysis/                 # R analysis scripts
│   ├── e1_ha_elisa_bridge.R
│   ├── e3_inhouse_vacczyme_bridge.R
│   └── end_to_end_bridge.R
├── intermediates/            # Saved model objects and data
│   ├── e2_search_results.md
│   └── [model fits, samples, etc.]
├── outputs/                  # Final outputs
│   └── parameter_estimates.csv
└── reviews/                  # Opus reviews and responses
    ├── opus_review_phase1.md
    ├── reviewer2_critique.md
    └── author_response.md
```

## Key Sources

1. **Barrett 1983** (Part 1 & 2): Original Vi HA vs ELISA comparison in typhoid carriers
2. **Lee 2020**: In-house ELISA vs VaccZyme comparison (r=0.991, 48 samples)
3. **Rijpkema 2018**: NIBSC 16/138 collaborative study
4. **Szu 2013**: Vi-IgGR1,2011 reference establishment

## Implementation Status

- [x] Phase 0: Web search for E2 bridge data (confirmed missing)
- [x] Phase 1: Infrastructure committed
- [x] Phase 2: Directory structure created
- [x] Phase 3: Data extraction
- [x] Phase 4: Bridge models (E1: R²≈0.45, E3: r=0.991)
- [x] Phase 5: End-to-end cascade with uncertainty propagation
- [x] Phase 6: Caveats documentation
- [ ] Phase 7: Opus review

## Key Results

### Bridge Conversion Table (Monte Carlo, n=5000)

| HA Titer | VaccZyme Median (EU/ml) | 95% Prediction Interval |
|----------|-------------------------|-------------------------|
| 20 | 50 | 5 - 457 |
| 40 | 86 | 9 - 848 |
| 80 | 144 | 14 - 1359 |
| 160 | 243 | 25 - 2404 |
| 320 | 412 | 43 - 3939 |
| 640 | 713 | 71 - 6252 |
| 1280 | 1142 | 119 - 10604 |

**WARNING**: 95% prediction intervals span ~2 orders of magnitude due to E2 structural uncertainty (τ=1.0 log unit). Use for rough comparisons only.

## Known Limitations

1. **E2 gap**: No direct bridge between Barrett ELISA and modern in-house ELISAs
2. **Population mismatch**: E1 data is acute typhoid carriers; E3 data is vaccinees
3. **No validation cohort**: Cannot perform true out-of-sample validation
4. **Barrett Figure 2**: Aggregated counts, not individual paired measurements

---
*Last updated: 2026-02-04*
