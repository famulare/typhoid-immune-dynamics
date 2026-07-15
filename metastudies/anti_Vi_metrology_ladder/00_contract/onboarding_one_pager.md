# Anti-Vi Metrology Bridge: Getting Your Bearings

**What are we doing?**

Building an empirical conversion cascade from 1960s-80s Vi hemagglutination (HA) titers to modern anti-Vi IgG ELISA (VaccZyme or IS-relative). This lets us use historical immunity data in contemporary dose-response models.

**Why is this hard?**

We have a 50-year gap between assay technologies with sparse bridging data:

| Era | Assay Type | Strengths | Gaps |
|-----|------------|-----------|------|
| **1960s-80s** | Vi HA/agglutination | Abundant historical data (Maryland challenge studies) | No direct link to modern ELISA; titer = endpoint dilution |
| **1980s-90s** | Early ELISA | Barrett 1983 compared to HA | Not standardized; "titer" still endpoint dilution |
| **2000s+** | VaccZyme/IS-calibrated | Oxford CHIM, vaccine trials; standardized units | No direct comparison to legacy HA |

The cascade must bridge: `HA → early ELISA → modern ELISA (VaccZyme/IS)`

Each arrow requires paired measurements on the same sera—which may not exist.

**The bridge structure**

Log-linear conversion with uncertainty:
```
log(titer_B) = β₀ + β₁ * log(titer_A) + ε
ε ~ Normal(0, σ²)
```

Per edge:
- **β₀**: Systematic offset (unit/scale differences)
- **β₁**: Proportionality (ideally ~1 if assays measure same thing)
- **σ²**: Combined measurement error + biological noise

End-to-end: propagate distributions through the cascade, not point estimates.

**Key documents to read**

1. `notes/metrology_bridge_specification.md` — Reference model and bridge graph
2. `notes/cross_cutting_observations.md` — Patterns across the literature
3. `notes/assay_mapping.md` — How reported assays map to canonical nodes

**Key papers**

- **Barrett 1983**: Developed Vi ELISA; compared to HA on shared sera
- **Lee 2020 (PLOS NTD)**: Paired in-house ELISA and VaccZyme; Figshare dataset
- **NIBSC 16/138 IFU**: Establishes ELISA commutability; IS calibration framework
- **Rijpkema 2018 (Vaccine)**: Published IS work; VaccZyme as reference

**The big open questions**

1. **HA → ELISA grading**: Does Barrett 1983 provide continuous paired data or only positive/negative concordance?

2. **Early → modern ELISA link**: Is there any paper with paired early-style ELISA and VaccZyme? Or must we assume "generic ELISA" commutability?

3. **Population dependence**: Do carriers, vaccinees, and acute cases have the same HA-ELISA relationship, or do we need stratification?

**What success looks like**

A calibrated cascade that:
- Converts historical HA titers to VaccZyme-equivalent or IS-relative distributions
- Has quantified uncertainty at each step
- Explicitly flags where bridging degrades to ordinal bins
- Can serve as a prior on immunity in dose-response calibration

The uncertainty will be large. That's honest. The goal is a model that's *usable* for historical data integration, with transparent limitations.

---

*Start with Barrett 1983 extraction, then the Lee 2020 Figshare data, then the NIBSC documentation to understand commutability.*
