# Cross-Cutting Observations

**Purpose**: Document patterns observed across the extraction corpus that inform bridging model design.

**Status**: Placeholder. To be populated during Phase 2 extraction.

---

## Assay Lineage

*Which ELISAs descended from Barrett 1983? Which use NIBSC standards?*

- To be filled

---

## Commutability Evidence

*What does the IS documentation say about assay equivalence?*

- NIBSC 16/138 IFU states commutability for "indirect in-house ELISAs of native Vi" with VaccZyme
- Specific conditions: [to be extracted]
- Non-commutable formats: [to be extracted]

---

## Population Effects

*Do bridging relationships differ by population type?*

| Population | Expected Effect | Evidence |
|------------|-----------------|----------|
| Acute typhoid | Mixed isotype; may inflate HA relative to IgG ELISA | TBD |
| Carriers | High avidity IgG; may be concordant across assays | TBD |
| Vaccinated | Predominantly IgG; should align well with ELISA | TBD |
| Endemic baseline | Variable; depends on exposure history | TBD |

---

## Era Patterns

*How did titer definitions and reporting change over time?*

| Era | Typical Reporting | Notes |
|-----|-------------------|-------|
| 1960s-70s | HA titer (1:X) | Endpoint dilution; 2-fold series |
| 1980s | ELISA titer (1:X) or OD ratio | Transition period |
| 1990s-2000s | ELISA concentration (µg/mL, EU/mL) | Move toward standardization |
| 2010s+ | ELISA (EU/mL, IU/mL) | IS-calibrated increasingly common |

---

## Antigen/Coating Effects

*Native Vi vs protein-coated vs conjugated*

- Native Vi: Most legacy and early ELISA
- Protein pre-/co-coating: Some modern in-house ELISAs; noted in NIBSC commutability
- Conjugate vaccines: May elicit different epitope profiles

---

## Recurring Data Gaps

*What's consistently missing across papers?*

- Individual-level paired data (often only summary statistics)
- Explicit LOD/LOQ values
- Assay CV estimates
- Protocol details for legacy methods

---

## Key Uncertainties for Bridge Model

1. **HA → ELISA granularity**: Is Barrett 1983 data ordinal or continuous?
2. **Early ELISA → modern ELISA link**: Does any direct bridge exist?
3. **VaccZyme kit version stability**: How much has the assay drifted?
4. **Population stratification feasibility**: Enough data per population?
