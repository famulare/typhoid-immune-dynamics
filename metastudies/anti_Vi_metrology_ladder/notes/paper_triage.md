# Paper Triage

**Purpose**: Track inclusion/exclusion decisions and bridge roles for all papers in the corpus.

**Status**: Placeholder. To be populated during Phase 2 extraction.

---

## Triage Table

| Paper | Status | Paired Data? | Assays Compared | Population | N pairs | Data Format | Bridge Role | Notes |
|-------|--------|--------------|-----------------|------------|---------|-------------|-------------|-------|
| Barrett 1983a (ELISA development) | pending | TBD | HA ↔ early ELISA | acute/carrier | ? | ? | core? | In input_studies/ |
| Barrett 1983b (carrier identification) | pending | TBD | HA ↔ early ELISA | carrier | ? | ? | core? | In input_studies/ |
| Lee 2020 (PLOS NTD) | pending | yes | in-house ↔ VaccZyme | vaccinees | ~300 | Figshare | core | In input_studies/ |
| NIBSC 16/138 IFU | pending | commutability | multiple ELISAs | standards panel | ? | IFU text | structural | In input_studies/ |
| Rijpkema 2018 (Vaccine) | pending | ? | multiple ELISAs | standards panel | ? | ? | structural? | In input_studies/ |
| WHO BS.2017.2307 | pending | ? | multiple ELISAs | standards panel | ? | ? | structural? | May need to access |

---

## Bridge Role Definitions

- **Core**: Direct paired measurements on shared samples → primary calibration
- **Structural**: Provides commutability/equivalence constraints but not direct pairs
- **Support**: Informs heterogeneity priors, sensitivity analyses, or validation
- **Exclude**: No usable bridging data after extraction

---

## Exclusion Log

*Document reasons for exclusion*

| Paper | Reason for Exclusion | Reviewed By | Date |
|-------|---------------------|-------------|------|
| (none yet) | | | |

---

## Papers to Seek

*Additional papers that may provide bridging data*

| Target | Rationale | Status |
|--------|-----------|--------|
| Any paper with early ELISA ↔ modern ELISA paired data | Critical missing link | Not yet searched |
| WHO/NIBSC collaborative study raw data | May have per-lab potencies | Check if available |
