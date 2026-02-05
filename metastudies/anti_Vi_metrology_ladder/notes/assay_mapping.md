# Assay Mapping: Reported Assays → Canonical Nodes

**Purpose**: Define how specific assays reported in the literature map to canonical nodes in the bridge graph.

**Status**: Initial draft. To be refined during Phase 2 extraction.

---

## Canonical Node Definitions

### LEGACY ASSAYS (pre-1990)

| Canonical Node | Description | Units | Typical LOD |
|----------------|-------------|-------|-------------|
| **Vi_HA** | Vi hemagglutination or passive hemagglutination | Endpoint dilution titer (1:X) | 1:10 or 1:20 |

**Mapping rules:**
- `Vi hemagglutination` → Vi_HA
- `Vi agglutination` → Vi_HA [ASSISTANT-PROPOSED: verify same technique]
- `Vi passive hemagglutination` → Vi_HA
- `indirect hemagglutination for Vi` → Vi_HA

### EARLY ELISA (1980s-1990s)

| Canonical Node | Description | Units | Typical LOD |
|----------------|-------------|-------|-------------|
| **early_ELISA** | Pre-standardization ELISA, endpoint dilution | Endpoint dilution titer (1:X) | Varies |

**Mapping rules:**
- `Barrett ELISA 1983` → early_ELISA (reference implementation)
- `in-house Vi ELISA, endpoint titer` → early_ELISA [requires protocol review]
- `Vi IgG ELISA, titer` → early_ELISA [if pre-2000 and endpoint dilution]

### MODERN ELISA (2000s+)

| Canonical Node | Description | Units | Typical LOD |
|----------------|-------------|-------|-------------|
| **VaccZyme** | Commercial VaccZyme anti-Vi IgG ELISA | EU/mL (kit units) | ~7.4 EU/mL |
| **IS_ELISA** | Any ELISA calibrated to IS or Vi-IgGR1 | IU/mL or relative potency | Varies |
| **modern_ELISA** | In-house ELISA, concentration-based, not IS-calibrated | µg/mL or arbitrary units | Varies |

**Mapping rules:**
- `VaccZyme anti-Vi IgG` → VaccZyme
- `Binding Site VaccZyme` → VaccZyme
- `in-house ELISA calibrated to NIBSC 16/138` → IS_ELISA
- `in-house ELISA calibrated to Vi-IgGR1, 2011` → IS_ELISA
- `Oxford in-house ELISA` → modern_ELISA [then check if VaccZyme-bridged]
- `ELISA, µg/mL or EU/mL` → modern_ELISA [then verify calibration]

### REFERENCE STANDARDS

| Node | Description | Role |
|------|-------------|------|
| **IS_16_138** | NIBSC 1st International Standard | Defines IU scale |
| **Vi_IgGR1_2011** | US reference reagent (FDA/CBER) | Alternative reference anchor |

---

## Decision Log

*Record mapping decisions here during extraction*

| Paper | Reported Assay | Assigned Node | Decision Tag | Notes |
|-------|----------------|---------------|--------------|-------|
| Barrett 1983 | "ELISA for Vi antibody" | early_ELISA | [ASSISTANT-PROPOSED] | Reference implementation |
| Barrett 1983 | "hemagglutination assay" | Vi_HA | [ASSISTANT-PROPOSED] | Cites prior method |
| Lee 2020 | "in-house ELISA" | modern_ELISA | [ASSISTANT-PROPOSED] | Bridged to VaccZyme |
| Lee 2020 | "VaccZyme" | VaccZyme | [USER-LOCKED] | Explicit |
| NIBSC 16/138 | "VaccZyme" | VaccZyme | [USER-LOCKED] | Explicit in IFU |
| ... | ... | ... | ... | ... |

---

## Ambiguous Cases

### When assay type is unclear

If a paper reports "Vi antibody" or "anti-Vi titer" without method details:
1. Check methods section for protocol
2. Check citations for source method
3. If still unclear, flag as [OPEN] and assign to most likely node with wide uncertainty

### When ELISA calibration is unclear

If a paper reports ELISA in concentration units but doesn't state calibration:
1. Check if reference standard mentioned anywhere
2. Check era (pre-2010 unlikely to use IS)
3. If concentration-based but no IS → modern_ELISA
4. If endpoint titer → early_ELISA

### When multiple assays share a name

E.g., "VaccZyme" may refer to different kit versions over time:
1. Record publication year
2. Note lot/version if stated
3. Flag for sensitivity analysis if spanning major kit revisions
