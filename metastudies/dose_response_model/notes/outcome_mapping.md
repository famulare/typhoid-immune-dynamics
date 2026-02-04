# Outcome Mapping: Observed Data → Model Variables

This document provides decision rules for mapping observed outcomes in the literature to the model variables defined in `dose_response_model_specification.md` (Section 1.4).

**Purpose**: Practical extraction guide for handling ambiguous cases during data coding.

**Related documents**:
- `dose_response_model_specification.md` - Formal definitions of latent states, observables, and observations
- `cross_cutting_observations.md` - Patterns across the extraction corpus

---

## Model Variables (Summary)

See `dose_response_model_specification.md` Section 1.4 for full definitions.

### Core Latent States
| Variable | Description |
|----------|-------------|
| Colonization | Successful establishment in the gut |
| Systemic invasion | Bacteria entering the bloodstream |
| Acute disease | Clinical illness (binary/graded) |
| Stool shedding | Fecal excretion of bacteria |
| Chronic carriage | Persistent gallbladder colonization (>1 year) |

### Observables
| Variable | Description |
|----------|-------------|
| Bacteremia | Bacteria in bloodstream (observable state) |
| Fever | Elevated temperature (ordinal, threshold-dependent) |
| Death | Fatal outcome |
| Immunity | Pre-challenge immune state |

### Observations (Measurements)
| Variable | Description |
|----------|-------------|
| Blood culture | Culture-based detection of bacteremia |
| Blood PCR | PCR-based detection of bacteremia |
| Stool culture | Culture-based detection of shedding |
| Stool PCR | PCR-based detection of shedding |
| Temperature | Thermometry reading |
| Serology titer | Antibody measurement |
| Seroconversion | Change in titer (derived) |
| Clinical typhoid dx | Composite endpoint (study-dependent) |

---

## Mapping Observed Outcomes to Latent States

### Colonization / Infection Proxies

| Observed Outcome | Maps To | Confidence | Notes |
|------------------|---------|------------|-------|
| Blood culture positive | Systemic invasion | High | Gold standard for invasion; implies colonization |
| Blood PCR positive | Systemic invasion | High | Higher sensitivity than culture |
| Stool culture positive | Stool shedding | High | Direct measure of shedding |
| Stool PCR positive | Stool shedding | High | Higher sensitivity than culture |
| Any positive (culture OR PCR) | Colonization | High | Any detection implies colonization occurred |
| Seroconversion only | Colonization | Moderate | Delayed; threshold-dependent; may miss mild cases |

**Decision rules**:
- **[USER-LOCKED]** Blood culture+ OR blood PCR+ → systemic invasion (and therefore colonization)
- **[USER-LOCKED]** Stool culture+ OR stool PCR+ → stool shedding (and therefore colonization)
- **[ASSISTANT-PROPOSED]** Seroconversion alone: flag as lower-confidence colonization proxy
- **[OPEN]** Studies reporting only seroconversion: include with sensitivity analysis flag?

### Fever / Acute Disease Proxies

| Observed Outcome | Maps To | Definition Clarity | Notes |
|------------------|---------|-------------------|-------|
| Oral temp ≥100.0°F (37.8°C) | Fever | Clear | Common older threshold |
| Oral temp ≥100.4°F (38.0°C) | Fever | Clear | CDC/modern threshold |
| Oral temp ≥101°F (38.3°C) | Fever | Clear | More conservative |
| Oral temp ≥103°F (39.4°C) sustained | Fever (severe) | Clear | Maryland "disease" threshold |
| "Clinical typhoid" | Acute disease | Variable | Composite; check definition |
| "Typhoid fever diagnosis" | Acute disease | Variable | May require fever + other criteria |
| "Typhoid fever" (unspecified) | Fever | Unclear | Verify paper definition |

**Decision rules**:
- **[USER-LOCKED]** Temperature threshold-based definitions map directly to fever
- **[USER-LOCKED]** "Clinical typhoid" or "typhoid diagnosis" maps to acute disease (composite)
- **[ASSISTANT-PROPOSED]** When fever duration specified, use as severity indicator
- **[OPEN]** Composite endpoints combining fever + bacteremia: stratify by component or treat as acute disease?

---

## Era-Specific Defaults

### Maryland Era (1960s-1970s)

| Reported Term | Default Mapping | Notes |
|---------------|-----------------|-------|
| "Disease" / "Ill" | Acute disease | Oral temp ≥103°F for 24-36hr + symptoms |
| "Typhoid fever" | Acute disease | Same as above; triggered chloramphenicol |
| "Infected" | Colonization | Usually means any positive culture |
| Attack rate | Fever (unless specified) | Check each paper |

### Oxford Era (2010s-present)

| Reported Term | Default Mapping | Notes |
|---------------|-----------------|-------|
| "Typhoid diagnosis" (TD) | Acute disease | Temp ≥38°C for ≥12hr OR blood culture+ |
| "Infection" | Colonization | TD OR PCR+ (includes subclinical) |
| "Definite typhoid" | Systemic invasion | Blood culture confirmed |
| Attack rate | Acute disease (TD) | Standard endpoint |

---

## Measurement Considerations

### Culture vs PCR

| Aspect | Culture | PCR |
|--------|---------|-----|
| Specificity | High | High |
| Sensitivity | Moderate (~50% for blood) | Higher |
| Affected by immunity | Yes (bactericidal Ab) | No |
| Affected by treatment | Yes | Less so |
| Available in | All eras | Oxford era only |

**Decision rule**: When both available, use "culture+ OR PCR+" as detection criterion.

### Temperature Thresholds

Different studies use different fever thresholds. For pooled analysis:

- **[ASSISTANT-PROPOSED]** Stratify by threshold when possible
- **[ASSISTANT-PROPOSED]** If pooling required, document threshold heterogeneity as limitation
- **[OPEN]** Develop measurement model for threshold differences?

---

## Mapping Ambiguities Log

*Record paper-specific decisions during Phase 4 joint review*

| Paper | Outcome as Reported | Proposed Mapping | Decision | Status |
|-------|---------------------|------------------|----------|--------|
| (to be filled during Phase 4) | | | | |

---

## Notes

- This document is updated as new outcome types are encountered during extraction
- Final mappings locked during Phase 4 joint review
- Ambiguous mappings may require sensitivity analysis (flag in Section 8 of specification)
- Cross-reference `dose_response_model_specification.md` Section 1.4 for formal probability statements
