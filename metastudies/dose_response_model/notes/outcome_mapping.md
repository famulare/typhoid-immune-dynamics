# Outcome Mapping: Observed Outcomes → Model Variables

This document maps the various outcome measures reported in the literature to the two primary model endpoints: **infection** and **fever**.

---

## Model Endpoints

### Infection
The model's infection endpoint represents successful colonization/invasion by S. Typhi following oral challenge. This is a necessary precursor to symptomatic disease.

### Fever
The model's fever endpoint represents symptomatic typhoid fever. In the model, P(fever|dose) ≤ P(infection|dose) always, since fever requires infection.

---

## Infection Proxies

| Observed Outcome | Maps To | Specificity | Sensitivity | Notes |
|------------------|---------|-------------|-------------|-------|
| Blood culture positive | infection | High | Low-moderate | Gold standard but misses ~50% of infections |
| PCR positive (blood) | infection | High | Higher than culture | Modern studies only |
| Stool shedding | infection | Moderate | Variable | Depends on sampling schedule; transient shedding possible |
| Seroconversion | infection | Moderate | Moderate | Delayed response; threshold-dependent |

### Decision Rules (to be refined during extraction)

- **[ASSISTANT-PROPOSED]** Blood culture OR PCR positive = infection
- **[ASSISTANT-PROPOSED]** Stool shedding alone: flag as lower-confidence infection proxy
- **[OPEN]** How to handle studies that only report seroconversion?

---

## Fever Proxies

| Observed Outcome | Maps To | Definition Clarity | Notes |
|------------------|---------|-------------------|-------|
| Oral temp ≥100.0°F (37.8°C) | fever | Clear | Common older threshold |
| Oral temp ≥100.4°F (38.0°C) | fever | Clear | CDC/modern threshold |
| Oral temp ≥101°F (38.3°C) | fever | Clear | More conservative |
| "Clinical typhoid" | fever | Variable | Often composite; definition varies |
| "Typhoid fever diagnosis" | fever | Variable | May include fever + other criteria |
| "Typhoid fever" (unspecified) | fever | Unclear | Need to check paper definition |

### Decision Rules (to be refined during extraction)

- **[ASSISTANT-PROPOSED]** Temperature threshold-based definitions preferred
- **[ASSISTANT-PROPOSED]** "Clinical typhoid" or "typhoid fever" mapped to fever endpoint unless definition clearly excludes fever
- **[OPEN]** How to handle composite endpoints that combine fever + bacteremia?

---

## Outcome Definitions by Study Era

### Maryland Era (1960s-1970s)
- Typical definition: Oral temperature ≥100°F for specified duration
- Often reported as "typhoid fever" = fever + systemic symptoms
- Need to verify each paper's exact definition

### Oxford Era (2010s-present)
- Typically uses "typhoid diagnosis" composite:
  - Oral temperature ≥38°C sustained for ≥12 hours, OR
  - Blood culture positive for S. Typhi
- PCR increasingly used alongside culture
- More standardized definitions across papers

---

## Mapping Ambiguities Log

| Paper | Outcome as Reported | Proposed Mapping | Status |
|-------|---------------------|------------------|--------|
| (to be filled during extraction) | | | |

---

## Notes

- This document will be updated as we encounter new outcome types during extraction
- Final mappings will be locked during Phase 3 joint review
- Ambiguous mappings may require sensitivity analysis in calibration
