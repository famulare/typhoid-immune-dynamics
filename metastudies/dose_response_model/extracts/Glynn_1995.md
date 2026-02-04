# Glynn et al. 1995 - Infecting dose and severity of typhoid: analysis of volunteer data

## Citation
Glynn JR, Hornick RB, Levine MM, Bradley DJ. Infecting dose and severity of typhoid: analysis of volunteer data and examination of the influence of the definition of illness used. Epidemiol Infect. 1995;115(1):23-30.

## Study Overview
- **Primary question**: Explore relationship between challenge dose and severity of disease; examine how illness definition affects results
- **Study type**: Retrospective reanalysis of Maryland CHIM data using original patient charts
- **Data source**: Chart review from studies conducted 1959-1970 at University of Maryland
- **Population**: Unvaccinated volunteers, excluding those with prior S. typhi challenge

**[CRITICAL]** This paper provides the most detailed individual-level data from the Maryland studies, with multiple illness definitions analyzed.

## Subject Characteristics

### Demographics (p.25)
- **N total**: 278 volunteers challenged with 10^3 to 10^9 organisms
- **Age**: 20-54 years (median 28)
- **Race**: 116 Black, 153 White (9 unrecorded)
- **Military service**: 102 had military service (unknown for 70); service ended mean 12 years (SD 5.9) previously
- **[OPEN]** Military service implies prior typhoid vaccination - treated as confounder

### Risk Factors for Attack (using Definition I, p.25)
- **Age**: <30 years higher attack rate (82/145 vs 37/117; RR 1.79, 95% CI 1.32-2.42)
- **Military service**: Without service higher attack rate (57/105 vs 30/96; RR 1.74, 95% CI 1.23-2.45)
- **Race**: Black higher attack rate (59/114 vs 56/142; RR 1.31, 95% CI 1.00-1.72)
- These associations were independent of each other

## Challenge Protocol

### Strain & Delivery
- **Strain**: Quailes strain S. typhi (phase type D-1)
- **Delivery**: "Administered in milk and gargled before swallowing"
- **Treatment trigger**: "Antibiotics (usually chloramphenicol) were started when the temperature was at least 103°F (39.4°C) for 24-36 h"

### Study Period Changes
- "From 1965, most volunteers were seen as outpatients unless and until they became ill"
- "All those challenged after 1965 received 10^5 organisms, compared to only 12 earlier"

## Illness Definitions (p.24)

Four definitions compared:

| Definition | Criteria |
|------------|----------|
| **(I)** | Temperature ≥100°F (37.8°C) for ≥12 h continuously AND peak ≥101°F (38.3°C) |
| **(II)** | Illness that resulted in antibiotic therapy |
| **(III)** | Temperature ≥101°F (38.3°C) for ≥12 h continuously |
| **(IV)** | Temperature spiking to ≥103°F (39.4°C) for ≥36 h |

**[CRITICAL]** Definition IV matches the original Hornick papers' criterion. Definition I is more sensitive.

## Attack Rate Data (Table 1, p.26)

**[KEY DATA - LARGEST N, MULTIPLE DEFINITIONS]**

| Dose | Def I (sensitive) | Def II (Abx given) | Def III | Def IV (strict, Hornick) |
|------|-------------------|-------------------|---------|--------------------------|
| 10^3 | 0/13 (0%) | 0/13 (0%) | 0/13 (0%) | 0/13 (0%) |
| 10^5 | 83/200 (41.5%) | 72/204 (35.3%) | 77/200 (38.5%) | 52/200 (26.0%) |
| 10^7 | 13/27 (48.1%) | 13/27 (48.1%) | 13/27 (48.1%) | 8/27 (29.6%) |
| 10^8-10^9 | 24/25 (96.0%) | 30/34 (88.2%) | 24/25 (96.0%) | 18/25 (72.0%) |

**Notes**:
- 10^5 group: N=200-204 (largest by far)
- Definition I attack rate at 10^5 (41.5%) is HIGHER than Hornick papers reported (28%)
- **[OPEN]** The 10^5 Hornick papers showed 32/116 = 28%; here it's 83/200 = 41.5%. Different denominators suggest additional subjects or different inclusion criteria.

## Incubation Period and Severity (Table 2, p.26)

For ill volunteers using Definition I:

| Dose | N ill | Incubation (geom. mean, days) | Peak temp (°F) | Time >103°F (geom. mean 4h periods) | Symptom score |
|------|-------|-------------------------------|----------------|-------------------------------------|---------------|
| 10^5 | 83 | 9.3 (8.4-10.4) | 103.8 (103.6-104.0) | 3.2 (2.5-4.0) | 7.7 (7.2-8.2) |
| 10^7 | 13 | 7.4 (4.9-11.2) | 103.6 (103.1-104.1) | 2.6 (1.4-5.0) | 8.9 (7.8-10.0) |
| 10^8-10^9 | 24 | 4.7 (4.1-5.4) | 104.4 (104.1-104.7) | 4.8 (3.4-6.7) | 9.2 (8.3-10.0) |

**Key observation**: Incubation period inversely related to dose (log-normally distributed).

## Dose-Severity Correlations (Table 3, p.27)

Correlations weaken with stricter illness definitions:

| Outcome | Def I (n=120) | Def III (n=112) | Def IV (n=78) |
|---------|---------------|-----------------|---------------|
| Peak temperature r | 0.22 (0.04-0.39) | 0.18 (0.00-0.35) | 0.17 (-0.05-0.38) |
| Log time >103°F r | 0.13 (-0.05-0.31) | 0.09 (-0.09-0.28) | 0.03 (-0.20-0.25) |
| Symptom score r | 0.27 (0.09-0.43) | 0.24 (0.06-0.41) | 0.22 (0.00-0.42) |
| Symptom score r (adjusted for year) | 0.27 (-0.12-0.65) | 0.09 (-0.29-0.46) | -0.12 (-0.56-0.32) |

**[ASSISTANT-PROPOSED]** Weak dose-severity correlation; association with symptom score lost after adjusting for year. Definition used strongly affects conclusions.

## Incubation Period as Severity Predictor (Table 4, p.28)

Using Definition I:

| Incubation | N | Peak temp (°F) | Time >103°F | Symptom score | Abx (%) | Relapse (%) |
|------------|---|----------------|-------------|---------------|---------|-------------|
| <7 days | 41 | 104.2 (103.9-104.4) | 4.2 (3.2-5.5) | 9.0 (8.4-9.6) | 97.6% | 32.5% |
| 7-9 days | 47 | 103.9 (103.6-104.2) | 3.2 (2.3-4.5) | 8.2 (7.5-8.9) | 91.5% | 21.3% |
| ≥10 days | 31 | 103.6 (103.3-103.9) | 2.6 (1.8-3.8) | 6.9 (6.2-7.7) | 77.4% | 6.5% |

Shorter incubation → more severe disease, more likely to need treatment, more relapses.

## Data Quality Notes

### Strengths
- Individual-level reanalysis from original charts
- Multiple illness definitions compared
- Largest N for Maryland data
- Adjusts for confounders (age, race, military service, year)

### Limitations
- Retrospective chart review
- Temporal confounding (year) hard to separate from dose (all post-1965 received 10^5)
- Outpatient vs inpatient ascertainment differences

### Page References
- Methods and definitions: p.24
- Table 1 (attack rates): p.26
- Table 2 (incubation/severity by dose): p.26
- Table 3 (correlations): p.27
- Table 4 (incubation as predictor): p.28

## Cross-References
- Uses same underlying data as Hornick papers
- Authors include Hornick and Levine - authoritative re-analysis
- Cites Hornick 1970 as primary source

## Fit Role Assessment
**[ASSISTANT-PROPOSED]** **CORE** - Most detailed attack rate data with multiple illness definitions; essential for understanding how outcome definition affects dose-response estimates. Table 1 provides the definitive dose-response curve data.

## Key Extractions for Calibration

### Primary dose-response data (Definition I - most sensitive):
| Dose | N | N ill | Attack rate |
|------|---|-------|-------------|
| 10^3 | 13 | 0 | 0% |
| 10^5 | 200 | 83 | 41.5% |
| 10^7 | 27 | 13 | 48.1% |
| 10^8-10^9 | 25 | 24 | 96.0% |

### For comparison: Definition IV (strict, matches Hornick):
| Dose | N | N ill | Attack rate |
|------|---|-------|-------------|
| 10^3 | 13 | 0 | 0% |
| 10^5 | 200 | 52 | 26.0% |
| 10^7 | 27 | 8 | 29.6% |
| 10^8-10^9 | 25 | 18 | 72.0% |

**[OPEN]** Should we use Definition I (sensitive, captures more true infections) or Definition IV (strict, matches original Hornick fever definition) for calibration? This affects N50 estimate substantially.
