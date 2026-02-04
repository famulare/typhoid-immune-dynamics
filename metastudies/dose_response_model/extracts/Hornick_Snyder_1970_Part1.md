# Hornick et al. 1970 - Typhoid Fever: Pathogenesis and Immunologic Control (Part 1)

## Citation
Hornick RB, Greisman SE, Woodward TE, DuPont HL, Dawkins AT, Snyder MJ. Typhoid Fever: Pathogenesis and Immunologic Control (First of Two Parts). N Engl J Med. 1970 Sep 24;283(13):686-691.

**Note**: File labeled as 2008 but is actually 1970 Part 1

## Study Overview
- **Primary question**: Comprehensive review of typhoid pathogenesis from volunteer studies
- **Study type**: Review of Maryland CHIM program with detailed pathogenesis data
- **Institution**: Division of Infectious Diseases, University of Maryland School of Medicine

## Dose-Response Data (Table 1, p.687)

**CRITICAL - Most complete dose-response table with incubation periods**

| Dose | N Volunteers | N with Disease | Attack Rate | Incubation (Median) | Incubation (Range) |
|------|--------------|----------------|-------------|---------------------|-------------------|
| 10^9 | 42 | 40 (95%) | 95% | 5 days | 3-32 days |
| 10^8 | 9 | 8 (89%) | 89% | - | - |
| 10^7 | 32 | 16 (50%) | 50% | 7.5 days | 4-56 days |
| 10^5 | 116 | 32 (28%) | 28% | 9 days | 6-33 days |
| 10^3 | 14 | 0 (-) | 0% | - | - |

**Key**: Disease = "temperature of 103°F or higher by mouth for at least 24 to 36 hours" (p.687)

## Strain Virulence Comparison (Table 2, p.690)

**CRITICAL - Shows INFECTION vs DISEASE distinction and Vi antigen effect**

Challenge dose: ~10^7 organisms

### Vi-containing strains:
| Strain | Disease | Infection‡ | No Infection§ |
|--------|---------|------------|---------------|
| Quailes | 16 of 30 | 12 of 30 | 2 of 30 |
| Zermatt | 6 of 11 | 4 of 11 | 1 of 11 |
| Ty2V | 2 of 6 | 3 of 6 | 1 of 6 |
| **Totals (Vi)** | **24 of 47 (51%)** | **19 of 47 (40%)** | **4 of 47 (9%)** |

### Non-Vi strains:
| Strain | Disease | Infection‡ | No Infection§ |
|--------|---------|------------|---------------|
| 0-901 | 6 of 20 | 6 of 20 | 8 of 20 |
| Ty2W | 4 of 19 | 10 of 19 | 5 of 19 |
| **Totals (non-Vi)** | **10 of 39 (26%)** | **16 of 39 (41%)** | **13 of 39 (33%)** |

*Mouse virulence: Vi strains LD50 ~10^2-10^4; non-Vi strains LD50 ~10^5-10^6*

**Definitions**:
- †Disease: Temperature ≥103°F for >36 hr & treatment with antibiotic
- ‡Infection: Low-grade fever or significant serological response, or positive blood culture or excretion of S. typhosa in stools for >5 days & no specific therapy
- §No Infection: No clinical, cultural or serologic evidence

**[CRITICAL]** This table separates DISEASE from INFECTION outcomes - essential for model that distinguishes infection from fever!

### Vi Antigen Effect
- "Disease rates were significantly higher in volunteers who ingested Vi-containing strains than with non-Vi strains: 51 versus 26 per cent (p less than 0.05)" (p.690)

## Clinical Description

### Disease Definition (p.687)
- "A standard protocol was established to estimate vaccine effect. Thus, antibiotic treatment was begun when the temperature reached 103°F or higher, by mouth, and persisted at that level for 24 to 36 hours."

### Treatment
- "250 men received antibiotic treatment (215 with chloramphenicol) for typhoid fever"
- "Two cases of moderate hemolytic anemia occurred before chloramphenicol was given"
- "No permanent typhoid carriers"

## Pathogenesis Observations

### Intestinal Phase (p.688)
- "Stools yielded S. typhosa within 24 hours after ingestion and were often positive during the incubation period"
- "In volunteers given 10^5 to 10^7 organisms, the presence of S. typhi in the stools for the first five days after challenge did not necessarily mean that they would become ill"
- Initial stool shedding "is not a definite indication of failing host defense"

### Asymptomatic Bacteremia (p.690)
- "In two of 64 patients this condition was documented"
- "In one volunteer, bacteremia occurred for 10 days; there was no overt clinical disease and no fever"

## Strain Information

### Quailes Strain
- "Obtained from a carrier"
- "Well studied, and its virulence in relation to other classic typhoid strains has been defined"
- Contains Vi antigen
- Mouse LD50: 2.8 × 10^2

### Other Strains Tested
- **Zermatt**: Isolated from 1963 epidemic; contains Vi; mouse LD50 3.0 × 10^4
- **Ty2V**: Classic strain; contains Vi; mouse LD50 3.0 × 10^6
- **0-901**: From Russian village 1915; lacks Vi; mouse LD50 3.11 × 10^5
- **Ty2W**: Derived from Ty2V; lacks Vi; mouse LD50 1.9 × 10^5

## Cross-References
- Part 1 of 2-part NEJM review
- Part 2 (same issue/following week) contains vaccine efficacy data
- Data overlap with earlier Hornick publications but this has the cleanest tables

## Fit Role Assessment
**[ASSISTANT-PROPOSED]** **CORE** - Contains the most complete dose-response data AND the critical infection vs disease distinction needed for model. Table 2 is unique in showing this separation.

## Key Extractions for Calibration

### For Infection Model (Quailes strain, 10^7 dose):
- 30 challenged
- 28 infected (disease + infection categories) = 93%
- 2 no infection = 7%

### For Fever|Infection Model (Quailes strain, 10^7 dose):
- 28 infected
- 16 disease = 57% of infected developed fever

**[OPEN]** How to reconcile Table 1 (32 challenged at 10^7, 16 disease) with Table 2 (30 challenged at 10^7, 16 disease)? Small discrepancy likely due to different analysis subsets.
