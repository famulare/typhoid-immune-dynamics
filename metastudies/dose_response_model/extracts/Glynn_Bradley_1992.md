# Glynn & Bradley 1992 - Dose-Severity Relationship in Salmonella Outbreaks

## Citation
Glynn JR, Bradley DJ. The relationship between infecting dose and severity of disease in reported outbreaks of salmonella infections. *Epidemiol Infect*. 1992;109:371-388. DOI: 10.1017/S0950268800050366

## Study Overview

**Study Type:** Meta-analysis of outbreak data (not a challenge study)

**Objective:** Investigate the relationship between infecting dose and severity of resulting disease for salmonella infections using proxy measures of dose from outbreak reports.

**Key Finding for S. typhi:** No evidence for a dose-severity relationship. Attack rate and incubation period are related to dose, but there is no evidence that they are related to severity (case fatality rate). This contrasts with non-typhoid salmonellas where dose-severity relationships were found.

**Organisms Studied:**
- *Salmonella typhi* (typhoid fever) - 68 epidemics
- *S. typhimurium* - 16 epidemics
- *S. enteritidis* - 11 epidemics
- *S. infantis* - 7 epidemics
- *S. thompson* - 8 epidemics
- *S. newport* - mentioned in within-epidemic analyses

## Methods

### Proxy Measures Used
**Dose proxies:**
- Attack rate (higher attack rate = higher dose)
- Incubation period (shorter incubation = higher dose)
- Amount of infected food consumed
- Type of vehicle (water vs food; water-borne assumed lower dose)

**Severity measures:**
- Case fatality rate (CFR) for typhoid
- Hospitalization rate for non-typhoid salmonellas

### Study Selection Criteria

**For typhoid (between-epidemic comparison):**
- Common source epidemics
- Pre-1945 for CFR analysis (no antibiotic therapy)
- At least 8 cases
- Sufficient information on CFR and attack rate or incubation period
- Excluded outbreaks among hospital patients

**For food-poisoning salmonellas:**
- CDC Salmonella Surveillance Reports 1964-1974
- Common source outbreaks
- More than 8 people involved
- Excluded hospital patients and mixed infections

### Statistical Analysis
- Unweighted linear regression (each epidemic equal weight)
- Logistic regression for circumscribed epidemics
- Models: Logit(outcome) = Constant + beta(proxy measure)

## Subject Characteristics

**Typhoid epidemics (n=69):**
- Mostly from Britain and United States
- 35 water-borne, remainder food-borne
- Pre-war circumscribed subgroup: n=31 (well-defined exposed population)
- Age, sex, immune status: insufficient information to control for these variables

**Food-poisoning salmonella epidemics (n=49):**
- All from USA (CDC reports)
- 1964-1974 period
- Most food-borne (one water-borne S. typhimurium epidemic with 14,000 cases)

## Outcomes with Definitions

### Typhoid Fever
| Outcome | Definition | Notes |
|---------|------------|-------|
| Case fatality rate (CFR) | Deaths / Total cases | Pre-1945 epidemics only |
| Attack rate | Cases / Exposed population | Proxy for dose |
| Incubation period | Time from exposure to illness onset | Median used; proxy for dose |

### Food-Poisoning Salmonellas
| Outcome | Definition | Notes |
|---------|------------|-------|
| Hospitalization rate | Hospitalized cases / Total cases | Time/culture dependent measure |
| Attack rate | Cases / Exposed population | Proxy for dose |
| Incubation period | Time from exposure to illness onset | Median in hours |

## Key Results for Salmonella typhi

### Attack Rate vs Incubation Period (p. 375-376, Fig. 2)
- Correlation coefficient = -0.55 (95% CI: -0.78 to -0.17)
- Regression coefficient = -4.01 (95% CI: -6.64 to -1.38)
- p < 0.001
- **Confirms expected relationship:** Higher attack rate associated with shorter incubation period

### Water vs Food-borne Epidemics (p. 377, Table 1)
| Parameter | Water-borne | Food-borne | P-value |
|-----------|-------------|------------|---------|
| Incubation (days), mean | 18.5 (n=8) | 13.2 (n=18) | 0.007 |
| Attack rate (%), geometric mean | 3.8 (n=34) | 16.7 (n=28) | <0.001 |
| CFR (%), mean (pre-1945) | 10.3 (n=32) | 11.4 (n=26) | 0.8 |

### Dose Proxies vs Case Fatality Rate
**Unweighted analysis (all epidemics):**
- Attack rate vs CFR: No significant correlation (p. 377, Fig. 3a)
- Incubation period vs CFR: No significant correlation (p. 377, Fig. 3b)

**Logistic regression (circumscribed pre-war epidemics):**
- Attack rate vs CFR: No association
- Incubation period vs CFR: Marginal positive association (LRS=4.4, 1 df, P=0.04)
  - [ASSISTANT-PROPOSED] This paradoxical finding (longer incubation associated with higher CFR) is opposite to what dose-severity would predict; authors note it was "only just significant"

### Key Conclusion for S. typhi (p. 381-382)
> "Overall, the comparison of data between typhoid epidemics provides no evidence of an association between dose (as measured by incubation period, attack rate or vehicle) and severity (as measured by the case fatality rate)."

This is consistent with Hornick volunteer studies: "In volunteer experiments with typhoid, vaccines gave protection against low but not high challenge doses, but once clinical disease occurred the severity of the disease and the number of relapses were not altered by vaccination."

## Dose-Outcome Data Tables

### Appendix 1: Typhoid Epidemics (p. 386-387)
Complete data for 69 epidemics. Selected examples showing range:

| Year | Place | Vehicle | Incubation (days) | Attack Rate % (n) | CFR % (n) |
|------|-------|---------|-------------------|-------------------|-----------|
| 1893 | Worthing | Water | - | 8.33 (1298/15579) | 12.48 (162/1298) |
| 1902 | Winchester | Oysters | 14 | 8.20 (10/122) | 40.00 (4/10) |
| 1914 | Hanford, CA | Spaghetti | 7 | 56.67 (85/150) | 3.23 (3/93) |
| 1916 | Helm, CA | Icecream | 6 | 95.83 (23/24) | 13.04 (3/23) |
| 1963 | Zermatt | Water | 17 | 4.37 (437/10000) | 0.69 (3/437) |

**Circumscribed pre-1945 outbreaks (marked with * in Appendix):** 31 epidemics identified with well-defined exposed populations (e.g., specific meals, defined milk rounds).

### Appendix 2: Non-Typhoid Salmonella Epidemics (p. 387-388)
Data from CDC Surveillance Reports 1964-1974.

**S. typhimurium (16 epidemics):**
| Year | Place | Vehicle | Attack Rate % | Hosp Rate % | Incubation (h) |
|------|-------|---------|---------------|-------------|----------------|
| 1965 | California | Water | 12.7 (14000/110000) | 0.5 | - |
| 1968 | N. Carolina | Icecream | 89.5 (17/19) | 0.0 | - |
| 1970 | Missouri | Icecream | 100.0 (11/11) | 0.0 | - |
| 1971 | Michigan | Smoked fish | 75.7 (28/37) | 10.7 | - |

**Logistic regression result (Table 2, p. 379):**
- Coefficient = 0.0406, SE = 0.0028, P < 0.001
- Model: Logit HR = -5.6 + 0.0406(AR)
- Deviance explained: 58.7%

**S. enteritidis (11 epidemics):**
- Coefficient = 0.0223, SE = 0.0063, P < 0.001
- Deviance explained: 10.7%
- Association lost after excluding one outlier with high hospitalization rate

**S. infantis (7 epidemics):**
- Coefficient = 0.0415, SE = 0.0072, P < 0.001
- Deviance explained: 43.6%

**S. thompson (8 epidemics):**
- Coefficient = -0.0069, SE = 0.0043, P = 0.1 (not significant)
- No dose-severity relationship found in between-epidemic analysis
- [OPEN] Authors note within-epidemic evidence exists for S. thompson (Tennessee outbreak, ref 45)

## Within-Epidemic Analyses (Single Outbreak Studies)

### S. typhi - Amount of Food (p. 374)
| Outbreak | Finding | P-value |
|----------|---------|---------|
| School picnic (ice cream) | CFR 3/17 (whole portion) vs 0/6 (half portion) | P = 0.5 (Fisher's) |
| Milk-borne (68 cases) | "Milk in tea/coffee only" had very mild attacks | Anecdotal |

### S. typhi - Incubation Period vs Outcome (p. 375)
Five outbreaks examined; incubation periods for fatal cases slightly shorter on average, but no significant differences at 5% level.

### Food-Poisoning Salmonellas - Incubation vs Severity (p. 375)
| Organism | Finding | Reference |
|----------|---------|-----------|
| S. newport (Sweden) | Shorter incubation = more severe illness (n=161) | [44] |
| S. thompson (Tennessee) | Earlier onset in hospitalized (n=51) vs others (n=72) | [45] |
| S. typhimurium (Italy) | Mean incubation 21h overall, 14h for 5 deaths (83 cases) | [26] |
| S. newport | No relationship (median 29h overall, 30h for severe) | [43] |

## Data Quality Notes

### Limitations Acknowledged (p. 379-381)
1. **Attack rate accuracy** (p. 379-380):
   - Depends on complete case ascertainment AND correct identification of exposed population
   - Water-borne outbreaks: population estimates approximate; some avoid exposure by boiling water
   - Most accurate for circumscribed epidemics (specific meals, camps)

2. **Incubation period** (p. 380):
   - Only estimable for point-source epidemics
   - Requires onset dates (not notification dates)
   - Late cases may be missed; secondary cases may contaminate data
   - Median used (epidemic curve approximately log-normal)

3. **Case fatality rate** (p. 380):
   - Deaths more completely identified than cases
   - Leads to CFR overestimation, varying by epidemic
   - Crude measure; numbers small in some epidemics

4. **Confounders** (p. 380):
   - Age, sex, immune status: insufficient data to control
   - Year of epidemic: controlled for in multiple regression; did not affect results

5. **Potential biases** (p. 380):
   - Water-borne: attack rates unduly low, CFR disproportionately high
   - Large epidemics often more inaccurate and biased than small ones

### Data Sources (p. 372-373)
- Bulletin of Hygiene / Abstracts of Hygiene (from 1926)
- The Lancet (1920-45)
- American Journal of Hygiene (1921-64)
- British Local Government Board Medical Officer's Reports (from 1900)
- CDC Salmonella Surveillance Reports (1964-76)
- Morbidity and Mortality Weekly Reports (from 1976)

## Cross-References

### Related Challenge Studies Cited
- **Hornick et al. 1970** (ref [2]): S. typhi volunteer studies showing dose-attack rate and dose-incubation relationships, but no dose-severity relationship
- **McCullough & Eisele 1951** (refs [3-6]): Non-typhoid salmonella volunteer experiments

### Related Outbreak Reviews
- Blaser & Newman 1982 (ref [22]): Review of human salmonellosis, infective dose
- Naylor 1983 (ref [23]): Incubation period features of typhoid outbreaks

### Specific Outbreaks of Interest
- **Zermatt 1963** (ref [19]): Large typhoid outbreak with detailed data
- **Mintz et al. 1991** (ref [21]): S. enteritidis outbreak with detailed dose-response by amount of Hollandaise sauce consumed

## Fit Role

**Fit Role: Support**

**Rationale:** [ASSISTANT-PROPOSED]
- This is a meta-analysis of outbreak data, not a controlled challenge study
- No direct dose measurements; relies entirely on proxy measures
- Primary value is the conclusion that for S. typhi, there is no evidence of dose-severity relationship, which is consistent with Hornick volunteer data
- Provides useful epidemiological context for interpreting dose-response relationships
- The finding that dose affects infection probability but not disease severity (given infection) is mechanistically important for modeling

**Specific Uses:**
1. Supports the assumption that typhoid disease severity is independent of infecting dose
2. Confirms dose-attack rate and dose-incubation relationships from outbreak data
3. Provides contrast with non-typhoid salmonellas where dose-severity relationships exist

**Limitations for Dose-Response Modeling:**
- No quantitative dose data
- Proxy measures have significant measurement error
- Cannot distinguish dose effects from strain virulence variation between outbreaks
- Pre-antibiotic era CFR data may not reflect modern clinical outcomes

## Open Questions

[OPEN] The marginal positive association between incubation period and CFR in circumscribed typhoid epidemics (P=0.04) runs counter to the dose-severity hypothesis. Authors suggest this is likely a chance finding, but alternative explanations should be considered.

[OPEN] The Appendix 1 typhoid epidemic data could potentially be used for more detailed analysis of attack rate distributions across outbreaks, but vehicle heterogeneity and data quality issues would need careful consideration.

[OPEN] Several within-epidemic studies cited (refs 44, 45, 21) may contain additional quantitative data on dose-response that could be extracted separately.
