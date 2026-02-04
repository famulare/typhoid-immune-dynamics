# Cross-Cutting Observations: Typhoid Challenge Literature

This document captures patterns, overlaps, and insights that span multiple papers in our extraction corpus.

---

## Cohort Overlaps

Papers that share subjects or reanalyze the same underlying data:

| Cohort | Papers | Notes |
|--------|--------|-------|
| Maryland volunteers (1960s-70s) | Hornick 1966, 1967, Appraisal, Snyder 1970; Dupont 1971; Gilman 1977; Levine 1976 | ~1,700 total volunteers over program |
| Maryland reanalyses | Glynn 1995; Woodward 1980; Waddington 2014 (review) | Secondary analyses of same underlying data |
| Oxford CHIM (2011-2017) | Waddington 2014; Darton 2012, 2016, 2017; Gibani 2019, 2020 | 430+ participants across studies |
| VAST trial | Jin 2017; Dahora 2019 | Same 103 participants |
| M01ZH09/Ty21a trial | Darton 2016; Juel 2018 | Same 91 participants |

### Known Relationships
- **Darton 2012** is preliminary data from **Waddington 2014** cohort
- **Darton 2017** reanalyzes **Waddington 2014** with PCR
- **Gibani 2019** pools 6 Oxford studies including Waddington, Darton, Jin
- **Glynn 1995** reanalyzes Hornick Maryland data
- **Woodward 1980** summarizes entire Maryland program

---

## Strain Lineage

### Quailes Strain
- **Origin**: Isolated 1958 from chronic carrier (female, gallbladder source) at University of Maryland Hospital following cholecystectomy
- **Studies using**: ALL Maryland studies; ALL Oxford studies
- **Key characteristics**: Wild-type, Vi-positive, virulent; stored at -70°C; same lineage used for >60 years
- **Genotype**: 3.1.0 (per Gibani 2019)

### Ty2 Strain
- **Origin**: Not used in challenge studies in this corpus
- **Studies using**: N/A (vaccine strain parent only)

### Other Strains
- **Vi-negative strains**: Tested in Hornick 1970 Part 2; lower attack rate (26% vs 51% at 10^7)
- **S. Paratyphi A NVGH308**: Used in Gibani 2019, 2020 for comparison

### Strain Comparability
- **[ASSISTANT-PROPOSED]** All human challenge data in this corpus uses Quailes strain; directly comparable
- **[OPEN]** Vi-positive vs Vi-negative comparison in Hornick 1970 suggests virulence difference

---

## Methodological Evolution

### Maryland Era (1960s-1970s)
- **Challenge setting**: Inpatient; Maryland House of Correction (prison volunteers)
- **Dose preparation**: Frozen milk specimens at -70°C; thawed for challenge
- **Delivery medium**: Skimmed milk (45 mL or 30 cc)
- **Buffering**: Not consistently documented; some studies mention NaHCO3
- **Outcome ascertainment**: Daily blood cultures days 4-10; all stools cultured; oral temp monitoring
- **Typical subject population**: Adult male inmates; median age not reported; mostly non-veterans
- **Disease threshold**: Temp >103°F for 24-36 hours + clinical symptoms → chloramphenicol treatment

### Oxford Era (2010s-present)
- **Challenge setting**: Outpatient with daily visits; hospitalization if diagnosis met
- **Dose preparation**: Fresh inoculum preparation from frozen stocks
- **Delivery medium**: Sodium bicarbonate solution (120 mL, 0.53g NaHCO3)
- **Buffering**: Standardized NaHCO3 pre-treatment
- **Fasting**: 90-minute fast before challenge
- **Outcome ascertainment**: Daily blood cultures; PCR available; standardized symptom scoring
- **Typical subject population**: Healthy UK adults 18-60; typhoid-naive; predominantly white
- **Disease threshold**: Temp ≥38°C for ≥12 hours OR positive blood culture

### Key Differences
1. **Delivery medium**: Milk (Maryland) vs bicarbonate solution (Oxford)
2. **Disease definition**: >103°F/36hr (Maryland) vs ≥38°C/12hr OR bacteremia (Oxford)
3. **Setting**: Inpatient prison (Maryland) vs outpatient community (Oxford)
4. **Diagnostics**: Culture only (Maryland) vs culture + PCR (Oxford)
5. **Dose range**: 10^3-10^9 (Maryland) vs 10^3-10^4 (Oxford)

---

## Outcome Definition Drift

How "typhoid fever" and related terms changed over time:

| Era/Paper | Term | Definition | Notes |
|-----------|------|------------|-------|
| Hornick 1966 | "disease" / "ill" | Temp >103°F for 24-36hr + clinical symptoms | Strict; requires sustained high fever |
| Hornick 1970 | "typhoid fever" | Same as 1966; chloramphenicol trigger | Disease, not just infection |
| Waddington 2014 (review) | compared | 26% (strict) vs 41.5% (relaxed) at 10^5 | Different definitions yield different rates |
| Oxford 2014+ | "typhoid diagnosis" (TD) | Temp ≥38°C for ≥12hr OR blood culture+ | More inclusive; lower threshold |
| Darton 2017 | "infection" | TD OR PCR+ | Captures subclinical infection |

**[OPEN]** Attack rate differences between Maryland (~28% at 10^5) and Oxford (~60% at 10^4) may partly reflect definition differences, not just dose.

---

## Immunity Measurement Evolution

### Assays Used by Era

| Era | Assays | Notes |
|-----|--------|-------|
| 1960s-70s | O agglutinins, H agglutinins | Not protective per Hornick 1966 |
| 1970s | H antibody titers | Protective in Gilman 1977 (≥1:20 threshold) |
| 2010s-present | Anti-Vi IgG, IgA, IgG subclasses, SBA | Anti-Vi IgG protective per Darton 2016 |

### Anti-Vi IgG
- **When first measured**: Modern Oxford studies (2014+)
- **Key finding**: 1 log increase → 71% reduction in hazard (Darton 2016)
- **Baseline prevalence**: ~29% of "naive" UK subjects have detectable anti-Vi (Darton 2016)

### Bactericidal Assays (SBA)
- **Papers reporting**: Juel 2018
- **Key finding**: Higher SBA associated with delayed onset, lower bacterial burden, reduced severity, but NOT protection against diagnosis
- **[OPEN]** SBA may measure disease-modifying vs protection from infection

### H Agglutinins
- **Papers reporting**: Gilman 1977
- **Key finding**: Pre-existing H antibody ≥1:20 protective (24% vs 61% attack rate, P=0.02)

---

## Dose Reporting Conventions

### Units
- **CFU (colony-forming units)**: Standard
- **"Organisms"**: Used interchangeably with CFU in older papers
- **Log notation**: 10^5, 10^7, etc.

### Dose Measurement Methods
- **Maryland**: Plate counts verified after preparation; frozen milk aliquots
- **Oxford**: Fresh preparation; ranges reported (e.g., "1-5 × 10^4")

### Calibration Issues
- **[OPEN]** Oxford reports dose ranges; exact CFU per subject varies within arm
- **[ASSISTANT-PROPOSED]** Maryland doses more precise due to verification; Oxford doses are target ranges

---

## Delivery Medium Effects

### Milk
- **Studies using**: ALL Maryland studies (Hornick, Dupont, Gilman, Levine)
- **Rationale**: Masks taste; may provide some gastric acid buffering

### Bicarbonate Solution
- **Studies using**: ALL Oxford studies (Waddington, Darton, Jin, Gibani)
- **Rationale**: Standardized gastric acid neutralization
- **NaHCO3 concentrations**: 0.53g in 120mL (Darton 2017)

### Evidence for Medium Effects
- **[OPEN]** No direct comparison in this corpus
- **[ASSISTANT-PROPOSED]** Bicarbonate may increase effective dose by reducing gastric acid killing
- Attack rate differences (Maryland ~28% at 10^5 vs Oxford ~60% at 10^4) could partly reflect medium effects

---

## Recurring Data Gaps

What's consistently missing across papers:

| Missing Element | Frequency | Impact on Analysis |
|-----------------|-----------|-------------------|
| Pre-challenge serology | High (esp. Maryland) | Requires latent CoP inference |
| Individual-level data | High | Limits hierarchical modeling |
| Exact dose per subject | Moderate (Oxford ranges) | Dose uncertainty in likelihood |
| Age distributions | High (esp. Maryland) | Cannot stratify by age |
| Prior typhoid history details | Moderate | Uncertainty in naive vs immune |
| Fasting compliance | High | Unquantified variability |

---

## Synthesis Notes

### Core Dose-Response Shape
- ID50 consistently estimated at ~10^7 CFU (Maryland studies)
- Attack rate at 10^5: ~28-40% (Maryland) vs ~60-77% at 10^4 (Oxford)
- Maximum attack rate: ~95% at 10^9 (Maryland)
- Minimum detectable dose: 10^3 shows 0% (Maryland) to low % (Oxford)

### Key Uncertainties for Modeling
1. **Outcome definition heterogeneity**: Maryland "disease" vs Oxford "diagnosis"
2. **Delivery medium effects**: Milk vs bicarbonate
3. **Baseline immunity**: Largely unmeasured in Maryland; variable in Oxford
4. **Cohort overlap**: Must avoid double-counting subjects

### Promising Data for Immunity Parameter (γ)
- Darton 2016: anti-Vi IgG vs attack rate
- Gibani 2020: prior exposure vs rechallenge attack rate
- Gilman 1977: H antibody vs attack rate
- Dupont 1971: prior infection vs attack rate
- Military vs non-veteran: Woodward 1980 (20% vs 48% at 10^5)

