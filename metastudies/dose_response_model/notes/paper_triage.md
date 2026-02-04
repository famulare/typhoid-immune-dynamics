# Paper Triage: Dose-Response Extraction

## Triage Criteria

### Inclusion Requirements (must have at least one)
- Oral challenge dose of S. Typhi (exact or range)
- Group-level denominators (N challenged per dose arm)
- At least one relevant outcome (infection or fever proxy)

### Fit Role Categories
- **Core**: Multi-dose studies with clear dose-stratified outcomes; best for estimating dose-response shape (N₅₀, α)
- **Support**: Single-dose or limited dose range but good immunity characterization; useful for γ estimation
- **Exclude**: No original dose-outcome data (pure reviews, duplicates, or data already captured elsewhere)

---

## Triage Summary Table

| # | Paper | Status | Dose Info | Multi-Dose? | Pre-Challenge Immunity? | Key Outcomes | Fit Role | Notes |
|---|-------|--------|-----------|-------------|------------------------|--------------|----------|-------|
| 1 | Dahora et al. 2019 | include | ~10^4 CFU | N | Y (Vi IgG, IgA) | 37% attack rate (vaccinees) | Support | Secondary immunological analysis of VAST trial; no unvacc controls |
| 2 | Darton et al. 2012 | exclude | 10^3-10^4 | Y | N | TD in 40 subjects | Exclude | Conference abstract; data superseded by Waddington 2014 |
| 3 | Darton et al. 2016 | include | ~1.8×10^4 | N | Y (anti-Vi IgG) | 67% (20/30) control attack | Support | Vaccine trial; 30 placebo controls valuable; individual data in supplement |
| 4 | Darton et al. 2017 | include | 10^3-10^4 | Y | N | 59% TD, 73% with PCR | Support | Same cohort as Waddington 2014; key for infection vs fever distinction |
| 5 | Dupont et al. 1971 | include | 10^5 | N | partial | 23% (prior inf) vs 30% (naive) | Support | Prior infection immunity data; single dose |
| 6 | Gibani et al. 2019 | include | 1-5×10^4 | N | Y (anti-Vi) | 56% TD rate pooled | Support | Pooled analysis of 6 Oxford studies; shedding focus |
| 7 | Gibani et al. 2020 | include | 1-5×10^4 | N | Y (prior challenge) | 63% naive, 26-68% rechallenge | Support | Re-challenge study; key for immunity modeling |
| 8 | Gilman et al. 1977 | include | 10^5 | N | partial (H antibody) | 48% (64 controls) | Support | Good control data; H antibody protective |
| 9 | Glynn & Bradley 1992 | include | outbreak proxies | N/A | N | CFR analysis | Support | Outbreak meta-analysis; no dose-severity relationship for typhoid |
| 10 | Glynn et al. 1995 | include | 10^3-10^9 | Y | N | reanalysis of Hornick | Support | Reanalysis of Maryland data; adds incubation period analysis |
| 11 | Hornick et al. 1967 | include | 10^5-10^9 | Y | N | vaccine efficacy | Support | Overlaps with Hornick 1966 |
| 12 | Hornick (Appraisal) | include | 10^5-10^7 | Y | N | vaccine efficacy | Support | Overlaps with other Hornick papers |
| 13 | Hornick 1966 | include | 10^3-10^9 | Y | N | 0%-95% by dose | **Core** | Multi-dose dose-response curve; foundational |
| 14 | Hornick & Snyder 1970 | include | 10^3-10^9 | Y | N | infection vs disease distinction | **Core** | Part 1+2; foundational; distinguishes infection from fever |
| 15 | Jin et al. 2017 | include | 1-5×10^4 | N | Y (postvax serology) | 77% (24/31) control | Support | VAST trial primary report; well-characterized single dose |
| 16 | Juel et al. 2018 | include | 10^4 | N | Y (SBA titers) | SBA correlates | Support | Antibody-severity correlations; single dose |
| 17 | Levine et al. 1976 | include | 10^5 | N | partial | ~40% (97 controls pooled) | Support | Vaccine trial; good control group data |
| 18 | Waddington 2014 (review) | exclude | summary | Y | N | synthesized | Exclude | Review paper; no new primary data |
| 19 | Waddington 2014 (outpatient) | include | 10^3-10^4 | Y | N | dose-escalation data | **Core** | Oxford dose-escalation study; modern multi-dose |
| 20 | Woodward 1980 | include | 10^3-10^9 | Y | partial (military) | 40% at 10^5 (N>300) | Support | Review/summary; military stratification unique |

---

## Detailed Triage Notes

### Tier 1 Candidates (Multi-Dose, Core)

1. **Hornick 1966** - Original Maryland dose-response data spanning 10^3 to 10^9 CFU
2. **Hornick & Snyder 1970** - Comprehensive Part 1+2; infection vs disease distinction
3. **Waddington 2014 (outpatient)** - Modern Oxford dose-escalation (10^3 to 10^4)

### Tier 2 Candidates (Support)

**Immunity-informative:**
- Darton 2016 - anti-Vi IgG correlates; individual-level data available
- Gibani 2020 - rechallenge data; prior exposure effects
- Gilman 1977 - H antibody protective in controls
- Dupont 1971 - prior infection immunity
- Jin 2017 - well-characterized single-dose vaccine trial
- Juel 2018 - bactericidal antibody correlates

**Attack rate validation:**
- Levine 1976 - 97 pooled controls at 10^5
- Glynn 1995 - reanalysis with incubation periods

**Methodological context:**
- Darton 2017 - infection vs fever; PCR sensitivity
- Gibani 2019 - shedding outcomes; pooled Oxford data
- Glynn & Bradley 1992 - outbreak data; no dose-severity relationship
- Woodward 1980 - military service stratification

### Exclusions

1. **Darton 2012** - Conference abstract superseded by Waddington 2014
3. **Waddington 2014 (review)** - Pure review; no new primary data

---

## Review Order for Phase 3

Based on triage, papers will be reviewed in this order:

1. **Core (multi-dose) papers**: Hornick 1966, Hornick & Snyder 1970, Waddington 2014 (outpatient)
2. **Support papers with immunity data**: Darton 2016, Jin 2017, Gibani 2020, Gilman 1977, Dupont 1971, Juel 2018
3. **Attack rate validation papers**: Levine 1976, Glynn 1995, Woodward 1980
4. **Methodological context**: Darton 2017, Gibani 2019, Glynn & Bradley 1992, Dahora 2019
5. **Exclusions to confirm**: Darton 2012, Waddington 2014 (review)

---

## Cohort Overlap Summary

### Maryland Studies (1960s-1970s) - Likely Overlapping Subjects
- Hornick 1966, 1967, Appraisal, Snyder 1970
- Dupont 1971, Gilman 1977, Levine 1976
- Glynn 1995, Woodward 1980 (reanalyses)

### Oxford Studies (2010s) - Known Overlapping Subjects
- Waddington 2014 (outpatient) - base cohort
- Darton 2012, 2016, 2017 - overlapping/extensions
- Gibani 2019, 2020 - pooled analyses
- Jin 2017, Dahora 2019 - VAST trial
- Juel 2018 - M01ZH09/Ty21a trial

