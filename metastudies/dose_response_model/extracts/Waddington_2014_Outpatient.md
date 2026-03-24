# Waddington et al. 2014 - An Outpatient, Ambulant-Design, Controlled Human Infection Model

## Citation
Waddington CS, Darton TC, Jones C, Sherwood E, Sherwood C, et al. An Outpatient, Ambulant-Design, Controlled Human Infection Model Using Escalating Doses of Salmonella Typhi Challenge Delivered in Sodium Bicarbonate Solution. Clin Infect Dis. 2014;58(9):1230-1240.

## Study Overview
- **Primary question**: Develop safe outpatient typhoid CHIM model with dose escalation
- **Study type**: Prospective dose-escalation controlled human infection study
- **Institution**: Oxford Vaccine Group, University of Oxford
- **Population**: Healthy adults 18-60 years, typhoid-naive
- **Registration**: NCT01405521

**[CRITICAL]** First modern (Oxford) CHIM study - establishes that sodium bicarbonate delivery dramatically increases effective dose compared to Maryland milk delivery.

## Subject Characteristics

### Demographics (Table 1, p.1232)
- **N total**: 41 enrolled (21 at 10^3, 20 at 10^4); 40 per-protocol (1 withdrawn from 10^3 group)
- **Age**: Mean 29.4 years (range 20-50)
- **Sex**: 24 female (59%), 17 male (41%)
- **Ethnicity**: 38 White, 2 Asian, 1 Mixed

### Exclusion Criteria
- Prior typhoid vaccination or infection
- Travel to endemic area in past 6 months
- Antibiotic use in past 7 days
- Chronic illness, immunosuppression

## Challenge Protocol

### Strain & Delivery
- **Strain**: Quailes strain S. Typhi (same as Maryland studies)
- **Delivery**: "120 mL of sodium bicarbonate solution (NaHCO3 2.1 g in sterile water)" (p.1231)
- **Fasting**: "at least 90 minutes" before challenge
- **Rationale**: Bicarbonate neutralizes gastric acid, increasing viable organisms reaching intestine

### Dose Escalation Design (p.1231-1232)
- Adaptive dose-escalation algorithm (Figure 1, p.1231): protocol starts by challenging 5 participants at initial dose; if ≥2/5 develop infection, challenges further 5 for n=10, then expands cohort
- **Two dose levels tested** (Table 1, p.1233; Figure 2 CONSORT flowchart, p.1234):
  - 10^3 CFU: 21 enrolled, 1 withdrawn, 20 per-protocol
  - 10^4 CFU: 20 enrolled, 20 per-protocol
- Escalation halted after 10^4 because target attack rates were achieved at both dose levels
- Abstract: "Two dose levels (10^3 or 10^4 colony-forming units) were required to achieve the primary objective."

### Actual Assayed Challenge Doses (Table 1, p.1233)

Inocula freshly prepared and back-titrated by triplicate plating on tryptone soya agar (p.1232).

| Characteristic | Dose 1 (nominal 10³) | Dose 2 (nominal 10⁴) |
|---|---|---|
| **Target range** | 1-5 × 10³ CFU | 10-50 × 10³ CFU |
| **Actual dose, all participants, median [IQR]** | 1.34 × 10³ [0.98-1.69] | 19.8 × 10³ [18.8-21.6] |
| **Actual dose, diagnosed, median [IQR]** | 1.05 × 10³ [0.97-1.58] | 20.3 × 10³ [18.8-20.3] |
| **Actual dose, not diagnosed, median [IQR]** | 1.39 × 10³ [1.00-1.79] | 19.4 × 10³ [18.8-22.8] |

**Key point for modeling**: Actual median doses (1.34 × 10³ and 1.98 × 10⁴) are close to nominal (10³ and 10⁴) but with narrow IQRs, confirming good dosing precision. The ~2× difference between nominal and actual at the higher dose level is typical for CHIM inocula preparation.

## Illness Definition (p.1232)

**Typhoid Diagnosis** (composite endpoint):
1. Temperature ≥38°C sustained for ≥12 hours, OR
2. S. Typhi bacteremia (positive blood culture)

**[CRITICAL]** More sensitive fever threshold (38°C = 100.4°F) than Maryland (103°F = 39.4°C for strict definition)

### Antibiotic Treatment Trigger
- "participants commenced a 14-day course of ciprofloxacin (500 mg twice daily)" upon meeting diagnosis criteria (p.1232)

## Attack Rate Data (Table 2, p.1234)

**[KEY DATA - OXFORD BICARBONATE DELIVERY]**

| Dose | N Challenged | N Diagnosed | Attack Rate |
|------|--------------|-------------|-------------|
| 10^3 CFU | 20 | 11 | 55% |
| 10^4 CFU | 20 | 13 | 65% |

**By Diagnosis Criterion** (Table 2):
| Dose | Temp + BC confirmation | Temp only | BC + clinical signs | BC only | Not diagnosed |
|------|------------------------|-----------|---------------------|---------|---------------|
| 10^3 | 6 | 1 | 1 | 3 | 9 |
| 10^4 | 5 | 2 | 5 | 1 | 7 |

**[CRITICAL]** At 10^3-10^4 CFU in bicarbonate: 55-65% attack rate
Compare to Maryland: 0% at 10^3 in milk, 28% at 10^5 in milk
**This represents ~10,000-fold (4 log) difference in effective dose!**

## Incubation Periods (Table 2, p.1234)

### All participants (n=20 per dose):
| Parameter | 10^3 (mean ± SD, median [IQR]) | 10^4 (mean ± SD, median [IQR]) |
|-----------|-------------------------------|-------------------------------|
| Time to antibiotic initiation (days) | 11.65 ± 3.27 (n=20), 14 [8.75-14] | 9.9 ± 3.35 (n=20), 9 [4-14] |

### Diagnosed participants only:
| Parameter | 10^3 (mean ± SD (n), median [IQR]) | 10^4 (mean ± SD (n), median [IQR]) |
|-----------|-------------------------------------|-------------------------------------|
| Time to diagnosis (any criterion) | 9.73 ± 3.35 (11), 9 [6.5-13] | 7.69 ± 1.65 (13), 8 [6-9] |
| Time to clinical diagnosis | 7.57 ± 1.74 (7), 7.5 [6-8.75] | 7.57 ± 1.90 (7), 7 [6-8.5] |
| Time to microbiological diagnosis | 13.5 ± 0.58 (4), 13.5 [13-14] | 7.83 ± 1.72 (6), 8 [6.25-9] |

*Source: Table 2, p.1234. Text (p.1233) summarizes: "the median incubation period from challenge to diagnosis was 8 days (10³ CFU, 9 days [IQR, 6.5–13]; 10⁴ CFU, 8 days [IQR, 6–9])."*

**Correction note (2026-03-24)**: Prior extraction listed fabricated values for "Fever duration", "Peak temperature", and "Time to diagnosis" at 10³ (2.1, 39.0, 8.1 respectively). These do not appear in any table in the main paper or supplement. Replaced with verified Table 2 values.

## Clinical Findings (Table 3, p.1237)

| Parameter | 10^3 (n=20), median [IQR] | 10^4 (n=20), median [IQR] |
|-----------|---------------------------|---------------------------|
| Oral temp at diagnosis (°C) | 37.9 [37-38.6] | 37.8 [37.4-38.1] |
| Oral temp at baseline (°C) | 36.3 [36.1-36.4] | 36.3 [36.1-36.6] |
| Radial pulse at diagnosis (bpm) | 91 [80.5-103.5] | 95 [87-102] |
| Rash, n (%) | 4 (20) | 4 (20) |

### Severe Typhoid (Table 2):
| Criterion | 10^3 | 10^4 |
|-----------|------|------|
| Oral temp ≥40°C | 0 (0%) | 2 (10%) |
| Grade 3+ laboratory abnormality | 1 (5%) | 1 (5%) |
| All severe typhoid | 1 (5%) | 3 (15%) |

## Blood Culture & PCR Results (p.1235)

### Bacteremia Detection
- Blood cultures taken every 12 hours
- PCR also performed (results in Darton et al. 2017)

| Dose | Blood culture positive | Rate |
|------|------------------------|------|
| 10^3 | 10/20 | 50% |
| 10^4 | 11/20 | 55% |

Note: At 10^3, 10 of 11 (90.9%) diagnosed participants were bacteremic. At 10^4, 11 of 13 (84.6%) diagnosed participants were bacteremic.

**Blood culture quantification at diagnosis** (Table 4, p.1237):
- 10³: median 0.47 CFU/mL (from text p.1235; Table 4 diagnosed: 0.5 [0-1.2])
- 10⁴: median 1.10 CFU/mL (from text p.1235; Table 4 diagnosed: 1.1 [0.4-2.1])
- Positive cultures detected: median 1.7 per bacteremic participant at 10³, 2.6 at 10⁴ (Table 4)
- Prior to antibiotics, 6.2% (10³) and 9.6% (10⁴) had positive blood cultures (p.1235)

## Stool Shedding (p.1235, 1237)

| Dose | N shedding at any time | Rate |
|------|------------------------|------|
| 10^3 | 13/20 | 65% |
| 10^4 | **not directly reported** | **bounded below** |

**Source for 10³**: Text (p.1237 top): "Stool quantification data for participants challenged with the 10³ CFU dose demonstrated that the median number of bacteria excreted was 32 (IQR, 18–41) CFU/g feces (n = 13)."

**Bounding W-I-4 (shedding at 10⁴)**: The paper does NOT report a participant-level shedding count at 10⁴. However:
- Text (p.1237): "Salmonella Typhi was subsequently shed by 18 of 24 (75%) of typhoid-diagnosed participants" (combined across both doses)
- Text (p.1237): "One otherwise asymptomatic participant (10³ CFU dose) had S. Typhi cultured from stool" — so 1 undiagnosed shedder at 10³
- Table 2: "≥38°C OR bacteremia OR positive stool" = 14/20 at 10⁴ = same as "≥38°C OR bacteremia" = 14/20. Therefore stool positivity adds zero undiagnosed participants at 10⁴ — all shedders at 10⁴ also had fever or bacteremia.
- If all 11 diagnosed at 10³ shed + 2 undiagnosed = 13 (matching n=13), then 18 - 11 = 7 diagnosed at 10⁴ shed, + 0 undiagnosed = **7/20 (35%) minimum**. But this assumes all 11 diagnosed at 10³ shed, which may overcount.
- Upper bound: ≤14/20 (70%) since that is the "fever OR bacteremia OR stool" endpoint.
- **Best estimate for modeling: use combined shedding rate 19/40 (47.5%) and note dose-specific is unavailable, OR use the Table 2 "≥38°C OR bacteremia OR positive stool" endpoint (14/20 = 70%) as an upper bound.**

## Alternative Fever Thresholds (Table 2, p.1234)

All participants (n=20 per dose):

| Endpoint | 10³ n (%) | 10⁴ n (%) |
|----------|-----------|-----------|
| ≥37.5°C (any duration) | 11 (52) | 14 (70) |
| ≥38.0°C (any duration) | 10 (48) | 13 (65) |
| ≥38.5°C (any duration) | 8 (38) | 11 (55) |
| ≥38.0°C OR bacteremia | 12 (57) | 14 (70) |
| ≥38.0°C AND bacteremia | 7 (33) | 10 (50) |
| ≥38.0°C with subsequent bacteremia | 5 (24) | 8 (40) |
| ≥38.0°C + subsequent bacteremia OR positive stool | 7 (33) | 10 (50) |
| ≥38.0°C OR bacteremia OR positive stool | 14 (67) | 14 (70) |

*These are critical for φ calibration (definition sensitivity) and cross-study comparison.*

## Serology Data

### Pre-challenge (p.1236)
- Anti-H and anti-Vi antibody measured at baseline
- "Low levels of anti-H and anti-Vi antibody were measured in all participants at baseline" (p.1236)
- "did not correlate with subsequent risk of infection" (p.1236)

### Post-challenge Antibody Responses (Supplementary Table 2)

**Anti-H (flagellin) IgG — typhoid-diagnosed participants:**

| Group | Baseline EU/mL (95% CI) | Day 14 EU/mL (95% CI) | Mean fold rise (% >4-fold) | Day 28 EU/mL (95% CI) | Mean fold rise (% >4-fold) |
|-------|-------------------------|------------------------|----------------------------|------------------------|----------------------------|
| 10³ TD | 27.8 (15.5-49.8) | 117.8 (53.9-257.6) | 4.2 (45%) | 150.6 (82.6-274.7) | 5.4 (64%) |
| 10⁴ TD | 38.8 (24.9-60.5) | 335.3 (128.4-875.4) | 8.6 (77%) | 286.5 (123.2-665.8) | 7.4 (77%) |
| No TD | 40.1 (29.4-54.6) | 47.5 (35.2-64.1) | 1.2 (0%) | 53.7 (37.9-75.9) | 1.3 (0%) |

**Anti-Vi (capsular) IgG — NO response in any group:**

| Group | Baseline EU/mL (95% CI) | Day 14 EU/mL (95% CI) | Mean fold rise (% >4-fold) |
|-------|-------------------------|------------------------|----------------------------|
| 10³ TD | 3.0 (1.6-5.6) | 3.2 (1.4-7.1) | 1.1 (0%) |
| 10⁴ TD | 4.2 (2.0-8.8) | 4.7 (2.5-9.0) | 1.1 (0%) |
| No TD | 6.9 (3.2-14.8) | 5.2 (2.7-9.7) | 0.8 (0%) |

**Key finding**: "anti-Vi levels remained unchanged throughout" (p.1236). No Vi seroconversion after challenge in any participant. This confirms that natural infection with Quailes strain does NOT induce anti-Vi antibody, consistent with the Maryland findings and supporting the model assumption that anti-Vi is vaccine-derived only.

**Anti-LPS (somatic O) IgG — strong response in diagnosed only:**

| Group | Baseline EU/mL (95% CI) | Day 14 EU/mL (95% CI) | Mean fold rise (% >4-fold) |
|-------|-------------------------|------------------------|----------------------------|
| 10³ TD | 33.7 (13.6-83.7) | 235.4 (65.8-842.4) | 7.0 (55%) |
| 10⁴ TD | 52.4 (27.3-100.4) | 590.0 (277.6-1254.0) | 11.3 (77%) |
| No TD | 38.9 (24.9-60.7) | 37.4 (23.4-59.7) | 1.0 (0%) |

*Source: Supplementary Table 2. TD = typhoid diagnosed; No TD = cumulative non-diagnosed across both doses.*

## Supplementary Symptom Data (Supplementary Table 1)

Solicited symptoms by dose and diagnosis status (n affected / n in group):

| Symptom | 10³ TD (n=11) | 10³ No TD (n=9) | 10⁴ TD (n=13) | 10⁴ No TD (n=7) | DOR (95% CI) |
|---------|---------------|-----------------|---------------|-----------------|--------------|
| Headache | 11/11 (100%) | 6/9 (67%) | 13/13 (100%) | 5/7 (71%) | 23.4 (1.2-460.7) |
| Generally unwell | 10/11 (91%) | 4/9 (44%) | 13/13 (100%) | 2/7 (29%) | 38.3 (4.1-361.3) |
| Loss of appetite | 9/11 (82%) | 4/9 (44%) | 13/13 (100%) | 1/7 (14%) | 19.8 (3.2-121.5) |
| Abdominal pain | 8/11 (73%) | 4/9 (44%) | 12/13 (92%) | 2/7 (29%) | 8.3 (1.9-36.4) |
| Nausea/vomiting | 8/11 (73%) | 4/9 (44%) | 12/13 (92%) | 1/7 (14%) | 11.0 (2.4-49.6) |
| Myalgia | 9/11 (82%) | 2/9 (22%) | 12/13 (92%) | 5/7 (71%) | 9.0 (1.9-42.9) |
| Arthralgia | 7/11 (64%) | 1/9 (11%) | 12/13 (92%) | 0/7 (0%) | 57.0 (6.0-541.5) |
| Cough | 8/11 (73%) | 2/9 (22%) | 8/13 (62%) | 2/7 (29%) | 9.0 (2.1-38.8) |
| Diarrhoea | 3/11 (27%) | 2/9 (22%) | 6/13 (46%) | 1/7 (14%) | 2.6 (0.6-11.7) |
| Constipation | 8/11 (73%) | 2/9 (22%) | 8/13 (62%) | 2/7 (29%) | 5.5 (1.3-22.9) |

*DOR = Diagnostic Odds Ratio (pooled across doses, correction applied per Glas et al. 2003). Arthralgia was the most discriminatory symptom (DOR 57.0).*

## Safety (p.1237)

- No serious adverse events
- All participants cleared infection with treatment
- No chronic carriage
- "Challenge model was well tolerated"

## Data Quality Notes

### Strengths
- Modern prospective study with systematic sampling
- Detailed individual-level data available
- Blood cultures q12h - high ascertainment
- Both bacteremia and fever endpoints captured
- PCR complement to blood culture

### Limitations
- Small sample sizes at each dose (n=20 per dose)
- Only two dose levels tested (10^3 and 10^4)
- Different fever threshold than Maryland (38°C vs 39.4°C)
- Bicarbonate delivery not directly comparable to milk

### Page References
- Study design: p.1231-1232
- Table 2 (attack rates): p.1234
- Table 3 (clinical outcomes): p.1235
- Bacteremia/shedding data: p.1235
- Safety: p.1237

### Correction History
- **2026-03-19**: Major corrections per `notes/verification/waddington_10e5_source_investigation.md`. Removed entirely fabricated 10^5 dose group (n=5, 100% attack rate) — AI confabulation from dose-escalation design diagram. Corrected 10^4 denominator from n=16 to n=20. Corrected 10^4 attack rate from 10/16 (63%) to 13/20 (65%). Corrected bacteremia counts from Table 2 diagnostic breakdown. Restructured diagnostic criterion table to match paper's actual Table 2 categories. Several values marked [NEEDS RE-EXTRACTION] pending verification against primary source.

## Cross-References
- Same strain (Quailes) as Maryland studies
- Darton et al. 2016 - detailed immunology from this cohort
- Darton et al. 2017 - PCR optimization using these samples
- Gibani et al. 2019 - later studies using same model

## Fit Role Assessment
**[ASSISTANT-PROPOSED]** **CORE** - Establishes modern CHIM dose-response with bicarbonate delivery. Critical for understanding delivery medium effect on N50. Provides both bacteremia and fever endpoints.

## Key Extractions for Calibration

### Primary dose-response data (fever OR bacteremia):
| Nominal Dose | Actual Median Dose | N | N diagnosed | Attack rate |
|--------------|--------------------|---|-------------|-------------|
| 10^3 | 1.34 × 10³ | 20 | 11 | 55% |
| 10^4 | 1.98 × 10⁴ | 20 | 13 | 65% |

### Bacteremia-only endpoint:
| Nominal Dose | Actual Median Dose | N | N bacteremic | Rate |
|--------------|--------------------|---|--------------|------|
| 10^3 | 1.34 × 10³ | 20 | 10 | 50% |
| 10^4 | 1.98 × 10⁴ | 20 | 11 | 55% |

### Stool shedding (infection proxy):
| Nominal Dose | Actual Median Dose | N | N shedding | Rate | Notes |
|--------------|--------------------|---|------------|------|-------|
| 10^3 | 1.34 × 10³ | 20 | 13 | 65% | From stool quantification text (n=13) |
| 10^4 | 1.98 × 10⁴ | 20 | 7-14 | 35-70% | **NOT DIRECTLY REPORTED.** Lower bound: 18 total diagnosed shedders - 11 at 10³ = 7. Upper bound: Table 2 "fever OR bacteremia OR stool" = 14. |

**[OPEN]** Delivery medium (bicarbonate vs milk) requires explicit modeling - cannot directly pool Oxford and Maryland data without accounting for ~4 log effective dose difference.
