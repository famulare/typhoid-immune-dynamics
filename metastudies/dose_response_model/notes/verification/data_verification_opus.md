# Data Extraction Verification Report

**Verifier**: Claude (Opus 4.6, independent verification agent)
**Date**: 2026-03-17
**Documents verified**: Joint Inference Plan Section 4 (Data Catalog) and Section 11 (Summary of All Observations), against 11 source extraction documents.

---

## 1. Executive Summary

**Total observations in Section 11**: 35 binomial terms (19 Oxford + 16 Maryland).

**Discrepancies found**: 8 items requiring attention, of which 3 are consequential for the model:

1. **CONSEQUENTIAL**: The Gilman "Gil-F-rest" observation (n~37, y~19) is fabricated from arithmetic -- the extraction reports H-antibody data for only 27 of 64 controls (Table 5), but it is unclear whether the remaining 37 subjects had no H-antibody measurement or simply were not reported in that table. The plan assumes the 37 are a separate group with the Maryland-mixture CoP, but the extraction does not support this decomposition.

2. **CONSEQUENTIAL**: The Darton 2016 infection definition used in the plan ("bacteremia OR shedding") yields D-I-plac = 26/30 (87%), but this comes from the sensitivity analysis row ("Bacteraemia or stool positive"), not the primary outcome. The Waddington 2014 infection definition is "shedding" (stool positive). This is a genuine definition inconsistency between the two Oxford infection observations that the plan acknowledges in Section 4.1.2 but does not resolve in the likelihood.

3. **CONSEQUENTIAL**: The Hornick 10^3 discrepancy (Hornick 1966: 9/14 vs Hornick 1970 Part 1: 0/14) remains unresolved. The plan uses 0/14 from the 1970 paper. The extractions do not provide enough information to determine whether these are the same 14 subjects reclassified or different subjects.

**Minor discrepancies**: 5 items that are either approximations acknowledged in the plan or small metadata issues (detailed below).

---

## 2. Observation-by-Observation Verification Table

### Oxford Observations

| Obs ID | Plan (n, y, dose) | Extraction (n, y, dose) | Match? | Notes |
|--------|-------------------|------------------------|--------|-------|
| W-F-3 | 20, 11, 10^3 | 20, 11, 10^3 (Waddington Table 2) | **Yes** | Exact match. TD composite (fever OR bacteremia). |
| W-F-4 | 16, 10, 10^4 | 16, 10, 10^4 (Waddington Table 2) | **Yes** | Exact match. |
| W-F-5 | 5, 5, 10^5 | 5, 5, 10^5 (Waddington Table 2) | **Yes** | Exact match. |
| W-I-3 | 20, 13, 10^3 | 20, 13, 10^3 (Waddington shedding data) | **Yes** | Outcome = stool shedding at any time. |
| W-I-4 | 16, 10, 10^4 | 16, 10, 10^4 (Waddington shedding data) | **Yes** | Exact match. |
| W-I-5 | 5, 4, 10^5 | 5, 4, 10^5 (Waddington shedding data) | **Yes** | Exact match. |
| D-F-plac | 30, 20, 1.82e4 | 30, 20, 1.82e4 (Darton Table 2) | **Yes** | Per-protocol placebo. Dose is median actual. |
| D-F-Ty21a | 30, 13, 1.82e4 | 30, 13, 1.82e4 (Darton Table 2) | **Yes** | Exact match. |
| D-F-M01 | 31, 18, 1.82e4 | 31, 18, 1.82e4 (Darton Table 2) | **Yes** | Exact match. |
| D-I-plac | 30, 26, 1.82e4 | 30, 26 (Darton Table 2 sensitivity: "Bacteraemia or stool positive") | **Yes, but see note** | The value 26/30 = 87% comes from the sensitivity analysis row, not the primary TD outcome. The infection definition is "bacteraemia OR stool positive", which is broader than Waddington's shedding-only definition. See Section 4 below. |
| D-I-M01 | 31, 21, 1.82e4 | 31, 21 (Darton Table 2: 68% "Bacteraemia or stool positive") | **Yes** | Same definition note as D-I-plac. |
| D-I-Ty21a | 30, 16, 1.82e4 | 30, 16 (Darton Table 2: 53% "Bacteraemia or stool positive") | **Yes** | Same definition note. |
| J-F-ctrl | 31, 24, ~1e4 | 31, 24, 1-5x10^4 (Jin Table 2) | **Yes** | Per-protocol control. Note: 32 enrolled, 31 per-protocol. Plan correctly uses 31. |
| J-F-ViTT | 37, 13, ~1e4 | 37, 13, 1-5x10^4 (Jin Table 2) | **Yes** | Exact match. |
| J-F-ViPS | 35, 13, ~1e4 | 35, 13, 1-5x10^4 (Jin Table 2) | **Yes** | Exact match. |
| J-I-ctrl | 31, 22, ~1e4 | 31, 22 (Jin Table 2: "All participants with positive stool" 22/31 = 71%) | **Yes** | Infection = stool culture positive. |
| J-I-ViTT | 37, 22, ~1e4 | 37, 22 (Jin Table 2: 22/37 = 59%) | **Yes** | Exact match. |
| J-I-ViPS | 35, 21, ~1e4 | 35, 21 (Jin Table 2: 21/35 = 60%) | **Yes** | Exact match. |
| G20-F-naive | 19, 12, ~2.5e4 | 19, 12, 1-5x10^4 (Gibani 2020 Table 2) | **Yes** | ST Naive arm. Actual median dose 22.9-26.2 x 10^3 CFU. Plan's "~2.5e4" is reasonable. |

### Maryland Observations

| Obs ID | Plan (n, y, dose) | Extraction (n, y, dose) | Match? | Notes |
|--------|-------------------|------------------------|--------|-------|
| H-F-3 | 14, 0, 10^3 | 14, 0, 10^3 (Hornick 1970 Part 1, Table 1) | **Yes** | "The smallest dose employed, 1000 cells, failed to induce disease in any of the 14 volunteers." But see Hornick 1966 discrepancy (Section 4). |
| H-F-5 | 116, 32, 10^5 | 116, 32, 10^5 (Hornick 1970 Part 1, Table 1) | **Yes** | Exact match. |
| H-F-7 | REMOVED | 32, 16, 10^7 (Hornick 1970 Part 1, Table 1) | **N/A** | Correctly removed to avoid double-counting with H-I-7 + H-FgI-7. See Section 3. |
| H-F-8 | 9, 8, 10^8 | 9, 8, 10^8 (Hornick 1970 Part 1, Table 1) | **Yes** | Exact match. |
| H-F-9 | 42, 40, 10^9 | 42, 40, 10^9 (Hornick 1970 Part 1, Table 1) | **Yes** | Exact match. |
| H-I-7 | 30, 28, 10^7 | 30, 28, ~10^7 (Hornick 1970 Part 1, Table 2, Quailes only: 16 disease + 12 infection = 28 out of 30) | **Yes** | Infection = disease + infection-without-disease categories. Correctly uses Quailes-only N=30 (not Table 1's N=32). |
| H-FgI-7 | 28, 16, 10^7 | 28, 16, ~10^7 (Hornick 1970 Part 1, Table 2: 16/30 with disease, conditional on 28 infected) | **Yes** | Fever given infection. Exact match. |
| Gil-F-Hlo | 14, ~9, 10^5 | 14, ~61%, 10^5 (Gilman Table 5: H Ab <1:20, 61% attack rate) | **Approximate** | 61% of 14 = 8.54, rounded to ~9. Plan correctly marks as approximate (~9). Could also be 8. The extraction gives only percentages, not exact counts. |
| Gil-F-Hhi | 13, ~3, 10^5 | 13, ~24%, 10^5 (Gilman Table 5: H Ab >=1:20, 24% attack rate) | **Approximate** | 24% of 13 = 3.12, rounded to ~3. Plan correctly marks as approximate. |
| Gil-F-rest | ~37, ~19, 10^5 | **NOT IN EXTRACTION** | **Cannot verify** | This observation is derived by subtraction: 64 total - 14 - 13 = 37 remaining controls. The y value ~19 is derived from 31 total - 9 - 3 = 19. However, the Gilman extraction reports H-antibody data for only 27 subjects (Table 5, which has n=14 and n=13). It is unclear whether the remaining 37 had H-antibody measured but unreported, or were never tested. **This is a fabricated observation.** See Recommendation 1. |
| Gil-I-ctrl | 43, 26, 10^5 | 43, 26, 10^5 (Gilman Table 4: Trials 1&3 controls, late shedding 4-30 days = 26/43 = 60%) | **Yes, but definition note** | The plan uses "late shedding (4-30 days)" as the infection proxy. This is a deliberate choice -- early shedding (0-3 days) may represent transit rather than true infection. Extraction confirms 26/43 for late shedding. However, the plan labels this as n=43, which is Trials 1&3 controls only, not the full 64 controls. |
| Lev-F-1 | 26, 13, 10^5 | 26, 13, 10^5 (Levine Table 2: Trial 1 controls, 13/26 = 50%) | **Yes** | Exact match. |
| Lev-F-2 | 33, 10, 10^5 | 33, 10, 10^5 (Levine Table 2: Trial 2 controls, 10/33 = 30%) | **Yes** | Exact match. |
| Lev-F-3 | 22, 12, 10^5 | 22, 12, 10^5 (Levine Table 2: Trial 3 controls, 12/22 = 55%) | **Yes** | Exact match. |
| Lev-F-4 | 16, 4, 10^5 | 16, 4, 10^5 (Levine Table 2: Trial 4 controls, 4/16 = 25%) | **Yes** | Exact match. |
| Lev-I-1 | 26, 19, 10^5 | 26, 19, 10^5 (Levine Table 2: Trial 1 controls, stool positive anytime = 19/26 = 73%) | **Yes** | Infection = any positive stool culture at any time after challenge. |

### Summary Counts

- **Exact match**: 27/35
- **Approximate (marked ~)**: 2/35 (Gil-F-Hlo, Gil-F-Hhi)
- **Cannot verify**: 1/35 (Gil-F-rest)
- **Correct but with definition note**: 4/35 (D-I-plac, D-I-M01, D-I-Ty21a, Gil-I-ctrl)
- **Correctly excluded**: 1/35 (H-F-7)

---

## 3. Cohort Independence Assessment

### 3.1 Hornick 1970 Table 1 (Part 1) vs Table 5 (Part 2) overlap

**Plan claim**: The "None" row in Table 5 (Part 2: 4/4 at 10^9, 30/15=50% at 10^7, 104/28=27% at 10^5) overlaps substantially with Table 1 (Part 1). The plan uses Table 1 as canonical for the multi-dose curve and excludes Table 5's unvaccinated rows.

**Assessment**: **Correct.** The Table 5 "None" row numbers exactly match Table 1 at the corresponding doses (N=4/4 at 10^9, N=30/15 at 10^7, N=104/28 at 10^5). However, Table 1 has N=32 at 10^7 while Table 5 has N=30. The extraction notes this: Table 2 (Part 1) uses Quailes-only N=30, while Table 1 may include 2 additional subjects from other batches or time periods. Table 5 likely uses the Quailes-only denominator at 10^7. Regardless, the plan's resolution (use Table 1 for the multi-dose curve; use only vaccinated arms from Table 5) is sound.

### 3.2 Gilman 1977 controls vs Hornick 1970

**Plan claim**: Independent. Gilman trials conducted 1973-1976, after Hornick 1970 enrollment.

**Assessment**: **Plausible but not certain.** The Gilman extraction states trials were conducted at the Maryland House of Correction, Jessup, Maryland -- the same institution as Hornick. However, Gilman 1977 was published later and the trials were in 1973-1976 (based on Levine 1976 Trial 3 = 1972, Trial 4 = 1973, and Gilman being a subsequent study). The Hornick 1970 review paper covers studies through ~1970. Given the temporal gap (Hornick enrollment ended ~1970, Gilman started ~1973), different volunteers are likely but cannot be proven from the extractions alone. The prison population turned over. **Treating as independent is reasonable.**

### 3.3 Levine 1976 Trial 1 (1970) vs Hornick 1970

**Plan claim**: Possibly overlapping for Trial 1 (1970); Trials 2-4 likely independent.

**Assessment**: **Correct concern.** Levine Trial 1 was conducted in 1970 at the same institution. If Hornick's multi-dose enrollment extended into 1970, some subjects could overlap. However, Levine Trial 1 controls (n=26) were challenged at 10^5, and Hornick Table 1 already has N=116 at 10^5. If these 26 are a subset of the 116, including both would double-count. The plan flags this appropriately. **Recommend: treat Levine Trial 1 as potentially overlapping with Hornick; run sensitivity excluding it.**

### 3.4 Gibani 2020 naive arms vs earlier Oxford studies

**Plan claim**: Independent (new naive recruits).

**Assessment**: **Correct.** The Gibani 2020 extraction (Table 1) describes the ST Naive arm (n=19) as new volunteers recruited for the PATCH study (OVG2014/01, conducted 2015-2017). These were screened for no prior typhoid vaccination or infection. The rechallenge arms overlap with earlier studies (Waddington 2014, Darton 2016, Jin 2017), but the naive arm is independent.

### 3.5 Darton 2017 as same cohort as Waddington 2014

**Plan claim**: Darton 2017 is a secondary analysis of the Waddington 2014 cohort; excluded from primary analysis.

**Assessment**: **Correct.** The Darton 2017 extraction explicitly states: "This is a secondary analysis of the Waddington et al. 2014 CHIM cohort" and "all participants are the same 41 individuals from Waddington 2014." Exclusion from the primary analysis is appropriate.

### 3.6 Darton 2016 vs Waddington 2014

**Plan claim**: Independent (different study: OVG2011/02 vs OVG2009/10).

**Assessment**: **Correct.** Different trial registrations, different enrollment periods (Darton: Nov 2011 - Jun 2012; Waddington: earlier), different participants. The Darton extraction explicitly references Waddington as a prior study, not a related analysis.

### 3.7 Gibani 2019 pooled analysis

**Plan claim**: DO NOT USE -- not independent; reanalysis of primary study data.

**Assessment**: **Correct.** The Gibani 2019 extraction explicitly lists all 6 included studies (OVG2009/10 = Waddington, OVG2011/02 = Darton, OVG2014/08 = Jin, OVG2014/01 = Gibani 2020, plus two others). Using this alongside the primary papers would double-count.

---

## 4. Cross-Study Consistency Findings

### 4.1 Hornick 10^3 Discrepancy (Item 8)

**Hornick 1966 extraction**: Figure 2 reports 10^3 dose: N=14, 9 ill, 64.3% attack rate. The extraction notes this is "anomalously high -- higher than 10^5 (27.6%)."

**Hornick 1970 Part 1 extraction**: Table 1 reports 10^3 dose: N=14, 0 with disease, 0% attack rate. The text states: "The smallest dose employed, 1000 cells, failed to induce disease in any of the 14 volunteers."

**Is there an explanation?** The extractions do not provide a definitive explanation. The same N=14 appears in both papers. Possible explanations:
- Hornick 1966 Figure 2 aggregated data from "cumulative totals across multiple experiments" that may have included non-Quailes strains or subjects treated with streptomycin pretreatment (Hornick 1970 Part 1 mentions 14 volunteers given 1000 bacteria where streptomycin pretreatment caused 1/4 treated volunteers to become ill).
- Hornick 1970 may represent a reclassification using stricter criteria.
- Hornick 1966 Figure 2 numbers were estimated from a bar chart and may be less reliable.

**Does Hornick 1970 Part 1 explicitly report 0/14?** Yes: "1000 cells, failed to induce disease in any of the 14 volunteers" (text) and Table 1 shows 0 disease at 10^3 dose.

**Recommendation**: The plan's choice to use the 1970 value (0/14) is defensible since it is from a later, more detailed publication with explicit textual confirmation. However, the 1966 value (9/14) should be tested as a sensitivity analysis (the plan already specifies this in Section 9.3).

### 4.2 Hornick Table 1 vs Table 2 at 10^7 (Item 9)

**Table 1 (Part 1)**: N=32, 16 disease at 10^7.
**Table 2 (Part 1, Quailes only)**: N=30, 16 disease, 12 infection-without-disease, 2 no infection.

**Are the 16 disease cases the same people?** Almost certainly yes -- same dose, same study, same strain. The extraction notes: "Slight discrepancy in N for 10^7 dose (Table 1: N=32, 16 disease; Table 2 Quailes only: N=30, 16 disease). This suggests Table 1 may include 2 additional volunteers from other Vi strains, or different analysis subsets."

**Where do the extra 2 subjects in Table 1 come from?** The extraction suggests they may be from non-Quailes Vi-positive strains (Zermatt or Ty2V) that were lumped into the 10^7 summary in Table 1 but excluded in the Quailes-only Table 2. Alternatively, Table 1 may include subjects from a slightly different time period.

**Is H-F-7 correctly REMOVED?** Yes. The plan removes H-F-7 (32, 16, 10^7 from Table 1) because it would double-count with H-I-7 (30, 28, 10^7 from Table 2) + H-FgI-7 (28, 16 from Table 2). The Table 2 data is more informative because it separates infection from disease. Using Table 1's 32/16 alongside Table 2's 30/28 and 28/16 would count the disease cases twice. The plan's resolution is correct.

**Minor note**: The plan uses N=30 (Table 2, Quailes only) for the infection observation, not N=32 (Table 1). This is the right choice for strain consistency but means 2 subjects at 10^7 are excluded from the model entirely.

### 4.3 Oxford Infection Definition Consistency (Item 10)

**Waddington 2014**: Infection = stool shedding at any time. Rates: 13/20 (65%) at 10^3, 10/16 (63%) at 10^4, 4/5 (80%) at 10^5.

**Darton 2016**: The plan uses "bacteremia OR stool positive" as the infection definition. From the Darton extraction (Table 2 sensitivity analysis): Placebo 26/30 (87%), M01ZH09 21/31 (68%), Ty21a 16/30 (53%).

**Jin 2017**: The plan uses "All participants with positive stool" (stool culture positivity). From the Jin extraction (Table 2): Control 22/31 (71%), Vi-TT 22/37 (59%), Vi-PS 21/35 (60%).

**Inconsistency**: Waddington and Jin use shedding-only definitions, while Darton uses a broader "bacteremia OR shedding" definition. This is acknowledged in the plan (Section 4.1.2) but not resolved. The broader Darton definition will yield higher infection rates than a shedding-only definition would.

**What are the shedding-only rates for Darton?** The extraction does not provide shedding-only rates per group. The stool shedding data (Table 4) reports aggregate early shedding across all groups (44/91 = 49%) but does not break this out by arm for the "any stool positive" definition.

**How large is the gap?** For Darton Placebo: 26/30 (87%) with bacteremia OR stool+ vs 20/30 (67%) with TD (fever OR bacteremia). The "any bacteremia" rate for placebo was 20/30 (67%). So shedding contributes additional cases beyond bacteremia alone. Without per-group shedding-only data, we cannot compute the exact gap, but the broader definition adds roughly 6/30 = 20% more infected subjects in the placebo arm.

**Recommendation**: For consistency, the model should ideally use the same infection definition across all Oxford studies. Shedding-only data for Darton 2016 by arm should be extracted if possible (from S1 Dataset). Alternatively, acknowledge that the Darton infection observations use a slightly broader definition, which will pull the infection curve upward at the Darton dose point.

### 4.4 Maryland Outcome Definitions Across Studies (Item 11)

Three Maryland studies use overlapping but non-identical fever/disease definitions:

**Hornick 1966/1970**: "Temperature exceeded 103F for over 24 to 36 hours" + "clinical manifestations of typhoid fever were present." Treatment with chloramphenicol initiated at that point.

**Gilman 1977**: "Presence of fever AND blood or stool culture positive for S. typhi." Treatment thresholds: fever >103F for 1 day, OR >101F for 3 days, OR >100F for 5 days.

**Levine 1976**: "Acute illness with oral temperature >= 101F accompanied by isolation of S. typhi from blood or stool." Treatment thresholds: >=103F for 24 hr OR >=101F for 48 hr.

**Are these the same outcome?**

No -- there are meaningful differences:

1. **Hornick is the strictest**: Requires >103F sustained for 24-36 hours. Fever alone suffices (no culture confirmation required for the definition, though cultures were done).

2. **Gilman is intermediate**: Requires fever (any level apparently, given the graduated treatment thresholds) PLUS culture confirmation. The treatment thresholds (103F/1d, 101F/3d, 100F/5d) are more liberal than Hornick's single 103F threshold.

3. **Levine is the most liberal**: Requires only 101F plus culture isolation. This is a substantially lower fever threshold than Hornick's 103F.

**Impact on the model**: Gilman and Levine would capture more "disease" cases than Hornick at the same dose and population, because their fever thresholds are lower. This could partially explain why Gilman controls (48% at 10^5) and some Levine controls (50-55%) have higher attack rates than Hornick (28% at 10^5). The plan treats all Maryland studies as using the same definition and applies a single phi correction, which may be an oversimplification.

**Recommendation**: The definition heterogeneity within the Maryland era should be noted as a limitation. Consider either: (a) applying a study-specific phi for Gilman and Levine (capturing the lower fever threshold), or (b) absorbing this into the study-level random effect (which the plan already includes). Option (b) is simpler and probably sufficient.

---

## 5. Missing Data Inventory

### 5.1 Data in Extractions NOT Used in the Plan (Item 12)

**Incubation period data:**
- Hornick 1970 Part 1, Table 1: Median and range incubation periods at each dose (5 days at 10^9, 7.5 at 10^7, 9 at 10^5). The plan notes this in Appendix C as "not included but could be added."
- Waddington 2014, Table 3: Time to diagnosis by dose (8.1 days at 10^3, 6.4 days at 10^4/10^5).
- Darton 2016: Median time to TD (265h M01ZH09, 172h placebo).
- Jin 2017: Time to diagnosis (6.0 days control, 6.5 Vi-TT, 7.2 Vi-PS).
- Gibani 2020: Time to diagnosis by group.

**Quantitative bacteremia (CFU/mL):**
- Darton 2016: Median 1.30 CFU/mL placebo, 0.13 M01ZH09, 0.05 Ty21a.
- Jin 2017: Median 0.4 CFU/mL control, 0.075 Vi-TT, 0.1 Vi-PS.

**Other data not used:**
- Darton 2017 evidence for primary bacteremia (S. Typhi DNA in blood within 48h).
- Darton 2017 estimate of asymptomatic infection (30/41 = 73% total infection rate vs 24/41 = 59% by standard criteria).
- Hornick 1970 Part 1 streptomycin pretreatment data (1/4 treated became ill at 10^3 dose; possibly explains the 1966 discrepancy).
- Gibani 2020 rechallenge data by prior disease status (10/38 = 26% vs 25/37 = 68%). The plan includes this as a holdout test but not in the likelihood.
- Hornick 1970 Part 2 re-challenge data ("illness rate was now about 25 per cent" but exact N not given).
- Levine 1976 vaccine arms (fresh vaccine: 5/30 = 17%, 3/26 = 12%; lyophilized: 16/49 = 33%, 3/17 = 18%).
- Gibani 2020 SPT-ST cross-challenge data (7/10 = 70%).
- Jin 2017 fever threshold sensitivity data for controls (>=37.5: 20/31 = 65%; >=38.0: 17/31 = 55%; >=38.5: 14/31 = 45%). **Note**: Jin does NOT report >=39.0°C; the >=39.0°C data (9/30 = 30%) is from Darton 2016 placebo Table 2. This was previously misattributed to Jin in the plan (corrected 2026-03-19).
- Darton 2016 fever threshold sensitivity data for all arms.
- Gilman 1977 shedding data for Trial 2 controls (8/21 late shedding = 38%).

### 5.2 Metadata Needed for Model but Missing from Extractions (Item 13)

**Exact actual doses administered:**
- Darton 2016 reports median 1.82x10^4 (range 1.46-2.66x10^4). The plan uses 1.82e4.
- Jin 2017 reports "1-5 x 10^4" as a range, with no median actual dose. The plan uses "~1e4" which is the lower bound of the range.
- Gibani 2020 reports actual median 22.9-26.2 x 10^3 CFU. The plan uses "~2.5e4".
- Waddington 2014 reports target doses of 10^3, 10^4, 10^5 without actual measured doses per subject (though Darton 2017 mentions the dose was verified).
- **All Maryland studies**: Doses are target doses from plate counts of the inoculum. Actual delivered dose depends on how much the volunteer ingested. No per-subject dose verification.

**Individual-level titer data for Oxford vaccinated arms:**
- Darton 2016: S1 Dataset exists but is not extracted. Only group-level attack rates and the anti-Vi hazard ratio are available.
- Jin 2017: Only GMT reported (563 for Vi-TT, 141 for Vi-PS). Individual titer distributions not in extraction.

**Cross-tabulated infection x fever counts:**
- For the cascaded model, the ideal data is: (infected+fever, infected+no fever, not infected) counts per group.
- Available only for Hornick Table 2 at 10^7 (16 disease, 12 infection-no-disease, 2 no infection) and partially for Waddington (bacteremia only, fever only, both).
- Not available for Darton 2016, Jin 2017, or Gibani 2020 in per-group cross-tabulated form.

**Pre-challenge anti-Vi IgG distribution parameters:**
- The plan uses the fraction with detectable anti-Vi (e.g., 40% for Darton placebo, 38% for Jin control). But the distribution of titers among the Vi+ subgroup is not extracted (no GMT or percentiles for the Vi+ subgroup of the naive/control arms).

**Waddington 2014 baseline anti-Vi**: The extraction says "Anti-Vi IgG measured but no statistically significant differences between those who met diagnostic criteria and those who did not." The fraction with detectable anti-Vi is not reported in the Waddington extraction. The plan uses "~29%" from comparable Oxford cohort rates, which is a reasonable approximation but not directly confirmed from this study.

---

## 6. Recommendations (Ranked by Consequence)

### Critical (resolve before Stan implementation)

**1. Resolve the Gil-F-rest observation.**
The plan creates a "remaining 37 controls" observation for Gilman that is not supported by the extraction. Options:
- (a) Drop Gil-F-rest entirely and use only the H-antibody-stratified observations (Gil-F-Hlo and Gil-F-Hhi, total n=27). This loses 37 subjects but avoids fabrication.
- (b) Use Gil-F-ctrl (n=64, y=31) as a single observation with the Maryland mixture model, and drop the stratified observations. This loses the within-study immunity contrast.
- (c) Check the primary paper to determine whether H-antibody was measured in all 64 controls or only 27. If all 64, the remaining 37 can be assigned by subtraction. If only 27, use option (a).
**Recommendation**: Option (c) first. If not possible, option (a) is safest.

**2. Harmonize the Oxford infection definition.**
Waddington and Jin use shedding-only; Darton uses bacteremia OR shedding. This affects the D-I-* observations. Options:
- (a) Extract shedding-only rates for Darton 2016 from S1 Dataset.
- (b) Accept the broader definition for Darton and note it as a limitation.
- (c) Use "any bacteremia" as the infection definition for all Oxford studies (available for Waddington: 8/20, 9/16, 3/5; and Darton: 16/31, 20/30, 11/30; and Jin: available from diagnosis type breakdown).
**Recommendation**: Option (b) is pragmatic. The plan should note that D-I-* uses a broader infection definition than W-I-* and J-I-*.

**3. Verify the Hornick 10^3 data against primary PDFs.**
The 1966 vs 1970 discrepancy (9/14 vs 0/14) is the single highest-influence data point for the low-dose behavior of the Maryland curve. The plan already flags this and requires a sensitivity analysis, which is appropriate. However, the primary publication should be checked to determine if the 1966 value includes non-Quailes subjects, streptomycin-pretreated subjects, or uses a different disease definition.

### Important (should be addressed but model can proceed)

**4. Address Maryland fever definition heterogeneity.**
Gilman (fever + culture) and Levine (101F + culture) use more liberal definitions than Hornick (103F for 24-36h). This could explain part of the higher attack rates in Gilman/Levine controls compared to Hornick. The study-level random effect will absorb some of this, but a systematic bias in one direction is not well-captured by a symmetric random effect. Consider noting this as a structural limitation.

**5. Verify Gilman H-antibody subsample representativeness.**
The plan lists this as Prerequisite 2 (Appendix). The extraction shows n=14 (H<1:20) + n=13 (H>=1:20) = 27 with H-antibody data, out of 64 total controls. If these 27 are a biased subsample (e.g., only Trial 2 controls, or only subjects who consented to additional blood draws), the stratified attack rates may not represent the full population.

**6. Check the approximate y values for Gil-F-Hlo and Gil-F-Hhi.**
The extraction gives only percentages (61% and 24%). Back-calculation: 61% of 14 = 8.54 (could be 8 or 9); 24% of 13 = 3.12 (could be 3). The plan uses ~9 and ~3. If these are wrong by 1, the impact is modest given small n, but the primary paper should confirm the exact counts.

### Minor (note but no action needed)

**7. Jin 2017 dose uncertainty.**
The plan uses "~1e4" for Jin 2017, but the actual dose is reported as "1-5 x 10^4." The Darton 2016 median at the same Oxford CHIM protocol was 1.82x10^4. Using 1e4 vs 2e4 affects the dose-response curve position. Consider using 2e4 as a point estimate (consistent with Darton) or modeling dose uncertainty.

**8. Levine 1976 infection definition for Lev-I-1.**
The plan uses "any stool+" (19/26 = 73%) from Trial 1. This is "stool positive anytime after challenge." The Gilman infection proxy Gil-I-ctrl uses "late shedding 4-30 days" (26/43 = 60%). These are different shedding definitions. The "anytime" definition will yield higher rates. For consistency, either use "anytime" for both or "late shedding" for both. The plan does not note this inconsistency.

---

## Appendix: Detailed Notes on Specific Observations

### A.1 Gilman H-Antibody Data (Table 5)

The Gilman extraction reports Table 5 as:

| Baseline H Antibody | N | Attack Rate |
|--------------------|---|-------------|
| < 1:20 | 14 | 61% |
| >= 1:20 | 13 | 24% |

Total in Table 5: 27. Total controls: 64 (43 from Trials 1&3 + 21 from Trial 2).

The plan assumes 14 + 13 = 27 with H-antibody data, leaving 64 - 27 = 37 without. This arithmetic is correct, but the extraction does not state whether H-antibody was measured in all controls or only a subsample. The Table 5 header in the extraction is "H Antibody and Protection in Controls" -- the Gilman paper's Table 5 title and scope should be checked.

Furthermore, if the 27 are from a single trial arm (e.g., Trials 1&3 only, n=43), then 43 - 27 = 16 would remain from Trials 1&3 and all 21 from Trial 2. But this is speculative.

### A.2 Levine H-Antibody Data

The Levine extraction reports: "Control volunteers with preexistent H-antibody titers >=20 had lower attack rate than those without (22% vs. 48%; P = 0.005)." This is stated as a general finding across all controls (n=97 pooled). The exact stratification (n per stratum) is not given in the extraction. The plan uses this only as confirmation of the Gilman finding and the pi_susc prior, not as a separate observation. This is appropriate.

### A.3 Darton 2016 Dose

The plan uses 1.82e4 CFU (median actual dose). The extraction confirms: "Actual median dose administered: 1.82 x 10^4 CFU (range: 1.46-2.66 x 10^4 CFU)." This is more precise than other Oxford studies where only the target range (1-5 x 10^4) is reported.

### A.4 Gibani 2020 Dose

The plan uses "~2.5e4". The extraction reports actual median doses of 22.9-26.2 x 10^3 CFU (varying by arm, from IQR data). The ~2.5e4 point estimate is reasonable.

### A.5 Waddington Infection Definition

The extraction reports stool shedding at any time: 13/20 (65%) at 10^3, 10/16 (63%) at 10^4, 4/5 (80%) at 10^5. The plan uses these as the W-I-* observations. This is a shedding-at-any-time definition, distinct from Gilman's "late shedding (4-30 days)" definition and Levine's "stool positive anytime." All three are slightly different, though all attempt to capture infection.

### A.6 Hornick Table 5 Vaccine Data

The plan includes vaccinated arms from Table 5 (Part 2):
- K at 10^5: 43, 4 -- matches extraction (43, 4, 9%)
- L at 10^5: 45, 3 -- matches extraction (45, 3, 7%)
- K at 10^7: 28, 12 -- matches extraction (28, 12, 43%)
- L at 10^7: 24, 13 -- matches extraction (24, 13, 54%)

All four match exactly.

The Vi vaccine arms (Vi at 10^5: 17, 3; Vi at 10^7: 14, 10) are not included in the plan's primary analysis but are noted in the extraction.
