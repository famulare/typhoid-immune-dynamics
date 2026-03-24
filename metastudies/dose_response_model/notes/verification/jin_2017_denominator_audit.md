# Jin et al. 2017 -- Table 2 Denominator Audit

**Auditor:** Claude Opus 4.6 (1M context)
**Date:** 2026-03-18
**Source:** Jin C et al. Lancet 2017;390:2472-2480. PDF read directly.

---

## 1. Per-Protocol Population (Reference Denominators)

From Figure 1 (trial profile) and Table 2:

| Group | Enrolled | Challenged | Per-Protocol (analyzed) |
|-------|----------|------------|------------------------|
| Control | 34 | 32 | **31** (1 excluded: antibiotics commenced before diagnosis or completion) |
| Vi-TT | 41 | 37 | **37** |
| Vi-PS | 37 | 35 | **35** |

The per-protocol Ns are **31, 37, 35**. These are the full-cohort denominators.

Number diagnosed (primary composite endpoint = bacteraemia OR persistent fever >= 38C for >= 12h):

| Group | Number Diagnosed | Denominator |
|-------|------------------|-------------|
| Control | 24 | /31 (77%) |
| Vi-TT | 13 | /37 (35%) |
| Vi-PS | 13 | /35 (37%) |

---

## 2. Blood Culture Sampling Protocol

### What the Methods say (p.2474, "Procedures" section):

> "Following challenge, participants were seen daily for vital sign measurement, **blood collection**, and general assessment in an outpatient clinic for a 2 week period"

Key details:
- **Blood (10 mL) samples** were collected for culture, haematological and biochemical testing
- Blood was collected **daily** during the 2-week challenge period
- This was done for **ALL participants**, not just symptomatic ones
- Stool culture samples were also collected and processed

### Quantitative blood culture:

> "Quantitative blood culture (10 mL) was done **at the time of diagnosis** and processed according to the manufacturer's instructions"

This is a critical distinction:
- **Daily blood culture**: performed on ALL participants throughout the 14-day challenge period (surveillance cultures)
- **Quantitative blood culture** (10 mL with CFU/mL quantification): performed only **at the time of diagnosis**

### Diagnosis criteria (p.2474):

> "Typhoid fever was diagnosed if pre-determined criteria were met: a positive blood culture with S Typhi collected more than 72 h post-challenge **or** a fever of 38C or higher persisting for 12 h or longer"

So bacteraemia (positive blood culture > 72h post-challenge) is itself one of the two diagnostic criteria in the composite endpoint. A participant could be diagnosed EITHER by:
1. Positive blood culture (= "microbiological diagnosis"), OR
2. Persistent fever >= 38C for >= 12h (= "clinical diagnosis")

---

## 3. Row-by-Row Audit of Table 2

### PRIMARY OUTCOME section

| Row | Control | Vi-TT | Vi-PS | Denominator Type | Notes |
|-----|---------|-------|-------|------------------|-------|
| Completed challenged | 31 | 37 | 35 | N/A (counts) | These ARE the per-protocol Ns |
| Total diagnosed | 24/31 (77%) | 13/37 (35%) | 13/35 (37%) | **Full cohort** | Denom = per-protocol N |

**Verification:** Matches PDF exactly. The denominators 31, 37, 35 are the full per-protocol populations.

### SECONDARY OUTCOMES section

| Row | Control | Vi-TT | Vi-PS | Denominator Type | Notes |
|-----|---------|-------|-------|------------------|-------|
| Microbiological diagnosis | 16/31 (52%) | 12/37 (32%) | 9/35 (26%) | **Full cohort** | Denom = per-protocol N. These are participants diagnosed via positive blood culture. |
| Clinical diagnosis | 8/31 (26%) | 1/37 (3%) | 4/35 (11%) | **Full cohort** | Denom = per-protocol N. These are participants diagnosed via persistent fever only (no preceding positive blood culture). |

**Cross-check:** Microbiological + Clinical should equal Total Diagnosed.
- Control: 16 + 8 = 24. Matches.
- Vi-TT: 12 + 1 = 13. Matches.
- Vi-PS: 9 + 4 = 13. Matches.

### CLINICAL OUTCOMES section

| Row | Control | Vi-TT | Vi-PS | Denominator Type |
|-----|---------|-------|-------|------------------|
| Fever >= 37.5C (any duration) | 20/31 (65%) | 13/37 (35%) | 18/35 (51%) | **Full cohort** |
| Fever >= 38.0C (any duration) | 17/31 (55%) | 6/37 (16%) | 11/35 (31%) | **Full cohort** |
| Fever >= 38.5C (any duration) | 14/31 (45%) | 4/37 (11%) | 9/35 (25%) | **Full cohort** |

**Verification:** All denominators are 31, 37, 35. These are full-cohort measures -- fever was monitored for ALL participants via twice-daily self-reported temperature measurements for 21 days.

### MICROBIOLOGICAL OUTCOMES section

| Row | Control | Vi-TT | Vi-PS | Denominator Type | CRITICAL NOTES |
|-----|---------|-------|-------|------------------|----------------|
| **S Typhi bacteraemia** | **24/24 (100%)** | **13/13 (100%)** | **11/13 (85%)** | **DIAGNOSED SUBSET** | Denominators are 24, 13, 13 = number diagnosed, NOT 31, 37, 35 |
| Time to first positive blood culture (days) | 6.1 (5.0-7.6) | 6.5 (6.1-8.6) | 6.1 (5.0-10.2) | Among diagnosed | Median (IQR) |
| Participants with positive S Typhi stool culture | **22/31 (71%)** | **22/37 (59%)** | **21/35 (60%)** | **Full cohort** | Denom = per-protocol N |
| Diagnosed participants with positive S Typhi stool culture | **19/24 (79%)** | **12/13 (92%)** | **10/13 (77%)** | **DIAGNOSED SUBSET** | Denom = number diagnosed |
| Median quantitative blood culture (CFU/mL; range) | 0.4 (0.05-22.7) | 0.075 (0.05-1.2) | 0.1 (0.05-5.6) | Among diagnosed | Quantitative culture done at time of diagnosis only |

---

## 4. The Critical Bacteraemia Denominator Issue

### The problem

The "S Typhi bacteraemia" row reports:
- Control: 24/24 (100%)
- Vi-TT: 13/13 (100%)
- Vi-PS: 11/13 (85%)

The denominators (24, 13, 13) are exactly the number of diagnosed participants, NOT the full per-protocol cohort (31, 37, 35).

### But does this mean blood cultures were only taken from diagnosed participants?

**No.** The Methods clearly state that blood was collected **daily from ALL participants** during the 2-week challenge period. The daily blood cultures were the surveillance mechanism used to detect asymptomatic bacteraemia.

### Then why is the bacteraemia row restricted to diagnosed participants?

This is the crux of the issue. There are two possible interpretations:

**Interpretation A: Reporting convention.** The authors chose to report bacteraemia prevalence only among diagnosed participants. This would mean that **among the undiagnosed participants, some may have had transient positive blood cultures that were not counted in this row** (or none did, but the denominator was still restricted). This interpretation is supported by the fact that the row is labeled "S Typhi bacteraemia" without qualification, appearing under "Microbiological outcomes" alongside other rows that mix full-cohort and diagnosed-subset denominators.

**Interpretation B: Bacteraemia was genuinely zero among undiagnosed.** Since daily blood cultures were taken from everyone, if any undiagnosed participant had a positive blood culture > 72h post-challenge, they would have met the microbiological diagnosis criterion and been diagnosed. Therefore, by definition, all bacteraemic participants ARE diagnosed participants. The denominator restriction to diagnosed participants is tautological -- bacteraemia is a subset of diagnosis.

**Interpretation B is almost certainly correct**, and here is why:

1. Positive blood culture > 72h post-challenge is one of the two criteria for diagnosis. Any participant with bacteraemia is automatically diagnosed (microbiological diagnosis).
2. The "Microbiological diagnosis" row shows: Control 16/31, Vi-TT 12/37, Vi-PS 9/35. These are people diagnosed BY positive blood culture.
3. The "S Typhi bacteraemia" row (24/24, 13/13, 11/13) is larger than the microbiological diagnosis row. The difference represents participants who were diagnosed clinically (by fever) AND ALSO had bacteraemia detected at some point.
4. Cross-check for Control: 16 diagnosed microbiologically + 8 diagnosed clinically = 24. Of the 24, all 24 had bacteraemia. This means all 8 clinically-diagnosed control participants ALSO had positive blood cultures (but were diagnosed by fever first, or had bacteraemia detected simultaneously with or after the fever criterion was met).
5. Cross-check for Vi-PS: 9 diagnosed microbiologically + 4 diagnosed clinically = 13. Of the 13, only 11 had bacteraemia. This means 2 of the 13 diagnosed Vi-PS participants were diagnosed by fever alone and never had a positive blood culture during the entire challenge period.

### The key logical point

**No undiagnosed participant can have bacteraemia**, because a positive blood culture automatically triggers diagnosis. Therefore:

- The TRUE full-cohort bacteraemia rates are: Control 24/31, Vi-TT 13/37, Vi-PS 11/35
- These are **identical to the diagnosis rates** for Control and Vi-TT (where 100% of diagnosed had bacteraemia)
- For Vi-PS, the full-cohort bacteraemia rate would be 11/35 (31%), slightly lower than the diagnosis rate of 13/35 (37%) because 2 diagnosed Vi-PS participants never had bacteraemia

### Can we reconstruct full-cohort bacteraemia rates?

Yes:

| Group | Bacteraemia among diagnosed | Undiagnosed (all bacteraemia-negative by definition) | **Full-cohort bacteraemia** |
|-------|---------------------------|-----------------------------------------------------|--------------------------|
| Control | 24/24 | 7/7 negative | **24/31 = 77.4%** |
| Vi-TT | 13/13 | 24/24 negative | **13/37 = 35.1%** |
| Vi-PS | 11/13 | 22/22 negative | **11/35 = 31.4%** |

---

## 5. Classification of Every Table 2 Endpoint

| Endpoint | Denominator Type | Control Denom | Vi-TT Denom | Vi-PS Denom | Usable for dose-response? |
|----------|------------------|---------------|-------------|-------------|--------------------------|
| Total diagnosed (composite) | Full cohort | 31 | 37 | 35 | Yes -- this is the primary endpoint |
| Microbiological diagnosis | Full cohort | 31 | 37 | 35 | Yes |
| Clinical diagnosis | Full cohort | 31 | 37 | 35 | Yes (but note: this is diagnosis by fever only, excluding those caught by blood culture first) |
| Fever >= 37.5C | Full cohort | 31 | 37 | 35 | Yes |
| Fever >= 38.0C | Full cohort | 31 | 37 | 35 | Yes |
| Fever >= 38.5C | Full cohort | 31 | 37 | 35 | Yes |
| S Typhi bacteraemia | **Diagnosed subset** | 24 | 13 | 13 | **See Section 6** |
| Stool culture (all participants) | Full cohort | 31 | 37 | 35 | Yes |
| Stool culture (diagnosed) | Diagnosed subset | 24 | 13 | 13 | No -- confounded by subset selection |
| Quantitative blood culture | Diagnosed subset | (diagnosed) | (diagnosed) | (diagnosed) | No -- only measured at diagnosis |

---

## 6. Assessment of Specific Endpoints for Dose-Response Modeling

### 6a. Stool shedding (22/31, 22/37, 21/35) -- FULL COHORT

**Usable: YES, with caveats.**

These are full-cohort denominators. S Typhi was detected from stool cultures in both diagnosed and undiagnosed participants ("S Typhi was detected from stool cultures in diagnosed and undiagnosed participants irrespective of vaccine assignment" -- p.2477).

Caveats:
- Stool shedding is a measure of GI colonization/infection, not systemic disease
- The high rates across all groups (59-71%) with modest group differences suggest stool shedding may reflect mucosal infection that Vi-vaccination does not strongly prevent
- This is consistent with Vi being a systemic (anti-capsular) immune mechanism, not a mucosal one
- Stool shedding may be a more complete measure of "infection" (any bacterial colonization) vs. "disease" (systemic infection with symptoms)

**For dose-response modeling:** Stool positivity could serve as an infection endpoint (as distinct from disease), but the narrow range of attack rates (59-71%) across vaccinated and unvaccinated groups means it has low discriminating power for vaccine efficacy. It may be most useful as evidence that Vi-vaccination does not prevent mucosal colonization.

### 6b. Bacteraemia (24/24, 13/13, 11/13) -- DIAGNOSED SUBSET

**As reported: NOT directly usable** (denominator is the diagnosed subset).

**However, full-cohort bacteraemia CAN be reconstructed** (see Section 4):
- Control: 24/31 = 77.4%
- Vi-TT: 13/37 = 35.1%
- Vi-PS: 11/35 = 31.4%

This reconstruction is valid because:
1. Daily blood cultures were taken from ALL participants
2. Any positive blood culture > 72h automatically triggers diagnosis
3. Therefore, NO undiagnosed participant can be bacteraemic
4. The full-cohort bacteraemia numerator = the diagnosed-subset bacteraemia numerator

**Important:** For Control and Vi-TT groups, the reconstructed full-cohort bacteraemia rate is IDENTICAL to the composite diagnosis rate (because 100% of diagnosed had bacteraemia). For Vi-PS, it is slightly lower (11/35 = 31.4% vs 13/35 = 37.1%) because 2 Vi-PS diagnosed participants were diagnosed by fever alone without bacteraemia.

**For dose-response modeling:** The reconstructed full-cohort bacteraemia rates are valid and can be used. However, they are nearly redundant with the composite diagnosis endpoint for Control and Vi-TT. The Vi-PS group is the only one where bacteraemia and diagnosis diverge (11 vs 13 diagnosed).

### 6c. Fever threshold rows -- FULL COHORT

**Usable: YES.**

All fever threshold rows (>= 37.5C, >= 38.0C, >= 38.5C) use per-protocol denominators (31, 37, 35). Fever was monitored for all participants via twice-daily self-reported temperature measurements.

These are particularly valuable because:
- They provide graded severity measures of the same underlying process
- The >= 38.0C threshold is used in the post-hoc "field trial approximation" endpoint
- They are independent of blood culture results

---

## 7. Errors in Current Extraction (Jin_2017.md)

### 7a. Transcription accuracy of Table 2 values

All numerical values in the extraction match the PDF. Specifically verified:

- Primary outcome: 24/31, 13/37, 13/35 -- **Correct**
- Microbiological diagnosis: 16/31, 12/37, 9/35 -- **Correct**
- Clinical diagnosis: 8/31, 1/37, 4/35 -- **Correct**
- Fever thresholds: All three rows match PDF -- **Correct**
- S Typhi bacteraemia: 24/24, 13/13, 11/13 -- **Correct**
- Stool culture (all): 22/31, 22/37, 21/35 -- **Correct**
- Stool culture (diagnosed): 19/24, 12/13, 10/13 -- **Correct**
- Quantitative blood culture: 0.4, 0.075, 0.1 -- **Correct**
- Time-to-event values: All match PDF -- **Correct**

### 7b. Labeling/interpretation issues

**CRITICAL ERROR in Summary Statistics table (line 348 of extraction):**

> "Control bacteraemia rate: 24/31 = 77.4% -- All diagnosed had bacteraemia"

This line presents 24/31 as the "control bacteraemia rate," which implies a full-cohort denominator. While the reconstructed full-cohort rate IS 24/31 (as shown in Section 4 above), the original Table 2 reports bacteraemia as 24/24 (100%) with a diagnosed-subset denominator. The extraction's summary is technically correct in its reconstructed value but misleading because:
1. It does not flag that the raw Table 2 value is 24/24 (100%) with a diagnosed-subset denominator
2. It silently performs the reconstruction without noting this
3. The casual reader would think Jin reported 24/31 directly

**The extraction should clearly distinguish between the raw Table 2 value (24/24) and any reconstructed full-cohort value (24/31).**

### 7c. The "Breakdown by Diagnosis Type" table in the extraction

The extraction (lines 136-141) presents three rows:
- Microbiological diagnosis: 16/31, 12/37, 9/35
- Clinical diagnosis only: 8/31, 1/37, 4/35
- S. Typhi bacteraemia among diagnosed: 24/24, 13/13, 11/13

The label "S. Typhi bacteraemia among diagnosed" is an editorial addition by the extractor -- the PDF Table 2 simply says "S Typhi bacteraemia." The extractor's label is actually MORE informative and accurate than the original, correctly flagging the diagnosed-subset denominator. This is good practice but should be noted as an interpretive addition.

---

## 8. Search for Supplementary Data on Full-Cohort Bacteraemia

### What the paper references

The paper references an appendix multiple times ("appendix p 1", "appendix pp 2-6", "appendix pp 4-5", "appendix pp 8-10"). The appendix is described as "See Online for appendix" (p.2473 margin note).

### What the appendix likely contains

Based on references in the main text:
- Appendix p 1: Alternative endpoint VE estimates (fever >= 38C preceding bacteraemia)
- Appendix p 2: IgG titre data, immunogenicity figures (Figure 3B)
- Appendix pp 4-5: Solicited symptom severity comparisons
- Appendix p 6: CRP concentrations in diagnosed participants
- Appendix p 7: Unadjusted logistic regression results
- Appendix pp 8-10: Full inclusion and exclusion criteria

### Is the appendix available?

**UPDATE 2026-03-19:** The supplementary appendix (11 pages) has been obtained and is available at `input_papers/Jin et al._2017_...mmc1.pdf`. Full data from Tables S1-S4, Figures S1-S3, and inclusion/exclusion criteria have been extracted into Jin_2017.md.

### Appendix verification of estimated numerators

The audit (below) estimated fever-preceding-bacteraemia numerators from VE values as ~13/31, ~2/37, ~5/35. Supplement Table S1 confirms: **13/31, 2/37, 7/35**. The Vi-PS estimate was incorrect (~5 vs actual 7). All other estimates were correct. The extract now uses the verified supplement values.

---

## 9. Summary of Findings

### Key conclusions

1. **Blood cultures were taken daily from ALL participants** during the 14-day challenge period. The bacteraemia data is NOT limited by selective sampling.

2. **The bacteraemia row denominators (24, 13, 13) ARE the diagnosed subset**, not the full cohort. This is confirmed by matching the numbers exactly to the "Total diagnosed" row.

3. **This is not a data quality problem -- it is a logical tautology.** Bacteraemia (positive blood culture > 72h) is one of the two criteria for diagnosis. Therefore, by definition, every bacteraemic participant is diagnosed. No undiagnosed participant can be bacteraemic. The denominator restriction to diagnosed participants, while potentially confusing, does not represent loss of information.

4. **Full-cohort bacteraemia rates can be reconstructed:** Control 24/31, Vi-TT 13/37, Vi-PS 11/35. These are valid and nearly identical to the composite diagnosis rates (exactly identical for Control and Vi-TT; slightly lower for Vi-PS where 2 diagnosed participants had fever without bacteraemia).

5. **Stool shedding (22/31, 22/37, 21/35) is a valid full-cohort infection endpoint** but has low discriminating power for vaccine efficacy (59-71% across all groups).

6. **Fever threshold rows are all full-cohort** and valid for dose-response modeling.

7. **The current extraction has one significant issue:** the Summary Statistics table presents "Control bacteraemia rate: 24/31" without flagging that this is a reconstruction from the raw Table 2 value of 24/24 (diagnosed subset denominator).

### Recommended data points for dose-response modeling

| Endpoint | Control | Vi-TT | Vi-PS | Type | Recommended? |
|----------|---------|-------|-------|------|-------------|
| Composite diagnosis (primary) | 24/31 | 13/37 | 13/35 | Disease | Yes -- primary endpoint |
| Reconstructed full-cohort bacteraemia | 24/31 | 13/37 | 11/35 | Infection (systemic) | Yes, but nearly redundant with composite for Control/Vi-TT |
| Fever >= 38.0C (any duration) | 17/31 | 6/37 | 11/35 | Disease (clinical) | Yes -- approximates field definitions |
| Stool culture positive (all participants) | 22/31 | 22/37 | 21/35 | Infection (mucosal) | Conditionally -- useful as colonization endpoint, poor VE discriminator |
| Fever >= 38.0C preceding bacteraemia (post-hoc) | 13/31 | 2/37 | 7/35 | Disease (field approximation) | Yes -- best field-trial analog; **verified from Supplementary Table S1 (2026-03-19)** |

---

*Audit completed 2026-03-18*
