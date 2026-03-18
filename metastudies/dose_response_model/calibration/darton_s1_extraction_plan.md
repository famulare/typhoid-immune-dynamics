# Darton S1 Individual-Level Data Extraction Plan

**Date**: 2026-03-18
**Source**: `input_papers/Darton et al._2016_...S1 data.xlsx`
**Purpose**: Extract individual-level data from Darton 2016 supplementary dataset to resolve definition inconsistencies and enable richer modeling.

---

## What the S1 Data Contains

### Endpoints sheet (n=92, one row per participant)
| Column | Meaning | Use |
|--------|---------|-----|
| `Group` | Vaccine arm: M01ZH09, Placebo, Ty21a | Stratification |
| `PPP Dx` | Per-protocol diagnosis (1/0) | TD composite |
| `TDCriteria` | How diagnosed: HiTemp, BCafterDay7, BCplusSympbeforeDay7, Day14 (=not diagnosed), Symptoms (=withdrawn) | Diagnosis pathway |
| `TTPosBC` | Time to positive blood culture (hours; 0 = never) | Bacteremia timing |
| `T37`...`T39` | Binary: ever reached fever threshold ≥37°C...≥39°C | Fever at multiple thresholds |
| `PosBCorStool` | Positive blood culture OR stool (1/0) | Broad infection composite |

### Microbiology results sheet (n=2276 samples)
| Column | Meaning | Use |
|--------|---------|-----|
| `SubjectID` | Participant number | Link to other sheets |
| `Specimen` | FAECES or BLOOD | Filter for stool vs blood |
| `STyphi` | S. Typhi isolated (1/0) | Positive result |
| `Group` | Vaccine arm | Stratification |
| `SampleDayNo` | Day relative to challenge | Timing |

### ELISA data sheet (n=92, one row per participant)
| Column | Meaning | Use |
|--------|---------|-----|
| `Group` | Vaccine arm | Stratification |
| `Vi IgG Day-28` | Anti-Vi IgG pre-vaccination (EU/mL) | Baseline CoP |
| `Vi IgG Day 0` | Anti-Vi IgG pre-challenge (EU/mL) | Challenge-time CoP |

---

## Extraction Script Specification

### Output 1: `darton_individual_endpoints.csv`

One row per participant. Columns:

| Column | Source | Logic |
|--------|--------|-------|
| `subject_id` | Microbiology `SubjectID` (need to map row order from Endpoints to SubjectID) | See note below |
| `group` | Endpoints `Group` | Direct |
| `ppp_dx` | Endpoints `PPP Dx` | 1 = per-protocol diagnosed |
| `td_criteria` | Endpoints `TDCriteria` | Diagnosis pathway |
| `fever_td` | Endpoints `PPP Dx` | 1 = composite TD (≥38°C/12h OR bacteremia) |
| `fever_38` | Endpoints `T38` | 1 = ever ≥38.0°C |
| `fever_385` | Endpoints `T38.5` | 1 = ever ≥38.5°C |
| `fever_39` | Endpoints `T39` | 1 = ever ≥39.0°C |
| `bacteremia` | Endpoints: `TTPosBC > 0` | 1 = any positive blood culture |
| `stool_positive` | Microbiology: any FAECES row with `STyphi=1` for this subject | 1 = any positive stool culture (shedding) |
| `bact_or_stool` | Endpoints `PosBCorStool` | 1 = bacteremia OR stool positive |
| `vi_igg_baseline` | ELISA `Vi IgG Day-28` | Pre-vaccination anti-Vi (EU/mL) |
| `vi_igg_prechallenge` | ELISA `Vi IgG Day 0` | Pre-challenge anti-Vi (EU/mL) |

**Subject ID mapping issue**: The Endpoints sheet does NOT have SubjectID — it's ordered by group with one row per participant. The Microbiology sheet has SubjectID. We need to determine whether the Endpoints rows are in the same order as subjects appear in other sheets. Strategy: count subjects per group in both sheets; if they match, the row order within each group is the mapping. If not, we need a cross-reference.

**Alternative simpler approach**: Since we mostly need group-level summaries, we can compute shedding per subject from Microbiology (grouping by SubjectID and Group), then merge group-level counts. Individual-level matching is only needed for cross-tabulation and individual CoP analysis.

### Output 2: `darton_group_summaries.csv`

One row per (group, endpoint) combination. This is what directly enters the plan's likelihood.

| Column | Values |
|--------|--------|
| `group` | M01ZH09, Placebo, Ty21a |
| `n_ppp` | Per-protocol N |
| `n_td` | Composite TD diagnosed |
| `n_fever_38` | Fever ≥38°C |
| `n_fever_39` | Fever ≥39°C |
| `n_bacteremia` | Any positive blood culture |
| `n_stool_positive` | Any positive stool culture (SHEDDING-ONLY) |
| `n_bact_or_stool` | Bacteremia OR stool positive |
| `mean_vi_igg_day0` | Geometric mean anti-Vi IgG at challenge |
| `frac_vi_detectable` | Fraction with Vi IgG > 7.4 EU/mL at baseline |

### Output 3: `darton_cross_tabulation.csv`

For each group, the 2×2 cross-tabulation of (stool_positive × fever_td). This enables bivariate modeling per Reviewer 2's suggestion.

| Group | Stool+ & Fever+ | Stool+ & Fever- | Stool- & Fever+ | Stool- & Fever- |
|-------|-----------------|-----------------|-----------------|-----------------|

---

## What This Resolves

### 1. Darton shedding-only rates (CRITICAL)
Currently the plan uses "bacteremia OR stool positive" (26/30 = 87% for placebo) because shedding-only wasn't available from the publication tables. The Microbiology sheet allows direct computation of shedding-only by group, resolving the Oxford infection definition inconsistency.

**Expected result**: Placebo shedding-only should be ~70-80% (lower than the 87% bact-OR-shed), bringing it closer to Jin's 71% at a similar dose.

### 2. Individual-level anti-Vi IgG for CoP mapping
The ELISA sheet has per-subject Vi IgG at both Day -28 (pre-vaccination) and Day 0 (pre-challenge). This enables:
- Direct computation of the CoP distribution within each group
- Individual-level logistic regression of outcome ~ CoP (replacing the aggregate HR = 0.29)
- Proper integration over the CoP distribution (addressing Jensen's inequality concern from Reviewer 2)

**Key values to extract**:
- Pre-challenge (Day 0) Vi IgG: this is the CoP at time of challenge
- Fraction with Day -28 Vi IgG > 7.4 EU/mL: validates the "40% of placebo had detectable baseline anti-Vi" claim

### 3. Cross-tabulated infection × fever (for bivariate modeling)
Individual-level data enables the 2×2 table that Reviewer 2 requested. We can compute P(fever AND infected), P(fever AND not infected), etc. directly.

### 4. Fever at multiple thresholds (for φ calibration)
Individual-level fever data at 37°C, 37.5°C, 38°C, 38.5°C, 39°C enables direct calibration of the definition sensitivity parameter φ — not just from group-level percentages but from individual-level threshold crossing.

---

## Implementation Notes

### Language: R (preferred per CLAUDE.md)

The script should:
1. Read the xlsx using `readxl::read_excel()`
2. Process each sheet
3. Handle the Subject ID mapping carefully
4. Output clean CSVs to `metastudies/dose_response_model/analysis_data/`
5. Print verification checks: do the computed group-level summaries match the published Table 2 values?

### Verification checks (MUST PASS before trusting outputs)

| Check | Expected | Source |
|-------|----------|--------|
| Placebo per-protocol N | 30 | Darton Table 2 |
| M01ZH09 per-protocol N | 31 | Darton Table 2 |
| Ty21a per-protocol N | 30 | Darton Table 2 |
| Placebo TD diagnosed | 20 | Darton Table 2 |
| Placebo bact-or-stool | 26 | Darton Table 2 |
| Placebo fever ≥38°C | 18 | Darton Table 2 |
| Placebo fever ≥39°C | 9 | Darton Table 2 |
| M01ZH09 TD diagnosed | 18 | Darton Table 2 |
| Ty21a TD diagnosed | 13 | Darton Table 2 |
| Placebo baseline Vi detectable (>7.4) | ~40% (12/30) | Darton p.13 |

### Edge cases to handle
- One M01ZH09 participant was "withdrawn at investigator's discretion due to unmanageable symptoms" (TDCriteria = 'Symptoms', PPP Dx = None). Exclude from per-protocol analysis.
- ELISA Vi IgG values of 3.7 appear to be below-detection-limit (the LOD is 7.4 EU/mL per Darton p.13). Treat 3.7 as censored at LOD/2.
- Some `TTT38` values are Excel formulas (`=F2-12`) rather than computed values. The script should use computed values from the xlsx, not parse formulas.

---

## Priority

This extraction should be done BEFORE Stan implementation. The shedding-only rates and individual-level anti-Vi data directly affect the likelihood specification. The current plan has placeholder values for Darton infection (flagged as definition mismatch) and aggregate-level CoP mapping — both are superseded by this individual-level data.
