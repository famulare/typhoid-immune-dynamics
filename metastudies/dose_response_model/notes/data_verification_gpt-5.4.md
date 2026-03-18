# Data Verification Report

## Executive Summary

I found 7 consequential issues and several smaller metadata problems in the current plan.

Most row-level `n` and `y` values in Section 11 do match the extraction notes, but the main exceptions are important:
- `H-F-3` has a direct source conflict: `Hornick_1966.md` says `9/14` at `10^3`, while `Hornick_Snyder_1970_Part1.md` says `0/14` and states explicitly that `10^3` failed to induce disease.
- The Maryland fever endpoint is not harmonized across studies. Hornick, Gilman, and Levine do not use the same disease definition, despite the plan saying they do.
- The Oxford infection endpoint is not harmonized either. Waddington uses shedding, Darton uses `bacteremia OR stool+`, and Section 11 uses Jin stool positivity only even though the plan text elsewhere groups Jin with the broader definition.
- Section 11 bookkeeping is inconsistent: it shows 19 Oxford rows and 15 active Maryland rows, not 16, and it omits the Hornick vaccinated rows that are listed in Section 4.

In short: the likelihood table is close, but it is not yet literal enough for implementation without resolving the Hornick `10^3` conflict, endpoint-definition heterogeneity, and the Section 11 row inventory.

## Observation-by-Observation Verification Table

Notes:
- `Plan values` and `Extraction values` are shown as `dose; n; y`.
- `Approximate` means the row is directionally supported but at least one component is approximate, inferred, or not literally tabulated in the extraction note.
- Section 11 currently contains 34 active rows, not 35. Two struck rows are shown in the plan but do not enter the likelihood.

| Obs ID | Plan values | Extraction values | Match? | Notes |
|---|---|---|---|---|
| W-F-3 | `10^3; 20; 11` | `10^3; 20; 11` | Yes | `Waddington_2014_Outpatient.md` Table 2: diagnosed typhoid. |
| W-F-4 | `10^4; 16; 10` | `10^4; 16; 10` | Yes | Exact match. |
| W-F-5 | `10^5; 5; 5` | `10^5; 5; 5` | Yes | Exact match. |
| W-I-3 | `10^3; 20; 13` | `10^3; 20; 13` | Yes | Exact match for stool shedding, not bacteremia. |
| W-I-4 | `10^4; 16; 10` | `10^4; 16; 10` | Yes | Exact match for stool shedding. |
| W-I-5 | `10^5; 5; 4` | `10^5; 5; 4` | Yes | Exact match for stool shedding. |
| D-F-plac | `1.82e4; 30; 20` | `1.82e4 median; 30; 20` | Yes | `Darton_2016.md` Table 2 uses per-protocol `n=30`; median actual dose `1.82 x 10^4` across study. |
| D-F-Ty21a | `1.82e4; 30; 13` | `1.82e4 median; 30; 13` | Yes | Exact count match. |
| D-F-M01 | `1.82e4; 31; 18` | `1.82e4 median; 31; 18` | Yes | Exact count match. |
| D-I-plac | `1.82e4; 30; 26` | `1.82e4 median; 30; 26` | Yes | Exact match for `bacteraemia or stool positive`. |
| D-I-Ty21a | `1.82e4; 30; 16` | `1.82e4 median; 30; 16` | Yes | Exact match for `bacteraemia or stool positive`. |
| D-I-M01 | `1.82e4; 31; 21` | `1.82e4 median; 31; 21` | Yes | Exact match for `bacteraemia or stool positive`. |
| J-F-ctrl | `~1e4; 31; 24` | `1-5 x 10^4; 31; 24` | Approximate | Counts match; dose is only given as a range in `Jin_2017.md`. |
| J-F-ViTT | `~1e4; 37; 13` | `1-5 x 10^4; 37; 13` | Approximate | Counts match; dose approximate. GMT `562.9 EU/mL` supports plan's `563`. |
| J-F-ViPS | `~1e4; 35; 13` | `1-5 x 10^4; 35; 13` | Approximate | Counts match; dose approximate. GMT `140.5 EU/mL` supports plan's `141`. |
| J-I-ctrl | `~1e4; 31; 22` | `1-5 x 10^4; 31; 22` | Approximate | Counts match for stool positivity, but this is shedding-only, not `bacteremia OR shedding`. |
| J-I-ViTT | `~1e4; 37; 22` | `1-5 x 10^4; 37; 22` | Approximate | Same issue: exact stool-positive count, but endpoint differs from Darton's infection endpoint. |
| J-I-ViPS | `~1e4; 35; 21` | `1-5 x 10^4; 35; 21` | Approximate | Same issue. |
| G20-F-naive | `~2.5e4; 19; 12` | `1-5 x 10^4; 19; 12` | Approximate | `Gibani_2020.md` gives `12/19`; actual S. Typhi median doses vary by arm from `22.9-26.2 x 10^3`, so `~2.5e4` is reasonable but not literal for the naive arm. |
| H-F-3 | `10^3; 14; 0` | `Hornick 1970 Part 1: 10^3; 14; 0` / `Hornick 1966: 10^3; 14; 9` | No | Direct source conflict. Plan follows `Hornick_Snyder_1970_Part1.md`, not `Hornick_1966.md`. This must be resolved before fitting. |
| H-F-5 | `10^5; 116; 32` | `10^5; 116; 32` | Yes | Exact match to `Hornick_Snyder_1970_Part1.md` Table 1. |
| H-F-8 | `10^8; 9; 8` | `10^8; 9; 8` | Yes | Exact match. |
| H-F-9 | `10^9; 42; 40` | `10^9; 42; 40` | Yes | Exact match. |
| H-I-7 | `10^7; 30; 28` | `10^7; 30; 28` | Approximate | Counts match exactly, but the extraction's outcome is broader than shedding: `low-grade fever or serology or blood culture or stool >5 days`. |
| H-FgI-7 | `10^7; 28; 16` | `10^7; 28; 16` | Yes | Exact match to disease-among-infected in `Hornick_Snyder_1970_Part1.md` Table 2. |
| Gil-F-Hlo | `10^5; 14; ~9` | `10^5; 14; 61%` | Approximate | `Gilman_1977.md` gives `14` and `61%`; `~9` is back-calculated, not tabulated exactly. |
| Gil-F-Hhi | `10^5; 13; ~3` | `10^5; 13; 24%` | Approximate | `~3` is back-calculated from `24% of 13`; not literal in the extraction. |
| Gil-F-rest | `10^5; ~37; ~19` | not directly reported | Approximate | Entire row is inferred from combined controls `64/31` minus the antibody subset `14 + 13`. This row is not present in the extraction and should be labeled derived. |
| Gil-I-ctrl | `10^5; 43; 26` | `10^5; 43; 26` | Approximate | Counts match exactly, but the endpoint is `late shedding 4-30 days`, not any shedding and not Hornick's broader infection definition. |
| Lev-F-1 | `10^5; 26; 13` | `10^5; 26; 13` | Yes | Exact count match to `Levine_1976.md` Table 2. |
| Lev-F-2 | `10^5; 33; 10` | `10^5; 33; 10` | Yes | Exact match. |
| Lev-F-3 | `10^5; 22; 12` | `10^5; 22; 12` | Yes | Exact match. |
| Lev-F-4 | `10^5; 16; 4` | `10^5; 16; 4` | Yes | Exact match. |
| Lev-I-1 | `10^5; 26; 19` | `10^5; 26; 19` | Approximate | Counts match exactly, but the outcome is `stool anytime`, not Hornick's infection composite and not Gilman's `late shedding` endpoint. |

### Struck Rows Shown in Section 11

| Obs ID | Plan values | Extraction values | Match? | Notes |
|---|---|---|---|---|
| H-F-7 | removed | `10^7; 32; 16` | N/A | The extraction supports `32/16` in Table 1, but including it with `H-I-7` and `H-FgI-7` would double-count at least 30 of the same challenged volunteers. Removal is reasonable. |
| Gil-F-ctrl | replaced | `10^5; 64; 31` | N/A | The extraction supports `64/31` combined controls. Replacing it with stratified rows is reasonable, but `Gil-F-rest` is derived rather than extracted. |

## Cohort Independence Assessment

| Relationship in plan | Assessment | Evidence from extraction notes |
|---|---|---|
| Gilman 1977 controls independent of Hornick 1970 | Unclear / not directly verified | `Gilman_1977.md` confirms same Maryland prison setting and same program, but does not give challenge years. The plan's `1973-1976` timing is not documented in the listed extraction note. Independence is plausible, not proven. |
| Levine 1976 controls overlap with Hornick 1970 | Partially supported as a risk, not resolved | `Levine_1976.md` dates Trials 1-4 to `1970-1973`. That makes Trial 1 a real overlap risk with late Hornick-era work, but the extraction does not identify individual volunteers or months. |
| Hornick Table 1 and Table 5 `None` rows overlap | Likely yes, but not identical | Table 1 has `116` controls at `10^5` and `32` at `10^7`; Table 5 `None` has `104` and `30`. These look like overlapping subsets rather than independent cohorts. |
| Hornick Table 1 and Table 2 share subjects | Almost certainly largely overlapping, not perfectly identical | At `10^7`, Table 1 has `32` volunteers with `16` disease; Table 2 Quailes has `30` total with the same `16` disease. The same disease cases are probably being re-used, but 2 volunteers are unaccounted for. |
| Gibani 2020 naive arm overlaps earlier Oxford studies | Supported as independent | `Gibani_2020.md` distinguishes naive controls from re-challenge participants and says prior-study recruitment applied to re-challenge participants. Nothing in the extraction suggests the naive `19` were previously challenged. |
| Gibani 2020 re-challenge arms overlap earlier Oxford studies | Supported | `Gibani_2020.md` explicitly says re-challenge participants were recruited from earlier Oxford studies, including Waddington, Darton, and Jin. |
| Darton 2017 is the same cohort as Waddington 2014 | Strongly supported | `Darton_2017.md` explicitly says it is the same cohort and should not be used as independent dose-response data. |
| Gilman H-antibody strata cover all controls | No | `Gilman_1977.md` gives `14 + 13 = 27` controls in Table 5 versus `64` combined controls. The stratified data are a subsample, not all controls. |
| Levine H-antibody strata available at row level | No | `Levine_1976.md` reports only the pooled contrast `22% vs 48%` with `P=0.005`; no `n/y` stratified rows are extracted. |

## Cross-Study Consistency Findings

### 8. Hornick `10^3` discrepancy

- `Hornick_1966.md` records `9/14` ill at `10^3`, taken from Figure 2.
- `Hornick_Snyder_1970_Part1.md` records `0/14` in Table 1 and also states in text that the smallest dose failed to induce disease.
- Supplemental evidence from `Hornick_1967.md` also gives `0/14` and explicitly says it clarifies the ambiguity in the 1966 figure.
- The listed extraction notes do not provide a definitive explanation for why `Hornick_1966.md` says `9/14`.

Assessment:
- The plan's use of `0/14` is better supported by the later tabular summaries, but the conflict is real and should not be hand-waved away.
- This is high impact because the point is influential at the low-dose end of the Maryland curve.

### 9. Hornick Table 1 vs Table 2 at `10^7`

- Table 1 (`Hornick_Snyder_1970_Part1.md`) reports Quailes `10^7: 32 challenged, 16 disease`.
- Table 2 reports Quailes `10^7: 30 total = 16 disease + 12 infection without disease + 2 no infection`.
- The identical `16` disease count strongly suggests the disease cases are the same people in both tables.
- The extra 2 subjects in Table 1 are not explained in the extraction note.

Assessment:
- Removing `H-F-7` from the likelihood is the right bookkeeping choice if `H-I-7` and `H-FgI-7` are retained.
- However, the report should say explicitly that this choice avoids double-counting at the cost of dropping 2 subjects from the Table 1 row.

### 10. Oxford infection definition consistency

- Waddington infection row uses stool shedding only: `13/20`, `10/16`, `4/5`.
- Darton infection row uses `bacteraemia or stool positive`: `26/30`, `21/31`, `16/30`.
- Jin Section 11 infection row uses stool positivity only: `22/31`, `22/37`, `21/35`.
- `Gibani_2019.md` and `Darton_2016.md` support the broader Darton infection endpoint.
- `Jin_2017.md` supports stool-positive counts, not a directly tabulated `bacteremia OR shedding` union count for all arms.

Assessment:
- The plan is internally inconsistent here. Section 11 uses Jin shedding-only rows, but the narrative elsewhere says Darton and Jin use the broader definition.
- This matters because the infection curve is being asked to absorb different observables across Oxford studies.
- The extraction set is sufficient to keep Waddington as shedding-only, Darton as `bacteraemia OR stool+`, and Jin as shedding-only. It is not sufficient to pretend those are identical without an observation model.

### 11. Maryland fever/disease definitions across studies

- Hornick: disease is sustained oral temperature `>=103 F` for `24-36` or `36-48` hours with treatment initiation.
- Gilman: typhoid fever is presence of fever plus positive blood or stool culture, with treatment triggers that can be much lower (`>101 F for 3 days` or `>100 F for 5 days`).
- Levine: typhoid fever is acute illness with oral temperature `>=101 F` plus isolation of `S. typhi` from blood or stool.

Assessment:
- These are not the same outcome.
- The plan's statement in Section 4.2 that all Maryland studies use the same fever definition is false as written.
- This is probably the single biggest cross-study comparability problem after the Hornick `10^3` conflict.

## Missing Data Inventory

### Mentioned in extractions but not used in the plan

- Hornick incubation-period data from `Hornick_Snyder_1970_Part1.md` Table 1.
- Waddington incubation/time-to-diagnosis summaries from `Waddington_2014_Outpatient.md` Table 3.
- Darton quantitative blood culture CFU/mL at diagnosis.
- Jin quantitative blood culture CFU/mL at diagnosis.
- Darton, Jin, and Gibani time-to-event summaries for diagnosis, fever, and bacteremia.
- Waddington criterion breakdown (`bacteremia only`, `fever only`, `both`, `neither`), which is useful partial joint outcome structure.
- Gilman combined-control row `64/31`, which is currently replaced but still useful as a cross-check on the back-calculated stratified rows.
- Levine pooled H-antibody protective contrast (`22% vs 48%`), which is discussed but not operationalized.

### Metadata still needed for the model but missing or incomplete in the extraction set

- Exact administered CFU per participant or at least arm-specific actual medians for Jin and the Gibani 2020 naive S. Typhi arm.
- Literal post-vaccination anti-Vi GMT values for Darton M01ZH09 and Ty21a arms. The plan describes these qualitatively, but the listed extraction note does not tabulate them.
- Exact infection-by-fever cross-tabs for most studies if a bivariate model is to be fit.
- Exact `bacteremia OR stool+` union counts for Jin if the Oxford infection endpoint is to be harmonized on that broader definition.
- Full pre-challenge anti-Vi distributions, not just fractions above detection, for Oxford placebo/control groups.
- Selection mechanism for the Gilman H-antibody subsample (`27/64` controls measured in the extracted table).
- Stronger overlap metadata linking Hornick, Levine, and Gilman volunteers or calendar windows.

## Recommendations

1. Resolve the Hornick `10^3` row against the primary PDFs before coding the Maryland low-dose point. This is the only direct source contradiction in an active likelihood row and it is high leverage.
2. Stop treating Hornick, Gilman, and Levine fever outcomes as one literal endpoint. Either add an observation/threshold model or limit the primary fever likelihood to a more homogeneous subset.
3. Harmonize the Oxford infection endpoint explicitly. Right now the plan mixes shedding-only and `bacteremia OR stool+` without saying so.
4. Fix Section 11 bookkeeping. It currently shows 34 active rows, not 35, and it omits the Hornick vaccinated rows listed in Section 4.
5. Mark `Gil-F-Hlo`, `Gil-F-Hhi`, and especially `Gil-F-rest` as derived rows rather than exact extracted observations. If they stay in the model, use sensitivity analyses for rounding and subsample representativeness.
6. Keep `H-F-7` removed if `H-I-7` and `H-FgI-7` stay in, but document that this avoids double-counting by sacrificing 2 subjects from the Table 1 row.
7. Downgrade unsupported CoP metadata claims to notes or assumptions: Waddington's `~29%` detectable anti-Vi, Hornick `no baseline serology`, and Darton post-vaccination GMT descriptions not tabulated in the current extraction.
