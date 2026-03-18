# Data Extraction Verification Prompt

**Purpose**: You are an independent data verification agent. Your job is to meticulously check every numerical data value and cohort metadata that will enter the likelihood of a typhoid dose-response model. Errors in data extraction propagate directly into parameter estimates, so this review must be thorough and literal — check numbers against source documents, flag discrepancies, and verify cohort definitions.

---

## Context

A joint inference plan has been written for calibrating a two-outcome (infection, fever) beta-Poisson dose-response model to data from Oxford (2010s, bicarbonate delivery) and Maryland (1960s-70s, milk delivery) controlled human infection model (CHIM) studies of S. Typhi Quailes strain.

The plan specifies ~35 binomial observations that enter the likelihood. Each observation has:
- **y**: number of subjects with the outcome
- **n**: total subjects in the group
- **Dose**: challenge dose in CFU
- **Outcome**: infection (shedding/bacteremia) or fever (clinical typhoid diagnosis)
- **Medium**: bicarbonate or milk
- **CoP status**: naive, vaccinated (with GMT if available), or unknown/latent
- **Study**: which publication the data comes from
- **Cohort independence**: whether subjects overlap with another observation

Your task is to verify every one of these values against the source extraction documents, flag any discrepancy, and check that the cohort metadata is correct for the model's needs.

---

## Files to Read

### The plan (what claims to be true):
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/notes/joint_inference_plan.md`
  - Focus on **Section 4 (Data Catalog)** and **Section 11 (Summary of All Observations)**

### The source extractions (intermediate — may contain errors):
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Waddington_2014_Outpatient.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Darton_2016.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Darton_2017.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Jin_2017.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Gibani_2020.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Gibani_2019.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Hornick_1966.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Hornick_Snyder_1970.md` (Part 2)
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Hornick_Snyder_1970_Part1.md` (Part 1)
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Gilman_1977.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/extracts/Levine_1976.md`

### The primary source papers (GROUND TRUTH — check numbers against these):

**IMPORTANT**: The extractions above were AI-generated and may themselves contain errors. For every number that enters the model likelihood, you MUST verify it against the original PDF. Read the relevant tables and text in the PDFs directly.

- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Waddington et al._2014_An Outpatient, Ambulant-Design, Controlled Human Infection Model Using Escalating Doses of Salmonell.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Darton et al._2016_Using a Human Challenge Model of Infection to Measure Vaccine Efficacy A Randomised, Controlled Tri.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Darton et al._2017_Blood culture-PCR to optimise typhoid fever diagnosis after controlled human infection identifies fr.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Jin et al._2017_Efficacy and immunogenicity of a Vi-tetanus toxoid conjugate vaccine in the prevention of typhoid fe.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Gibani et al._2020_Homologous and heterologous re-challenge with Salmonella Typhi and Salmonella Paratyphi A in a rando.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Gibani et al._2019_The Impact of Vaccination and Prior Exposure on Stool Shedding of Salmonella Typhi and Salmonella Pa.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Hornick_1966_Study_of_induced_typhoid_fever_in_man._I._Evaluati.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Hornick_Snyder_1970_Typhoid_fever_pathogenesis_and_immunologic_Part_1.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Hornick_Snyder_1970_Typhoid_fever_pathogenesis_and_immunologic_Part_2.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Gilman et al._1977_Evaluation of a UDP-Glucose-4-Epimeraseless Mutant of Salmonella typhi as a Live Oral Vaccine.pdf`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/input_papers/Levine et al._1976_Attenuated, Streptomycin-Dependent Salmonella typhi Oral Vaccine Potential Deleterious Effects of L.pdf`

### Additional context:
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/notes/outcome_mapping.md`
- `/Users/famulare/git/famulare/typhoid-immune-dynamics/metastudies/dose_response_model/notes/yolo_working_model_notes.md`

---

## Verification Checklist

For EACH observation in Section 11 of the plan, verify the following. Report your findings in a structured table.

### A. Numerical Accuracy

For each (Obs ID, Study, Outcome, Dose, n, y) row:

1. **Does n match the source extraction?** Check the exact sample size against the extraction document. If the extraction reports the number differently (e.g., "n=41 total" vs "n=20 at 10³"), trace how the per-dose n was derived.

2. **Does y match the source extraction?** Check the exact number of events. Watch for:
   - Confusion between "diagnosed" (composite TD) vs "fever only" vs "bacteremia only"
   - Differences between per-protocol and intention-to-treat populations
   - Rounding or estimation (marked with ~ in the plan)

3. **Does the dose match?** Verify the exact dose. Watch for:
   - "1-5 × 10⁴" ranges vs point estimates
   - Whether the reported dose is the target or the actual measured dose (e.g., Darton reports median 1.82 × 10⁴)

4. **Is the outcome definition correctly mapped?** Verify that what the plan calls "infection" matches the extraction's shedding or broad infection definition, and what the plan calls "fever" matches the extraction's TD composite or disease definition.

### B. Cohort Definitions and Independence

For each observation:

5. **Is the CoP status correctly described?** Verify:
   - "Naive" groups: Were they actually screened for prior typhoid? What fraction had detectable baseline anti-Vi?
   - "Vaccinated" groups: What vaccine? What was the post-vaccination GMT? Is the GMT from the extraction?
   - "Maryland mixture": Is there any baseline serology reported that the plan ignores?

6. **Is the cohort independence claim correct?** Specifically check:
   - Do the Gilman 1977 controls overlap with the Hornick 1970 controls? (Different time periods? Same prison?)
   - Do the Levine 1976 controls overlap with Hornick 1970? (Levine Trial 1 was in 1970)
   - Do the Hornick Table 1 (Part 1) and Table 5 (Part 2) share subjects? (The plan claims the "None" row in Table 5 overlaps with Table 1)
   - Do the Gibani 2020 naive arms overlap with any earlier Oxford study?
   - Is Darton 2017 really the same cohort as Waddington 2014? (The plan excludes it on this basis)

7. **For the H-antibody stratified data (Gilman, Levine)**:
   - How many subjects had H-antibody measured? Was it all controls or a subsample?
   - Are the y values exact or back-calculated from percentages? (The plan marks some with ~)
   - Do the stratified counts sum to the pooled total?

### C. Cross-Study Consistency

8. **Hornick 10³ discrepancy**: The plan uses 0/14 from Hornick 1970 Part 1, but Hornick 1966 reports 9/14 at the same dose.
   - What exactly does each extraction say?
   - Is there any explanation in the extractions for the discrepancy?
   - Does Hornick 1970 Part 1 explicitly report 0/14, or is this inferred?

9. **Hornick Table 1 vs Table 2 at 10⁷**: Table 1 reports N=32 with 16 disease. Table 2 reports N=30 (Quailes only) with 16 disease and 12 infection-without-disease.
   - Are the 16 disease cases the same people in both tables?
   - Where do the extra 2 subjects in Table 1 come from?
   - The plan uses H-F-7 as REMOVED (double-counting). Is this correct?

10. **Oxford infection definition consistency**: The plan notes that Waddington uses "shedding" while Darton/Jin use "bacteremia OR shedding."
    - Verify the exact infection definition for each Oxford study from the extractions
    - What are the shedding-only rates for Darton and Jin (if reported)?
    - How large is the gap between shedding-only and the broader definition?

11. **Maryland outcome definitions across studies**: Hornick, Gilman, and Levine all used the Maryland CHIM protocol, but the fever/disease definitions may differ slightly.
    - Hornick: "temperature exceeded 103°F for over 24 to 36 hours"
    - Gilman: "presence of fever AND blood or stool culture positive"
    - Levine: "acute illness with oral temperature ≥101°F accompanied by isolation of S. typhi"
    - Are these really the same outcome? Flag any discrepancies.

### D. Missing Data and Metadata

12. **What data mentioned in the extractions is NOT used in the plan?** Flag anything potentially useful that was left out. In particular:
    - Incubation period data (Hornick Table 1, Waddington Table 3)
    - Quantitative bacteremia (CFU/mL) from Darton and Jin
    - Time-to-event data from any study
    - Any other immunity correlate data beyond anti-Vi IgG

13. **What metadata is needed for the model but missing from the extractions?**
    - Exact actual doses administered (vs target ranges)
    - Individual-level titer data for Oxford vaccinated arms
    - Cross-tabulated infection × fever counts (for bivariate modeling)
    - Pre-challenge anti-Vi IgG distribution parameters (not just fraction detectable)

---

## Output Format

Please produce your report as a single markdown document with:

1. **Executive Summary**: How many discrepancies found? How many are consequential for the model?

2. **Observation-by-Observation Verification Table**: One row per observation in Section 11, with columns:
   - Obs ID
   - Plan values (n, y, dose)
   - Extraction values (n, y, dose)
   - Match? (Yes / No / Approximate)
   - Notes on any discrepancy

3. **Cohort Independence Assessment**: For each claimed independence/overlap relationship, your assessment with evidence.

4. **Cross-Study Consistency Findings**: Results for items 8-11 above.

5. **Missing Data Inventory**: Items 12-13.

6. **Recommendations**: Ranked list of issues that should be resolved before Stan implementation, from most to least consequential.

Be pedantic. Be literal. Check every number. If a value is marked with ~ in the plan, verify what the actual value should be and whether the approximation is reasonable. If you cannot verify a value because the extraction doesn't contain it, say so explicitly.
