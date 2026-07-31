# Gilman et al. 1977 - Evaluation of a UDP-Glucose-4-Epimeraseless Mutant of Salmonella typhi as a Live Oral Vaccine

## Citation
Gilman RH, Hornick RB, Woodward WE, DuPont HL, Snyder MJ, Levine MM, Libonati JP. Evaluation of a UDP-Glucose-4-Epimeraseless Mutant of Salmonella typhi as a Live Oral Vaccine. J Infect Dis. 1977 Dec;136(6):717-723.

**DOI/Source**: The Journal of Infectious Diseases, Vol. 136, No. 6, December 1977

## Study Overview
- **Primary question**: Evaluate safety and protective efficacy of the gal E mutant Ty21a oral typhoid vaccine strain grown with and without galactose
- **Study type**: Vaccine efficacy trial with subsequent virulent challenge (CHIM)
- **Institution**: Division of Infectious Diseases, University of Maryland School of Medicine
- **Funding**: Army grant no. 49-193-MD-2867 from the U.S. Department of Defense
- **Population**: Adult male inmates at Maryland House of Correction, Jessup, Maryland

**[ASSISTANT-PROPOSED]** This is primarily a vaccine efficacy study rather than a dose-response study. Challenge was performed at a single dose (10^5), so limited utility for dose-response modeling, but provides valuable control group data at that dose level.

## Challenge Protocol

### Strain
- **Name**: Quailes strain of *S. typhi*
- **Virulence**: Virulent wild-type strain (same as other Maryland CHIM studies)

### Dose & Delivery
- **Challenge dose**: 10^5 organisms (single dose level)
- **Delivery vehicle**: 45 ml of milk
- **Acid neutralization**: Not explicitly stated for challenge (NaHCO3 used for vaccination)
- **Timing**: 5-9 weeks following completion of vaccination series

### Vaccination Protocol (for vaccine groups)
- **Vaccine strain**: Ty21a (gal E mutant of *S. typhi*)
- **Vaccine doses**: 5-8 doses containing 3-10 x 10^10 viable bacteria per dose
- **Dosing schedule**: 3-4 day intervals over a 4-week period
- **Delivery**: Orally in 45 ml of milk, 5 min after 2 g NaHCO3 in 60 ml distilled water
- **Trials 1 & 3**: Vaccine grown in brain-heart infusion broth (BHIB) with 0.1% galactose
- **Trial 2**: Vaccine grown on BHI agar without galactose

## Subject Characteristics

### Inclusion/Exclusion
- Healthy adult male inmates
- Accepted for challenge if medical history, physical examination, CBC, urinalysis, blood chemistries (including liver function tests), ECG, chest X-ray, and oral cholecystogram were normal
- Informed consent obtained; no coercion; free to withdraw
- Study protocols reviewed and approved by University of Maryland Committee on Human Experimentation

### Sample Sizes
| Group | Vaccinated (N) | Challenged (N) | Controls Challenged (N) |
|-------|----------------|----------------|------------------------|
| Trials 1 & 3 (with galactose) | 100 (35 + 65) | 28 | 43 |
| Trial 2 (without galactose) | 56 | 27 | 21 |

**[ASSISTANT-PROPOSED]** Not all vaccinated individuals were subsequently challenged; numbers differ between vaccination phase (155 total vaccinated) and challenge phase (55 vaccinees + 64 controls = 119 challenged).

## Baseline Immunity

### Serologic Testing
- Sera collected prior to vaccination and biweekly for 8 weeks
- Antibodies to O, H, and Vi antigens measured
- Additional sera collected before challenge and biweekly for 8 weeks post-challenge

### Pre-existing Antibody Observations
- Titer of antibody to O antigen was rarely elevated prior to vaccination (p.720)
- H antibody base-line titer >= 1:20 associated with reduced seroconversion to vaccination (Table 3, p.720)
- In unvaccinated controls, pre-existing H antibody (>= 1:20) correlated with lower attack rate: 24% vs 61% in those without H antibody (P = 0.02) (Table 5, p.721)
- **O and Vi antibody did NOT correlate with protection, in either controls or vaccinees** (p.721,
  verbatim): "in vaccinees the presence of H antibody was not associated with protection from typhoid
  fever, nor was the presence of O or Vi antibody correlated with protection of vaccinees or controls."
  Added 2026-07-31 -- this sentence was previously absent from the extract. It is the most
  decision-relevant line in the paper for the calibration, because the model's CoP axis IS anti-Vi.
  Levine 1976 p.427 reports the same null for O and Vi independently. Neither paper publishes the
  O/Vi cross-tabulations, so the null cannot be quantified and is plausibly small-N plus a compressed
  baseline titre range rather than a biological null.

**[ASSISTANT-PROPOSED]** The H antibody-protection correlation in controls is notable but was not observed in vaccinees.

**[AUTHORS, p.722 -- do not attribute the above reading to them]** Gilman et al. draw the OPPOSITE
inference from the same discordance, and use it to doubt the effect: "Prior antibody to H antigen in
unvaccinated control volunteers was directly correlated with subsequent protection... This effect,
however, may not have been specifically due to H antibody since vaccinated volunteers with high H
titers were not preferentially protected. Indeed, humoral antibodies may be quite irrelevant to
immunity to typhoid fever, which could be due to local or cell-mediated phenomena."

**[SCOPE NOTE, Mike 2026-07-31]** Gilman's objection is about humoral antibody as a MECHANISTIC
correlate, and on that he is right. It does not bear on whether a titre is a useful PREDICTIVE
correlate, which is the only sense the dose-response model uses. Keep the two claims separate:
the authors' mechanistic scepticism is not an argument against a predictive CoP.

## Outcome Definitions

### Typhoid Fever (Clinical Disease)
- **Definition**: Presence of fever AND blood or stool culture positive for *S. typhi*
- **Treatment threshold** (for chloramphenicol initiation):
  - Fever >103 F (39.4 C) for 1 day, OR
  - Fever >101 F (38.3 C) for 3 consecutive days, OR
  - Fever >100 F (37.8 C) for 5 consecutive days

### Fecal Excretion
- Stool and rectal swabs cultured using Selenite F enrichment broth
- Plated onto salmonella-shigella agar and MacConkey's agar
- *S. typhi* identified by standard methods
- Excretion categorized as: 0-3 days and 4-30 days post-challenge

## Dose-Outcome Data Tables

### Table 4 (p.720): Challenge Outcomes - Typhoid Fever

| Group | N | Typhoid Fever (%) | P-value vs Controls |
|-------|---|-------------------|---------------------|
| **Trials 1 & 3 (vaccine + galactose)** | | | |
| Vaccinees | 28 | 2 (7%) | P = 0.0002 |
| Controls | 43 | 23 (53%) | - |
| **Trial 2 (vaccine - galactose)** | | | |
| Vaccinees | 27 | 5 (19%) | NS |
| Controls | 21 | 8 (38%) | - |

**Vaccine Efficacy (Trials 1 & 3)**: 87% (calculated as 1 - 7%/53%)

### Table 4 (p.720): Challenge Outcomes - Fecal Excretion of Virulent S. typhi

| Group | N | Excretion 0-3 days (%) | Excretion 4-30 days (%) |
|-------|---|------------------------|------------------------|
| **Trials 1 & 3** | | | |
| Vaccinees | 28 | 10 (36%) NS | 3 (11%) P = 0.00009 |
| Controls | 43 | 21 (49%) | 26 (60%) |
| **Trial 2** | | | |
| Vaccinees | 27 | 5 (19%) NS | 5 (19%) NS |
| Controls | 21 | 5 (24%) | 8 (38%) |

### Control Group Data for Dose-Response Modeling

**[CRITICAL DATA]** Unvaccinated control data at 10^5 dose:

| Source | N | Typhoid Fever | Attack Rate |
|--------|---|---------------|-------------|
| Trials 1 & 3 controls | 43 | 23 | 53% |
| Trial 2 controls | 21 | 8 | 38% |
| **Combined controls** | 64 | 31 | 48% |

**[ASSISTANT-PROPOSED]** The combined control attack rate of 48% (31/64) at 10^5 dose is higher than the 24-27% reported in Hornick 1970 (which had N=104, 28 cases = 27%). This may reflect natural variability, differences in control populations, or batch effects. Should be compared with other 10^5 dose data points.

### Table 5 (p.721): H Antibody and Protection in Controls

Table 5 cross-tabulates clinical outcome (Well/Ill) by baseline H-antibody status in unvaccinated controls:

| Clinical Outcome | H Ab <1:20 | H Ab ≥1:20 |
|------------------|------------|------------|
| Well | 14 | 13 |
| Ill | 22 (61%)* | 4 (24%)* |
| **Total** | **36** | **17** |

*Difference between attack rates significant at P = 0.02

**Denominators**: 36 controls had H Ab <1:20; 17 had H Ab ≥1:20. Total with H-antibody data: 53 of 64 controls (83%). The 11 controls without H-antibody data are unaccounted for in this table (likely missing baseline sera).

### Table 6 (p.721): Seroconversion After Challenge

| Group | N | O Antibody Seroconversion | H Antibody Seroconversion |
|-------|---|---------------------------|---------------------------|
| **With clinical disease** | | | |
| All vaccinees | 7 | 0% (P = 0.003 vs controls) | 14% (P = 0.05 vs controls) |
| All controls | 30 | 70% | 63% |
| **Without clinical disease** | | | |
| All vaccinees | 48 | 0% (NS) | 8% (P = 0.009 vs controls) |
| All controls | 32 | 6% | 34% |

**[ASSISTANT-PROPOSED]** The suppressed O antibody response in vaccinees who developed disease is an unexpected finding that authors attribute to possible oral tolerance mechanisms.

## Vaccine-Related Data (Not Directly Relevant to Dose-Response)

### Table 2 (p.720): Humoral Antibody Response to Vaccination

| Vaccine Trial | N | O Seroconversion | H Seroconversion | Vi Seroconversion |
|---------------|---|------------------|------------------|-------------------|
| Trials 1 & 3 (+ galactose) | 96 | 17% (P = 0.04) | 51% (NS) | 2% (NS) |
| Trial 2 (- galactose) | 53 | 4% | 60% | 0% |

### Table 3 (p.720): Effect of Prior H Antibody on Vaccine Seroconversion

| Vaccine Trial | Seroconversion if H < 1:20 | Seroconversion if H >= 1:20 |
|---------------|---------------------------|----------------------------|
| Trials 1 & 3 | 73% (P = 0.0001) | 31% |
| Trial 2 | 74% (P = 0.01) | 33% |

## Data Quality Notes

### Strengths
- Large control group (N=64 combined) at 10^5 dose
- Well-documented disease definition with explicit temperature thresholds
- Pre-challenge antibody status recorded
- Concurrent unvaccinated controls challenged alongside vaccinees
  - NOT randomized: the word "random" does not appear in the paper. p.719 says only
    "all available vaccinees and a comparable number of unvaccinated men serving as
    controls were fed 10^5 virulent S. typhi". Corrected 2026-07-31 (was "Prospective
    randomized design", which was ASSISTANT-INFERRED and unsupported).
- Same research group/methods as other Maryland CHIM studies (enables cross-study comparison)

### Limitations
- Single challenge dose (10^5) limits dose-response modeling utility
- Two control groups from different trials with different attack rates (53% vs 38%)
- Not all vaccinated individuals challenged (selection bias possible)
- Challenge timing varied (5-9 weeks post-vaccination)
- Volunteer population (healthy adult males) not representative of endemic populations

### Page References
- Study design and vaccine preparation: p.717-718
- Challenge protocol and disease definition: p.719
- Adverse reactions: p.719
- Table 1 (fecal excretion of vaccine): p.719
- Table 2 (antibody response to vaccination): p.720
- Table 3 (H antibody effect on seroconversion): p.720
- Table 4 (challenge outcomes): p.720
- Table 5 (H antibody and protection): p.721
- Table 6 (seroconversion after challenge): p.721
- Discussion of protective mechanisms: p.721-722

## Cross-References
- **Maryland CHIM program**: This is part of the larger University of Maryland human challenge studies; same Quailes strain used
- **Hornick 1970 NEJM**: Reference [3] in this paper; provides earlier dose-response data at multiple dose levels
- **DuPont 1970**: Reference [4]; streptomycin-dependent vaccine evaluation
- **Germanier & Furer 1975**: Reference [6]; original characterization of Ty21a strain
- **Control comparison**: Attack rate of 53% (Trials 1&3) and 38% (Trial 2) vs 24-27% in Hornick 1970 at same 10^5 dose

## Fit Role Assessment

**[ASSISTANT-PROPOSED]** **SUPPORT** - Limited utility for primary dose-response modeling due to single challenge dose. However, provides valuable data points:
1. Additional control group data at 10^5 dose (N=64)
2. An H-antibody/protection ASSOCIATION in unvaccinated controls (not "H antibody-mediated
   protection": post-hoc, RR 2.60 with 95% CI 1.06-6.36, Fisher p = 0.018, and the authors
   explicitly decline the causal reading. The analysis was prompted externally -- p.717 thanks
   "Dr. R. Edelman ... for drawing our attention to H antibody and disease protection.")
3. Suppressed antibody response in vaccinees suggests oral tolerance
4. Cross-validation data for Maryland CHIM program consistency

**[OPEN]** The higher attack rate in controls (48% combined) compared to Hornick 1970 (27% at 10^5) needs investigation. Possible explanations:
- Random variation (Hornick N=104 vs Gilman N=64)
- Temporal changes in challenge strain
- Different volunteer cohort characteristics
- Methodological differences

## Key Extractions for Calibration

### Unvaccinated control dose-response (10^5 only):
| Study Arm | Dose | N | Disease | Rate |
|-----------|------|---|---------|------|
| Trials 1&3 controls | 10^5 | 43 | 23 | 53% |
| Trial 2 controls | 10^5 | 21 | 8 | 38% |
| Combined | 10^5 | 64 | 31 | 48% |

### Antibody-stratified attack rates (controls only, 53/64 with H Ab data):
| H Antibody Status | N | Disease | Rate |
|-------------------|---|---------|------|
| < 1:20 | 36 | 22 | 61% |
| >= 1:20 | 17 | 4 | 24% |

**Note**: Case counts are exact from Table 5 (Well + Ill rows give both numerator and denominator). 53/64 controls had baseline H-antibody measured; 11 are missing.
