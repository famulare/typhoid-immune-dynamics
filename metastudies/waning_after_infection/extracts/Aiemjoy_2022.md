# Aiemjoy et al. 2022 - Estimating typhoid incidence from community-based serosurveys: a multicohort study

## Citation
Aiemjoy K, Seidman JC, Saha S, et al. Estimating typhoid incidence from community-based serosurveys: a multicohort study. Lancet Microbe 2022; 3: e578–87. https://doi.org/10.1016/S2666-5247(22)00114-8

## Study Overview
- **Primary question**: Can longitudinal antibody kinetics following blood culture-confirmed enteric fever be used to estimate population-level seroincidence from cross-sectional serosurveys?
- **Study type**: Longitudinal cohort of confirmed enteric fever cases + cross-sectional population serosurveys
- **Institutions**: Stanford, UC Davis, Sabin Vaccine Institute, Child Health Research Foundation (Bangladesh), Dhulikhel Hospital (Nepal), Aga Khan University (Pakistan), Massachusetts General Hospital, International Vaccine Institute
- **Population**: Blood culture-confirmed enteric fever patients (ages 0–60+) and population-based serosurvey participants (ages 0–25)
- **Countries**: Bangladesh, Nepal, Pakistan, Ghana
- **Enrollment period**: May 2016 – Feb 2021
- **Funding**: Bill & Melinda Gates Foundation (INV-000572); Vi antibody work supported by NIH/NIAID R01AI134814

**[CRITICAL]** This paper's primary purpose is seroincidence estimation, not waning characterization per se. However, the longitudinal antibody kinetics model and fitted parameters are directly relevant to understanding post-infection antibody boosting and waning. The paper provides fitted waning parameters for HlyE IgA, HlyE IgG, LPS IgA, LPS IgG, and (in supplement, Nepal only) Vi IgG.

## Subject Characteristics

### Longitudinal Case Cohort (Table 1, p.e581)
- **N total**: 1420 blood culture-confirmed enteric fever patients
  - Bangladesh: 407
  - Nepal: 543
  - Pakistan: 399
  - Ghana: 71
- **Age**: Median 9.0 years (IQR 4.5–19.7)
- **Sex**: 784 male (55.2%), 635 female (44.7%)
- **Serovar**: 1267 S. Typhi (89.2%), 153 S. Paratyphi A (10.8%)
- **Vaccine status**: 1395 (98.2%) not vaccinated; 6 (0.4%) received TCV; 10 (0.7%) received other typhoid vaccine
- **Follow-up**: Median 382.0 days (IQR 94.0–696.0)
- **Total longitudinal samples**: 4126
- **Hospitalized**: 458 (32.3%) yes, 939 (66.1%) no

### Age Distribution by Country
| Country | <5 years | 5–15 years | >16 years | Median age |
|---------|----------|------------|-----------|------------|
| Bangladesh | 179 (44.0%) | 227 (55.8%) | 1 (0.2%) | 5.4 (3.2–8.0) |
| Nepal | 12 (2.2%) | 127 (23.4%) | 404 (74.4%) | 20.9 (15.7–26.4) |
| Pakistan | 181 (45.4%) | 152 (38.1%) | 66 (16.5%) | 5.3 (3.0–12.0) |
| Ghana | 11 (15.5%) | 52 (73.2%) | 8 (11.3%) | 8.0 (6.0–12.0) |

**[NOTE]** Age distribution differs markedly across countries — Bangladesh and Pakistan are predominantly children <15 years; Nepal is predominantly adults >16 years. This affects interpretation of age-stratified parameters.

### Cross-Sectional Serosurvey Participants (Table 2, p.e584)
- **N total**: 1808 individuals aged 0–25 years
  - Dhaka, Bangladesh: 401
  - Kathmandu, Nepal: 353
  - Kavrepalanchok, Nepal: 481
  - AKU, Karachi, Pakistan: 294
  - KGH, Karachi, Pakistan: 200
  - Agogo, Ghana: 79

## Sampling Schedule

### Case Follow-up (p.e579)
- **SEAP sites (Bangladesh, Nepal, Pakistan)**: Plasma at enrollment; dried blood spot (DBS) at 28 days, 3 months, 6 months, 12 months, 18 months post-enrollment
- **SETA-Ghana**: Plasma at days 3–7, 28–30, 90, 180, 270, 360
- Retrospective cases enrolled with DBS at enrollment then followed on same schedule
- Median samples per participant: 3.0 (IQR 2.0–4.0)

## Antibody Targets Measured

### Antigens and Isotypes
| Antigen | Isotypes | ELISA dilutions | Notes |
|---------|----------|-----------------|-------|
| HlyE (Hemolysin E) | IgA, IgG | 1:500, 1:5000 | Cytotoxic pore-forming toxin; present in S. Typhi and S. Paratyphi A but rare in other serovars |
| LPS (S. Typhi lipopolysaccharide) | IgA, IgG | 1:1000, 1:5000 | Major outer membrane component; potential cross-reactivity with non-typhoidal Salmonella via shared O12 antigen |
| Vi (capsular polysaccharide) | IgG only | 1:100 | Nepal and Ghana only; cannot distinguish natural infection from Vi-based vaccination |

**[CRITICAL]** Anti-Vi was only measured in the Nepal sub-cohort (and Ghana). Vi results are reported only in the supplement (Fig S11, Table S1 overall row). The main seroincidence analysis used only HlyE IgA and IgG because LPS showed cross-reactivity with non-typhoidal Salmonella and lower responses in S. Paratyphi A cases.

### Laboratory Methods (Supplement p.2)
- Kinetic ELISAs on plasma and eluted DBS
- DBS: Two filled filter paper protrusions submerged in 133 µL 1xPBS 0.05% Tween, eluate assumed ~1:10 dilution
- Plates coated with: S. Typhi LPS from strain Ty21a (2.5 µg/mL), purified HlyE (1 µg/mL), Vi antigen (Sanofi Pasteur, 2 µg/mL)
- Results normalized to plate standard (plasma pool for LPS/Vi; chimeric monoclonal for HlyE), multiplied by 100, reported as ELISA units

## Antibody Kinetics Model

### Model Structure (p.e580; refs 20, 21; see Teunis et al. 2016 for derivation)

**Two-phase model with exponential rise, peak, then power-function decay.**

The model fits five parameters per antigen-isotype-age stratum:
1. **y0** — baseline antibody response (pre-infection level)
2. **y1** — peak antibody response
3. **t1** — time to peak antibody response (days since fever onset)
4. **ν** (nu) — decay rate parameter (1/days in Teunis notation; reported as α in Aiemjoy Table S1, in EU/year)
5. **r** — decay shape factor (power-function exponent)

#### Exact Functional Forms (from Teunis et al. 2016, Epidemics 16:33–39)

The underlying within-host model (Teunis Eq. 13) has two phases:

**Infection/colonization episode** (t ≤ t1):
- Pathogen: b'(t) = μ₀b(t) − cy(t)
- Antibody: y'(t) = μy(t)  (exponential growth)
- Solution: y(t) = y0·exp(μt), where μ = (1/t1)·log(y1/y0)

**Waning immunity episode** (t > t1):
- Pathogen cleared: b(t) = 0
- Antibody decay ODE: **y'(t) = −ν·y(t)^r**  (Teunis Eq. 13, right side)

This ODE has the closed-form solution (Teunis Eq. 14):

> **y(t > t1) = y1 · (1 + (r−1)·y1^(r−1)·ν·(t − t1))^(−1/(r−1))**

This is equivalent to a power-function decay (Teunis Eq. 11):

> **y(τ) = y1·(1 + β·τ)^(−α)**

where τ = t − t1, with the mapping: shape α = 1/(r−1), scale β = (r−1)·y1^(r−1)·ν.

#### Physical Interpretation (Teunis 2016, Sections 2.2–2.3)
The power-function decay arises from assuming the **distribution of decay rates across many antibody production sites follows a Gamma distribution** (shape α, scale β). This is biologically motivated: there are many populations of antibody-producing cells with heterogeneous lifetimes, and their aggregate decay is non-exponential. Key implications:
- When **r > 1**: initial decay is faster than exponential, followed by a prolonged period of slower-than-exponential decay (a "long tail")
- When **r → 1**: recovers simple exponential decay (single decay rate for all production sites)
- The model predicts **much more persistent antibody** at long times than exponential models — Teunis Fig. 4 shows that for pertussis IgG-PT, 80% of subjects still have concentrations >5 IU/mL at **30 years** post-infection under power-function decay, vs. ~0% under exponential

**[CRITICAL]** All shape factors r in the Aiemjoy fits are >1 (range 1.25–2.89), confirming non-exponential decay. This means exponential decay models would systematically underestimate long-term antibody persistence for all targets.

**[NOTE — Unit mapping]** Aiemjoy Table S1 reports the decay rate as "α" in units that appear to be ELISA units per year, while Teunis uses "ν" in units of 1/days. The Aiemjoy notation appears to use α and r as direct model parameters (where Teunis uses ν and r), with the correspondence: Aiemjoy's α ≈ Teunis's ν after unit conversion. The exact parameterization used in the JAGS implementation should be verified from the analysis code at https://github.com/kaiemjoy/TyphoidSeroIncidence.git before re-implementing.

### Fitting Method
- Bayesian hierarchical framework using MCMC (JAGS 4.3.0 via rjags)
- Fitted in three age strata: <5 years, 5–15 years, >15 years
- Results reported as median and 95% credible intervals (CrI) from posterior distributions
- Suspected reinfections (≥3-fold increase in ≥2 antigen-isotypes at ≥3 months post-onset) excluded from decay estimation

### Reinfection Handling (p.e580, Supplement p.2)
- 37 participants met reinfection definition
- Reinfection incidence per 100 person-years: Bangladesh 5.2 (95% CI 1.8–8.5), Pakistan 4.9 (4.0–6.8), Nepal 0.7 (0.0–1.6), Ghana 5.6 (0.0–13.3)
- Median time to reinfection detection: 13.8 months (IQR 10.3–18.4)
- Observations at and after suspected reinfection excluded from longitudinal decay models

## Fitted Model Parameters (Table S1, Supplement p.15)

### Overall (all countries pooled) — Median (95% CrI)

| Parameter | HlyE IgA | HlyE IgG | LPS IgA | LPS IgG | Vi IgG |
|-----------|----------|----------|---------|---------|--------|
| **Peak response y1** (EU) | 94.0 (14.1–892.2) | 226.3 (56.5–1096.1) | 271.1 (31.3–3949.3) | 200.5 (47.9–1550) | 542.7 (143.5–3015.6) |
| **Decay rate α** (EU/yr) | 0.03 (0–1.8) | 0.02 (0–0.6) | 0.03 (0–3.6) | 0.07 (0–1.3) | 0 (0–0.4) |
| **Shape factor r** | 2.81 (1.7–5.7) | 2.18 (1.5–3.9) | 2.89 (1.5–7.9) | 2.16 (1.4–4.1) | 1.25 (1–3.1) |
| **Time to peak t1** (days) | 19.63 (5.2–70.8) | 15.6 (3.8–63.7) | 11.07 (2.5–36.9) | 5.48 (0.7–37.1) | 2.92 (0.6–11.3) |
| **Baseline y0** (EU) | 13.81 (0.8–254.5) | 14.47 (0.6–323.9) | 32.67 (1.4–655.9) | 13.31 (0.5–302.1) | 126.72 (7.6–1806.3) |

### Age-Stratified Parameters

#### Peak antibody response (y1, ELISA units) — Median (95% CrI)

| Age | HlyE IgA | HlyE IgG | LPS IgA | LPS IgG | Vi IgG |
|-----|----------|----------|---------|---------|--------|
| <5 | 37.6 (7–469.8) | 159.7 (20.7–1686.8) | 81.7 (25.8–254.6) | 157.1 (23.3–1601.8) | 325.4 (7.4–22290.5) |
| 5–15 | 103.3 (13–1385.9) | 222.7 (49.6–1541.7) | 338.6 (56.7–8402.8) | 220.0 (48.9–4235.2) | 399.6 (147.4–6077.6) |
| 16+ | 130.9 (10.1–2569.6) | 248.6 (30–2726.3) | 88.7 (4–2788.1) | 112.8 (8.9–1625.4) | 650.9 (218.5–3107.3) |

**[NOTE]** Peak HlyE responses increase with age. LPS and Vi peaks are more variable.

#### Decay rate (α, ELISA units per year) — Median (95% CrI)

| Age | HlyE IgA | HlyE IgG | LPS IgA | LPS IgG | Vi IgG |
|-----|----------|----------|---------|---------|--------|
| <5 | 0.03 (0–1.5) | 0.04 (0–0.6) | 0.05 (0–4.4) | 0.11 (0–1.5) | 0 (0–0) |
| 5–15 | 0.03 (0–2.3) | 0.03 (0–0.8) | 0.06 (0–4.3) | 0.08 (0–1.6) | 0.02 (0–3.9) |
| 16+ | 0.02 (0–0.5) | 0.04 (0–0.3) | 0.01 (0–1.9) | 0.05 (0–0.5) | 0.01 (0–0.1) |

**[CRITICAL — Vi IgG]** Anti-Vi IgG decay rate is estimated as **0 (0–0) for age <5** and **0 (0–0.4) overall** — the model finds essentially **no decay** in Vi IgG after infection. This is confirmed in main text (p.e583): "anti-Vi antibody responses plateaued with no subsequent decay (rate=0)." Median anti-Vi IgG levels remained elevated above baseline for **at least 32 months** (Supplement p.13, Fig S11).

#### Shape factor (r) — Median (95% CrI)

| Age | HlyE IgA | HlyE IgG | LPS IgA | LPS IgG | Vi IgG |
|-----|----------|----------|---------|---------|--------|
| <5 | 3.24 (1.8–7) | 2.2 (1.6–3.2) | 3.26 (1.6–10.3) | 2.19 (1.3–5.1) | 2.5 (1.2–4.7) |
| 5–15 | 2.58 (1.7–4.5) | 2.08 (1.4–4) | 2.48 (1.3–8.6) | 2.07 (1.4–3.6) | 1.28 (1–2.8) |
| 16+ | 2.64 (1.9–4.1) | 2.01 (1.4–3.4) | 3.26 (2.4–4.9) | 2.26 (1.8–3) | 1.35 (1.2–1.8) |

**[NOTE]** All shape factors r > 1, indicating power-law decay that is faster-than-exponential initially but slower-than-exponential at long times. Vi IgG has the lowest shape factors (closest to 1), consistent with the near-zero decay rate — the shape parameter is less informative when there is essentially no decay.

#### Time to peak (t1, days since fever onset) — Median (95% CrI)

| Age | HlyE IgA | HlyE IgG | LPS IgA | LPS IgG | Vi IgG |
|-----|----------|----------|---------|---------|--------|
| <5 | 28.53 (8–85.1) | 16.15 (2.3–104.5) | 9.3 (2.3–33.7) | 6.2 (0.6–56.3) | 1.4 (0.6–3.3) |
| 5–15 | 19.95 (6.4–54) | 12.34 (3.5–41.4) | 6.72 (1.2–32.4) | 6.19 (1.1–32.9) | 4.71 (2.1–10.4) |
| 16+ | 20.76 (6–66.4) | 19.06 (6.4–46.8) | 26.52 (13.7–52.1) | 10.74 (1.9–43.2) | 9.57 (2.8–25.7) |

**[NOTE]** Antibody responses reach peak within ~3 weeks of fever onset (p.e581). HlyE peaks later than LPS. Vi IgG peaks very rapidly (median 1.4–9.6 days) but note this may reflect that Vi responses are modest and plateau early.

#### Baseline antibody response (y0, ELISA units) — Median (95% CrI)

| Age | HlyE IgA | HlyE IgG | LPS IgA | LPS IgG | Vi IgG |
|-----|----------|----------|---------|---------|--------|
| <5 | 20.77 (2.4–181.4) | 14.67 (1.2–181.2) | 77.47 (24–216.8) | 10.6 (0.5–228.8) | 0.6 (0.2–1.7) |
| 5–15 | 14.93 (0.7–295.4) | 15.06 (0.7–313.2) | 19.67 (1.1–298.6) | 20.78 (1.2–344.1) | 164.55 (64.2–354.4) |
| 16+ | 13.99 (0.7–290.1) | 15.59 (0.7–381.4) | 68.8 (3.2–1441.8) | 32.72 (1.8–539.2) | 201.78 (22.3–1762.1) |

**[CRITICAL — Vi IgG baseline]** Vi IgG baseline (y0) is very low in <5 year olds (0.6 EU) but very high in 5–15 (164.6 EU) and 16+ (201.8 EU) groups. This reflects cumulative prior exposure in endemic settings. The ratio of 28-day case response to population mean is only 1.5 for Vi IgG versus 30.9 for HlyE IgG (p.e583), making Vi a poor discriminator of acute infection from prior exposure.

## Duration of Elevated Responses (p.e582)

Median time above baseline, from model fits:
| Antigen-Isotype | Months above baseline |
|----------------|-----------------------|
| HlyE IgG | 29 |
| LPS IgG | 14 |
| HlyE IgA | 11 |
| LPS IgA | 3 |
| Vi IgG | ≥32 (no observed decay) |

**[CRITICAL]** This ordering — IgG longer than IgA for both antigens, HlyE longer than LPS — is consistent across the study. Vi IgG is the most durable, with no measurable decay over the ~2-year observation window.

## Key Findings for Waning/Boosting Model

### Consistency Across Settings
- Model-fitted antibody trajectories were **similar across all four countries** (Fig 2A, p.e583)
- Peak antibody response distributions similar across countries; differences between country distributions centered near 0
- **Some inter-country variation in decay rate**: Bangladesh had slowest decay, Ghana fastest
  - Proposed explanations: (1) undetected reinfections more common in high-transmission settings; (2) secondary responses may wane more slowly in frequently-exposed populations

### Severity Independence
- Antibody kinetics did **not differ between hospitalized and non-hospitalized** patients (Fig S6, Supplement p.8)
- This suggests antibody response magnitude and decay are not dependent on clinical disease severity

**[OPEN]** Whether truly asymptomatic infections produce similar antibody kinetics is unknown — the study could only compare hospitalized vs. non-hospitalized clinical cases, not asymptomatic infections.

### S. Typhi vs. S. Paratyphi A (Figure 1C, p.e582)
- HlyE IgG and IgA trajectories were **similar** between S. Typhi and S. Paratyphi A
- **LPS IgA and IgG peak responses were lower for S. Paratyphi A** — expected because S. Paratyphi A only shares the O12 antigen with S. Typhi (S. Typhi has O9 and O12)

### Anti-Vi IgG Sub-analysis (p.e583; Fig S11, Supplement p.13)
- Conducted in Nepal sub-cohort only
- After a modest initial rise, anti-Vi IgG **plateaued with no subsequent decay** (decay rate α = 0)
- Median anti-Vi IgG remained elevated above baseline for **≥32 months**
- Ratio of 28-day case response to population mean:
  - Anti-Vi IgG: **1.5** (IQR 1.0–2.1) — lowest among all antigen-isotypes
  - Anti-HlyE IgG: **30.9** (IQR 16.1–54.6) — highest
- Anti-Vi IgG responses **did not increase with age** in cross-sectional serosurveys, echoing Pulickal et al. 2009 from Kathmandu (ref 26)
- **Vi cannot distinguish natural infection from Vi-based vaccination**

**[CRITICAL — for waning modeling]** The near-zero decay rate for anti-Vi IgG after natural infection is a key finding. However, the very modest boosting (only 1.5× above population mean at 28 days) means that Vi IgG levels in endemic populations are predominantly determined by cumulative lifetime exposure, not recent infection. This is qualitatively different from the sharp boost-then-decay seen with HlyE. Any waning model for anti-Vi must account for this plateau behavior.

## Measurement Parameters (Table S3, Supplement p.17)

| Antigen-Isotype | Biologic Noise* | Measurement Error (CV) | Lower Censoring Limit |
|-----------------|-----------------|------------------------|-----------------------|
| **HlyE IgA** | | | |
| Bangladesh | 2.87 | 0.28 | 0.38 |
| Ghana | 2.87 | 0.24 | 0.18 |
| Nepal | 2.87 | 0.24 | 0.85 |
| Pakistan | 2.87 | 0.28 | 0.51 |
| **HlyE IgG** | | | |
| Bangladesh | 1.08 | 0.31 | 0.79 |
| Ghana | 1.08 | 0.16 | 0.65 |
| Nepal | 1.08 | 0.13 | 1.89 |
| Pakistan | 1.08 | 0.15 | 1.59 |
| **LPS IgA** | | | |
| Bangladesh | 2.72 | 0.30 | 0.66 |
| Ghana | 2.72 | 0.16 | 0.86 |
| Nepal | 2.72 | 0.11 | 1.79 |
| Pakistan | 2.72 | 0.25 | 5.13 |
| **LPS IgG** | | | |
| Bangladesh | 0.96 | 0.30 | 0.99 |
| Ghana | 0.96 | 0.20 | 0.88 |
| Nepal | 0.96 | 0.18 | 0.65 |
| Pakistan | 0.96 | 0.27 | 4.84 |

\* Biologic noise = upper 95th percentile from North American controls (ELISA units, natural scale)

## Data Quality Notes

### Strengths
- Very large longitudinal cohort (n=1420) across four countries with diverse epidemiological settings
- Up to 2 years follow-up with median ~1 year
- Blood culture-confirmed cases — gold standard diagnosis
- Bayesian hierarchical modeling accounts for inter-individual heterogeneity
- DBS-plasma correlation validated (Fig S2: r = 0.84–0.98 across antigen-isotypes)
- Reinfections explicitly identified and handled
- Analysis code publicly available (GitHub)

### Limitations
- **Vi IgG measured only in Nepal (and Ghana)**; not available for Bangladesh or Pakistan sub-cohorts
- Vi IgG boosting is very modest (1.5× population mean), limiting power to characterize Vi waning precisely
- Power-function decay model chosen a priori; other functional forms not compared
- Waning parameters reflect primary infection kinetics — secondary/tertiary infections may differ (ref 35: Diekmann et al. 2018)
- Cannot determine antibody kinetics of truly asymptomatic infections
- Age strata are broad (<5, 5–15, 16+); continuous age-dependence not formally modeled
- This is a corrected publication (correction published Feb 6, 2023)

### Key Page References
- Study design: p.e579–e580
- Table 1 (case demographics): p.e581
- Figure 1 (longitudinal kinetics, age, serovar): p.e582
- Figure 2 (country comparisons): p.e583
- Anti-Vi sub-analysis: p.e583
- Table S1 (all fitted parameters): Supplement p.15
- Figure S11 (Vi IgG kinetics): Supplement p.13
- Table S3 (measurement noise): Supplement p.17

## Cross-References
- **Teunis PFM, van Eijkeren JCH, de Graaf WF, Marinović AB, Kretzschmar MEE. Linking the seroresponse to infection to within-host heterogeneity in antibody production. Epidemics 2016; 16:33–39.** (ref 21) — Derives the power-function decay model used here. Key equations: ODE y'(t) = −νy(t)^r (Eq. 13); solution y(t>t1) = y1·(1+(r−1)·y1^(r−1)·ν·(t−t1))^(−1/(r−1)) (Eq. 14). Shows power-function decay arises from Gamma-distributed heterogeneity in antibody production site decay rates. Fitted to pertussis IgG-PT with r ≈ 2.2 (95% PI 1.7–2.8). Predicts 80% of subjects still above 5 IU/mL at 30 years post-infection under power-function vs. ~0% under exponential.
- de Graaf WF, Kretzschmar MEE, Teunis PFM, Diekmann O. A two-phase within-host model for immune response and its application to serological profiles of pertussis. Epidemics 2014; 9:1–7. (ref 20) — original two-phase within-host model (exponential decay version) that Teunis 2016 extends
- Diekmann et al. 2018, J Math Biol (ref 35) — theoretical framework for waning and boosting dynamics
- Pulickal et al. 2009, Clin Vaccine Immunol (ref 26) — earlier study of natural humoral immune response to S. Typhi in Kathmandu; found no increase in Vi IgG with age
- Watson et al. 2017, PLoS NTD (ref 8) — cross-sectional Vi seroepidemiological survey in Fiji
- Garrett et al. 2022, Lancet Glob Health (ref 17) — SEAP clinical incidence estimates from same catchment areas
- Meiring et al. 2021, Lancet Glob Health (ref 27) — used Vi IgG seroconversion to estimate incidence; similar order of magnitude results
- Analysis code: https://github.com/kaiemjoy/TyphoidSeroIncidence.git

## Fit Role Assessment
**[ASSISTANT-PROPOSED]** **CORE** for waning characterization — provides the most comprehensive published dataset on post-infection antibody waning kinetics for typhoidal Salmonella across multiple antigen-isotypes. The anti-Vi IgG result (essentially no decay over ≥32 months) is particularly important for modeling natural immunity dynamics, even though Vi is poorly represented relative to HlyE and LPS. The power-function decay model with fitted parameters provides a quantitative framework that can be directly incorporated into immunity models.

## Key Extractions for Waning/Boosting Model

### Summary of waning behavior by antigen-isotype (overall, all ages)

| Antigen-Isotype | Peak y1 (EU) | Decay rate α | Shape r | Months above baseline | Boosting ratio (28d case / pop mean) |
|-----------------|-------------|-------------|---------|----------------------|--------------------------------------|
| HlyE IgG | 226.3 (56.5–1096.1) | 0.02 (0–0.6) | 2.18 (1.5–3.9) | 29 | 30.9 (16.1–54.6) |
| HlyE IgA | 94.0 (14.1–892.2) | 0.03 (0–1.8) | 2.81 (1.7–5.7) | 11 | — |
| LPS IgG | 200.5 (47.9–1550) | 0.07 (0–1.3) | 2.16 (1.4–4.1) | 14 | — |
| LPS IgA | 271.1 (31.3–3949.3) | 0.03 (0–3.6) | 2.89 (1.5–7.9) | 3 | — |
| **Vi IgG** | **542.7 (143.5–3015.6)** | **0 (0–0.4)** | **1.25 (1–3.1)** | **≥32** | **1.5 (1.0–2.1)** |

### Anti-Vi IgG — Key quantitative findings
- **Decay rate: 0 (95% CrI 0–0.4)** — no measurable waning over observation period
- **Duration above baseline: ≥32 months** (study observation limit)
- **Peak: 542.7 EU (143.5–3015.6)** — but high baseline in older ages (164.6 in 5–15, 201.8 in 16+)
- **Boosting ratio at 28 days: 1.5× (IQR 1.0–2.1)** — modest boost relative to endemic population background
- **No increase in Vi IgG with age** in cross-sectional serosurveys
- Measured only in Nepal sub-cohort (n=543 cases, though Nepal had very few <5 year olds: n=12)

**[CRITICAL — Effective censoring of Vi waning signal]** The near-zero Vi IgG decay rate almost certainly reflects an inability to detect waning rather than true absence of waning. The argument:

1. **Signal-to-noise ratio is fundamentally limiting.** The boosting ratio for Vi IgG at 28 days is only **1.5× (IQR 1.0–2.1)** the population mean. Compare HlyE IgG at **30.9× (16.1–54.6)**. The Vi "signal" from a single infection is barely distinguishable from background.

2. **Baseline dominates the trajectory.** Vi IgG baseline y0 in the 5–15 age group is **164.6 EU** and the peak y1 is **399.6 EU** — the peak is only ~2.4× baseline. In adults 16+, y0 = **201.8 EU** vs. y1 = **650.9 EU** (~3.2× baseline). By contrast, HlyE IgG has y0 ~15 EU vs. y1 ~225 EU (~15× baseline). When the peak-to-baseline ratio is small, any waning quickly brings the trajectory into the noise floor of the baseline, making it undetectable.

3. **The model structure itself creates a censoring effect.** In the Teunis power-function model, y(t) = y1·(1 + (r−1)·y1^(r−1)·ν·(t−t1))^(−1/(r−1)), the decay is relative to y1, but the observed trajectory also includes y0 as a floor. If y1/y0 is small, the model cannot distinguish "slow decay toward y0" from "no decay" because the trajectory never moves far enough from y0 to constrain ν. The MCMC will pile posterior mass at ν ≈ 0 simply because values of ν that produce observable decay over the ~2-year follow-up are inconsistent with trajectories that hover near y0. This is effective **right-censoring of the decay rate parameter**.

4. **The <5 age group is the most informative but has the least data.** For <5 year olds, Vi IgG baseline is very low (y0 = 0.6 EU) with peak 325.4 EU — a huge ratio (~540×) that in principle would allow clean waning detection. But Nepal contributed only **12 cases** aged <5 (Table 1), and Vi was measured only in Nepal. The model correctly returns α = 0 (0–0) for this stratum — not because there's no waning, but because there's essentially no data.

5. **Cross-sectional Vi IgG shows no age trend** — the paper notes this (p.e583) and Pulickal et al. 2009 found the same in Kathmandu. If Vi truly didn't wane, one would expect population Vi IgG to accumulate with age (more infections → higher levels). The flat age profile is more consistent with a steady-state where boosting and waning roughly balance, or with Vi responses being driven primarily by something other than recent acute infection (e.g., chronic low-level exposure, carriage contacts, cross-reactive exposures).

6. **The HlyE age-trend contrast confirms the test works.** This is the model-free clincher: in the same cross-sectional serosurveys, anti-HlyE IgG and IgA **do** show significant increases with age at all sites (Fig S8, Supplement p.10; Kruskal-Wallis p < 0.001 at nearly every site for both HlyE isotypes). This is exactly what one expects for an antigen with slow waning (HlyE IgG duration above baseline ~29 months) in an endemic setting — each successive infection adds an increment that hasn't fully decayed before the next. LPS IgG and IgA also show significant age trends. **Vi IgG is the only antigen-isotype that fails to show an age trend despite being measured in the same populations.** If the explanation were "Vi just doesn't boost much from natural infection," that would also mean the longitudinal Vi kinetics are uninformative about waning (you can't characterize the decay of a signal that barely exists). If the explanation were "Vi boosts but truly doesn't wane," then Vi should accumulate with age at least as strongly as HlyE — it doesn't. Either way, the "decay rate = 0" parameter estimate cannot be taken at face value as evidence of durable Vi immunity from natural infection.

**Bottom line**: The Aiemjoy Vi IgG "no waning" result should be interpreted as "waning is undetectable with this study design and signal-to-noise ratio," not as evidence that Vi IgG truly does not wane after natural infection. The combination of modest boosting, high endemic baseline, Vi measurement restricted to one country, and very few young children makes this parameter effectively unidentifiable. Independent data (e.g., from challenge studies, vaccine trials with controlled Vi exposure, or low-endemic settings where baseline is low) would be needed to characterize Vi IgG waning.

**[NOTE — Teunis model verification]** The exact functional forms from Teunis et al. 2016 are now documented above (Section "Antibody Kinetics Model"). The unit mapping between Aiemjoy's Table S1 parameterization (α in EU/year, r) and Teunis's original notation (ν in 1/days, r) should be verified from the analysis code before re-implementing.
