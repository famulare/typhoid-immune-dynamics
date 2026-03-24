# Mylona et al. 2024 - The Identification of Enteric Fever-Specific Antigens for Population-Based Serosurveillance

## Citation
Mylona E, Hefele L, Nga Tran Vu Thieu, Tan Trinh Van, Chau Nguyen Ngoc Minh, Anh Tran Tuan, Karkey A, Dongol S, Basnyat B, Phat Voong Vinh, Thanh Ho Ngoc Dan, Russell P, Charles RC, Parry CM, Baker S. The Identification of Enteric Fever-Specific Antigens for Population-Based Serosurveillance. J Infect Dis 2024; 229:833–44. https://doi.org/10.1093/infdis/jiad242

Clinical Trial Registration: ISRCTN63006567

## Study Overview
- **Primary question**: Can protein antigens expressed during S. Typhi/S. Paratyphi A infection serve as better markers of enteric fever exposure than conventional Vi, O-antigen, and H-antigen serological targets?
- **Study type**: Longitudinal observational study (descriptive) nested within a treatment trial; comparative serology across patient groups and community controls
- **Institutions**: Cambridge Institute of Therapeutic Immunology and Infectious Diseases (Cambridge, UK); Hospital for Tropical Diseases / Oxford University Clinical Research Unit (Ho Chi Minh City, Vietnam); Oxford University Clinical Research Unit (Kathmandu, Nepal); Harvard Medical School; Clinical Sciences (Liverpool); IAVI Human Immunology Laboratory (London)
- **Population**: Blood culture-confirmed S. Typhi and S. Paratyphi A patients, culture-negative febrile patients, and afebrile community controls in Kathmandu, Nepal
- **Country**: Nepal (Patan Hospital, Kathmandu)
- **Enrollment period**: From parent trial Koirala et al. 2013 (ISRCTN63006567); treatment trial comparing gatifloxacin vs. ofloxacin for uncomplicated enteric fever
- **Funding**: Wellcome Trust (grant 215515/Z/19/Z to S.B.); NIH/NIAID (R01 AI134814 to R.C.C.)
- **Follow-up duration**: 3 months (4 time points)

**[CRITICAL]** This paper does NOT contain fitted waning model parameters. It is a descriptive comparative study using LOESS smoothing, Friedman test, Wilcoxon signed-rank test, Kruskal-Wallis test, and Spearman correlations. Its value for the waning metastudy is qualitative: it characterizes the shape and discriminatory power of anti-Vi IgG responses relative to protein antigens over the first 3 months post-infection, providing independent evidence that Vi IgG is a poor marker of acute enteric fever exposure.

## Subject Characteristics

### Patient Groups (Table 1, p.834)

| Characteristic | Community Controls (CC) | All Patient Groups | FCN | SPA | ST |
|---|---|---|---|---|---|
| **N** | 49 | 486 | 322 | 61 | 103 |
| Female, n (%) | NA | 166 (34.2) | 111 (34.5) | 23 (37.7) | 32 (31.1) |
| Male, n (%) | NA | 320 (65.9) | 211 (65.5) | 38 (62.3) | 71 (68.9) |
| Age, mean (SD) | NA | 17.3 (9.8) | 17.9 (10.6) | 16.8 (8.1) | 15.6 (7.9) |
| Age NA | ... | 3 | 1 | 0 | 2 |

Abbreviations: CC = afebrile community controls; FCN = febrile culture-negative patients; SPA = S. Paratyphi A-confirmed; ST = S. Typhi-confirmed.

**[NOTE]** CC group had only a single time-point sample (baseline). Sex data not available for CC. All patients were from Kathmandu, a known endemic area for enteric fever.

**[NOTE]** Only patients presenting with 72-96 hours of fever (day 1) were included. This is a symptomatic, blood culture-confirmed cohort -- no asymptomatic or subclinical infections.

### Parent Trial Context
- Open-label, randomized, controlled trial comparing gatifloxacin vs. ofloxacin for uncomplicated enteric fever (Koirala et al. 2013, ref 17)
- Patients randomly assigned to 7-day course of treatment
- Longitudinal samples collected from 103 S. Typhi and 61 S. Paratyphi A patients
- Comparative samples from 322 culture-negative febrile patients (FCN)
- 49 afebrile community controls (CC) -- single time-point only

## Sampling Schedule (p.834)

| Time Point | Label | Day |
|---|---|---|
| Day 1 (enrollment, day of fever) | D1 | 1 |
| Day 8 | D8 | 8 |
| Month 1 | M1 | 28 |
| Month 3 | M3 | 90 |

- Blood collected in EDTA tubes
- Immediately separated into cells and plasma
- CC group: single time-point plasma sample only

**[NOTE]** The 3-month follow-up window is short relative to waning timescales. This study captures only the early acute/convalescent phase, not the longer-term decay observable in Aiemjoy 2022 (up to ~2 years).

## Antigens Tested

### Panel (Supplementary Table 3, referenced on p.834-835)
A panel of 17 S. Typhi and S. Paratyphi A antigens was purified and tested by indirect ELISA. The paper focuses on 6 key antigens:

| Antigen | Gene/Protein | Description | Discriminatory? |
|---|---|---|---|
| **STY1479 (YncE)** | ATP-binding protein | Putative function; iron restriction-induced; potential chronic carrier marker | Yes -- elevated in ST and SPA vs. FCN |
| **STY1886 (CdtB)** | Typhoid toxin component | Damages host DNA, eliminates immune cells | Yes -- elevated in ST and SPA vs. FCN |
| **STY1498 (HlyE)** | Hemolysin E | Pore-forming toxin; absent from most NTS | Yes -- most strongly elevated in ST and SPA |
| **Vi** | Vi capsular polysaccharide | S. Typhi capsule; ELISA at 15 ug/mL coating | **No** -- poor discrimination |
| **O2** | O-antigen conformation | Serovar-specific O-antigen | Partially -- elevated but cross-reactive |
| **O9** | O-antigen conformation | Serovar-specific O-antigen | Partially -- elevated but cross-reactive |

Additional antigens tested but showing limited or no discrimination: STY1612, STY4539, STY1086, STY1703, STY1767, and others (Supplementary Figures 3-5).

### ELISA Methods (p.834)
- Indirect ELISA on flat-bottom 96-well plates (Nunc 2404, Thermo Scientific)
- Antigen coating: 1 ug/mL for proteins and Vi; 15 ug/mL for O2 and O9
- Coating buffer: 50 mM carbonate bicarbonate buffer
- Plasma dilution: 1:200
- Detection: alkaline phosphatase-conjugated anti-human IgG, 1 hour at ambient temperature
- Substrate: phosphate-buffered saline p-nitrophenyl phosphate (SigmaFAST N1891), 30 minutes
- Absorbance read at 405 and 490 nm
- Results expressed as ELISA units (EU) relative to standard curves (4-parameter modified Hill plot)
- Reference standard: pooled antigen-specific antibody human plasma from Nepalese plasma bank
- Positive cutoff: absorbance > blank wells mean + 2 SD
- Standard curves: 10-point 2-fold serial dilutions; 1 EU = reciprocal of dilution giving OD = 1
- All responses log-transformed for analysis

## Statistical Methods (p.835)

- **LOESS smooth functions** fitted to longitudinal antibody trajectories (Figure 1)
- **Friedman test** followed by **Wilcoxon signed-rank test with Bonferroni correction** for within-group time-point comparisons (Figures 2, 3)
- **Kruskal-Wallis test** followed by **Wilcoxon signed-rank test with Bonferroni correction** for between-group comparisons at each time point (Figures 5, 6)
- **Spearman correlation coefficients** for inter-antigen correlation at day 8 (Figure 4)
- Nonparametric methods chosen after Shapiro-Wilk (normality) and Levene (homogeneity of variances) tests
- All analysis in R (R Foundation for Statistical Computing, 2014; ref 18)
- Significance threshold: P < .05

**[CRITICAL]** No parametric waning models were fitted. No decay rate parameters were estimated. The study is purely descriptive. The LOESS smooths provide qualitative trajectory shapes but no quantitative waning parameters.

## Results

### Longitudinal Antibody Trajectories (Figure 1, p.836; text p.835)

#### Anti-Vi IgG (Figure 1D)
- **FCN group**: LOESS smooth essentially flat over 3 months; comparable to community control mean (dashed line)
- **SPA group**: antibody dynamics against Vi "were not substantially different" from FCN group (p.835)
- **ST group**: Vi IgG trajectory showed "an early peak" returning "rapidly to baseline" by month 3 (p.835)
- Anti-Vi IgG responses in the ST group were described as "stable" and "comparable to those measured at day 1" (p.835-836, in context of Figure 3A discussion)

**[CRITICAL]** The LOESS smooth for anti-Vi IgG in S. Typhi patients (Figure 1D, rightmost panel) shows a very modest rise from ~day 1 to ~day 25, with the LOESS curve remaining within approximately 2.0-2.5 log EU throughout the entire 90-day observation window. The curve barely exceeds the community control baseline (dashed line at ~2.3 log EU). This contrasts dramatically with O2, O9, and HlyE, which show clear rises of 1+ log units above baseline in ST patients.

#### Protein Antigens That DID Show Clear Dynamics
- **O2 (Figure 1A)**: Markedly elevated LOESS curve in SPA and ST groups vs. FCN; clear peak around day 25-50
- **O9 (Figure 1B)**: Similar to O2; markedly elevated in SPA and ST groups
- **STY1498/HlyE (Figure 1C)**: Most distinct fitted curves; strongly elevated in ST and SPA vs. FCN; remained elevated above baseline over 3 months
- **STY1479/YncE (Figure 1E)**: Higher in SPA and ST vs. FCN; remained above baseline over 3-month period
- **STY1886/CdtB (Figure 1F)**: Higher in SPA and ST vs. FCN; remained above baseline

### Within-Group Time-Point Comparisons (Figures 2-3, pp.837-838)

#### Anti-Vi IgG over time (Figure 3A)
- **FCN**: No significant differences between any time points
- **SPA**: Significant increase from D1 to D8 (shown as ** in Figure 3A SPA panel) and D1 to M1 (shown as *)
- **ST**: **No significant pairwise differences** between any time points for Vi IgG

**[CRITICAL]** In S. Typhi patients -- the primary group where Vi exposure is expected -- anti-Vi IgG showed NO statistically significant change at any time point vs. day 1 (Figure 3A, ST panel). This is in stark contrast to protein antigens (HlyE, CdtB, YncE) and O-antigens (O2, O9) which showed highly significant increases (P < .0001 at multiple comparisons).

#### Protein antigens over time
- **STY1479/YncE** (Figure 2A): Significant increases in SPA (* at D8, from Friedman) and ST (* at M1) groups
- **STY1886/CdtB** (Figure 2B): Significant increases in SPA (*** D1-D8, *** D1-M1) and ST (* D1-M1)
- **STY1498/HlyE** (Figure 2C): Significant increases across all groups; most prominent -- **** at D1-D8 and multiple comparisons in FCN, SPA, and ST groups

#### O-antigens over time (Figures 3B-C)
- **O2** (Figure 3B): Highly significant increases in SPA (**** D1-D8, **** D1-M1, **** D1-M3) and ST (** D1-D8, **** D1-M1, *** D1-M3)
- **O9** (Figure 3C): Highly significant increases in all three patient groups (FCN, SPA, ST) -- **** at most comparisons

### Peak Antibody Responses (p.835-836)

The paper describes peak timing qualitatively:
- **O2, O9, HlyE**: Increased significantly (P < .05) from D1 to D8 in ST and/or SPA patients; remained elevated at M1
- **STY1479 (YncE), STY1886 (CdtB)**: Similar rise pattern but "less prominent"
- **Vi**: "The anti-Vi IgG responses remained stable with time in the S. Typhi group and were comparable to those measured at day 1" (p.835-836)
- **Most IgG trajectories exhibited a notable decline by month 3** for protein antigens (p.836): STY1479, STY1498 (HlyE), O2, O9 all showed "a comparable decrease in both the S. Typhi and S. Paratyphi A groups"

**[NOTE]** The decline by M3 for protein antigens indicates that even within the 3-month window, waning is observable for HlyE, O2, O9, etc. The absence of any such decline for Vi IgG is consistent with either (a) no boosting occurred (so there is nothing to wane) or (b) the signal is too small to detect waning above noise.

### Correlation Between Antibody Responses (Figure 4, pp.837-839)

#### Day 8 Spearman Correlations

**Community Controls (Figure 4A):**
All antibody responses correlated moderately to strongly (rho > 0.6) with each other. Vi correlated with protein antigens at rho = 0.480-0.668 and with O2 at 0.571 and O9 at 0.456.

**FCN group (Figure 4B):**
Moderate correlations (rho > 0.5) among most antigens. Vi correlated moderately with protein antigens (rho = 0.174-0.790). Anti-O2 and O9 responses correlated poorly with most protein antigens.

**S. Paratyphi A group (Figure 4C):**
Vi correlated with STY1479 at rho = 0.825, STY1886 at 0.705, STY1498 at 0.091. O2 correlated poorly with Vi (rho = -0.078). STY1498 (HlyE) correlated poorly with Vi (rho = 0.091).

**[CRITICAL] S. Typhi group (Figure 4D):**
| Antigen Pair | Spearman rho |
|---|---|
| Vi vs. STY1479 (YncE) | 0.111 |
| Vi vs. STY1886 (CdtB) | 0.088 |
| Vi vs. STY1498 (HlyE) | -0.049 |
| Vi vs. O2 | 0.094 |
| Vi vs. O9 | 0.067 |

**All Vi correlations with other antigens in the S. Typhi group are rho < 0.2, and Vi-HlyE is actually slightly negative (-0.049).** This means that among confirmed S. Typhi patients at day 8 (peak of acute response), the anti-Vi IgG response provides essentially no information about the magnitude of the protein antigen or O-antigen response, and vice versa. Vi IgG is behaving independently of the acute infection response.

The text confirms (p.838): "STY1498 (HlyE)...correlated poorly to all antibody responses (rho < 0.2) due to its unique trajectory" in ST patients. And specifically for Vi in ST patients: correlation with protein antigens was poor (rho > 0.4 threshold not met for most).

**[NOTE]** In SPA patients, Vi had moderate correlations with STY1479 (0.825) and STY1886 (0.705), but these patients would not be expected to have strong Vi responses (Vi capsule is specific to S. Typhi). The SPA correlations may reflect shared endemic background exposure rather than acute infection.

### Between-Group Comparisons: Anti-Vi IgG vs. Controls (Figures 5-6, pp.840-841)

#### Anti-Vi IgG: Patient groups vs. CC and FCN (Figure 6A, p.841)

| Time Point | CC vs. ST | CC vs. SPA | FCN vs. ST | FCN vs. SPA | Notes |
|---|---|---|---|---|---|
| Day 1 | NS | NS | NS | NS | All groups comparable at baseline |
| Day 8 | * (SPA-CC) or ** (SPA-CC) | * or ** | * | NS | Only SPA-CC and ST-FCN borderline significant |
| Day 28 (M1) | * (SPA-CC) | * | NS | NS | Minimal significance |
| Day 90 (M3) | NS | NS | NS | NS | **All groups comparable at month 3** |

**[CRITICAL]** From Figure 6A: At day 8, there are modest significant differences between SPA and CC (* or **) and between ST and FCN (*). By day 28, only SPA-CC retains marginal significance. By day 90 (month 3), **there are NO significant differences between any patient group and controls for anti-Vi IgG**. This means that by 3 months after blood culture-confirmed S. Typhi infection, anti-Vi IgG levels have returned to population background levels.

**[NOTE]** Significance annotations in Figure 6A are read from the figure panels. Exact P-values are not reported in the text for Vi between-group comparisons. The significance symbols follow the legend: *P < .05, **P < .01, ***P < .001, ****P < .0001.

#### Contrast with protein antigens and O-antigens (Figures 5-6)
- **STY1498/HlyE** (Figure 5C): Strongly significant (**** P < .0001) separation between ST/SPA and CC/FCN at ALL time points (D1, D8, M1, M3). HlyE discriminates enteric fever patients from controls even at 3 months.
- **O2** (Figure 6B): Strongly significant (**** at most comparisons) separation between ST/SPA and CC/FCN at D1, D8, M1, and M3.
- **O9** (Figure 6C): Similarly strong separation at all time points.
- **STY1479/YncE** (Figure 5A): Significant at D8 and M1 but not at D1 or M3 for some comparisons.
- **STY1886/CdtB** (Figure 5B): Significant at D8 and M1 for SPA and ST vs. FCN, but weaker signal.

### Authors' Explicit Assessment of Vi (Discussion, p.842)

Direct quotes:

> "the baseline antibody responses against Vi in enteric fever patients suggest that anti-Vi antibodies are circulating in this community and the Vi is not suitable as a marker of enteric fever exposure due to high endemic background"

> "The use of Vi in this context is likely to be further hampered by ongoing TCV vaccination schemes and it is known that responses to Vi cannot distinguish between natural infection and vaccine-mediated exposure"

**[CRITICAL]** The authors themselves conclude Vi is unsuitable for serosurveillance of enteric fever. Their reasoning aligns with our analysis in the Aiemjoy extraction: high endemic background makes Vi uninformative as a marker of recent acute infection.

### Antigens Recommended for Serosurveillance (Discussion, p.842)

The paper identifies a subset of protein antigens as promising serosurveillance candidates:
1. **STY1498 (HlyE)** -- pore-forming toxin; absent from most non-typhoidal Salmonella; strongest discriminator
2. **STY1886 (CdtB)** -- typhoid toxin component; rarely found in NTS
3. **STY1479 (YncE)** -- ATP-binding protein; identified as potential chronic carrier marker

These are recommended for use "in combination" for multiplex serosurveillance approaches.

## Data Quality Notes

### Strengths
- Large cohort (n=486 patients + 49 community controls) with blood culture-confirmed diagnoses
- Comprehensive panel of 17 antigens, allowing direct head-to-head comparison of Vi vs. protein antigens
- Both S. Typhi and S. Paratyphi A patients included
- Culture-negative febrile controls provide a clinically relevant comparator (febrile patients without enteric fever)
- Afebrile community controls provide a population baseline
- Nonparametric statistical methods appropriate given the data distributions
- Same endemic setting (Kathmandu) as the Nepal sub-cohort in Aiemjoy 2022, facilitating cross-study comparison

### Limitations
- **Short follow-up (3 months only)** -- cannot characterize waning beyond the early convalescent phase
- **No fitted waning parameters** -- purely descriptive study; no quantitative decay rates extractable
- **Single site** (Kathmandu, Nepal) -- findings may not generalize to all endemic settings
- **Community controls were single time-point only** -- no longitudinal baseline data
- **ELISA units are relative** (normalized to Nepalese plasma pool); not directly comparable to Aiemjoy EU without calibration
- **No individual-level quantitative data tables** -- results presented only as boxplots and LOESS curves; median/IQR values must be estimated from figures
- **Cross-reactivity concerns** acknowledged: O2 and O9 cross-react between S. Typhi and S. Paratyphi A (shared O12 antigen), and O-antigens are present in invasive non-typhoidal Salmonella
- **Treatment effects not analyzed** -- patients received either gatifloxacin or ofloxacin; potential treatment effects on immune kinetics not assessed
- **Age range relatively narrow**: mean ~16-18 years; very few young children, so results may not apply to the <5 age group most relevant for TCV policy

### Key Page References
- Study design: pp.833-835
- Table 1 (demographics): p.834
- ELISA methods: p.834
- Statistical methods: p.835
- Figure 1 (longitudinal LOESS trajectories): p.836
- Figure 2 (boxplots -- protein antigens over time): p.837
- Figure 3 (boxplots -- Vi, O2, O9 over time): p.838
- Figure 4 (Spearman correlation matrices): p.839
- Figure 5 (between-group comparisons -- protein antigens): p.840
- Figure 6 (between-group comparisons -- Vi, O2, O9): p.841
- Discussion (Vi assessment): p.842
- Discussion (recommended antigens): p.842

## Cross-References

- **Aiemjoy et al. 2022** (ref 10 in this paper) -- Cited for seroincidence estimation using anti-HlyE. Uses the same Nepal population setting (Kathmandu) and overlapping antigen panel (HlyE, Vi). Aiemjoy provides quantitative waning parameters; Mylona provides complementary qualitative evidence about Vi's poor performance.
- **Koirala et al. 2013** (ref 17) -- Parent treatment trial from which all plasma samples originate. Open-label RCT of gatifloxacin vs. ofloxacin for uncomplicated enteric fever at Patan Hospital, Kathmandu, Nepal. ISRCTN63006567.
- **Tran Vu Thieu et al. 2017** (ref 16) -- Evaluation of purified S. Typhi protein antigens for serological diagnosis of acute typhoid fever; describes the ELISA methodology used here.
- **Watson et al. 2017** (ref 19) -- Cross-sectional Vi seroepidemiological survey in Fiji showing poor seroconversion against Vi; cited in discussion of Vi limitations.
- **Charles et al. 2013** (ref 15) -- Identification of immunogenic S. Typhi antigens in chronic biliary carriers; YncE identified as potential carrier marker.
- **Darton et al. 2017** (ref 31) -- Novel serodiagnostic signatures using S. Typhi proteome array; identified HlyE and CdtB as discriminatory.
- **Andrews et al. 2019** (ref 32) -- Plasma IgA responses against 2 S. Typhi antigens identify typhoid patients.

## Fit Role Assessment

**[ASSISTANT-PROPOSED]** **SUPPORTING** for waning characterization -- this paper does not provide quantitative waning parameters but provides critical qualitative evidence about the behavior of anti-Vi IgG following natural S. Typhi infection. Its primary contribution to the metastudy is threefold:

1. **Independent confirmation that anti-Vi IgG shows negligible boosting from natural S. Typhi infection** relative to protein antigen responses (HlyE, CdtB, YncE) and O-antigen responses (O2, O9) in the same patients.

2. **Direct evidence that anti-Vi IgG cannot discriminate enteric fever patients from controls by month 3** (Figure 6A), consistent with either rapid return to baseline or no meaningful boosting in the first place.

3. **Quantitative Spearman correlations (Figure 4D)** showing anti-Vi IgG is essentially uncorrelated (rho < 0.2) with all other antigen responses at day 8 in S. Typhi patients, indicating Vi IgG is measuring something fundamentally different from the acute infection response.

## Key Extractions for Waning/Boosting Model

### Anti-Vi IgG -- Critical Qualitative Findings

| Finding | Source | Interpretation |
|---|---|---|
| Vi LOESS trajectory essentially flat in ST patients | Figure 1D, p.836 | Minimal boosting from acute infection |
| No significant Vi IgG change at any time point in ST group | Figure 3A (ST panel), p.838 | Friedman/Wilcoxon non-significant |
| Vi IgG NOT different between ST patients and controls at M3 | Figure 6A (Day 90 panel), p.841 | Response returned to population background by 3 months |
| Vi IgG modest significance at D8 only (ST vs. FCN) | Figure 6A (Day 8 panel), p.841 | Very transient, minimal differentiation |
| Vi correlations rho < 0.2 with all antigens in ST at D8 | Figure 4D, p.839 | Vi response decoupled from acute infection response |
| Authors state Vi "not suitable as a marker of enteric fever exposure" | Discussion, p.842 | Explicit conclusion about Vi unsuitability |

### Contrast: HlyE (STY1498) IgG as Positive Control

| Finding | Source |
|---|---|
| HlyE LOESS trajectory markedly elevated in ST/SPA vs. FCN | Figure 1C |
| HlyE IgG significantly elevated at ALL time points (D1-M3) vs. controls | Figure 5C |
| HlyE IgG **** (P < .0001) separation from CC and FCN at D1, D8, M1, M3 | Figure 5C |
| HlyE shows clear peak at D8 and decline by M3 (i.e., waning is detectable within 3 months) | Figure 2C, p.837 |

### Implications for Aiemjoy 2022 Vi IgG Decay Rate = 0

**[CRITICAL]** The Mylona findings provide strong independent support for interpreting the Aiemjoy Vi IgG decay rate = 0 as an artifact of signal-to-noise limitations rather than evidence of truly durable Vi immunity after natural infection. The argument chains as follows:

1. **Mylona shows Vi IgG barely responds to acute S. Typhi infection** (Figure 1D, Figure 3A ST panel): if there is no meaningful boost, there is nothing to wane, and any waning model will estimate decay rate = 0.

2. **Mylona shows Vi IgG is at population background by 3 months** (Figure 6A Day 90): if the signal has already merged with background noise by month 3, a model fitting over months 3-24+ (as in Aiemjoy) will see a flat line and estimate zero decay.

3. **Mylona shows Vi IgG is uncorrelated with the acute infection response at peak** (Figure 4D, rho < 0.2 in ST): the Vi IgG being measured in endemic populations appears to reflect cumulative lifetime exposure or cross-reactive background, not recent infection. A waning model fit to this signal is characterizing the waning of something other than the acute infection response.

4. **The same study, same patients, same assay platform shows robust and waning-detectable responses for HlyE** (Figure 1C, 2C, 5C): this rules out methodological artifact. The problem is specific to Vi, not to the study design.

5. **Both studies are from Kathmandu, Nepal**, making the populations directly comparable. The Mylona cohort predates Aiemjoy but the epidemiological setting is the same.

**Bottom line**: Mylona 2024 demonstrates that anti-Vi IgG is a fundamentally poor marker of acute enteric fever exposure in endemic settings -- not because Vi IgG "doesn't wane" (which would imply durable immunity), but because the signal-to-noise ratio is too low to characterize any waning that might occur. The Aiemjoy decay rate = 0 for Vi IgG should be treated as effectively censored/uninformative rather than as a biological parameter indicating permanent immunity.
