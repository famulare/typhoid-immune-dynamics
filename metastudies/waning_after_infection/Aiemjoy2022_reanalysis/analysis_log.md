# Aiemjoy 2022 Reanalysis — Analysis Log

## 2026-03-24: Session 1 — Extractions, PDF digitization, and initial reanalysis

### Step 1: Aiemjoy 2022 extraction for waning metastudy

Created `metastudies/waning_after_infection/extracts/Aiemjoy_2022.md` — full extraction of Aiemjoy et al. 2022 (Lancet Microbe 3:e578-87) focused on antibody waning parameters.

Key content extracted:
- All fitted parameters from Table S1 (y0, y1, t1, decay rate α, shape r) for 5 antigen-isotypes (HlyE IgA/IgG, LPS IgA/IgG, Vi IgG), 3 age strata + overall
- Anti-Vi IgG: decay rate = 0 (95% CrI 0–0.4), duration above baseline ≥32 months, boosting ratio only 1.5× population mean
- Duration of elevation ordering: HlyE IgG (29 mo) > LPS IgG (14 mo) > HlyE IgA (11 mo) > LPS IgA (3 mo) > Vi IgG (≥32 mo, reported as no decay)

### Step 2: Documented Teunis 2016 decay model functional forms

Read Teunis et al. 2016 (Epidemics 16:33-39) and added exact equations to the Aiemjoy extraction:
- Decay ODE: y'(t) = −ν·y(t)^r
- Closed-form: y(t > t1) = y1·(1 + (r−1)·y1^(r−1)·ν·(t−t1))^(−1/(r−1))
- Equivalent power-function: y(τ) = y1·(1 + βτ)^(−α) with Gamma distribution interpretation
- Physical interpretation: power-function decay from Gamma-distributed heterogeneity in antibody production site decay rates

### Step 3: Developed censoring argument for Vi IgG "no waning"

Added 6-point argument to the Aiemjoy extraction that the Vi IgG decay rate = 0 is a signal-to-noise artifact:

1. Boosting ratio only 1.5× (vs 30.9× for HlyE IgG)
2. Baseline dominates trajectory (y1/y0 only ~2.4× in 5-15 age group)
3. Model structure right-censors the decay rate when signal-to-noise is poor
4. <5 age group (where y0 is low) has only n=12 from Nepal
5. Flat cross-sectional Vi age profile contradicts "no waning" — levels should accumulate
6. HlyE IgG DOES show age trends in same populations — positive control confirming the test works

Point 6 is the model-free clincher: works as a logical disjunction without requiring trust in any model.

### Step 4: Mylona 2024 extraction

Created `metastudies/waning_after_infection/extracts/Mylona_2024.md` — extraction of Mylona et al. 2024 (JID 229:833-44).

Key findings relevant to Vi waning:
- Vi IgG LOESS trajectory essentially flat in S. Typhi patients (Figure 1D)
- No significant within-group change for Vi at any timepoint (Friedman test)
- Vi IgG NOT different from community controls by month 3 (Figure 6A)
- Vi correlated poorly with all protein antigens (Spearman rho < 0.2 at day 8)
- Authors explicitly: "Vi is not suitable as a marker of enteric fever exposure"
- HlyE shows robust waning-detectable responses in same patients/assay — rules out methodological artifact

### Step 5: Figure S11 vector extraction from PDF

Created extraction pipeline in `Aiemjoy2022_reanalysis/fig_s11_extraction/`.

**Phase 1** (`phase1_extract.py`): Raw primitive extraction from supplement page 13 using PyMuPDF.
- 1977 drawing entries extracted (all vector, no rasterization)
- 62 multi-segment line paths identified as potential trajectories
- Axis labels rendered as paths (not text) — required hard-coded calibration

**Phase 2** (`phase2_extract_data.py`): Panel segmentation and coordinate conversion.
- Panel A: x-axis is LINEAR (not log), verified by 5-tick residuals < 0.04 pts. Y-axis log10.
- Panel B: age group color mapping identified: #f38efc = <5, #9133be = 5-15, #000080 = 16+
- Panel C data points are 2-bezier half-ellipses (not 4-bezier circles as initially expected)

Extraction results:
- Panel A: 141 individual trajectories + 1 overall model line
- Panel B: 261 data points (108 Kathmandu, 153 Kavre)
- Panel C: 500 data points (81 Vi IgG)
- **Validation**: Vi IgG boosting ratio from Panel C = median 1.46 (IQR 1.03–2.07) vs published 1.5 (1.0–2.1). Near-perfect match.

### Step 6: Initial trajectory slope analysis

Quick analysis of extracted Panel A trajectories:
- 69 declining, 71 rising, 1 flat — perfect 50/50 split
- Median log10 slope/day ≈ 0 (0.000004)
- Median fold change = 1.000
- **Conclusion**: no detectable directional waning signal. Pure noise.

Long-duration subset (>200 days, n=39): 12 declining, 27 rising — slight upward bias, possibly reinfections.

### Step 7: Jin 2017 extraction for waning metastudy

Created `metastudies/waning_after_infection/extracts/Jin_2017.md` — re-evaluation of Jin et al. 2017 (Lancet 390:2472-80) Vi-TCV CHIM trial through the lens of Vi waning.

Key findings:
- Vi-TT produces ~100-fold rise in anti-Vi IgG (5.6 → 562.9 EU/mL) vs natural infection's 1.5× boost
- Anti-Vi IgG is a continuous correlate of protection (OR 0.37 per log10 increase, adjusted)
- No direct waning data (only pre-vax and day 28 measured)
- Jin uses VaccZyme commercial assay — DIFFERENT from Aiemjoy in-house ELISA

Modeling framework (Mike confirmed):
1. Vi IgG wanes; Aiemjoy decay rate = 0 is a non-informative upper bound
2. Vaccine immunogenicity data should constrain waning rates
3. Jin Fig S3 correlate-of-protection curve bridges antibody levels to protection probability
4. Population-level significance of natural-infection Vi boosting is an open question requiring dynamic modeling

### Step 8: Panel A reanalysis — faceted trajectory plots

Created `fig_s11_reanalysis/panel_a_reanalysis.py` with:
- Reproduction of Panel A from extracted data
- Slope distribution analysis (histogram, fold-change, slope vs duration, slope vs starting EU)
- Spaghetti plot colored by observation duration
- Faceted grid: columns = # time points (2/3/4), rows = declining/rising

Trajectory composition: 98 with 2 points, 35 with 3 points, 8 with 4 points. Most are just two-point line segments.

Clarified that "time 0" = days since fever onset, NOT enrollment. No patient has a pre-infection baseline. First sample is at first clinical visit (prospective) or later (retrospective enrollees).

### Step 9: Measurement noise analysis

Created `fig_s11_reanalysis/measurement_noise_analysis.py`.

Vi IgG CV not reported in Aiemjoy Table S3. Used plausible range 0.15–0.30 based on same-platform HlyE/LPS CVs (0.13–0.31).

Key results:
- Observed SD of log2(FC) = 0.87, exceeds pure-noise prediction at all CVs (0.30–0.60)
- **Measurement noise alone cannot explain all observed variation** — real biological heterogeneity exists
- KS test rejects noise-only null at all CVs (p < 0.02)
- But excess variation is SYMMETRIC around zero — mixture of waning + reinfection + heterogeneity, not a clean directional signal

Population-level minimum detectable Vi IgG half-life (80% power):
- CV=0.15: 4.7 years
- CV=0.20: 3.5 years
- CV=0.25: 2.8 years
- CV=0.30: 2.4 years

**Conclusion**: This study design can only detect Vi waning if half-life < ~2.5–5 years. Anything longer is invisible.

### Step 10: Discovery of Aiemjoy GitHub repository

Scanned https://github.com/kaiemjoy/TyphoidSeroIncidence — discovered the actual individual-level data CSV.

`data/TypoidCaseData_github_09.30.21.csv`: 1,667 subjects, **341 with Vi IgG data** (all Nepal), **611 total Vi IgG measurements**. This is far richer than our 141 trajectories extracted from the figure.

JAGS model file confirms decay functional form and reveals that Vi IgG measurement precision (`prec.logy[7]`) is estimated from data, not hardcoded.

**Decision**: Crossover to using the actual CSV data. PDF extraction preserved in git history.

---

## 2026-03-24: Session 1 (continued) — Crossover to CSV data

### Step 11: Discovery and download of Aiemjoy GitHub data

Discovered https://github.com/kaiemjoy/TyphoidSeroIncidence contains the actual individual-level CSV data. Downloaded:
- `data/TypoidCaseData_github_09.30.21.csv` (1,666 subjects, 477 KB)
- `longitudinal_model/v9na.model.jags` (JAGS model definition)
- `longitudinal_model/graph-func.r` (ab() and serocourse() implementations)
- `longitudinal_model/v9na.data.r` (data loading, column mappings, zero handling)

**Key data facts**: 341 subjects with Vi IgG (all Nepal), 194 longitudinal (≥2 obs), 611 total Vi IgG measurements. Vi IgG is test index 7 in their model. Measurement precision estimated from data (not hardcoded). Repo has no LICENSE file; data provenance documented in `data/README.md`.

### Step 12: Data preparation (01_prepare_vi_data.py)

Melted wide-format CSV to long format for Vi IgG. Key findings:
- 159 S. Typhi + 35 S. Paratyphi longitudinal subjects
- Observation counts: 147 single, 130 with 2, 52 with 3, 12 with 4
- Time range: 0–1133 days (extended beyond figure's 540-day range)
- 3 zero values in dataset
- Age: only 8 subjects <5 years (Nepal cohort is predominantly adults)
- All 152 longitudinal subjects with reinf_obs data have reinf_obs=1 — meaning uncertain (needs investigation)

### Step 13: Figure S11A reproduction from CSV (02_reproduce_fig_s11a.py)

Successfully reproduced Figure S11A from the CSV. The S. Typhi-only spaghetti plot (n=159) matches the published figure well. The overall model curve (Table S1 parameters, decay rate=0) is flat, as published. Trajectory count discrepancy: 159 in CSV vs 141 in published figure — likely data freeze or exclusion criteria differences.

### Step 14: PDF extraction validation (03_validate_pdf_extraction.py)

Matched 74/141 PDF-extracted trajectories to CSV subjects (52% match rate, strict tolerance).

Point-level accuracy (n=165 matched points):
- **ELISA units**: median error 0.0%, IQR [0.0%, 0.1%] — essentially perfect
- **Time**: median error -3.0 days, IQR [-4.0, -2.0] — small systematic offset from x-axis calibration

**Conclusion**: PDF vector extraction methodology validated. The coordinate calibration was accurate to within the expected precision limits.

### Step 15: Slope analysis on full CSV (04_slope_analysis.py)

Results with 191 trajectories (vs 141 from PDF):
- 92 declining, 99 rising (48/52 split — consistent with PDF result of 69/71)
- Median fold change: 1.053 (PDF: 1.000)
- Median log10 slope/day: 0.000063 (essentially zero)
- Long-duration (>200 days, n=57): 19 declining, 38 rising — upward bias persists, likely reinfections
- S. Typhi and S. Paratyphi show similar patterns

Generated: slope distributions, faceted trajectories (npoints × direction), spaghetti by serovar.

### Step 16: Measurement noise analysis on full CSV (05_measurement_noise.py)

Observed SD of log2(FC) = 0.888 (PDF: 0.871 — consistent).

Minimum detectable Vi IgG half-life (population, 80% power):
- CV=0.15: 6.2 years
- CV=0.20: 4.6 years
- CV=0.25: 3.7 years
- CV=0.30: 3.1 years

Slightly longer than PDF-based estimates due to longer median duration (139 vs 123 days). Core conclusion unchanged: Vi waning with half-life < ~3-6 years cannot be detected with this study design.

### Step 17: Descriptive summary (06_descriptive_summary.py)

Dataset characterization: observation patterns, time coverage, Vi IgG distributions by visit/age/serovar.

### Step 18: Cleanup

Removed `fig_s11_extraction/` and `fig_s11_reanalysis/` from working tree. These are preserved in git history at commit 720ea47 and were validated against CSV ground truth in Step 14.

---

### Step 19: JAGS model formal extraction (their_code/model_extraction.md)

Full reverse-engineering of the Aiemjoy/Teunis JAGS model following the scientific model extraction template. 420 lines covering: mathematical specification, all priors with exact hyperparameters, observation model, hierarchical structure, evidence flow analysis, and red flags.

Key findings from model code analysis:

#### (1) Estimating the CV

The JAGS model (line 10): `logy[subj,obs,test] ~ dnorm(mu.logy[subj,obs,test], prec.logy[test])`

`prec.logy[test]` is a **single precision parameter per antigen-isotype** (not per subject, not per lab). It's estimated from data with prior `Gamma(4.0, 1.0)` — prior mean precision = 4, implying prior SD(log y) = 0.5, which corresponds to CV ≈ 53%.

**Critically: this is NOT just the assay CV.** It captures everything the per-subject trajectory model fails to explain:
- Assay measurement noise
- Within-subject biological fluctuation
- Model misspecification (the trajectory shape doesn't perfectly fit)

So `prec.logy[7]` for Vi IgG absorbs all residual variance on the log scale. If the model poorly fits Vi trajectories (because the true dynamics are a mixture of flat and waning), the "measurement noise" estimate inflates to compensate. **This means the model conflates genuine Vi heterogeneity with measurement noise.**

#### (2) Reinfection handling

**There is NO reinfection model in the JAGS code.** The model assumes one infection → rise → decay per subject, period.

The paper says reinfections were excluded from the data before fitting. But `v9na.data.r` doesn't do this exclusion — it loads the full CSV and sends everything to JAGS. So either:
- The CSV already has post-reinfection observations removed (the `reinf_obs` column flags this somehow)
- The exclusion was done in a preprocessing step not in this repo

The `reinf_obs` column is suspicious: 152/194 longitudinal Vi subjects have `reinf_obs = 1`. That can't mean "reinfection was detected" (the paper reports only 37/1420 total). It more likely means "this subject was evaluable for reinfection" (followed ≥3 months). The actual post-reinfection data trimming probably happened upstream.

#### (3) Latent structure / grouping

**No mixture model, no latent classes.** The hierarchical structure is:

`par[subj,test,1:5] ~ dmnorm(mu.par[test,], prec.par[test,,])`

All subjects for a given test are drawn from a **single multivariate normal** on the 5 log-transformed parameters. The full 5×5 covariance matrix (via Wishart prior, df=20) captures continuous correlations between parameters (e.g., high-peak subjects may have different decay rates). But there is **no mechanism for discrete subpopulations**.

This means: if Vi IgG trajectories are a mixture of (a) genuine responders who boost and then wane, and (b) subjects where the measured values are just endemic background fluctuation, the model averages over both. The posterior for decay rate will be wide, centered near zero, because the mixed-population average is dominated by the non-responders. The model literally **cannot distinguish "Vi doesn't wane" from "Vi wanes but most of these trajectories are noise."**

Also notable: **each antigen-isotype is modeled independently.** `par[subj,test,1:5]` has separate population distributions per test. Vi IgG borrows no strength from HlyE or LPS for the same subjects. The only sharing is observation times.

Age dependence: **not in the model.** The commented-out filter lines in `v9na.data.r` show they run the model separately on age-filtered subsets. No continuous age covariate.

#### (4) Total evidence flow

For Vi IgG (test 7):
- 611 observations from 341 subjects feed the likelihood
- These inform: 341 subject-level parameter vectors, 1 population mean vector (5D), 1 population precision matrix (5×5), and 1 residual precision scalar
- The 147 single-observation subjects contribute very weakly — with only 1 point, all 5 parameters are essentially drawn from the prior
- The 194 longitudinal subjects (2-4 points each) constrain their individual parameters, but with 2-4 points for 5 parameters per subject, each subject is underdetermined. The hierarchical prior does the heavy lifting.

One cute trick: line 95-97 adds a **ghost subject** with 3 observations at days 5, 30, 90 and NA antibody values. JAGS imputes these from the posterior predictive — it's a standard trick for getting model predictions at specific time points.

#### Why it's "not bad" — it's just not designed for this question

The model is well-crafted for its intended purpose: estimating **population-average antibody kinetics** for seroincidence calculation. For that, you want the marginal kinetic parameters, and the single-MVN hierarchical model gives you exactly that with proper uncertainty propagation.

But for asking "does Vi wane?", the model has a structural blind spot: it cannot separate the signal (genuine post-infection waning) from the noise floor (endemic background fluctuation) when the two overlap as heavily as they do for Vi. A mixture model or a model with an explicit "responder/non-responder" latent class would be better suited, but would need substantially more data per subject to identify — and with 2-4 points per subject, that's not feasible.

So the authors aren't wrong or careless — they're using an appropriate model for their question (seroincidence) that happens to give a misleading answer for our question (waning).

#### Additional model mechanics
- y1 = y0 + exp(par[2]) guarantees y1 > y0 (peak above baseline)
- shape = exp(par[5]) + 1 guarantees shape > 1 (faster-than-exponential initial decay)
- The inner term of the power-function stays positive when shape > 1 (no negativity bug)
- Ghost subject (nsubj+1) at days 5/30/90 with NA data: prior predictive device, zero impact on inference
- Zero replacement (0 → 0.01 before log) is pragmatic but arbitrary; 3 Vi IgG zeros affected

---

## Infrastructure notes

- `pyproject.toml` created for `uv sync` with extraction extras (PyMuPDF, matplotlib, scipy)
- All Python scripts run from project root via `uv run python <script>`
- Directory renamed from `waning` → `waning_after_infection`
- Reanalysis scripts numbered sequentially (01–06) in `reanalysis/`
- CSV data from Aiemjoy GitHub repo is the authoritative data source; PDF extraction preserved in git history as methodology validation
