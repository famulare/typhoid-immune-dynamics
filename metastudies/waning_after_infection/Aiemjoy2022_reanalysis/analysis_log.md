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

## 2026-03-24: Session 2 — reinf_obs investigation and mixture EDA

### Step 20: reinf_obs variable investigation

Added `04_faceted_by_reinf_obs.png` — trajectories faceted by n_obs (columns) × reinf_obs status (rows). Immediately revealed a structural pattern: ALL reinf_obs=NaN longitudinal subjects have exactly 2 observations. All 3- and 4-obs subjects have reinf_obs=1.

Investigated thoroughly to determine what `reinf_obs` means. The `data/README.md` says: "Suspected reinfection observed during follow-up."

#### What reinf_obs is NOT

1. **NOT "reinfection was detected."** 716/1667 subjects have reinf_obs=1, but the paper reports only 37/1420 suspected reinfections total. The base rate is 2.6%, not 43%.

2. **NOT a country/cohort indicator.** reinf_obs=1 is present across all three countries at varying rates (Bangladesh 44%, Nepal 33%, Pakistan 59%). The paper describes one SEAP protocol with identical visit schedules for all Asian sites.

3. **NOT "included in the paper's analysis."** The paper reports Bangladesh=407, Nepal=543, Pakistan=399 subjects. The CSV has Bangladesh=407, Nepal=798, Pakistan=462 — 255 extra Nepal and 63 extra Pakistan subjects, no Ghana. reinf_obs=1 counts (181, 263, 272) don't match the paper Ns for any country.

4. **NOT purely "had ≥3 months follow-up."** 296 reinf_obs=NaN subjects DO have time data at visit 3 (≥3 months). The correlation is strong (93.6% of reinf=1 vs 31.1% of NaN have visit 3 time) but not clean.

5. **NOT used in the JAGS model.** The column is never referenced in v9na.data.r or v9na.model.jags.

#### Best-supported hypothesis

**reinf_obs=1 ≈ "subject had sufficient multi-antigen longitudinal data to be screened for reinfection"** using the paper's criterion (≥3-fold rise in ≥2 antigen-isotypes at ≥3 months post-onset).

Evidence: 90.5% of reinf=1 subjects have ≥2 antigen-isotype pairs spanning ≥3 months, vs only 39.1% of NaN subjects. This is the strongest single predictor found, though still imperfect. The remaining 9.5% discrepancy suggests an additional criterion in the upstream preprocessing code not included in the GitHub repo.

#### Alternative hypothesis

reinf_obs=1 could mean "this subject was actively enrolled in the systematic longitudinal follow-up arm" (vs subjects with incidental or incomplete follow-up). The structural correlation with having richer data would follow naturally.

#### Practical implication for reanalysis

The reinf_obs facet is essentially a **data-quality stratification**: reinf=1 subjects have richer multi-antigen follow-up data, reinf=NaN subjects have sparser data (≤2 Vi observations, limited antigen coverage). For the waning analysis, the reinf=1 group is the better-characterized set. The NaN group's trajectories appear flatter and more compressed, which could reflect either genuinely different biology (non-responders?) or simply measurement limitations from shorter follow-up.

#### CSV vs paper subject count discrepancy

The CSV contains 1667 subjects; the paper reports 1420 (407 Bangladesh, 543 Nepal, 399 Pakistan, 71 Ghana). The CSV has no Ghana subjects but 255 extra Nepal and 63 extra Pakistan subjects. These extras may be subjects enrolled after the paper's analysis cutoff or subjects excluded from the paper for reasons not documented in the repo. The upstream preprocessing code is not in the published GitHub repository.

### Step 21: Reinfection algorithm implementation (07_reinf_obs_investigation.py)

Applied the Aiemjoy reinfection detection algorithm ourselves to definitively test the hypothesis that reinf_obs flags algo output. Algorithm: ≥3-fold increase in ≥2 antigen-isotype combinations at visits ≥3 months from fever onset, with absolute difference ≥1 EU.

**Results:**
- Our algo flags 43/1667 suspected reinfections (paper reports 37 — close, minor implementation differences likely in "consecutive" vs "all" visit pairs)
- Among 191 longitudinal Vi subjects: 10 algo-reinfected, 181 algo-clean

**reinf_obs vs algo cross-tab (longitudinal Vi):**

|  | algo=clean | algo=reinfected |
|---|---|---|
| reinf_obs=1 | 145 | 4 |
| reinf_obs=NaN | 36 | 6 |

**reinf_obs is definitively NOT the algo output.** 6/10 algo-reinfected Vi subjects have reinf_obs=NaN — the opposite of what you'd expect. The cross-tab shows essentially no relationship.

**The algo does identify a real biological signal:**
- Algo-reinfected Vi subjects: 9/10 rising, median FC=1.53
- Algo-clean Vi subjects: 91 declining / 90 rising, median FC=0.99

These 10 subjects' rising Vi trajectories are likely reinfection-driven boosts, useful for the mixture analysis.

**Conclusion:** reinf_obs is an opaque data-provenance flag from upstream preprocessing not in the repo. We cannot determine its exact meaning. Moving on — it will not be used in the reanalysis. The reinfection algo classification is independently useful for understanding latent structure.

Moved reinf_obs faceted plot from 04 to 07. Script 04 now only contains slope analysis, direction facet, and spaghetti plots.

### Step 22: Folder reorganization

Refactored `reanalysis/` into logical subfolders:
- `data_prep/` — 01_prepare_vi_data.py
- `validation/` — 02, 03 (figure reproduction, PDF validation)
- `eda/` — 04, 05, 06 (slopes, noise, descriptive summary)
- `reinf_obs/` — 07 (reinf_obs investigation)
- `mixture_eda/` — 08+ (fold-change mixture models)

Each script's OUT_DIR updated to its subfolder. Scripts and their output PNGs live together.

### Step 23: Fold-change mixture model EDA (mixture_eda/08_fold_change_mixture.py)

Fitted four mixture models to the 191 longitudinal Vi IgG log2(fold change) values. All models enforce σ₂ = √(σ₁² + σ_bio²) so the signal component inherits measurement noise and σ_bio captures excess biological variance.

#### Model comparison

| Model | k | LL | AIC | BIC | CV | σ_bio |
|---|---|---|---|---|---|---|
| 1. Two Gaussians | 5 | -237.3 | 484.5 | 500.8 | 0.201 | 0.919 |
| 2. Gauss + Skew-Normal | 6 | -237.3 | 486.5 | 506.0 | 0.201 | 0.919 |
| 3. Signal time-dep (log2 Δt) | 6 | -237.1 | 486.2 | 505.7 | 0.199 | 0.916 |
| **4. Teunis power-law (τ-ratio)** | **6** | **-234.6** | **481.2** | **500.7** | **0.211** | **0.856** |

**Step 4 wins** (ΔAIC = 3.3 vs Step 1, same parameter count).

#### Step 1: Two-Gaussian mixture (MLE)

Noise: μ₁=-0.218, σ₁=0.406, π=0.425. Signal: μ₂=0.537, σ₂=1.005.

**Implied assay CV ≈ 0.20** — first data-driven estimate for Vi IgG. Derived from the noise component σ₁ via σ_per_meas = σ₁/√2 → CV. Consistent across all four models (0.199-0.211).

89/191 (47%) classified as responders. Responders have lower starting EU (330 vs 589 for noise). The noise group's higher starting EU and negative μ₁ reflect waning from elevated levels — not a statistical artifact (there is no selection on the first measurement), but actual Vi IgG decay. This observation is consistent with Vi as a correlate of protection: lower baseline Vi → more susceptible → infection → fold rise.

#### Step 2: Gaussian + Skew-Normal

Skew α = 0.000 — converged at the Gaussian boundary. Right-tail asymmetry undetectable at n=191. AIC worse (486.5).

#### Step 3: Signal time-dependent (log2 Δt)

Signal: N(μ₂ − β·log2(Δt), σ₂). β = -0.052 (wrong sign for waning). The log2(Δt) covariate doesn't capture the waning physics because it ignores when the observations occurred relative to peak antibody. AIC worse (486.2).

#### Step 4: Teunis power-law mixture (τ-ratio covariate)

Uses the Teunis power-law decay covariate z = log2(τ₁/τ₀) where τ = max(t − t_peak, 1 day). t_peak = 15 days from fever onset, based on HlyE IgG peak timing (15.6d in Aiemjoy Table S1 — the best-identified antigen in the same subjects/platform). Vi IgG's own t1=2.9d is unreliable.

**Noise component:** N(−α₁·z, σ₁) — no intercept, FC=1 at τ₁=τ₀ by construction.
- α₁ = 0.025 (power-law waning exponent)
- Predicted FC at t₀=14→t₁=365d: 0.86 (14% decline)

**Signal component:** N(μ_bio − α₂·z, σ₂)
- α₂ = -0.111 (negative: signal FC increases with time ratio — mixture of boosting/waning)
- μ_bio = -0.093 (near zero)

**Development notes on the Teunis covariate:**
- First tried the full Teunis form with β estimated: (1+β·τ)^(−α). β diverged to ∞, collapsing to the pure power-law τ^(−α). The data cannot separate β from α in fold-change data — this requires absolute antibody levels.
- The τ-ratio covariate z = log2(τ₁/τ₀) is the natural fold-change prediction from the power-law limit and has no free rate parameter.
- τ is floored at 1 day (not 0) to handle subjects observed before t_peak without singularities.

#### Key findings

1. **Vi IgG assay CV ≈ 0.20** — robust across all four models. First data-driven estimate, consistent with the range assumed in Step 9/16.

2. **Two-population structure confirmed** — ~40-46% noise (waning-dominated), ~54-60% signal (boost-dominated). Split is consistent across models.

3. **The Teunis power-law covariate matters** — using log2(τ₁/τ₀) instead of log2(Δt) improves AIC by 5 points (comparing Step 4 vs Step 3 at same k=6). Accounting for when observations fall relative to peak antibody is important.

4. **α₁ = 0.025 is a detectable waning exponent** — small but captures the systematic decline in the noise component. The noise group's negative fold change IS waning, not a statistical artifact. Higher starting EU → further from equilibrium → more waning.

5. **Vi as correlate of protection signal** — lower starting Vi IgG is associated with the responder component (median 394 vs 546 EU). Consistent with Jin's correlate-of-protection curve.

6. **σ_bio ≈ 0.86** — biological heterogeneity in the signal component is ~2× measurement noise (0.86 vs 0.43). This is the variation in boost magnitude, waning stage, and reinfection across subjects.

7. **Fold-change analysis reaches its limit** — per-subject trajectory modeling (Stan) is needed to properly estimate waning rates and separate boost timing from decay.

### Step 24: Threshold rules vs mixture model (mixture_eda/09_threshold_vs_mixture.py)

Compared the Step 4 mixture P(responder) to simple fold-change threshold rules. All 191 subjects are blood-culture-confirmed enteric fever: 156 S. Typhi, 35 S. Paratyphi A.

#### Serovar-blind mixture assignments

The mixture model assigns Typhi and Paratyphi subjects to the responder component at **identical rates** (~54% each). The model sees no serovar information — it operates on unlabeled Vi IgG fold changes. The equal assignment means the Vi IgG fold-change distributions of confirmed Typhi and Paratyphi patients are indistinguishable at the population level.

#### Threshold sensitivity on confirmed infections

Every Typhi subject had a confirmed infection that should in principle boost Vi IgG. Detection rates:

| Method | Typhi detected (of 156) | Paratyphi "detected" (of 35) |
|---|---|---|
| Mixture P>0.5 | 85 (54%) | 19 (54%) |
| Vi FC ≥ 2× | 29 (19%) | 9 (26%) |
| Vi FC ≥ 3× | 12 (8%) | 5 (14%) |
| Vi FC ≥ 4× | 6 (4%) | 2 (6%) |

The 3× threshold catches only **8% of confirmed S. Typhi infections** via Vi IgG alone. The threshold rules are strictly nested inside the mixture classification — every subject above any threshold is also classified as a mixture responder, but 84% of mixture-responders fall below the 3× threshold. The sensitivity curve (Panel D) shows the 3× rule catches ~16% of mixture-identified responders; the threshold must drop to ~1.3× to catch 80%.

#### Paratyphi subjects are NOT clean negative controls

Initial framing treated Paratyphi subjects as Vi-negative controls (S. Paratyphi A doesn't express Vi antigen). However, in a high-endemic setting like Nepal, Paratyphi-confirmed subjects may have had prior or concurrent **asymptomatic or paucisymptomatic S. Typhi infections**. The equal Vi boosting rate in Typhi and Paratyphi groups could reflect genuine serologic evidence of unconfirmed Typhi exposure, not false positives. This is consistent with high rates of subclinical typhoid transmission and is relevant to the cohort incidence modeling question.

#### Implications

1. **Vi IgG is nearly useless as a single-antigen seroincidence marker at the individual level.** The 1.5× median boost is buried in noise. Aiemjoy's actual seroincidence method uses all 7 antigen-isotypes jointly, not Vi alone — and for good reason.

2. **The mixture model captures population-level signal** (CV estimate, two-component structure, waning exponent) but **cannot identify Vi-specific boosting** at the individual level. Serovar is invisible to it.

3. **Asymptomatic Typhi infections may be common** in this cohort. The equal Vi boosting rates across serovars, combined with Nepal's known typhoid endemicity, suggest ongoing Typhi exposure regardless of which serovar caused the index blood-culture-confirmed episode. This is directly relevant to the dynamic modeling question about population-level significance of natural Vi boosting.

### Step 25: Truncated signal component and limits of fold-change analysis

Attempted to fix the non-monotonic P(responder) in the left tail (where the broad signal component's tail exceeds the noise component, causing P(responder) to rise for very negative fold changes).

#### Truncation approach

The biological boost B₀ at t_peak should always be ≥ 0 on the log2(FC) scale (a genuine boost produces FC ≥ 1). So the signal component is:

- B₀ ~ TruncNorm(μ_bio, σ_bio, lower=0) — initial boost, truncated at 0
- Observed: X = B₀ − α₂·z + ε — after waning and measurement noise
- Signal PDF: truncnorm_conv_pdf(x + α₂·z, μ_bio, σ₁, σ_bio)

The convolution of TruncNorm with Normal has a closed form: N(x; μ, σ₂)·Φ(m/v) / (1−Φ(−μ/σ_bio)). The truncation is on the **initial boost at t_peak**, not on the observed fold change — subjects who boosted then waned past baseline (FC < 1) are allowed.

Verified: analytical formula matches numerical convolution exactly. Left tail decays like N(0, σ₁) (measurement noise shape), as expected.

#### Key finding: truncation on initial boost requires knowing the baseline

The truncated model consistently fits worse than the untruncated Gaussian (AIC ~490 vs 484.5). The optimizer finds a different regime: noise component broadens (σ₁=0.71, CV≈0.36) to absorb left-tail subjects, π_noise increases to ~0.79, and the signal becomes a narrow population of clearly-boosted subjects. Constraining μ_bio ≥ 0 pushes it to the boundary at 0.

**The fundamental issue:** B₀ is the boost above *baseline*, but we don't observe baseline separately from the first measurement. The first Vi IgG measurement conflates (baseline + residual boost at the time of sampling). Subjects with high starting EU (median 1225 for the most negative FC bin) could have high baseline OR high residual boost — we can't distinguish these without a pre-infection measurement or a multi-timepoint trajectory model.

The strong negative correlation between starting EU and fold change (Spearman ρ = −0.62, p < 10⁻²¹) confirms this: high-start subjects decline, low-start subjects rise, and the fold-change summary cannot separate "waning from elevated baseline" from "waning from acute boost."

#### Also tested: full Teunis form with β estimated

Attempted (1+β·τ)^(−α) with β as a free parameter. β diverges to ∞ — the data can't separately identify β and α from fold changes alone (β requires absolute antibody levels). The pure power-law τ^(−α) is the correct limit for fold-change data.

#### Conclusion: definitive wall for fold-change analysis

The fold-change mixture models have extracted what they can:
1. **CV ≈ 0.20** (robust across all models)
2. **Two-population structure** (~40-55% noise/waning, ~45-60% signal/boosting)
3. **α₁ ≈ 0.025** waning exponent (from the untruncated Teunis model)
4. **Vi as correlate of protection signal** (lower starting EU → more likely responder)
5. **Threshold rules catch only ~8% of confirmed infections** (3× rule)
6. **Paratyphi subjects indistinguishable from Typhi** in Vi fold changes — possible evidence of widespread asymptomatic Typhi co-exposure in this endemic population

To go further requires per-subject trajectory modeling (Stan) that jointly estimates:
- Individual baseline (pre-infection level)
- Boost magnitude and timing
- Power-law waning rate
- Measurement noise (informed by CV ≈ 0.20)
- Latent class (responder vs non-responder)

The fold-change EDA provides strong priors for this model: CV, the population mixture proportions, the waning exponent range, and the correlate-of-protection relationship.

### Step 26: Fold-rise-informed mixture model (Step 5 in 08_fold_change_mixture.py)

Replaced the free μ_bio boost parameter with the fold-rise model from `scratch/cohort_incidence_model_proof_of_concept.R`:

fold_rise(eu_start) = 10^(μ₀ · (1 − log10(eu_start) / log10(CoP_max)))

The signal component's expected boost is now a **decreasing function of starting titer** — subjects with high pre-challenge Vi IgG get smaller boosts (ceiling effect). This mechanistically explains the ρ=-0.62 correlation we observed.

Signal: B₀ ~ TruncNorm(log2(fold_rise(eu_start)), σ_bio, ≥0), then X = B₀ − α₂·z + ε

#### Results

| Model | k | AIC | ΔAIC vs Step 1 |
|---|---|---|---|
| Steps 1-4 | 5-6 | 484-490 | baseline |
| **Step 5: Fold-rise** | **7** | **419.8** | **+64.7** |

**ΔAIC = 65 over Step 1** — overwhelmingly better with only 2 extra parameters.

Fitted parameters:
- **μ₀ = 3.32** — boost intensity. Predicted fold rises: 23× at eu_start=50, 4.6× at 200, 1.6× at 500, 1× at 1000 EU
- **CoP_max = 753 EU** — titer ceiling where boost → 1×. Lower than the cohort model default (3162 EU), suggesting Vi IgG boost saturates earlier than assumed
- **α₁ = -0.023** — noise waning exponent (near zero, noise component essentially flat with time)
- **α₂ = 0.107** — signal waning exponent (positive, correct sign: responders' fold changes decline with time ratio)
- **σ_bio = 0.000** — collapsed to zero! The fold-rise model explains ALL biological heterogeneity through the eu_start dependence
- **CV ≈ 0.30** — higher than Steps 1-3 (0.20). The fold-rise model absorbs what was previously called "noise"
- **π_noise = 0.45** — 55% noise, 45% signal; but component assignments shifted (64% classified as responders at P>0.5 because the signal component now correctly predicts small boosts for high-titer subjects)

#### Interpretation

The fold-rise model does two critical things:
1. **Explains the eu_start → FC correlation mechanistically** — it's not regression to the mean, it's the immunological ceiling effect where high-titer subjects boost less
2. **Enables correct classification of high-titer subjects with small FC** — these were misclassified as noise in Steps 1-4 but are now correctly identified as responders who had small expected boosts

The σ_bio = 0 finding means individual variation in boost magnitude is entirely captured by the fold-rise model's dependence on starting titer. No additional subject-level heterogeneity is needed beyond measurement noise.

CoP_max = 753 EU is a calibration target for the cohort model. It represents the Vi IgG level at which natural infection produces no additional boost — effectively the immune saturation point for Vi.

#### Caution: error-in-covariates and the Step 5 AIC

Despite the ΔAIC=65, Step 5's parameter estimates are biased and should not be used directly:

1. **Error-in-covariates.** eu_start is a noisy measurement of something that isn't even the right quantity. The fold-rise model requires the *pre-infection* titer, but eu_start is a *post-infection* observation taken days to months after fever onset. Conditioning on a noisy covariate inflates apparent explanatory power because the model overfits to measurement error in eu_start.

2. **σ_bio = 0 is a red flag.** The model claims ALL individual variation in boost magnitude is explained by eu_start. This is too good — it's absorbing eu_start measurement error into the fit rather than identifying genuine biological heterogeneity.

3. **CoP_max = 753 EU is biased downward.** Fitting from noisy, post-infection first-observations rather than true pre-challenge titers systematically underestimates the ceiling.

**Decision: Step 4 (untruncated Teunis τ-ratio) is the working model** for fold-change classification and parameter estimation. It uses only timing information (measured precisely), gives interpretable parameters (α₁, CV), and its limitations are well-characterized.

Step 5 is retained as a proof-of-concept that the fold-rise model structure explains the eu_start→FC correlation mechanistically. Its parameters (μ₀, CoP_max) should be re-estimated in the Stan individual-level trajectory model where a proper observation model handles measurement timing and error-in-covariates.

---

## Infrastructure notes

- `pyproject.toml` created for `uv sync` with extraction extras (PyMuPDF, matplotlib, scipy)
- All Python scripts run from project root via `python3 <script>` (or `uv run python`)
- Directory renamed from `waning` → `waning_after_infection`
- Reanalysis scripts organized into subfolders: `data_prep/`, `validation/`, `eda/`, `reinf_obs/`, `mixture_eda/`
- CSV data from Aiemjoy GitHub repo is the authoritative data source; PDF extraction preserved in git history as methodology validation
