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

## Infrastructure notes

- `pyproject.toml` created for `uv sync` with extraction extras (PyMuPDF, matplotlib, scipy)
- All Python scripts run from project root via `uv run python <script>`
- Directory renamed from `waning` → `waning_after_infection` and `fig_s11_extraction` → `Aiemjoy2022_fig_s11_extraction` → current structure for self-documentation
