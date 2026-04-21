# Joint Inference Plan: Typhoid Dose-Response Model Calibration

**Date**: 2026-03-17
**Authors**: Mike Famulare, Claude (Opus 4.6)
**Status**: Draft revised after Reviewer 2 feedback (see `reviewer2_response.md`)
**Purpose**: Complete specification of a Bayesian joint inference framework for calibrating a two-outcome typhoid dose-response model to all available CHIM data, with explicit treatment of nuisance parameters, exchangeability assumptions, and cross-era bridging.

---

## 1. Executive Summary

We propose a joint Bayesian model that fits two dose-response curves — **infection** and **fever** — to data from all Oxford (2011-2017) and Maryland (1960s-1970s) controlled human infection studies. The model uses the modified beta-Poisson form with immunity scaling, a delivery medium offset, and latent immunity for the Maryland cohort. We work in **bicarbonate-equivalent dose units** as the reference frame, with the milk-to-bicarbonate conversion as an explicit parameter.

The inference has three layers, each contributing information that the others cannot:
1. **Oxford data** bounds N50_fev from above (fever >50% at 10³ bicarb) and estimates immunity scaling (γ_fev) from Vi-vaccine contrasts (Jin Vi-TT/Vi-PS). Oxford alone is **underconstrained** for α — with only 2 dose levels 1 log apart, both in the flat upper tail, the dose-response shape is not identifiable from Oxford alone. **Oxford shedding is excluded from the infection likelihood** (see Section 2.6) because early antibiotic treatment truncates shedding detection, producing systematic downward bias that worsens at higher doses.
2. **Maryland data** is where α gets pinned: 6 orders of magnitude of dose range (10³-10⁹) spanning the full S-curve. Maryland also uniquely identifies the infection-fever conditional (Hornick Table 2: P(fever|infected) = 57% at 10⁷) and constrains the medium offset δ.
3. **Cross-era bridging** requires explicit assumptions about what biology is portable across 50 years (strain virulence, α, γ, functional form) and what is not (delivery medium, background immunity, outcome definitions). This is the load-bearing structure of the entire inference.

Total parameters: 6 biological + 4 nuisance + 1 overdispersion = 11 free parameters
Active calibration observations: ~25 binomial terms across ~8 non-overlapping study arms
Validation observations: 4 (Darton M01ZH09/Ty21a, non-Vi mechanisms)

---

## 2. Model Structure

### 2.1 Cascaded Beta-Poisson Model

The model has two stages, each a beta-Poisson:

**Stage 1 — Infection:**

$$
P_{\text{inf}}(D, \text{CoP}) = 1 - \left(1 + \frac{D}{\delta_{\text{medium}}} \cdot \frac{2^{1/\alpha_{\text{inf}}} - 1}{N_{50,\text{inf}}}\right)^{-\alpha_{\text{inf}} / \text{CoP}^{\gamma_{\text{inf}}}}
$$

**Stage 2 — Fever given infection:**

$$
P_{\text{fev|inf}}(D, \text{CoP}) = 1 - \left(1 + \frac{D}{\delta_{\text{medium}}} \cdot \frac{2^{1/\alpha_{\text{fev|inf}}} - 1}{N_{50,\text{fev|inf}}}\right)^{-\alpha_{\text{fev|inf}} / \text{CoP}^{\gamma_{\text{fev|inf}}}}
$$

**Observable marginal — Fever:**

$$
P_{\text{fev}}(D, \text{CoP}) = P_{\text{inf}}(D, \text{CoP}) \times P_{\text{fev|inf}}(D, \text{CoP})
$$

This is the product of two beta-Poissons — **not itself a beta-Poisson**. The marginal P(fever) has a richer shape family than a single beta-Poisson can produce.

**Parameters:**
- $D$ = administered dose in CFU (as reported by the study)
- $\delta_{\text{medium}}$ = delivery medium offset (milk-to-bicarb conversion; bicarb = 1, milk = δ > 1)
- $N_{50,k}$ = dose for 50% outcome probability in naive individuals, in bicarb-equivalent CFU
- $\alpha_k$ = heterogeneity parameter
- $\gamma_k$ = immunity scaling exponent
- $\text{CoP}$ = correlate of protection (1 = naive, >1 = some immunity)
- $k \in \{\text{inf}, \text{fev|inf}\}$ — 6 biological parameters total

**Convention**: All N50 values are in bicarb-equivalent units. For Oxford data, δ = 1. For Maryland data, δ > 1.

### 2.2 Why the Cascaded Form

The cascade parameterization — where $P(\text{fev|inf})$ is an explicit beta-Poisson rather than a derived ratio of two marginals — has important behavioral advantages:

**1. Better behavior at high immunity (low CoP).** In the marginal-ratio formulation ($P_{\text{fev}} / P_{\text{inf}}$, both beta-Poisson), the conditional approaches an ill-defined 0/0 limit as both curves approach zero at high CoP. In the cascaded form, $P(\text{fev|inf})$ is its own well-behaved function with its own asymptotic behavior.

**2. Can mimic stable protection not captured by CoP.** If cellular immunity (T cells, mucosal IgA) provides protection against fever that is not captured by the anti-Vi IgG CoP, this manifests as a *lower* $P(\text{fev|inf})$ that persists even at CoP ≈ 1. The cascaded form accommodates this naturally through $\gamma_{\text{fev|inf}}$ operating on a different mechanism than $\gamma_{\text{inf}}$. The marginal-ratio form would require the fever curve to track the infection curve's CoP dependence, which it may not.

**3. The marginal P(fever) is a richer function.** The product of two beta-Poissons can produce shapes that no single beta-Poisson can: floors, shoulders, and differential curvature at different dose ranges. This matters for fitting the Maryland fever data, where the curve shape may reflect both the infection process and the fever-given-infection process operating at different effective thresholds.

**4. Identifiability is preserved.** Both $P_{\text{inf}}$ and $P_{\text{fev}} = P_{\text{inf}} \times P_{\text{fev|inf}}$ are observable as group-level marginal attack rates (shedding and fever respectively). The conditional $P_{\text{fev|inf}}$ is directly observable at Hornick Table 2 (16/28 at 10⁷ milk). The likelihood structure is the same as the two-marginal formulation: binomial observations on shedding rates, fever rates, and the one conditional rate.

**Prior art**: This cascaded beta-Poisson-given-beta-Poisson structure was used successfully for COVID-19 dose-response modeling (infection vs severe disease), where it captured the empirical observation that vaccines provided stronger protection against severe disease than against infection — analogous to the "leaky protection" pattern expected here.

**Hornick Table 2 likelihood**: The conditional observation enters directly:

$$y_{\text{fev|inf}} = 16 \sim \text{Binomial}(28, P_{\text{fev|inf}}(10^7/\delta, \text{CoP}_{\text{md}}))$$

The infection observation: $y_{\text{inf}} = 28 \sim \text{Binomial}(30, P_{\text{inf}}(10^7/\delta, \text{CoP}_{\text{md}}))$. No double-counting.

**Darton cross-tabulation (from S1 individual-level data)**: The Darton S1 extraction provides individual-level infection × fever crossover for all 91 per-protocol subjects. For the Placebo group (n=30): 13 stool+/TD+, 6 stool+/TD-, 7 stool-/TD+, 4 stool-/TD-. This enables a joint bivariate likelihood for Darton (per Reviewer 2 Major Concern 2), rather than two independent marginal binomials. Cross-tabulations at fever thresholds ≥38°C and ≥39°C are also available. See `analysis_data/darton_cross_tabulation.csv`.

For all other data where only marginal rates are observed (no individual-level crossover), the likelihood uses the product:

$$y_{\text{fev},j} \sim \text{Binomial}(n_j, P_{\text{inf}}(D_j, \text{CoP}_j) \times P_{\text{fev|inf}}(D_j, \text{CoP}_j))$$

### 2.3 CoP Mapping

The correlate of protection CoP is a scalar summary of the immune state at the time of challenge.

**For Oxford studies with measured anti-Vi IgG**: CoP is derived from the titer model (fit separately to longitudinal data). For this calibration, we treat the mapping CoP = $g(\text{anti-Vi IgG})$ as given from the titer model, propagating uncertainty if available.

**For Oxford naive groups**: CoP = 1 for subjects with anti-Vi below detection. For the ~30% of Oxford "naive" subjects with detectable baseline anti-Vi (>7.4 EU/mL), CoP > 1 as determined by the titer model. In practice, because we use group-level attack rates, we model the group-average probability:

$$
\bar{p}_j = f_{\text{naive},j} \cdot P_k(D, 1) + (1 - f_{\text{naive},j}) \cdot P_k(D, \overline{\text{CoP}}_{\text{Vi+},j})
$$

where $f_{\text{naive},j}$ is the fraction with undetectable anti-Vi and $\overline{\text{CoP}}_{\text{Vi+},j}$ is the average CoP of the anti-Vi-positive subgroup.

**For Oxford vaccinated groups**: CoP is determined by post-vaccination anti-Vi IgG GMT through the titer model.

**For Maryland studies**: CoP is a **latent variable** (anti-Vi not measured). See Section 2.4.

### 2.4 Maryland Immunity: Mixture Model

The Maryland prison population had heterogeneous background immunity. Direct evidence:
- Gilman 1977: H antibody ≥1:20 → 24% attack rate vs 61% without (P=0.02) at 10⁵
- Levine 1976: H antibody ≥1:20 → 22% attack rate vs 48% without (P=0.005) at 10⁵
- The Maryland dose-response curve is asymmetric (gentle lower tail, steep upper portion), consistent with a mixture of susceptible and immune subjects

We model the Maryland population as a two-component mixture:

$$
P_{\text{maryland},k}(D) = \pi_{\text{susc}} \cdot P_k(D/\delta, \text{CoP}_{\text{susc}}) + (1 - \pi_{\text{susc}}) \cdot P_k(D/\delta, \text{CoP}_{\text{imm}})
$$

Where:
- $\pi_{\text{susc}}$ = fraction of susceptible (near-naive) subjects
- $\text{CoP}_{\text{susc}}$ = CoP of susceptible fraction (near 1, may be slightly >1)
- $\text{CoP}_{\text{imm}}$ = CoP of immune fraction (estimated)

The Gilman/Levine H-antibody data directly informs $\pi_{\text{susc}}$: in Gilman, 36/53 controls (68%) had H antibody <1:20 (susceptible) and 17/53 (32%) had H antibody ≥1:20 (partially immune). Note: 53 of 64 total controls had baseline H-antibody measured.

**Simplification option**: Set $\text{CoP}_{\text{susc}} = 1$ (truly naive), reducing the mixture to 2 parameters ($\pi_{\text{susc}}$, $\text{CoP}_{\text{imm}}$).

### 2.5 Outcome Definition Adjustment

Maryland and Oxford use different outcome definitions, and within Maryland, different studies use different fever thresholds:

| Study | Fever Definition | Treatment Trigger | Culture Required? |
|---|---|---|---|
| **Oxford (Waddington, Darton, Jin)** | ≥38°C (100.4°F) for ≥12h | At diagnosis | No (OR bacteremia) |
| **Hornick 1966** | ≥103°F (39.4°C) for 24-36h | ≈ definition | Implied |
| **Gilman 1977** | Fever + culture (any threshold) | ≥103°F/1d, ≥101°F/3d, ≥100°F/5d | Yes |
| **Levine 1976** | ≥101°F (38.3°C) + culture | ≥103°F/24h, ≥101°F/48h | Yes |

These definitions observe different thresholds on the **same** underlying biological process: the distribution of peak fever severity given infection. We model this as a **temperature-dependent definition sensitivity function** $\phi(T)$ calibrated from the Oxford multi-threshold data.

#### Empirical φ calibration from Oxford threshold ladder

The Oxford studies report attack rates at multiple fever thresholds on the same participants, providing a direct calibration of how the observed fever rate scales with temperature cutoff:

| Threshold | Waddington 10³ | Waddington 10⁴ | Jin control | Darton placebo |
|---|---|---|---|---|
| ≥37.5°C (99.5°F) | 11/20 = 55% | 14/20 = 70% | 20/31 = 65% | — |
| ≥38.0°C (100.4°F) | 10/20 = 50% | 13/20 = 65% | 17/31 = 55% | — |
| ≥38.5°C (101.3°F) | 8/20 = 40% | 11/20 = 55% | 14/31 = 45% | — |
| ≥39.0°C (102.2°F) | — | — | — | 9/30 = 30% |

Normalizing to the Oxford composite TD rate (fever ≥38°C/12h OR bacteremia) to get $\phi(T)$:

| Threshold | φ (ratio to composite TD) | Notes |
|---|---|---|
| ≥38.0°C (Oxford reference) | 0.71-0.85 | Varies by study (Jin: 0.71; Waddington 10⁴: 1.0) |
| ≥38.3°C (Levine ≈101°F) | **~0.60-0.75** | Interpolated between ≥38.0 and ≥38.5 |
| ≥38.5°C (≈101.3°F) | 0.58-0.85 | Jin: 0.58; Waddington 10⁴: 0.85 |
| ≥39.0°C (≈102.2°F) | **~0.45** | Darton placebo |
| ≥39.4°C (Hornick ≈103°F) | **~0.30-0.37** | Extrapolated; Hornick also requires 24-36h sustained duration |

**Key insight (2026-03-24)**: The same φ(T) function that bridges Oxford↔Maryland also resolves the Levine↔Hornick definition gap *without any extra parameters*. Each Maryland study simply evaluates φ at its own temperature threshold:

$$P_{\text{study},j}^{\text{obs}}(D) = P_{\text{fev}}(D_{\text{eff}}, \text{CoP}_j) \cdot \phi(T_j) \cdot g_j(\text{duration})$$

where:
- $\phi(T_j)$ is the temperature-dependent definition sensitivity, calibrated from Oxford
- $T_{\text{Hornick}} = 39.4°C$, $T_{\text{Levine}} = 38.3°C$, $T_{\text{Gilman}} \approx 38.3°C$ (Gilman uses tiered triggers but counts any fever + culture)
- $g_j(\text{duration})$ is an additional penalty for duration requirements (Hornick requires 24-36h sustained; Levine and Gilman have tiered duration thresholds; this is ≤1 and likely ~0.7-0.9 for Hornick)

**Study-specific φ values** (estimated from Oxford ladder, no free parameters):
- **Hornick**: $\phi_H \approx 0.30 \times g(\text{24-36h duration}) \approx 0.25$
- **Levine**: $\phi_L \approx 0.65$ (≥101°F + culture; no strict duration requirement)
- **Gilman**: $\phi_G \approx 0.65$ (similar to Levine: fever + culture, tiered treatment triggers)

**Predicted Levine/Hornick ratio**: $\phi_L / \phi_H \approx 0.65/0.25 = 2.6\times$. This predicts Levine control attack rates should be ~2.6× Hornick at the same dose and immunity level. Observed: Levine pooled controls at 10⁵ ≈ 40% vs Hornick at 10⁵ ≈ 28%, ratio ≈ 1.4×. The predicted ratio is higher than observed, suggesting either (a) Levine's culture-confirmation requirement filters out some low-fever cases, or (b) the φ extrapolation overpredicts, or (c) real between-study immunity differences partially offset the definition gap. This calibration check can be performed during model fitting.

#### Dose-dependence of φ

$\phi$ is likely dose-dependent: at high doses, infections are more severe, so a higher fraction meets strict criteria. We parameterize:

$$
\phi(T, D_{\text{eff}}) = \phi_0(T) + (1 - \phi_0(T)) \cdot \frac{P_{\text{fev}}(D_{\text{eff}}, 1)^{\beta_\phi}}{1}
$$

where $\phi_0(T)$ is the minimum definition sensitivity at threshold T (at low effective doses where infections are mild) and $\beta_\phi$ controls how quickly it approaches 1 at high doses.

**Pragmatic recommendation (revised per Reviewer 2)**: Use dose-dependent $\phi$ as the default, since the Maryland data span a 6-log dose range while the Oxford threshold calibration comes from a single dose. A constant $\phi$ would systematically bias the fit at extreme doses, with the bias absorbed into $\alpha$ and $\delta$.

**Simplification option for initial exploration**: Fix $\phi$ as study-specific constants from the Oxford ladder (Hornick: 0.25, Levine: 0.65, Gilman: 0.65). Zero extra parameters. Check residual patterns for evidence of dose-dependence.

**Implementation note**: The Waddington alternative fever threshold data (Section 4.1.1, Table 2) provides the φ calibration at two dose levels (10³ and 10⁴). If φ values differ between these doses for the same threshold, that directly constrains $\beta_\phi$.

### 2.6 Oxford Shedding Exclusion (Treatment-Truncation Bias)

**Discovery (2026-03-24)**: Oxford stool shedding systematically underestimates true infection rate due to early antibiotic treatment. Participants were treated at typhoid diagnosis (median day 8-9), and ciprofloxacin rapidly clears gut colonization. Shedding that would have developed after treatment is never detected. The bias worsens at higher doses where diagnosis is faster:

| Study/Dose | Shedding | Diagnosis (fever/BC) | Shedding ≥ Diagnosis? |
|---|---|---|---|
| Waddington 10³ | 65% | 55% | Yes |
| Waddington 10⁴ | ~40% | 65% | **No — inverted** |
| Darton placebo ~10⁴ | 63% | 67% | **No** |
| Jin control ~10⁴ | 71% | 77% | **No** |

At every Oxford dose ≥10⁴, measured shedding < diagnosis. This violates the biological ordering P(infection) ≥ P(fever) and would produce systematic bias in N50_inf and P(fever|infection) estimates.

**Maryland shedding is unaffected.** Treatment was delayed until strict fever criteria (≥103°F/1 day, ≥101°F/3 days, ≥100°F/5 days), and shedding was measured 4-30 days post-challenge — ample time for detection before and after treatment.

**Decision**: Two model tiers:

**Tier 1 (primary model)**: Exclude all Oxford shedding from the infection likelihood. Oxford contributes **fever data only** plus vaccine contrasts on fever. The infection curve is identified entirely from Maryland shedding data through the δ bridge. This is conservative and avoids introducing bias. **However**, γ_inf is then unconstrained by Oxford vaccine data — only the Maryland mixture (H-antibody, not Vi) informs it. Since there are biological reasons to expect γ_inf < γ_fev (Vi immunity reduces fever severity more than colonization probability), and only the Jin Vi-vaccine shedding contrasts can distinguish them, the Tier 1 model may need to share γ = γ_inf = γ_fev, losing this biological nuance.

**Tier 2 (η-correction model — RECOMMENDED)**: Restore Oxford shedding with an explicit shedding detection probability η that corrects for treatment-truncation bias. See Section 2.7.

### 2.7 η-Correction for Oxford Shedding Detection (Tier 2 Model)

**Rationale**: The treatment-truncation bias has a specific causal structure. Shedding is missed only if it would have *started* after antibiotic treatment. The bias is therefore:

$$P(\text{detected shedding}_j) = \eta_j \cdot P_{\text{inf}}(D_j, \text{CoP}_j)$$

where $\eta_j = P(\tau_{\text{shed}} < \tau_{\text{treat},j})$ is the probability that shedding onset precedes treatment in study $j$. This is NOT a free-floating nuisance parameter — it has an observable anchor.

**Empirical constraint on η**: Waddington reports "Salmonella Typhi was subsequently shed by 18 of 24 (75%) of typhoid-diagnosed participants" and we established that ~100% of diagnosed participants at 10³ shed (11/11) vs ~54% at 10⁴ (7/13). Since true infection ≥ diagnosis (by definition), the diagnosed-who-shed fraction is a lower bound on η:

| Study | Dose | Diagnosed who shed / Diagnosed | η lower bound | Time to treatment, median [IQR] |
|-------|------|-------------------------------|---------------|--------------------------------|
| Waddington | 10³ | ~11/11 (100%) | ~1.0 | 14 [8.75-14] days |
| Waddington | 10⁴ | ~7/13 (54%) | ~0.54 | 9 [4-14] days |
| Darton | ~10⁴ | not separately reported | ~0.63/0.67 = 0.94? | similar to Jin |
| Jin | ~10⁴ | not separately reported | ~0.71/0.77 = 0.92? | 9.9 ± 3.35 days |

**Note**: The Darton and Jin "η" estimates from shedding/diagnosis ratios are NOT the same as the Waddington diagnosed-who-shed ratio. The shedding/diagnosis ratio mixes two effects: (1) some diagnosed people don't shed before treatment (η < 1), and (2) some shedders are undiagnosed (η correction applied to denominator). The Waddington disaggregation is the cleanest anchor.

**Implementation**:

*Option A — Single η per dose level (simplest)*:

$$y_{\text{shed},j} \sim \text{Binomial}(n_j, \, \eta(D_j) \cdot P_{\text{inf}}(D_j, \text{CoP}_j))$$

Model η as dose-dependent (since higher dose → faster diagnosis → less time for shedding):

$$\eta(D) = \eta_{\text{lo}} + (1 - \eta_{\text{lo}}) \cdot \exp(-\kappa \cdot D / N_{50,\text{inf}})$$

where $\eta_{\text{lo}}$ is the minimum detection probability at very high doses (≈0.5 based on Waddington 10⁴) and κ controls the dose-dependence. At low doses (D << N50), η → 1 (plenty of time before diagnosis). At high doses, η → η_lo.

This adds 2 parameters (η_lo, κ) but restores 8 observations (net gain of 6 degrees of freedom). This is better than dropping the observations.

*Option B — Study-level constant η (simpler)*:

Assign a single η per study-dose group, with a shared prior:

$$\eta_j \sim \text{Beta}(a, b)$$

Prior centered on ~0.75 (Waddington overall: 18/24 diagnosed shed) with moderate variance. Each study-dose group gets its own η_j drawn from this population distribution. This is a hierarchical nuisance parameter with ~3 effective parameters (hyperprior a, b + partial pooling).

*Option C — Fixed η from observed ratios (simplest, no extra parameters)*:

Fix η at observed shedding/diagnosis ratios per study-dose:
- Waddington 10³: η = 1.0 (shedding > diagnosis)
- Waddington 10⁴: η = 0.62 (bounded estimate ~8/13)
- Darton ~10⁴: η = 0.94 (19/20 shedders per diagnosed)
- Jin ~10⁴: η = 0.92 (22/24)

This adds zero parameters but bakes in assumptions. Good for initial exploration; too rigid for final inference.

**Recommendation**: Start with **Option C** (fixed η) for initial model fitting to verify the correction restores coherent results (shedding > fever after correction). Then upgrade to **Option A** (dose-dependent η) for the final model. Option B is the fallback if Option A is poorly identified.

**What this restores**:
- All 8 Oxford shedding observations (W-I-3, W-I-4, D-I-plac, J-I-ctrl, J-I-ViTT, J-I-ViPS, plus D-I-Ty21a/M01 for validation)
- The Jin Vi-TT and Vi-PS shedding contrasts, which are the **only data** that separately constrain γ_inf from γ_fev
- The biological prior γ_inf < γ_fev becomes testable
- Active observations increase from 19 to 25 (Tier 2 Option A) with +2 parameters, net improvement in data-to-parameter ratio

**Key insight**: The η-correction is not merely rescuing biased data — it is *modeling a real process* (time-to-shedding vs time-to-treatment) that we have partial information about. This is more principled than either dropping the data or using it uncorrected.

For the **infection** outcome, the definition difference is smaller. Maryland's "infection" (Hornick Table 2: low-grade fever, serology, blood culture, or shedding >5 days) is broadly comparable to Oxford's shedding definition. We apply no definition correction for infection, but note that Maryland's broader infection criteria may yield slightly higher rates than Oxford's shedding-only ascertainment.

---

## 3. Parameter Table

### 3.1 Biological Parameters (shared across all studies)

| Parameter | Symbol | Meaning | Units |
|-----------|--------|---------|-------|
| Infection N50 | $N_{50,\text{inf}}$ | Dose for 50% infection in naive, bicarb delivery | bicarb-equiv CFU |
| Fever-given-infection N50 | $N_{50,\text{fev\|inf}}$ | Dose for 50% fever conditional on infection, naive | bicarb-equiv CFU |
| Infection heterogeneity | $\alpha_{\text{inf}}$ | Shape/heterogeneity for infection | dimensionless |
| Fever-given-infection heterogeneity | $\alpha_{\text{fev\|inf}}$ | Shape/heterogeneity for fever given infection | dimensionless |
| Infection immunity scaling | $\gamma_{\text{inf}}$ | How steeply CoP reduces infection risk | dimensionless |
| Fever-given-infection immunity scaling | $\gamma_{\text{fev\|inf}}$ | How steeply CoP reduces fever risk given infected | dimensionless |

Note: The marginal P(fever) is not a beta-Poisson — it is $P_{\text{inf}} \times P_{\text{fev|inf}}$, the product of two beta-Poissons. The "N50 for fever" in the marginal sense is not a model parameter; it is a derived quantity that depends on all 6 biological parameters.

### 3.2 Nuisance Parameters

| Parameter | Symbol | Meaning | Constrained by |
|-----------|--------|---------|----------------|
| Medium offset | $\delta$ | Milk-to-bicarb dose conversion | ID50 position + CoP_md constraint |
| Maryland susceptible fraction | $\pi_{\text{susc}}$ | Fraction near-naive in Maryland | Gilman/Levine H-antibody (~0.5) |
| Maryland immune CoP | $\text{CoP}_{\text{imm}}$ | Immunity level of immune fraction | Curve shape asymmetry |
| Definition sensitivity | $\phi$ | Maryland/Oxford fever definition ratio | Oxford threshold analysis |

### 3.3 Design Choices

| Choice | Default | Alternative | Rationale |
|--------|---------|-------------|-----------|
| Shared vs specific α | Specific ($\alpha_{\text{inf}}$, $\alpha_{\text{fev}}$) | Shared with ratio | Data may not support 2 independent α |
| Shared vs specific γ | Specific ($\gamma_{\text{inf}}$, $\gamma_{\text{fev}}$) | Shared | Evidence for differential immunity effect on infection vs fever |
| Maryland mixture | 2-component | Single CoP_md | Gilman data supports mixture |
| Definition correction | Constant φ | Dose-dependent | Start simple, check residuals |

---

## 4. Data Catalog

### 4.1 Oxford Studies (Bicarbonate, Quailes strain, δ = 1)

All Oxford studies use the same strain, same delivery, same fever definition (TD = ≥38°C/12h OR bacteremia), and screen for typhoid naivety.

#### 4.1.1 Waddington 2014 — Multi-dose naive (INDEPENDENT)

Study ID: OVG2009/10. N=40 per-protocol (20 at 10³, 20 at 10⁴). **Two dose levels only** — the extraction erroneously included a 10⁵ arm that does not exist in this paper.

**Fever (TD composite) — PDF-verified 2026-03-18:**

| Obs ID | Dose (CFU) | n | y | p | Outcome |
|--------|-----------|---|---|---|---------|
| W-F-3 | 10³ | 20 | 11 | 0.55 | TD |
| W-F-4 | 10⁴ | 20 | 13 | 0.65 | TD |

**Infection (stool shedding) — EXCLUDED from primary likelihood (Section 2.6: treatment-truncation bias):**

| Obs ID | Dose (CFU) | n | y | p | Outcome |
|--------|-----------|---|---|---|---------|
| ~~W-I-3~~ | 10³ | 20 | 13 | 0.65 | Shedding. Excluded: treatment-truncation bias. |
| ~~W-I-4~~ | 10⁴ | 20 | ~8 | ~0.40 | Shedding. Excluded: treatment-truncation bias. Bounded 7-14 (not directly reported in paper). |

**Immunity status**: Screened naive. Baseline anti-Vi not reported per subject; Waddington notes "no statistically significant differences" in pre-challenge anti-Vi between diagnosed and not-diagnosed. Treat as CoP ≈ 1 with a note that ~29% may have had detectable anti-Vi (based on comparable Oxford cohort rates).

**Extraction errors corrected**: (1) 10⁵ arm removed — does not exist in paper. (2) 10⁴ denominator corrected from n=16 to n=20, fever from y=10 to y=13. Source of errant 10⁵ data under investigation (see `waddington_10e5_source_investigation.md`).

#### 4.1.2 Darton 2016 — Vaccine trial at ~1.8×10⁴ (INDEPENDENT)

Study ID: OVG2011/02. N=91 per-protocol. Dose: median 1.82×10⁴ CFU.

**Fever (TD composite):**

| Obs ID | Group | Dose (CFU) | n | y | p | CoP source |
|--------|-------|-----------|---|---|---|------------|
| D-F-plac | Placebo | 1.82e4 | 30 | 20 | 0.67 | Baseline anti-Vi measured; 37% >7.4 EU/mL (S1 data). Individual-level titers available. |
| D-F-M01 | M01ZH09 | 1.82e4 | 31 | 18 | 0.58 | VALIDATION ONLY. Did not raise anti-Vi; 19% had detectable baseline. |
| D-F-Ty21a | Ty21a | 1.82e4 | 30 | 13 | 0.43 | VALIDATION ONLY. Lacks Vi gene; 28% had detectable baseline anti-Vi. |

**Infection (shedding-only, from S1 individual-level data) — EXCLUDED from primary likelihood (Section 2.6: treatment-truncation bias):**

| Obs ID | Group | n | y | p | Outcome |
|--------|-------|---|---|---|---------|
| ~~D-I-plac~~ | Placebo | 30 | 19 | 0.63 | Stool shedding (S1 extraction). Excluded: treatment-truncation bias. |
| ~~D-I-M01~~ | M01ZH09 | 31 | 20 | 0.65 | Stool shedding (S1). Excluded. |
| ~~D-I-Ty21a~~ | Ty21a | 30 | 18 | 0.60 | Stool shedding (S1). Excluded. |

*Note: Oxford shedding rates (Waddington 65% at 10³, Darton 63%, Jin 71%) are internally consistent but systematically biased downward by early antibiotic treatment — shedding < diagnosis at all doses ≥10⁴. See Section 2.6.*

**Key immunity data**: HR = 0.29 per log₁₀ anti-Vi IgG for TD (p=0.006). This enters the likelihood as a regression constraint on γ_fev (see Section 5.3).

**Anti-Vi GMT at challenge (from Darton Figure/Supplementary data)**:
- Placebo: baseline levels (low, ~7-15 EU/mL for Vi+ subjects)
- Ty21a: modest rise
- M01ZH09: no significant rise

#### 4.1.3 Jin 2017 (VAST) — Vi vaccine trial at ~10⁴ (INDEPENDENT)

Study ID: OVG2014/08. N=103 per-protocol. Dose: 1-5×10⁴ CFU.

**Fever (TD composite):**

| Obs ID | Group | Dose | n | y | p | Post-vax anti-Vi GMT |
|--------|-------|------|---|---|---|---------------------|
| J-F-ctrl | Control | ~1e4 | 31 | 24 | 0.77 | Baseline only; 38% >7.4 |
| J-F-ViTT | Vi-TT | ~1e4 | 37 | 13 | 0.35 | 563 EU/mL |
| J-F-ViPS | Vi-PS | ~1e4 | 35 | 13 | 0.37 | 141 EU/mL |

**Infection (shedding) — EXCLUDED from primary likelihood (Section 2.6: treatment-truncation bias):**

| Obs ID | Group | n | y | p | Notes |
|--------|-------|---|---|---|-------|
| ~~J-I-ctrl~~ | Control | 31 | 22 | 0.71 | Excluded: shedding (71%) < diagnosis (77%). |
| ~~J-I-ViTT~~ | Vi-TT | 37 | 22 | 0.59 | Excluded. VE ratio (0.59/0.71 = 0.83) may still inform γ_inf as sensitivity check. |
| ~~J-I-ViPS~~ | Vi-PS | 35 | 21 | 0.60 | Excluded. VE ratio (0.60/0.71 = 0.85) available for sensitivity. |

**Additional fever thresholds (controls only)**: Used to inform φ.
- ≥37.5°C: 20/31 = 65%
- ≥38.0°C: 17/31 = 55%
- ≥38.5°C: 14/31 = 45%
- (Jin does NOT report ≥39.0°C; use Darton 2016 placebo for that threshold: 9/30 = 30%)

**Alternative diagnostic criteria (Supplementary Table S1, controls only)**: Provide additional validation targets.
- Fever ≥37.5°C + bacteraemia (any order): 18/31 = 58%
- Fever ≥38.0°C + bacteraemia (any order): 16/31 = 52%
- Fever ≥38.0°C with *subsequent* bacteraemia: 13/31 = 42% — best field-trial analog

**Pre-vaccination anti-Vi GMTs (Supplementary Table S2)**: Control 7.7, Vi-TT 5.6, Vi-PS 6.3 EU/mL (all NS). Confirms groups balanced at baseline. Control Day 28 GMT = 8.0 (no drift).

**Anti-Vi IgG by diagnosis status (Supplementary Table S2)**:
- Vi-TT: Diagnosed GMT 522, Undiagnosed 586 (p=0.84, NS). Within-group anti-Vi does NOT predict outcome → group-GMT approximation adequate.
- Vi-PS: Diagnosed GMT 73, Undiagnosed 207 (p=0.007). Within-group anti-Vi DOES predict outcome → group-GMT approximation introduces Jensen's inequality bias; integration over CoP distribution warranted.
- Control: Diagnosed 7.1, Undiagnosed 12.0 (p=0.19, NS).

**Anti-Vi IgG → TD**: OR = 0.37 (0.15-0.88), p<0.0001. Enters as regression constraint on γ_fev.

#### 4.1.4 Gibani 2020 — Rechallenge at ~2.5×10⁴ (PARTIALLY INDEPENDENT)

Study ID: OVG2014/01 (PATCH). Naive arms are independent; rechallenge arms overlap with earlier studies.

**Naive arms only (INDEPENDENT):**

| Obs ID | Group | Dose | n | y | p | Outcome |
|--------|-------|------|---|---|---|---------|
| G20-F-naive | ST Naive | ~2.5e4 | 19 | 12 | 0.63 | TD |

**Rechallenge arms (OVERLAP — use for γ estimation with caution):**

| Obs ID | Group | n | y | p | Notes |
|--------|-------|---|---|---|-------|
| G20-F-rech | ST-ST Rechallenge | 27 | 12 | 0.44 | Mixed prior vaccine + challenge history |

**Critical subgroup analysis** (informs susceptibility heterogeneity):
- No prior disease: 10/38 = 26% on rechallenge
- Prior disease: 25/37 = 68% on rechallenge

This enters as a constraint on the interpretation of γ (see Section 6).

### 4.2 Maryland Studies (Milk, Quailes strain, δ > 1)

All Maryland studies use the same strain, same delivery (milk), and Maryland fever definition (≥103°F for 24-36h = "disease").

#### 4.2.1 Hornick 1970 — Multi-dose controls (CANONICAL)

Source: Part 1, Table 1. This is the primary multi-dose dataset for the milk dose-response curve.

**Fever (Maryland definition):**

| Obs ID | Dose (CFU) | n | y | p | Notes |
|--------|-----------|---|---|---|-------|
| H-F-3 | 10³ | 14 | 0 | 0.00 | Upper 95% CI: 23% |
| H-F-5 | 10⁵ | 116 | 32 | 0.28 | Largest N; anchor point |
| H-F-7 | 10⁷ | 32 | 16 | 0.50 | ID50 in milk units |
| H-F-8 | 10⁸ | 9 | 8 | 0.89 | Small N |
| H-F-9 | 10⁹ | 42 | 40 | 0.95 | Near saturation |

**Immunity status**: Unknown. No baseline serology. Prison population, likely some prior Salmonella exposure. Modeled via the mixture (Section 2.4).

#### 4.2.2 Hornick 1970 — Infection vs disease at 10⁷ (UNIQUE)

Source: Part 1, Table 2 (Quailes strain only).

| Obs ID | Outcome | n | y | p |
|--------|---------|---|---|---|
| H-I-7 | Infected (any evidence) | 30 | 28 | 0.93 |
| H-FgI-7 | Disease given infected | 28 | 16 | 0.57 |

This is the **only data point** where P(fever|infected) is substantially below 1. It enters as a joint constraint on the infection and fever curves evaluated at the same effective dose.

#### 4.2.3 Hornick 1970 — Vaccine trial (PARTIALLY OVERLAPPING with 4.2.1)

Source: Part 2, Table 5. The "None" (unvaccinated) row likely overlaps substantially with Table 1. **Use vaccinated arms only as additional immunity data**:

| Obs ID | Vaccine | Dose | n | y | p |
|--------|---------|------|---|---|---|
| H-V-K5 | K (acetone) | 10⁵ | 43 | 4 | 0.09 |
| H-V-L5 | L (phenol) | 10⁵ | 45 | 3 | 0.07 |
| H-V-K7 | K | 10⁷ | 28 | 12 | 0.43 |
| H-V-L7 | L | 10⁷ | 24 | 13 | 0.54 |

**CoP for K/L vaccines**: Unknown in anti-Vi units. These are whole-cell killed vaccines with complex immune profiles. **Cannot be mapped to the anti-Vi CoP framework without additional assumptions.** Recommend: include as informative but with a separate, weakly-constrained CoP_KL parameter, or exclude from primary analysis and use for validation.

#### 4.2.4 Gilman 1977 — Ty21a trial controls at 10⁵ (INDEPENDENT of 4.2.1)

Trials conducted 1973-1976, after Hornick 1970 enrollment. Different volunteers.

**Fever (Maryland definition):**

| Obs ID | Group | n | y | p |
|--------|-------|---|---|---|
| Gil-F-ctrl | Combined controls | 64 | 31 | 0.48 |

**Fever by H-antibody stratum:**

| Obs ID | H Ab status | n | y | p |
|--------|-------------|---|---|---|
| Gil-F-Hlo | H < 1:20 | 36 | 22 | 0.61 |
| Gil-F-Hhi | H ≥ 1:20 | 17 | 4 | 0.24 |

**Shedding in controls (Trials 1&3 only):**

| Obs ID | n | y (late shedding 4-30d) | p |
|--------|---|-------------------------|---|
| Gil-I-ctrl | 43 | 26 | 0.60 |

These H-antibody-stratified data directly constrain $\pi_{\text{susc}}$ and provide information on the within-Maryland immunity heterogeneity.

#### 4.2.5 Levine 1976 — S27v trial controls at 10⁵ (LIKELY INDEPENDENT of 4.2.1)

Trials conducted 1970-1973. Some early trials may overlap with late Hornick enrollment, but later trials are independent.

**Fever (Maryland definition):**

| Obs ID | Trial | Year | n | y | p |
|--------|-------|------|---|---|---|
| Lev-F-1 | 1 | 1970 | 26 | 13 | 0.50 |
| Lev-F-2 | 2 | 1971 | 33 | 10 | 0.30 |
| Lev-F-3 | 3 | 1972 | 22 | 12 | 0.55 |
| Lev-F-4 | 4 | 1973 | 16 | 4 | 0.25 |

**H-antibody confirmation**: H ≥1:20 → 22% vs 48% (P=0.005). Consistent with Gilman.

**Note**: The trial-to-trial variability (25% to 55%) at the same dose suggests either batch effects, temporal shifts in cohort immunity, or random variation. This informs the study-level random effect (Section 5.5).

**Shedding in controls (Trial 1 only):**

| Obs ID | n | y (any stool+) | p |
|--------|---|----------------|---|
| Lev-I-1 | 26 | 19 | 0.73 |

### 4.3 Cohort Overlap Map

| Dataset | Time Period | Overlaps With | Resolution |
|---------|-----------|---------------|------------|
| Hornick 1970 Table 1 (multi-dose) | 1960-1970 | Hornick Table 5 unvax | Use Table 1 as canonical; Table 5 unvax excluded |
| Hornick 1970 Table 5 (vaccines) | 1960-1970 | Table 1 unvax rows | Use only vaccinated arms from Table 5 |
| Gilman 1977 | 1973-1976 | Minimal with Hornick | Treat as independent |
| Levine 1976 | 1970-1973 | Possibly Trial 1 with late Hornick | Flag Trial 1; Trials 2-4 likely independent |
| Waddington 2014 | 2011-2012 | None | Independent |
| Darton 2016 | 2011-2012 | None (different study) | Independent |
| Jin 2017 | 2015-2016 | None | Independent |
| Gibani 2020 naive | 2015-2017 | None (new naive recruits) | Independent |
| Gibani 2020 rechallenge | 2015-2017 | Waddington, Darton | Use with caution; model prior exposure explicitly |
| Gibani 2019 (pooled) | 2011-2017 | All Oxford | **DO NOT USE** — not independent; use primary papers |

### 4.4 Data Exclusions

| Excluded | Reason |
|----------|--------|
| ~~Hornick 1966 10³ data (9/14 = 64%)~~ | ~~Anomalous~~ — **RESOLVED**: The 9/14 was an extraction error. Both papers show 0/14. No exclusion needed. |
| Hornick 10⁸ data (8/9 = 89%) | Include but note n=9 is very small |
| Gibani 2019 pooled analysis | Reanalysis of primary study data; use primaries instead |
| Non-Quailes strain data (Hornick Table 2 non-Vi strains) | Different virulence; not comparable |
| Darton 2017 (PCR analysis) | Same cohort as Waddington; no new dose-response data |

---

## 5. Likelihood Specification

### 5.1 General Form

Each observation $j$ contributes a binomial likelihood:

$$
y_j \sim \text{Binomial}(n_j, p_j)
$$

The predicted probability $p_j$ depends on the data source. Below we specify $p_j$ for each category.

### 5.2 Oxford Naive Data (Waddington + Gibani naive)

For Oxford naive groups with δ = 1 and CoP ≈ 1 (screened typhoid-naive, low baseline anti-Vi):

$$
p_j = P_k(D_j, \text{CoP}=1; N_{50,k}, \alpha_k, \gamma_k)
$$

where $k$ is the outcome (infection or fever).

**Refinement**: If we account for the ~30% of Oxford controls with detectable anti-Vi:

$$
p_j = 0.7 \cdot P_k(D_j, 1) + 0.3 \cdot P_k(D_j, \overline{\text{CoP}}_{\text{Vi+}})
$$

This matters most for the Jin 2017 control arm (38% had baseline anti-Vi) and Darton 2016 placebo (37% had baseline anti-Vi per S1 data; published paper claims 40%). Individual-level anti-Vi titers are available from the Darton S1 Dataset for all 91 per-protocol subjects. For Waddington, the baseline anti-Vi was not significantly associated with outcome; treat as all-naive for the initial model.

### 5.3 Oxford Immunity Variation Data

For vaccinated groups at fixed dose $D$ with known post-vaccination CoP distribution:

$$
p_j = \int P_k(D, \text{CoP}) \cdot f_j(\text{CoP}) \, d\text{CoP}
$$

In practice, approximate with the group-level GMT:

$$
p_j \approx P_k(D, g(\text{GMT}_j))
$$

where $g(\cdot)$ is the CoP mapping from the titer model. This applies to:
- Jin 2017 Vi-TT (GMT = 563 EU/mL), Vi-PS (GMT = 141 EU/mL)

**Empirical evidence on the group-GMT approximation (Jin Supplementary Table S2)**: The adequacy of the group-GMT approximation differs by vaccine group. Vi-TT diagnosed vs undiagnosed GMT = 522 vs 586 (p=0.84): within-group anti-Vi variation does NOT predict outcome, so the group-GMT approximation is empirically justified. Vi-PS diagnosed vs undiagnosed GMT = 73 vs 207 (p=0.007): within-group variation DOES predict outcome, so the group-GMT approximation introduces Jensen's inequality bias for this group. For Vi-PS, integrating over the within-group CoP distribution is warranted. Individual-level anti-Vi data are not tabulated in Jin (only group GMTs and CI), but the Darton S1 individual-level data provides an empirical CoP distribution that can inform the integration for both studies.

**Darton vaccine arms are EXCLUDED from γ estimation via anti-Vi CoP.** Neither M01ZH09 nor Ty21a raised anti-Vi IgG (Darton p.13). Ty21a genetically lacks the Vi capsule gene. M01ZH09 is derived from Vi-expressing Ty2 (ΔaroC ΔssaV) but a single oral dose does not generate measurable anti-Vi IgG; its protection is primarily through anti-LPS. Both vaccines operate through non-Vi mechanisms and cannot be mapped to CoP via anti-Vi GMT. They are retained in Section 11 for **validation** (does the model correctly predict their attack rates are between placebo and Vi-vaccinated?) but do not enter the γ likelihood.

**The Darton anti-Vi regression — upgraded with individual-level data**: The published HR = 0.29 per log₁₀ anti-Vi can now be replaced by direct individual-level modeling. The S1 Dataset provides pre-challenge anti-Vi IgG for all 91 per-protocol subjects. This enables:
- Individual-level logistic regression: TD ~ anti-Vi IgG (replacing the aggregate HR)
- Proper integration over the within-group CoP distribution (addressing Jensen's inequality, per Reviewer 2 Minor Concern 4)
- Direct calibration of the CoP mapping g(anti-Vi) → CoP for use in the dose-response model

For the primary likelihood, the Darton placebo enters as a group-level **fever-only** observation (D-F-plac: 20/30) with the group's CoP distribution characterized by the individual-level anti-Vi titers from S1. D-I-plac (shedding) is excluded per Section 2.6. The individual-level data is available in `analysis_data/darton_individual_endpoints.csv`.

### 5.4 Maryland Multi-Dose Data (Hornick)

For Hornick multi-dose controls, the model probability accounts for (a) the medium offset δ, (b) the Maryland mixture, and (c) the definition correction φ:

**For infection** (Hornick Table 2 at 10⁷ only):

$$
p_j = \pi_{\text{susc}} \cdot P_{\text{inf}}(D_j/\delta, \text{CoP}_{\text{susc}}) + (1-\pi_{\text{susc}}) \cdot P_{\text{inf}}(D_j/\delta, \text{CoP}_{\text{imm}})
$$

**For fever** (Hornick Table 1 at all doses):

$$
p_j = \phi \cdot \left[\pi_{\text{susc}} \cdot P_{\text{fev}}(D_j/\delta, \text{CoP}_{\text{susc}}) + (1-\pi_{\text{susc}}) \cdot P_{\text{fev}}(D_j/\delta, \text{CoP}_{\text{imm}})\right]
$$

The $\phi$ factor converts from the model's Oxford-calibrated fever probability to the Maryland-observed disease rate.

### 5.5 Maryland Single-Dose Studies (Gilman, Levine)

These provide additional 10⁵ control data with the same likelihood structure as Section 5.4. The trial-to-trial variability in Levine controls (25%-55% at the same dose) motivates a **study-level random effect**:

$$
p_{j,s} = \phi \cdot \left[\pi_{\text{susc}} \cdot P_{\text{fev}}(D_j/\delta, \text{CoP}_{\text{susc}} \cdot e^{\epsilon_s}) + (1-\pi_{\text{susc}}) \cdot P_{\text{fev}}(D_j/\delta, \text{CoP}_{\text{imm}} \cdot e^{\epsilon_s})\right]
$$

where $\epsilon_s \sim \text{Normal}(0, \sigma_{\text{study}})$ is a study-level random effect that perturbs the effective immunity level.

**Revised per Reviewer 2**: A logit-scale random effect is more standard and better-behaved near 0 and 1: $\text{logit}(p_{j,s}) = \text{logit}(p_{\text{model},j}) + \epsilon_s$ where $\epsilon_s \sim \text{Normal}(0, \sigma_{\text{study}})$. This absorbs batch effects, temporal cohort shifts, and other unmodeled variability.

**Alternative**: Instead of a random effect, use a beta-binomial likelihood (overdispersed binomial) for Maryland studies:

$$
y_j \sim \text{BetaBinomial}(n_j, p_j \cdot \kappa, (1-p_j) \cdot \kappa)
$$

where $\kappa$ is a concentration parameter. This is simpler and achieves the same goal of accommodating extra-binomial variation.

**Recommendation**: Start with beta-binomial overdispersion for Maryland studies; switch to random effects if you need the CoP interpretation.

### 5.6 Maryland Infection-Disease Split (Hornick Table 2)

At 10⁷ Quailes in milk, n=30:
- 28 infected, 16 with disease, 2 not infected

This enters as two binomial terms:

$$
y_{\text{inf}} = 28 \sim \text{Binomial}(30, p_{\text{inf}})
$$
$$
y_{\text{dis|inf}} = 16 \sim \text{Binomial}(28, p_{\text{dis|inf}})
$$

where:

$$
p_{\text{inf}} = \pi_{\text{susc}} \cdot P_{\text{inf}}(10^7/\delta, 1) + (1-\pi_{\text{susc}}) \cdot P_{\text{inf}}(10^7/\delta, \text{CoP}_{\text{imm}})
$$

$$
p_{\text{dis|inf}} = \frac{\phi \cdot [\pi_{\text{susc}} \cdot P_{\text{fev}}(10^7/\delta, 1) + (1-\pi_{\text{susc}}) \cdot P_{\text{fev}}(10^7/\delta, \text{CoP}_{\text{imm}})]}{p_{\text{inf}}}
$$

This joint likelihood term connects the infection and fever curves through the conditional, providing the strongest constraint on their relative positions.

### 5.7 Gilman H-Antibody Stratified Data

The H-antibody-stratified attack rates enter as separate binomial observations for each stratum:

$$
y_{\text{H-lo}} \sim \text{Binomial}(36, \phi \cdot P_{\text{fev}}(10^5/\delta, \text{CoP}_{\text{susc}}))
$$
$$
y_{\text{H-hi}} \sim \text{Binomial}(17, \phi \cdot P_{\text{fev}}(10^5/\delta, \text{CoP}_{\text{imm}}))
$$

This directly constrains the contrast between the mixture components, providing the most direct information about $\text{CoP}_{\text{imm}}$ relative to $\text{CoP}_{\text{susc}}$.

---

## 6. Exchangeability and Portability Assumptions

### 6.1 What Is Shared Across All Studies (50-year span)

| Shared quantity | Justification | Risk |
|----------------|---------------|------|
| N50_inf, N50_fev | Same strain (Quailes), same species, same biological process | Medium: host population genetics may differ |
| α_inf, α_fev | Heterogeneity in host biology; should be similar in healthy adults | Medium: protocol differences may inflate apparent α |
| γ_inf, γ_fev | How immunity scales dose-response; fundamental host-pathogen interaction | Low-medium: assumes same immune mechanisms operate |
| β-Poisson functional form | Mechanistically derived from independent single-hit infection | Low: well-established for enteric pathogens |

**Core assumption**: The biological dose-response curve (in bicarb-equivalent units, for a naive host) is the same in 1965 as in 2015. What changes is the delivery medium, the outcome definition, and the host population's immunity distribution.

### 6.2 What Differs Between Eras

| Differing quantity | Parameterized as | Justification |
|-------------------|-----------------|---------------|
| Delivery medium | δ_medium | Bicarbonate vs milk: different gastric acid neutralization |
| Background immunity | CoP distribution | Prison population (1960s) vs screened volunteers (2010s) |
| Outcome definition | φ (definition sensitivity) | 103°F/36h vs 38°C/12h OR bacteremia |
| Ascertainment methods | Not modeled (absorbed) | Manual vs automated blood culture; culture vs PCR |

### 6.3 What Differs Between Studies Within an Era

**Within Maryland**:
- Batch-to-batch variation in challenge preparation (Levine controls: 25-55% at same dose)
- Temporal shifts in prison population immunity
- Modeled as: study-level random effect or beta-binomial overdispersion

**Within Oxford**:
- Dose varies (10³-10⁵ in Waddington; ~10⁴ in others)
- Baseline anti-Vi fraction varies (19-40% across studies)
- Modeled as: explicit CoP variation per group

### 6.4 Explicit Assumptions About the 50-Year Gap

These are the load-bearing assumptions of the cross-era bridge. Each should be testable:

**Assumption 1: Strain virulence is constant.** The Quailes strain was propagated and stored frozen. Laboratory passage could alter virulence, but Oxford explicitly uses the same strain lineage. *Test*: Compare Oxford and Maryland attack rates at equivalent effective doses (after adjusting for δ and CoP). If the adjusted rates are consistent, the assumption holds.

**Assumption 2: The β-Poisson shape parameter α is portable.** This assumes that the distribution of host susceptibility factors (gastric pH variation, mucosal immunity, genetic susceptibility) is similar in healthy adult males in 1960s Maryland prisons and healthy adults in 2010s Oxford. *Risk*: Diet, microbiome, and genetic composition may differ. *Mitigation*: α is estimated jointly from both datasets; if it's truly different, the posterior will show tension (detectable via leave-one-era-out cross-validation).

**Assumption 3: The immunity scaling γ is portable.** This assumes the mechanism by which immunity reduces susceptibility is the same across eras. Since γ is estimated primarily from Oxford data (where immunity is measured), this assumption is mainly invoked when applying the Oxford-derived γ to predict the Maryland immunity effect. *Test*: The Maryland vaccine data (K, L vaccines) should show consistent VE predictions if γ is correct.

**Assumption 4: The medium offset δ is a pure biological parameter.** Bicarbonate neutralizes gastric acid → more bacteria survive → higher effective dose. This is a biochemical fact that should be the same in any human gut. *Risk*: The effective offset may depend on stomach contents, fasting protocol, milk fat content — factors that differ across protocols. *Mitigation*: Treat δ as uncertain with a wide prior spanning 2-5 logs.

**Assumption 5: The CoP mapping (anti-Vi IgG → protection) is applicable to the Maryland era.** Anti-Vi was measured in Oxford but not in Maryland. We assume that whatever anti-Vi levels the Maryland subjects had would predict their susceptibility through the same mechanism observed in Oxford. *Risk*: Maryland studies found NO correlation between O/H/Vi agglutinins and protection. This may reflect (a) different assays, (b) different immune mechanisms, or (c) genuinely different biology. *Mitigation*: The Maryland CoP is latent and estimated from the data (via the mixture model), not from assumed anti-Vi levels. The assumption is softer than it appears: we only assume that SOME scalar immunity measure exists and that the β-Poisson with CoP^γ scaling is the right functional form.

---

## 7. Prior Specification

### 7.1 Biological Parameters

| Parameter | Prior | Justification |
|-----------|-------|---------------|
| $\log_{10}(N_{50,\text{inf}})$ | Normal(2.5, 1.0) | Expect ~300 bicarb CFU; constrained by Maryland shedding through δ bridge (Oxford shedding excluded per Section 2.6) |
| $\log_{10}(N_{50,\text{fev}})$ | Normal(2.8, 1.0) | Expect ~600 bicarb CFU; slightly higher than infection |
| $\alpha_{\text{inf}}$ | LogNormal(-1.5, 0.8) | Expect ~0.1-0.5; must be positive; Oxford data suggests small α |
| $\alpha_{\text{fev}}$ | LogNormal(-1.5, 0.8) | Same reasoning as infection |
| $\gamma_{\text{inf}}$ | LogNormal(-0.5, 0.7) | Expect ~0.2-1.0; positive; immunity reduces infection somewhat |
| $\gamma_{\text{fev}}$ | LogNormal(0.0, 0.7) | Expect ~0.3-2.0; positive; immunity reduces fever more than infection |

**Ordering constraint**: Apply a soft prior penalty if $N_{50,\text{inf}} > N_{50,\text{fev}}$ (infection threshold should be lower than fever threshold). Specifically: $\log_{10}(N_{50,\text{fev}}) - \log_{10}(N_{50,\text{inf}}) \sim \text{HalfNormal}(0, 1)$ (positive, allowing 0-2 log difference).

### 7.2 Nuisance Parameters

| Parameter | Prior | Justification |
|-----------|-------|---------------|
| $\log_{10}(\delta)$ | Normal(3.5, 0.7) | Expect ~1000-30,000× (3-4.5 logs); biochemistry of gastric acid neutralization |
| $\pi_{\text{susc}}$ | Beta(7, 4) | Centered at ~0.65; Gilman shows 36/53 (68%) with H antibody <1:20 |
| $\text{CoP}_{\text{imm}}$ | LogNormal(1.0, 0.5) | Expect CoP ~2-5 for partially immune; must be >1 |
| $\phi$ | Beta(5, 5) | Centered at 0.5; Oxford threshold analysis suggests ~0.3-0.7 |
| $\sigma_{\text{study}}$ | HalfNormal(0, 0.5) | Study-level variability in effective immunity; modest |
| $\text{CoP}_{\text{susc}}$ | LogNormal(0, 0.2) | Near 1 (near-naive); allows slight deviation |

### 7.3 Prior Predictive Check

Before fitting, simulate from the prior to verify:
1. The prior-predicted Oxford naive attack rates span the observed range (~40-80% at 10³-10⁴)
2. The prior-predicted Maryland curve spans 0-100% over 10³-10⁹
3. The infection curve is above the fever curve at all doses
4. The medium offset δ is in the 2-5 log range
5. The Maryland mixture produces plausible asymmetric S-curves

---

## 8. Identifiability Analysis

### 8.1 Well-Identified Parameters

| Parameter | Primary data source | Why well-identified |
|-----------|-------------------|---------------------|
| $N_{50,\text{fev&#124;inf}}$ | Waddington multi-dose + pooled 10⁴ data | Near N50 at 10³ bicarb; 5+ data groups |
| $\phi$ | Jin + Darton fever threshold analyses | Direct measurement of definition gap (Jin: ≥38.0, ≥38.5; Darton: ≥39.0) |
| $\pi_{\text{susc}}$ | Gilman + Levine H-antibody stratification | Direct measurement of immunity prevalence |

### 8.2 Moderately Identified Parameters

| Parameter | Primary data source | Why moderate (not strong) |
|-----------|-------------------|--------------------------|
| $\gamma_{\text{fev&#124;inf}}$ | Jin Vi-TT/Vi-PS vaccine VE + Darton baseline anti-Vi HR | **Downgraded from "well-identified"**: Darton M01ZH09 and Ty21a excluded (non-Vi mechanisms). Remaining constraints: Jin Vi-TT (35% vs 77%) and Vi-PS (37% vs 77%) = 2 Vi-vaccine contrasts at known GMT, plus Darton baseline anti-Vi HR (continuous, across all arms including placebo). Moderate but thinner than previously claimed. |
| $\gamma_{\text{inf}}$ | Jin Vi-TT/Vi-PS shedding + Gibani 2019 shedding ORs | Jin Vi-vaccine arms show modest shedding reduction (59-60% vs 71%); Gibani 2019 pooled anti-Vi → shedding P<.0001 |

### 8.3 Weakly Identified Parameters

| Parameter | Why weak | What helps |
|-----------|----------|------------|
| $\alpha_{\text{inf}}$, $\alpha_{\text{fev&#124;inf}}$ | Oxford has only ~1 log of informative dose range | Maryland 6-log range through δ bridge |
| $N_{50,\text{inf}}$ | Shedding rates ~63-65% at 10³-10⁴ (flat); N50 well below observed range | Bounded above by 10³; lower bound from model coherence |
| $\text{CoP}_{\text{imm}}$ | Confounded with δ through the ID50 position | Gilman stratified data + curve shape asymmetry |

### 8.3 Structural Confounds

**Confound 1: δ × CoP_imm on the Maryland ID50 position.**
The product $\delta \cdot f(\text{CoP}_{\text{imm}})$ is constrained by where the Maryland curve crosses 50%, but the individual factors are not. Breaking this requires either: (a) external information on δ (biochemistry), (b) external information on CoP_imm (serology we don't have), or (c) the curve shape constraint from Maryland's 6-log range (partial).

**Resolution**: The prior on δ (informed by gastric acid biochemistry: ~3-4.5 logs) and the prior on CoP_imm (modest: ~2-5) jointly constrain the ridge. The posterior will show a negative correlation between δ and CoP_imm, but the marginals should be interpretable.

**Confound 2: α versus CoP_imm in the Maryland curve shape.**
A flatter α (more heterogeneity) and a lower CoP_imm can both explain a gentle lower tail. Breaking this requires the steepness at HIGH doses (where even immune subjects get infected), which constrains α_eff = α/CoP^γ. With γ estimated from Oxford, and the high-dose Maryland data (89% at 10⁸, 95% at 10⁹), α and CoP_imm are partially separable.

**Confound 3: φ × δ × CoP_imm three-way degeneracy (added per Reviewer 2).**
The Maryland fever likelihood at any dose is $p = \phi \cdot [\text{mixture probability}]$, where the mixture probability depends on $\delta$ and $\text{CoP}_{\text{imm}}$. Any decrease in $\phi$ can be compensated by an increase in the mixture probability (via lower $\delta$ or lower $\text{CoP}_{\text{imm}}$). The three-way surface should be mapped explicitly via the $\phi_{\min} \times \delta$ grid analysis (Section 9.3). The dose-dependent $\phi$ and the Oxford threshold calibration partially constrain this surface, but the residual degeneracy should be reported transparently.

**Confound 4: α × N50 correlation from Oxford.**
Standard β-Poisson identifiability issue. With only ~1 log of dose range from Oxford, N50 and α are correlated. The Maryland data (with much wider dose range, after the δ bridge) breaks this correlation substantially.

**Confound 5: γ identification after Darton vaccine exclusion.**
With Darton M01ZH09 and Ty21a excluded (non-Vi mechanisms), the γ parameters rest on: (a) Jin Vi-TT and Vi-PS contrasts (2 groups with known anti-Vi GMT vs 1 control group, at a single dose), and (b) the Darton baseline anti-Vi HR (continuous variation within the full Darton cohort, but confounded with other baseline differences). This is thinner than previously assessed. If γ is weakly identified from Oxford alone, the Maryland data — which should be constraining δ and CoP_imm — will also absorb γ uncertainty, creating a four-way interaction (δ, CoP_imm, γ, φ). The staged fitting diagnostic (Section 10.4) should explicitly report the γ marginal posterior width at Stage 1 to assess whether this is a problem.

**Mitigating factor**: The Gibani 2020 rechallenge data (naive 63% vs rechallenge 44%) provides an additional γ constraint through natural immunity, independent of the Vi pathway. However, this is confounded by the susceptibility paradox (Section 9.2) and by mixed prior vaccination in the rechallenge arm.

**Diagnostic**: Compute and report the expected posterior correlation matrix for all parameters from a Laplace approximation at a plausible parameter point. Verify that the pairwise confounds discussed above are the dominant structure and that no unexpected higher-order interactions exist.

### 8.4 Information Flow Summary

```
LAYER 1: Oxford naive (Waddington — 2 doses only: 10³ and 10⁴)
  → N50_inf, N50_fev|inf [bounded above by 10³; lower bound very weak]
  → α_inf, α_fev|inf [essentially unconstrained — 2 points in the flat tail]
  NOTE: Oxford alone is underconstrained for N50 and α. This is expected
  and is why the cross-era bridge exists. Stage 1 fit will show wide
  posteriors; Stage 2 (adding Maryland) is where α gets pinned.

LAYER 2: Oxford immunity variation (Jin Vi-TT/Vi-PS only for γ calibration)
  → γ_inf, γ_fev|inf [moderate — only 2 Vi vaccine contrasts + baseline anti-Vi HR]
  → Darton M01ZH09/Ty21a EXCLUDED from γ (non-Vi mechanisms); used for validation
  → Validates CoP mapping from titer model

LAYER 3: Maryland multi-dose (Hornick) via δ bridge
  → Sharpens α_inf, α_fev [strong — 6 log dose range]
  → Constrains δ × f(CoP_imm) [position]
  → Partially separates δ from CoP_imm [shape]

LAYER 4: Maryland immunity heterogeneity (Gilman, Levine)
  → π_susc, CoP_imm [direct measurement via H-antibody]
  → Further separates δ from CoP_imm

LAYER 5: Maryland infection-fever split (Hornick Table 2)
  → P(fever|infected) = 57% at 10⁷ milk
  → Constrains relative position of infection and fever curves
  → Uniquely identifies the infection-fever gap at sub-saturated levels

LAYER 6: Cross-checks and overdispersion
  → Levine trial-to-trial variability → σ_study
  → φ from Oxford threshold analysis → definition correction
  → Maryland vaccine data → qualitative γ cross-check
```

---

## 9. Model Checking Plan

### 9.1 Posterior Predictive Checks

For each observation in the data:
1. Simulate $y_j^{\text{rep}} \sim \text{Binomial}(n_j, p_j^{\text{rep}})$ from posterior draws
2. Compare distribution of $y_j^{\text{rep}}$ to observed $y_j$
3. Flag observations where observed value is in the extreme tails (<2.5% or >97.5%)

### 9.2 Key Predictions to Verify

| Prediction | Expected | Why it matters |
|------------|----------|---------------|
| Oxford attack rates at 10³-10⁴ | 50-70% fever, 60-75% shedding | Core fit to best data |
| Maryland curve shape (0% → 95% over 10³-10⁹) | Monotonic S-curve | Tests δ + mixture model |
| Gilman H-antibody contrast (24% vs 61%) | Captured by mixture | Tests CoP_imm identification |
| Hornick Table 2 conditional (57%) | Consistent with infection/fever curve ratio | Tests two-curve coherence |
| Oxford vaccine VE (~55%) | From predicted CoP shift | Tests γ extrapolation |
| Jin Vi-TT VE at fever+subsequent bacteraemia (~87%) | Predicted from differential γ_inf vs γ_fev\|inf | Tests whether cascade model captures stronger VE at stricter endpoints |
| Jin control fever+subsequent bacteraemia (42%) | Predicted from P_inf × P_fev\|inf at ~10⁴ | Tests cascade model internal consistency (Table S1 data not in likelihood) |
| Levine trial-to-trial variability | Within overdispersion bounds | Tests σ_study |
| **Gibani rechallenge subgroups (HOLD-OUT)** | **Model should predict prior-disease → higher rechallenge risk** | **Critical test of γ interpretation** |

**Gibani paradox hold-out (added per Reviewer 2)**: Fit the model excluding the Gibani 2020 rechallenge arms. Then predict the split: prior disease (25/37 = 68%) vs no prior disease (10/38 = 26%). Under the adaptive immunity model (γ > 0), the model should predict that subjects with prior disease are MORE protected (lower attack rate). The data show the opposite. If the model predicts the wrong direction, this is strong evidence for innate susceptibility classes and a fundamental limitation of the single-γ framework. This test should be run at Stage 1 (Oxford only) and again at Stage 3 (full model).

### 9.3 Sensitivity Analyses

| What to vary | How | Purpose |
|-------------|-----|---------|
| Exclude Maryland entirely | Oxford-only fit | Assess Maryland's marginal contribution |
| Exclude Oxford entirely | Maryland-only fit (if identifiable) | Assess Oxford's marginal contribution |
| Fix δ at 3, 3.5, 4, 4.5 logs | Grid of fixed values | Map the δ-CoP_imm ridge |
| Share α across outcomes | α_inf = α_fev | Test whether data supports separate α |
| Share γ across outcomes | γ_inf = γ_fev | Test whether data supports separate γ |
| Drop definition correction (φ = 1) | Treat definitions as equivalent | Assess φ's impact on δ and CoP |
| Replace mixture with single CoP_md | Simpler Maryland model | Test whether heterogeneity matters |
| Hornick 10³ influence check | Fit with 0/14 vs excluded | High-influence zero; 1966-1970 discrepancy resolved (both show 0/14, extraction was wrong) |
| Add susceptibility classes | π_innate_resistant fraction | Address Gibani paradox |

### 9.4 Leave-One-Study-Out Cross-Validation

For each study, refit the model excluding that study and predict its outcomes. Particularly informative:
- Leave out Jin 2017 → predict 77% control attack rate
- Leave out Gilman 1977 → predict H-antibody stratification
- Leave out Hornick 10⁸ → predict 89% attack rate
- Leave out Hornick Table 2 → predict infection-fever split

---

## 10. Implementation Notes

### 10.1 Stan Parameterization

Work in **log space** for positive parameters: $\log_{10}(N_{50})$, $\log(\alpha)$, $\log(\gamma)$, $\log_{10}(\delta)$, $\log(\text{CoP}_{\text{imm}})$. This ensures positivity and improves sampler geometry.

For the β-Poisson, precompute the constant $c_k = (2^{1/\alpha_k} - 1) / N_{50,k}$ and evaluate:

```
p = 1 - (1 + D_eff * c_k) ^ (-alpha_k / CoP ^ gamma_k)
```

where $D_{\text{eff}} = D / \delta$ for milk studies and $D_{\text{eff}} = D$ for bicarb studies.

### 10.2 Numerical Stability

The expression $(1 + x)^{-a}$ for very large $x$ and very small $a$ (the Maryland regime with immunity) may lose floating-point precision. Use `log1p` and `exp`:

```
log_survival = -alpha_eff * log1p(D_eff * c_k)
p = 1 - exp(log_survival)
```

For $p$ very close to 0 or 1, use `log1m_exp` for the log-likelihood to avoid catastrophic cancellation.

### 10.3 Expected Sampling Challenges

1. **δ-CoP_imm ridge**: The posterior will have a narrow ridge in (log δ, log CoP_imm) space. Use informative priors to regularize. Consider reparameterizing to (sum, difference) coordinates.

2. **α-N50 correlation**: Standard for β-Poisson models. The Maryland data should help, but check the bivariate posterior.

3. **Mixture model multimodality**: The two-component Maryland mixture could have label-switching issues. Enforce ordering: CoP_susc < CoP_imm.

4. **Small-n data**: Hornick 10⁸ (n=9) contributes weak likelihood but can create posterior irregularities. The beta-binomial overdispersion helps.

### 10.4 Staged Fitting Strategy

To diagnose issues, fit in stages:

**Stage 1**: Oxford only (infection + fever). Fix δ = 1, CoP_susc = 1. Estimate (N50_inf, N50_fev|inf, α_inf, α_fev|inf, γ_inf, γ_fev|inf). **Expected behavior**: N50 bounded above (~10³) but with very wide lower bound. α essentially unconstrained (only 2 dose points, 1 log apart, both in the flat upper tail). γ moderately constrained from Jin Vi-TT/Vi-PS contrasts + Darton baseline anti-Vi variation. This stage confirms that the β-Poisson machinery works and that γ is estimable, but it will NOT produce useful N50 or α estimates. That is by design — α requires the 6-log Maryland dose range.

**Stage 2**: Add Maryland multi-dose fever data. Add δ, π_susc, CoP_imm, φ. **This is where α gets pinned.** The Maryland curve spanning 10³-10⁹ (0% to 95%) constrains the dose-response shape that Oxford alone cannot. Verify that (a) the cross-era bridge produces a coherent joint fit, (b) the α posterior narrows dramatically relative to Stage 1, and (c) the δ-CoP_imm ridge is present but manageable.

**Stage 3**: Add Maryland infection data (Hornick Table 2) and single-dose studies (Gilman, Levine). Add σ_study. Full model.

**Stage 4**: Sensitivity analyses from Section 9.3.

---

## 11. Summary of All Observations Entering the Likelihood

For quick reference, every binomial observation:

### Oxford (δ = 1)

| ID | Study | Outcome | Dose | n | y | CoP status |
|----|-------|---------|------|---|---|------------|
| W-F-3 | Waddington | Fever | 10³ | 20 | 11 | Naive |
| W-F-4 | Waddington | Fever | 10⁴ | 20 | 13 | Naive. **CORRECTED**: n=20 not 16, y=13 not 10 (PDF-verified 2026-03-18). |
| ~~W-F-5~~ | | | | | | **REMOVED**: 10⁵ arm does not exist in Waddington 2014. Extraction hallucination. |
| ~~W-I-3~~ | ~~Waddington~~ | ~~Infection (shedding)~~ | ~~10³~~ | ~~20~~ | ~~13~~ | **EXCLUDED (Section 2.6)**: Treatment-truncation bias. Oxford shedding underestimates true infection. |
| ~~W-I-4~~ | ~~Waddington~~ | ~~Infection (shedding)~~ | ~~10⁴~~ | ~~20~~ | ~~7-14~~ | **EXCLUDED (Section 2.6)**: Same bias + not directly reported. Bounded 7-14. |
| ~~W-I-5~~ | | | | | | **REMOVED**: 10⁵ arm does not exist. |
| D-F-plac | Darton | Fever | 1.82e4 | 30 | 20 | Mixed (40% Vi+) |
| D-F-Ty21a | Darton | Fever | 1.82e4 | 30 | 13 | Ty21a. **VALIDATION ONLY** — lacks Vi gene; non-Vi protection mechanism. Cannot map to anti-Vi CoP. |
| D-F-M01 | Darton | Fever | 1.82e4 | 31 | 18 | M01ZH09. **VALIDATION ONLY** — did not raise anti-Vi IgG; protection via anti-LPS. Cannot map to anti-Vi CoP. |
| ~~D-I-plac~~ | ~~Darton~~ | ~~Infection (shedding)~~ | ~~1.82e4~~ | ~~30~~ | ~~19~~ | **EXCLUDED (Section 2.6)**: Treatment-truncation bias. Shedding 63% < diagnosis 67%. |
| ~~D-I-Ty21a~~ | ~~Darton~~ | ~~Infection (shedding)~~ | ~~1.82e4~~ | ~~30~~ | ~~18~~ | **EXCLUDED (Section 2.6)**. |
| ~~D-I-M01~~ | ~~Darton~~ | ~~Infection (shedding)~~ | ~~1.82e4~~ | ~~31~~ | ~~20~~ | **EXCLUDED (Section 2.6)**. |
| J-F-ctrl | Jin | Fever | ~2e4 | 31 | 24 | Mixed (38% Vi+). Dose is target range 1-5×10⁴; no median reported; ~2e4 is geometric mean of range. |
| J-F-ViTT | Jin | Fever | ~2e4 | 37 | 13 | Vi-TT (GMT 563) |
| J-F-ViPS | Jin | Fever | ~2e4 | 35 | 13 | Vi-PS (GMT 141) |
| ~~J-I-ctrl~~ | ~~Jin~~ | ~~Infection (shedding)~~ | ~~~2e4~~ | ~~31~~ | ~~22~~ | **EXCLUDED (Section 2.6)**: Shedding 71% < diagnosis 77%. VE ratios available for γ_inf sensitivity. |
| ~~J-I-ViTT~~ | ~~Jin~~ | ~~Infection (shedding)~~ | ~~~2e4~~ | ~~37~~ | ~~22~~ | **EXCLUDED (Section 2.6)**. Shedding VE ratio = 0.59/0.71 = 0.83 for sensitivity. |
| ~~J-I-ViPS~~ | ~~Jin~~ | ~~Infection (shedding)~~ | ~~~2e4~~ | ~~35~~ | ~~21~~ | **EXCLUDED (Section 2.6)**. Shedding VE ratio = 0.60/0.71 = 0.85 for sensitivity. |
| G20-F-naive | Gibani | Fever | ~2.5e4 | 19 | 12 | Naive |

### Maryland (δ > 1)

| ID | Study | Outcome | Dose | n | y | CoP status |
|----|-------|---------|------|---|---|------------|
| H-F-3 | Hornick | Fever | 10³ | 14 | 0 | Maryland mixture |
| H-F-5 | Hornick | Fever | 10⁵ | 116 | 32 | Maryland mixture |
| ~~H-F-7~~ | ~~Hornick~~ | ~~Fever~~ | ~~10⁷~~ | ~~32~~ | ~~16~~ | **REMOVED: double-counts with H-I-7 + H-FgI-7** |
| H-F-8 | Hornick | Fever | 10⁸ | 9 | 8 | Maryland mixture |
| H-F-9 | Hornick | Fever | 10⁹ | 42 | 40 | Maryland mixture |
| *H-V-K5* | *Hornick* | *Fever* | *10⁵* | *43* | *4* | *K vaccine. CoP unknown (whole-cell killed). Excluded from primary; available for validation.* |
| *H-V-L5* | *Hornick* | *Fever* | *10⁵* | *45* | *3* | *L vaccine. CoP unknown. Excluded from primary; available for validation.* |
| *H-V-K7* | *Hornick* | *Fever* | *10⁷* | *28* | *12* | *K vaccine at 10⁷. Excluded from primary.* |
| *H-V-L7* | *Hornick* | *Fever* | *10⁷* | *24* | *13* | *L vaccine at 10⁷. Excluded from primary.* |
| H-I-7 | Hornick | Infection | 10⁷ | 30 | 28 | Maryland mixture |
| H-FgI-7 | Hornick | Fever\|Inf | 10⁷ | 28 | 16 | Maryland mixture |
| ~~Gil-F-ctrl~~ | ~~Gilman~~ | ~~Fever~~ | ~~10⁵~~ | ~~64~~ | ~~31~~ | **REPLACED by stratified observations below** |
| Gil-F-Hlo | Gilman | Fever | 10⁵ | 36 | 22 | Susceptible stratum (H Ab <1:20). Exact counts from Table 5. |
| Gil-F-Hhi | Gilman | Fever | 10⁵ | 17 | 4 | Immune stratum (H Ab ≥1:20). Exact counts from Table 5. |
| Gil-F-rest | Gilman | Fever | 10⁵ | 11 | 5 | **DERIVED** by subtraction (64-36-17=11, 31-22-4=5). Controls without baseline H-antibody data. Sensitivity: fit with and without this row. |
| Gil-I-ctrl | Gilman | Infection | 10⁵ | 43 | 26 | Maryland mixture. **NOTE**: endpoint is late shedding (4-30 days post-challenge), not any-time shedding. Trials 1&3 controls only. |
| Lev-F-1 | Levine | Fever | 10⁵ | 26 | 13 | Maryland mixture. Levine def = ≥101°F + culture → φ_L ≈ 0.65 (Section 2.5). |
| Lev-F-2 | Levine | Fever | 10⁵ | 33 | 10 | Maryland mixture. Same φ_L ≈ 0.65. |
| Lev-F-3 | Levine | Fever | 10⁵ | 22 | 12 | Maryland mixture. Same φ_L ≈ 0.65. |
| Lev-F-4 | Levine | Fever | 10⁵ | 16 | 4 | Maryland mixture. Same φ_L ≈ 0.65. |
| Lev-I-1 | Levine | Infection | 10⁵ | 26 | 19 | Maryland mixture. **NOTE**: endpoint is any-time stool positive (differs from Gilman's 4-30 day definition). |
| Lev-I-2 | Levine | Infection | 10⁵ | 33 | 15 | Maryland mixture. Trial 2, 1971. Any-time stool positive. Same definition caveat as Lev-I-1. |
| Lev-I-3 | Levine | Infection | 10⁵ | 22 | 17 | Maryland mixture. Trial 3, 1972. Any-time stool positive. |
| Lev-I-4 | Levine | Infection | 10⁵ | 16 | 6 | Maryland mixture. Trial 4, 1973. Any-time stool positive. |

**Active observations (Tier 1 / primary likelihood)**: 7 Oxford fever (W-F-3, W-F-4, D-F-plac, J-F-ctrl, J-F-ViTT, J-F-ViPS, G20-F-naive) + 18 Maryland (4 Hornick fever + 2 Hornick infection/conditional + 3 Gilman fever strata + 1 Gilman infection + 4 Levine fever + 4 Levine infection) = **25 active calibration observations**. Struck rows (H-F-7, W-F-5, W-I-5, Gil-F-ctrl) and validation-only rows (D-F-Ty21a, D-F-M01, H-V-*) are already excluded from this count. Plus 6 validation observations.
**Active observations (Tier 2 / η-correction)**: Restores 6 Oxford shedding rows (W-I-3, W-I-4, D-I-plac, J-I-ctrl, J-I-ViTT, J-I-ViPS) with η correction = **31 active calibration observations** + 2 parameters (η_lo, κ) or 0 (Option C fixed η).
**Available but excluded**: 4 Hornick vaccine rows (CoP unmappable), 8 Oxford shedding rows (Tier 1 only; restored in Tier 2).
**Effective independent observations** (per Reviewer 2): Tier 1 ~27-29; Tier 2 ~31-33, after correcting for within-group correlation.
**Total parameters**: Tier 1: 6 biological + 4 nuisance + 1 overdispersion = 11. Tier 2: +2 (η_lo, κ) = 13.
**Data-to-parameter ratio**: Tier 1 ~2.3:1; Tier 2 ~2.4:1. Both require informative priors. Tier 2 is preferred for γ_inf identification.
**Canonical data source**: `dose_response_data.csv` in this directory. All observation counts verified against this file.
**Definition warnings**: (1) Oxford shedding **EXCLUDED** (Section 2.6): treatment-truncation bias makes shedding < diagnosis at all doses ≥10⁴. Oxford contributes fever only. (2) Levine/Gilman/Hornick fever definitions **RESOLVED** (Section 2.5): study-specific φ(T) values derived from Oxford threshold ladder — Hornick φ≈0.25, Levine φ≈0.65, Gilman φ≈0.65. No extra parameters needed. (3) Darton S1 fever thresholds (T38, T39) are 1-2 subjects lower than published Table 2 — likely a minor definition difference in how the S1 encodes temperature events vs the publication's analysis.

---

## Appendix: Prerequisites Before Implementation

**~~PREREQUISITE 1~~** (**RESOLVED 2026-03-18**): The Hornick 10³ discrepancy does not exist. Independent PDF verification confirmed both Hornick 1966 Figure 2 and Hornick 1970 Table 1 show 0/14 at 10³. The extraction's 9/14 was an AI hallucination. The Hornick_1966.md extraction has been corrected.

**~~PREREQUISITE 2~~** (**RESOLVED 2026-03-24**): The Gilman Table 5 cross-tabulates Well/Ill × H-antibody status. The original extraction misread the "Well" row (14, 13) as denominators. Correct denominators: H<1:20 = 36 (14 well + 22 ill), H≥1:20 = 17 (13 well + 4 ill). Total with H-antibody data: 53 of 64 controls (83%). The 11 missing are likely controls without baseline sera. Strata are well-powered (not a thin subsample). π_susc prior updated from Beta(5,5) to Beta(7,4) to reflect 68% susceptible.

**~~PREREQUISITE 3~~** (**RESOLVED 2026-03-18**): Oxford infection definitions now harmonized on shedding-only across all three studies. Darton shedding-only extracted from S1 individual-level data: Placebo 19/30 (63.3%), M01ZH09 20/31 (64.5%), Ty21a 18/30 (60.0%). Jin confirmed as shedding-only (22/31, 22/37, 21/35). Waddington shedding (13/20 at 10³, ?/20 at 10⁴). All consistent.

---

## Appendix A: Gamma Portability Across Immunity Sources (added per Reviewer 2)

The model applies a single γ per outcome across all immunity sources: Vi vaccination (Oxford), whole-cell killed vaccination (Maryland), natural infection (Gibani rechallenge), and unknown background immunity (Maryland mixture). This assumes the mechanism by which immunity reduces susceptibility is the same regardless of how immunity was acquired.

**Evidence for portability**: Both Vi vaccination (Darton 2016, anti-Vi HR = 0.29) and prior natural infection (Gibani 2020, RR = 0.64) show protective effects. Both reduce shedding (Gibani 2019, anti-Vi: P<.0001; prior exposure: OR 0.33). The direction and rough magnitude are consistent.

**Evidence against portability**: Vi vaccination induces primarily anti-Vi IgG. Natural infection induces a broader response (T cells, mucosal IgA, anti-O, anti-H). The Gibani 2020 paradox (non-sick subjects more protected on rechallenge) suggests innate susceptibility, not adaptive immunity, dominates in the rechallenge setting. Hornick 1970 found NO correlation between O/H/Vi agglutinins and protection in Maryland.

**Implication**: The estimated γ is an effective parameter that averages over multiple immunity mechanisms. It may not generalize to immunity sources not represented in the calibration data (e.g., endemic natural exposure in children). This limitation should be stated explicitly when using the model for field predictions.

---

## Appendix B: Alternative Model Configurations (was Appendix A)

### A1: Shared α Model (8 free parameters)

Set $\alpha_{\text{inf}} = \alpha_{\text{fev}} = \alpha$. Removes 1 parameter. Justified if host heterogeneity in susceptibility is the same for infection and fever. Test by comparing WAIC/LOO-CV.

### A2: Shared γ Model (8 free parameters)

Set $\gamma_{\text{inf}} = \gamma_{\text{fev}} = \gamma$. Removes 1 parameter. Justified if immunity affects infection and fever equally. Contradicted by prior expectation (immunity should affect fever more than infection — "leaky protection"). Test by comparing to separate-γ model.

### A3: Fixed-Ratio Model (7 free parameters)

Following the prototype: set $N_{50,\text{inf}} = N_{50,\text{fev}} / r_{N50}$, $\alpha_{\text{inf}} = \alpha_{\text{fev}} \cdot r_\alpha$, $\gamma_{\text{inf}} = \gamma_{\text{fev}} / r_\gamma$, with $r_{N50}$, $r_\alpha$, $r_\gamma$ either fixed or estimated. Reduces to 3 biological parameters + 3 ratios. More parsimonious but less flexible.

### A4: Susceptibility Class Extension (+2 parameters)

Add a fraction $\pi_{\text{innate}}$ of innately resistant individuals with a much lower dose-response (or a maximum susceptibility ceiling <1). Motivated by Gibani 2020 paradox. Adds $\pi_{\text{innate}}$ and $P_{\max,\text{resistant}}$ (or a resistant-class N50 multiplier). Consider only if the base model shows systematic misfit for rechallenge cohorts.

---

## Appendix C: What This Plan Does NOT Cover

1. **Titer dynamics**: The CoP = g(anti-Vi IgG) mapping is taken as given from a separate model. Uncertainty propagation from the titer model into this inference is desirable but not specified here.

2. **Field predictions**: Converting bicarb-unit predictions to field setting predictions requires additional parameters (exposure rate, field dose distribution) that are outside the scope of this CHIM-data calibration.

3. **Incubation period modeling**: Time-to-event data is available (Hornick Table 1, Waddington Table 3) and could provide additional constraints on α and the dose-response shape. Not included here but could be added as an additional likelihood layer.

4. **Seroconversion modeling**: Per the YOLO notes, seroconversion is too conditional on disease status to model as a direct dose-response. Reserved for validation.

5. **Non-Quailes strains**: Hornick Table 2 includes Zermatt, Ty2V, O-901, and Ty2W strains. These are excluded from calibration (different virulence) but could inform strain-specific random effects in a future extension.
