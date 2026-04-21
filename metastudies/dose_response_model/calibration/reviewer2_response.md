# Author Response to Reviewer 2

**Manuscript**: Joint Inference Plan: Typhoid Dose-Response Model Calibration
**Date**: 2026-03-17

We thank Reviewer 2 for a thorough and constructive review. The review identified several genuine issues that strengthen the plan. We address each concern below and indicate where the plan has been revised.

---

## Responses to Major Concerns

### MC1: The phi-delta-CoP_imm three-way confound

**The reviewer is correct.** We understated the severity of this confound. The constant-phi approximation is particularly dangerous because phi is calibrated at ~10⁴ bicarb (one Oxford dose), then applied to Maryland data spanning 10³-10⁹ milk (which, after δ conversion, maps to ~10⁻¹ to 10⁵ bicarb-equivalent — a 6-log range).

**Revised approach**: We accept the reviewer's recommendation to make dose-dependent phi the default. Specifically:

$$\phi(D_{\text{eff}}) = 1 - (1 - \phi_{\min}) \cdot e^{-D_{\text{eff}} / D_{\phi}}$$

where $\phi_{\min}$ is the baseline definition sensitivity at low effective doses and $D_{\phi}$ is the effective dose at which $\phi$ approaches 1. This adds 1 parameter (2 total for the definition model: $\phi_{\min}$ and $D_{\phi}$), but the Oxford multi-threshold data directly constrains both: the decline from 77% (composite) to ~30% (≥39°C, from Darton placebo) pins $\phi_{\min}$, while the known convergence at high doses pins $D_{\phi}$.

However, we note that even with dose-dependent phi, the three-way confound remains at the level of the Maryland ID50 position. **We will add** to Section 8.3 an explicit discussion of the phi-delta-CoP_imm surface, and we will present the 2D grid sensitivity analysis (phi_min × delta) as recommended.

An alternative approach, which we now discuss in the plan, is to **fit the Maryland data to the strict-fever model directly** — i.e., fit a third dose-response curve for "strict fever" (≥39.4°C) that is constrained to lie below the Oxford-definition fever curve. This avoids phi entirely by treating the two definitions as two different outcomes, each with their own N50. The Maryland-Oxford bridge then happens through the infection curve (which has minimal definition differences) and the shared biological parameters. This is arguably cleaner than a phi correction.

### MC2: Double-counting in independent curves and Hornick Table 2

**The reviewer correctly identifies a double-counting issue.** At 10⁷ milk, the Hornick Table 1 entry (H-F-7: 16/32 fever) and the Table 2 data (H-I-7: 28/30 infected, H-FgI-7: 16/28 fever|infected) partially overlap. The 30 in Table 2 are a subset of the 32 in Table 1 (the 2 extra subjects in Table 1 are likely from other Vi-containing strains contributing to the total count).

**Resolution**: We will use Table 2 as the canonical source at 10⁷ and **remove H-F-7** from the marginal fever likelihood. The Table 2 data enters as a joint observation:

- H-I-7: 28/30 infected (constrains infection curve)
- H-FgI-7: 16/28 disease given infected (constrains the ratio P_fev/P_inf)

This eliminates the double-counting.

**On the cascade vs independent formulation**: We accept that the cascade formulation P(fever|infected) is more principled. However, implementing it properly requires individual-level crossover data, which we have only at one dose. We propose a **compromise**: fit two independent curves for the Oxford data (where we lack crossover data), but use the **cascade constraint** specifically at the Maryland 10⁷ data point. This amounts to the original plan but with the double-counting fix above.

**On within-group infection-fever correlation**: The reviewer is correct that the infection and fever observations within the same group (e.g., D-I-plac and D-F-plac) come from the same subjects. A joint bivariate model (2×2 table: infected+fever, infected+no-fever, not-infected+no-fever) would be more appropriate. However, we lack the individual-level crossover counts for most Oxford groups — the publications report only marginal rates. Where cross-tabulation data exists (Waddington diagnosis breakdown), we should use it. **We will add a note** to the plan recommending bivariate modeling where cross-tabulated data is available, and acknowledging the independence approximation as a simplification elsewhere.

### MC3: Gamma identification and the interaction web

**This concern is well-taken.** The reviewer correctly notes that gamma_fev rests on thin data if the Darton HR is used only for validation. Let us clarify the intended approach:

Gamma is identified from the **group-level attack rate contrasts** across immunity levels at fixed dose:
- Jin control (77%) vs Vi-TT (35%) vs Vi-PS (37%) → 3 groups spanning CoP ~1 to ~500 EU/mL
- Darton placebo (67%) vs Ty21a (43%) → 2 groups at different CoP
- Gibani rechallenge (44%) vs naive (63%) → natural immunity contrast

That's 7-8 group-level observations across 3-4 immunity levels at ~10⁴ bicarb, constraining one parameter (gamma_fev). This is not lavish but is more than the "3 binomial contrasts" the reviewer suggests. The Darton HR is used as a **cross-check**: if the model-predicted dose-response slope with respect to CoP matches the observed HR, we have confidence in gamma.

**We will add** to the plan: (a) the expected posterior correlation matrix from a Laplace approximation (Section 8.3, new subsection), and (b) the gamma marginal from Stage 1 vs Stages 2-3 comparison (Section 10.4).

### MC4: The Hornick 10³ data point (0/14) and the 1966 discrepancy

**We agree this is more consequential than originally acknowledged.** The 0/14 is highly influential because the binomial likelihood strongly penalizes p > 0.05 for this observation.

**On the 1966 vs 1970 discrepancy**: We have NOT independently verified this against the primary publications. The extracts note that the 1966 paper reports 9/14 and the 1970 review reports 0/14 at the same dose. Possible explanations:
1. The 1966 paper's 10³ dose group may have included streptomycin-pretreated subjects (Hornick 1970 Part 1, p.689 notes that 1/4 streptomycin-pretreated subjects at 10³ became ill)
2. The 1966 figure may be misread or mislabeled
3. The 1970 paper may use a different criterion for counting cases

**We will add**: (a) a note that the discrepancy must be verified against the primary PDFs, (b) a three-way sensitivity analysis: 0/14, excluded, and 9/14, as the reviewer recommends, and (c) the beta-binomial treatment of the 10³ observation as an option if the discrepancy cannot be resolved.

### MC5: Gibani susceptibility paradox integration

**We agree this should be a specific hold-out prediction, not just a deferred sensitivity analysis.**

**Revised plan**: Add the Gibani rechallenge subgroup data (prior disease: 25/37 = 68%, no prior disease: 10/38 = 26%) as explicit **hold-out observations** in the model checking plan. The model is fit without these observations, then asked to predict the split. If the model predicts the wrong direction (prior disease should be more protected under adaptive immunity, but the data show the opposite), this is strong evidence for susceptibility classes.

We also accept the suggestion to discuss whether gamma_from_vaccination and gamma_from_natural_exposure should be the same. The biological argument: Vi vaccination induces anti-Vi IgG specifically. Natural infection induces a broader immune response (T cells, mucosal IgA, anti-O, anti-H). These may protect through different mechanisms, and gamma may not be portable across immunity sources. This is a genuine limitation that the plan should acknowledge.

### MC6: Effective number of independent observations

**We accept this concern and will revise the count.**

Corrections:
- Gilman: Gil-F-ctrl (64) partially overlaps with Gil-F-Hlo (14) + Gil-F-Hhi (13). We should use EITHER the pooled observation OR the stratified observations, not both. **Resolution**: Use the stratified observations (Gil-F-Hlo, Gil-F-Hhi) plus a residual observation for the remaining ~37 subjects without H-antibody measurement: Gil-F-rest (n≈37, y≈19, p≈51%). Alternatively, if H-antibody was measured in all subjects, replace the pooled with the two strata.
- Hornick at 10⁷: Remove H-F-7 and use H-I-7 + H-FgI-7 (as discussed in MC2).
- Oxford infection+fever from same group: Acknowledge as an approximation; recommend bivariate modeling where cross-tabulated data is available.

**Revised count**: ~28-30 effectively independent observations for 11-13 parameters. The ratio is tight (2.5:1) but comparable to other CHIM dose-response calibrations in the literature (e.g., Haas 1999 β-Poisson fits for Salmonella). The Bayesian framework with informative priors can handle this ratio, but we should be transparent about it.

---

## Responses to Minor Concerns

**MC-m1 (alpha prior)**: Agreed. We will add a specific prior predictive check: "verify that P(D=1, CoP=1) < 10% under the prior" to exclude pathological near-flat curves.

**MC-m2 (random effect parameterization)**: Agreed that logit-normal is more standard. We will switch to a logit-scale random effect: logit(p_j) = logit(p_model) + epsilon_s, epsilon_s ~ Normal(0, sigma_study). This also has better behavior near 0 and 1.

**MC-m3 (Waddington 10⁵ anomaly)**: Noted. At n=5, this is noise. We will add a footnote acknowledging this and noting it as expected from the imperfect sensitivity of stool culture for detecting infection.

**MC-m4 (Jensen's inequality for GMT)**: Valid point. For the initial model, we'll use the GMT approximation and verify the bias is small post-hoc by computing E[P(D, CoP)] under a log-normal titer distribution with the reported variance. If the bias exceeds 5 percentage points, we'll switch to numerical integration.

**MC-m5 (notation consistency)**: We will audit the notation throughout and ensure D/δ is used consistently. The prototype's absorbed-offset convention (no explicit δ) will be clearly distinguished from the plan's explicit-δ convention.

**MC-m6 (Levine shedding data)**: Good catch. Lev-I-1 (19/26 shedding at 10⁵ milk) provides the only Maryland infection data outside Hornick Table 2 at 10⁷. Combined with Gil-I-ctrl (26/43 shedding at 10⁵) and H-I-7 (28/30 infection at 10⁷), this gives 3 infection data points at 2 milk doses — enough to contribute to N50_inf estimation through the δ bridge. We will acknowledge this in the identifiability analysis.

---

## Responses to Questions

**Q1 (phi in the conditional)**: Good catch. Phi should apply to the *observation* of fever, not to the underlying biology. In the conditional P(disease|infected), the numerator is the Maryland-observed disease rate (which has phi applied) and the denominator is the infection rate (no phi). So the conditional IS phi * P_fev / P_inf. This is correct: we're modeling the probability that an infected person is *observed to have disease under the Maryland criterion*, which is the underlying fever probability times the probability that the fever is severe enough to be captured. So yes, phi appears in the conditional intentionally.

**Q2 (Oxford infection definition heterogeneity)**: The Waddington infection definition is shedding-only; Darton/Jin use bacteremia OR shedding (broader). The difference is modest: at ~10⁴, Waddington shows 63% shedding while Darton shows 87% by the broader definition. This ~24-point gap is real and systematic. **We should model this.** Options: (a) use shedding consistently (requires extracting shedding-only rates from Darton/Jin — the Gibani 2019 pooled analysis may help), or (b) add a definition indicator for infection analogous to phi for fever. We will add this as a note in the plan. For the initial model, using shedding consistently across studies is the cleaner approach.

**Q3 (CoP mapping)**: The titer model gives CoP as a transformed anti-Vi IgG titer. The simplest mapping is CoP = max(1, titer / LOD) where LOD is the limit of detection (~7.4 EU/mL). More sophisticated mappings (logistic transform, power transform) are possible but introduce additional parameters. **For this calibration, we treat the CoP mapping as fixed** from the titer model and estimate gamma conditional on that mapping. If the titer model is updated, gamma must be re-estimated. This is analogous to a two-stage estimation approach — valid if the titer model has much more data than the dose-response model (which it does, from longitudinal vaccine trial serology).

**Q4 (Gilman H-antibody subsample)**: The 27 H-antibody-measured subjects appear to be a subset. We do not know if they are representative. **We should assume they are a convenience subsample** and note this uncertainty. The prior on pi_susc (Beta(5,5)) is centered at 0.5 but allows wide deviation, which accommodates the possibility that the subsample is biased.

**Q5 (1966 vs 1970 verification)**: Not verified against primary sources. **This must be done before implementation.** We will flag this as a prerequisite in the plan.

---

## Summary of Plan Revisions

Based on this review, the following changes will be made to the joint inference plan:

1. **Section 2.5**: Dose-dependent phi becomes the default; constant phi becomes the simplification option. Add discussion of "strict fever as a third outcome" as an alternative to phi correction.
2. **Section 5.6**: Remove H-F-7 from the marginal fever likelihood to eliminate double-counting; use H-I-7 + H-FgI-7 as the canonical 10⁷ data.
3. **Section 5.7**: Clarify that Gil-F-Hlo + Gil-F-Hhi replace Gil-F-ctrl (not supplement it).
4. **Section 8.3**: Add three-way confound discussion (phi × delta × CoP_imm); add expected posterior correlation matrix.
5. **Section 9.2**: Add Gibani rechallenge subgroup as explicit hold-out prediction.
6. **Section 9.3**: Promote Hornick 10³ sensitivity to a required analysis (3-way: 0/14, excluded, 9/14).
7. **Section 10.4**: Add Stage 1 → Stage 2 gamma shift diagnostic.
8. **Section 3.2**: Add 10³ verification against primary PDFs as a prerequisite.
9. **Minor**: Switch study-level random effect to logit scale; add prior predictive check for alpha; acknowledge Jensen's inequality for GMT; use consistent shedding definition across Oxford studies.
10. **New Appendix**: Discussion of gamma portability across immunity sources (vaccine vs natural vs innate).
