# Peer Review: Anti-Vi Metrology Bridge YOLO Implementation

**Reviewer**: Reviewer 2 (Adversarial but Fair)
**Date**: 2026-02-04
**Manuscript**: Anti-Vi Metrology Bridge from Historical HA Titers to VaccZyme ELISA

---

## Summary Recommendation

**MAJOR REVISION**

The authors present an ambitious attempt to bridge historical Vi hemagglutination titers (1980s) to modern VaccZyme ELISA values through a multi-step cascade. While the documentation is commendably transparent about limitations, and the statistical implementation is technically sound, **the E2 structural assumption is so consequential and so weakly justified that the entire bridge rests on quicksand**. The work is publishable with major revisions, but the current framing overstates what this bridge can deliver.

---

## Major Concerns

### 1. The E2 Structural Assumption is Arbitrary and Dominates the Result

**Severity: Critical**

The authors acknowledge that E2 (Barrett ELISA to Modern In-house ELISA) has "NO direct data support" and accounts for ~70% of total variance. However, the chosen tau = 1.0 log unit is not justified by any empirical data or formal elicitation.

The decision log states tau is composed of:
- Inter-assay CV: ~20-30% (reasonable)
- Coating method differences: ~0.5 log (where does this come from?)
- Population differences: ~0.5 log (why this value?)

**Problem**: These components are not additive on the log scale in the way implied. If they were independent multiplicative errors, you would add variances, not standard deviations. The "conservative rounding up" to 1.0 is hand-waving.

**Required change**: Either (a) provide a formal uncertainty budget with documented sources for each component, or (b) present results across a wide range of tau values (0.5, 1.0, 1.5, 2.0) as equally plausible scenarios, not as sensitivity analysis around a "default."

### 2. The Calibration Assumption "Barrett titer 160 ~ 10 ug/ml" Has No Foundation

**Severity: Critical**

The end-to-end bridge code contains this comment:

```r
# Simplification: Assume Barrett titer 160 ~ 10 ug/ml (midrange responder)
```

This is a critical parameter that anchors the entire bridge. The stated rationale is:

> "Vi-IgGR1,2011 at 33 ug/ml has titer ~160-320 based on context"

What context? The Vi-IgGR1,2011 reference did not exist in 1983 and has never been run on the Barrett 1983 assay. This "calibration" is fabricated.

**Impact**: An error of 2-5x in this calibration (acknowledged as possible in the caveats) would shift all predictions by the same factor, on top of the E2 uncertainty.

**Required change**: Either (a) remove this calibration and present results in Barrett ELISA units only (breaking the bridge at E2), or (b) quantify this calibration uncertainty as a separate term in the cascade, or (c) find actual literature values for Barrett ELISA sensitivity in ug/ml terms.

### 3. The Barrett 1983 Data Extraction Contains Potential Errors

**Severity: Moderate**

The raw counts table in `barrett_1983_extraction.md` shows:

| HA | <20 | 20 | ... | Row Total |
|----|-----|----| ... |-----------|
| <20 | 37 | 4 | ... | 42 |

But the CSV file shows (ELISA, HA, count):
- (10, 10, 37) - <20 ELISA, <20 HA
- (20, 10, 4) - 20 ELISA, <20 HA

This looks like the axes may be transposed. The extraction markdown shows HA as rows and ELISA as columns, but the description says "HA, Y-axis" and "ELISA, X-axis." In Barrett 1983 Figure 2, **HA is on the Y-axis (vertical)**.

Let me verify with the data:
- Extraction says "Positive by ELISA (>=20): 40" and "Positive by HA (>=20): 35"
- From the CSV, summing count where ELISA_titer >= 20: 4+2+1+1+2+1+1+2+2+2+4+5+1+1+2+1+1+1+1+1 = 36
- Summing count where HA_titer >= 20: 1+1+1+1+2+1+1+2+2+2+4+5+1+1+2+1+1+1+1+1 = 32

**This does not match the paper's reported values of 40 and 35.**

**Required change**: Re-verify the data extraction against the original figure. Provide a clear mapping showing which CSV column corresponds to which axis in the original figure. The sum of positives must match the paper's reported values.

### 4. The E1 Model Specification Ignores the Data Structure

**Severity: Moderate**

The Barrett 1983 data are aggregated cell counts, yet the analysis uses `uncount()` to expand them to individual observations and then fits OLS regression. This treats all observations as independent and ignores:

1. **Left-censoring**: Values coded as <20 (set to 10) are left-censored, not point measurements. A Tobit model or survival regression would be more appropriate.

2. **Discreteness**: Both axes are 2-fold dilution series. The data are doubly discrete, yet OLS assumes continuous outcomes.

3. **Measurement error in predictor**: The HA titer has measurement error too, yet standard regression assumes the predictor is known exactly.

The impact is likely modest given the large scatter, but the authors should acknowledge these limitations.

**Required change**: Either (a) fit a model appropriate for censored, discrete data (e.g., ordinal regression, Tobit), or (b) justify why OLS is adequate and discuss how these issues might bias estimates.

### 5. The Lee 2020 E3 Bridge Parameters Cannot Be Verified

**Severity: Moderate**

The authors state they attempted to fetch individual-level data from Figshare but used published parameters as a fallback. However:

- The residual SD (sigma = 0.30 log units) is "estimated from r," not from actual residuals
- The confidence intervals on intercept and slope are marked "Not reported"
- No validation that the published regression was fit correctly is possible

The derivation of sigma from r requires assuming the residuals are homoscedastic and the relationship is truly linear. Given r = 0.991, small deviations from these assumptions could substantially change sigma.

**Required change**: (a) Document the Figshare attempt - was data truly unavailable or just not sought? (b) Present sensitivity to sigma_E3 values, or (c) clearly state that E3 uncertainty may be underestimated.

### 6. Population Mismatch is More Serious Than Acknowledged

**Severity: Moderate**

The bridge combines:
- E1: Acute typhoid carriers (natural infection, diagnostic setting, 1980s El Salvador)
- E3: Healthy vaccinees (controlled trial, 2010s Philippines)

The caveats mention "avidity maturation differs" and "IgG subclass distribution differs," but this understates the problem:

1. **Natural infection produces IgM early**: Barrett's HA assay may detect IgM preferentially. His ELISA is IgG-specific. The HA-ELISA relationship from carriers may not transfer to vaccinees at all.

2. **Antibody levels differ dramatically**: Carrier titers of 640-5120 were common in Barrett's data. Lee 2020 vaccinees had median ~30 ug/ml. The bridge extrapolates from a carrier-dominated calibration to a vaccinee application.

3. **Epitope differences**: Natural infection antibodies may recognize different Vi epitopes than vaccine-induced antibodies, affecting assay comparability.

**Required change**: Add explicit caveat that the E1 data from carriers may be fundamentally non-transferable to vaccinee samples, not just systematically offset.

---

## Minor Concerns

### 7. Citation Inconsistency in Barrett Reference

The extraction file cites "Barrett TJ et al. J Clin Microbiol 1983;17(4):625-627" but the README cites "J Clin Microbiol. 1983;18(3):625-632."

The correct citation appears to be J Clin Microbiol 1983;18(3):625-632 (September 1983). Please verify and correct.

### 8. Log Scale Confusion in E3

The Lee 2020 extraction states the regression uses "natural log" but also provides "0.13 log10 units." The conversion factor used (divide by ln(10)) is correct, but the presentation is confusing. State consistently whether the cascade operates on natural log or log10.

In fact, in `end_to_end_bridge.R`, the E2 tau is specified in "log units" but the conversion chain uses different bases:
- E1: log2 scale
- E2: "log10" per bridge_graph.md, but "log" (ambiguous) in code
- E3: natural log

**Required change**: Add a unit conversion table that explicitly shows which log base is used at each step and the conversion factors applied.

### 9. Incomplete Stop Criteria Evaluation

The decision log states "Posterior predictive check fails: TBD" but the final caveats document says "N/A - No validation data available."

The posterior predictive check should evaluate whether the model's predictions are consistent with the observed data used to fit it, not external validation data. You could check:
- Does the E1 model reproduce the observed scatter?
- Are prediction intervals calibrated (e.g., 95% PI contains 95% of data)?

**Required change**: Perform a posterior predictive check on the E1 model at minimum, and document results.

### 10. Sensitivity Analysis is Incomplete

The tau sensitivity analysis shows distributions for HA = 160 only. This hides the fact that uncertainty propagation may be non-linear across the HA range.

**Required change**: Show uncertainty as a function of input HA titer, not just at a single point.

### 11. Missing Variance Decomposition Details

The uncertainty decomposition shows "SD_log" and "fold_95" for three scenarios but does not break down variance into E1/E2/E3 components directly. The claim that "E2 contributes ~70% of total variance" should be derived from actual variance arithmetic, not eyeballed.

**Required change**: Provide explicit variance decomposition calculation: Var(total) = Var(E1) + Var(E2) + Var(E3) + 2*Cov(E1,E2) + ...

### 12. Output File Naming

`parameter_estimates.csv` is misleading - it contains bridge prediction summaries, not parameter estimates. The actual E1 and E3 parameters are in `intermediates/`.

### 13. Missing Date Context in Paper Citations

Barrett 1983 was published before ELISA standardization efforts. Lee 2020 was published during the NIBSC 16/138 era. This 37-year gap is mentioned but its implications for assay comparability deserve more emphasis.

### 14. Typo in README

"E1 data is acute typhoid carriers" should be "E1 data are from acute typhoid patients/carriers."

---

## Specific Required Changes

1. **Re-verify Barrett 1983 data extraction** against original figure, with explicit axis labeling
2. **Justify or abandon the "titer 160 = 10 ug/ml" calibration**
3. **Present tau as a scenario parameter**, not a default with sensitivity analysis
4. **Add calibration uncertainty as a separate cascade term**
5. **Clarify log base conventions** throughout the cascade
6. **Perform posterior predictive check** on E1 model
7. **Expand sensitivity analysis** to show uncertainty across HA titer range
8. **Quantify variance decomposition** with explicit arithmetic
9. **Strengthen population non-transferability caveat**
10. **Correct citation inconsistencies**

---

## Questions for Authors

1. Why was tau = 1.0 chosen rather than, say, tau = 0.7 or tau = 1.5? What would change your choice?

2. Have you attempted to contact the original Barrett 1983 authors or their institutions to locate archived samples or assay protocols?

3. The E3 data from Lee 2020 may be available in supplementary materials or by author request. Did you attempt this?

4. Could you present an alternative "two-bridge" formulation where E1 (HA to Barrett ELISA) and E3 (In-house to VaccZyme) are kept separate, acknowledging that E2 cannot be crossed?

5. What would be the minimum sample size needed for a proper E2 bridging study, given the current uncertainty estimates?

6. The Rijpkema 2018 collaborative study included 7 in-house ELISAs. Were any of these descended from or comparable to the Barrett 1983 method? Could that study inform E2?

---

## Positive Aspects

Despite the major concerns, this work has several strengths worth acknowledging:

1. **Exceptional transparency**: The caveats document is brutally honest about limitations
2. **Well-documented search strategy**: The E2 gap was systematically confirmed, not assumed
3. **Appropriate use of Monte Carlo**: The cascade propagates uncertainty correctly
4. **Clear visual presentation**: The bridge graph and conversion table are helpful
5. **Reproducible code**: The R scripts are well-documented and readable
6. **Thoughtful decision logging**: The rationale for each choice is documented

This is a creditable attempt at an essentially impossible bridging problem. With revisions addressing the major concerns, it could serve as a useful "rough guide" with appropriately massive uncertainty bounds.

---

*Reviewer 2 recommends Major Revision. The work is not rejectable because the authors have been transparent about limitations, but the current framing implies more precision and certainty than the data support.*
