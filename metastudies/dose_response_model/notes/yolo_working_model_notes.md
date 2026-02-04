# YOLO Working Model Notes

**Date**: 2026-02-04
**Purpose**: Claude's first-pass intuition on collapsing the reference model to something fittable, and which papers power what.

---

## My Proposed Reduced Working Model

### The Core Simplification: One Curve, Medium Offset, Immunity Scaling

After reading through everything, here's my intuition: the minimal model that captures the important dynamics is:

$$
P(\text{typhoid fever} | D, \text{medium}, \text{CoP}) = 1 - \left(1 + D \cdot \delta_{\text{medium}} \cdot \frac{2^{1/\alpha} - 1}{N_{50}}\right)^{-\alpha / \text{CoP}^{\gamma}}
$$

Where:
- **$N_{50}$**: ID50 in "milk-equivalent" CFU (pick one medium as reference)
- **$\alpha$**: heterogeneity parameter (aggregates everything: gastric survival variance, host factors, etc.)
- **$\delta_{\text{medium}}$**: multiplicative offset for delivery medium (milk = 1, bicarbonate = $\delta$)
- **$\text{CoP}$**: scalar correlate of protection (1 = naive; >1 = some immunity)
- **$\gamma$**: how steeply immunity scales the effective dose

### Why This and Not the Full Cascade?

The reference model has this beautiful chain: challenge dose → gastric dose → colonization → systemic invasion → acute disease → clinical typhoid diagnosis. But:

1. **We can't see the intermediate states.** Colonization and systemic invasion are latent. Stool shedding is a noisy proxy for colonization, bacteremia for systemic invasion, but they're not the states themselves.

2. **The product of sigmoids looks like a sigmoid.** Per Section 6.4 of the spec, cascaded monotonic dose-response curves collapse to something that *looks* like a single sigmoid with finite data. Unless we have specific reasons to break them apart, we can't identify them separately.

3. **The one place we might separate the cascade**: if infection (bacteremia/shedding) and disease (fever) have different dose-responses. The Oxford data hints at this—shedding rates exceed typhoid diagnosis rates at low doses. But...

**My call**: Start with a single curve to "typhoid fever" (the composite endpoint). If there's signal for bifurcating infection vs disease, we can revisit. The Maryland studies defined their endpoint as clinical typhoid (the thing that got treated), and that's the outcome we care about for public health anyway.

---

## The Delivery Medium Problem

This is the single biggest issue in the data. The Maryland and Oxford data look incompatible at face value:

| Study Era | Dose | Attack Rate | Medium |
|-----------|------|-------------|--------|
| Maryland | 10^5 | ~28% | Milk |
| Oxford | 10^4 | ~60-65% | Bicarbonate |

A naive pooling would be nonsense. The options are:

1. **Fit eras separately, compare.** Safe but loses power.
2. **Treat medium as a fixed offset ($\delta_{\text{medium}}$).** Pool with a multiplicative correction.
3. **Treat medium as a random effect.** Let era-level differences absorb medium + everything else.

**My intuition**: Go with (2), a fixed offset. Bicarbonate neutralizes gastric acid, so more bacteria reach the gut. This is biologically grounded. The offset appears to be ~3-4 log10—i.e., 10^4 in bicarbonate behaves like 10^7-10^8 in milk. That's consistent with:
- Oxford 10^4 bicarb → ~60-65% attack = Maryland 10^7-10^8 milk
- Oxford 10^3 bicarb → ~55% attack ≈ Maryland 10^7 milk

**Concern**: We're estimating a 3-4 log offset from era-level comparisons with different populations, different fever thresholds, 50 years apart. This is a huge extrapolation. But it's the only way to use both datasets.

---

## Outcome Definition Heterogeneity

Maryland: ">103F for 24-36hr + symptoms" (strict)
Oxford: ">=38C for 12hr OR bacteremia" (permissive)

These are different outcomes! Oxford catches more subclinical disease.

**Options**:
1. **Treat as different outcome variables.** Fit "strict disease" for Maryland, "any typhoid diagnosis" for Oxford.
2. **Assume they're observing the same latent state with different sensitivity.** Add measurement model.
3. **Use bacteremia as common denominator.** Both eras measured it.

**My intuition**: Bacteremia is the cleanest common endpoint. Both eras did daily blood cultures. The Oxford rates:
- 10^3: 40% bacteremic
- 10^4: 56% bacteremic

Compare to Maryland 10^5 milk: ~28% "ill" (but bacteremia rates not always separated from clinical disease).

Actually, looking at Hornick 1966 and related—they often equate "blood culture positive" with early disease detection. The "illness" rate *is* approximately the bacteremia rate for those who progressed to symptoms.

**Revised intuition**: The outcome definitions may be closer than they look. Maryland's strict fever threshold meant they only treated (and thus "diagnosed") severe cases—but those are the bacteremic ones. Oxford's permissive threshold catches more people, but ~95% of their diagnosed cases were bacteremic anyway (per Gibani 2020).

**Pragmatic call**: Fit to "typhoid fever/diagnosis" as defined by each study. Accept that there's ~10-20% definition slop. The medium offset will absorb some of this anyway.

---

## Immunity Representation

The reference model imagines CoP as an array collapsed to a scalar. What do we actually have?

### Maryland era
- Almost nothing. Hornick 1966 notes O/H agglutinins didn't correlate with protection.
- Gilman 1977: H antibody >=1:20 was protective in controls (61% → 24% attack rate).
- Woodward 1980: Prior military service (a proxy for "probably vaccinated or exposed") reduced attack rate 20% vs 48%.
- Dupont 1971: Prior infection (23% vs 30%) showed modest protection.

### Oxford era
- Darton 2016: Anti-Vi IgG → 1 log increase = 71% reduction in hazard. This is the best signal.
- Gibani 2020: Prior challenge reduced attack rate 63% → 44% (homologous). But confounded by vaccine history and prior disease status.
- Juel 2018: SBA correlates with severity, not protection.

### The weird Gibani 2020 finding
People who *didn't* get sick on first challenge were more protected on rechallenge (26% vs 68% attack rate). This is backwards from "adaptive immunity from infection" and suggests innate susceptibility factors dominate.

**My intuition on CoP**:

For the working model, I'd treat CoP as binary or ordinal:
- **Naive**: CoP = 1
- **Some immunity** (prior vaccination OR prior exposure without disease): CoP = ? (estimate)
- **Prior disease**: CoP = ? (may actually be ~1, per Gibani finding)

Or maybe even simpler: just have a "naive" vs "not naive" indicator and estimate the risk ratio.

**What I'm uncertain about**: The Darton 2016 anti-Vi IgG result is the strongest signal, but ~29% of their "naive" UK subjects had detectable anti-Vi. That's weird. Cross-reactivity? Subclinical exposure? If we take it at face value, we could use anti-Vi titer as a continuous CoP predictor—but we don't have titers for Maryland, so it's only useful for Oxford calibration.

---

## Which Papers Power Which Parameters

### $N_{50}$ (ID50)

**Primary data**:
- **Hornick 1966/1970**: Multi-dose (10^3-10^9). The only data spanning the full dose-response curve. ~1700 subjects total across the Maryland program. This is the anchor.
- **Waddington 2014 (outpatient)**: Dose escalation 10^3-10^5 in bicarbonate. N=41. Provides modern anchor for bicarbonate.

**Validation**:
- Glynn 1995: Reanalysis of Maryland data, adds incubation period.
- Control arms from vaccine trials: Levine 1976 (N=97 at 10^5), Gilman 1977 (N=64 at 10^5).

**My concern**: The Hornick 10^3 data point is weird (64% attack rate at 10^3, higher than 10^5). Small N=14. Probably ignore it or treat as outlier.

### $\alpha$ (heterogeneity)

Estimated jointly with $N_{50}$ from the shape of the dose-response curve. Requires data across multiple doses.

**Primary data**: Hornick 1966/1970 (shape at high doses), Waddington 2014 (shape at low doses with bicarbonate offset).

**Identifiability concern**: $\alpha$ and $N_{50}$ will be correlated. Need sufficient dose range to separate them. Hornick spans 6 logs—that should be enough. Waddington only spans 2 logs, which is thin.

### $\delta_{\text{medium}}$ (medium offset)

**Primary data**: Cross-era comparison.
- Maryland at 10^5: ~28%
- Oxford at 10^4: ~60%

If these are the same underlying curve, $\delta$ is what makes them line up.

**My concern**: This is comparing apples to oranges (different eras, populations, fever thresholds). We're inferring a biological parameter from confounded observational data. But I don't see an alternative.

**Sanity check**: The biological story (bicarbonate neutralizes acid → more bacteria survive) suggests $\delta > 1$ (bicarbonate is more potent). A ~10,000x offset means gastric acid kills 99.99% of bacteria in milk delivery. Is that plausible? Literature on gastric killing of Salmonella suggests pH-dependent log-reductions in that range, so... maybe?

### $\gamma$ (immunity scaling)

This is the hard one.

**Best data**:
- **Darton 2016**: Anti-Vi IgG vs attack rate. HR = 0.29 per log10 increase. This is the cleanest immunity-outcome relationship.
- **Gibani 2020**: Prior challenge vs rechallenge. ~36% relative risk reduction for S. Typhi. But confounded.

**Weaker data**:
- **Gilman 1977**: H antibody >= 1:20 → 61% vs 24% attack rate. N=27. Old assay.
- **Dupont 1971**: Prior infection → 23% vs 30% attack rate. Weak effect.
- **Woodward 1980**: Military (probably vaccinated) → 20% vs 48%. Good effect size but coarse immunity measure.

**My concern**: The best $\gamma$ data (Darton 2016) is at a single dose. We need to extrapolate the immunity-dose interaction. If immunity just shifts the curve left (like medium offset does), we can estimate $\gamma$ from single-dose immunity contrasts. But if immunity changes the curve shape, we need immunity variation across doses—and we don't have that.

**Identifiability question**: Can we actually estimate $\gamma$ separately from a baseline immunity random effect? I'm not sure. The model has $\text{CoP}^{-\gamma}$ scaling $\alpha$. If everyone has a different CoP (which they do—immunity varies), how do we separate "$\gamma$ is 1 and CoP variance is high" from "$\gamma$ is 2 and CoP variance is lower"?

We might need to fix $\gamma$ based on external information (e.g., assume it's 1 = linear immunity scaling) and estimate CoP distributions.

---

## The Proposed Working Model (Specific)

### Likelihood

For observation $i$ in study $s$:

$$
y_i \sim \text{Binomial}(n_i, p_i)
$$

$$
p_i = 1 - \left(1 + D_i \cdot \delta_s \cdot \frac{2^{1/\alpha} - 1}{N_{50}}\right)^{-\alpha \cdot \text{RR}_s}
$$

Where:
- $D_i$ = challenge dose for observation $i$
- $\delta_s$ = medium offset (1 for Maryland/milk, estimated for Oxford/bicarbonate)
- $\text{RR}_s$ = relative risk for immunity stratum $s$ (1 for naive, estimated for "immune")

**Simplifications made**:
1. Collapsed CoP to binary (naive vs immune)
2. Replaced $\text{CoP}^{-\gamma}$ with simple relative risk multiplier $\text{RR}$
3. No cascade (single outcome: typhoid fever)
4. Study-level effects absorbed into $\delta$ and $\text{RR}$

### Parameters to estimate

| Parameter | Prior | Rationale |
|-----------|-------|-----------|
| $N_{50}$ | Log-normal(7, 1) | ~10^7 CFU, order-of-magnitude uncertainty |
| $\alpha$ | Half-normal(0, 1) | Positive, weakly informative |
| $\delta_{\text{bicarb}}$ | Log-normal(3, 1) | ~1000-10000x offset expected |
| $\text{RR}_{\text{immune}}$ | Beta(2, 5)? | Expect 0.3-0.7 |

### Data to include

**Tier 1 (core dose-response)**:
- Hornick 1966/1970: Doses 10^5, 10^7, 10^8, 10^9 (exclude 10^3 outlier)
- Waddington 2014: Doses 10^3, 10^4, 10^5

**Tier 2 (single-dose validation)**:
- Levine 1976 controls: 10^5 milk, N=97
- Gilman 1977 controls: 10^5 milk, N=64
- Jin 2017 controls: 10^4 bicarb, N=31
- Darton 2016 placebo: 10^4 bicarb, N=30

**Tier 3 (immunity contrast)**:
- Darton 2016: Placebo (67%) vs Ty21a (43%) at same dose
- Gibani 2020: Naive (63%) vs rechallenge (44%) for S. Typhi
- Gilman 1977: H antibody low (61%) vs high (24%) in controls

---

## Things I Don't Know How to Do

### 1. Dealing with cohort overlap

Many papers share subjects:
- Hornick 1966, 1967, Appraisal, Snyder 1970, Woodward 1980, Glynn 1995 all draw from the same ~1700 Maryland volunteers.
- Waddington 2014, Darton 2012/2016/2017, Gibani 2019/2020 share Oxford cohorts.

**Problem**: If I naively pool them, I'm double-counting subjects.

**Possible solutions**:
- Pick one paper per cohort as the "canonical" source.
- Try to reconstruct non-overlapping subsets.
- Model the overlap explicitly (hierarchical structure with shared subjects).

I'd probably go with: Use Hornick 1970 as the Maryland canonical source (most complete), Waddington 2014 for Oxford dose-response, and carefully vet the vaccine trial controls for non-overlap.

### 2. The 10^3 milk anomaly

Hornick 1966 reports 64% attack rate at 10^3 (N=14). This is higher than 10^5 (28%, N=116). Either:
- Small sample fluke (most likely)
- Different subject population for that arm
- Dose measurement error
- Something biological we don't understand

**My call**: Exclude it from primary analysis, note as sensitivity.

### 3. Bayesian vs frequentist

The model has maybe 4-5 parameters and ~10-15 data points (dose-by-study combinations). That's a pretty weak data-to-parameter ratio. Bayesian inference with informative priors is probably necessary.

But how informative? If priors dominate, we're not really learning from data.

### 4. What is "immune" anyway?

The binary naive/immune split is crude. In reality:
- Some "naive" subjects have cross-reactive antibodies (Darton 2016: 29% with anti-Vi)
- Some "previously exposed" subjects never mounted a response
- Vaccine-induced immunity differs from infection-induced immunity

A latent immunity model (where observed antibodies are noisy observations of true CoP) would be cleaner but probably unidentifiable with this data.

### 5. Should we model bacteremia separately?

The Oxford data distinguishes "bacteremia only" from "fever" from "both". We could fit two curves:
- P(bacteremia | dose)
- P(fever | bacteremia)

This would let us separate "infection" from "disease". But:
- Maryland didn't always report bacteremia separately
- It adds parameters
- The curves will be highly correlated

Maybe worth trying as sensitivity.

---

## Summary: My YOLO Recommendation

**Minimal viable model**:
1. Single beta-Poisson dose-response to "typhoid fever"
2. Fixed medium offset (milk=1, estimate $\delta_{\text{bicarb}}$)
3. Binary immunity (naive vs immune, estimate relative risk)
4. Fit Hornick + Waddington for dose-response shape
5. Use vaccine trial controls and Darton/Gibani for immunity contrast

**Parameters**: $N_{50}$, $\alpha$, $\delta_{\text{bicarb}}$, $\text{RR}_{\text{immune}}$

**Biggest uncertainties**:
1. Medium offset is confounded with era/population/definition differences
2. Immunity effect is confounded with prior disease susceptibility (Gibani 2020)
3. Cohort overlap makes sample sizes unclear

**What to do next**: Actually extract the numbers, build the likelihood, fit it, see what breaks.

---

## Appendix: Paper-Parameter Mapping Table

| Parameter | Primary Papers | Supporting Papers | Confounders/Concerns |
|-----------|----------------|-------------------|---------------------|
| $N_{50}$ | Hornick 1966/1970 | Glynn 1995, Waddington 2014 | Era differences, medium |
| $\alpha$ | Hornick 1966/1970 | Waddington 2014 | Correlated with $N_{50}$ |
| $\delta_{\text{medium}}$ | Cross-era comparison | - | Everything is confounded |
| $\text{RR}_{\text{immune}}$ | Darton 2016, Gibani 2020 | Gilman 1977, Dupont 1971, Woodward 1980 | Definition of "immune", Gibani susceptibility finding |

---

## Appendix B: The Three-Outcome Model (Shedding, Seroconversion, Fever)

*Added after user push: "what if we want shedding, seroconversion, AND fever?"*

### Assumed Simplifications (Per User)

1. **Ignore cohort overlap** — accept double-counting, YOLO
2. **Prior on medium offset** — from biochemistry/animal models; attribute some era difference to Maryland's unknown baseline immunity
3. **The Gibani susceptibility paradox** — *hmmm* (we'll deal with this somehow)
4. **Gamma identifiability** — get clever with similar cohorts across studies
5. **10^3 outlier** — fine, exclude it
6. **Pool gamma across outcomes** — one immunity scaling parameter for all three

---

### The Three-Outcome Structure

The outcomes form a partial causal ordering:

```
                                    ┌─→ Stool Shedding
Dose → Colonization → Infection ────┤
                                    └─→ Systemic Invasion → Fever
                                              ↓
                                        Seroconversion
```

But we can't observe colonization or systemic invasion directly. So pragmatically:

- **Shedding** ≈ proxy for colonization (can happen without fever)
- **Fever** ≈ proxy for clinical disease (the downstream outcome)
- **Seroconversion** ≈ proxy for immune response to infection (requires some level of antigen exposure)

The key insight: **these outcomes have different thresholds on the same latent infection process**. Shedding happens first and at lower doses; fever requires more severe infection; seroconversion requires enough antigen exposure to trigger measurable antibody response.

---

### Proposed Model: Shared Shape, Outcome-Specific Thresholds

$$
P(\text{outcome}_k | D, \text{medium}, \text{CoP}) = 1 - \left(1 + D \cdot \delta_{\text{medium}} \cdot \frac{2^{1/\alpha} - 1}{N_{50,k}}\right)^{-\alpha / \text{CoP}^{\gamma}}
$$

Where:
- $N_{50,k}$ = outcome-specific ID50 (shedding < seroconversion < fever, presumably)
- $\alpha$ = **shared** heterogeneity (same underlying host/dose variability)
- $\delta_{\text{medium}}$ = **shared** medium offset
- $\gamma$ = **shared** immunity scaling (per user request)

**Parameter count**:
- 3 outcome-specific $N_{50,k}$ values
- 1 shared $\alpha$
- 1 shared $\delta_{\text{medium}}$
- 1 shared $\gamma$
- Total: 6 parameters

---

### Data Inventory by Outcome

#### Shedding

| Source | Dose | Medium | N | Shedding Rate | Notes |
|--------|------|--------|---|---------------|-------|
| Waddington 2014 | 10^3 | bicarb | 20 | 65% | Oxford dose-escalation |
| Waddington 2014 | 10^4 | bicarb | 16 | 63% | |
| Waddington 2014 | 10^5 | bicarb | 5 | 80% | |
| Gilman 1977 | 10^5 | milk | 43 | 49% (early), 60% (late) | Maryland controls |
| Gibani 2019 | 10^4 | bicarb | 158 | ~65%* | Pooled unvaccinated |
| Darton 2016 | 10^4 | bicarb | 30 | ~50%* | Placebo arm |

*Estimated from aggregate percentages; individual counts need verification.

**Observation**: Shedding rates are ~10-15% higher than fever rates at the same dose (Oxford). This is expected if shedding captures subclinical infection.

**Maryland gap**: Hornick 1966/1970 say "all stools were cultured" but I don't see dose-stratified shedding rates reported separately from disease. Gilman 1977 has the clearest data (49-60% in controls). **[OPEN]** Need to check if Hornick papers have dose-stratified shedding hidden somewhere.

#### Seroconversion

This is the messiest outcome. Problems:
- Multiple antigens (O, H, Vi)
- Different assays across eras
- Threshold-dependent (what counts as "conversion"?)
- Confounded with disease status (sick people seroconvert more)

| Source | Dose | Medium | N | Seroconversion Rate | Antigen | Notes |
|--------|------|--------|---|---------------------|---------|-------|
| Gilman 1977 | 10^5 | milk | 30 (with disease) | 70% O, 63% H | O, H | Maryland controls who got sick |
| Gilman 1977 | 10^5 | milk | 32 (no disease) | 6% O, 34% H | O, H | Maryland controls who didn't get sick |
| Darton 2016 | 10^4 | bicarb | 30 | varies by antigen | Vi, LPS, H | Oxford placebo |
| Hornick 1970 | various | milk | ? | "no correlation with protection" | O, H, Vi | But rates not tabulated by dose |

**Key finding (Gilman 1977)**: Seroconversion is highly conditional on disease status. Among those who got sick: 70% seroconverted (O). Among those who didn't: 6%. This means seroconversion is downstream of infection/disease, not just exposure.

**Implication**: We can't fit a simple dose→seroconversion curve. We'd need:
$$P(\text{seroconversion} | D) = P(\text{seroconversion} | \text{disease}) \cdot P(\text{disease} | D) + P(\text{seroconversion} | \text{no disease}) \cdot P(\text{no disease} | D)$$

Or equivalently, model seroconversion as conditional on infection/disease status rather than directly on dose.

**This is different from shedding and fever**, which can be modeled as direct dose-response curves.

#### Fever (Clinical Typhoid)

This is what we have the most data for (see main document).

| Source | Dose | Medium | N | Fever Rate | Notes |
|--------|------|--------|---|------------|-------|
| Hornick 1970 | 10^5 | milk | 104 | 27% | Maryland unvaccinated |
| Hornick 1970 | 10^7 | milk | 30 | 50% | |
| Hornick 1970 | 10^9 | milk | 4 | 100% | Small N |
| Waddington 2014 | 10^3 | bicarb | 20 | 55% | Oxford |
| Waddington 2014 | 10^4 | bicarb | 16 | 63% | |
| Waddington 2014 | 10^5 | bicarb | 5 | 100% | |
| Darton 2016 | 10^4 | bicarb | 30 | 67% | Placebo |
| Jin 2017 | 10^4 | bicarb | 31 | 77% | Control |
| Gilman 1977 | 10^5 | milk | 64 | 48% | Maryland |
| Levine 1976 | 10^5 | milk | 97 | ~40% | Maryland |

---

### What Would Break?

#### 1. Seroconversion isn't a simple dose-response

The Gilman 1977 data shows seroconversion is conditional on disease, not a direct function of dose. This breaks the "three parallel curves" model.

**Possible fix**: Model seroconversion as:
$$P(\text{seroconversion}) = P(\text{seroconversion} | \text{infected}) \cdot P(\text{infected})$$

Where P(infected) comes from the dose-response, and P(seroconversion | infected) is a separate parameter (~70% based on Gilman). But this adds complexity and we lose the elegant "shared gamma" structure.

**Alternative**: Redefine the outcome as "seroconversion OR disease" (i.e., any evidence of infection). But this is awkward because seroconversion is measured after the fact.

**My suspicion**: Seroconversion may not be worth modeling as a separate outcome. It's too downstream and too assay-dependent. Better to use it as a validation check (do our infection predictions match observed seroconversion rates?) rather than a calibration target.

#### 2. Maryland shedding data is sparse

We have Gilman 1977 (one dose, 10^5) and... that's kind of it for dose-stratified shedding in Maryland. The multi-dose curve comes entirely from Oxford.

**Implication**: We can't independently estimate the medium offset for shedding. We'd have to assume it's the same as for fever (which is reasonable biologically—same gastric transit—but unverified).

#### 3. The outcome ordering constraint

If shedding ≤ seroconversion ≤ fever in terms of threshold, we should have $N_{50,\text{shed}} < N_{50,\text{sero}} < N_{50,\text{fever}}$. But this isn't obviously true from the data:

- Waddington 2014 at 10^3: 65% shedding, 55% fever (shedding > fever ✓)
- Waddington 2014 at 10^4: 63% shedding, 63% fever (equal)
- Gilman 1977 controls at 10^5: 49-60% shedding, 48% fever (about equal)

The shedding and fever rates are remarkably similar. Maybe the cascade collapses faster than I thought—once you're colonized enough to shed, you're also likely to get sick.

**Implication**: Shedding and fever may not be separable. A single "infection" outcome might be sufficient.

#### 4. Pooling gamma across outcomes assumes immunity affects all outcomes equally

If anti-Vi IgG protects against fever but not shedding (or vice versa), the shared-gamma assumption breaks.

**Evidence**: Gibani 2019 shows anti-Vi IgG protects against shedding (P < .0001) AND Darton 2016 shows it protects against fever (HR = 0.29/log). So the direction is consistent. But the *magnitude* may differ.

**Test**: Compare vaccine efficacy against shedding vs fever. From Gibani 2019:
- Vi-PS vs shedding: OR 0.34 (~66% reduction)
- Vi-TT vs shedding: OR 0.41 (~59% reduction)

From Darton 2016:
- Ty21a vs fever: VE 35%
- M01ZH09 vs fever: VE 13%

These aren't directly comparable (different vaccines, different studies), but the magnitudes are in the same ballpark. Shared gamma seems defensible.

---

### Revised Three-Outcome Proposal

Given the above, here's what I'd actually try:

**Model Structure**:
1. **Fever**: Full dose-response with $N_{50,\text{fever}}$, $\alpha$, $\delta$, $\gamma$
2. **Shedding**: Same model with outcome-specific $N_{50,\text{shed}}$, shared $\alpha$, $\delta$, $\gamma$
3. **Seroconversion**: Don't model directly. Use as validation.

**Constraint**: $N_{50,\text{shed}} \leq N_{50,\text{fever}}$ (shedding threshold ≤ fever threshold)

**Parameters**: 5 total (vs. 6 for full three-outcome)
- $N_{50,\text{fever}}$
- $N_{50,\text{shed}}$ (or: ratio $N_{50,\text{shed}}/N_{50,\text{fever}}$)
- $\alpha$ (shared)
- $\delta$ (shared)
- $\gamma$ (shared)

---

### Paper-Parameter Mapping (Three-Outcome Version)

| Parameter | Outcome | Primary Papers | Supporting | Concerns |
|-----------|---------|----------------|------------|----------|
| $N_{50,\text{fever}}$ | Fever | Hornick 1970, Waddington 2014 | Glynn 1995, Gilman 1977, Levine 1976 | Era confounding |
| $N_{50,\text{shed}}$ | Shedding | Waddington 2014, Gibani 2019 | Gilman 1977 (one dose only) | Maryland gap |
| $\alpha$ | Both | Hornick 1970 (fever), Waddington 2014 (both) | — | Correlated with N50s |
| $\delta$ | Both | Cross-era comparison | — | Mega-confounded |
| $\gamma$ | Both | Darton 2016 (fever), Gibani 2019 (shedding) | Gilman 1977, Gibani 2020 | Assumes equal effect |

---

### What About the Gibani Susceptibility Paradox? (Item 3: "hmmm")

The finding that people who *didn't* get sick on first challenge are *more* protected on rechallenge suggests innate susceptibility heterogeneity. This is a problem for any model that assumes:
- Prior exposure → immunity → protection

Instead it suggests:
- Some people are innately resistant
- Resistant people don't get sick on first challenge AND remain resistant on rechallenge
- Susceptible people get sick on first challenge AND remain susceptible on rechallenge

**Implications for the model**:

1. **Host heterogeneity ($\alpha$) may be doing double duty**: capturing both dose variability AND stable susceptibility differences.

2. **Prior disease is a selection marker, not a treatment**: People with prior disease aren't "trained" — they were susceptible all along. People without prior disease aren't "naive who got lucky" — they may be innately resistant.

3. **This could break the immunity model**: If "prior exposure" protection is actually "resistant phenotype revelation", then gamma isn't measuring immunity scaling — it's measuring selection.

**Possible approaches**:

a) **Ignore it** (purest YOLO). Assume traditional immunity model is approximately correct. Accept that gamma is an effective parameter mixing true immunity with susceptibility selection.

b) **Add a susceptibility class**. Model subjects as mixture: fraction $\pi$ are "susceptible" (get the normal dose-response), fraction $1-\pi$ are "resistant" (much flatter dose-response or lower max attack rate). Prior disease status reveals which class you're in.

c) **Reframe prior exposure effects**. Don't estimate "immunity from exposure" directly. Instead, note that cohorts with prior exposure have enriched resistant phenotypes, and adjust expectations accordingly.

**What I'd try first**: Option (a). Fit the standard model, see if it works. If posterior predictive checks show systematic misfit for rechallenge cohorts, revisit with option (b).

The key test: Does the model predict that prior-disease subjects have *higher* attack rates on rechallenge than prior-no-disease subjects (as Gibani found)? If the standard model can't accommodate this, we need susceptibility classes.

---

### Summary: Three-Outcome YOLO Plan

1. **Fit fever + shedding** as parallel dose-response curves with shared $\alpha$, $\delta$, $\gamma$ and outcome-specific $N_{50}$

2. **Use seroconversion for validation** rather than calibration (too conditional on disease status)

3. **Accept** that gamma pools immunity effects with susceptibility selection (revisit if it breaks)

4. **Prior on $\delta$** from biochemistry/animal models + attribute some era difference to Maryland baseline immunity (reduces effective $\delta$ from naive comparison)

5. **Constrain** $N_{50,\text{shed}} \leq N_{50,\text{fever}}$

6. **Test** with posterior predictive checks whether model can handle:
   - Shedding > fever at low doses
   - Rechallenge dynamics (including Gibani paradox)

7. **If it breaks**: add susceptibility classes or separate immunity parameters by outcome

---

## Appendix C: Comparison with `cohort_incidence_model_proof_of_concept.R`

*Added after reviewing Mike's existing prototype implementation.*

### What the Prototype Does

The first 200 lines of `cohort_incidence_model_proof_of_concept.R` implement:

1. **Titer dynamics model** (`titer_vs_time`, lines 59-77): Power-law decay of anti-Vi IgG titers after immunizing events, with age-dependent decay parameters.

2. **Fold-rise model** (`fold_rise_model`, lines 91-96): How CoP changes after an immunizing event, as a function of pre-challenge CoP.

3. **Dose-response model** (`p_outcome_given_dose`, lines 135-161): Beta-Poisson model for two outcomes—infection and fever—as a function of dose and CoP.

4. **Protective efficacy** (`protective_efficacy`, lines 183-187): Relative risk reduction for immune vs naive individuals.

---

### How the Prototype Reflects What I Noticed

#### The dose-response functional form is the same

The prototype uses the beta-Poisson form from the model specification:

```r
p = 1 - (1+dose*(2^(1/alpha)-1)/n50)^(-alpha/(CoP_pre^gamma))
```

(This is expected—the model specification I read today was written based on this prototype, so the functional form match is by construction.)

#### Two outcomes, not three

The prototype models:
- `infection_given_dose` — analogous to my "shedding" outcome
- `fever_given_dose` — analogous to my "fever" outcome
- `fever_given_infection` — derived as the ratio

No explicit seroconversion outcome, consistent with my recommendation to use seroconversion for validation rather than calibration.

#### CoP is a continuous scalar (anti-Vi IgG)

The prototype treats CoP as a continuous variable (anti-Vi IgG titer in EU/mL), not the binary naive/immune indicator I proposed for simplicity. This is more realistic but requires the titer dynamics machinery.

Baseline CoP = 1 (below assay LOD), maximum CoP ~ 10^3.5. The model is parameterized in terms of this titer scale.

---

### Where the Prototype Differs from My YOLO Proposal

#### Outcome-specific parameters vs pooled

I proposed pooling $\gamma$ across outcomes to reduce parameters. The prototype takes a different approach—outcome-specific values for all three parameters:

| Parameter | Infection | Fever | Relationship |
|-----------|-----------|-------|--------------|
| $N_{50}$  | 2,780 | 27,800 | `n50_infection = n50_fever / 10` |
| $\alpha$  | 1.68 | 0.84 | `alpha_infection = alpha_fever * 2` |
| $\gamma$  | 0.2 | 0.4 | `gamma_infection = gamma_fever / 2` |

**Implications**:
- Infection has a ~10x lower threshold than fever (makes sense—you can be infected without fever)
- Infection has more heterogeneity ($\alpha$ = 1.68 vs 0.84)
- Immunity has a **weaker effect** on infection than fever ($\gamma$ = 0.2 vs 0.4)

That last point is notable: the prototype encodes that immunity reduces fever more than it reduces infection. This is consistent with the "leaky protection" concept—immunity may not prevent colonization but does prevent progression to clinical disease.

**Two reasonable approaches**: My pooling proposal reduces parameters at the cost of flexibility. The prototype's approach is more flexible and can capture differential effects of immunity on infection vs disease. Which is better depends on what the data can support.

#### No explicit medium offset

The prototype doesn't have a `delta_medium` parameter. The N50 values (2,780 and 27,800 CFU) appear to be pre-calibrated to a reference condition.

For context: my YOLO estimated N50 ~ 10^7 for fever in milk delivery. The prototype's N50 = 27,800 is ~360x lower. This is consistent with a ~2.5 log medium offset baked into the parameters—i.e., these N50 values are probably in "bicarbonate-equivalent" CFU.

**This means**: The medium offset I was worried about may already be absorbed into the prototype's parameter values. If we want to fit Maryland milk data directly, we'd need to either (a) add the offset back, or (b) convert Maryland doses to bicarbonate-equivalents before fitting.

#### Continuous titer dynamics vs static immunity

The prototype has a full time-dynamic model for how titers evolve:
- Rise phase (logistic rise over ~30 days)
- Decay phase (power-law decay with age-dependent exponent)
- Fold-rise model for boosting from subsequent exposures

My YOLO treated immunity as static (binary or categorical at time of challenge). The prototype is much richer.

**Implication**: For calibration to challenge study data (single time point), my static approach is probably fine. But for downstream use in a transmission model, the prototype's dynamics are necessary.

---

### Open Questions in Both Approaches

#### The Gibani susceptibility paradox

Not addressed in the prototype. The model assumes everyone has the same dose-response conditional on their CoP—no innate susceptibility classes.

If the Gibani finding (non-sick-on-first-challenge people are more protected on rechallenge) is real, this model can't capture it. The prototype would predict that people who didn't get sick had higher CoP, and that higher CoP explains their subsequent protection. But Gibani found baseline antibodies were NOT different between groups.

**Potential fix**: Add a latent susceptibility class (fraction $\pi$ are "high susceptibility", $1-\pi$ are "low susceptibility") with class-specific dose-response parameters.

#### Era/study heterogeneity

The prototype has fixed parameters, no study-level random effects. For use in a transmission model (single simulated world), this is fine. For calibration to heterogeneous challenge study data across eras, we'd need to add study-level effects.

#### Measurement model for outcomes

The prototype treats outcomes as binary (infected or not, fever or not). No explicit measurement error or threshold-dependent observation model.

For challenge studies where outcomes are measured intensively (daily cultures, continuous temperature), this is probably fine. For field data with imperfect ascertainment, we'd need more machinery.

---

### How Consistent is My Approach with the Prototype?

| Aspect | My YOLO | Prototype | Comparison |
|--------|---------|-----------|------------|
| Functional form | Beta-Poisson | Beta-Poisson | Same (by construction) |
| Outcomes modeled | Shedding, Fever | Infection, Fever | Same concept (infection ≈ shedding) |
| Seroconversion | Validation only | Not modeled | Aligned |
| CoP representation | Binary/categorical | Continuous titer | Prototype more detailed |
| $N_{50}$ | Outcome-specific | Outcome-specific | Same approach |
| $\alpha$ | Shared | Outcome-specific | Different choices |
| $\gamma$ | Shared (proposed) | Outcome-specific | Different choices |
| Medium offset | Explicit parameter | Absorbed into N50 | Different approaches |
| Titer dynamics | Static | Time-varying | Prototype more detailed |
| Susceptibility classes | Flagged as potential need | Not included | Open question for both |

**Bottom line**: The functional form and overall structure align. The key difference is that **the prototype has outcome-specific $\alpha$ and $\gamma$**, while I proposed pooling them. These represent different tradeoffs between parsimony and flexibility—what's right depends on what the data can identify.

---

### Synthesis: Options Going Forward

Based on comparing the two approaches, here are the decision points for the calibration work:

1. **Pooled vs outcome-specific $\gamma$**: My YOLO proposed pooling; the prototype separates them. The prototype's `gamma_fever = 2 * gamma_infection` encodes "leaky protection" (immunity affects disease more than infection). This is a modeling choice—pooling is more parsimonious, separating is more flexible.

2. **Pooled vs outcome-specific $\alpha$**: Same tradeoff. The prototype has `alpha_infection = 2 * alpha_fever`, encoding more heterogeneity in infection outcomes. Could also use a ratio constraint as a middle ground.

3. **Explicit medium offset vs absorbed into N50**: My YOLO proposed an explicit $\delta_{\text{medium}}$ parameter. The prototype absorbs the offset into the N50 values (calibrated to bicarbonate-equivalent doses). Both work—explicit offset is clearer for interpretation, absorbed is simpler for the likelihood.

4. **Static vs dynamic immunity**: For calibration to challenge study data (single time point), the static approach is sufficient. The prototype's titer dynamics become important for transmission model applications.

5. **Susceptibility classes**: Neither approach addresses the Gibani finding. This remains an open modeling question that may require a structural extension.

---

### Parameter Translation Table

For reference, mapping between my YOLO notation and the prototype:

| My YOLO | Prototype Variable | Prototype Default |
|---------|-------------------|-------------------|
| $N_{50,\text{fever}}$ | `n50_fever_given_dose` | 27,800 |
| $N_{50,\text{shed}}$ | `n50_infection_given_dose` | 2,780 |
| $\alpha$ (shared) | — | — |
| $\alpha_{\text{fever}}$ | `alpha_fever_given_dose` | 0.84 |
| $\alpha_{\text{infection}}$ | `alpha_infection_given_dose` | 1.68 |
| $\gamma$ (shared) | — | — |
| $\gamma_{\text{fever}}$ | `gamma_fever_given_dose` | 0.4 |
| $\gamma_{\text{infection}}$ | `gamma_infection_given_dose` | 0.2 |
| $\delta_{\text{medium}}$ | (absorbed into N50) | — |
| CoP | `CoP_pre` | 1 (naive baseline) |

---

## Appendix D: Is Anything Hopeless?

*Reflection on whether the prototype's implicit choices have fatal identifiability problems.*

**The parameter structure isn't as heavy as it first appears**

Looking again at the prototype, the "6 parameters" (N50, alpha, gamma for each of two outcomes) are actually 3 free parameters + 3 fixed ratios:

```r
n50_infection = n50_fever / 10
alpha_infection = alpha_fever * 2
gamma_infection = gamma_fever / 2
```

So you're really estimating `n50_fever`, `alpha_fever`, `gamma_fever`, and the ratios encode biological assumptions (infection is easier than fever, more heterogeneous, and less affected by immunity). That's quite parsimonious.

**The hard part: separating what from what**

The challenge I see is the cross-era identification problem:

| Era | What it gives us | What's missing |
|-----|------------------|----------------|
| Maryland | Multi-dose data (10^5 to 10^9) → can identify N50 and alpha jointly | No CoP measurements → can't identify gamma |
| Oxford | CoP variation (anti-Vi titers) → can identify gamma | Narrow dose range (10^3-10^4) → N50 and alpha less constrained |

The implicit assumption bridging them: **alpha is the same across eras**, and all the era difference goes into either medium offset or baseline immunity differences.

Is that assumption defensible? Maybe. The heterogeneity parameter alpha captures host variability in the dose-response. If that's mostly biological (genetics, microbiome, innate immunity), it shouldn't change much between healthy US adults in the 1960s and healthy UK adults in the 2010s.

But if alpha also absorbs protocol differences (fasting, timing of cultures, fever ascertainment), it might differ.

**What doesn't seem hopeless**

1. **The functional form**: Beta-Poisson is well-motivated for single-hit infection processes. Even if the parameters are uncertain, the shape family is right.

2. **The two-outcome structure**: Infection and fever as separate dose-response curves with related parameters—this makes biological sense and the prototype's ratios (10x threshold difference, 2x heterogeneity difference, 0.5x immunity effect difference) are reasonable starting points.

3. **Using continuous CoP instead of binary immunity**: This is actually *better* for identifiability than my binary proposal. With continuous CoP, you can estimate gamma from the slope of attack rate vs titer.

**What's genuinely uncertain (but not hopeless)**

1. **The medium offset / baseline immunity confound**: Your assumption (2) in the YOLO—that some era difference is due to Maryland's unknown baseline immunity—is probably necessary. A pure medium offset of ~4 logs seems too large to be gastric acid alone. But how to partition it is genuinely uncertain.

2. **Whether the Gibani paradox matters**: If susceptibility classes are real and large, the current model structure will give biased gamma estimates (confusing selection with immunity). But it might still *fit* the data acceptably, just with a mechanistically muddled interpretation.

**One thing that gives me pause**

The prototype's N50 = 27,800 for fever implies an ID50 of ~28,000 CFU in bicarbonate. My back-of-envelope from Maryland milk data suggested ID50 ~ 10^7 in milk. That's a ~360x difference, or about 2.5 logs.

But if bicarbonate really increases effective dose by 2.5-4 logs (which the Waddington vs Hornick comparison suggests), then the prototype's N50 is in the right ballpark for bicarbonate-equivalent doses.

The question is whether that offset is stable. If gastric acid killing varies a lot between individuals or conditions, the "effective dose" framing might be leaky.

**Bottom line**

Nothing seems *hopeless*. The hardest part is the cross-era bridge—you need Maryland for the dose-response shape and Oxford for the immunity effect, and connecting them requires assumptions about what's shared vs what differs. Your prototype's structure handles this by absorbing the offset into N50 (calibrated to one reference condition) and using fixed ratios to constrain the outcome-specific parameters. That's a reasonable set of choices given the data constraints.
