# Typhoid dose-response from human challenge: what the Tier 2 calibration says

**Source:** `t2-indiv-vax`, stage `phi-rho-eta-psi-vax`, 182 observations. 1
divergent transition in 4,000 draws, max R-hat 1.003. All intervals are 90%
posterior credible intervals, median first. Everything is conditional on this
model and this corpus.

**This document is the interpretation.** It selects the results that carry a
scientific or decision consequence and says what we think they mean. It is not the
complete output of the fit. For that, read the generated
[Tier 2 posterior summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md):
every fitted parameter with its interval and sampler diagnostics, the pairwise
correlations above |r| = 0.7, the full priorsense power-scaling table, and an index
of all 33 figures — including the ones not discussed here (`calibration_targets`,
`delta_bridge`, `dose_cop_surface`, `maryland_mixture`, `titre_protection`, the
per-cohort `grid/` panels, and the sampler-diagnostic set). That file is written by
`diagnose_fit()` and is authoritative on numbers; where it and this document
disagree, it wins. The corresponding
[prior-predictive summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax-prior/summary.md)
is the matching prior run.

## 1. What the model is fitted to

One dose-response model, fitted jointly to seven Quailes-strain challenge studies
from two eras.

| Era | Studies | Vehicle | Dose range | Serology |
|---|---|---|---|---|
| Maryland, 1960s–70s | Hornick, Levine, Gilman | milk | 10³–10⁹ CFU | none |
| Oxford, 2010s | Waddington, Darton, Jin, Gibani | bicarbonate | 10³–2.5×10⁴ CFU | per-subject anti-Vi IgG |

The biology is a cascade: infection first, then fever given infection. Each step
is a modified beta-Poisson curve.

    P = 1 − (1 + D·(2^(1/α) − 1)/N50)^(−α / (CoP^γ · V))

`D` is bicarbonate-equivalent dose. `CoP` is anti-Vi IgG relative to naive. `V` is
a per-vaccine factor for the protection that the subject's own anti-Vi titer does
not explain; `V = 1` for placebo.

Three terms describe the observation process, not the biology. These three make
the studies combinable:

- **φ(T, D)** — fever-threshold sensitivity. The probability that a typhoid case
  crosses the study's temperature threshold `T`.
- **ψ** — infection-definition sensitivity, relative to blood-or-stool culture.
- **η(D)** — shedding-detection probability. The probability that Oxford stool
  shedding is seen before antibiotics truncate it.

## 2. The dose axis

Milk and bicarbonate are not the same dose. The bridge parameter δ = 101 [31–384]:
one bacterium delivered in bicarbonate is worth about a hundred delivered in milk.

δ trades off against both N50 values (r = −0.61 against `log10_N50_inf`, −0.69
against `log10_N50_fevginf`). The trade-off is moderate: each marginal contracts
about 2× from its prior, the milk-frame product N50×δ contracts 3.3×, and the
split between them still contracts 1.7×. Both directions are informed. The caution
is interpretive rather than statistical — δ absorbs whatever the milk did,
including anything else that differs between the eras, so 101× is a fitted vehicle
scalar and not a measurement of gastric survival.

δ is also what makes the rest of this document possible. Without it, the 1960s
milk ladder and the modern bicarbonate trials cannot share a dose axis, and the
modern trials alone span barely one decade of dose. The wide dose range in every
curve below is borrowed from Maryland.

On the milk scale used in the figures:

- N50 for infection: 4,200 CFU milk [800–13,000], which is 40 CFU
  bicarbonate-equivalent [5.8–184].
- N50 for fever given infection: 14,000 CFU milk [4,300–38,000], which is 139 CFU
  bicarbonate-equivalent [28–519].

Both curves are shallow. α_inf = 0.21 [0.15–0.31] and α_fev|inf = 0.27
[0.18–0.40]. In the beta-Poisson, small α means wide person-to-person variation in
susceptibility. The curve climbs slowly across four orders of magnitude, so low
doses are not safe and high doses do not guarantee disease.

## 3. Infection and fever against dose

![dose_response_fit](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/dose_response_fit.png)

Naive subject (CoP = 1), milk scale:

| milk dose | P(infection) | P(fever \| infection) | P(fever), latent | P(diagnosed at ≥39.4 °C) |
|---|---|---|---|---|
| 10³ | 0.34 [0.16–0.51] | 0.15 [0.05–0.33] | 0.05 [0.01–0.15] | 0.01 [0.00–0.05] |
| 10⁴ | 0.58 [0.47–0.68] | 0.46 [0.30–0.58] | 0.26 [0.16–0.37] | 0.11 [0.05–0.20] |
| 10⁵ | 0.74 [0.68–0.81] | 0.70 [0.61–0.77] | 0.52 [0.44–0.60] | 0.32 [0.25–0.41] |
| 10⁷ | 0.90 [0.85–0.95] | 0.91 [0.85–0.96] | 0.82 [0.76–0.88] | 0.71 [0.62–0.80] |
| 10⁹ | 0.96 [0.93–0.99] | 0.97 [0.94–0.99] | 0.94 [0.89–0.97] | 0.89 [0.82–0.95] |

The latent fever curve is the biology. The diagnosed curve is what a 1960s
Maryland protocol would have written down at Hornick's ≥39.4 °C threshold. At 10⁴
CFU milk the two differ by a factor of two.

Infection and fever separate most strongly at low dose. At 10³ CFU milk the model
puts a third of subjects into infection and one in twenty into fever. The infection
curve saturates about two logs earlier than the fever curve.

## 4. Anti-Vi titer changes progression more than acquisition

![cop_response_milk_doses](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/cop_response_milk_doses.png)

The immunity exponents differ between the two cascade steps:

- γ_inf = 0.081 [0.037–0.151]
- γ_fev|inf = 0.192 [0.089–0.307]

P(γ_fev|inf > γ_inf) = 0.91. Anti-Vi does comparatively little to stop infection.
It does more to stop an infection from becoming fever. That is the shape of the
figure: the green P_inf line is nearly flat across four decades of titer, and the
purple P_fev|inf line falls steeply.

At 10⁴ CFU milk:

| anti-Vi (CoP) | P(infection) | P(fever \| infection) | P(fever), latent |
|---|---|---|---|
| 1 (naive) | 0.58 | 0.46 | 0.26 |
| 10 | 0.51 | 0.32 | 0.16 |
| 100 | 0.44 | 0.22 | 0.10 |

**This is the result with the clearest decision consequence.** If anti-Vi protects
against disease more than against acquisition, then Vi conjugate vaccines will cut
clinical typhoid faster than they cut transmission. Disease endpoints and
transmission endpoints then need different vaccine-efficacy inputs, and a single
efficacy number applied to both will overstate the transmission benefit.

Two things qualify it.

**The separation is not stable across the tier ladder.** At `t1-indiv` the two
exponents are indistinguishable. The split emerges when the vaccine arms enter and
widens at Tier 2:

| configuration | individual rows | γ_inf | γ_fev\|inf |
|---|---|---|---|
| `t1-indiv` | 56 | 0.166 [0.068–0.303] | 0.163 [0.066–0.298] |
| `t1-indiv-vax` | 153 | 0.111 [0.046–0.215] | 0.164 [0.069–0.284] |
| `t2-indiv-vax` | 153 | 0.081 [0.037–0.151] | 0.192 [0.089–0.307] |

**The overall slope is shallower than any direct contrast in the corpus** — see
§5. That affects the level of both exponents, though not obviously their ordering.

## 5. How steep is titer → protection? The corpus disagrees

Three sources speak to this, and they do not agree. All figures below are
converted to the same scale: percentage reduction in cumulative risk per ten-fold
rise in anti-Vi IgG.

| source | endpoint | risk reduction per 10× titer |
|---|---|---|
| this model, posterior | fever | 22% [14–31%] |
| Jin ViTT arm, raw (CoP 2.2 → 152) | fever | 34% [20–48%] |
| Jin ViPS arm, raw (CoP 2.2 → 38) | fever | 44% [25–60%] |
| Darton Table 7, hazard model | typhoid diagnosis | 59% [20–81%] |

Darton's paper states this as "a 1-log increase in anti-Vi IgG resulted in a 71%
decrease in hazard ratio for typhoid diagnosis [95% CI 30–88%]". Two conversions
sit between that sentence and the table above.

**First, hazard is not risk.** Their 0.29 is a hazard ratio on time to diagnosis
inside a 14-day challenge window. Ours is a ratio of cumulative attack rates at
the end of that window. At Darton's placebo attack rate of 20/30, a hazard ratio
of 0.29 implies a cumulative risk ratio of 0.41 — a 59% risk reduction, not 71%.
Roughly a fifth of the apparent gap is this conversion alone.

**Second, the endpoints differ.** Darton's endpoint is typhoid diagnosis, a
composite of clinical and microbiological criteria. Ours is fever. The model's
corresponding infection-only slope is much shallower still: 6% [3–12%] per 10×.

After both corrections the model is at 22% and Darton is at 59%, with intervals
that overlap only at the edges. That is a smaller disagreement than 22-versus-71
but it is a real one.

**Why our number is the low one.** The model is shallower than *every* direct
contrast, including both Jin arms, and the reason is likelihood weight rather than
evidence. Individualizing Darton's three arms puts 153 per-subject rows into the
likelihood. Those subjects' anti-Vi titers are mostly at or near the naive floor —
a thin range with a weak association. Jin's Vi-vaccine arms, which carry nearly all
the high-titer information, are three grouped rows. γ_inf falls monotonically as
the individual rows accumulate (table in §4).

It shows in the fit. The model compresses Jin's contrast:

| Jin arm | CoP | observed fever | fitted |
|---|---|---|---|
| control | 2.2 | 24/31 = 0.77 | 0.69 [0.63–0.75] |
| ViPS | 38 | 13/35 = 0.37 | 0.48 [0.39–0.57] |
| ViTT | 152 | 13/37 = 0.35 | 0.38 [0.27–0.50] |

We do not read Darton's hazard ratio as contradicting the model. We do read the
tier ladder as saying γ_inf is set more by row count than by titer range, and that
is a consequence of the individualization design rather than of the data. The
check is a sensitivity fit with the Darton individual rows grouped or
down-weighted: if γ_inf recovers toward 0.15, the effect is weighting.
`[open, not run]`

## 6. Fever is a threshold, not a state

![phi_severity](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/phi_severity.png)

Every study in the corpus defines "typhoid fever" at a different temperature.
φ(T, D) is the map between them, pinned by Darton's own threshold ladder (16, 10,
and 8 of 20 typhoid diagnoses crossing 38.0, 38.5, and 39.0 °C).

At low dose, the definition dominates:

- φ₀(38.0 °C) = 0.77 [0.65–0.88]
- φ₀(38.3 °C) = 0.67 [0.56–0.77] — Levine and Gilman
- φ₀(39.4 °C) = 0.21 [0.09–0.40] — Hornick

Hornick's endpoint is roughly a third as sensitive as Levine's at the same
underlying disease. Correct for the threshold before comparing historical attack
rates across these studies.

Dose then lifts the whole map. φ(39.4 °C) rises from 0.26 [0.13–0.44] at 10³ CFU
milk to 0.86 [0.81–0.91] at 10⁷ and 0.95 [0.92–0.98] at 10⁹. **Severity rises with
dose**, and the model gets this without a separate severity parameter: φ is a
survival function in temperature, so −dφ/dT is the implied peak-temperature
distribution among diagnosed cases (panel D).

**Caveat on panel D.** The dose-lift exponent `beta_phi` is pinned to 1, so the
lift has no temperature dependence. The implied severity density therefore
collapses at high dose rather than shifting upward, and φ(41 °C) = 0.90 at 10⁶
CFU. The direction of the severity effect is real. The high-dose shape is an
artifact of the pin. `[known limitation, documented in the figure]`

## 7. Infection is a definition too — treatment, shedding, and case criteria

![eta_detection](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/eta_detection.png)

Oxford treats at diagnosis. Shedding that would have started after treatment never
appears in the record. η(D) carries that loss.

Fitted: η_lo = 0.71 [0.61–0.81], κ = 0.98 [0.17–5.09]. Across the actual Oxford
dose range η is flat — 0.716 at 10³ CFU, 0.711 at 10⁴ CFU. **About 29% of true
Oxford infections are invisible in the stool record.**

Read that as a constant detection loss, not a dose effect. `tier2_plan.md` stated a
falsifiable read: if κ came back prior-dominated, that would be the cross-study
protocol difference showing up in the fit rather than physiology. It did —
priorsense puts κ's prior sensitivity at 0.63 against a likelihood sensitivity of
0.05. The modeled decline happens below about 100 CFU, under every observed
Oxford dose.

It also does not explain the one dose trend it was reached for. Waddington's
shedding counts *fall* with dose, 13/20 at 10³ and 8/20 at 10⁴. With η flat and
P_inf rising, the model predicts the opposite — fitted 0.53 [0.45–0.61] and 0.60
[0.52–0.68] — and absorbs the discrepancy as sampling noise. At n=20 per arm it
can: both observations sit inside their posterior-predictive intervals.
`[the model does not reproduce this trend; it declines to]`

The Maryland studies used different infection criteria, and those are not
equivalent either:

- Hornick, stool or blood culture: ψ = 1 by definition
- Levine, any-time stool: ψ = 0.82 [0.71–0.92]
- Gilman, late stool only (days 4–30): ψ = 0.71 [0.54–0.84]

Levine's stool-only criterion misses about one infection in five that Hornick's
would have caught.

## 8. Ty21a and M01ZH09 protect beyond anti-Vi

Ty21a is Vi-negative. M01ZH09 did not raise anti-Vi IgG. Neither vaccine's
protection can run through the anti-Vi correlate. The model therefore gives each
vaccine a residual factor `V` on the same exponent as immunity, after accounting
for each subject's own measured titer.

- V_Ty21a = 1.88 [1.35–2.68]
- V_M01ZH09 = 1.44 [1.06–1.96]
- P(V_Ty21a > V_M01ZH09) = 0.86

At the Darton challenge dose of 18,200 CFU, holding anti-Vi fixed:

| | efficacy against infection | efficacy against fever |
|---|---|---|
| Ty21a | 0.24 [0.11–0.38] | 0.43 [0.21–0.62] |
| M01ZH09 | 0.13 [0.02–0.25] | 0.24 [0.04–0.44] |

These are close to what Darton reported after adjusting for baseline anti-Vi
titer: Ty21a 48% [4–72] against fever ≥38.0 °C and 41% [2–64] against bacteremia;
M01ZH09 19% [−27–48] and 28% [−7–52]. The cascade model reaches these numbers
without fitting a per-arm regression.

**Two caveats.** First, independence from anti-Vi is structural, not tested: `V`
multiplies `CoP^γ` in the exponent, so the model cannot detect an interaction it
has no term for. Second, `V_Ty21a` rests on one open-label arm of 29 subjects. The
historical Gilman and Hornick vaccine arms are in the corpus but are not in this
likelihood.

## 9. How well it fits

![ppc](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/ppc.png)

All 29 grouped observations fall inside their 90% posterior-predictive intervals.
Note that the error bars in the figure are credible intervals on the *fitted
probability*, not predictive intervals for the count, so points sitting outside
their bars are expected at these sample sizes and are not misfits.

The shared beta-binomial overdispersion is small: ICC ρ = 0.024 [0.005–0.055].
After the observation-process terms are in, these studies are more consistent with
each other than the raw attack rates suggest.

Two structural misses are directional rather than large: the Waddington shedding
pair (§7), where the model's dose trend runs opposite to the data, and the Jin
titer contrast (§5), which the model compresses.

## 10. The observation-process terms are a result, not plumbing

φ, ψ, and η were built to make seven studies commensurable. They came back large
enough to be findings in their own right.

Each one is a factor the raw literature does not report. The fever-threshold map
spans a factor of three: at low dose, 67% of typhoid cases cross Levine's 38.3 °C
and only 21% cross Hornick's 39.4 °C, so the same underlying disease produces
attack rates that differ threefold by definition alone. The infection-definition
map costs about one infection in five: Levine's any-time-stool criterion has
sensitivity 0.82 [0.71–0.92] against Hornick's stool-or-blood culture, and
Gilman's late-shedding window 0.71 [0.54–0.84]. The treatment-truncation
correction costs three in ten: about 29% of true Oxford infections never appear in
the stool record because antibiotics start before shedding would have.

Those three numbers multiply. A Hornick-definition fever count and a Levine
stool-positive count, taken from studies of the same organism at the same dose,
are measuring quantities that differ by more than the dose effects people usually
argue about. Any meta-analysis that treats these endpoints as interchangeable is
pooling different quantities, and any transmission model that takes a historical
attack rate at face value inherits whichever definition that study happened to use.

The corollary matters for study design as much as for synthesis. Reporting the
threshold, the culture criterion, and the treatment protocol alongside an attack
rate is not bookkeeping — it is the difference between a number that can be pooled
and one that cannot. The modern Oxford trials do report all three, which is why
they can anchor φ and η at all. The Maryland studies mostly do not, and every
observation-process parameter here is an attempt to reconstruct from the outside
what those protocols would have recorded.

## 11. What is weakly identified

Do not quote these as measurements.

- **δ and the N50 values trade off** (r ≈ −0.6 to −0.7) but are not jointly
  unidentified. Marginals contract ~2× from prior; the sloppy direction still
  contracts 1.7×. δ is the weakest of the three on priorsense and the one most
  worth a sensitivity check, but it is not prior-carried.
- **κ** (η's dose dependence) is prior-dominated: priorsense prior sensitivity
  0.63, likelihood 0.05.
- **The Maryland immunity mixture.** π_susc = 0.67 [0.42–0.86] and CoP_imm = 9.2
  [1.0–37] are latent and prior-dominated. A third of the 1960s Maryland
  volunteers behave as if they carried roughly a nine-fold anti-Vi-equivalent
  titer — plausible for the period, but inferred from attack rates, not measured.
- **ψ_stool and η are confounded by construction.** The Darton cross-tab (15 of
  26) pins their product. Only η's dose dependence separates them, and η's dose
  dependence is not identified.
- **γ_inf may be weighting rather than evidence** (§5). Open.
- **`VE_fev_ViTT` and `VE_fev_ViPS`** in the generated quantities use placeholder
  CoP values. They are wiring, not results.
- One divergent transition remains unexplained.
- `frac_late` was expected to be prior-dominated and is not — priorsense
  likelihood sensitivity 0.20 against prior 0.08. Gilman's late-shedding row
  constrains it more than the plan assumed. `eta_lo` meets `tier2_plan.md` §5's
  stated trigger for a second look at the η confound, though priorsense flags 16
  of 26 parameters at similar magnitudes, so the signal is not specific. `[open]`
