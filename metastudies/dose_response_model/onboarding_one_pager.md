# Dose-Response Metastudy: Getting Your Bearings

**What are we doing?**

Calibrating a mechanistic dose-response model for typhoid fever. The model predicts probability of infection and disease as a function of bacterial dose and prior immunity. This is better than standard compartmental models because it generalizes across settings where exposure intensity varies by orders of magnitude.

**Why is this hard?**

We have two eras of human challenge data that don't directly compare:

| Era | Strengths | Gaps |
|-----|-----------|------|
| **Maryland (1960s-70s)** | Multi-dose data spanning 10^3 to 10^9 CFU; large N (~1700 total) | No immunity measurements; milk delivery |
| **Oxford (2010s)** | Measured anti-Vi IgG titers; modern protocols | Narrow dose range (10^3-10^4); bicarbonate delivery |

The delivery medium matters enormously—bicarbonate neutralizes stomach acid, so ~10^4 CFU in bicarbonate produces attack rates similar to ~10^7 in milk. This ~1000x "effective dose" difference must be handled somehow.

**The model structure**

Beta-Poisson dose-response with immunity scaling:

```
P(outcome | dose, CoP) = 1 - (1 + dose * (2^(1/α) - 1) / N50)^(-α / CoP^γ)
```

- **N50**: dose for 50% probability in naive individuals
- **α**: heterogeneity (host variability in susceptibility)
- **γ**: how strongly immunity scales effective dose
- **CoP**: correlate of protection (anti-Vi IgG titer; 1 = naive)

Two outcomes: **infection** (colonization/shedding) and **fever** (clinical disease). Same functional form, different parameters. The prototype uses fixed ratios: infection has 10x lower N50, 2x higher α, and 0.5x γ compared to fever.

**Key documents to read**

1. `dose_response_model_specification.md` — The reference model and simplification principles
2. `cross_cutting_observations.md` — Patterns across the literature (strain lineage, outcome definitions, era differences)
3. `yolo_working_model_notes.md` — First-pass thinking on what's fittable and what powers what
4. `cohort_incidence_model_proof_of_concept.R` (first 200 lines) — The working prototype

**Key papers**

- **Hornick 1966/1970**: Multi-dose data; anchors the dose-response shape
- **Waddington 2014**: Modern dose-escalation; establishes bicarbonate protocol
- **Darton 2016**: Anti-Vi IgG as correlate of protection; best immunity data
- **Gibani 2020**: Rechallenge study; raises questions about immunity vs innate susceptibility

**The big open questions**

1. **Medium offset**: How much of the Maryland-Oxford difference is gastric acid (biological) vs unknown baseline immunity in Maryland (confounding)?

2. **The Gibani paradox**: People who *didn't* get sick on first challenge are *more* protected on rechallenge. This suggests stable innate susceptibility differences, not just adaptive immunity. The current model can't capture this.

3. **Cross-era bridging**: We need Maryland for dose-response shape and Oxford for immunity effects. Connecting them requires assuming the heterogeneity parameter (α) is the same across eras.

**What success looks like**

A calibrated model that:
- Reproduces observed attack rates across doses and immunity levels
- Has interpretable parameters with quantified uncertainty
- Can be embedded in a transmission model for vaccine policy analysis

The parameters will be uncertain. That's okay. The goal is a model that's *less wrong* than ignoring dose-response entirely, with honest uncertainty quantification.

---

*Start with the model specification doc, then skim the extracts for papers that interest you, then read the YOLO notes to see one path through the calibration problem.*
