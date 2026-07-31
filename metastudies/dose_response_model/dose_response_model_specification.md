# Typhoid Dose-Response Model Specification

**Purpose**: Define the reference model—the most elaborate model the literature supports, what we would fit with optimal data and unlimited precision. This serves as a benchmark for reasoning clearly about which simplifications are forced by data limitations versus chosen for parsimony.

**Status**: Phase 3 complete. Sections 6-8 to be filled during Phases 4-5.

**Related documents**:
- `outcome_mapping.md` - Decision rules for mapping observed data to model variables
- `cross_cutting_observations.md` - Patterns across the extraction corpus

---

## 1. Latent Biological Processes

The reference model distinguishes these underlying processes:

### 1.1 Infection Cascade

```
Ingestion → Gastric survival → Colonization → {Systemic invasion, Stool shedding} → {Clinical disease, Chronic carriage}
```

- Bacterial survival through gastric barrier → colonization
- Colonization → systemic invasion (bacteremia)
- Colonization → stool shedding
- Systemic invasion → clinical disease
- Systemic invasion → chronic carriage (independent of clinical disease)

### 1.2 Outcomes to Model

Each outcome is a distinct random variable:

| Outcome | Latent Process | Observable Via |
|---------|---------------|----------------|
| **Infection** | Successful colonization | Latent, partially observed through downstream markers |
| **Bacteremia** | Bacteria in bloodstream | Blood culture, PCR |
| **Stool shedding** | Intestinal colonization with fecal excretion | Stool culture |
| **Fever** | Elevated temperature | Thermometry (multiple threshold definitions) |
| **Clinical typhoid** | Composite syndrome | Clinical diagnosis (definition varies by study/era) |
| **Chronic carriage** | Persistent gallbladder colonization | Prolonged persistent or intermittent stool culture positivity (>1 year) |
| **Seroconversion** | Antibody response | Serology (multiple assays, threshold-dependent) |

### 1.3 Severity Gradations

- Fever severity (temperature magnitude, duration)
- Disease severity (mild/moderate/severe; hospitalization; complications)
- Death

### 1.4 Joint/Conditional Relationships

#### Core Latent States

The causal chain from challenge to outcomes flows through these latent states:

**Gastric survival** (count distribution):
$$P(D_{\text{gastric}} | D_{\text{challenge}}, \text{strain}, \text{medium}, \text{host})$$
The dose surviving the gastric barrier depends on challenge dose, strain, delivery medium (milk, bicarbonate, etc.), and individual host factors.

**Colonization** (binary):
$$P(\text{colonization} | D_{\text{gastric}}, \text{strain}, \text{immunity}, \text{host})$$
Successful establishment in the gut.

**Systemic invasion** (binary):
$$P(\text{systemic invasion} | \text{colonization}, D_{\text{gastric}}, \text{strain}, \text{immunity}, \text{host})$$
Bacteria entering the bloodstream. *Note: Whether colonization is a necessary cause is uncertain; included in conditioning set as modeling choice.*

**Acute disease** (binary/graded):
$$P(\text{acute disease} | \text{systemic invasion}, \text{immunity}, \text{host})$$

**Chronic carriage** (binary):
$$P(\text{chronic carriage} | \text{systemic invasion}, \text{immunity}, \text{host})$$

**Stool shedding** (binary/duration):
$$P(\text{stool shedding} | \text{colonization}, D_{\text{gastric}}, \text{strain}, \text{immunity}, \text{host})$$

*Design note*: $D_{\text{gastric}}$ drops out of acute disease and chronic carriage conditioning. This encodes a memorylessness hypothesis: from systemic invasion forward, outcomes depend on the invasion state, not the original dose.

#### Observables

Biological states accessible from outside the body (without histopathology):

**Bacteremia**:
$$P(\text{bacteremia} | \text{systemic invasion}, \text{immunity}, \text{host})$$

**Fever** (ordinal, threshold-dependent):
$$P(\text{fever} | \text{acute disease}, \text{immunity}, \text{host})$$

**Death**:
$$P(\text{death} | \text{acute disease}, \text{immunity}, \text{host})$$

**Immunity** (pre-challenge):
$$P(\text{immunity} | \text{host})$$
The immune state at the time of challenge. This is the "immunity" appearing in conditioning throughout this section. See Section 2 for representation as a Correlate of Protection (CoP).

#### Observations

Measurements of observables, dependent on assay characteristics:

**Stool culture**:
$$P(\text{stool culture}^+ | \text{stool shedding}, \text{immunity}, \text{assay})$$
*Note: Immunity (bactericidal antibodies) affects culture viability.*

**Stool PCR**:
$$P(\text{stool PCR}^+ | \text{stool shedding}, \text{assay})$$

**Blood culture**:
$$P(\text{blood culture}^+ | \text{bacteremia}, \text{immunity}, \text{assay})$$
*Note: Immunity (bactericidal antibodies) affects culture viability.*

**Blood PCR**:
$$P(\text{blood PCR}^+ | \text{bacteremia}, \text{assay})$$

**Temperature recorded**:
$$P(T > \theta | \text{fever}, \text{measurement protocol})$$

**Serology titer** (observation of immunity):
$$P(\text{titer} | \text{immunity}, \text{assay})$$

**Seroconversion** (derived comparison of two titer observations):
$$P(\text{seroconversion} | \text{titer}_{\text{pre}}, \text{titer}_{\text{post}}, \text{assay})$$
*Note: Assay variance determines whether observed change is called "real."*

**Clinical typhoid diagnosis** (composite observation, study-dependent):
$$P(\text{clinical typhoid dx} | \text{systemic invasion}, \text{acute disease}, \text{immunity}, \text{host}) = \sum_{o \in \mathcal{O}_{\text{clinical}}} P(o | \cdot)$$
The definition varies by study/era—e.g., Maryland: ">103°F for 24-36hr + symptoms"; Oxford: "≥38°C for 12hr OR blood culture+". This composite aggregates over whichever observations constitute the clinical endpoint in a given protocol.

#### Correlation Structure

Outcomes are not fully independent given shared conditioning. Key correlations:
- Bacteremia duration and fever severity
- Shedding duration and seroconversion magnitude
- (To be elaborated based on data availability)

---

## 2. Immunity Representation

### 2.1 Immune Compartments

The true immune state comprises multiple compartments:

- **Humoral immunity**: Circulating antibodies (IgG, IgA) against various antigens (Vi, O, H, LPS)
- **Mucosal immunity**: Secretory IgA, gut-resident memory
- **Cellular immunity**: T cell responses (CD4+, CD8+), memory T cells

### 2.2 Correlate of Protection (CoP)

The **Correlate of Protection (CoP)** is a summary measure of the immune state sufficient to predict protection. In principle, it is an array of cellular and molecular responses; in practice, it collapses to a scalar.

This collapse to a scalar is implicit in the concept of a "correlate of protection"—the CoP framework assumes a low-dimensional summary sufficient to predict protection.

The pre-challenge CoP is the "immunity" appearing in conditioning throughout Section 1.

### 2.3 Scope Note: Post-Challenge Dynamics

The immune state evolves in response to infection. Pre-post dynamics are important to typhoid modeling but are treated as a separate data and calibration project in this repository. For this project, if seroconversion is used, we assume a lognormally distributed difference with standard error of approximately ±1 log2 unit.

*See Section 5 for mathematical formalization of CoP.*

---

## 3. Heterogeneity Sources

### 3.1 Study-Level

| Source | Impact |
|--------|--------|
| **Strain** | Quailes vs Ty2 vs wild-type may differ in virulence |
| **Delivery medium** | Milk vs bicarbonate affects gastric survival |
| **Fasting protocol** | Fasting duration before challenge affects gastric survival |
| **Subject population** | Naive volunteers vs endemic-area residents |
| **Era** | 1960s vs 2010s protocols, ascertainment methods |
| **Outcome definition** | Clinical typhoid criteria vary (e.g., ">103°F/36hr" vs "≥38°C/12hr OR culture+") |
| **Diagnostic methods** | Culture only vs culture + PCR affects infection ascertainment |

### 3.2 Individual-Level

- Age, sex
- Prior exposure history
- Baseline immune state
- Genetic susceptibility (e.g., HLA type)

---

## 4. DAG: Latent Processes → Observables → Observations

```
═══════════════════════════════════════════════════════════════════════════════
                              INPUTS & CONTEXT
═══════════════════════════════════════════════════════════════════════════════

    ┌───────────┐   ┌───────────┐   ┌───────────┐   ┌───────────┐
    │D_challenge│   │  Strain   │   │  Medium   │   │   Host    │
    └─────┬─────┘   └─────┬─────┘   └─────┬─────┘   └─────┬─────┘
          │               │               │               │
          └───────────────┴───────┬───────┴───────────────┘
                                  │
═══════════════════════════════════════════════════════════════════════════════
                            CORE LATENT STATES
═══════════════════════════════════════════════════════════════════════════════
                                  │
                                  ▼
    ┌───────────┐          ┌─────────────┐
    │    CoP    │─────────▶│  D_gastric  │
    │ (immunity)│          └──────┬──────┘
    └─────┬─────┘                 │
          │                       ▼
          │               ┌─────────────┐
          ├──────────────▶│ Colonization│
          │               └──────┬──────┘
          │                      │
          │         ┌────────────┴────────────┐
          │         │                         │
          │         ▼                         ▼
          │  ┌─────────────┐          ┌─────────────┐
          ├─▶│  Systemic   │          │    Stool    │
          │  │  Invasion   │          │  Shedding   │
          │  └──────┬──────┘          └──────┬──────┘
          │         │                        │
          │    ┌────┴────┐                   │
          │    │         │                   │
          │    ▼         ▼                   │
          │ ┌────── ┐ ┌─────────┐            │
          └▶│ Acute │ │ Chronic │            │
            │Disease│ │Carriage │            │
            └───┬───┘ └─────────┘            │
                │                            │
═══════════════════════════════════════════════════════════════════════════════
                              OBSERVABLES
═══════════════════════════════════════════════════════════════════════════════
                │                            │
       ┌────────┼────────┐                   │
       │        │        │                   │
       ▼        ▼        ▼                   ▼
  ┌────────┐┌───────┐┌───────┐         ┌───────────┐
  │Bactere-││ Fever ││ Death │         │   Stool   │
  │  mia   │└───┬───┘└───────┘         │  Shedding │
  └────┬───┘    │                      └─────┬─────┘
       │        │                            │
═══════════════════════════════════════════════════════════════════════════════
                              OBSERVATIONS
═══════════════════════════════════════════════════════════════════════════════
       │        │                            │
   ┌───┴───┐    │                       ┌────┴────┐
   │       │    │                       │         │
   ▼       ▼    ▼                       ▼         ▼
┌──────┐┌─────┐┌─────┐              ┌───────┐   ┌─────┐
│Blood ││Blood││Temp │              │Stool  │   │Stool│
│Culture││PCR ││ >θ  │              │Culture│   │ PCR │
└──────┘└─────┘└──┬──┘              └───────┘   └─────┘
                  │
                  ▼
          ┌────────────────┐
          │Clinical Typhoid│
          │   Diagnosis    │
          └────────────────┘

───────────────────────────────────────────────────────────────────────────────
                         IMMUNITY OBSERVATIONS
───────────────────────────────────────────────────────────────────────────────

    ┌───────────┐
    │    CoP    │
    │ (immunity)│
    └─────┬─────┘
          │
          ▼
    ┌───────────┐         ┌───────────────┐
    │  Titer    │────────▶│ Seroconversion│
    │  (assay)  │         │  (derived)    │
    └───────────┘         └───────────────┘
```

**Legend:**
- Boxes = random variables
- Arrows = conditional dependencies
- CoP influences all latent states (arrows simplified for readability)
- Strain conditions D_gastric, colonization, systemic invasion, shedding
- Host conditions all latent states

---

## 5. Mathematical Framework

This section formalizes the relationships described qualitatively above.

### 5.1 General Form

For each outcome $Y$:

$$
P(Y | D, \text{CoP}, \vec{\theta}_Y) = f_Y(D, \text{CoP}; \vec{\theta}_Y)
$$

Where:
- $D$ = bacterial dose (CFU)
- $\text{CoP}$ = scalar correlate of protection (see Section 2)
- $\vec{\theta}_Y$ = outcome-specific parameters

### 5.2 Outcome Hierarchy

Outcomes decompose along the causal chain defined in Section 1.4:

$$P(\text{acute disease}) = P(\text{acute disease} | \text{systemic invasion}) \cdot P(\text{systemic invasion} | \text{colonization}, D_{\text{gastric}}) \cdot P(\text{colonization} | D_{\text{gastric}}) \cdot P(D_{\text{gastric}} | D_{\text{challenge}})$$

(Conditioning on strain, CoP, host suppressed for readability.)

Similar decompositions for:
- Bacteremia given systemic invasion
- Fever given acute disease
- Shedding given colonization and $D_{\text{gastric}}$
- Chronic carriage given systemic invasion
- Clinical typhoid given acute disease and observations

### 5.3 Mechanistic Basis: Modified Beta-Poisson

The working model uses a modified beta-Poisson framework:

$$
P(\text{infection} | D, \text{CoP}) = 1 - \left(1 + D \cdot \frac{2^{1/\alpha} - 1}{N_{50}}\right)^{-\alpha / \text{CoP}^{\gamma}}
$$

This form has mechanistic interpretation:
- Each bacterium has independent probability of initiating infection
- $N_{50}$ = dose for 50% probability of outcome in naive individuals
- $\alpha$ = unexplained heterogeneity of outcome given dose (aggregates all sources of variation in the cascade from challenge dose to outcome)
- $\gamma$ = how strongly immunity scales effective dose

### 5.4 CoP Formalization

The CoP is in principle an array:

$$\vec{\text{CoP}} = (\text{CoP}_1, \text{CoP}_2, \ldots, \text{CoP}_n)$$

For the working model, this collapses to a scalar:

$$\text{CoP} = g(\vec{\text{CoP}})$$

**Key simplification assumptions** (to be documented in Section 6.5):
- Which immune compartments dominate protection?
- How do different compartments combine (additive? multiplicative? minimum?)?
- What is lost by ignoring dynamics within a challenge study timescale?

### 5.5 Titers as Observations of CoP

Each serological assay measures one (or a combination of) CoP components:

$$P(\text{titer} | \vec{\text{CoP}}, \text{assay})$$

For example, anti-Vi IgG ELISA (Vacczyme assay) observes a specific humoral component.

**Seroconversion** is a derived quantity comparing two titer observations:
$$P(\text{seroconversion} | \text{titer}_{\text{pre}}, \text{titer}_{\text{post}}, \text{assay})$$

Assay variance determines whether an observed change is classified as "real" seroconversion.

### 5.6 Post-Challenge CoP Dynamics

The immune state evolves in response to infection:
$$\vec{\text{CoP}}_{\text{post}} = f(\vec{\text{CoP}}_{\text{pre}}, \text{colonization}, \text{systemic invasion}, \text{acute disease}, \ldots)$$

*Scope note*: This is treated as a separate project. See Section 2.3.

---

## 6. Simplifications Needed for Practical Model

### 6.1 Principles of Simplification

Building a practical working model requires integrating out degrees of freedom in the reference model. To see this in action, consider the first simplification one encounters when building a dose-response model.

In our framework (which is already simplifying interesting biology), we assume the probability of colonization depends on the strain, pre-challenge immunity, host factors, and the gastric dose (CFU that reach the gut). But the gastric survival process and the array of host factors that influence an individual's colonization outcome are unobservable. A model tied only to observables must integrate them out. Formally:

$$
P(\text{colonization} | D_{\text{chall}}, \text{strain}, \text{medium}) = \int \mathcal{D}D_{\text{gastric}} \, \mathcal{D}\text{host} \; P(\text{colonization} | D_{\text{gastric}}, \text{host}, \text{strain}) \, P(D_{\text{gastric}} | D_{\text{chall}}, \text{strain}, \text{medium}, \text{host}) \, P(\text{host})
$$

This is where dose-response modeling usually starts, and why all realistic dose-response models need to accommodate generic heterogeneity.

### 6.2 Why Beta-Poisson

The beta-Poisson model is natural in this context. It derives from the assumption that colonization is a "single-hit" phenomenon: under zero heterogeneity, the probability of colonization equals the probability that at least one "colony-forming unit" (CFU) successfully establishes a colony. The beta distribution then accounts for how each CFU has a different probability of success due to all the vagaries of gastric passage, delivery medium, and host factors.

### 6.3 Non-Dose Heterogeneity

For variation not involving dose, the natural families of simplifying models depend on the nature of the outcomes and are more familiar: normally-distributed host effects, multinomial/1-of-N factor choices, ordinal fever scales (without concern for fever trajectories), and all the standard exponential family regressions.

### 6.4 Cascaded Outcomes

Done correctly, conditional outcomes require nonlinear regression. For example, the relationship between fever and challenge dose is mediated by acute disease and systemic invasion. Formally:

$$
P(\text{fever} | D_{\text{chall}}, \text{immunity}, \text{strain}, \text{host}, \text{medium}) = P(\text{fever} | \text{acute disease}, \text{immunity}, \text{host}) \cdot P(\text{acute disease} | \text{systemic invasion}, \text{immunity}, \text{host}) \cdot P(\text{systemic invasion} | \text{colonization}, D_{\text{chall}}, \text{strain}, \text{immunity}, \text{medium}) \cdot P(\text{colonization} | D_{\text{chall}}, \text{strain}, \text{medium})
$$

Systemic invasion and colonization independently depend on $D_{\text{chall}}$ (by assumption). Biologically, this reflects a hypothesis that systemic invasion is possible without mounting a stable enteric colony and thus has its own dose-response dynamics. The probability of fever is therefore a product of multiple conditional processes, some of which are dose-dependent.

But done practically, the product of monotonic sigmoids is monotonic and definitely not a standard sigmoid—yet it looks like one with any finite data. It may therefore be reasonable to assume the product collapses to a single dose-response curve, although this is only fully justified if all but one of the stages in the cascade are not dose-dependent.

It thus only makes sense practically to keep nonlinear cascades when they are necessary to capture dynamics required by your reasons for modeling.

### 6.5 Specific Simplifications

Status 2026-07-31: the practical model is `calibration/typhoid_dose_response.stan`;
Tier 2 design is `joint_inference_plan.md`. This table records what was
actually adopted, so the reference model and the fitted model stay reconcilable.

| Reference Model Component | Simplification adopted | Justification / where |
|---|---|---|
| Multiple immune compartments | Single scalar CoP (anti-Vi IgG, VaccZyme EU/mL, naive = 1) | plan §2.4; axis caveat below |
| Cascaded dose-response | Retained as a two-stage cascade `P_inf x P_fev|inf` | Per §6.4, collapsing is justified only if all but one stage is dose-independent; both stages here are dose-dependent. The product is not itself a beta-Poisson. |
| Gastric survival / delivery medium | One scalar `delta`, milk -> bicarbonate-equivalent dose | plan §2.3. Strongly confounded with both N50s (r = -0.72 / -0.78). |
| Joint outcome distribution | Independent outcomes **except where one cohort supplies two endpoints**, which are factorized marginal + conditional | See the no-double-counting rule below. |
| Fever definition heterogeneity | Fitted definition map `phi(T,D)` | plan §2.5; §3.1 "Outcome definition" |
| Infection definition heterogeneity | Fitted definition map `psi_d` (Tier 2) | plan §2.8; §3.1 "Diagnostic methods" |
| Shedding ascertainment (Oxford) | `eta(D)` treatment-truncation correction (Tier 2) | plan §2.7 |
| Unobserved Maryland immunity | Two-component latent mixture (`pi_susc`, `CoP_susc`, `CoP_imm`) | plan §2.4 |
| Extra-binomial variation | One shared parameter, beta-binomial `grand_overdispersion_rho` | plan §5.5 |
| Chronic carriage, seroconversion, death, incubation period | Not modelled | No likelihood term exists; §1.2/§1.3 are aspirational for these |

**No-double-counting rule (package-wide).** Where a single cohort of volunteers
supplies two endpoints, they are NOT two independent marginals — they are factorized
into a marginal and a conditional on the same men:

- Hornick 10^7: `H-I-7` (28/30 infected) + `H-FgI-7` (16/28 fever|infected). The
  marginal fever row `H-F-7` was **deleted** from the dataset for double-counting.
- Darton placebo: groups 6/7, infection over all 30 + fever|infection over the infected.
- Levine is the outstanding exception — `Lev-F-k` and `Lev-I-k` are two marginals over
  the same men, treated as independent. **Known limitation, tolerated:** the clean fix
  needs fever nested inside the infection endpoint, but Levine's endpoint is stool-only
  and fever is NOT nested in it (Darton: 7 of 20 TD+ stool-negative), and Levine
  publishes no cross-tabulation. Mitigated only in LOO, where the pair is one unit.

**CoP caveat — the axis is known to be wrong for the Maryland era.** The scalar CoP is
anti-Vi IgG. In the Maryland volunteers, every paper that looked found baseline anti-Vi
did **not** predict outcome (Gilman 1977 p.721; Levine 1976 p.427), and purified Vi
vaccine gave a 13-fold antibody rise with zero protection (Woodward 1980 p.554), while
**H agglutinin did** stratify. The Maryland latent mixture is therefore an
anti-Vi-*equivalent* device standing in for a non-Vi immunity axis, not a measurement
of one. Documented at `calibration/priors.yaml` (`pi_susc`).

**Fever severity.** `phi(T,D)` is the implemented severity model (§1.3). Known
limitation: the dose-lift term carries no T dependence, so `phi -> 1` uniformly in T as
dose rises — at Hornick's top dose the model implies `phi(41C) = 0.95`. It represents
"threshold choice stops mattering at high dose", not "higher dose gives hotter fevers".
Consequence of pinning `beta_phi = 1`.

---

## 7. Mapping to Available Data

Verified against `calibration/dose_response_data.csv` and
`analysis_data/darton_individual_endpoints.csv`, 2026-07-31.

**Scale.** 37 CSV rows: 25 `tier1_active`, 31 `tier2_active`, 6 `validation_only`.
After Darton individualization, Tier 1 presents **80 observations** to Stan — 24
grouped rows plus 56 per-subject n=1 rows (30 infection + 26 fever|infection).

| Reference component (§1.2–1.3) | Data available | Quality | Notes |
|---|---|---|---|
| **Dose** | All 37 rows | Variable | 10^3–10^9 CFU overall, but the eras barely overlap: Maryland 10^3–10^9 (milk), Oxford 10^3–2x10^4 (bicarbonate). The `delta` bridge is what joins them, and it is confounded with both N50s. Published dose *ranges* are collapsed to point estimates. |
| **Pre-challenge immunity** | 13 of 37 rows | Oxford only | `CoP` is NA for **all 24 Maryland rows** — no 1960s serology exists. Maryland immunity is latent (the mixture). Darton adds 30 per-subject anti-Vi titres, but only **12 of 30 are above detection**, which is why `gamma_inf` and `gamma_fevginf` do not separate. |
| **Infection** | 12 rows (6 Maryland, 6 Oxford) | Definition-heterogeneous | The endpoint differs materially by study — Hornick stool-or-blood culture, Levine any-time stool, Gilman late shedding 4–30 d, Oxford stool shedding. This is what `psi_d` (§2.8 of the plan) exists to reconcile. All 6 Oxford infection rows are Tier 2 only. |
| **Bacteremia** | Darton per-subject only | Good | 20 of 30 placebo bacteremic; 20 blood-positive. Not a separate likelihood layer — folded into `bact_or_stool` (26 of 30) as the broad infection marker. No Maryland bacteremia counts. |
| **Stool shedding** | Oxford (Tier 2) + Levine + Gilman + Darton | Good where present | Darton measures it against the broad marker in the same men: 19 stool+ of 26 `bact_or_stool` — the anchor for `psi_stool`. Note fever is **not** nested in stool-positivity: 7 of 20 TD+ were stool-negative. |
| **Fever** | 24 rows (15 Maryland, 9 Oxford) | Good | Threshold definitions vary and are modelled by `phi(T,D)`: Hornick 39.4C, Levine 38.3C, Gilman assumed 38.3C (imported, see `data_prep.R`), Oxford composite TD at T_ref = 38.0C. |
| **Fever severity** | Darton ladder only | Good but thin | Among 20 TD+ placebo subjects: 20/19/16/10/8 crossing >=37/37.5/38/38.5/39 C. Pins `phi0(T)`; **39.0 C is the maximum**, so Hornick's 39.4 C is extrapolation. |
| **Fever given infection** | 1 Maryland row + 26 Darton subjects | Sparse | `H-FgI-7` (16/28) is the only grouped conditional. This layer is the thinnest in the whole dataset. |
| **Clinical typhoid (composite)** | Oxford TD | Good | Taken as the reference definition, so `phi == 1` for Oxford by construction. |
| **Chronic carriage** | None | — | No likelihood term. |
| **Seroconversion** | Post-challenge only | Not usable as an outcome | Darton has baseline/pre-challenge titres; post-challenge seroconversion is not modelled. |
| **Death / severe disease** | None | — | No deaths in the Maryland program (Woodward 1980). No likelihood term. |
| **Incubation period** | Available, unused | Good | Hornick Table 1 and Waddington Table 3 report time-to-event. Not modelled; could constrain `alpha` and the dose-response shape. |

**What the mapping implies.** Three reference components carry most of the model's
weight (dose, fever, infection) and three are effectively single-source: fever severity
and fever-given-infection rest almost entirely on one Darton cohort of 30, and all
Maryland immunity is latent. Carriage, seroconversion, death and incubation have no
likelihood term at all — §1.2 and §1.3 remain aspirational for those.

---

## 8. Calibration Setup

Implemented in `calibration/typhoid_dose_response.stan`; Tier 2 design and locked
decisions in `joint_inference_plan.md`. This section is kept in sync with
that plan — if they disagree, the plan is authoritative and this is stale.

### 8.1 Likelihood Structure

| Component | Form | Notes |
|---|---|---|
| Dose-outcome counts | Beta-binomial | `y ~ beta_binomial(n, p*k, (1-p)*k)`, `k = (1-rho)/rho`, `rho = grand_overdispersion_rho`. `rho = 0` is exactly the binomial. Grouped rows only; n=1 rows stay binomial. |
| Fever definition | Fitted map `phi(T,D)` | Not stratification — a measurement model, pinned by the Darton temperature ladder |
| Infection definition | Fitted map `psi_d` (Tier 2) | Pinned by the Darton stool-vs-`bact_or_stool` cross-tab, 19/26 |
| Shedding ascertainment | `eta(D)` (Tier 2, Oxford) | Treatment truncation; distinct mechanism from `psi` |
| Dose uncertainty | **Point estimates** | Oxford reports ranges; NOT integrated over. Open. |
| Missing immunity | Latent 2-component mixture | Maryland only |

### 8.2 Study-Level Effects

| Effect | Structure | Rationale |
|---|---|---|
| Strain | Excluded | Quailes throughout the calibration set |
| Delivery medium | Fixed scalar `delta` | Milk -> bicarbonate-equivalent dose |
| Era | Not modelled separately | Constraint: nothing may be given per-cohort freedom on `N50`. Hornick's cohorts ARE the dose ladder, so such a term flattens the dose-response while improving fit. |
| Cohort / study | Single shared overdispersion, `grand_overdispersion_rho` | Tier 1 has 16 cohorts over 24 grouped observations, 10 of them singletons, so per-cohort terms are ~one parameter per datum. Reasoning: `calibration/cohort_random_effects_design.md`. |
| Outcome definition | Measurement model, not stratification | `phi(T,D)` for fever, `psi_d` for infection |
| Immunology across eras | **Assumed invariant** | `gamma_inf` / `gamma_fevginf` are SHARED between Maryland and Oxford by design: human immunology is taken to be the same in both. This is why an era-specific protection scale is not an available fix for Maryland misfit. |

### 8.3 Individual-Level Nuisances

| Nuisance | Treatment | Notes |
|---|---|---|
| Host heterogeneity | Absorbed into `alpha` | Per §6.2 |
| Age / sex | Not used | Unavailable for most Maryland rows. Glynn 1995 finds age <30 RR 1.79 in this population, so this is a real unmodelled confound. |
| Prior exposure | Latent mixture (Maryland); measured CoP (Oxford) | Maryland `pi_susc ~ 0.62`; concordant with Gilman's H-negative fraction 36/53 = 0.68 and Woodward's non-veteran fraction 200/305 = 0.66 |
| Cohort membership | Recorded as `cohort_id` | Not passed to Stan; no likelihood term reads it. Post-fit it keys the LOO units (80 rows -> 46), so model comparison does not count the same volunteers twice. |

### 8.4 Prior Specification

Single source of truth: `calibration/priors.yaml` (hyperparameters are Stan *data*).

| Parameter | Prior | Rationale |
|---|---|---|
| `log10_N50_inf` | Normal(2.5, 1.0) | Order-of-magnitude uncertainty |
| `alpha_inf`, `alpha_fevginf` | LogNormal(-1.5, 0.8) | Positive, weakly informative |
| `gamma_inf`, `gamma_fevginf` | LogNormal(-1.6, 0.9) | Anchored to Darton HR 0.29 per log10 anti-Vi |
| `log10_delta` | Normal(3.5, 0.7) | Milk-to-bicarb bridge |
| `pi_susc` | Beta(7, 4) | ~0.65; provenance caveat in priors.yaml |
| `CoP_imm` | Exponential(0.074) | **Prior-carried, unidentified** — no 1960s serology exists |
| `grand_overdispersion_rho` | Beta(1, 49) | ICC; `rho = 0` is the binomial. Anchored to sigma ~ 0.26 measured across the five Maryland 10^5 cohorts |

### 8.5 Identifiability Concerns

Measured on `results/tier1`, not hypothetical:

| Issue | Status | Diagnostic / resolution |
|---|---|---|
| `N50` vs `delta` | **Confounded**: r = -0.72 (inf), -0.78 (fev\|inf) | Structural; `delta_bridge.png` shows it |
| `CoP_imm` | **Unidentified**, prior-carried (priorsense prior 0.386 vs lik 0.071) | No Maryland serology; report as an assumption |
| `gamma_inf` vs `gamma_fevginf` | **Do not separate** (both ~0.16-0.18) | Darton's titre range is thin: 12/30 above detection |
| `psi_stool` vs `eta` | Confounded at Oxford (both multiply `P_inf`) | Darton cross-tab pins `psi*eta(18200)`; only eta's dose-dependence separates them |
| `psi_late` (Gilman 4-30 d) | Will be prior-carried | No cross-tab exists; constrain below `psi_stool` |
| Cohort effects vs dose-response | Held by using one shared overdispersion parameter | Per-cohort offsets would compete with `N50_inf`/`alpha_inf` |
| `sigma_study` | Declared in the .stan, given a prior, used in zero likelihood terms | Dead; delete when Step 2 lands |

### 8.6 Gates

Before adopting any Tier 2 increment:
1. `log10_N50_inf` and `alpha_inf` must not move materially (the signal-absorption gate).
2. priorsense on the new parameter; if prior-dominated, report it as such.
3. LOO must improve, using `compute_loo_units()` grouping (which merges the Hornick
   marginal + conditional into one unit).
4. The `Gil-F-Hlo` residual is **not** a gate — it is a within-cohort stratum contrast
   and is expected to persist. If it vanishes, something else moved; explain it.
