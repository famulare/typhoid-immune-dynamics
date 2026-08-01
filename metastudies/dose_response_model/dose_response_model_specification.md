# Typhoid Dose-Response Model Specification

**Purpose**: Define the reference model—the most elaborate model the literature supports, what we would fit with optimal data and unlimited precision. This serves as a benchmark for reasoning clearly about which simplifications are forced by data limitations versus chosen for parsimony.

**Status (2026-08-01)**: Reference model and practical mapping are complete.
**§5.7 gives the implemented model.** It was checked against
`calibration/typhoid_dose_response.stan`, `data_prep.R`, `tier_specs.R`,
`priors.yaml`, and the 2026-08-01 fit. Sections 6-8 record the simplifications,
data mapping, Tier 2 decisions, and remaining identifiability limits. The
canonical fit is `t2-indiv-vax`, stage `phi-rho-eta-psi-vax`, N_obs = 182.

Sections 1-4 and 6.1-6.4 describe the *reference* model; sections 5.7, 7, and 8
describe the *fitted* model. If they differ, §5.7 describes the fitted model.
Program dataflow and Stan-block dispatch are in `stan_model_structure.md`. The
generated summary is authoritative for numbers; where it and this document
disagree, it wins.

**Related documents**:
- `stan_model_structure.md` - Stan program dataflow and likelihood dispatch
- `notes/outcome_mapping.md` - Decision rules for mapping observed data to model variables
- `notes/cross_cutting_observations.md` - Patterns across the extraction corpus
- `joint_inference_plan.md` - Locked joint likelihood and identifiability assumptions
- `calibration/cohort_random_effects_design.md` - Locked design for one shared overdispersion parameter
- `calibration/jin_within_arm_titre_plan.md` - Draft plan for within-arm titre quadrature
- `calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md` - Canonical Tier 2 fit; authoritative for numbers

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
- These correlations are not separately parameterized in the practical fit; the
  retained infection -> fever cascade and the observation-specific measurement
  maps carry the dependency that the available data can support.

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

The implemented kernel adds a protection factor. See §5.7; that section
supersedes this equation for the fitted model.

### 5.4 CoP Formalization

The CoP is in principle an array:

$$\vec{\text{CoP}} = (\text{CoP}_1, \text{CoP}_2, \ldots, \text{CoP}_n)$$

For the working model, this collapses to a scalar:

$$\text{CoP} = g(\vec{\text{CoP}})$$

**Key simplification assumptions** (documented in Section 6.5 and the joint plan):
- A scalar CoP is sufficient for the fitted Oxford protection channel.
- Maryland protection is represented by a latent two-component mixture because
  baseline anti-Vi is unavailable.
- Within-challenge immune dynamics are outside the current likelihood.

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

### 5.7 The implemented model

Sections 5.1–5.6 are the reference model. Where they differ from this section,
this section describes the fitted model. It was checked against
`calibration/typhoid_dose_response.stan`, `calibration/data_prep.R`,
`calibration/tier_specs.R`, and `calibration/priors.yaml` on 2026-08-01. For the
mapping from these expressions to Stan blocks, see `stan_model_structure.md`.

#### 5.7.1 The kernel

The model uses one kernel for all dose-response probabilities. For non-anti-Vi
vaccine protection, the factor $V$ divides the exponent:

$$
P(D_{\text{eff}}; N_{50}, \alpha, \text{CoP}, \gamma, V) = 1 - \left(1 + D_{\text{eff}} \cdot \frac{2^{1/\alpha} - 1}{N_{50}}\right)^{-\alpha / (\text{CoP}^{\gamma} \cdot V)}
$$

$V = 1$ recovers §5.3. $V$ is not a factor on dose or probability, and the model
has no interaction term between vaccine protection and the subject's titre.

$D_{\text{eff}} = D$ for Oxford and $D / \delta$ for Maryland (§5.7.4).

#### 5.7.2 The two-stage cascade and the N50 ordering constraint

$$P_{\text{inf}} = P(D_{\text{eff}}; N_{50}^{\text{inf}}, \alpha_{\text{inf}}, \text{CoP}, \gamma_{\text{inf}}, V) \qquad P_{\text{fev}|\text{inf}} = P(D_{\text{eff}}; N_{50}^{\text{fev}|\text{inf}}, \alpha_{\text{fev}|\text{inf}}, \text{CoP}, \gamma_{\text{fev}|\text{inf}}, V)$$

The two $N_{50}$ values are **not** independent. The model samples
$\log_{10} N_{50}^{\text{inf}}$ and a non-negative offset $d_{\text{fev}}$:

$$\log_{10} N_{50}^{\text{fev}|\text{inf}} = \log_{10} N_{50}^{\text{inf}} + d_{\text{fev}}, \qquad d_{\text{fev}} \ge 0$$

This constrains $N_{50}^{\text{fev}|\text{inf}} \ge N_{50}^{\text{inf}}$: fever
requires at least as much dose as infection. It is a biological assumption. All
three original prior terms are retained, and the unit-Jacobian change of variables
leaves the prior unchanged.

#### 5.7.3 Observation-process maps

**Fever definition, $\varphi(T, D_{\text{eff}})$:** Maryland rows only.

$$\varphi_0(T) = \text{logit}^{-1}\!\left(\varphi_{0a} - \varphi_{0b}(T - T_{\text{ref}})\right), \qquad T_{\text{ref}} = 38.0\ ^\circ\text{C}$$

$$\varphi(T, D_{\text{eff}}) = \varphi_0(T) + (1 - \varphi_0(T)) \cdot P^{\text{naive}}_{\text{fev}}(D_{\text{eff}})$$

$P^{\text{naive}}_{\text{fev}}$ is the hard-coded cascade product at $\text{CoP} = 1$.
Because immunity already acts through $\text{CoP}^\gamma$, this avoids counting it
twice. The dose-lift exponent $\beta_\varphi$ is **pinned to 1** and does not
appear in the code. See §6.5.

**Infection definition, $\psi$:** applies only to Maryland infection rows
(5 of 182 rows in the canonical fit).

$$\psi(\text{def}) = \begin{cases} 1 & \text{Hornick, stool-or-blood (broad reference)} \\ \psi_{\text{stool}} & \text{Levine, any-time stool} \\ \psi_{\text{stool}} \cdot f_{\text{late}} & \text{Gilman, late shedding (days 4–30)} \end{cases}$$

The multiplicative form gives $\psi_{\text{late}} \le \psi_{\text{stool}}$: a
narrower window cannot detect more. A `psi_active` flag sets $\psi \equiv 1$
when the infection-definition map is inactive.

**Shedding detection, $\eta(D_{\text{eff}})$** — grouped Oxford shedding rows only:

$$\eta(D_{\text{eff}}) = \eta_{\text{lo}} + (1 - \eta_{\text{lo}}) \exp\!\left(-\kappa \, D_{\text{eff}} / N_{50}^{\text{inf}}\right)$$

$\eta$ decreases from 1 toward $\eta_{\text{lo}}$ on a dose scale tied to
$N_{50}^{\text{inf}}$. $\kappa$ is the most prior-dominated parameter in the
retained fit (§8.5).

**Maryland immunity mixture:**

$$M(D_{\text{eff}}, \cdot) = \pi_{\text{susc}} P(D_{\text{eff}}; \cdot, \text{CoP}_{\text{susc}}) + (1 - \pi_{\text{susc}}) P(D_{\text{eff}}; \cdot, \text{CoP}_{\text{imm}})$$

For Maryland **fever**, the mixture is over the cascade **product**:
$\mathbb{E}[P_{\text{inf}} P_{\text{fev}|\text{inf}}]$, not
$\mathbb{E}[P_{\text{inf}}]\,\mathbb{E}[P_{\text{fev}|\text{inf}}]$. These are not
algebraically equal. The Hornick conditional row is a **ratio of two
differently weighted mixtures**.

#### 5.7.4 The seven observation groups

Every likelihood row is one of seven kinds. `obs_prob()` has no other branches.

| # | group | success probability | δ | φ | ψ | η | V |
|---|---|---|---|---|---|---|---|
| 1 | `ox_fev` | $P_{\text{inf}}(D,\text{CoP}) \cdot P_{\text{fev}\mid\text{inf}}(D,\text{CoP})$ | — | — | — | — | — |
| 2 | `ox_inf` | $\eta(D) \cdot P_{\text{inf}}(D,\text{CoP})$ | — | — | — | ✓ | — |
| 3 | `md_fev` | $\varphi(T, D/\delta) \cdot M(\text{cascade product})$, or a Gilman-stratum branch at $\text{CoP}_{\text{susc}}$ / $\text{CoP}_{\text{imm}}$ | ✓ | ✓ | — | — | — |
| 4 | `md_inf` | $\psi(\text{def}) \cdot M(P_{\text{inf}}, D/\delta)$ | ✓ | — | ✓ | — | — |
| 5 | `hornick_cond` | $\left[\varphi \cdot M(\text{cascade product})\right] / M(P_{\text{inf}})$, clamped to $[10^{-12}, 1-10^{-12}]$ | ✓ | ✓ | — | — | — |
| 6 | `ox_inf_indiv` | $P_{\text{inf}}(D, \text{CoP}_{\text{subject}}, V)$ | — | — | — | — | ✓ |
| 7 | `ox_fevginf_indiv` | $P_{\text{fev}\mid\text{inf}}(D, \text{CoP}_{\text{subject}}, V)$ | — | — | — | — | ✓ |

The dispatch implies:

- **Groups 6 and 7 get no $\eta$.** Darton's per-subject endpoint is
  `bact_or_stool`, the broad blood-or-stool marker, so there is no stool
  truncation to correct. $\eta$ applies only to the five grouped Oxford shedding
  rows and to the $\psi$ cross-tab term.
- **A vaccine's $V$ enters the composite fever probability twice**, once through
  group 6 and once through group 7, because both route through the kernel's $V$.
- **The Gilman H-strata bypass the mixture.** A nonzero `gilman_stratum` selects a
  single component ($\text{CoP}_{\text{susc}}$ or $\text{CoP}_{\text{imm}}$)
  instead of mixing, still multiplied by $\varphi$. Two rows use this branch.

#### 5.7.5 CoP construction

$\text{CoP} = \text{anti-Vi IgG} / 3.7$, where 3.7 EU/mL is the imputed value at
the VaccZyme assay's limit of detection, so **CoP = 1 is the LLD, not an interior
point**. The distribution is left-censored with an atom at 1: 64 of the 90 Darton
per-subject rows sit at exactly 1.0.

Two constructions are used; this distinction affects identification of the titre
slope (see §8.5):

- **Darton** (groups 6/7) enters **per subject**, each with their own measured
  pre-challenge titre.
- **Jin, Waddington, Gibani** (groups 1/2) enter as **arm-level rows at the arm's
  geometric-mean titre** — Jin control 2.16, Vi-PS 38.11, Vi-TT 152.16; Waddington
  and Gibani naive at 1. Plugging a GMT into a nonlinear response is a Jensen
  approximation. Its measured size at the retained posterior is small
  (`joint_inference_plan.md` §5.3). A draft quadrature plan is in
  `calibration/jin_within_arm_titre_plan.md`.
- Maryland rows have no serology at all; immunity is the latent mixture.

#### 5.7.6 Observation distribution and the decoupled sub-likelihoods

Every one of the $N_{\text{obs}}$ rows uses the same beta-binomial:

$$y_i \sim \text{BetaBinomial}\!\left(n_i,\; p_i k,\; (1-p_i)k\right), \qquad k = \frac{1 - \rho}{\rho}$$

$\rho = 0$ is exactly the binomial. There is **no branch** for $n = 1$: the code
relies on $\text{BetaBinomial}(1, pk, (1-p)k) \equiv \text{Bernoulli}(p)$ for any
$k$. The 153 per-subject rows do not inform $\rho$ or receive its discount. The 29
grouped rows estimate $\rho$ and are discounted by it (§8.5). $p_i$ is clamped to
$[10^{-12}, 1-10^{-12}]$ before use, because the beta-binomial requires strictly
positive shape parameters where the binomial tolerates $p \in \{0,1\}$.

Two further likelihood terms are not dose-response rows and are not overdispersed:

1. **Darton temperature ladder** → identifies $\varphi_0(T)$. Three binomials:
   16, 10, and 8 of 20 TD+ placebo subjects crossing 38.0, 38.5, 39.0 °C. The
   extract also records 20/20 at 37.0 and 19/20 at 37.5; those two rungs are
   **deliberately excluded** [Mike, 2026-08-01] — they sit below $T_{\text{ref}}$
   where $\varphi_0$ should be ≈1 by construction, they are near-saturated, and
   including them would inflate $\varphi_{0b}$ without adding slope information.
2. **Darton stool-vs-broad cross-tab** → pins $\psi_{\text{stool}} \cdot \eta(18200)$
   jointly. One binomial, **15 of 26** — the nested intersection, not the marginal
   19 of 26. Gated by `psi_active`.

`log_lik` has length $N_{\text{obs}} + N_{\text{ladder}} + \texttt{psi\_active}$
= 182 + 3 + 1 = 186 for the canonical fit, so
$\text{target} = \text{lprior} + \sum \text{log\_lik}$ holds exactly and priorsense
power-scales the whole likelihood. LOO deliberately uses only the first
$N_{\text{obs}}$ entries.

#### 5.7.7 Generated quantities

Posterior-predictive `p_pred` / `y_rep` per row, plus reference-dose curves
($P_{\text{inf}}$, $P_{\text{fev}}$ at $10^3$/$10^4$ naive; Maryland fever at
$10^3$/$10^5$/$10^7$; the Hornick conditional; $\varphi_0$ at 38.3 and 39.4;
$\varphi$ at Hornick 39.4 across doses; $\eta$ at $10^3$/$10^4$).

Generated quantities are not included in `summary.csv` or `summary.md`, which cover
parameters and transformed parameters only (26 variables). GQ are stored in
`fit.rds` and are read by the figures.

The PPC rows provide the fitted arm probabilities. The generated quantities do not
affect sampling.

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

Status 2026-08-01: the practical model is `calibration/typhoid_dose_response.stan`;
Tier 2 design is `joint_inference_plan.md`. This table records what was
actually adopted, so the reference model and the fitted model stay reconcilable.

| Reference Model Component | Simplification adopted | Justification / where |
|---|---|---|
| Multiple immune compartments | Single scalar CoP (anti-Vi IgG, VaccZyme EU/mL, naive = 1) | plan §2.3; §5.7.5; axis caveat below |
| Non-anti-Vi vaccine protection | Per-vaccine scalar `V` dividing the exponent (`V_M01ZH09`, `V_Ty21a`; `V = 1` at placebo) | Ty21a is Vi-negative and M01ZH09 raised no anti-Vi, so neither protection can run through the CoP channel. §5.7.1. **Independence from anti-Vi is structural, not tested** — `V` multiplies `CoP^gamma`, so no interaction term exists to detect. |
| Cascaded dose-response | Retained as a two-stage cascade `P_inf x P_fev\|inf` | Per §6.4, collapsing is justified only if all but one stage is dose-independent; both stages here are dose-dependent. The product is not itself a beta-Poisson. |
| Relative position of the two N50s | **Constrained**, not estimated freely: `log10_N50_fevginf = log10_N50_inf + d_fev`, `d_fev >= 0` | §5.7.2. A biological ordering assumption the data cannot overturn; also removed the density cliff causing ~99% of divergences. |
| Gastric survival / delivery medium | One shared scalar `delta` (fitted, not fixed), milk -> bicarbonate-equivalent dose | plan §2.1–2.2. Trades off with both N50s: r = **-0.54** (inf), **-0.59** (fev\|inf) in the retained fit — moderate, not a hard confound. See §8.5. |
| Joint outcome distribution | Independent outcomes **except where one cohort supplies two endpoints**, which are factorized marginal + conditional | Partly violated; see the no-double-counting rule below. |
| Fever definition heterogeneity | Fitted definition map `phi(T,D)`, dose-lift exponent pinned to 1 | plan §2.5; §5.7.3; §3.1 "Outcome definition" |
| Infection definition heterogeneity | Fitted definition map `psi` (`psi_stool`, `frac_late`; Tier 2) | plan §2.8; §5.7.3; §3.1 "Diagnostic methods". Touches 5 of 182 rows. |
| Shedding ascertainment (Oxford) | `eta(D)` treatment-truncation correction (Tier 2) | plan §2.7; §5.7.3. Grouped Oxford rows only — **not** the Darton per-subject rows. |
| Unobserved Maryland immunity | Two-component latent mixture (`pi_susc`, `CoP_susc`, `CoP_imm`), taken over the cascade **product** | plan §2.4; §5.7.3 |
| Within-arm titre variation (grouped Oxford rows) | Arm GMT plugged into the nonlinear response | §5.7.5. A Jensen approximation, small at the retained fit. A quadrature plan is in `calibration/jin_within_arm_titre_plan.md` (DRAFT). |
| Extra-binomial variation | One shared parameter, beta-binomial `grand_overdispersion_rho` | `calibration/cohort_random_effects_design.md` (LOCKED). |
| Dose uncertainty | Published ranges collapsed to point estimates | Not integrated over. Two Waddington rows carry doses that differ from their own reported medians (§7). Open. |
| Chronic carriage, seroconversion, death, incubation period | Not modelled | No likelihood term exists; §1.2/§1.3 are aspirational for these |

**No-double-counting rule (package-wide).** Where a single cohort of volunteers
supplies two endpoints, they are NOT two independent marginals — they are factorized
into a marginal and a conditional on the same men:

- Hornick 10^7: `H-I-7` (28/30 infected) + `H-FgI-7` (16/28 fever|infected). The
  marginal fever row `H-F-7` is excluded to avoid double-counting.
- Darton, all three arms: groups 6/7, infection over all 90 subjects +
  fever|infection over the 63 infected.

**Where the rule is not applied.** The canonical `t2-indiv-vax` fit contains **ten**
cohorts supplying two independent marginals over the same volunteers, not one:

- **Levine, 4 pairs** — `Lev-F-k` and `Lev-I-k`, k = 1..4.
- **Oxford, 5 pairs** — grouped shedding rows for
  `W-I-3`/`W-F-3`, `W-I-4`/`W-F-4`, and Jin's `J-I-*`/`J-F-*` for all three arms.
  Each pair is the same men counted twice.
- **Gilman, 1 partial overlap** — `Gil-I-ctrl` (43 of 64) partly overlaps the three
  `Gil-F-*` stratum rows.

**Known limitation, tolerated:** the clean fix needs fever nested inside the
infection endpoint, but these infection endpoints are stool-based and fever is NOT
nested in stool-positivity (Darton: 7 of 20 TD+ were stool-negative), and none of
these studies publishes a cross-tabulation. Mitigated only in LOO, where
`compute_loo_units()` keys grouped rows by `cohort_id` and so merges every one of
these pairs into a single unit.

**CoP caveat — the axis is known to be wrong for the Maryland era.** The scalar CoP is
anti-Vi IgG. In the Maryland volunteers, every paper that looked found baseline anti-Vi
did **not** predict outcome (Gilman 1977 p.721; Levine 1976 p.427), and purified Vi
vaccine gave a 13-fold antibody rise with zero protection (Woodward 1980 p.554), while
**H agglutinin did** stratify. The Maryland latent mixture is therefore an
anti-Vi-*equivalent* device standing in for a non-Vi immunity axis, not a measurement
of one. Documented at `calibration/priors.yaml` (`pi_susc`).

**Fever severity.** `phi(T,D)` is the implemented severity model (§1.3). Known
limitation: the dose-lift term carries no T dependence, so `phi -> 1` uniformly in T as
dose rises — at Hornick's top dose (10^9 milk) the model implies `phi(41C) = 0.94`.
It represents "threshold choice stops mattering at high dose", not "higher dose gives
hotter fevers". Consequence of pinning `beta_phi = 1`.

---

## 7. Mapping to Available Data

Verified against `calibration/dose_response_data.csv` and
`analysis_data/darton_individual_endpoints.csv`, 2026-07-31.

**Scale of the canonical fit.** 37 CSV rows: 25 `tier1_active`, 31 `tier2_active`,
6 `validation_only`. After Darton individualization, `t2-indiv-vax` presents
**182 observations** to Stan:

| | group | rows |
|---|---|---|
| 1 | `ox_fev` | 6 |
| 2 | `ox_inf` | 5 |
| 3 | `md_fev` | 11 |
| 4 | `md_inf` | 6 |
| 5 | `hornick_cond` | 1 |
| | **grouped subtotal** | **29** |
| 6 | `ox_inf_indiv` | 90 |
| 7 | `ox_fevginf_indiv` | 63 |
| | **per-subject subtotal** | **153** |

Plus 3 ladder binomials and 1 cross-tab binomial in `log_lik` (§5.7.6). For
reference, `t1-indiv` presented 80 (24 grouped + 56 per-subject); the older figure
appears in Tier 1 documents and is not this fit.

| Reference component (§1.2–1.3) | Data available | Quality | Notes |
|---|---|---|---|
| **Dose** | All 37 rows | Variable | 10^3–10^9 CFU overall, but the eras barely overlap: Maryland 10^3–10^9 (milk), Oxford {10^3, 10^4, 1.82x10^4, 2x10^4, 2.5x10^4} (bicarbonate). The `delta` bridge is what joins them. Published dose *ranges* are collapsed to point estimates, and two Waddington rows use round numbers against their own reported medians (`W-F-3`: 1000 vs median 1.34x10^3; `W-F-4`: 10000 vs median 1.98x10^4). |
| **Pre-challenge immunity** | 13 of 37 rows carry a `CoP`; 24 are NA | Oxford only | The 24 NA rows are **22 Maryland** (no 1960s serology exists) **plus 2 `validation_only` Oxford rows** (`D-F-Ty21a`, `D-F-M01`). Maryland immunity is latent (the mixture). Darton contributes **90** per-subject titres across all three arms, but 64 of the 90 sit at exactly the detection floor (12 of 30 Placebo are above it) — a thin titre range, which is the core identifiability problem for `gamma_inf` (§8.5). |
| **Infection** | 12 rows (6 Maryland, 6 Oxford) | Definition-heterogeneous | The endpoint differs materially by study — Hornick stool-or-blood culture, Levine any-time stool, Gilman late shedding 4–30 d, Oxford stool shedding. The fitted `psi` map (§5.7.3) accounts for these differences. Of the 6 Oxford infection rows, **5 enter the fitted likelihood**: `D-I-plac` is dropped by individualization as a double count of the 30 Darton placebo `ox_inf_indiv` rows. |
| **Bacteremia** | Darton per-subject only | Good | Placebo: 20 of 30 bacteremic, 20 blood-positive. Not a separate likelihood layer — folded into `bact_or_stool` as the broad infection marker (Placebo 26/30, M01ZH09 21/31, Ty21a 16/29). No Maryland bacteremia counts. |
| **Stool shedding** | Oxford (Tier 2) + Levine + Gilman + Darton | Good where present | Darton measures it against the broad marker in the same men: 15 nested stool+ of 26 `bact_or_stool+` — the fitted anchor for `psi_stool`. The marginal total is 19 stool+ of 30; fever is **not** nested in stool-positivity: 7 of 20 TD+ were stool-negative. |
| **Fever** | 24 rows (15 Maryland, 9 Oxford) | Good | Threshold definitions vary and are modelled by `phi(T,D)`: Hornick 39.4C, Levine 38.3C, Gilman assumed 38.3C (imported, see `data_prep.R`), Oxford composite TD at T_ref = 38.0C. |
| **Fever severity** | Darton ladder only | Good but thin | Among 20 TD+ placebo subjects the extract records 20/19/16/10/8 crossing >=37/37.5/38/38.5/39 C. **Only the top three rungs are in the likelihood** (16/10/8 at 38.0/38.5/39.0) — see §5.7.6 for why. Pins `phi0(T)`; **39.0 C is the maximum**, so Hornick's 39.4 C is extrapolation. |
| **Fever given infection** | 1 Maryland row + 63 Darton subjects | Improved at Tier 2 | `H-FgI-7` (16/28) is the only grouped conditional. The per-subject layer is now 26 Placebo + 21 M01ZH09 + 16 Ty21a = 63 rows, the second-largest group in the fit — it is no longer the thinnest layer, though it rests on one cohort. |
| **Clinical typhoid (composite)** | Oxford TD | Good | Taken as the reference definition, so `phi == 1` for Oxford by construction. |
| **Chronic carriage** | None | — | No likelihood term. |
| **Seroconversion** | Post-challenge only | Not usable as an outcome | Darton has baseline/pre-challenge titres; post-challenge seroconversion is not modelled. |
| **Death / severe disease** | None | — | No deaths in the Maryland program (Woodward 1980). No likelihood term. |
| **Incubation period** | Available, unused | Good | Hornick Table 1 and Waddington Table 3 report time-to-event. Not modelled; could constrain `alpha` and the dose-response shape. |

**What the mapping implies.** Three reference components carry most of the model's
weight (dose, fever, infection). Fever severity rests entirely on one Darton cohort
of 30 via three ladder rungs, all Maryland immunity is latent, and the high-titre
end of the CoP axis is carried by three grouped Jin arms against 26
titre-informative Darton subjects. Carriage, seroconversion, death and incubation
have no likelihood term at all — §1.2 and §1.3 remain aspirational for those.

---

## 8. Calibration Setup

Implemented in `calibration/typhoid_dose_response.stan`; Tier 2 design and locked
decisions in `joint_inference_plan.md`. This section is kept in sync with
that plan — if they disagree, the plan is authoritative and this is stale.

### 8.1 Likelihood Structure

Full expressions in §5.7. This table is the index.

| Component | Form | Notes |
|---|---|---|
| Dose-outcome counts | Beta-binomial, **every row** | `y ~ beta_binomial(n, p*k, (1-p)*k)`, `k = (1-rho)/rho`. `rho = 0` is exactly the binomial. There is **no `n = 1` branch** — the code relies on `beta_binomial(1, pk, (1-p)k) == bernoulli(p)` for any `k`. A reader looking for the special case will not find one. |
| Observation dispatch | Seven groups | §5.7.4. Codes 1–7; `obs_prob()` has no other branches and is the single source shared by the model block and generated quantities. |
| Non-anti-Vi vaccine protection | `V_M01ZH09`, `V_Ty21a` dividing the exponent | Groups 6/7 only. A vaccine arm's `V` therefore enters the composite fever probability **twice**, once per cascade stage. |
| Fever definition | Fitted map `phi(T,D)` | Not stratification — a measurement model. `phi0(T)` is pinned by the Darton ladder; the **dose-lift half** is pinned by `beta_phi = 1` plus the naive cascade curve, which is fit to Hornick's multi-dose fever rows. "Pinned by the ladder" covers only `phi0`. |
| Infection definition | Fitted map `psi` (Tier 2) | Pinned by the nested Darton stool-vs-`bact_or_stool` cross-tab, **15/26**; the marginal 19/26 is not the fitted likelihood. Applies to Maryland infection rows only — 5 of 182. |
| Shedding ascertainment | `eta(D)` (Tier 2, grouped Oxford) | Treatment truncation; distinct mechanism from `psi`. **Not** applied to the Darton per-subject rows, whose endpoint is the broad `bact_or_stool` marker. |
| Decoupled sub-likelihoods | 3 ladder binomials + 1 cross-tab binomial | Deliberately **not** overdispersed. `log_lik` is `N_obs + N_ladder + psi_active` = 186 so that `target == lprior + sum(log_lik)` exactly; LOO uses only the first `N_obs`. |
| Numerical guards | `p` clamped to `[1e-12, 1-1e-12]` | Required by the beta-binomial's positive-shape constraint. The group-5 conditional ratio is clamped separately inside `obs_prob()`. |
| Prior-predictive mode | `prior_only` data flag | Skips the whole likelihood; the paired run is `..._phi-rho-eta-psi-vax-prior/`. |
| Dose uncertainty | **Point estimates** | Oxford reports ranges; NOT integrated over. Open. |
| Within-arm titre variation | **Arm GMT plug-in** | Grouped Oxford rows. A Jensen approximation, <1% at the retained posterior; quadrature replacement is DRAFT. §5.7.5. |
| Missing immunity | Latent 2-component mixture | Maryland only; taken over the cascade **product**, not per stage. |

### 8.2 Study-Level Effects

| Effect | Structure | Rationale |
|---|---|---|
| Strain | Excluded | Quailes throughout the calibration set |
| Delivery medium | Fixed scalar `delta` | Milk -> bicarbonate-equivalent dose |
| Era | Not modelled separately | Constraint: nothing may be given per-cohort freedom on `N50`. Hornick's cohorts ARE the dose ladder, so such a term flattens the dose-response while improving fit. |
| Cohort / study | Single shared overdispersion, `grand_overdispersion_rho` | The fit has 29 grouped observations over 16 cohorts, with **1** singleton. See `calibration/cohort_random_effects_design.md` (LOCKED). |
| Outcome definition | Measurement model, not stratification | `phi(T,D)` for fever, `psi` for infection |
| Immunology across eras | **Assumed invariant** | `gamma_inf` / `gamma_fevginf` are SHARED between Maryland and Oxford by design: human immunology is taken to be the same in both. This is why an era-specific protection scale is not an available fix for Maryland misfit. |
| What else is shared | Nearly everything | `alpha_inf`, `alpha_fevginf`, `N50_inf`, `N50_fevginf`, `phi0_a`, `phi0_b` and `rho` are also fully shared across eras. The **only** era-specific machinery is `delta`, the Maryland mixture (`pi_susc`, `CoP_susc`, `CoP_imm`), and `phi`. |

### 8.3 Individual-Level Nuisances

| Nuisance | Treatment | Notes |
|---|---|---|
| Host heterogeneity | Absorbed into `alpha` | Per §6.2 |
| Age / sex | Not used | Unavailable for most Maryland rows. Glynn 1995 finds age <30 RR 1.79 in this population, so this is a real unmodelled confound. |
| Prior exposure | Latent mixture (Maryland); measured CoP per subject (Darton) or per arm GMT (Jin/Waddington/Gibani) | Maryland `pi_susc` = 0.670 [0.42–0.86] in the retained fit. It is numerically close to Gilman's H-negative fraction 36/53 = 0.68 and Woodward's non-veteran fraction 200/305 = 0.66 — but `priors.yaml` warns this is **near-circular** and is not evidence the mixture weight is identified. Do not cite the agreement as corroboration. |
| Cohort membership | Recorded as `cohort_id` | Not passed to Stan; no likelihood term reads it. Post-fit it keys the LOO units, so model comparison does not count the same volunteers twice. `compute_loo_units()` merges **all** same-cohort grouped rows — the four Levine pairs, the Gilman control cohort, and (at Tier 2) the five Oxford fever/shedding pairs — not just the Hornick marginal + conditional. |

### 8.4 Prior Specification

Single source of truth: `calibration/priors.yaml` (hyperparameters are Stan *data*).

All 20 priors passed as Stan data, complete as of 2026-08-01:

| Parameter | Prior | Rationale |
|---|---|---|
| `log10_N50_inf` | Normal(2.5, 1.0) | Order-of-magnitude uncertainty (~300 bicarb CFU) |
| `log10_N50_fevginf` | Normal(2.8, 1.0) | ~600 bicarb CFU; **applied to the derived quantity**, see `d_fev` |
| `d_fev` | half-Normal(0, 1) | The ordering offset (§5.7.2). The three N50 terms retain their joint prior because the change of variables has unit Jacobian |
| `alpha_inf`, `alpha_fevginf` | LogNormal(-1.5, 0.8) | Positive, weakly informative (~0.1–0.5) |
| `gamma_inf`, `gamma_fevginf` | LogNormal(-1.6, **0.7**) | Median ~0.2. Anchored to **two concordant slopes** — Jin 2017 OR 0.37 and Darton HR 0.29, both per log10 anti-Vi. `sdlog` was tightened 0.9 → 0.7 when the second anchor was added |
| `log10_delta` | Normal(3.5, 0.7) | Milk-to-bicarb bridge, ~1000–30000x |
| `pi_susc` | Beta(7, 4) | ~0.65 (Gilman 36/53, H-agglutinin). Circularity caveat in `priors.yaml` and §8.3 |
| `CoP_susc` | LogNormal(0, 0.3) | ~1x naive. **A free parameter, not fixed at 1** |
| `CoP_imm` | Exponential(0.074) | Mean 13.5 (= 50 EU/mL absolute). **Prior-carried** — no 1960s serology exists |
| `phi0_a` | Normal(1.0, 1.5) | logit `phi0` at `T_ref` = 38.0 C |
| `phi0_b` | half-Normal(0, 2.0) | logit decay per degC; the half-normal makes `phi0` monotone-decreasing in threshold **structurally** |
| `eta_lo` | Beta(5, 5) | ~0.5; the high-dose asymptote of shedding detection |
| `kappa` | LogNormal(0, 1.0) | `eta`'s dose scaling. **The most prior-dominated parameter in the retained fit** (§8.5) |
| `psi_stool` | Beta(3, 2) | ~0.60, weakly informative |
| `frac_late` | Beta(4, 1.5) | ~0.73, vague. Makes `psi_late <= psi_stool` structural |
| `grand_overdispersion_rho` | Beta(1, 49) | ICC; `rho = 0` is the binomial. Anchored to sigma ~ 0.26 measured across the five Maryland 10^5 cohorts |
| `log_V_M01ZH09`, `log_V_Ty21a` | Normal(0, 1) | Log scale, so `V = 1` (no residual effect beyond the subject's own titre) is the prior centre and the data can move it either direction |

### 8.5 Identifiability Concerns

Measured in the retained Tier 2 result, not hypothetical; the corresponding
figures and full diagnostics are in
`calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/`:

| Issue | Status | Diagnostic / resolution |
|---|---|---|
| `N50` vs `delta` | **Trade off, not confounded**: r = **-0.54** (inf), **-0.59** (fev\|inf) | Neither pair reaches the \|r\| > 0.7 threshold, so neither appears in the generated summary's correlation table. Marginals contract ~2x from prior, and the weakly identified direction contracts 1.7x. `delta_bridge.png`. |
| `CoP_imm` | **Prior-sensitive / weakly identified** (priorsense prior 0.425 vs likelihood 0.082) | No Maryland serology; report as an assumption |
| `kappa` | **The most prior-dominated parameter in the fit** (prior 0.631 vs likelihood 0.050) | `eta`'s modelled decline happens below ~100 CFU, under every observed Oxford dose, so the data cannot see the shape. Read `eta` as a constant detection loss |
| `pi_susc`, `CoP_susc` | Also flagged "strong prior / weak likelihood" (0.098 vs 0.014; 0.159 vs 0.035) | Both are latent Maryland mixture quantities with no serology behind them |
| `gamma_inf` vs `gamma_fevginf` | Separated in the retained fit (medians 0.081 and 0.192; r = -0.27) but reflect a **contradiction**, not merely thin data | 64 of the 90 Darton per-subject rows sit at exactly CoP = 1 and provide no titre-slope information. The other 26 give an arm-adjusted infection slope of the **opposite sign** (+0.71 per log10, z = +1.25) and a fever\|inf slope near zero (-0.08, z = -0.14). Jin's grouped arms provide the high-titre end but give ~9%/10x on shedding. See `tier2_findings_report.md` §5 |
| Grouped vs individualized weighting | **Asymmetric by construction** | At n = 1 the beta-binomial is exactly binomial, so the 153 per-subject rows neither inform `rho` nor receive its discount, while the 29 grouped rows do both. At `rho` = 0.024 the design effect discounts Jin's arms ~1.8x. |
| `psi_stool` vs `eta` | Confounded, but **not** because both multiply `P_inf` | `psi` applies only to Maryland infection rows; `eta` only to grouped Oxford rows. They never touch the same row. The confound lives entirely in the cross-tab sub-likelihood, where `psi_stool * eta(18200)` is *itself* the binomial probability. Only `eta`'s dose dependence separates them — and that is `kappa`, which is prior-dominated |
| `frac_late` (Gilman 4-30 d) | Expected prior-dominated; **it is not** (prior 0.08 vs likelihood 0.20) | Gilman's single `psi_def = 2` row constrains it more than the plan assumed. Structurally constrains `psi_late = psi_stool * frac_late` |
| Cohort effects vs dose-response | Held by using one shared overdispersion parameter | Per-cohort offsets would compete with `N50_inf`/`alpha_inf` |
| `sigma_study` | Not included | The shared beta-binomial `grand_overdispersion_rho` models extra-binomial variation. |

### 8.6 Gates

The Tier 2 fit uses the following gates. It reports parameter movement, priorsense,
and diagnostics; the single remaining divergence is a residual risk:
1. `log10_N50_inf` and `alpha_inf` must not move materially (the signal-absorption gate).
2. priorsense on the new parameter; if prior-dominated, report it as such.
3. LOO must improve, using `compute_loo_units()` grouping, which keys grouped rows
   by `cohort_id` and so merges every same-cohort pair (Hornick marginal +
   conditional, the four Levine pairs, the Gilman control cohort, and at Tier 2 the
   five Oxford fever/shedding pairs). **Not evidenced for the Tier 2 increment:** no
   `loo.json` exists for any `t2-*` run. The only LOO artifacts in the repo are
   under `results/scenarios/t1-indiv/`. `[open — gate asserted, evidence absent]`
4. The `Gil-F-Hlo` residual is **not** a gate — it is a within-cohort stratum contrast
   and is expected to persist. If it vanishes, something else moved; explain it.
