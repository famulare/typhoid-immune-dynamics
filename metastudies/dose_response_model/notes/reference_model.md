# Reference Model for Typhoid Dose-Response

**Purpose**: Define the most elaborate model the literature supports—what we would fit with optimal data and unlimited precision. This serves as a benchmark for reasoning clearly about which simplifications are forced by data limitations versus chosen for parsimony.

**Status**: Draft - to be refined during Phase 3

---

## 1. Latent Biological Processes

The reference model distinguishes these underlying processes:

### 1.1 Infection Cascade

```
Ingestion → Gastric survival → Colonization → {Systemic invasion, Intestinal shedding} → {Clinical disease, Chronic carriage}
```

- Bacterial survival through gastric barrier → colonization
- Colonization → systemic invasion (bacteremia)
- Colonization → intestinal shedding
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
| **Seroconversion** | Antibody response | Serology (multiple assays) |

### 1.3 Time-to-Event Structure

For each outcome, model:
- Time to bacteremia onset
- Time to fever onset
- Duration of fever
- Duration of shedding
- Time to seroconversion
- Time to chronic carriage establishment

### 1.4 Severity Gradations

- Fever severity (temperature magnitude, duration)
- Disease severity (mild/moderate/severe; hospitalization; complications)
- Death

### 1.5 Joint/Conditional Relationships

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

**Intestinal shedding** (binary/duration):
$$P(\text{shedding} | \text{colonization}, D_{\text{gastric}}, \text{strain}, \text{immunity}, \text{host})$$

*Design note*: $D_{\text{gastric}}$ drops out of acute disease and chronic carriage conditioning. This encodes a memorylessness hypothesis: from systemic invasion forward, outcomes depend on the invasion state, not the original dose. 

#### Observables

Biological states accessible from outside the body (without histopathology):

**Stool shedding**:
$$P(\text{stool shedding} | \text{colonization}, D_{\text{gastric}}, \text{strain}, \text{immunity}, \text{host})$$

**Bacteremia**:
$$P(\text{bacteremia} | \text{systemic invasion}, \text{immunity}, \text{host})$$

**Fever** (ordinal, threshold-dependent):
$$P(\text{fever} | \text{acute disease}, \text{immunity}, \text{host})$$

**Death**:
$$P(\text{death} | \text{acute disease}, \text{immunity}, \text{host})$$

**CoP (Correlate of Protection)**: The true immune state—an array of cellular/molecular responses. See Section 2 for full treatment. Pre-challenge CoP is the "immunity" appearing in conditioning throughout this document.

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

**Serology titer** (observation of CoP component):
$$P(\text{titer} | \text{CoP}, \text{assay})$$

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

### 2.1 CoP as the True Immune State (Observable)

The **Correlate of Protection (CoP)** is the true immune state—in principle, an array of cellular and molecular responses:

$$\vec{\text{CoP}} = (\text{CoP}_1, \text{CoP}_2, \ldots, \text{CoP}_n)$$

This array includes components from:

**Immune compartments**:
- **Humoral immunity**: Circulating antibodies (IgG, IgA) against various antigens (Vi, O, H, LPS)
- **Mucosal immunity**: Secretory IgA, gut-resident memory
- **Cellular immunity**: T cell responses (CD4+, CD8+), memory T cells

The pre-challenge CoP array is the "immunity" appearing in conditioning throughout Section 1.

### 2.2 Titers as Observations of CoP

Each serological assay measures one (or a combination of) CoP components:

$$P(\text{titer} | \vec{\text{CoP}}, \text{assay})$$

For example, anti-Vi IgG ELISA (Vacczyme assay) observes a specific humoral component.

**Seroconversion** is a derived quantity comparing two titer observations:
$$P(\text{seroconversion} | \text{titer}_{\text{pre}}, \text{titer}_{\text{post}}, \text{assay})$$

Assay variance determines whether an observed change is classified as "real" seroconversion.

### 2.3 Post-Challenge CoP Dynamics

The immune state evolves in response to infection:
$$\vec{\text{CoP}}_{\text{post}} = f(\vec{\text{CoP}}_{\text{pre}}, \text{colonization}, \text{systemic invasion}, \text{acute disease}, \ldots)$$

*Scope note*: Pre-post dynamics are important to typhoid modeling but are treated as a separate data and calibration project in this repository. For this project, if seroconversion is used, we assume a lognormally distributed difference with standard error of approximately ±1 log2 unit.

### 2.4 Simplification Pathway

For the working model, the CoP array collapses to a scalar:

$$\text{CoP} = g(\vec{\text{CoP}})$$

This collapse to a scalar is implicit in the concept of a "correlate of protection"—the CoP framework assumes a low-dimensional summary sufficient to predict protection.

**Key simplification assumptions** (to be documented in Phase 5):
- Which immune compartments dominate protection?
- How do different compartments combine (additive? multiplicative? minimum?)?
- What is lost by ignoring dynamics within a challenge study timescale?

---

## 3. Dose-Response Framework

### 3.1 General Form

For each outcome $Y$:

$$
P(Y | D, \text{CoP}, \vec{\theta}_Y) = f_Y(D, \text{CoP}; \vec{\theta}_Y)
$$

Where:
- $D$ = bacterial dose (CFU)
- $\text{CoP}$ = scalar correlate of protection (see Section 2)
- $\vec{\theta}_Y$ = outcome-specific parameters

### 3.2 Hierarchy of Outcomes

Outcomes decompose along the causal chain defined in Section 1.5:

$$P(\text{acute disease}) = P(\text{acute disease} | \text{systemic invasion}) \cdot P(\text{systemic invasion} | \text{colonization}, D_{\text{gastric}}) \cdot P(\text{colonization} | D_{\text{gastric}}) \cdot P(D_{\text{gastric}} | D_{\text{challenge}})$$

(Conditioning on strain, CoP, host suppressed for readability.)

Similar decompositions for:
- Bacteremia given systemic invasion
- Fever given acute disease
- Shedding given colonization and $D_{\text{gastric}}$
- Chronic carriage given systemic invasion
- Clinical typhoid given acute disease and observations

### 3.3 Mechanistic Basis

The working model uses a modified beta-Poisson framework:

$$
P(\text{infection} | D, \text{CoP}) = 1 - \left(1 + D \cdot \frac{2^{1/\alpha} - 1}{N_{50}}\right)^{-\alpha / \text{CoP}^{\gamma}}
$$

This form has mechanistic interpretation:
- Each bacterium has independent probability of initiating infection
- $N_{50}$ = dose for 50% probability of outcome in naive individuals
- $\alpha$ = unexplained heterogeneity of outcome given dose (aggregates all sources of variation in the cascade from challenge dose to outcome)
- $\gamma$ = how strongly immunity scales effective dose

---

## 4. Observational Model

For each latent process, specify the probability of observing a positive result given the true state:

$$
P(\text{observed} | \text{latent state}, \text{ascertainment protocol})
$$

### 4.1 Measurement Processes

| Observable | Model |
|------------|-------|
| **Blood culture** | $P(\text{culture}^+ \mid \text{bacteremia}, \text{timing}, \text{antibiotic treatment})$ |
| **PCR** | $P(\text{PCR}^+ \mid \text{bacteremia}, \text{timing})$ |
| **Stool culture** | $P(\text{stool}^+ \mid \text{shedding}, \text{sampling schedule})$ |
| **Temperature** | $P(T > \theta \mid \text{fever state}, \text{measurement frequency})$ |
| **Serology** | $P(\text{titer} > \text{threshold} \mid \text{seroconversion}, \text{assay}, \text{timing})$ |

### 4.2 Key Considerations

- **Sensitivity vs specificity** for each assay
- **Timing relative to infection** affects detection probability
- **Threshold definitions** vary across studies (e.g., fever = 100°F vs 100.4°F)
- **Sampling frequency** affects time-to-event precision

Note: Functional forms left unspecified at this stage; the reference model establishes *what* needs modeling, not *how*.

---

## 5. Heterogeneity Sources

### 5.1 Study-Level

| Source | Impact |
|--------|--------|
| **Strain** | Quailes vs Ty2 vs wild-type may differ in virulence |
| **Delivery medium** | Milk vs bicarbonate affects gastric survival |
| **Fasting protocol** | Fasting duration before challenge affects gastric survival |
| **Subject population** | Naive volunteers vs endemic-area residents |
| **Era** | 1960s vs 2010s protocols, ascertainment methods |
| **Outcome definition** | Clinical typhoid criteria vary (e.g., ">103°F/36hr" vs "≥38°C/12hr OR culture+") |
| **Diagnostic methods** | Culture only vs culture + PCR affects infection ascertainment |

### 5.2 Individual-Level

- Age, sex
- Prior exposure history
- Baseline immune state
- Genetic susceptibility (e.g., HLA type)

---

## 6. DAG: Latent Processes → Observables → Observations

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
          ├─▶│  Systemic   │          │  Intestinal │
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

## 7. Simplifications Needed for Practical Model

*To be completed after reviewing available data*

| Reference Model Component | Simplification | Justification |
|--------------------------|----------------|---------------|
| Multiple immune compartments | Single scalar CoP | TBD |
| Time-to-event | Binary outcome | TBD |
| Joint outcome distribution | Independent outcomes | TBD |
| ... | ... | ... |

---

## 8. Mapping to Available Data

*To be completed after Phase 4 joint review*

| Reference Model Component | Data Available | Data Quality | Notes |
|--------------------------|----------------|--------------|-------|
| Dose | Yes (most studies) | Variable | Units/calibration issues |
| Pre-challenge immunity | Partial | Variable | Often missing or indirect |
| Bacteremia | Yes (most studies) | Good | Culture timing varies |
| Fever | Yes | Good | Threshold definitions vary |
| ... | ... | ... | ... |
