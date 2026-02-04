# Reference Model for Typhoid Dose-Response

**Purpose**: Define the most elaborate model the literature supports—what we would fit with optimal data and unlimited precision. This serves as a benchmark for reasoning clearly about which simplifications are forced by data limitations versus chosen for parsimony.

**Status**: Draft - to be refined during Phase 3

---

## 1. Latent Biological Processes

The reference model distinguishes these underlying processes:

### 1.1 Infection Cascade

```
Ingestion → Gastric survival → Colonization → {Systemic invasion, Intestinal shedding} → Clinical disease
```

- Bacterial survival through gastric barrier → colonization
- Colonization → systemic invasion (bacteremia)
- Colonization → intestinal shedding
- Systemic invasion → clinical disease

### 1.2 Outcomes to Model

Each outcome is a distinct random variable:

| Outcome | Latent Process | Observable Via |
|---------|---------------|----------------|
| **Infection** | Successful colonization | Latent, partially observed through downstream markers |
| **Bacteremia** | Bacteria in bloodstream | Blood culture, PCR |
| **Stool shedding** | Intestinal colonization with fecal excretion | Stool culture |
| **Fever** | Elevated temperature | Thermometry (multiple threshold definitions) |
| **Clinical typhoid** | Composite syndrome | Clinical diagnosis (definition varies by study/era) |
| **Seroconversion** | Antibody response | Serology (multiple assays) |

### 1.3 Time-to-Event Structure

For each outcome, model:
- Time to bacteremia onset
- Time to fever onset
- Duration of fever
- Duration of shedding
- Time to seroconversion

### 1.4 Severity Gradations

- Fever severity (temperature magnitude, duration)
- Disease severity (mild/moderate/severe; hospitalization; complications)

### 1.5 Joint/Conditional Relationships

Outcomes are not independent. Key relationships:

- $P(\text{fever} | \text{bacteremia}, D, \text{immunity})$
- $P(\text{shedding} | \text{infection}, D, \text{immunity})$
- $P(\text{seroconversion} | \text{infection}, \text{bacteremia}, \text{fever})$
- Correlation structure among outcomes (e.g., bacteremia duration and fever severity)

---

## 2. Immunity Representation

### 2.1 Full Mechanistic Representation

What we're simplifying from:

**Immune compartments**:
- **Humoral immunity**: Circulating antibodies (IgG, IgA) against various antigens (Vi, O, H, LPS)
- **Mucosal immunity**: Secretory IgA, gut-resident memory
- **Cellular immunity**: T cell responses (CD4+, CD8+), memory T cells

**Dynamics**:
- Pre-existing immunity from prior exposure/vaccination
- Waning over time
- Boosting upon re-exposure
- Cross-reactivity with related organisms

### 2.2 Simplification Pathway

These compartments collapse to a scalar Correlate of Protection (CoP) that modulates susceptibility:

$$
\text{CoP} = g(\vec{I})
$$

Where $\vec{I}$ is the full vector of immune state variables. The reference model makes explicit what biological complexity is being summarized by this scalar.

**Key simplification assumptions** (to be documented):
- Which immune compartments dominate protection?
- How do different compartments combine (additive? multiplicative? minimum?)?
- What is lost by ignoring dynamics within a challenge study timescale?

---

## 3. Dose-Response Framework

### 3.1 General Form

For each outcome $Y$ (infection, fever, bacteremia, etc.):

$$
P(Y | D, \vec{I}, \vec{\theta}_Y) = f_Y(D, g(\vec{I}); \vec{\theta}_Y)
$$

Where:
- $D$ = bacterial dose (CFU)
- $\vec{I}$ = vector of immune state variables
- $g(\vec{I})$ = function mapping immune state to effective protection (simplifies to scalar CoP)
- $\vec{\theta}_Y$ = outcome-specific parameters

### 3.2 Hierarchy of Outcomes

Some outcomes are conditional on others:

$$
P(\text{bacteremia} | D, \vec{I}) = P(\text{bacteremia} | \text{infection}) \cdot P(\text{infection} | D, \vec{I})
$$

Similar decompositions for:
- Fever given bacteremia
- Shedding given infection
- Seroconversion given infection/bacteremia/fever
- Clinical typhoid given component symptoms

### 3.3 Mechanistic Basis

The working model uses a modified beta-Poisson framework:

$$
P(\text{infection} | D, \text{CoP}) = 1 - \left(1 + D \cdot \frac{2^{1/\alpha} - 1}{N_{50}}\right)^{-\alpha / \text{CoP}^{\gamma}}
$$

This form has mechanistic interpretation:
- Each bacterium has independent probability of initiating infection
- $N_{50}$ = dose for 50% infection in naive individuals
- $\alpha$ = heterogeneity in bacterial "hit" probability
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
| **Subject population** | Naive volunteers vs endemic-area residents |
| **Era** | 1960s vs 2010s protocols, ascertainment, definitions |

### 5.2 Individual-Level

- Age, sex
- Prior exposure history
- Baseline immune state
- Genetic susceptibility (e.g., HLA type)

---

## 6. DAG: Latent Processes → Observables

```
                        ┌─────────────┐
                        │    Dose     │
                        └──────┬──────┘
                               │
                               ▼
┌─────────────┐         ┌─────────────┐
│  Immunity   │────────▶│  Infection  │
│   (vec I)   │         │  (latent)   │
└─────────────┘         └──────┬──────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
       ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
       │ Bacteremia  │  │  Shedding   │  │Seroconversion│
       │  (latent)   │  │  (latent)   │  │   (latent)   │
       └──────┬──────┘  └──────┬──────┘  └──────┬──────┘
              │                │                │
              ▼                ▼                ▼
       ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
       │Blood culture│  │Stool culture│  │  Serology   │
       │    PCR      │  │             │  │   titer     │
       └──────┬──────┘  └─────────────┘  └─────────────┘
              │
              ▼
       ┌─────────────┐
       │   Fever     │
       │  (latent)   │
       └──────┬──────┘
              │
              ▼
       ┌─────────────┐
       │ Temperature │
       │ measurement │
       └─────────────┘
```

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
