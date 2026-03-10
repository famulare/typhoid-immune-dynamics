# Anti-Vi Metrology Bridge Specification

**Purpose**: Define the reference bridging model—the most elaborate model the literature supports, what we would fit with optimal data and unlimited precision. This serves as a benchmark for reasoning clearly about which simplifications are forced by data limitations versus chosen for parsimony.

**Status**: Phase 3 pending. Sections 5-7 to be filled during Phases 4-5.

**Related documents**:
- `assay_mapping.md` - Decision rules for mapping reported assays to canonical nodes
- `cross_cutting_observations.md` - Patterns across the extraction corpus

---

## 1. Latent Immunity Construct

The reference model posits a latent "true immunity signal" that all assays attempt to measure:

### 1.1 What Are We Trying to Measure?

All Vi serology assays attempt to quantify the immune response to the Vi capsular polysaccharide antigen. The underlying biological signal includes:

- **Circulating anti-Vi IgG**: Primary target of modern ELISAs
- **Circulating anti-Vi IgM/IgA**: May contribute to HA, less to IgG-specific ELISA
- **Antibody avidity**: Affects agglutination efficiency
- **Total antibody vs specific isotype**: HA measures functional agglutination; ELISA measures specific binding

### 1.2 Latent Variable Representation

$$\text{Immunity}_{\text{latent}} = f(\text{IgG}, \text{IgM}, \text{IgA}, \text{avidity}, \ldots)$$

In practice, this collapses to a scalar for modeling purposes. The key question is whether different assays observe the same scalar or systematically different projections.

### 1.3 Population Heterogeneity in the Latent Signal

The relationship between latent immunity and observed titers may vary by:

| Population | Expected Characteristics |
|------------|-------------------------|
| **Acute typhoid** | Recent exposure; rising titers; mixed isotype response |
| **Chronic carriers** | Long-term Vi exposure; potentially high avidity; stable IgG |
| **Vaccinated** | Defined antigen exposure; predominantly IgG response |
| **Endemic baseline** | Repeated environmental exposure; variable |
| **Non-endemic naive** | Low or undetectable; near assay LOD |

---

## 2. Assay Graph: Nodes and Edges

### 2.1 Canonical Assay Nodes

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           REFERENCE SCALE                                    │
│                                                                              │
│                        ┌──────────────────┐                                  │
│                        │   IS-relative    │                                  │
│                        │  (NIBSC 16/138)  │                                  │
│                        └────────┬─────────┘                                  │
│                                 │                                            │
│                    ┌────────────┼────────────┐                              │
│                    │            │            │                              │
│                    ▼            ▼            ▼                              │
│            ┌───────────┐ ┌───────────┐ ┌───────────┐                        │
│            │ VaccZyme  │ │ In-house  │ │ Vi-IgGR1  │                        │
│            │  ELISA    │ │ IS-ELISA  │ │   2011    │                        │
│            └─────┬─────┘ └───────────┘ └───────────┘                        │
│                  │                                                           │
├──────────────────┼───────────────────────────────────────────────────────────┤
│                  │          MODERN ELISA (2000s+)                            │
├──────────────────┼───────────────────────────────────────────────────────────┤
│                  │                                                           │
│                  │  ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┐                     │
│                  │       MISSING LINK ZONE                                   │
│                  │  │   (early ELISA → modern ELISA)  │                     │
│                  │   ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─                      │
│                  │                                                           │
├──────────────────┼───────────────────────────────────────────────────────────┤
│                  │          EARLY ELISA (1980s-1990s)                        │
│                  │                                                           │
│                  │            ┌───────────┐                                  │
│                  └ ─ ─ ─ ? ─ ▶│  Barrett  │                                  │
│                               │ ELISA 1983│                                  │
│                               └─────┬─────┘                                  │
│                                     │                                        │
├─────────────────────────────────────┼────────────────────────────────────────┤
│                                     │          LEGACY (1960s-1980s)          │
│                                     │                                        │
│                                     ▼                                        │
│                               ┌───────────┐                                  │
│                               │    Vi     │                                  │
│                               │    HA     │                                  │
│                               └───────────┘                                  │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘

Legend:
  ─────▶  = Empirical bridging edge (paired data exists)
  ─ ─ ─▶  = Hypothesized/uncertain edge (may not have paired data)
  [box]   = Canonical assay node
```

### 2.2 Edge Inventory

| Edge | Source | Target | Data Type | Status |
|------|--------|--------|-----------|--------|
| E1 | Vi HA | Barrett ELISA | Paired sera | **To verify** (Barrett 1983) |
| E2 | Barrett ELISA | Modern in-house | ? | **Unknown** (gap in cascade) |
| E3 | In-house ELISA | VaccZyme | Paired sera | **Available** (Lee 2020 Figshare) |
| E4 | VaccZyme | IS-relative | Calibration | **Structural** (NIBSC 16/138) |
| E5 | In-house ELISAs | IS-relative | Commutability | **Structural** (NIBSC 16/138) |

---

## 3. Measurement Model

### 3.1 General Form

Each assay observes the latent immunity with assay-specific transformation and noise:

$$\log(\text{titer}_{\text{assay}}) = \beta_0^{(\text{assay})} + \beta_1^{(\text{assay})} \cdot \log(\text{Immunity}_{\text{latent}}) + \epsilon^{(\text{assay})}$$

Where:
- $\beta_0^{(\text{assay})}$: Assay-specific offset (units, calibration)
- $\beta_1^{(\text{assay})}$: Assay-specific scaling (ideally = 1 if linear in immunity)
- $\epsilon^{(\text{assay})} \sim N(0, \sigma^2_{\text{assay}})$: Measurement noise

### 3.2 Censoring

Most assays have detection limits:

$$\text{titer}_{\text{observed}} = \begin{cases}
\text{LOD} & \text{if titer}_{\text{true}} < \text{LOD} \text{ (left-censored)} \\
\text{titer}_{\text{true}} & \text{if LOD} \le \text{titer}_{\text{true}} \le \text{UOD} \\
\text{UOD} & \text{if titer}_{\text{true}} > \text{UOD} \text{ (right-censored)}
\end{cases}$$

### 3.3 Discreteness

Endpoint dilution titers are discrete (2-fold series): 1:20, 1:40, 1:80, ...

This introduces:
- Interval censoring (true value between adjacent dilutions)
- Heaping at common values

---

## 4. Bridging Model

### 4.1 Edge-by-Edge Conversion

For each edge E connecting assay A to assay B:

$$\log(\text{titer}_B) = \alpha_E + \beta_E \cdot \log(\text{titer}_A) + \epsilon_E$$

**Parameters:**
- $\alpha_E$: Intercept (accounts for unit differences)
- $\beta_E$: Slope (ideally ~1 if assays measure same thing proportionally)
- $\sigma_E$: Residual SD (combined measurement noise from both assays + true variability)

### 4.2 Population Effects

Bridging may depend on population:

$$\alpha_E^{(p)}, \beta_E^{(p)}, \sigma_E^{(p)} \text{ for population } p \in \{\text{acute}, \text{carrier}, \text{vaccinated}, \ldots\}$$

### 4.3 Lab/Protocol Effects

If multiple labs contribute data for an edge:

$$\alpha_E^{(\ell)}, \beta_E^{(\ell)} \text{ for lab } \ell$$

Or as random effects:
$$\alpha_E^{(\ell)} \sim N(\alpha_E, \tau_\alpha^2)$$

---

## 5. Simplifications for Working Model

*To be completed after Phase 4 data review*

| Reference Component | Simplification | Justification |
|---------------------|----------------|---------------|
| Population-specific edges | Pool if N too small | TBD |
| Slope parameter β | Fix to 1 if poorly identified | TBD |
| Early→modern ELISA link | Assume commutability or stop | TBD |
| ... | ... | ... |

---

## 6. Mapping to Available Data

*To be completed after Phase 4 joint review*

| Edge | Papers | N pairs | Population | Data Format | Estimable Parameters |
|------|--------|---------|------------|-------------|---------------------|
| HA → early ELISA | Barrett 1983 | ? | acute/carrier | ? | TBD |
| In-house → VaccZyme | Lee 2020 | ~300 | vaccinees | Individual | α, β, σ |
| ... | ... | ... | ... | ... | ... |

---

## 7. Calibration Setup

*To be completed after data mapping*

### 7.1 Likelihood Structure

| Edge | Likelihood | Notes |
|------|------------|-------|
| Continuous pairs | Censored normal | Handle LOD/LOQ |
| Binned/ordinal | Ordered logit | If only bins available |
| Concordance only | Binomial sens/spec | Minimal information |

### 7.2 Prior Specification

| Parameter | Prior | Rationale |
|-----------|-------|-----------|
| α (intercept) | N(0, 2) on log scale | Weakly informative; allows 100x offset |
| β (slope) | N(1, 0.5) | Centered on proportional, allows deviation |
| σ (residual) | Half-N(0, 1) | Weakly informative; ~1 log unit typical |
| τ (random effect SD) | Half-N(0, 0.5) | Modest lab-to-lab variation |

### 7.3 End-to-End Propagation

Given historical HA titer $t_{\text{HA}}$:

1. Sample from posterior of E1: $\log(t_{\text{early}}) \sim N(\hat\alpha_1 + \hat\beta_1 \cdot \log(t_{\text{HA}}), \hat\sigma_1^2)$
2. [If E2 exists] Sample from posterior of E2: $\log(t_{\text{modern}}) \sim N(\hat\alpha_2 + \hat\beta_2 \cdot \log(t_{\text{early}}), \hat\sigma_2^2)$
3. Sample from posterior of E3: $\log(t_{\text{VaccZyme}}) \sim N(\hat\alpha_3 + \hat\beta_3 \cdot \log(t_{\text{modern}}), \hat\sigma_3^2)$
4. (Optional) Convert to IS-relative via E4 calibration

**Output**: Posterior predictive distribution over VaccZyme-equivalent or IS-relative values, not a point estimate.
