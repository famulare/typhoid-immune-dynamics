# Dose-Response Model Data Extraction Contract

## Problem Statement

We are calibrating a modified beta-Poisson dose-response model for typhoid fever that incorporates pre-existing immunity. The model predicts probability of infection and fever given bacterial challenge dose and a correlate of protection (CoP), with immunity scaling the effective "hit rate" in a mechanistically interpretable way.

### The Model

**Infection given dose:**
$$
P(\text{infection} | D, \text{CoP}) = 1 - \left(1 + D \cdot \frac{2^{1/\alpha_I} - 1}{N_{50,I}}\right)^{-\alpha_I / \text{CoP}^{\gamma_I}}
$$

**Fever given dose:**
$$
P(\text{fever} | D, \text{CoP}) = 1 - \left(1 + D \cdot \frac{2^{1/\alpha_F} - 1}{N_{50,F}}\right)^{-\alpha_F / \text{CoP}^{\gamma_F}}
$$

**Parameters to estimate:**
- $N_{50}$: Dose for 50% probability in naive individuals
- $\alpha$: Shape parameter controlling dose-response steepness
- $\gamma$: Immunity scaling exponent (how strongly CoP reduces risk)

### The Challenge

The historical and contemporary literature on typhoid challenge studies is highly heterogeneous:
- **Non-standardized reporting**: Dose units, outcome definitions, and subject characteristics vary widely
- **Scattered information**: Key data spread across tables, text, and supplementary materials
- **Overlapping cohorts**: Some subjects appear in multiple publications
- **Missing immunity data**: Pre-challenge serology often unavailable, requiring latent variable inference
- **Outcome definition drift**: "Typhoid fever" meant different things across decades

### Goals

1. **Systematic extraction**: Create structured, traceable records from all relevant papers
2. **Preserve uncertainty**: Store ambiguous values as intervals or alternative interpretations
3. **Enable calibration**: Produce analysis-ready data that supports Bayesian inference with:
   - Binomial likelihoods for dose-outcome observations
   - Study-level random effects for protocol heterogeneity
   - Latent immunity inference where serology is missing
4. **Document decisions**: Maintain clear audit trail of interpretive choices

---

## Integrated Plan: Typhoid Dose-Response Model Data Extraction & Calibration

### Decision Bookkeeping Convention

Throughout extraction, all interpretive decisions are tagged:
- **[USER-LOCKED]**: Decisions explicitly approved during joint review (binding)
- **[ASSISTANT-PROPOSED]**: Default interpretation pending review
- **[OPEN]**: Requires judgment call or later sensitivity analysis

---

## Phase 1: Setup & Infrastructure

**1.1 Create folder structure:**
```
metastudies/dose_response_model/
├── input_papers/              # (exists) PDFs
├── extracts/                  # markdown extracts per paper
├── schemas/                   # YAML schemas (created in Phase 4, not upfront)
├── analysis_data/             # final CSV(s) for calibration
├── notes/
│   ├── cross_cutting_observations.md
│   ├── paper_triage.md
│   ├── outcome_mapping.md     # how observed outcomes → model variables
│   └── identifiability_memo.md # which params constrained by which data
└── calibration/
    ├── likelihood_design.md   # statistical model choices + rationale
    └── priors.yaml            # prior specifications
```

**1.2 Create extraction template (markdown):**

Each extract will use decision bookkeeping throughout with [USER-LOCKED], [ASSISTANT-PROPOSED], and [OPEN] tags.

**1.3 Create outcome mapping draft** (`notes/outcome_mapping.md`):

Initial mapping (to refine as we extract):
```
INFECTION proxies:
- blood_culture_positive → infection (high specificity, imperfect sensitivity)
- PCR_positive → infection (higher sensitivity than culture)
- stool_shedding → infection (depends on sampling schedule)
- seroconversion → infection (delayed, may miss some)

FEVER proxies:
- oral_temp_≥100F → fever (common threshold)
- oral_temp_≥100.4F → fever (CDC definition)
- "clinical typhoid" → fever (composite, definition varies)
- "typhoid fever diagnosis" → fever (often composite)
```

---

## Phase 2: Batch First-Pass Extraction (Claude solo)

**2.1 For each paper, create `extracts/{author}_{year}.md`:**

#### Study Metadata
- Citation, DOI if available
- Year(s) conducted (not just published)
- Location, institution
- Primary study question (challenge study, vaccine trial, re-analysis?)

#### Challenge Protocol
- **Strain**: Quailes, Ty2, other (with notes on provenance)
- **Dose(s)**: Exact CFU or range; if ambiguous, store as `[low, high]` interval
- **Dose measurement basis**: Plated CFU / estimated / stated target / unknown
- **Delivery medium**: Milk, bicarbonate solution, water, capsule, other
- **Buffering**: Yes/no/unknown; NaHCO₃ concentration if stated
- **Fasting window**: Hours pre-challenge if stated
- **Inoculum volume**: mL if stated
- **Administration timing/procedure**: Any relevant details

#### Subject Characteristics (per cohort/arm)
- **N**: Total and by dose arm
- **Age**: Mean/median/range
- **Sex**: Distribution if reported
- **Geography/origin**: Where recruited, nationality
- **Year**: When study conducted
- **SES proxies**: Occupation, student status, institutionalization
- **Inclusion/exclusion criteria**: Especially prior typhoid, vaccination status

#### Baseline Immunity (CRITICAL - capture everything)
- **Prior typhoid history**: Yes/no/unknown; how assessed (self-report, records, serology)
- **Prior vaccination**: Type, timing
- **Baseline serology** (for each assay reported):
  - Assay type (Vi IgG, H agglutinin, O agglutinin, bactericidal, other)
  - Units
  - Timing relative to challenge
  - Summary statistics (mean, median, range, distribution)
- **[OPEN]** flags for missing immunity data that will need latent inference

#### Outcomes (with ascertainment details)
For each outcome reported:
- **Raw definition**: Quote exactly how paper defines it
- **Ascertainment schedule**: How often checked, for how long
- **Diagnostic method**: Blood culture, PCR, stool culture, temperature measurement
- **Time window**: e.g., "within 14 days of challenge"
- **Numerator/denominator by dose arm**: Preserve exact reporting
- **Sensitivity context**: Any discussion of assay sensitivity/specificity

**Key rule**: If counts are ambiguous, store **alternative interpretations** explicitly, don't pick one silently.

#### Post-Challenge Measures (for later use)
- Serology changes
- Symptom timing
- Treatment details

#### Data Quality & Uncertainty Notes
- What's explicit vs inferred?
- Page/table references for every key number
- Ambiguities flagged with [OPEN] or [ASSISTANT-PROPOSED]

#### Cross-References
- Does this paper cite/reanalyze data from other papers in our set?
- Suspected cohort overlaps with other papers

---

**2.2 Create cross-cutting observations** (`notes/cross_cutting_observations.md`):

- **Cohort overlaps**: Which papers share subjects?
- **Strain lineage**: Quailes history, Ty2, comparability
- **Methodological evolution**: Maryland 1960s-70s → Oxford 2010s
- **Outcome definition drift**: Did definitions change over time?
- **Immunity measurement evolution**: What assays when?
- **Dose reporting conventions**: Units, calibration issues
- **Delivery medium effects**: Evidence for milk vs bicarbonate differences?
- **Recurring data gaps**: What's consistently missing?

---

**2.3 Create triage document** (`notes/paper_triage.md`):

| Paper | Status | Dose Info | Multi-Dose? | Pre-Challenge Immunity? | Key Outcomes | Fit Role | Notes |
|-------|--------|-----------|-------------|------------------------|--------------|----------|-------|
| ... | include/exclude/maybe | y/n/partial | y/n | y/n/partial | list | core/support/exclude | ... |

**Fit Role:**
- **Core**: Best-reported dose arms + clear outcomes + some immunity context → primary calibration
- **Support**: Informs heterogeneity priors, sensitivity analyses, or validation
- **Exclude**: No usable dose-outcome data after extraction

---

**2.4 Definition of Done for Phase 2:**

A paper extraction is complete when:
- [ ] Status determined (include/exclude)
- [ ] If included: markdown extract exists with all sections
- [ ] Page/table references for every key number
- [ ] All ambiguities flagged as [ASSISTANT-PROPOSED] or [OPEN]
- [ ] Cross-references to other papers noted
- [ ] Triage table updated

---

## Phase 3: Paper-by-Paper Joint Review

**For each paper (you reading alongside):**

1. I present my extract summary
2. You verify against your reading
3. We resolve [ASSISTANT-PROPOSED] items → convert to [USER-LOCKED]
4. We discuss [OPEN] items → either lock a decision or flag for sensitivity analysis
5. Finalize the markdown extract
6. Update outcome mapping if new outcome types appear
7. Update cross-cutting observations if new patterns emerge

**Review order:**
1. **Core multi-dose papers first** (Hornick series likely) - establish dose-response shape
2. **Support papers with immunity data** - immunity dependence
3. **Remaining papers** - quick pass for anything useful
4. **Exclusions** - confirm and document reasons

---

## Phase 4: Normalization & Schema Design

**4.1 Design YAML schema** (AFTER seeing all data)

Based on what we actually extracted, design schema covering:
- Study-level metadata
- Cohort/arm characteristics
- Dose-outcome observations
- Immunity measures
- Decision provenance ([USER-LOCKED] vs sensitivity analysis needed)

**4.2 Convert extracts to YAML**

Each finalized markdown → structured YAML

**4.3 Compile analysis-ready CSV**

Flatten to modeling-ready format:
- One row per dose-arm-outcome combination
- Include outcome_type, outcome_proxy_type, definition_notes
- Include immunity summary measures where available
- Flag rows needing latent immunity inference

---

## Phase 5: Calibration Problem Design

**5.1 Outcome mapping finalization** (`notes/outcome_mapping.md`)

Lock decisions on how each observed outcome maps to model variables.

**5.2 Likelihood structure** (`calibration/likelihood_design.md`)

- Binomial: `y ~ Binomial(n, p(dose, CoP; θ))`
- Separate likelihoods for infection vs fever
- How to handle:
  - Dose intervals → integrate over range or sample latent true dose
  - Multiple outcome interpretations → mixture likelihood or sensitivity analysis
  - Outcome proxy differences → measurement model or stratification

**5.3 Heterogeneity structure**

- Study random effects on which parameters? (N₅₀? multiplicative dose factor?)
- Outcome-definition random effects if using mixed proxies

**5.4 Latent immunity model**

- For cohorts without baseline serology: latent CoP with hierarchical prior
- Prior informed by: geography, era, recruitment criteria, any partial info
- If baseline serology exists: measurement model linking titer → CoP

**5.5 Identifiability memo** (`notes/identifiability_memo.md`)

Pre-fitting diagnostic:
- Can we identify N₅₀ and α from available dose spread?
- Is γ identifiable without baseline immunity variation?
- Which parameters are constrained by which studies?
- What requires strong priors vs data-driven?

---

## Phase 6: Prior Specification (`calibration/priors.yaml`)

- Core parameters: Weakly informative centered at current defaults
- Study random effects: Half-normal on SDs
- Latent CoP: Informative based on cohort characteristics
- Document rationale for each prior choice

---

## Phase 7: Fit, Validate, Document

**7.1 Implementation** (Stan/brms)

**7.2 Diagnostics**
- Convergence (R-hat, ESS, trace plots)
- Posterior predictive checks

**7.3 Sensitivity analyses**
- Prior sensitivity
- Include/exclude specific studies
- Alternative outcome mappings
- With/without latent immunity

**7.4 Documentation**
- Update model description with fitted values
- Clear data provenance
- Archive all materials

---

## Working Agreement

This document serves as the contract for the dose-response model calibration work. Progress will be tracked in the accompanying `progress_checklist.md` file. The work will proceed across multiple sessions as needed, with clear handoff points documented in the checklist.
