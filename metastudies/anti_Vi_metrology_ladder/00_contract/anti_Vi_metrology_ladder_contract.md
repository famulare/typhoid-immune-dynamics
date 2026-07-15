# Anti-Vi Metrology Bridge Data Extraction Contract

## Problem Statement

We are constructing an empirical bridge from **legacy Vi hemagglutination/agglutination (HA)** measurements (1960s–80s) to **modern anti-Vi IgG ELISA** (VaccZyme-compatible or IS-relative scale). This enables integration of historical typhoid immunity data with contemporary vaccine trial results for dose-response modeling.

### The Bridge

We are fitting assay-to-assay conversion models that propagate:
- Legacy HA titer → Early ELISA titer → Modern ELISA (VaccZyme/IS-relative)

**Parameters to estimate per edge:**
- Location (intercept on log scale)
- Scale (slope on log scale)
- Residual variance (measurement + biological)
- Population-specific effects (if estimable)

*See `notes/metrology_bridge_specification.md` for full mathematical formulation and graph structure.*

### The Challenge

The historical and contemporary anti-Vi serology literature is highly heterogeneous:
- **Assay evolution**: HA → early ELISA → modern standardized ELISA over 50+ years
- **Unit drift**: "Titer" definitions changed (endpoint dilution → concentration → IS-relative)
- **Scattered bridging data**: Few papers directly compare assays on the same sera
- **Population dependence**: Bridging relationships may differ for carriers vs vaccinees vs acute cases
- **Antigen/coating differences**: Native Vi vs modified formats affect commutability

### Goals

1. **Systematic extraction**: Create structured records from all papers with paired assay data
2. **Preserve uncertainty**: Store censored values, alternative interpretations, assay-specific quirks
3. **Enable bridging**: Produce analysis-ready paired datasets that support Bayesian inference with:
   - Log-normal measurement models with censoring
   - Lab/protocol random effects
   - Population stratification where data allow
4. **Document decisions**: Maintain clear audit trail of bridging assumptions and stop conditions

---

## Integrated Plan: Anti-Vi Metrology Bridge Extraction & Calibration

### Decision Bookkeeping Convention

Throughout extraction, all interpretive decisions are tagged:
- **[USER-LOCKED]**: Decisions explicitly approved during joint review (binding)
- **[ASSISTANT-PROPOSED]**: Default interpretation pending review
- **[OPEN]**: Requires judgment call or later sensitivity analysis

---

## Phase 1: Setup & Infrastructure

**1.1 Create folder structure:**
```
metastudies/anti_Vi_metrology_ladder/
├── input_studies/              # (exists) PDFs and source files
├── extracts/                   # markdown extracts per paper
├── schemas/                    # YAML schemas (created in Phase 5, not upfront)
├── analysis_data/              # final paired datasets for bridging
├── notes/
│   ├── metrology_bridge_specification.md  # reference bridging model
│   ├── cross_cutting_observations.md
│   ├── paper_triage.md
│   ├── assay_mapping.md        # how reported assays → canonical assay nodes
│   └── identifiability_memo.md # which edges constrained by which data
└── calibration/
    ├── likelihood_design.md    # statistical model choices + rationale
    └── priors.yaml             # prior specifications
```

**1.2 Create extraction template (markdown):**

Each extract will use decision bookkeeping throughout with [USER-LOCKED], [ASSISTANT-PROPOSED], and [OPEN] tags.

**1.3 Create assay mapping draft** (`notes/assay_mapping.md`):

Initial mapping (to refine as we extract):
```
LEGACY ASSAYS (pre-1990):
- Vi_hemagglutination → HA (endpoint dilution, 2-fold series)
- Vi_agglutination → HA (may be same or different technique)
- Vi_passive_hemagglutination → HA (similar to above)

EARLY ELISA (1980s-1990s):
- Barrett_ELISA_1983 → early_ELISA (endpoint dilution titer)
- in_house_ELISA_unspecified → early_ELISA (requires protocol details)

MODERN ELISA (2000s+):
- VaccZyme_anti_Vi_IgG → VaccZyme (EU/mL, kit-calibrated)
- in_house_ELISA_IS_calibrated → IS_ELISA (IU/mL or relative potency)
- in_house_ELISA_Vi_IgGR1_calibrated → IS_ELISA (linked to US reference)

REFERENCE STANDARDS:
- NIBSC_16_138 → IS (1st International Standard)
- Vi_IgGR1_2011 → US_reference (FDA reference reagent)
```

---

## Phase 2: Batch First-Pass Extraction (Claude solo)

**2.1 For each paper, create `extracts/{author}_{year}.md`:**

#### Study Metadata
- Citation, DOI if available
- Year(s) conducted (not just published)
- Location, institution
- Primary study purpose (assay development, outbreak investigation, vaccine trial, standards work?)

#### Assays Used
For each assay reported:
- **Assay type**: HA, ELISA, other
- **Antigen**: Native Vi, modified Vi, coating chemistry
- **Detection**: IgG, total Ig, agglutination endpoint
- **Units**: Titer (dilution), concentration (EU/mL, IU/mL), OD, other
- **Calibration**: Reference standard used, lot if stated
- **LOD/LOQ**: Detection and quantitation limits
- **Protocol details**: Dilution series, threshold definitions, incubation conditions

#### Paired Measurements (CRITICAL)
For each assay-assay comparison:
- **Assay A**: Type, units, values reported
- **Assay B**: Type, units, values reported
- **Sample pairing**: Same sera? Same individuals different timepoints?
- **N**: Number of paired observations
- **Population**: Acute typhoid, carriers, vaccinees, endemic controls, standards panel
- **Format**: Individual values, grouped summary, correlation only, positivity concordance only
- **Page/table reference**: Exact location in paper

#### Sample Characteristics
- **Population types**: Acute, carrier, vaccinated, endemic baseline, non-endemic control
- **Geography**: Endemic vs non-endemic regions
- **N per population**: Sample sizes
- **Selection criteria**: How samples were chosen (consecutive, convenience, designed panel)

#### Assay Comparison Results
- **Correlation**: r, r², ρ if reported
- **Regression**: Equation if fitted
- **Concordance**: Sensitivity/specificity at thresholds
- **Bias/systematic differences**: Any noted discrepancies

#### Data Quality & Uncertainty Notes
- What's explicit vs inferred?
- Page/table references for every key number
- Ambiguities flagged with [OPEN] or [ASSISTANT-PROPOSED]
- Digitization required? (figure only, no table)

#### Cross-References
- Does this paper cite/reanalyze data from other papers in our set?
- Is this assay protocol used in other papers?
- Link to NIBSC/WHO standards work?

---

**2.2 Create cross-cutting observations** (`notes/cross_cutting_observations.md`):

- **Assay lineage**: Which ELISAs descended from Barrett 1983? Which use NIBSC standards?
- **Commutability evidence**: What does the IS documentation say about assay equivalence?
- **Population effects**: Do bridging relationships differ by population type?
- **Era patterns**: How did titer definitions change over time?
- **Antigen/coating effects**: Native Vi vs protein-coated vs conjugated
- **Recurring data gaps**: What's consistently missing?

---

**2.3 Create triage document** (`notes/paper_triage.md`):

| Paper | Status | Paired Data? | Assays Compared | Population | N pairs | Data Format | Bridge Role | Notes |
|-------|--------|--------------|-----------------|------------|---------|-------------|-------------|-------|
| Barrett 1983a | include | yes | HA ↔ early ELISA | acute/carrier | ? | table/figure | core | ... |
| Lee 2020 | include | yes | in-house ↔ VaccZyme | vaccinees | ~300 | Figshare | core | ... |
| NIBSC 16/138 | include | commutability | multiple ELISAs | standards panel | ? | IFU | structural | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |

**Bridge Role:**
- **Core**: Direct paired measurements on shared samples → primary calibration
- **Structural**: Provides commutability/equivalence constraints but not direct pairs
- **Support**: Informs heterogeneity priors, sensitivity analyses, or validation
- **Exclude**: No usable bridging data after extraction

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

## Phase 3: Reference Bridge Specification

**Purpose**: Define the most elaborate bridging model the literature supports—what we would fit with optimal data. This serves as a benchmark for reasoning about which simplifications are forced by data limitations versus chosen for parsimony.

**3.1 Write reference bridge document** (`notes/metrology_bridge_specification.md`):

- Assay graph: Nodes = canonical assay types, Edges = bridging relationships
- Latent immunity construct: The "true signal" all assays attempt to measure
- Measurement model: Each assay as noisy observation of latent + assay-specific effects
- Bridging model: Edge-by-edge conversion functions
- Population effects: How bridging relationships vary by serum provenance
- Commutability constraints: What IS documentation tells us about ELISA equivalence
- Observational model: $P(\text{observed titer} | \text{latent immunity}, \text{assay}, \text{lab})$

**3.2 Create bridge graph** showing assay nodes and empirical edges

**3.3 Document simplifications** needed to reach practical model (to be filled after Phase 4)

**3.4 Map reference model to available data** (to be filled after Phase 4)

---

## Phase 4: Paper-by-Paper Joint Review

**For each paper (you reading alongside):**

1. I present my extract summary
2. You verify against your reading
3. We resolve [ASSISTANT-PROPOSED] items → convert to [USER-LOCKED]
4. We discuss [OPEN] items → either lock a decision or flag for sensitivity analysis
5. Finalize the markdown extract
6. Update assay mapping if new assay types appear
7. Update cross-cutting observations if new patterns emerge

**Review order:**
1. **Core paired-data papers first** (Barrett 1983, Lee 2020) - establish bridging edges
2. **Structural papers** (NIBSC/WHO) - establish commutability constraints
3. **Support papers** - additional population contexts
4. **Exclusions** - confirm and document reasons

---

## Phase 5: Normalization & Schema Design

**5.1 Collapse reference model to working model** (`notes/working_model.md`)

Based on Phase 4 data review, document:
- Which reference model edges have data support
- Specific simplifications chosen and justification
- Working model specification (edges we'll actually calibrate)
- Stop conditions triggered (ordinal-only bridges, population-specific limits)

**5.2 Design YAML schema** (informed by working model)

Schema covering:
- Paper-level metadata
- Assay protocol descriptions
- Paired observations (with censoring flags)
- Population characteristics
- Decision provenance ([USER-LOCKED] vs sensitivity analysis needed)

**5.3 Convert extracts to YAML**

Each finalized markdown → structured YAML

**5.4 Compile analysis-ready paired datasets**

Format per bridge edge:
- One row per paired observation (or summary if individual data unavailable)
- Include assay_A_type, assay_A_value, assay_A_units, assay_A_censoring
- Include assay_B_type, assay_B_value, assay_B_units, assay_B_censoring
- Include population, lab, timepoint (if applicable)
- Flag rows needing special handling

---

## Phase 6: Bridge Model Design

**6.1 Assay mapping finalization** (`notes/assay_mapping.md`)

Lock decisions on canonical assay node assignments.

**6.2 Likelihood structure** (`calibration/likelihood_design.md`)

For each edge:
- **If continuous paired data**: Log-linear regression with censored normal errors
  ```
  log(titer_B) ~ Normal(β₀ + β₁ * log(titer_A), σ²)
  ```
- **If ordinal/binned data**: Ordered logit mapping
  ```
  P(bin_B = j | titer_A) = ordinal_regression(titer_A, θ)
  ```
- **If concordance only**: Sensitivity/specificity model
  ```
  P(positive_B | titer_A) = logistic(β₀ + β₁ * titer_A)
  ```

How to handle:
- Left-censoring at LOD
- Right-censoring at saturation
- Heteroscedasticity (variance increasing with titer)

**6.3 Heterogeneity structure**

- Lab random effects on intercept and/or slope
- Population random effects (carriers vs vaccinees vs acute)
- Era effects if spanning decades

**6.4 End-to-end propagation**

- Cascade through multiple edges with uncertainty propagation
- Monte Carlo sampling from posterior predictive at each edge
- Final output: distribution over VaccZyme-equivalent or IS-relative values

**6.5 Identifiability memo** (`notes/identifiability_memo.md`)

Pre-fitting diagnostic:
- Which edges have enough data for slope estimation vs intercept-only?
- Can population effects be identified or must they be pooled?
- What requires strong priors vs data-driven?

---

## Phase 7: Prior Specification (`calibration/priors.yaml`)

- Edge parameters: Weakly informative (log-linear slopes near 1, intercepts near 0)
- Random effect SDs: Half-normal
- Residual variance: Half-normal, informed by known assay CVs
- Document rationale for each prior choice

---

## Phase 8: Fit, Validate, Document

**8.1 Implementation** (Stan/brms or R)

**8.2 Diagnostics**
- Convergence (R-hat, ESS, trace plots)
- Posterior predictive checks per edge
- Leave-one-out cross-validation

**8.3 Sensitivity analyses**
- Prior sensitivity
- Include/exclude specific studies
- With/without population stratification
- Ordinal vs continuous treatment of binned data

**8.4 Stop condition evaluation**
- If posterior predictive checks fail, document which bridges are not supportable
- Downgrade to ordinal bins or "latent immunity signal" interpretation as needed

**8.5 Documentation**
- Update bridge specification with fitted values
- Clear data provenance
- Archive all materials

---

## Working Agreement

This document serves as the contract for the anti-Vi metrology bridge calibration work. Progress will be tracked in the accompanying `progress_checklist.md` file. The work will proceed across multiple sessions as needed, with clear handoff points documented in the checklist.

---

## Evidence Anchors (from original contract)

### A1) Legacy method → early ELISA bridge (same era, same antigen target)

1. **Barrett et al., J Clin Microbiol (1983)**: Developed anti-Vi ELISA and compared results to a previously described hemagglutination assay (HA) on human sera. ([ASM Journals][1])
2. **Barrett et al., J Clin Microbiol (1983)** (related outbreak/carrier context): Vi antibody by ELISA and hemagglutination techniques in a small investigation. ([ASM Journals][2])

*Role in cascade:* establishes at least one empirically documented relationship between HA and an early ELISA readout on overlapping sera.

### A2) ELISA ↔ VaccZyme bridge (modern era)

3. **Lee et al., PLOS NTD (2020)**: Measures a subset of sera in both an in-house ELISA and VaccZyme, fits a statistical model, and converts in-house titers to VaccZyme-scale values. ([PLOS][3])

*Role in cascade:* explicit, data-driven mapping from a non-VaccZyme ELISA to VaccZyme.

### A3) Standardization / commutability bridge across ELISAs (VaccZyme included)

4. **NIBSC 16/138 IFU**: States that indirect in-house ELISAs of native Vi are commutable with the VaccZyme ELISA under the collaborative study; describes comparison against U.S. reference reagent Vi-IgGR1, 2011, multiple ELISA formats including VaccZyme. ([nibsc.org][4])
5. **WHO collaborative study report / ECBS document** (candidate IS evaluation): Again states candidate IS 16/138 suitability and commutability in VaccZyme and in-house ELISAs. ([WHO][5])
6. **Rijpkema et al., Vaccine (2018)**: Published version of the IS work, including VaccZyme as the "commercial ELISA." ([ScienceDirect][6])

*Role in cascade:* gives a "hub" that supports ELISA comparability and provides a reference scale (IS-relative).

---

## Caveats (first-class concerns)

### D1) HA and ELISA may not be monotone-equivalent across populations

Reasons:
- HA depends on isotype mix, avidity, and assay conditions; ELISA depends on antigen presentation, conjugate, and cutoff definition.

Operationalization:
- Fit separate mappings by serum provenance class (acute infection vs vaccinated vs carriers), if data allow.
- If not, widen uncertainty and prohibit fine-grained conversion.

### D2) "ELISA titers" are not concentrations

Early ELISAs often report endpoint dilution titers, not standardized IgG concentration.

Operationalization:
- Keep "titer" as its own measurement type with censoring and dilution-step discreteness.
- Avoid pretending titers are linear in IgG.

### D3) Antigen differences: native Vi vs modified/coated formats

NIBSC notes commutability for certain ELISA formats (e.g., protein pre-/co-coating). ([nibsc.org][4])

Operationalization:
- Annotate each ELISA by coating chemistry; treat non-commutable formats as separate nodes.

### D4) Vendor kit drift / lot effects

VaccZyme has evolved; calibrator linkage matters.

Operationalization:
- Record kit version/calibrator strategy for every dataset; if unknown, treat as additional multiplicative noise.

### D5) Unavoidable missing link risk: early ELISA ↔ modern ELISA may be absent

Operationalization:
- Predefine stop criteria: if no paired datasets or credible intermediate bridges exist, you will only support:
  - HA → "latent Vi seropositivity / coarse bins"
  - modern ELISA ↔ VaccZyme mapping separately
  - but not a continuous HA→VaccZyme conversion

---

## Stop Criteria (to avoid forcing coherence)

You stop claiming a continuous cascade if any of the following occur:

1. **E1 (HA→ELISA)** provides only positivity concordance and no graded mapping, and no other paper provides graded paired data. ([Europe PMC][7])
2. No paper provides any plausible intermediate that links "endpoint dilution ELISA titer" style outputs to "IgG concentration/potency" style outputs.
3. Mappings differ qualitatively (non-overlapping) across populations and you lack metadata to stratify.
4. Predictive checks show that converted values systematically mis-predict known vaccine response distributions (e.g., implausible baselines).

In those cases, the output becomes:
- a qualitative/ordinal prior (e.g., "low/medium/high Vi immunity")
- or a latent factor model where historical HA contributes weakly rather than deterministically.

---

[1]: https://journals.asm.org/doi/10.1128/jcm.17.4.625-627.1983 "Enzyme-linked immunosorbent assay for detection of human antibodies to ..."
[2]: https://journals.asm.org/doi/10.1128/jcm.18.6.1320-1322.1983 "Identification of a carrier by using Vi enzyme-linked immunosorbent ..."
[3]: https://journals.plos.org/plosntds/article?id=10.1371%2Fjournal.pntd.0008171 "Comparison of anti-Vi IgG responses between two clinical studies ... - PLOS"
[4]: https://nibsc.org/documents/ifu/16-138.pdf "1st IS for Anti-Typhoid capsular Vi polysaccharide IgG (Human) - NIBSC"
[5]: https://cdn.who.int/media/docs/default-source/biologicals/bs-documents-%28ecbs%29/2017-documents/bs.2017.2307_who_is_anti-typhoid_capsular_vi_polysaccharide_igg_%28human%29.pdf "WHO IS Anti-Typhoid Capsular Vi Polysaccharide IgG (Human)"
[6]: https://www.sciencedirect.com/science/article/pii/S1045105618302835 "Establishment of the first International Standard for human anti ..."
[7]: https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC272705&blobtype=pdf "Enzyme-Linked Human Antibodies to Salmonella typhi Vi Antigen"
