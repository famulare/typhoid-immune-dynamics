# Typhoid Dose-Response Model Calibration

This folder contains materials for calibrating a modified beta-Poisson dose-response model for typhoid fever from human challenge study literature.

## Current result

The calibration workflow has landed and works end to end: fitting, diagnostics,
summary generation, and the figure suite are all committed. The preferred
working model is **Tier 2 individualized vaccine** (`t2-indiv-vax`, stage
`phi-rho-eta-psi-vax`), which includes the individualized Darton placebo and
vaccine arms, Oxford shedding correction, and the infection-definition
correction.

The [Tier 2 posterior summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md)
is the entry point for the fitted results, parameter tables, diagnostics, and
the full set of linked figures. The [Tier 2 prior-predictive summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax-prior/summary.md)
contains the corresponding prior run. The posterior fit has one divergent
transition in 4,000 draws; see the summaries for the complete diagnostic and
prior-sensitivity record.

For a fresh clone, the discoverable rebuild ladder is exactly
`t1-grouped`, `t1-indiv`, `t1-indiv-vax`, `t2-indiv`, and `t2-indiv-vax`.
The retired `t2-grouped` entry remains in the registry for historical provenance
but is omitted from the CLI list and generated ladder.

## Purpose

Estimate how the probability of infection and fever depends on:
- Bacterial challenge dose
- Pre-existing immunity (Correlate of Protection)

## Key Documents

| Document | Purpose |
|----------|---------|
| `dose_response_extraction_contract.md` | Project contract: goals, workflow, decision conventions |
| `progress_checklist.md` | Phase-by-phase progress tracking |
| `dose_response_model_specification.md` | Model specification: biology, causal structure, equations |
| `notes/outcome_mapping.md` | Decision rules for mapping observed outcomes to model variables |
| `notes/paper_triage.md` | Paper-by-paper inclusion/exclusion decisions |
| `notes/cross_cutting_observations.md` | Patterns across the literature corpus |
| `calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md` | Preferred Tier 2 fit, diagnostics, parameter results, and figure index |

## Folder Structure

```
dose_response_model/
├── dose_response_model_specification.md   # Model specification
├── input_papers/       # Source PDFs (21 papers)
├── extracts/           # Markdown extracts per paper (Phase 2 output)
├── notes/              # Working notes
├── schemas/            # YAML schemas (Phase 5)
├── analysis_data/      # Final CSV for calibration (Phase 5)
└── calibration/        # Likelihood design, priors (Phases 6-7)
```

## Workflow Phases

1. **Setup** - Folder structure, templates
2. **Extraction** - First-pass data extraction from each paper
3. **Model Specification** - Define reference model (what we'd fit with perfect data)
4. **Joint Review** - Paper-by-paper validation with user
5. **Normalization** - Collapse to working model, create analysis-ready data
6. **Calibration Design** - Likelihood structure, heterogeneity model
7. **Priors** - Prior specification with rationale
8. **Fit & Validate** - Implementation, diagnostics, sensitivity analyses

Current status tracked in `progress_checklist.md`.

## Decision Conventions

Throughout extraction, interpretive decisions are tagged:
- **[USER-LOCKED]**: Approved during joint review (binding)
- **[ASSISTANT-PROPOSED]**: Default pending review
- **[OPEN]**: Requires judgment or sensitivity analysis
