# Typhoid Dose-Response Model Calibration

This folder contains materials for calibrating a modified beta-Poisson dose-response model for typhoid fever from human challenge study literature.

## Current result

**Project status (locked 2026-08-01):** the reproducible calibration deliverable
is complete. The canonical working model is **Tier 2 individualized vaccine**
(`t2-indiv-vax`, stage `phi-rho-eta-psi-vax`). It includes the individualized
Darton placebo and vaccine arms, the Oxford shedding correction, the infection-
definition correction, and the shared beta-binomial overdispersion term.

This is a lock on the default configuration and result artifact for downstream
work, not a claim that the model is fully identified or that all adequacy questions
are resolved.

Two entry points, depending on what you want:

- [**Tier 2 findings report**](tier2_findings_report.md) — what the calibration
  says and what it means. Dose-response for infection and fever, the anti-Vi
  titer effect, fever severity and detection, treatment truncation and stool
  shedding, Ty21a and M01ZH09 protection beyond anti-Vi, and where the corpus
  disagrees with itself. Interpretation, with the open questions named.
- [**Tier 2 posterior summary**](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md)
  — the generated output: every fitted parameter, sampler diagnostics, priorsense
  power-scaling, and the full figure index. Authoritative on numbers.

The [Tier 2 prior-predictive summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax-prior/summary.md)
contains the corresponding prior run. The posterior has **1 divergent transition
in 4,000 draws (0.025%)**, no max-treedepth hits, minimum E-BFMI 0.947, and
maximum reported R-hat 1.003; retain that residual as a review item rather than
calling the fit diagnostically perfect.

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
| `tier2_findings_report.md` | **What the calibration says.** Scientific findings and their decision consequences, written for a reader who will not open the fit |
| `dose_response_extraction_contract.md` | Project contract: goals, workflow, decision conventions |
| `progress_checklist.md` | Phase-by-phase progress tracking |
| `onboarding_one_pager.md` | Current project orientation and Tier 2 result lock |
| `dose_response_model_specification.md` | Model specification: biology, causal structure, equations |
| `joint_inference_plan.md` | Locked joint likelihood and identifiability assumptions |
| `stan_model_structure.md` | Current Stan program dataflow and likelihood dispatch |
| `notes/outcome_mapping.md` | Decision rules for mapping observed outcomes to model variables |
| `notes/paper_triage.md` | Paper-by-paper inclusion/exclusion decisions |
| `notes/cross_cutting_observations.md` | Patterns across the literature corpus |
| `calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md` | Generated: every parameter, sampler diagnostics, priorsense, full figure index. Authoritative on numbers |

## Folder Structure

```
dose_response_model/
├── dose_response_model_specification.md   # Model specification
├── input_papers/       # Source PDFs (21 papers)
├── extracts/           # Markdown extracts per paper (Phase 2 output)
├── notes/              # Working notes
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
