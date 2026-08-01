# Dose-response metastudy: current onboarding

**Status: Tier 2 working result locked 2026-08-01.** The project is complete as a
reproducible calibration deliverable. The lock fixes the default model and result
artifact for downstream work; it does not claim that every parameter is strongly
identified or that no sensitivity analysis is warranted.

## What this project does

This project calibrates a mechanistic typhoid dose-response model to human
challenge studies. It estimates how challenge dose and pre-existing immunity shape
two linked outcomes:

1. **Infection**: a broad latent infection process, observed through study-specific
   culture or shedding definitions.
2. **Fever**: a clinical outcome conditional on infection, observed with different
   temperature and diagnostic thresholds.

The two historical data regimes provide complementary information:

| Regime | Contribution | Main limitation |
|---|---|---|
| Maryland, 1960s-70s, milk delivery | Broad dose ladder, about 10^3-10^9 CFU | No baseline anti-Vi measurements; heterogeneous definitions |
| Oxford, 2010s, bicarbonate delivery | Modern protocols, anti-Vi titres, vaccine contrasts | Narrow dose range; treatment truncates later shedding |

The scalar `delta` bridges milk to bicarbonate-equivalent dose. It is useful but
strongly confounded with the infection and fever dose scales, so it is not a direct
measurement of gastric survival.

## The current model

The core dose-response is a modified beta-Poisson curve with immunity scaling:

```text
P(outcome | dose, CoP) = 1 - (1 + dose * (2^(1/alpha) - 1) / N50)^(-alpha / CoP^gamma)
```

The implemented model retains the cascade `P(infection) * P(fever | infection)`
and adds the measurement processes needed to combine the studies:

- `phi(T,D)`: fever-definition sensitivity;
- `psi_d`: infection-definition sensitivity, anchored by the nested Darton
  stool-vs-`bact_or_stool` intersection (15 of 26);
- `eta(D)`: Oxford shedding detection after treatment truncation;
- a two-component latent Maryland immunity mixture;
- `grand_overdispersion_rho`: one shared beta-binomial overdispersion parameter;
- individualized Darton, M01ZH09, and Ty21a vaccine arms.

CoP is denominated on the VaccZyme anti-Vi IgG scale for Oxford data, relative to a
naive reference. Maryland CoP is a latent anti-Vi-equivalent proxy, not a measured
Maryland serologic quantity.

## Locked Tier 2 result

The canonical result is **`t2-indiv-vax`**, stage
**`phi-rho-eta-psi-vax`**:

- 182 Stan observations: 29 grouped and 153 individual;
- posterior run: 4 chains x 1,000 warmup + 1,000 sampling draws;
- 1 divergent transition / 4,000 draws (0.025%);
- no max-treedepth hits, minimum E-BFMI 0.947, maximum reported R-hat 1.003;
- prior-predictive run and posterior-predictive figures retained alongside it.

Read the [posterior summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md)
first, then the [progress checklist](progress_checklist.md) for completion status,
residual risks, and provenance. The fit is a locked working baseline, not a
permission to present weakly identified quantities as precise biological facts.

## What remains a limitation

- `N50` and `delta` remain structurally confounded.
- The Maryland immunity mixture is latent and partly prior-carried.
- The Oxford `eta` correction and infection-definition terms are partly
  confounded; `frac_late` is expected to be weakly identified.
- The posterior's single divergence should be reviewed before treating the fit as
  diagnostically final.
- The Gibani rechallenge/susceptibility paradox and several unused outcomes remain
  outside this likelihood.

These are explicit limitations of the locked result, not unfinished wiring implied
by the project-completion status.

## Start here

1. [README](README.md) for the folder map and rebuild ladder.
2. [Model specification](dose_response_model_specification.md) for biology,
   equations, and data mapping.
3. [Joint inference plan](joint_inference_plan.md) for locked likelihood choices
   and identifiability assumptions.
4. [Tier 2 fit summary](calibration/results/t2-indiv-vax__phi-rho-eta-psi-vax/summary.md)
   for parameters, diagnostics, and figures.
