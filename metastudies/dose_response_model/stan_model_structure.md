# Stan model structure — `typhoid_dose_response.stan`

Updated 2026-07-31 from [typhoid_dose_response.stan](calibration/typhoid_dose_response.stan). Per-configuration observation counts: [calibration/TIER_LADDER.md](calibration/TIER_LADDER.md)
(Tier 2 machinery; C3 dose-dependent phi + issue-#15 cascade). Two views:

1. **Program dataflow** — how the Stan blocks feed each other.
2. **`obs_prob` dispatch** — the scientific core: how each observation
   group maps to a binomial success probability.

Core kernels:

- **beta-Poisson**: `P = 1 - (1 + D_eff·(2^(1/alpha) - 1)/N50)^(-alpha / CoP^gamma)`,
  with `D_eff = D/delta` (caller applies `/delta` where needed).
- **eta detection**: `eta_lo + (1 - eta_lo)·exp(-kappa·D_eff/N50_inf)`.
- **Maryland mixture**: `pi_susc·P(CoP_susc) + (1 - pi_susc)·P(CoP_imm)`.
- **phi0(T)** (low-dose definition sensitivity): `inv_logit(phi0_a - phi0_b·(T - T_ref))`.
- **phi(T,D)** (dose-dependent): `phi0(T) + (1 - phi0(T))·P_fev_naive(D_eff)`
  (beta_phi dose-lift exponent PINNED to 1; P_fev_naive = CoP=1 cascade fever curve).

## 1. Program dataflow

```mermaid
flowchart TD
    subgraph DATA["data"]
        OBS["Observations — N_obs flat rows<br/>group in 1..7 · n · y · dose · CoP · stratum · T_thresh"]
        LAD["Darton phi0 ladder<br/>ladder_T · ladder_count · ladder_N"]
        HYP["Prior hyperparameters<br/>(single source: priors.yaml)"]
        FLAG["prior_only flag · T_ref"]
    end

    subgraph PARS["parameters"]
        BIO["Biological (shared across studies)<br/>log10_N50_inf · d_fev<br/>alpha_inf · alpha_fevginf<br/>gamma_inf · gamma_fevginf"]
        NUIS["Nuisance<br/>log10_delta · pi_susc<br/>CoP_susc · CoP_imm (Exp prior)"]
        PHI["phi(T,D) shape<br/>phi0_a · phi0_b<br/>(beta_phi pinned = 1)"]
        ETAP["eta-correction<br/>eta_lo · kappa"]
        OD["(sigma_study DELETED 2026-07-31:<br/>inert at every config, no study index existed)"]
    end

    subgraph TP["transformed parameters"]
        REPARAM["log10_N50_fevginf = log10_N50_inf + d_fev<br/>(smooth lower-bound reparam: fever >= infection)<br/>N50_inf, N50_fevginf, delta = 10^(...)"]
        LPRIOR["lprior = sum of prior log-densities"]
    end

    subgraph MODEL["model"]
        TGT["target += lprior"]
        LIK["if prior_only == 0:<br/>y[i] ~ binomial(n[i], obs_prob(...))"]
        LADLIK["ladder_count[k] ~<br/>binomial(ladder_N, phi0(ladder_T[k]))"]
    end

    subgraph GQ["generated quantities"]
        PPC["p_pred · y_rep · log_lik (loo-ready)"]
        DERIV["Reference quantities:<br/>p_inf / p_fev at 1e3,1e4 naive ·<br/>Maryland fever curve · Hornick conditional ·<br/>phi0(38.3/39.4) · phi(39.4, dose) · eta · VE (placeholder)"]
    end

    subgraph FUN["functions (computational engine)"]
        OP["obs_prob() — unified per-obs success prob"]
        BP["beta_poisson()"]
        MM["maryland_mixture()"]
        ED["eta_detection()"]
        PHIF["phi_TD() / phi0_fn()"]
        OP --> BP
        OP --> MM
        OP --> ED
        OP --> PHIF
        MM --> BP
        PHIF --> BP
    end

    HYP --> LPRIOR
    BIO --> REPARAM
    NUIS --> REPARAM
    BIO --> LPRIOR
    NUIS --> LPRIOR
    PHI --> LPRIOR
    ETAP --> LPRIOR
    OD --> LPRIOR

    LPRIOR --> TGT
    TGT --> LIK
    FLAG --> LIK
    OBS --> LIK
    REPARAM --> LIK
    NUIS --> LIK
    PHI --> LIK
    ETAP --> LIK
    LAD --> LADLIK
    PHI --> LADLIK
    LIK -. calls .-> OP

    LIK --> PPC
    REPARAM --> DERIV
    PPC -. calls .-> OP
    DERIV -. calls .-> OP
```

## 2. `obs_prob` group dispatch (likelihood core)

```mermaid
flowchart TD
    START["obs_prob(group, dose, CoP, T_thresh, stratum, ..., T_ref, phi0_a, phi0_b)"]
    START --> G{"group?"}

    G -->|"1 · ox_fev"| OXF["D = dose (delta=1, no mixture)<br/>P_inf(CoP) x P_fev|inf(CoP)"]
    G -->|"2 · ox_inf"| OXI["D = dose<br/>eta_detection(D) x P_inf(CoP)"]
    G -->|"3 · md_fev"| MDF{"stratum?"}
    G -->|"4 · md_inf"| MDI["D = dose/delta<br/>maryland_mixture of P_inf<br/>pi·CoP_susc + (1-pi)·CoP_imm"]
    G -->|"5 · hornick_cond"| HC["D = dose/delta<br/>P(fever | infected) =<br/>phi(T,D) · mixture(P_fev) / mixture(P_inf)<br/>(guarded division)"]
    G -->|"6 · ox_inf_indiv"| OI6["D = dose (Oxford, no eta)<br/>P_inf(CoP) — individual infection"]
    G -->|"7 · ox_fevginf_indiv"| OF7["D = dose (Oxford)<br/>P_fev|inf(CoP) — individual fever|infection"]

    MDF -->|"1 · susceptible"| MS["phi(T,D) x P_inf(CoP_susc) x P_fev|inf(CoP_susc)"]
    MDF -->|"2 · immune"| MI["phi(T,D) x P_inf(CoP_imm) x P_fev|inf(CoP_imm)"]
    MDF -->|"0 · mixture"| MX["phi(T,D) x (pi·prod_susc + (1-pi)·prod_imm)"]
```

Groups 6/7 (the issue-#15 cascade) carry the Darton placebo per-subject endpoints:
6 = infection (`bact_or_stool`, all 30), 7 = fever|infection (`fever_td`, the 26
infected). Together they factorize the per-subject cascade and separate `gamma_inf`
from `gamma_fevginf` (weakly — bounded by Darton's thin anti-Vi titre range).
