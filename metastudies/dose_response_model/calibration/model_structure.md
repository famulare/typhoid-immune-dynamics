# Stan model structure — `typhoid_dose_response.stan`

Generated 2026-06-26 from [typhoid_dose_response.stan](typhoid_dose_response.stan)
(Tier 2: eta-correction, full complexity). Two views:

1. **Program dataflow** — how the Stan blocks feed each other.
2. **`obs_prob` dispatch** — the scientific core: how each observation
   group maps to a binomial success probability.

Core kernels:

- **beta-Poisson**: `P = 1 - (1 + D_eff·(2^(1/alpha) - 1)/N50)^(-alpha / CoP^gamma)`,
  with `D_eff = D/delta` (caller applies `/delta` where needed).
- **eta detection**: `eta_lo + (1 - eta_lo)·exp(-kappa·D_eff/N50_inf)`.
- **Maryland mixture**: `pi_susc·P(CoP_susc) + (1 - pi_susc)·P(CoP_imm)`.

## 1. Program dataflow

```mermaid
flowchart TD
    subgraph DATA["data"]
        OBS["Observations — N_obs flat rows<br/>group in 1..5 · n · y · dose · CoP · stratum"]
        HYP["Prior hyperparameters<br/>(single source: priors.yaml)"]
        FLAG["prior_only flag"]
    end

    subgraph PARS["parameters"]
        BIO["Biological (shared across studies)<br/>log10_N50_inf · d_fev<br/>alpha_inf · alpha_fevginf<br/>gamma_inf · gamma_fevginf"]
        NUIS["Nuisance<br/>log10_delta · pi_susc<br/>CoP_susc · CoP_imm · phi_md"]
        ETAP["eta-correction<br/>eta_lo · kappa"]
        OD["sigma_study (inert in Tier 1)"]
    end

    subgraph TP["transformed parameters"]
        REPARAM["log10_N50_fevginf = log10_N50_inf + d_fev<br/>(smooth lower-bound reparam: fever >= infection)<br/>N50_inf, N50_fevginf, delta = 10^(...)"]
        LPRIOR["lprior = sum of prior log-densities"]
    end

    subgraph MODEL["model"]
        TGT["target += lprior"]
        LIK["if prior_only == 0:<br/>y[i] ~ binomial(n[i], obs_prob(...))"]
    end

    subgraph GQ["generated quantities"]
        PPC["p_pred · y_rep · log_lik (loo-ready)"]
        DERIV["Reference-dose quantities:<br/>p_inf / p_fev at 1e3,1e4 naive ·<br/>Maryland fever curve · Hornick conditional ·<br/>eta(1e3,1e4) · VE (placeholder CoP)"]
    end

    subgraph FUN["functions (computational engine)"]
        OP["obs_prob() — unified per-obs success prob"]
        BP["beta_poisson()"]
        MM["maryland_mixture()"]
        ED["eta_detection()"]
        OP --> BP
        OP --> MM
        OP --> ED
        MM --> BP
    end

    HYP --> LPRIOR
    BIO --> REPARAM
    NUIS --> REPARAM
    BIO --> LPRIOR
    NUIS --> LPRIOR
    ETAP --> LPRIOR
    OD --> LPRIOR

    LPRIOR --> TGT
    TGT --> LIK
    FLAG --> LIK
    OBS --> LIK
    REPARAM --> LIK
    NUIS --> LIK
    ETAP --> LIK
    LIK -. calls .-> OP

    LIK --> PPC
    REPARAM --> DERIV
    PPC -. calls .-> OP
    DERIV -. calls .-> OP
```

## 2. `obs_prob` group dispatch (likelihood core)

```mermaid
flowchart TD
    START["obs_prob(group, dose, CoP, phi, stratum, ...)"]
    START --> G{"group?"}

    G -->|"1 · ox_fev"| OXF["D = dose (delta=1, no mixture)<br/>P_inf(CoP) x P_fev|inf(CoP)"]
    G -->|"2 · ox_inf"| OXI["D = dose<br/>eta_detection(D) x P_inf(CoP)"]
    G -->|"3 · md_fev"| MDF{"stratum?"}
    G -->|"4 · md_inf"| MDI["D = dose/delta<br/>maryland_mixture of P_inf<br/>pi·CoP_susc + (1-pi)·CoP_imm"]
    G -->|"5 · hornick_cond"| HC["D = dose/delta<br/>P(fever | infected) =<br/>phi · mixture(P_fev) / mixture(P_inf)<br/>(guarded division)"]

    MDF -->|"1 · susceptible"| MS["phi x P_inf(CoP_susc) x P_fev|inf(CoP_susc)"]
    MDF -->|"2 · immune"| MI["phi x P_inf(CoP_imm) x P_fev|inf(CoP_imm)"]
    MDF -->|"0 · mixture"| MX["phi x (pi·prod_susc + (1-pi)·prod_imm)"]
```
