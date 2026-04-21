// =============================================================================
// Typhoid Dose-Response Model: Tier 2 (η-correction, full complexity)
// =============================================================================
// UNTESTED DRAFT — skeleton for colleague review, not ready to run.
// Implements: cascaded beta-Poisson (infection × fever|infection),
//   cross-era δ bridge, Maryland mixture model, study-specific φ(T),
//   η-correction for Oxford shedding detection bias.
// Reference: joint_inference_plan.md Sections 2.1-2.7
// Authors: Mike Famulare, Claude (Opus 4.6)
// Date: 2026-03-24
// =============================================================================

functions {
  // Beta-Poisson dose-response: P(outcome | dose, N50, alpha, CoP, gamma)
  // Returns probability of outcome given effective dose and immunity.
  // D_eff = D / delta_medium (already converted by caller)
  real beta_poisson_lp(real D_eff, real N50, real alpha, real CoP, real gamma) {
    // Modified beta-Poisson with immunity scaling:
    //   P = 1 - (1 + D_eff * (2^(1/alpha) - 1) / N50) ^ (-alpha / CoP^gamma)
    real scale = (pow(2.0, 1.0 / alpha) - 1.0) / N50;
    real exponent = -alpha / pow(CoP, gamma);
    return 1.0 - pow(1.0 + D_eff * scale, exponent);
  }

  // Maryland mixture model: two-component (susceptible + immune)
  real maryland_mixture(real D_eff, real N50, real alpha, real gamma,
                        real pi_susc, real CoP_susc, real CoP_imm) {
    real p_susc = beta_poisson_lp(D_eff, N50, alpha, CoP_susc, gamma);
    real p_imm  = beta_poisson_lp(D_eff, N50, alpha, CoP_imm, gamma);
    return pi_susc * p_susc + (1.0 - pi_susc) * p_imm;
  }

  // η shedding detection probability (dose-dependent, Option A)
  // At low dose: η → 1 (plenty of time before diagnosis)
  // At high dose: η → η_lo (treatment truncates shedding detection)
  real eta_detection(real D_eff, real N50_inf, real eta_lo, real kappa) {
    return eta_lo + (1.0 - eta_lo) * exp(-kappa * D_eff / N50_inf);
  }
}

data {
  // ---- Oxford fever observations (δ = 1) ----
  int<lower=0> N_ox_fev;                    // number of Oxford fever groups
  array[N_ox_fev] int<lower=0> n_ox_fev;    // sample sizes
  array[N_ox_fev] int<lower=0> y_ox_fev;    // fever cases
  vector<lower=0>[N_ox_fev] dose_ox_fev;    // doses in CFU (bicarb frame)
  vector<lower=0>[N_ox_fev] CoP_ox_fev;     // group-average CoP (1=naive, >1=vaccinated)

  // ---- Oxford shedding observations (η-corrected, Tier 2) ----
  int<lower=0> N_ox_inf;                    // number of Oxford shedding groups
  array[N_ox_inf] int<lower=0> n_ox_inf;    // sample sizes
  array[N_ox_inf] int<lower=0> y_ox_inf;    // shedding cases
  vector<lower=0>[N_ox_inf] dose_ox_inf;    // doses in CFU
  vector<lower=0>[N_ox_inf] CoP_ox_inf;     // group-average CoP

  // ---- Maryland fever observations (δ > 1) ----
  // Hornick multi-dose + Gilman strata + Levine trials
  int<lower=0> N_md_fev;                    // number of Maryland fever groups
  array[N_md_fev] int<lower=0> n_md_fev;    // sample sizes
  array[N_md_fev] int<lower=0> y_md_fev;    // fever cases
  vector<lower=0>[N_md_fev] dose_md_fev;    // doses in CFU (milk frame, pre-δ)
  vector<lower=0>[N_md_fev] phi_md_fev;     // study-specific φ(T) for each obs
  // phi values: Hornick ≈ 0.25, Levine ≈ 0.65, Gilman ≈ 0.65
  // These are DATA, not parameters (fixed from Oxford threshold ladder)

  // Flag: is this a Gilman H-antibody stratified observation?
  // 0 = mixture model, 1 = susceptible stratum (CoP_susc), 2 = immune stratum (CoP_imm)
  array[N_md_fev] int<lower=0, upper=2> gilman_stratum;

  // ---- Maryland infection observations ----
  // Hornick Table 2 infection (10^7), Gilman shedding, Levine shedding
  int<lower=0> N_md_inf;
  array[N_md_inf] int<lower=0> n_md_inf;
  array[N_md_inf] int<lower=0> y_md_inf;
  vector<lower=0>[N_md_inf] dose_md_inf;

  // ---- Hornick Table 2 conditional: fever | infected at 10^7 ----
  int<lower=0> n_hornick_cond;   // = 28 (infected)
  int<lower=0> y_hornick_cond;   // = 16 (febrile among infected)
  real<lower=0> dose_hornick_7;  // = 1e7
  real<lower=0> phi_hornick;     // = 0.25 (Hornick φ for conditional)
}

parameters {
  // ---- Biological parameters (shared across all studies) ----
  // TO BE EXAMINED: all priors are defaults from plan Section 7

  real log10_N50_inf;      // log10 of infection N50 in bicarb-equivalent CFU
  real log10_N50_fevginf;  // log10 of fever|infection N50
  real<lower=0> alpha_inf;       // beta-Poisson heterogeneity (infection)
  real<lower=0> alpha_fevginf;   // beta-Poisson heterogeneity (fever|inf)
  real<lower=0> gamma_inf;       // immunity scaling exponent (infection)
  real<lower=0> gamma_fevginf;   // immunity scaling exponent (fever|inf)

  // ---- Nuisance parameters ----
  real log10_delta;               // log10 milk-to-bicarb dose offset
  real<lower=0, upper=1> pi_susc; // Maryland susceptible fraction
  real<lower=0> CoP_imm;          // Maryland immune fraction CoP (>1)
  real<lower=0> CoP_susc;         // Maryland susceptible fraction CoP (near 1)

  // ---- η-correction parameters (Tier 2, Option A) ----
  real<lower=0, upper=1> eta_lo;  // minimum shedding detection prob at high dose
  real<lower=0> kappa;            // dose-scaling for η

  // ---- Study-level overdispersion ----
  real<lower=0> sigma_study;      // study-level random effect SD
}

transformed parameters {
  real<lower=0> N50_inf = pow(10.0, log10_N50_inf);
  real<lower=0> N50_fevginf = pow(10.0, log10_N50_fevginf);
  real<lower=0> delta = pow(10.0, log10_delta);
}

model {
  // ===========================================================================
  // PRIORS — TO BE EXAMINED
  // All values from joint_inference_plan.md Section 7
  // ===========================================================================

  // -- Biological parameters --
  log10_N50_inf ~ normal(2.5, 1.0);       // ~300 bicarb CFU; TO BE EXAMINED
  log10_N50_fevginf ~ normal(2.8, 1.0);   // ~600 bicarb CFU; TO BE EXAMINED
  alpha_inf ~ lognormal(-1.5, 0.8);       // ~0.1-0.5; TO BE EXAMINED
  alpha_fevginf ~ lognormal(-1.5, 0.8);   // same; TO BE EXAMINED
  gamma_inf ~ lognormal(-0.5, 0.7);       // ~0.2-1.0; TO BE EXAMINED
  gamma_fevginf ~ lognormal(0.0, 0.7);    // ~0.3-2.0; TO BE EXAMINED

  // Soft ordering constraint: N50_fev|inf > N50_inf
  // (infection threshold lower than fever threshold)
  // Implemented as half-normal on the difference
  (log10_N50_fevginf - log10_N50_inf) ~ normal(0, 1) T[0, ];  // TO BE EXAMINED

  // -- Nuisance parameters --
  log10_delta ~ normal(3.5, 0.7);         // ~1000-30000x; TO BE EXAMINED
  pi_susc ~ beta(7, 4);                   // ~0.65; Gilman 36/53=68%; TO BE EXAMINED
  CoP_imm ~ lognormal(1.0, 0.5);          // ~2-5; TO BE EXAMINED
  CoP_susc ~ lognormal(0, 0.2);           // near 1; TO BE EXAMINED

  // -- η-correction parameters --
  eta_lo ~ beta(5, 5);                    // ~0.5; TO BE EXAMINED
  kappa ~ lognormal(0, 1);               // weakly informative; TO BE EXAMINED

  // -- Study-level overdispersion --
  sigma_study ~ normal(0, 0.5) T[0, ];   // half-normal; TO BE EXAMINED

  // ===========================================================================
  // LIKELIHOOD
  // ===========================================================================

  // ---- Oxford fever (δ = 1, no mixture) ----
  for (i in 1:N_ox_fev) {
    real D_eff = dose_ox_fev[i];  // bicarb frame, δ=1
    real CoP = CoP_ox_fev[i];

    // P(fever) = P(inf) × P(fev|inf)
    real p_inf = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP, gamma_inf);
    real p_fevginf = beta_poisson_lp(D_eff, N50_fevginf, alpha_fevginf,
                                      CoP, gamma_fevginf);
    real p_fev = p_inf * p_fevginf;

    y_ox_fev[i] ~ binomial(n_ox_fev[i], p_fev);
  }

  // ---- Oxford shedding with η-correction (Tier 2) ----
  for (i in 1:N_ox_inf) {
    real D_eff = dose_ox_inf[i];
    real CoP = CoP_ox_inf[i];

    real p_inf = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP, gamma_inf);
    real eta = eta_detection(D_eff, N50_inf, eta_lo, kappa);

    // Observed shedding = η × P(infection)
    y_ox_inf[i] ~ binomial(n_ox_inf[i], eta * p_inf);
  }

  // ---- Maryland fever (δ > 1, mixture model, study-specific φ) ----
  for (i in 1:N_md_fev) {
    real D_eff = dose_md_fev[i] / delta;  // convert milk dose to bicarb-equivalent
    real phi = phi_md_fev[i];             // study-specific definition sensitivity

    real p_fev;

    if (gilman_stratum[i] == 1) {
      // Gilman susceptible stratum: use CoP_susc directly (no mixture)
      real p_inf_s = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP_susc, gamma_inf);
      real p_fg_s  = beta_poisson_lp(D_eff, N50_fevginf, alpha_fevginf,
                                      CoP_susc, gamma_fevginf);
      p_fev = phi * p_inf_s * p_fg_s;
    } else if (gilman_stratum[i] == 2) {
      // Gilman immune stratum: use CoP_imm directly (no mixture)
      real p_inf_i = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP_imm, gamma_inf);
      real p_fg_i  = beta_poisson_lp(D_eff, N50_fevginf, alpha_fevginf,
                                      CoP_imm, gamma_fevginf);
      p_fev = phi * p_inf_i * p_fg_i;
    } else {
      // Standard Maryland mixture (Hornick, Levine, Gilman unstratified)
      real p_inf_mix = maryland_mixture(D_eff, N50_inf, alpha_inf, gamma_inf,
                                         pi_susc, CoP_susc, CoP_imm);
      // For fever, need mixture of the PRODUCT P_inf × P_fev|inf
      real p_fev_susc = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP_susc, gamma_inf)
                        * beta_poisson_lp(D_eff, N50_fevginf, alpha_fevginf,
                                           CoP_susc, gamma_fevginf);
      real p_fev_imm  = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP_imm, gamma_inf)
                        * beta_poisson_lp(D_eff, N50_fevginf, alpha_fevginf,
                                           CoP_imm, gamma_fevginf);
      p_fev = phi * (pi_susc * p_fev_susc + (1.0 - pi_susc) * p_fev_imm);
    }

    y_md_fev[i] ~ binomial(n_md_fev[i], p_fev);
  }

  // ---- Maryland infection (δ > 1, mixture model) ----
  for (i in 1:N_md_inf) {
    real D_eff = dose_md_inf[i] / delta;

    real p_inf = maryland_mixture(D_eff, N50_inf, alpha_inf, gamma_inf,
                                   pi_susc, CoP_susc, CoP_imm);

    y_md_inf[i] ~ binomial(n_md_inf[i], p_inf);
  }

  // ---- Hornick Table 2 conditional: P(fever | infected) at 10^7 ----
  {
    real D_eff = dose_hornick_7 / delta;

    // Infection probability (mixture)
    real p_inf = maryland_mixture(D_eff, N50_inf, alpha_inf, gamma_inf,
                                   pi_susc, CoP_susc, CoP_imm);

    // Fever probability (mixture of products)
    real p_fev_susc = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP_susc, gamma_inf)
                      * beta_poisson_lp(D_eff, N50_fevginf, alpha_fevginf,
                                         CoP_susc, gamma_fevginf);
    real p_fev_imm  = beta_poisson_lp(D_eff, N50_inf, alpha_inf, CoP_imm, gamma_inf)
                      * beta_poisson_lp(D_eff, N50_fevginf, alpha_fevginf,
                                         CoP_imm, gamma_fevginf);
    real p_fev_mix = phi_hornick * (pi_susc * p_fev_susc + (1.0 - pi_susc) * p_fev_imm);

    // Conditional: P(fever | infected) = P(fever) / P(infection)
    real p_cond = p_fev_mix / p_inf;

    // Infection marginal
    n_hornick_cond ~ binomial(30, p_inf);  // 28/30 infected (hard-coded n=30)
    y_hornick_cond ~ binomial(n_hornick_cond, p_cond);
  }
}

generated quantities {
  // Posterior predictive checks and derived quantities

  // Oxford predicted attack rates at reference doses (bicarb frame)
  real p_inf_1e3_naive = beta_poisson_lp(1e3, N50_inf, alpha_inf, 1.0, gamma_inf);
  real p_inf_1e4_naive = beta_poisson_lp(1e4, N50_inf, alpha_inf, 1.0, gamma_inf);
  real p_fev_1e3_naive = p_inf_1e3_naive
                         * beta_poisson_lp(1e3, N50_fevginf, alpha_fevginf, 1.0, gamma_fevginf);
  real p_fev_1e4_naive = p_inf_1e4_naive
                         * beta_poisson_lp(1e4, N50_fevginf, alpha_fevginf, 1.0, gamma_fevginf);

  // Maryland predicted curve at key doses (milk frame)
  real p_fev_md_1e3 = 0.25 * (pi_susc
    * beta_poisson_lp(1e3/delta, N50_inf, alpha_inf, CoP_susc, gamma_inf)
      * beta_poisson_lp(1e3/delta, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf)
    + (1-pi_susc)
    * beta_poisson_lp(1e3/delta, N50_inf, alpha_inf, CoP_imm, gamma_inf)
      * beta_poisson_lp(1e3/delta, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf));
  real p_fev_md_1e5 = 0.25 * (pi_susc
    * beta_poisson_lp(1e5/delta, N50_inf, alpha_inf, CoP_susc, gamma_inf)
      * beta_poisson_lp(1e5/delta, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf)
    + (1-pi_susc)
    * beta_poisson_lp(1e5/delta, N50_inf, alpha_inf, CoP_imm, gamma_inf)
      * beta_poisson_lp(1e5/delta, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf));
  real p_fev_md_1e7 = 0.25 * (pi_susc
    * beta_poisson_lp(1e7/delta, N50_inf, alpha_inf, CoP_susc, gamma_inf)
      * beta_poisson_lp(1e7/delta, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf)
    + (1-pi_susc)
    * beta_poisson_lp(1e7/delta, N50_inf, alpha_inf, CoP_imm, gamma_inf)
      * beta_poisson_lp(1e7/delta, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf));

  // Hornick Table 2 conditional prediction
  real p_cond_pred;
  {
    real D7 = dose_hornick_7 / delta;
    real p_i = maryland_mixture(D7, N50_inf, alpha_inf, gamma_inf,
                                 pi_susc, CoP_susc, CoP_imm);
    real pfs = beta_poisson_lp(D7, N50_inf, alpha_inf, CoP_susc, gamma_inf)
               * beta_poisson_lp(D7, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf);
    real pfi = beta_poisson_lp(D7, N50_inf, alpha_inf, CoP_imm, gamma_inf)
               * beta_poisson_lp(D7, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf);
    p_cond_pred = phi_hornick * (pi_susc * pfs + (1-pi_susc) * pfi) / p_i;
  }

  // η at reference Oxford doses
  real eta_1e3 = eta_detection(1e3, N50_inf, eta_lo, kappa);
  real eta_1e4 = eta_detection(1e4, N50_inf, eta_lo, kappa);

  // Vaccine efficacy predictions (Jin-like: CoP from Vi-TT GMT 563)
  // NOTE: requires CoP mapping g(anti-Vi) which is NOT in this model yet.
  // Placeholder: assume CoP_ViTT ~ 5, CoP_ViPS ~ 2 (TO BE REPLACED)
  real VE_fev_ViTT;
  real VE_fev_ViPS;
  {
    real CoP_ViTT = 5.0;  // PLACEHOLDER — TO BE REPLACED with titer model
    real CoP_ViPS = 2.0;  // PLACEHOLDER — TO BE REPLACED
    real p_fev_ctrl = p_fev_1e4_naive;  // approximate
    real p_fev_vitt = beta_poisson_lp(2e4, N50_inf, alpha_inf, CoP_ViTT, gamma_inf)
                      * beta_poisson_lp(2e4, N50_fevginf, alpha_fevginf,
                                         CoP_ViTT, gamma_fevginf);
    real p_fev_vips = beta_poisson_lp(2e4, N50_inf, alpha_inf, CoP_ViPS, gamma_inf)
                      * beta_poisson_lp(2e4, N50_fevginf, alpha_fevginf,
                                         CoP_ViPS, gamma_fevginf);
    VE_fev_ViTT = 1.0 - p_fev_vitt / p_fev_ctrl;
    VE_fev_ViPS = 1.0 - p_fev_vips / p_fev_ctrl;
  }

  // Infection-fever gap: P(inf) - P(fev) at reference doses
  real inf_fev_gap_1e4 = p_inf_1e4_naive - p_fev_1e4_naive;

  // Delta in human-interpretable units
  real delta_fold = delta;  // milk dose / bicarb dose equivalence
}
