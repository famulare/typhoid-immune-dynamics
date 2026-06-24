// =============================================================================
// Typhoid Dose-Response Model: Tier 2 (eta-correction, full complexity)
// =============================================================================
// Refactored 2026-06-23 (workflow upgrade; see calibration/CALIBRATION_WORKFLOW.md
// and plan please-plan-a-b-smooth-ritchie.md). Changes vs the prior draft:
//   - Unified obs_prob() helper: ONE source for every observation's success
//     probability, called by both the model likelihood and generated quantities.
//   - Flat per-observation data layout (N_obs rows + group code), replacing the
//     five per-group data blocks.
//   - N50 reparameterization: log10_N50_fevginf = log10_N50_inf + d_fev with
//     d_fev<lower=0>. This replaces the T[0,] density CLIFF on the difference
//     (the diagnosed cause of ~99% divergences) with Stan's smooth lower-bound
//     transform. The PRIOR is preserved exactly: all three original N50 prior
//     terms are retained (change of variables Jacobian = 1), so a drop in
//     divergences is attributable to geometry alone, not a changed prior.
//   - lprior accumulator (brms idiom): every prior written once; total prior in
//     `lprior` for priorsense power-scaling. Prior hyperparameters are DATA,
//     sourced from calibration/priors.yaml (single source of truth).
//   - generated quantities emit per-observation p_pred, y_rep, log_lik (PPC from
//     the model, loo-ready) instead of an R-side likelihood mirror.
// Update 2026-06-23: phi (Maryland fever definition-sensitivity) is now an
//   ESTIMATED scalar parameter phi_md ~ Beta (plan Section 7), NOT fixed per-obs
//   data. The old fixed phi=0.25 was the low-dose asymptote of phi(T) applied
//   dose-wide; it capped Maryland fitted fever at 0.25 while Hornick high-dose
//   data are 0.89-0.95 (structurally unfittable). A single floated scalar is
//   identified by the Hornick dose-range plateau and covaries with delta.
// Implements: cascaded beta-Poisson (infection x fever|infection),
//   cross-era delta bridge, Maryland mixture, estimated scalar phi_md,
//   eta-correction for Oxford shedding detection bias.
// Reference: joint_inference_plan.md Sections 2.1-2.7, Section 7 (priors)
// Authors: Mike Famulare, Claude (Opus 4.6 draft; Opus 4.8 refactor)
// =============================================================================

functions {
  // Beta-Poisson dose-response: P(outcome | dose, N50, alpha, CoP, gamma)
  //   P = 1 - (1 + D_eff * (2^(1/alpha) - 1) / N50) ^ (-alpha / CoP^gamma)
  // D_eff = D / delta_medium (already converted by caller).
  real beta_poisson(real D_eff, real N50, real alpha, real CoP, real gamma) {
    real scale = (pow(2.0, 1.0 / alpha) - 1.0) / N50;
    real exponent = -alpha / pow(CoP, gamma);
    return 1.0 - pow(1.0 + D_eff * scale, exponent);
  }

  // Maryland two-component (susceptible + immune) mixture.
  real maryland_mixture(real D_eff, real N50, real alpha, real gamma,
                        real pi_susc, real CoP_susc, real CoP_imm) {
    real p_susc = beta_poisson(D_eff, N50, alpha, CoP_susc, gamma);
    real p_imm  = beta_poisson(D_eff, N50, alpha, CoP_imm, gamma);
    return pi_susc * p_susc + (1.0 - pi_susc) * p_imm;
  }

  // eta shedding-detection probability (dose-dependent, Option A).
  real eta_detection(real D_eff, real N50_inf, real eta_lo, real kappa) {
    return eta_lo + (1.0 - eta_lo) * exp(-kappa * D_eff / N50_inf);
  }

  // ---- Unified observation probability -------------------------------------
  // ONE place that turns an observation into its binomial success probability.
  // group: 1=ox_fev  2=ox_inf  3=md_fev  4=md_inf  5=hornick_cond
  // Covariates: dose (raw CFU), CoP (group-average, used by ox only),
  //             phi (definition sensitivity, md_fev/hornick only),
  //             stratum (gilman: 0=mixture, 1=susceptible, 2=immune).
  // Unused covariates take harmless defaults (CoP=phi=1, stratum=0) from the
  // caller; the math below ignores them per-group.
  real obs_prob(int group, real dose, real CoP, real phi, int stratum,
                real N50_inf, real N50_fevginf,
                real alpha_inf, real alpha_fevginf,
                real gamma_inf, real gamma_fevginf,
                real delta, real pi_susc, real CoP_susc, real CoP_imm,
                real eta_lo, real kappa) {
    if (group == 1) {                                   // ox_fev (delta=1, no mixture)
      real D = dose;
      return beta_poisson(D, N50_inf, alpha_inf, CoP, gamma_inf)
           * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP, gamma_fevginf);
    } else if (group == 2) {                            // ox_inf (eta-corrected shedding)
      real D = dose;
      real p_inf = beta_poisson(D, N50_inf, alpha_inf, CoP, gamma_inf);
      return eta_detection(D, N50_inf, eta_lo, kappa) * p_inf;
    } else if (group == 3) {                            // md_fev (delta>1, mixture/strata, *phi)
      real D = dose / delta;
      if (stratum == 1) {                               // Gilman susceptible stratum
        return phi * beta_poisson(D, N50_inf, alpha_inf, CoP_susc, gamma_inf)
                   * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf);
      } else if (stratum == 2) {                        // Gilman immune stratum
        return phi * beta_poisson(D, N50_inf, alpha_inf, CoP_imm, gamma_inf)
                   * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf);
      } else {                                          // mixture of the P_inf x P_fev|inf product
        real p_fev_susc = beta_poisson(D, N50_inf, alpha_inf, CoP_susc, gamma_inf)
                        * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf);
        real p_fev_imm  = beta_poisson(D, N50_inf, alpha_inf, CoP_imm, gamma_inf)
                        * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf);
        return phi * (pi_susc * p_fev_susc + (1.0 - pi_susc) * p_fev_imm);
      }
    } else if (group == 4) {                            // md_inf (delta>1, mixture)
      return maryland_mixture(dose / delta, N50_inf, alpha_inf, gamma_inf,
                              pi_susc, CoP_susc, CoP_imm);
    } else {                                            // group == 5 hornick_cond: P(fever | infected)
      real D = dose / delta;
      real p_inf = maryland_mixture(D, N50_inf, alpha_inf, gamma_inf,
                                    pi_susc, CoP_susc, CoP_imm);
      real p_fev_susc = beta_poisson(D, N50_inf, alpha_inf, CoP_susc, gamma_inf)
                      * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf);
      real p_fev_imm  = beta_poisson(D, N50_inf, alpha_inf, CoP_imm, gamma_inf)
                      * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf);
      real p_fev_mix = phi * (pi_susc * p_fev_susc + (1.0 - pi_susc) * p_fev_imm);
      real p_cond = p_fev_mix / p_inf;
      return fmin(fmax(p_cond, 1e-12), 1.0 - 1e-12);    // guard the division
    }
  }
}

data {
  // ---- Flat observation layout --------------------------------------------
  int<lower=0> N_obs;
  array[N_obs] int<lower=1, upper=5> group;   // 1=ox_fev 2=ox_inf 3=md_fev 4=md_inf 5=hornick_cond
  array[N_obs] int<lower=0> n;                 // sample sizes
  array[N_obs] int<lower=0> y;                 // events
  vector<lower=0>[N_obs] dose;                 // raw dose in CFU (helper applies /delta where needed)
  // CoP = correlate of protection. DEFINITIONAL UNITS: anti-Vi IgG titre on the
  // commercial VaccZyme ELISA scale (The Binding Site; EU/mL, LLD 7.4), expressed
  // RELATIVE TO THE NAIVE REFERENCE so CoP=1 at naive (<LLD, imputed ~3.7 EU/mL).
  // Both modern Oxford inputs (Jin 2017, Darton 2016) use this assay (verified).
  // It enters the dose-response as CoP^gamma — a per-log10-titre power law — so
  // gamma is the protection slope anchorable to Darton's HR 0.29/log10 anti-Vi.
  // NOTE: current per-group values are INTERIM placeholders (Jin 5.0/2.0 etc.)
  // pending the EU/mL value-swap; see tier1_lab_notebook.md D1.
  vector<lower=0>[N_obs] CoP;                  // group-average CoP (anti-Vi/naive, VaccZyme EU/mL; 1 elsewhere)
  array[N_obs] int<lower=0, upper=2> stratum;  // gilman stratum (md_fev; 0 elsewhere)
  // NOTE: phi is no longer data — it is the estimated scalar parameter phi_md below.

  // ---- Control flag --------------------------------------------------------
  int<lower=0, upper=1> prior_only;            // 1 = skip likelihood (prior predictive)

  // ---- Prior hyperparameters (DATA; single source = calibration/priors.yaml)
  // Change a value here (via priors.yaml) and refit -- no recompile needed.
  real pr_log10_N50_inf_mu;      real<lower=0> pr_log10_N50_inf_sd;
  real pr_log10_N50_fevginf_mu;  real<lower=0> pr_log10_N50_fevginf_sd;
  real pr_d_fev_mu;              real<lower=0> pr_d_fev_sd;     // half-normal (d_fev>=0)
  real pr_alpha_inf_mu;          real<lower=0> pr_alpha_inf_sd;
  real pr_alpha_fevginf_mu;      real<lower=0> pr_alpha_fevginf_sd;
  real pr_gamma_inf_mu;          real<lower=0> pr_gamma_inf_sd;
  real pr_gamma_fevginf_mu;      real<lower=0> pr_gamma_fevginf_sd;
  real pr_log10_delta_mu;        real<lower=0> pr_log10_delta_sd;
  real<lower=0> pr_pi_susc_a;    real<lower=0> pr_pi_susc_b;    // beta
  real pr_CoP_imm_mu;            real<lower=0> pr_CoP_imm_sd;
  real pr_CoP_susc_mu;           real<lower=0> pr_CoP_susc_sd;
  real<lower=0> pr_phi_md_a;     real<lower=0> pr_phi_md_b;     // beta (Maryland definition-sensitivity)
  real<lower=0> pr_eta_lo_a;     real<lower=0> pr_eta_lo_b;     // beta
  real pr_kappa_mu;              real<lower=0> pr_kappa_sd;
  real pr_sigma_study_mu;        real<lower=0> pr_sigma_study_sd; // half-normal (sigma_study>=0)
}

parameters {
  // ---- Biological parameters (shared across all studies) ----
  real log10_N50_inf;             // log10 infection N50 (bicarb-equivalent CFU)
  real<lower=0> d_fev;            // log10 gap: fever threshold ABOVE infection (reparam offset)
  real<lower=0> alpha_inf;        // beta-Poisson heterogeneity (infection)
  real<lower=0> alpha_fevginf;    // beta-Poisson heterogeneity (fever|inf)
  real<lower=0> gamma_inf;        // immunity scaling exponent (infection)
  real<lower=0> gamma_fevginf;    // immunity scaling exponent (fever|inf)

  // ---- Nuisance parameters ----
  real log10_delta;               // log10 milk-to-bicarb dose offset
  real<lower=0, upper=1> pi_susc; // Maryland susceptible fraction
  // Latent Maryland CoP (anti-Vi not measured pre-VaccZyme); same units as CoP
  // above — anti-Vi-equivalent titre relative to naive (VaccZyme EU/mL).
  real<lower=0> CoP_imm;          // Maryland immune-component CoP (>1)
  real<lower=0> CoP_susc;         // Maryland susceptible-component CoP (near 1)
  real<lower=0, upper=1> phi_md;  // Maryland fever definition-sensitivity (estimated scalar; was fixed 0.25/0.65)

  // ---- eta-correction parameters (Tier 2, Option A) ----
  real<lower=0, upper=1> eta_lo;  // minimum shedding detection prob at high dose
  real<lower=0> kappa;            // dose-scaling for eta

  // ---- Study-level overdispersion ----
  real<lower=0> sigma_study;      // study-level random effect SD (INERT in Tier 1)
}

transformed parameters {
  real log10_N50_fevginf = log10_N50_inf + d_fev;   // reparam: fever threshold >= infection
  real<lower=0> N50_inf = pow(10.0, log10_N50_inf);
  real<lower=0> N50_fevginf = pow(10.0, log10_N50_fevginf);
  real<lower=0> delta = pow(10.0, log10_delta);

  // ---- lprior accumulator (priors written ONCE; hyperparameters from data) --
  // Reproduces the original three-term N50 prior exactly under the (Jacobian=1)
  // change of variables (inf, fevginf) -> (inf, d_fev).
  real lprior = 0;
  lprior += normal_lpdf(log10_N50_inf     | pr_log10_N50_inf_mu,     pr_log10_N50_inf_sd);
  lprior += normal_lpdf(log10_N50_fevginf | pr_log10_N50_fevginf_mu, pr_log10_N50_fevginf_sd);
  lprior += normal_lpdf(d_fev             | pr_d_fev_mu,             pr_d_fev_sd);       // half-normal via lower=0
  lprior += lognormal_lpdf(alpha_inf      | pr_alpha_inf_mu,         pr_alpha_inf_sd);
  lprior += lognormal_lpdf(alpha_fevginf  | pr_alpha_fevginf_mu,     pr_alpha_fevginf_sd);
  lprior += lognormal_lpdf(gamma_inf      | pr_gamma_inf_mu,         pr_gamma_inf_sd);
  lprior += lognormal_lpdf(gamma_fevginf  | pr_gamma_fevginf_mu,     pr_gamma_fevginf_sd);
  lprior += normal_lpdf(log10_delta       | pr_log10_delta_mu,       pr_log10_delta_sd);
  lprior += beta_lpdf(pi_susc             | pr_pi_susc_a,            pr_pi_susc_b);
  lprior += lognormal_lpdf(CoP_imm        | pr_CoP_imm_mu,           pr_CoP_imm_sd);
  lprior += lognormal_lpdf(CoP_susc       | pr_CoP_susc_mu,          pr_CoP_susc_sd);
  lprior += beta_lpdf(phi_md              | pr_phi_md_a,             pr_phi_md_b);
  lprior += beta_lpdf(eta_lo              | pr_eta_lo_a,             pr_eta_lo_b);
  lprior += lognormal_lpdf(kappa          | pr_kappa_mu,             pr_kappa_sd);
  lprior += normal_lpdf(sigma_study       | pr_sigma_study_mu,       pr_sigma_study_sd); // half-normal via lower=0
}

model {
  target += lprior;

  if (prior_only == 0) {
    for (i in 1:N_obs) {
      y[i] ~ binomial(n[i], obs_prob(group[i], dose[i], CoP[i], phi_md, stratum[i],
                                     N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                                     gamma_inf, gamma_fevginf, delta, pi_susc,
                                     CoP_susc, CoP_imm, eta_lo, kappa));
    }
  }
}

generated quantities {
  // ---- Per-observation posterior predictive + pointwise log-likelihood -----
  // p_pred: fitted success probability (same obs_prob as the likelihood -> no
  // R mirror); y_rep: replicate counts; log_lik: for loo (with correct unit
  // grouping applied downstream -- the Hornick marginal H-I-7 and conditional
  // hornick_cond are one table factorized; combine them into one loo unit).
  vector[N_obs] p_pred;
  array[N_obs] int y_rep;
  vector[N_obs] log_lik;
  for (i in 1:N_obs) {
    p_pred[i]  = obs_prob(group[i], dose[i], CoP[i], phi_md, stratum[i],
                          N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                          gamma_inf, gamma_fevginf, delta, pi_susc,
                          CoP_susc, CoP_imm, eta_lo, kappa);
    y_rep[i]   = binomial_rng(n[i], p_pred[i]);
    log_lik[i] = binomial_lpmf(y[i] | n[i], p_pred[i]);
  }

  // ---- Reference-dose derived quantities (naive Oxford, bicarb frame) -------
  real p_inf_1e3_naive = beta_poisson(1e3, N50_inf, alpha_inf, 1.0, gamma_inf);
  real p_inf_1e4_naive = beta_poisson(1e4, N50_inf, alpha_inf, 1.0, gamma_inf);
  real p_fev_1e3_naive = p_inf_1e3_naive
                         * beta_poisson(1e3, N50_fevginf, alpha_fevginf, 1.0, gamma_fevginf);
  real p_fev_1e4_naive = p_inf_1e4_naive
                         * beta_poisson(1e4, N50_fevginf, alpha_fevginf, 1.0, gamma_fevginf);

  // Maryland predicted fever curve at key doses (milk frame; estimated phi_md)
  real p_fev_md_1e3 = obs_prob(3, 1e3, 1.0, phi_md, 0, N50_inf, N50_fevginf,
                               alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                               delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa);
  real p_fev_md_1e5 = obs_prob(3, 1e5, 1.0, phi_md, 0, N50_inf, N50_fevginf,
                               alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                               delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa);
  real p_fev_md_1e7 = obs_prob(3, 1e7, 1.0, phi_md, 0, N50_inf, N50_fevginf,
                               alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                               delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa);

  // Hornick Table 2 conditional prediction (estimated phi_md)
  real p_cond_pred = obs_prob(5, 1e7, 1.0, phi_md, 0, N50_inf, N50_fevginf,
                              alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                              delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa);

  // eta at reference Oxford doses
  real eta_1e3 = eta_detection(1e3, N50_inf, eta_lo, kappa);
  real eta_1e4 = eta_detection(1e4, N50_inf, eta_lo, kappa);

  // Vaccine efficacy predictions (PLACEHOLDER CoP; requires titer model g(anti-Vi))
  real VE_fev_ViTT;
  real VE_fev_ViPS;
  {
    real CoP_ViTT = 5.0;  // PLACEHOLDER
    real CoP_ViPS = 2.0;  // PLACEHOLDER
    real p_fev_ctrl = p_fev_1e4_naive;
    real p_fev_vitt = beta_poisson(2e4, N50_inf, alpha_inf, CoP_ViTT, gamma_inf)
                      * beta_poisson(2e4, N50_fevginf, alpha_fevginf, CoP_ViTT, gamma_fevginf);
    real p_fev_vips = beta_poisson(2e4, N50_inf, alpha_inf, CoP_ViPS, gamma_inf)
                      * beta_poisson(2e4, N50_fevginf, alpha_fevginf, CoP_ViPS, gamma_fevginf);
    VE_fev_ViTT = 1.0 - p_fev_vitt / p_fev_ctrl;
    VE_fev_ViPS = 1.0 - p_fev_vips / p_fev_ctrl;
  }

  real inf_fev_gap_1e4 = p_inf_1e4_naive - p_fev_1e4_naive;
  real delta_fold = delta;
}
