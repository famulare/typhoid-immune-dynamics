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
// Update 2026-06-23: phi (Maryland fever definition-sensitivity) was promoted
//   from fixed per-obs data (0.25/0.65) to an ESTIMATED scalar phi_md ~ Beta. The
//   old fixed phi=0.25 was the low-dose asymptote of phi(T) applied dose-wide; it
//   capped Maryland fitted fever at 0.25 while Hornick high-dose data are 0.89-0.95.
// Update 2026-07-15 (C3): phi is now DOSE-DEPENDENT phi(T, D_eff). A dose-constant
//   phi(T) re-introduces the high-dose Hornick cap (plugging the Darton-ladder
//   phi0(39.4)~0.30 caps H-F-9 at ~0.27 vs observed 0.95). Form (plan Sec 2.5):
//     phi(T,D) = phi0(T) + (1 - phi0(T)) * P_fev_naive(D_eff)^beta_phi
//   phi0(T) = inv_logit(phi0_a - phi0_b*(T - T_ref)) is the LOW-DOSE definition
//   sensitivity (monotone-decreasing in threshold T), pinned by the Darton placebo
//   temperature ladder (>=38/38.5/39 among the 20 TD+ subjects = 16/10/8) as a
//   decoupled 2-param binomial sub-likelihood. The dose-lift exponent beta_phi is
//   PINNED to 1 (2026-07-15): the C3 fit found it prior-dominated (4 Hornick points
//   can't identify the transition shape), so phi rises exactly with the naive fever
//   curve -- the value is the correct limits + inherited scale, not a fitted shape.
//   Approximation: the ladder pins phi0 treating Darton's dose as ~low. Retires phi_md.
// Update 2026-07-15 (+cascade, issue #15): Darton placebo per-subject endpoints enter
//   as a proper cascade -- group 6 (ox_inf_indiv, P_inf at subject anti-Vi) for all
//   30, group 7 (ox_fevginf_indiv, P_fev|inf) for the infected -- replacing the C1
//   composite-fever rows. Splits gamma_inf vs gamma_fevginf via the titre spread.
// Update 2026-07-31 (Step 2, LOCKED cohort_random_effects_design.md): ONE beta-binomial
//   overdispersion parameter, grand_overdispersion_rho (the ICC; rho=0 IS the binomial),
//   with concentration grand_concentration_k = (1-rho)/rho derived. Replaces sigma_study
//   (deleted 2026-07-31: declared, used in ZERO likelihood terms at any tier). A per-cohort
//   random effect was evaluated and WITHDRAWN (16 cohorts / 24 group-level obs, 10
//   singletons; Hornick's 5 cohorts ARE the dose ladder, so a free per-cohort offset
//   competes with N50_inf/alpha_inf for the same variance) -- see the design doc. Applied
//   to every flat obs_prob() row (beta_binomial(n, p*k, (1-p)*k) is EXACTLY binomial(n,p)
//   at n=1, so this is a no-op for the Darton individual rows without a branch); the phi0
//   ladder binomial is deliberately NOT overdispersed (a different sub-model, not the
//   Maryland replication this parameter targets).
// Implements: cascaded beta-Poisson (infection x fever|infection),
//   cross-era delta bridge, Maryland mixture, dose-dependent phi(T,D),
//   individual-subject cascade endpoints, one beta-binomial overdispersion parameter,
//   eta-correction for Oxford shedding bias.
// Reference: ../joint_inference_plan.md Sections 2.1-2.7, Section 7 (priors)
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

  // Low-dose definition sensitivity phi0(T): fraction of TD cases crossing the
  // strict threshold T at low dose. Monotone-decreasing in T (phi0_b >= 0),
  // logit-linear, centered at T_ref. Pinned by the Darton placebo ladder.
  real phi0_fn(real T, real T_ref, real phi0_a, real phi0_b) {
    return inv_logit(phi0_a - phi0_b * (T - T_ref));
  }

  // Dose-dependent definition sensitivity phi(T, D_eff) (plan Sec 2.5). Rises
  // from phi0(T) toward 1 as the naive fever probability saturates with dose.
  // P_fev_naive is the CoP=1 (naive) cascade fever curve at D_eff, so immunity
  // is NOT double-counted here (it already acts upstream via CoP^gamma).
  // beta_phi PINNED to 1 [Mike 2026-07-15]: the C3 fit showed the dose-lift shape
  // exponent prior-dominated (data can't identify it from 4 Hornick points), so we
  // keep the simplest identified form -- phi rises exactly with the naive fever
  // curve. The value is the correct limits (phi0 at low dose -> 1 at saturation)
  // and the inherited dose scale, not a fitted transition shape.
  real phi_TD(real T, real D_eff, real T_ref, real phi0_a, real phi0_b,
              real N50_inf, real N50_fevginf, real alpha_inf, real alpha_fevginf,
              real gamma_inf, real gamma_fevginf) {
    real p_fev_naive = beta_poisson(D_eff, N50_inf, alpha_inf, 1.0, gamma_inf)
                     * beta_poisson(D_eff, N50_fevginf, alpha_fevginf, 1.0, gamma_fevginf);
    real phi0 = phi0_fn(T, T_ref, phi0_a, phi0_b);
    return phi0 + (1.0 - phi0) * p_fev_naive;
  }

  // ---- Unified observation probability -------------------------------------
  // ONE place that turns an observation into its binomial success probability.
  // group: 1=ox_fev  2=ox_inf  3=md_fev  4=md_inf  5=hornick_cond
  // Covariates: dose (raw CFU), CoP (group-average, used by ox only),
  //             T_thresh (fever threshold in degC, md_fev/hornick only),
  //             stratum (gilman: 0=mixture, 1=susceptible, 2=immune).
  // phi(T,D) is computed internally for md_fev/hornick from (phi0_a,phi0_b,
  // beta_phi,T_ref); unused covariates take harmless defaults from the caller.
  real obs_prob(int group, real dose, real CoP, real T_thresh, int stratum,
                real N50_inf, real N50_fevginf,
                real alpha_inf, real alpha_fevginf,
                real gamma_inf, real gamma_fevginf,
                real delta, real pi_susc, real CoP_susc, real CoP_imm,
                real eta_lo, real kappa,
                real T_ref, real phi0_a, real phi0_b) {
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
      real phi = phi_TD(T_thresh, D, T_ref, phi0_a, phi0_b,
                        N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                        gamma_inf, gamma_fevginf);
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
    } else if (group == 5) {                            // hornick_cond: P(fever | infected)
      real D = dose / delta;
      real phi = phi_TD(T_thresh, D, T_ref, phi0_a, phi0_b,
                        N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                        gamma_inf, gamma_fevginf);
      real p_inf = maryland_mixture(D, N50_inf, alpha_inf, gamma_inf,
                                    pi_susc, CoP_susc, CoP_imm);
      real p_fev_susc = beta_poisson(D, N50_inf, alpha_inf, CoP_susc, gamma_inf)
                      * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_susc, gamma_fevginf);
      real p_fev_imm  = beta_poisson(D, N50_inf, alpha_inf, CoP_imm, gamma_inf)
                      * beta_poisson(D, N50_fevginf, alpha_fevginf, CoP_imm, gamma_fevginf);
      real p_fev_mix = phi * (pi_susc * p_fev_susc + (1.0 - pi_susc) * p_fev_imm);
      real p_cond = p_fev_mix / p_inf;
      return fmin(fmax(p_cond, 1e-12), 1.0 - 1e-12);    // guard the division
    } else if (group == 6) {                            // ox_inf_indiv: individual Oxford infection
      return beta_poisson(dose, N50_inf, alpha_inf, CoP, gamma_inf);   // P_inf at subject CoP (no eta/delta)
    } else {                                            // group == 7 ox_fevginf_indiv: P(fever | infected)
      return beta_poisson(dose, N50_fevginf, alpha_fevginf, CoP, gamma_fevginf);
    }
  }
}

data {
  // ---- Flat observation layout --------------------------------------------
  int<lower=0> N_obs;
  array[N_obs] int<lower=1, upper=7> group;   // 1=ox_fev 2=ox_inf 3=md_fev 4=md_inf 5=hornick_cond 6=ox_inf_indiv 7=ox_fevginf_indiv
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
  // Per-obs strict fever threshold (degC) for the phi(T,D) definition map: Hornick
  // 39.4, Levine/Gilman 38.3, T_ref elsewhere (unused where phi is not applied).
  vector<lower=0>[N_obs] T_thresh;
  real T_ref;                                  // reference threshold (Oxford composite, 38.0 degC)

  // ---- Darton placebo temperature ladder (identifies phi0(T)) --------------
  // Among the N_TD TD+ placebo subjects, how many crossed each strict threshold.
  // A decoupled binomial sub-likelihood: count[k] ~ Binomial(N_TD, phi0(T[k])).
  int<lower=0> N_ladder;
  vector<lower=0>[N_ladder] ladder_T;          // thresholds (degC)
  array[N_ladder] int<lower=0> ladder_count;   // TD+ subjects crossing each threshold
  int<lower=0> ladder_N;                        // number of TD+ placebo subjects (denominator)

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
  real<lower=0> pr_CoP_imm_rate;                                // exponential (Maryland immune anti-Vi-equiv)
  real pr_CoP_susc_mu;           real<lower=0> pr_CoP_susc_sd;
  real pr_phi0_a_mu;             real<lower=0> pr_phi0_a_sd;    // normal (logit phi0 at T_ref)
  real pr_phi0_b_mu;             real<lower=0> pr_phi0_b_sd;    // half-normal (logit slope per degC, >=0)
  real<lower=0> pr_eta_lo_a;     real<lower=0> pr_eta_lo_b;     // beta
  real pr_kappa_mu;              real<lower=0> pr_kappa_sd;
  real<lower=0> pr_grand_overdispersion_rho_a;   // beta (ICC; rho=0 is the binomial)
  real<lower=0> pr_grand_overdispersion_rho_b;
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
  real<lower=0> CoP_imm;          // Maryland immune-component CoP (anti-Vi-equiv; Exp prior)
  real<lower=0> CoP_susc;         // Maryland susceptible-component CoP (near 1)
  // ---- Dose-dependent Maryland fever definition-sensitivity phi(T,D) ----
  // (beta_phi dose-lift exponent PINNED to 1: prior-dominated in the C3 fit.)
  real phi0_a;                    // logit phi0 at T_ref
  real<lower=0> phi0_b;           // logit decay per degC (monotone-decreasing in threshold)

  // ---- eta-correction parameters (Tier 2, Option A) ----
  real<lower=0, upper=1> eta_lo;  // minimum shedding detection prob at high dose
  real<lower=0> kappa;            // dose-scaling for eta

  // ---- Overdispersion (Step 2, LOCKED 2026-07-31) ----
  real<lower=0, upper=1> grand_overdispersion_rho;  // beta-binomial ICC (0 = binomial)
}

transformed parameters {
  real log10_N50_fevginf = log10_N50_inf + d_fev;   // reparam: fever threshold >= infection
  real<lower=0> N50_inf = pow(10.0, log10_N50_inf);
  real<lower=0> N50_fevginf = pow(10.0, log10_N50_fevginf);
  real<lower=0> delta = pow(10.0, log10_delta);
  // beta-binomial concentration (Step 2). rho -> 0 drives k -> inf, recovering the
  // binomial exactly; rho is the ICC / design-effect driver (1 + (n-1)*rho).
  real<lower=0> grand_concentration_k = (1.0 - grand_overdispersion_rho) / grand_overdispersion_rho;

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
  lprior += exponential_lpdf(CoP_imm      | pr_CoP_imm_rate);
  lprior += lognormal_lpdf(CoP_susc       | pr_CoP_susc_mu,          pr_CoP_susc_sd);
  lprior += normal_lpdf(phi0_a            | pr_phi0_a_mu,            pr_phi0_a_sd);
  lprior += normal_lpdf(phi0_b            | pr_phi0_b_mu,            pr_phi0_b_sd);   // half-normal via lower=0
  lprior += beta_lpdf(eta_lo              | pr_eta_lo_a,             pr_eta_lo_b);
  lprior += lognormal_lpdf(kappa          | pr_kappa_mu,             pr_kappa_sd);
  lprior += beta_lpdf(grand_overdispersion_rho | pr_grand_overdispersion_rho_a,
                                                  pr_grand_overdispersion_rho_b);
}

model {
  target += lprior;

  if (prior_only == 0) {
    for (i in 1:N_obs) {
      real p = obs_prob(group[i], dose[i], CoP[i], T_thresh[i], stratum[i],
                       N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                       gamma_inf, gamma_fevginf, delta, pi_susc,
                       CoP_susc, CoP_imm, eta_lo, kappa,
                       T_ref, phi0_a, phi0_b);
      // Guard the beta_binomial's alpha/beta > 0 requirement (binomial tolerates
      // p in {0,1}; beta_binomial does not). At n=1 this is exactly binomial(1,p).
      real pc = fmin(fmax(p, 1e-12), 1.0 - 1e-12);
      y[i] ~ beta_binomial(n[i], pc * grand_concentration_k,
                           (1.0 - pc) * grand_concentration_k);
    }
    // Darton placebo temperature ladder -> phi0(T) (decoupled from dose-response;
    // NOT overdispersed -- a different sub-model than the Maryland replication rho targets).
    for (k in 1:N_ladder) {
      ladder_count[k] ~ binomial(ladder_N, phi0_fn(ladder_T[k], T_ref, phi0_a, phi0_b));
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
  // log_lik covers the WHOLE likelihood: N_obs dose-response rows followed by the
  // N_ladder phi0(T) threshold binomials, in that order (the layout is recorded in
  // each run's run_manifest.json as data.log_lik_layout).
  //
  // Why the tail is here: the ladder terms contribute to `target` in the model block
  // but used to appear in NEITHER log_lik NOR lprior, so target != lprior +
  // sum(log_lik). Two consequences, both silent: priorsense power-scales exactly
  // those two objects, so the reported likelihood-sensitivity of phi0_a/phi0_b
  // omitted the three binomials that IDENTIFY them; and loo dropped them from every
  // unit. They are genuine observations (16/10/8 of 20 TD+ Darton placebo subjects
  // crossing >=38/38.5/39 degC) and belong in both.
  vector[N_obs + N_ladder] log_lik;
  for (i in 1:N_obs) {
    p_pred[i]  = obs_prob(group[i], dose[i], CoP[i], T_thresh[i], stratum[i],
                          N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                          gamma_inf, gamma_fevginf, delta, pi_susc,
                          CoP_susc, CoP_imm, eta_lo, kappa,
                          T_ref, phi0_a, phi0_b);
    {
      real pc = fmin(fmax(p_pred[i], 1e-12), 1.0 - 1e-12);
      real a = pc * grand_concentration_k;
      real b = (1.0 - pc) * grand_concentration_k;
      y_rep[i]   = beta_binomial_rng(n[i], a, b);
      log_lik[i] = beta_binomial_lpmf(y[i] | n[i], a, b);
    }
  }
  for (k in 1:N_ladder)
    log_lik[N_obs + k] = binomial_lpmf(ladder_count[k] | ladder_N,
                                       phi0_fn(ladder_T[k], T_ref, phi0_a, phi0_b));

  // ---- Reference-dose derived quantities (naive Oxford, bicarb frame) -------
  real p_inf_1e3_naive = beta_poisson(1e3, N50_inf, alpha_inf, 1.0, gamma_inf);
  real p_inf_1e4_naive = beta_poisson(1e4, N50_inf, alpha_inf, 1.0, gamma_inf);
  real p_fev_1e3_naive = p_inf_1e3_naive
                         * beta_poisson(1e3, N50_fevginf, alpha_fevginf, 1.0, gamma_fevginf);
  real p_fev_1e4_naive = p_inf_1e4_naive
                         * beta_poisson(1e4, N50_fevginf, alpha_fevginf, 1.0, gamma_fevginf);

  // Maryland predicted fever curve at key doses (milk frame; phi(T,D) at Hornick 39.4)
  real p_fev_md_1e3 = obs_prob(3, 1e3, 1.0, 39.4, 0, N50_inf, N50_fevginf,
                               alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                               delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa,
                               T_ref, phi0_a, phi0_b);
  real p_fev_md_1e5 = obs_prob(3, 1e5, 1.0, 39.4, 0, N50_inf, N50_fevginf,
                               alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                               delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa,
                               T_ref, phi0_a, phi0_b);
  real p_fev_md_1e7 = obs_prob(3, 1e7, 1.0, 39.4, 0, N50_inf, N50_fevginf,
                               alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                               delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa,
                               T_ref, phi0_a, phi0_b);

  // Hornick Table 2 conditional prediction (phi(T,D) at Hornick 39.4)
  real p_cond_pred = obs_prob(5, 1e7, 1.0, 39.4, 0, N50_inf, N50_fevginf,
                              alpha_inf, alpha_fevginf, gamma_inf, gamma_fevginf,
                              delta, pi_susc, CoP_susc, CoP_imm, eta_lo, kappa,
                              T_ref, phi0_a, phi0_b);

  // ---- phi(T,D) diagnostics: low-dose asymptote phi0(T) + dose-lift at Hornick 39.4
  real phi0_38_3 = phi0_fn(38.3, T_ref, phi0_a, phi0_b);   // Levine/Gilman threshold
  real phi0_39_4 = phi0_fn(39.4, T_ref, phi0_a, phi0_b);   // Hornick threshold
  real phi_hornick_1e3 = phi_TD(39.4, 1e3 / delta, T_ref, phi0_a, phi0_b,
                                N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                                gamma_inf, gamma_fevginf);
  real phi_hornick_1e5 = phi_TD(39.4, 1e5 / delta, T_ref, phi0_a, phi0_b,
                                N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                                gamma_inf, gamma_fevginf);
  real phi_hornick_1e9 = phi_TD(39.4, 1e9 / delta, T_ref, phi0_a, phi0_b,
                                N50_inf, N50_fevginf, alpha_inf, alpha_fevginf,
                                gamma_inf, gamma_fevginf);

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
