#' M4 parity gate (LOAD-BEARING) — see plan please-plan-a-b-smooth-ritchie.md.
#'
#' The recovery harness (simulate_recovery.R) shares obs_prob() between the data
#' simulator and the fitted model, so it is BLIND to an obs_prob() bug. This test
#' is the only guard that the unified obs_prob() faithfully reproduces the original
#' five-group likelihood. It compares the new Stan generated-quantities p_pred and
#' log_lik against an independent R transcription of the ORIGINAL per-group model
#' block, at several fixed parameter vectors, for every active Tier-1 row.
#'
#' Run:  Rscript test_obs_prob_parity.R   (exits non-zero on failure)

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})
source("priors.R"); source("data_prep.R")

# Tolerances at the floating-point / implementation-difference level (Stan pow &
# exact-lgamma binomial_lpmf vs R ^ & saddlepoint dbinom). A real logic error in
# obs_prob() would produce O(0.01-1) differences, far above these.
TOL_P  <- 1e-7   # probabilities (O(1) values)
TOL_LL <- 1e-3   # pointwise log-likelihood

# ---- Independent R reference (transcribed from the ORIGINAL model block) -------
bp <- function(D, N50, alpha, CoP, gamma) {
  scale <- (2^(1 / alpha) - 1) / N50
  1 - (1 + D * scale)^(-alpha / CoP^gamma)
}
md_mix <- function(D, N50, alpha, gamma, pi, CoPs, CoPi)
  pi * bp(D, N50, alpha, CoPs, gamma) + (1 - pi) * bp(D, N50, alpha, CoPi, gamma)

obs_prob_R <- function(row, p) {
  N50i <- 10^p$log10_N50_inf
  N50f <- 10^(p$log10_N50_inf + p$d_fev)         # reparam: log10_N50_fevginf = inf + d_fev
  delta <- 10^p$log10_delta
  g <- row$group
  if (g == 1L) {                                  # ox_fev
    D <- row$dose_cfu
    bp(D, N50i, p$alpha_inf, row$CoP, p$gamma_inf) *
      bp(D, N50f, p$alpha_fevginf, row$CoP, p$gamma_fevginf)
  } else if (g == 3L) {                           # md_fev
    D <- row$dose_cfu / delta; phi <- p$phi_md; st <- row$gilman_stratum  # phi now estimated scalar
    pf <- function(C) bp(D, N50i, p$alpha_inf, C, p$gamma_inf) *
                      bp(D, N50f, p$alpha_fevginf, C, p$gamma_fevginf)
    if (!is.na(st) && st == 1L) phi * pf(p$CoP_susc)
    else if (!is.na(st) && st == 2L) phi * pf(p$CoP_imm)
    else phi * (p$pi_susc * pf(p$CoP_susc) + (1 - p$pi_susc) * pf(p$CoP_imm))
  } else if (g == 4L) {                           # md_inf
    md_mix(row$dose_cfu / delta, N50i, p$alpha_inf, p$gamma_inf, p$pi_susc, p$CoP_susc, p$CoP_imm)
  } else if (g == 5L) {                           # hornick_cond
    D <- row$dose_cfu / delta; phi <- p$phi_md      # phi now estimated scalar
    pf <- function(C) bp(D, N50i, p$alpha_inf, C, p$gamma_inf) *
                      bp(D, N50f, p$alpha_fevginf, C, p$gamma_fevginf)
    p_inf <- md_mix(D, N50i, p$alpha_inf, p$gamma_inf, p$pi_susc, p$CoP_susc, p$CoP_imm)
    pc <- phi * (p$pi_susc * pf(p$CoP_susc) + (1 - p$pi_susc) * pf(p$CoP_imm)) / p_inf
    min(max(pc, 1e-12), 1 - 1e-12)
  } else stop("group 2 (ox_inf) not in Tier 1 parity set")
}

# ---- Build data + compile -----------------------------------------------------
priors <- load_priors()
sd <- build_stan_data("dose_response_data.csv", priors, tier_col = "tier1_active")
obs <- attr(sd, "obs")
mod <- cmdstan_model("typhoid_dose_response.stan")

# ---- Parameter vectors to test (constrained scale; must cover the model's params) ----
PARAM_NAMES <- c("log10_N50_inf","d_fev","alpha_inf","alpha_fevginf","gamma_inf",
                 "gamma_fevginf","log10_delta","pi_susc","CoP_imm","CoP_susc",
                 "phi_md","eta_lo","kappa","sigma_study")
vecs <- list(                                                  # phi_md added after CoP_susc
  c(2.5, 0.3, 0.30, 0.35, 0.60, 0.90, 3.5, 0.65, 3.0, 1.0, 0.90, 0.5, 1.0, 0.3),
  c(2.0, 0.0, 0.15, 0.50, 0.20, 1.50, 2.0, 0.40, 5.0, 1.1, 0.50, 0.4, 0.7, 0.1),  # d_fev=0 edge
  c(3.1, 1.2, 0.50, 0.20, 1.00, 0.30, 4.5, 0.80, 2.0, 0.9, 0.97, 0.6, 1.5, 0.5)
)
truth <- posterior::as_draws_matrix(do.call(rbind, lapply(vecs, function(v) setNames(v, PARAM_NAMES))))

# ---- New-model p_pred / log_lik via generate_quantities -----------------------
assert_fitted_params_match(mod, PARAM_NAMES)
gq <- mod$generate_quantities(fitted_params = truth, data = sd, seed = 1)
p_stan  <- posterior::as_draws_matrix(gq$draws("p_pred"))    # ndraws x N_obs
ll_stan <- posterior::as_draws_matrix(gq$draws("log_lik"))

# ---- Compare ------------------------------------------------------------------
max_dp <- 0; max_dl <- 0
for (di in seq_len(nrow(p_stan))) {
  p <- as.list(setNames(vecs[[di]], PARAM_NAMES))
  p_ref  <- vapply(seq_len(nrow(obs)), function(i) obs_prob_R(obs[i, ], p), numeric(1))
  ll_ref <- dbinom(obs$y, obs$n, p_ref, log = TRUE)
  max_dp <- max(max_dp, max(abs(p_stan[di, ] - p_ref)))
  max_dl <- max(max_dl, max(abs(ll_stan[di, ] - ll_ref)))
}
cat(sprintf("Parity over %d param vectors x %d rows:\n", nrow(p_stan), nrow(obs)))
cat(sprintf("  max |p_pred_stan - p_ref|   = %.3e\n", max_dp))
cat(sprintf("  max |log_lik_stan - ll_ref| = %.3e\n", max_dl))
if (max_dp < TOL_P && max_dl < TOL_LL) {
  cat("PARITY PASS\n")
} else {
  cat("PARITY FAIL — obs_prob() does not match the original per-group likelihood\n")
  quit(status = 1)
}
