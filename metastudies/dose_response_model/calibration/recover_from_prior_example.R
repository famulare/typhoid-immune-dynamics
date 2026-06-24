#' Tier-1 simulated-data parameter-recovery example (single prior draw).
#'
#' The standard synthetic-data check: draw ONE parameter vector from the prior,
#' simulate a synthetic dataset matching the real Tier-1 covariate table (same
#' group/dose/CoP/stratum/n; only y resampled from the drawn truth), refit, and
#' confirm the parameters are recovered within statistical expectations. One
#' point sample. Writes the recovery report to results/recovery/tier1/.
#'
#' Run from the calibration directory (relative paths, matching the harness):
#'   Rscript recover_from_prior_example.R
#'
#' Run test_obs_prob_parity.R first — recovery shares obs_prob() with the fit and
#' cannot catch a bug there.

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})
source("simulate_recovery.R")   # also sources priors.R / data_prep.R / diagnostics.R

priors    <- load_priors()
mod       <- cmdstan_model("typhoid_dose_response.stan")
stan_data <- build_stan_data("dose_response_data.csv", priors, tier_col = "tier1_active")

recover_from_prior(mod, stan_data, attr(stan_data, "obs"), priors, seed = 2026)
