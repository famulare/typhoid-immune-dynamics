#' Simulated-data parameter-recovery example (single prior draw), for any tier.
#'
#' The standard synthetic-data check: draw ONE parameter vector from the prior,
#' simulate a synthetic dataset matching a tier's real covariate table (same
#' group/dose/CoP/stratum/n; only y resampled from the drawn truth), refit, and
#' confirm the parameters are recovered within statistical expectations. One point
#' sample. Report goes to results/recovery/<tier>/prior_draw/.
#'
#' The tier comes from tier_specs.R, so which rows/representation are being recovered
#' is recorded rather than implied -- and report_pars/inert_pars are derived from that
#' spec instead of being a hardcoded Tier-1 pair.
#'
#' Run from the calibration directory (relative paths, matching the harness):
#'   Rscript recover_from_prior_example.R            # t1-indiv
#'   Rscript recover_from_prior_example.R t1-grouped
#'
#' Run test_obs_prob_parity.R first — recovery shares obs_prob() with the fit and
#' cannot catch a bug there.

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})
source("simulate_recovery.R")   # also sources priors.R / data_prep.R / diagnostics.R
if (!exists("tier_spec")) source("tier_specs.R")

args <- commandArgs(trailingOnly = TRUE)
key  <- if (length(args)) args[1] else "t1-indiv"

priors    <- load_priors()
mod       <- cmdstan_model("typhoid_dose_response.stan")
stan_data <- build_tier_data(key, "dose_response_data.csv", priors, mod = mod)
spec      <- tier_spec(key)

recover_from_prior(mod, stan_data, attr(stan_data, "obs"), priors, seed = 2026,
                   report_pars = tier_report_pars(spec),
                   inert_pars  = spec$inert_pars)
