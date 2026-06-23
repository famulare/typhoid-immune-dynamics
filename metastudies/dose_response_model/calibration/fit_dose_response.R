#' ---
#' title: "Typhoid dose-response: Tier 1 calibration driver (cmdstanr)"
#' ---
#'
#' Step 1 of the calibration ladder (see CALIBRATION_WORKFLOW.md): minimal Tier 1
#' fit of typhoid_dose_response.stan to the 25 tier1_active observations. Oxford
#' shedding / eta-correction and the study random effect are disabled (ox_inf rows
#' are tier1_active==0; sigma_study, eta_lo, kappa sample their priors only and are
#' INERT — do not interpret them).
#'
#' Workflow (Buffalo-style, R-native):
#'   - flat Stan data + prior hyperparameters-as-data via data_prep.R + priors.yaml
#'   - prior predictive (prior_only=1) and posterior fits
#'   - full diagnostics battery + PPC + priorsense via diagnostics.R
#'   - PPC fed by Stan generated-quantities p_pred (no R likelihood mirror)

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("^--file=", "", fa))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
calib_dir <- get_script_dir()
setwd(calib_dir)
source("priors.R"); source("data_prep.R"); source("diagnostics.R")

stan_file   <- "typhoid_dose_response.stan"
data_file   <- "dose_response_data.csv"
results_dir <- "results"

# Interpretable parameters (the inert sigma_study/eta_lo/kappa are excluded here).
interp_pars <- c("log10_N50_inf", "d_fev", "log10_N50_fevginf",
                 "alpha_inf", "alpha_fevginf", "gamma_inf", "gamma_fevginf",
                 "log10_delta", "pi_susc", "CoP_imm", "CoP_susc",
                 "N50_inf", "N50_fevginf", "delta")

priors <- load_priors()
stan_data <- build_stan_data(data_file, priors, tier_col = "tier1_active", prior_only = 0L)
obs <- attr(stan_data, "obs")
cat(sprintf("Tier 1 obs: %d (groups: %s)\n", stan_data$N_obs,
            paste(names(table(stan_data$group)), table(stan_data$group), sep = "x", collapse = ", ")))

mod <- cmdstan_model(stan_file)

# ---- Prior predictive --------------------------------------------------------
cat("\n=== PRIOR PREDICTIVE ===\n")
fit_prior <- mod$sample(
  data = modifyList(stan_data, list(prior_only = 1L)),
  chains = 2, parallel_chains = 2, iter_warmup = 500, iter_sampling = 500,
  seed = 1234, refresh = 0, show_messages = FALSE
)
diagnose_fit(fit_prior, file.path(results_dir, "tier1_prior"),
             pars = setdiff(interp_pars, c("N50_inf", "N50_fevginf", "delta")),
             obs = obs, priors = priors, model_name = "tier1_prior")

# ---- Posterior ---------------------------------------------------------------
cat("\n=== POSTERIOR ===\n")
t0 <- Sys.time()
fit <- mod$sample(
  data = stan_data,
  chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000,
  seed = 2024, refresh = 200, adapt_delta = 0.9
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("\n=== SAMPLER DIAGNOSTICS ===\n")
fit$cmdstan_diagnose()

diag <- diagnose_fit(fit, file.path(results_dir, "tier1"),
                     pars = interp_pars, obs = obs, priors = priors,
                     model_name = "tier1", elapsed_s = elapsed)

cat("\n=== INTERPRETABLE PARAMETERS ===\n")
print(diag$table, n = Inf)
cat("\nDone.\n")
