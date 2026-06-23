#' ---
#' title: "Typhoid dose-response: Tier 1 calibration driver (cmdstanr)"
#' ---
#'
#' Step 1 of the calibration ladder (see CALIBRATION_WORKFLOW.md):
#' Minimal Tier 1 fit of typhoid_dose_response.stan to the 25 tier1_active
#' observations in dose_response_data.csv. Oxford shedding / eta-correction and
#' the study random effect are disabled here (N_ox_inf = 0; sigma_study, eta_lo,
#' kappa sample their priors only and are INERT — do not interpret them).

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
})

# --- Paths (work whether sourced or run via Rscript) ---------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
calib_dir   <- get_script_dir()
stan_file   <- file.path(calib_dir, "typhoid_dose_response.stan")
data_file   <- file.path(calib_dir, "dose_response_data.csv")
results_dir <- file.path(calib_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- Assemble the Stan data list (Tier 1) --------------------------------------
d    <- read_csv(data_file, show_col_types = FALSE)
dat1 <- d %>% filter(tier1_active == 1)
grp  <- function(g) dat1 %>% filter(likelihood_group == g) %>% arrange(obs_id)

ox_fev <- grp("ox_fev")
md_fev <- grp("md_fev")
md_inf <- grp("md_inf")
hc     <- grp("hornick_cond")
stopifnot(nrow(hc) == 1L)
cat(sprintf("Tier 1 obs: ox_fev=%d, md_fev=%d, md_inf=%d, hornick_cond=%d (total=%d)\n",
            nrow(ox_fev), nrow(md_fev), nrow(md_inf), nrow(hc),
            nrow(ox_fev) + nrow(md_fev) + nrow(md_inf) + nrow(hc)))

stan_data <- list(
  N_ox_fev    = nrow(ox_fev),
  n_ox_fev    = as.integer(ox_fev$n),
  y_ox_fev    = as.integer(ox_fev$y),
  dose_ox_fev = ox_fev$dose_cfu,
  CoP_ox_fev  = ox_fev$CoP,

  # Oxford shedding disabled in Tier 1 (zero-length -> empty Stan loop)
  N_ox_inf    = 0L,
  n_ox_inf    = integer(0),
  y_ox_inf    = integer(0),
  dose_ox_inf = numeric(0),
  CoP_ox_inf  = numeric(0),

  N_md_fev       = nrow(md_fev),
  n_md_fev       = as.integer(md_fev$n),
  y_md_fev       = as.integer(md_fev$y),
  dose_md_fev    = md_fev$dose_cfu,
  phi_md_fev     = md_fev$phi,
  gilman_stratum = as.integer(md_fev$gilman_stratum),

  N_md_inf    = nrow(md_inf),
  n_md_inf    = as.integer(md_inf$n),
  y_md_inf    = as.integer(md_inf$y),
  dose_md_inf = md_inf$dose_cfu,

  n_hornick_cond = as.integer(hc$n),     # 28 infected (conditional denominator)
  y_hornick_cond = as.integer(hc$y),     # 16 febrile among infected
  dose_hornick_7 = hc$dose_cfu,          # 1e7
  phi_hornick    = hc$phi,               # 0.25

  prior_only = 0L
)
stopifnot(!anyNA(unlist(stan_data[c("CoP_ox_fev","phi_md_fev","gilman_stratum",
                                    "dose_md_inf","dose_md_fev")])))

# --- Compile -------------------------------------------------------------------
mod <- cmdstan_model(stan_file)

# --- Prior predictive check ----------------------------------------------------
gq_vars <- c("p_inf_1e3_naive","p_inf_1e4_naive","p_fev_1e3_naive","p_fev_1e4_naive",
             "p_fev_md_1e3","p_fev_md_1e5","p_fev_md_1e7","p_cond_pred")
fit_prior <- mod$sample(
  data = modifyList(stan_data, list(prior_only = 1L)),
  chains = 2, parallel_chains = 2, iter_warmup = 500, iter_sampling = 500,
  seed = 1234, refresh = 0
)
cat("\n=== PRIOR PREDICTIVE (generated-quantity attack rates) ===\n")
print(fit_prior$summary(variables = gq_vars,
                        ~quantile(.x, c(0.05, 0.5, 0.95), na.rm = TRUE)))

# --- Posterior fit -------------------------------------------------------------
fit <- mod$sample(
  data = stan_data,
  chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000,
  seed = 2024, refresh = 200, adapt_delta = 0.9
)

cat("\n=== SAMPLER DIAGNOSTICS ===\n")
fit$cmdstan_diagnose()

interp_pars <- c("log10_N50_inf","log10_N50_fevginf","alpha_inf","alpha_fevginf",
                 "gamma_inf","gamma_fevginf","log10_delta","pi_susc",
                 "CoP_imm","CoP_susc","N50_inf","N50_fevginf","delta")
inert_pars  <- c("sigma_study","eta_lo","kappa")  # Step 1: prior-only, do not interpret

summ <- fit$summary()
cat("\n=== INTERPRETABLE PARAMETERS ===\n")
print(fit$summary(variables = interp_pars))
cat("\n=== INERT PARAMETERS (Step 1: prior-only; R-hat shown only as a sampler check) ===\n")
print(fit$summary(variables = inert_pars))
write_csv(summ, file.path(results_dir, "tier1_summary.csv"))

# --- Posterior predictive: fitted attack rate per observation ------------------
# Mirror of the Stan likelihood, evaluated on posterior draws.
dr <- as_draws_df(fit$draws(c("N50_inf","N50_fevginf","alpha_inf","alpha_fevginf",
                              "gamma_inf","gamma_fevginf","delta",
                              "pi_susc","CoP_imm","CoP_susc")))
bp <- function(D_eff, N50, alpha, CoP, gamma) {
  scale <- (2^(1 / alpha) - 1) / N50
  1 - (1 + D_eff * scale)^(-alpha / CoP^gamma)
}
md_mix <- function(D_eff, N50, alpha, gamma, pi, CoPs, CoPi) {
  pi * bp(D_eff, N50, alpha, CoPs, gamma) + (1 - pi) * bp(D_eff, N50, alpha, CoPi, gamma)
}

fitted_p <- function(row) {
  with(dr, {
    if (row$likelihood_group == "ox_fev") {
      D <- row$dose_cfu; CoP <- row$CoP
      bp(D, N50_inf, alpha_inf, CoP, gamma_inf) *
        bp(D, N50_fevginf, alpha_fevginf, CoP, gamma_fevginf)
    } else if (row$likelihood_group == "md_inf") {
      D <- row$dose_cfu / delta
      md_mix(D, N50_inf, alpha_inf, gamma_inf, pi_susc, CoP_susc, CoP_imm)
    } else if (row$likelihood_group == "md_fev") {
      D <- row$dose_cfu / delta; phi <- row$phi; st <- row$gilman_stratum
      pf <- function(CoP) bp(D, N50_inf, alpha_inf, CoP, gamma_inf) *
                          bp(D, N50_fevginf, alpha_fevginf, CoP, gamma_fevginf)
      if (st == 1) phi * pf(CoP_susc)
      else if (st == 2) phi * pf(CoP_imm)
      else phi * (pi_susc * pf(CoP_susc) + (1 - pi_susc) * pf(CoP_imm))
    } else if (row$likelihood_group == "hornick_cond") {
      D <- row$dose_cfu / delta; phi <- row$phi
      pf <- function(CoP) bp(D, N50_inf, alpha_inf, CoP, gamma_inf) *
                          bp(D, N50_fevginf, alpha_fevginf, CoP, gamma_fevginf)
      p_inf <- md_mix(D, N50_inf, alpha_inf, gamma_inf, pi_susc, CoP_susc, CoP_imm)
      phi * (pi_susc * pf(CoP_susc) + (1 - pi_susc) * pf(CoP_imm)) / p_inf
    } else rep(NA_real_, length(N50_inf))
  })
}

ppc <- dat1 %>%
  rowwise() %>%
  mutate(
    p_draws = list(fitted_p(cur_data())),
    fit_med = median(p_draws),
    fit_lo  = quantile(p_draws, 0.05),
    fit_hi  = quantile(p_draws, 0.95),
    obs_rate = y / n
  ) %>%
  ungroup() %>%
  select(obs_id, study, likelihood_group, dose_cfu, n, y, obs_rate,
         fit_med, fit_lo, fit_hi)
write_csv(ppc, file.path(results_dir, "tier1_ppc.csv"))

p <- ggplot(ppc, aes(obs_rate, fit_med, color = likelihood_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_linerange(aes(ymin = fit_lo, ymax = fit_hi), alpha = 0.5) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = obs_id), size = 2.5, max.overlaps = 20,
                           show.legend = FALSE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Observed attack rate (y/n)", y = "Posterior fitted probability (median, 90% CI)",
       title = "Tier 1 posterior predictive check", color = "group") +
  theme_minimal()
ggsave(file.path(results_dir, "tier1_ppc.png"), p, width = 8, height = 6, dpi = 150)

# --- Persist -------------------------------------------------------------------
fit$save_object(file.path(results_dir, "fit_tier1.rds"))
cat("\nDone. Outputs in:", results_dir, "\n")
