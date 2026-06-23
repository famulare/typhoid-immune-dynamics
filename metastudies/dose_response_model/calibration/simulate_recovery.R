#' Synthetic-data recovery harness — the instrument for the divergence block.
#' See plan please-plan-a-b-smooth-ritchie.md (workstream B).
#'
#' Uses the SAME Stan code as a known-truth simulator (mod$generate_quantities at
#' fixed parameter values -> y_rep), then refits to check parameter recovery. This
#' separates "model/sampler broken" from "data don't inform the parameter".
#'
#' NOTE: recovery CANNOT catch an obs_prob() bug (the simulator and the fitted
#' model share it). The M4 parity gate (test_obs_prob_parity.R) is the guard for
#' that; run it first.
#'
#' Functions:
#'   make_truth()                  - 1-row draws_matrix of all 13 params (constrained)
#'   simulate_y()                  - draw synthetic y at a truth via generate_quantities
#'   recover_once()                - simulate -> refit -> diagnose vs truth
#'   recover_repeated()            - k truths from the prior -> SBC-lite coverage
#'   attribution_cliff_vs_reparam()- controlled cliff(T[0,]) vs reparam divergence test
#'
#' Quick demo:  Rscript simulate_recovery.R

suppressPackageStartupMessages({library(cmdstanr); library(posterior); library(dplyr)})
source("priors.R"); source("data_prep.R"); source("diagnostics.R")

PARAM_NAMES <- c("log10_N50_inf","d_fev","alpha_inf","alpha_fevginf","gamma_inf",
                 "gamma_fevginf","log10_delta","pi_susc","CoP_imm","CoP_susc",
                 "eta_lo","kappa","sigma_study")

# A realistic truth (near the Tier-1 posterior) for point recovery.
TRUTH_REALISTIC <- c(log10_N50_inf = 2.0, d_fev = 0.3, alpha_inf = 0.16,
                     alpha_fevginf = 0.8, gamma_inf = 0.6, gamma_fevginf = 0.7,
                     log10_delta = 1.6, pi_susc = 0.7, CoP_imm = 1.8, CoP_susc = 0.8,
                     eta_lo = 0.5, kappa = 1.0, sigma_study = 0.3)

#' Build a draws_matrix of true parameter values (constrained scale).
#' @param values named numeric over PARAM_NAMES (defaults to TRUTH_REALISTIC).
make_truth <- function(values = TRUTH_REALISTIC) {
  stopifnot(setequal(names(values), PARAM_NAMES))
  posterior::as_draws_matrix(t(as.matrix(values[PARAM_NAMES])))
}

#' Draw `reps` synthetic y vectors at the given truth via generate_quantities.
#' @return integer matrix reps x N_obs (one synthetic dataset per row).
simulate_y <- function(mod, truth, stan_data, seed = 99) {
  gq <- mod$generate_quantities(fitted_params = truth, data = stan_data, seed = seed)
  yr <- posterior::as_draws_matrix(gq$draws("y_rep"))
  matrix(as.integer(yr), nrow = nrow(yr), dimnames = list(NULL, NULL))
}

#' Simulate one dataset at `truth` and refit; diagnose recovery.
recover_once <- function(mod, stan_data, obs, priors, truth_values = TRUTH_REALISTIC,
                         out_dir = "results/recovery/point",
                         chains = 4, warmup = 1000, sampling = 1000,
                         adapt_delta = 0.9, sim_seed = 99, fit_seed = 7) {
  truth <- make_truth(truth_values)
  y_sim <- simulate_y(mod, truth, stan_data, seed = sim_seed)[1, ]
  rec_data <- modifyList(stan_data, list(y = as.integer(y_sim)))
  fit <- mod$sample(data = rec_data, chains = chains, parallel_chains = chains,
                    iter_warmup = warmup, iter_sampling = sampling,
                    adapt_delta = adapt_delta, seed = fit_seed, refresh = 0,
                    show_messages = FALSE)
  # derived truths so the table can mark recovery on transformed params too
  tv <- truth_values
  tv["log10_N50_fevginf"] <- tv["log10_N50_inf"] + tv["d_fev"]
  tv["N50_inf"] <- 10^tv["log10_N50_inf"]
  tv["N50_fevginf"] <- 10^tv["log10_N50_fevginf"]
  tv["delta"] <- 10^tv["log10_delta"]
  pars <- c("log10_N50_inf","d_fev","log10_N50_fevginf","alpha_inf","alpha_fevginf",
            "gamma_inf","gamma_fevginf","log10_delta","pi_susc","CoP_imm","CoP_susc")
  diagnose_fit(fit, out_dir, pars = pars, true_params = tv, obs = obs,
               priors = priors, model_name = "recovery_point")
}

#' SBC-lite: k truths drawn from the prior, one synthetic dataset each, refit each,
#' aggregate 90% CI coverage per parameter. Tests calibration on the real design.
recover_repeated <- function(mod, stan_data, priors, k = 20,
                             out_csv = "results/recovery/coverage.csv",
                             chains = 4, warmup = 600, sampling = 600,
                             adapt_delta = 0.9, seed0 = 1000) {
  pars <- c("log10_N50_inf","d_fev","alpha_inf","alpha_fevginf","gamma_inf",
            "gamma_fevginf","log10_delta","pi_susc","CoP_imm","CoP_susc")
  rng <- set.seed(seed0)
  rows <- list()
  for (rep in seq_len(k)) {
    tv <- vapply(PARAM_NAMES, function(p) sample_prior(priors, p, 1L), numeric(1))
    truth <- make_truth(tv)
    y_sim <- simulate_y(mod, truth, stan_data, seed = seed0 + rep)[1, ]
    rec_data <- modifyList(stan_data, list(y = as.integer(y_sim)))
    fit <- mod$sample(data = rec_data, chains = chains, parallel_chains = chains,
                      iter_warmup = warmup, iter_sampling = sampling,
                      adapt_delta = adapt_delta, seed = seed0 + rep,
                      refresh = 0, show_messages = FALSE)
    q <- fit$summary(pars, ~posterior::quantile2(.x, c(0.05, 0.95)))
    div <- sum(fit$sampler_diagnostics()[, , "divergent__"])
    for (p in pars) {
      qq <- q[q$variable == p, ]
      rows[[length(rows) + 1]] <- data.frame(
        rep = rep, param = p, true = tv[p],
        q5 = qq$q5, q95 = qq$q95,
        covered = tv[p] >= qq$q5 & tv[p] <= qq$q95,
        divergences = div)
    }
    cat(sprintf("  rep %2d/%d: %d divergences\n", rep, k, div))
  }
  df <- dplyr::bind_rows(rows)
  cov <- df %>% group_by(param) %>%
    summarise(coverage90 = mean(covered), n = n(), .groups = "drop")
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, out_csv)
  cat("\n90% CI coverage (target ~0.90):\n"); print(as.data.frame(cov))
  cat(sprintf("Total divergences across %d fits: %d\n", k, sum(unique(df[c("rep","divergences")])$divergences)))
  invisible(list(detail = df, coverage = cov))
}

#' Controlled attribution: same real data + priors + likelihood, ONLY the N50
#' parameterization differs (T[0,] cliff vs <lower=0> reparam). Reconstructs the
#' cliff model from git HEAD; runs both as short fits and reports divergence rate.
attribution_cliff_vs_reparam <- function(data_csv = "dose_response_data.csv",
                                         priors = load_priors(),
                                         chains = 4, warmup = 1000, sampling = 1000,
                                         adapt_delta = 0.9, seed = 2024) {
  base <- "results/_attribution"; dir.create(base, showWarnings = FALSE, recursive = TRUE)

  # --- reparam (current) model on flat data ---
  mod_new <- cmdstan_model("typhoid_dose_response.stan")
  sd_new <- build_stan_data(data_csv, priors, tier_col = "tier1_active")
  fit_new <- mod_new$sample(data = sd_new, chains = chains, parallel_chains = chains,
                            iter_warmup = warmup, iter_sampling = sampling,
                            adapt_delta = adapt_delta, seed = seed, refresh = 0,
                            show_messages = FALSE)
  div_new <- sum(fit_new$sampler_diagnostics()[, , "divergent__"])
  n_new <- chains * sampling

  # --- cliff (git HEAD) model on old per-group data ---
  cliff_stan <- file.path(base, "typhoid_dose_response_cliff.stan")
  ok <- tryCatch({
    src <- system2("git", c("show", "HEAD:metastudies/dose_response_model/calibration/typhoid_dose_response.stan"),
                   stdout = TRUE, stderr = FALSE)
    writeLines(src, cliff_stan); TRUE
  }, error = function(e) FALSE)
  if (!ok || !grepl("T\\[0, \\]", paste(readLines(cliff_stan), collapse = "\n"))) {
    cat("[skip] could not reconstruct cliff model from git HEAD\n")
    return(invisible(list(reparam = list(div = div_new, n = n_new))))
  }
  mod_cliff <- cmdstan_model(cliff_stan, dir = base)
  sd_old <- build_old_pergroup_data(data_csv)
  fit_cliff <- mod_cliff$sample(data = sd_old, chains = chains, parallel_chains = chains,
                                iter_warmup = warmup, iter_sampling = sampling,
                                adapt_delta = adapt_delta, seed = seed, refresh = 0,
                                show_messages = FALSE)
  div_cliff <- sum(fit_cliff$sampler_diagnostics()[, , "divergent__"])
  n_cliff <- chains * sampling

  cat("\n=== DIVERGENCE ATTRIBUTION (same data/priors/likelihood) ===\n")
  cat(sprintf("  cliff  (T[0,] on N50 difference): %d / %d (%.1f%%)\n",
              div_cliff, n_cliff, 100 * div_cliff / n_cliff))
  cat(sprintf("  reparam (d_fev <lower=0> offset): %d / %d (%.1f%%)\n",
              div_new, n_new, 100 * div_new / n_new))
  cat("  => divergences attributable to the N50 geometry alone.\n")
  invisible(list(cliff = list(div = div_cliff, n = n_cliff),
                 reparam = list(div = div_new, n = n_new)))
}

#' Old per-group Stan data for the cliff (git HEAD) model.
build_old_pergroup_data <- function(data_csv) {
  d <- readr::read_csv(data_csv, show_col_types = FALSE)
  d1 <- d %>% filter(tier1_active == 1)
  g <- function(grp) d1 %>% filter(likelihood_group == grp) %>% arrange(obs_id)
  ox_fev <- g("ox_fev"); md_fev <- g("md_fev"); md_inf <- g("md_inf"); hc <- g("hornick_cond")
  stopifnot(nrow(hc) == 1L)
  list(
    N_ox_fev = nrow(ox_fev), n_ox_fev = as.integer(ox_fev$n), y_ox_fev = as.integer(ox_fev$y),
    dose_ox_fev = ox_fev$dose_cfu, CoP_ox_fev = ox_fev$CoP,
    N_ox_inf = 0L, n_ox_inf = integer(0), y_ox_inf = integer(0),
    dose_ox_inf = numeric(0), CoP_ox_inf = numeric(0),
    N_md_fev = nrow(md_fev), n_md_fev = as.integer(md_fev$n), y_md_fev = as.integer(md_fev$y),
    dose_md_fev = md_fev$dose_cfu, phi_md_fev = md_fev$phi,
    gilman_stratum = as.integer(md_fev$gilman_stratum),
    N_md_inf = nrow(md_inf), n_md_inf = as.integer(md_inf$n), y_md_inf = as.integer(md_inf$y),
    dose_md_inf = md_inf$dose_cfu,
    n_hornick_cond = as.integer(hc$n), y_hornick_cond = as.integer(hc$y),
    dose_hornick_7 = hc$dose_cfu, phi_hornick = hc$phi,
    prior_only = 0L
  )
}

# ---- demo main (small/fast; scale up k and iters for real validation) ---------
if (sys.nframe() == 0) {
  priors <- load_priors()
  mod <- cmdstan_model("typhoid_dose_response.stan")
  stan_data <- build_stan_data("dose_response_data.csv", priors, tier_col = "tier1_active")
  obs <- attr(stan_data, "obs")

  cat("\n[1] smoke test: generate_quantities -> y_rep length\n")
  ys <- simulate_y(mod, make_truth(), stan_data, seed = 1)
  cat(sprintf("    y_rep dim = %d x %d (expect 1 x %d)\n", nrow(ys), ncol(ys), stan_data$N_obs))

  cat("\n[2] point recovery at a realistic truth\n")
  recover_once(mod, stan_data, obs, priors)

  cat("\n[3] cliff-vs-reparam divergence attribution\n")
  attribution_cliff_vs_reparam(priors = priors)

  cat("\n[4] SBC-lite coverage (k=6 demo; use k>=40 for real SBC)\n")
  recover_repeated(mod, stan_data, priors, k = 6)
}
