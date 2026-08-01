#' ---
#' title: "Typhoid dose-response: tier-parameterized calibration driver (cmdstanr)"
#' ---
#'
#' Fits typhoid_dose_response.stan for ONE spec from tier_specs.R. The spec -- not
#' this file -- declares the row set (`tier_col`), the Darton representation
#' (`individualize_darton`), the expected observation count, which parameters are
#' reported vs INERT, the model stage token, and the output directory. Default spec:
#' `t2-indiv-vax`. Counts and the naming rule for the discoverable rebuild ladder:
#' TIER_LADDER.md (generated from the registry).
#'
#' Apart from the preferred default below, this file names no tier as a literal and
#' hardcodes no parameter list. Its predecessor (fit_dose_response.R) did both, and its header claimed "the 25
#' tier1_active observations" while fitting 80 -- because build_stan_data() defaults
#' to individualize_darton = TRUE and nothing overrode it.
#'
#' `sigma_study` was DELETED 2026-07-31: declared with a prior, used in ZERO
#' likelihood terms at every tier, and there was no study index in the Stan data at
#' all. The cohort random effect it was for is LOCKED-WITHDRAWN
#' (cohort_random_effects_design.md), which already called for its deletion.
#' INERT until a group-2 (`ox_inf`) row is active, i.e. in every `t1-*` spec:
#' `eta_lo`, `kappa`. `beta_phi` is pinned to 1 in code. `phi_md` does NOT exist --
#' retired at C3 (d1880a8), replaced by `phi0_a`/`phi0_b`.
#'
#' Usage:
#'   Rscript fit_tier.R                      # default spec (t2-indiv-vax)
#'   Rscript fit_tier.R t1-grouped
#'   Rscript fit_tier.R --list               # the discoverable rebuild ladder
#'   Rscript fit_tier.R t1-indiv --quick     # short chains, for a plumbing check

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})
if (!exists("calib_dir")) source("utils.R")

DEFAULT_TIER <- "t2-indiv-vax"

#' Fit one tier: prior predictive + posterior + diagnostics + figure suite.
#'
#' @param key tier registry key (see tier_keys()).
#' @param allow_blocked fit a rung the registry declares blocked. Deliberate only --
#'   the blocked_reason explains what would silently be mis-specified.
fit_tier <- function(key = DEFAULT_TIER,
                     data_file = "dose_response_data.csv",
                     stan_file = "typhoid_dose_response.stan",
                     priors = load_priors(),
                     prior_predictive = TRUE, figures = TRUE,
                     chains = 4, warmup = 1000, sampling = 1000, adapt_delta = 0.9,
                     seed_post = 2024,
                     prior_chains = 2, prior_warmup = 500, prior_sampling = 500,
                     seed_prior = 1234,
                     out_root = "results", allow_blocked = FALSE) {
  spec <- tier_spec(key)

  # stanc only (no compilation) so the tier's parameter requirements and its stage
  # token are checked in ~1 s, before anything expensive happens.
  mod_chk <- tryCatch(cmdstan_model(stan_file, compile = FALSE), error = function(e) NULL)
  stan_data <- build_tier_data(key, data_file, priors, prior_only = 0L,
                               allow_blocked = allow_blocked, mod = mod_chk)
  obs  <- attr(stan_data, "obs")
  tier <- attr(stan_data, "tier")
  pars_post  <- tier_plot_pars(spec)
  pars_prior <- tier_report_pars(spec)

  cat(sprintf("\n=== %s (%s) ===\n", spec$key, spec$label))
  cat(sprintf("  %s | individualize_darton=%s | stage %s | status %s\n",
              spec$tier_col, spec$individualize_darton, spec$stage, tier$status))
  cat(sprintf("  N_obs %d (%d grouped, %d individual) | groups: %s\n",
              stan_data$N_obs, sum(obs$n > 1), sum(obs$n == 1),
              paste(sprintf("%s=%d", names(table(obs$likelihood_group)),
                            table(obs$likelihood_group)), collapse = ", ")))
  cat(sprintf("  reporting %d params; INERT (prior-only, do not interpret): %s\n",
              length(pars_post), paste(spec$inert_pars, collapse = ", ")))

  mod <- cmdstan_model(stan_file)

  if (isTRUE(prior_predictive)) {
    cat("\n--- prior predictive ---\n")
    sd_prior <- copy_stan_attrs(modifyList(stan_data, list(prior_only = 1L)), stan_data)
    fit_p <- mod$sample(data = sd_prior, chains = prior_chains,
                        parallel_chains = prior_chains, iter_warmup = prior_warmup,
                        iter_sampling = prior_sampling, seed = seed_prior,
                        refresh = 0, show_messages = FALSE)
    man_p <- run_manifest("prior_predictive", stan_data = sd_prior, priors = priors,
                          model_name = paste0(spec$key, "-prior"), data_csv = data_file,
                          stan_file = stan_file,
                          sampler = list(chains = prior_chains, iter_warmup = prior_warmup,
                                         iter_sampling = prior_sampling, seed = seed_prior,
                                         prior_only = 1L))
    diagnose_fit(fit_p, tier_out_dir(spec, out_root, prior = TRUE),
                 pars = pars_prior, obs = obs, priors = priors,
                 model_name = paste0(spec$key, "-prior"),
                 stan_data = sd_prior, manifest = man_p,
                 extra_plots = if (isTRUE(figures))
                   model_figures_hook(sd_prior, label = paste(spec$label, "(PRIOR)")))
  }

  cat("\n--- posterior ---\n")
  t0 <- Sys.time()
  fit <- mod$sample(data = stan_data, chains = chains, parallel_chains = chains,
                    iter_warmup = warmup, iter_sampling = sampling,
                    adapt_delta = adapt_delta, seed = seed_post, refresh = 200)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("\n--- sampler diagnostics ---\n"); fit$cmdstan_diagnose()

  man <- run_manifest("posterior", stan_data = stan_data, priors = priors,
                      model_name = spec$key, data_csv = data_file, stan_file = stan_file,
                      sampler = list(chains = chains, iter_warmup = warmup,
                                     iter_sampling = sampling, adapt_delta = adapt_delta,
                                     seed = seed_post, prior_only = 0L))
  diag <- diagnose_fit(fit, tier_out_dir(spec, out_root), pars = pars_post,
                       obs = obs, priors = priors, model_name = spec$key,
                       elapsed_s = elapsed, stan_data = stan_data, manifest = man,
                       extra_plots = if (isTRUE(figures))
                         model_figures_hook(stan_data, label = spec$label))

  cat("\n=== interpretable parameters ===\n"); print(diag$table, n = Inf)
  invisible(list(fit = fit, spec = spec, stan_data = stan_data, diag = diag))
}

#' The registry as a table -- the "declare a tier without running it" surface.
list_tiers <- function(mod = NULL) {
  rows <- do.call(rbind, lapply(discoverable_tier_keys(), function(k) {
    s <- tier_spec(k)
    data.frame(key = k, tier_col = s$tier_col,
               indiv = s$individualize_darton, N_obs = s$expect$N_obs,
               stage = s$stage, status = tier_status(s, mod),
               run_dir = tier_out_dir(s), stringsAsFactors = FALSE)
  }))
  cat("\n", paste(knitr_table(rows), collapse = "\n"), "\n", sep = "")
  cat("\nThis is the discoverable rebuild ladder. Retired/blocked registry entries are\n",
      "kept in tier_specs.R for history but omitted here.\n")
  invisible(rows)
}

if (sys.nframe() == 0) {
  setwd(calib_dir())
  source("priors.R"); source("data_prep.R"); source("diagnostics.R")
  source("figures.R")          # figure suite, attached via diagnose_fit(extra_plots=)
  source("tier_specs.R")

  args  <- commandArgs(trailingOnly = TRUE)
  flags <- grep("^--", args, value = TRUE)
  keys  <- setdiff(args, flags)
  quick <- "--quick" %in% flags

  if ("--list" %in% flags) {
    list_tiers(tryCatch(cmdstan_model("typhoid_dose_response.stan", compile = FALSE),
                        error = function(e) NULL))
  } else {
    key <- if (length(keys)) keys[1] else DEFAULT_TIER
    fit_tier(key, allow_blocked = "--allow-blocked" %in% flags,
             chains = if (quick) 2 else 4,
             warmup = if (quick) 300 else 1000,
             sampling = if (quick) 300 else 1000,
             prior_chains = if (quick) 1 else 2,
             prior_warmup = if (quick) 200 else 500,
             prior_sampling = if (quick) 200 else 500)
  }
}
