#' Light cross-run comparison harness (workstream C).
#' See plan please-plan-a-b-smooth-ritchie.md.
#'
#' Plain R (no Snakemake/targets): a named list of scenario specs + a runner +
#' an aggregator. Each fit drops results.json via diagnose_fit(); summarize reads
#' them into a comparison table + forest plot. Because priors are data (workstream
#' D), prior-override sensitivities (e.g. the delta-prior grid) are cheap refits
#' with NO recompile.
#'
#' Scenario spec fields: label, tier_col, drop_obs, keep_obs, prior_overrides.
#' Structural-Stan variants (share alpha/gamma, single CoP_md, drop phi) need new
#' .stan files; they plug into this same runner once written (deferred).
#'
#' Quick demo:  Rscript run_scenarios.R

suppressPackageStartupMessages({library(cmdstanr); library(posterior); library(dplyr)
                                library(ggplot2); library(tidyr)})
source("priors.R"); source("data_prep.R"); source("diagnostics.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

KEY_PARS <- c("log10_N50_inf","d_fev","gamma_inf","gamma_fevginf","log10_delta",
              "pi_susc","CoP_imm","alpha_inf","alpha_fevginf","CoP_susc")

# Example scenario set (all no-new-Stan-code: row filters + prior overrides).
SCENARIOS <- list(
  list(label = "tier1_base"),
  list(label = "hornick3_excluded", drop_obs = "H-F-3"),                 # MC4 sensitivity
  list(label = "gilrest_excluded",  drop_obs = "Gil-F-rest"),           # derived-by-subtraction row
  list(label = "delta_prior_lo",    prior_overrides = list(log10_delta = list(mu = 3.0))),
  list(label = "delta_prior_hi",    prior_overrides = list(log10_delta = list(mu = 4.0)))
)

#' Run one scenario -> results/scenarios/<label>/ (diagnostics + results.json + loo.json).
run_scenario <- function(spec, mod, data_csv = "dose_response_data.csv",
                         priors0 = load_priors(),
                         chains = 4, warmup = 800, sampling = 800,
                         adapt_delta = 0.9, seed = 2024) {
  priors <- apply_prior_overrides(priors0, spec$prior_overrides %||% list())
  stan_data <- build_stan_data(data_csv, priors,
                               tier_col = spec$tier_col %||% "tier1_active",
                               drop_obs = spec$drop_obs %||% character(),
                               keep_obs = spec$keep_obs)
  obs <- attr(stan_data, "obs")
  out_dir <- file.path("results", "scenarios", spec$label)

  t0 <- Sys.time()
  fit <- mod$sample(data = stan_data, chains = chains, parallel_chains = chains,
                    iter_warmup = warmup, iter_sampling = sampling,
                    adapt_delta = adapt_delta, seed = seed, refresh = 0,
                    show_messages = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  diagnose_fit(fit, out_dir,
               pars = c("log10_N50_inf","d_fev","log10_N50_fevginf","alpha_inf",
                        "alpha_fevginf","gamma_inf","gamma_fevginf","log10_delta",
                        "pi_susc","CoP_imm","CoP_susc"),
               obs = obs, priors = priors, model_name = spec$label, elapsed_s = elapsed)

  lj <- compute_loo_units(fit, obs)
  if (!is.null(lj)) jsonlite::write_json(lj, file.path(out_dir, "loo.json"),
                                         auto_unbox = TRUE, pretty = TRUE, digits = 6)
  invisible(out_dir)
}

#' PSIS-LOO on correctly-grouped observation UNITS: the Hornick infection marginal
#' (H-I-7) and fever conditional (H-FgI-7) are one table factorized -> one unit.
compute_loo_units <- function(fit, obs) {
  if (!requireNamespace("loo", quietly = TRUE)) return(NULL)
  if (!"log_lik" %in% fit$metadata()$stan_variables) return(NULL)
  ll_arr <- fit$draws("log_lik")                          # iter x chain x N_obs (col order = obs rows)
  ll <- posterior::as_draws_matrix(ll_arr)                 # (iter*chain) x N_obs
  unit <- obs$obs_id
  unit[unit %in% c("H-I-7", "H-FgI-7")] <- "hornick_1e7"   # combine the joint factorization
  u <- unique(unit)
  ll_u <- vapply(u, function(k) rowSums(ll[, unit == k, drop = FALSE]), numeric(nrow(ll)))
  r_eff <- tryCatch(loo::relative_eff(exp(ll_u),
                      chain_id = rep(seq_len(posterior::nchains(ll_arr)),
                                     each = posterior::niterations(ll_arr))),
                    error = function(e) NULL)
  lo <- tryCatch(loo::loo(ll_u, r_eff = r_eff), error = function(e) {
    message("  [skip] loo: ", conditionMessage(e)); NULL })
  if (is.null(lo)) return(NULL)
  list(n_units = length(u),
       elpd_loo = lo$estimates["elpd_loo", "Estimate"],
       elpd_loo_se = lo$estimates["elpd_loo", "SE"],
       p_loo = lo$estimates["p_loo", "Estimate"],
       n_pareto_k_gt_0.7 = sum(loo::pareto_k_values(lo) > 0.7))
}

#' Aggregate results.json (+ loo.json) across scenarios into a comparison table + forest plot.
summarize_scenarios <- function(labels, out_dir = "results/summaries") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  recs <- lapply(labels, function(lab) {
    rj <- file.path("results", "scenarios", lab, "results.json")
    if (!file.exists(rj)) return(NULL)
    r <- jsonlite::read_json(rj, simplifyVector = FALSE)
    lj <- file.path("results", "scenarios", lab, "loo.json")
    loo_r <- if (file.exists(lj)) jsonlite::read_json(lj, simplifyVector = FALSE) else NULL
    list(label = lab, r = r, loo = loo_r)
  })
  recs <- Filter(Negate(is.null), recs)

  # long table of all params
  rows <- list()
  for (e in recs) for (p in names(e$r$parameters)) {
    pd <- e$r$parameters[[p]]
    rows[[length(rows) + 1]] <- data.frame(
      scenario = e$label, param = p,
      mean = pd$mean %||% NA, lo90 = pd$ci90[[1]] %||% NA, hi90 = pd$ci90[[2]] %||% NA,
      rhat = pd$rhat %||% NA, ess_bulk = pd$ess_bulk %||% NA)
  }
  long <- dplyr::bind_rows(rows)

  # scenario-level summary (divergences, min ebfmi, elpd_loo)
  scen <- dplyr::bind_rows(lapply(recs, function(e) data.frame(
    scenario = e$label,
    n_divergent = e$r$sampler$n_divergent %||% NA,
    max_rhat = max(vapply(e$r$parameters, function(x) x$rhat %||% NA, numeric(1)), na.rm = TRUE),
    min_ess_bulk = min(vapply(e$r$parameters, function(x) x$ess_bulk %||% NA, numeric(1)), na.rm = TRUE),
    elpd_loo = if (!is.null(e$loo)) e$loo$elpd_loo else NA,
    elpd_loo_se = if (!is.null(e$loo)) e$loo$elpd_loo_se else NA,
    pareto_k_gt_0.7 = if (!is.null(e$loo)) e$loo$n_pareto_k_gt_0.7 else NA)))

  readr::write_csv(long, file.path(out_dir, "comparison.csv"))
  readr::write_csv(scen, file.path(out_dir, "comparison_scenarios.csv"))

  # comparison.md
  L <- c("# Scenario comparison", "",
         "## Sampler health + LOO (loo across scenarios with identical units only)", "",
         knitr_table(scen), "",
         "## Key parameters (posterior mean [90% CI])", "")
  wide <- long %>% filter(param %in% KEY_PARS) %>%
    mutate(cell = sprintf("%.3g [%.3g, %.3g]", mean, lo90, hi90)) %>%
    select(scenario, param, cell) %>% tidyr::pivot_wider(names_from = param, values_from = cell)
  L <- c(L, knitr_table(as.data.frame(wide)))
  writeLines(L, file.path(out_dir, "comparison.md"))

  # forest plot of key params across scenarios
  fp <- long %>% filter(param %in% KEY_PARS)
  p <- ggplot(fp, aes(mean, scenario)) +
    geom_pointrange(aes(xmin = lo90, xmax = hi90)) +
    facet_wrap(~param, scales = "free_x") +
    labs(title = "Posterior mean +/- 90% CI across scenarios", x = NULL, y = NULL)
  ggplot2::ggsave(file.path(out_dir, "forest.png"), p, width = 12, height = 7, dpi = 150, bg = "white")

  cat(sprintf("\nScenario comparison saved to: %s/\n", normalizePath(out_dir)))
  cat("  comparison.md / comparison.csv / comparison_scenarios.csv / forest.png\n")
  print(scen)
  invisible(list(long = long, scen = scen))
}

# minimal markdown table (avoid a knitr dependency)
knitr_table <- function(df) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1, function(r) paste0("| ", paste(format(r, trim = TRUE), collapse = " | "), " |"))
  c(hdr, sep, body)
}

# ---- demo main ----------------------------------------------------------------
if (sys.nframe() == 0) {
  mod <- cmdstan_model("typhoid_dose_response.stan")
  priors0 <- load_priors()
  for (spec in SCENARIOS) {
    cat(sprintf("\n=== scenario: %s ===\n", spec$label))
    run_scenario(spec, mod, priors0 = priors0)
  }
  summarize_scenarios(vapply(SCENARIOS, function(s) s$label, character(1)))
}
