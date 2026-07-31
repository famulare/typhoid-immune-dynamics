#' Light cross-run comparison harness (workstream C).
#' See plan please-plan-a-b-smooth-ritchie.md.
#'
#' Plain R (no Snakemake/targets): a named list of scenario specs + a runner +
#' an aggregator. Each fit drops results.json via diagnose_fit(); summarize reads
#' them into a comparison table + forest plot. Because priors are data (workstream
#' D), prior-override sensitivities (e.g. the delta-prior grid) are cheap refits
#' with NO recompile.
#'
#' Scenario spec fields: label, tier (a tier_specs.R key), drop_obs, keep_obs,
#' prior_overrides. Output goes to results/scenarios/<tier>/<label>/.
#' Structural-Stan variants (share alpha/gamma, single CoP_md, drop phi) need new
#' .stan files; they plug into this same runner once written (deferred).
#'
#' Quick demo:  Rscript run_scenarios.R

suppressPackageStartupMessages({library(cmdstanr); library(posterior); library(dplyr)
                                library(ggplot2); library(tidyr)})
source("priors.R"); source("data_prep.R"); source("diagnostics.R")
source("figures.R")   # bespoke model figure suite, via diagnose_fit(extra_plots=)


KEY_PARS <- c("log10_N50_inf","d_fev","gamma_inf","gamma_fevginf","log10_delta",
              "pi_susc","CoP_imm","alpha_inf","alpha_fevginf","CoP_susc")

# Example scenario set (all no-new-Stan-code: row filters + prior overrides).
# Every scenario names its tier. `tier` defaults to DEFAULT_SCENARIO_TIER rather
# than to a tier_col literal, so a scenario set can never silently describe a
# different configuration than the fit it is a sensitivity of.
DEFAULT_SCENARIO_TIER <- "t1-indiv"

SCENARIOS <- list(
  list(label = "base"),
  list(label = "hornick3_excluded", drop_obs = "H-F-3"),                 # MC4 sensitivity
  list(label = "gilrest_excluded",  drop_obs = "Gil-F-rest"),           # derived-by-subtraction row
  list(label = "delta_prior_lo",    prior_overrides = list(log10_delta = list(mu = 3.0))),
  list(label = "delta_prior_hi",    prior_overrides = list(log10_delta = list(mu = 4.0)))
)

#' Run one scenario -> results/scenarios/<label>/ (diagnostics + results.json + loo.json).
run_scenario <- function(spec, mod, data_csv = "dose_response_data.csv",
                         priors0 = load_priors(),
                         chains = 4, warmup = 800, sampling = 800,
                         adapt_delta = 0.9, seed = 2024, figures = TRUE,
                         allow_blocked = FALSE) {
  priors <- apply_prior_overrides(priors0, spec$prior_overrides %||% list())
  key    <- spec$tier %||% DEFAULT_SCENARIO_TIER
  tspec  <- tier_spec(key)
  # build_tier_data() asserts the UNFILTERED tier against the registry first, so a
  # scenario's drop_obs cannot mask a change in the tier's composition.
  stan_data <- build_tier_data(key, data_csv, priors,
                               drop_obs = spec$drop_obs %||% character(),
                               keep_obs = spec$keep_obs,
                               allow_blocked = allow_blocked, mod = mod)
  obs <- attr(stan_data, "obs")
  # tier in the PATH: three tiers of scenarios can no longer collide on one label.
  out_dir <- file.path("results", "scenarios", key, spec$label)

  t0 <- Sys.time()
  fit <- mod$sample(data = stan_data, chains = chains, parallel_chains = chains,
                    iter_warmup = warmup, iter_sampling = sampling,
                    adapt_delta = adapt_delta, seed = seed, refresh = 0,
                    show_messages = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  man <- run_manifest("scenario", stan_data = stan_data, priors = priors,
                      model_name = spec$label, data_csv = data_csv,
                      sampler = list(chains = chains, iter_warmup = warmup,
                                     iter_sampling = sampling, adapt_delta = adapt_delta,
                                     seed = seed, prior_only = 0L),
                      extra = list(scenario = spec$label,
                                   prior_overrides = spec$prior_overrides %||% list(),
                                   drop_obs = spec$drop_obs %||% character()))
  diagnose_fit(fit, out_dir,
               pars = tier_report_pars(tspec),   # derived; was a hardcoded 11
               obs = obs, priors = priors, model_name = spec$label, elapsed_s = elapsed,
               stan_data = stan_data, manifest = man,
               extra_plots = if (isTRUE(figures))
                 model_figures_hook(stan_data, label = sprintf("%s / %s", key, spec$label)))

  lj <- compute_loo_units(fit, obs)
  if (!is.null(lj)) jsonlite::write_json(lj, file.path(out_dir, "loo.json"),
                                         auto_unbox = TRUE, pretty = TRUE, digits = 6)
  invisible(out_dir)
}

#' PSIS-LOO on correctly-grouped observation UNITS. Rows that share volunteers are not
#' independent and must not be separate LOO units: grouped rows are keyed by cohort,
#' individual (n=1) rows by subject. Tier 1: 80 rows -> 46 units (was 79).
compute_loo_units <- function(fit, obs) {
  if (!requireNamespace("loo", quietly = TRUE)) return(NULL)
  if (!"log_lik" %in% fit$metadata()$stan_variables) return(NULL)
  ll_arr <- fit$draws("log_lik")
  ll_all <- posterior::as_draws_matrix(ll_arr)
  # log_lik is [N_obs dose-response rows, then N_ladder phi0 threshold binomials].
  # SLICE EXPLICITLY: `unit` below has length nrow(obs), and indexing a wider matrix
  # with a shorter logical RECYCLES it -- which would silently fold ladder columns
  # into the first few units rather than erroring.
  n_obs <- nrow(obs)
  if (ncol(ll_all) < n_obs)
    stop("log_lik has ", ncol(ll_all), " columns for ", n_obs, " observations",
         call. = FALSE)
  ll <- ll_all[, seq_len(n_obs), drop = FALSE]
  # The ladder tail is deliberately NOT a LOO unit. Those three binomials are counts
  # over the SAME 20 Darton placebo TD+ subjects who already appear as ox_fevginf_indiv
  # rows, so they are not independent of existing units, and an aggregate threshold
  # count cannot be split per subject to merge properly. They belong in log_lik (so
  # target == lprior + sum(log_lik) and priorsense power-scales the whole likelihood)
  # but not in model comparison.
  # LOO units must be INDEPENDENT, and rows that share volunteers are not. Group by
  # who the people are, not by which row they came from:
  #   - grouped rows (n > 1): the cohort. Lev-F-k and Lev-I-k are the same men, as are
  #     H-I-7 and H-FgI-7 (nested); previously only the Hornick pair was merged, so
  #     Levine's four pairs were each counted twice.
  #   - individual rows (n = 1): the SUBJECT, not the cohort. Darton's 56 rows are 30
  #     different men, so cohort-level merging would wrongly collapse them into one.
  #     D-I-plac-102 and D-FgI-plac-102 are one man and do merge.
  # Conservative where overlap is real but unresolvable: the Gilman control cohort
  # becomes one unit, because Gil-I-ctrl (43 of 64) partially overlaps all three
  # disjoint H-strata and no cross-tab exists to separate them.
  # NOTE this corrects MODEL COMPARISON only. The posterior itself still treats
  # Lev-F-k and Lev-I-k as independent binomials on the same men -- a known and
  # tolerated double-count. See ../joint_inference_plan.md Sec 6.5 / Sec 8.1.
  unit <- ifelse(obs$n == 1L,
                 paste0(obs$cohort_id, "::subj-", sub("^.*-", "", obs$obs_id)),
                 obs$cohort_id)
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
       # Comparability key. elpd is comparable only between fits over the IDENTICAL
       # observation set, and the unit set is too coarse a test: dropping Gil-F-rest
       # removes a row but not a unit (it shares a cohort_id), leaving the unit keys
       # identical while that unit now sums one fewer term -- a mechanically higher elpd
       # that passes a unit-set check. Hash the obs_ids.
       unit_keys_md5 = digest_chr(sort(u)),
       data_keys_md5 = digest_chr(sort(obs$obs_id)),
       elpd_loo = lo$estimates["elpd_loo", "Estimate"],
       elpd_loo_se = lo$estimates["elpd_loo", "SE"],
       p_loo = lo$estimates["p_loo", "Estimate"],
       n_pareto_k_gt_0.7 = sum(loo::pareto_k_values(lo) > 0.7))
}

#' Aggregate results.json (+ loo.json) across scenarios into a comparison table + forest plot.
summarize_scenarios <- function(labels, tier = DEFAULT_SCENARIO_TIER,
                                reference = "base",
                                out_dir = file.path("results", "scenarios", tier,
                                                    "cross_scenario_comparisons")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  recs <- lapply(labels, function(lab) {
    rj <- file.path("results", "scenarios", tier, lab, "results.json")
    if (!file.exists(rj)) return(NULL)
    r <- jsonlite::read_json(rj, simplifyVector = FALSE)
    lj <- file.path("results", "scenarios", tier, lab, "loo.json")
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
    n_units = if (!is.null(e$loo)) e$loo$n_units else NA,
    data_keys_md5 = if (!is.null(e$loo)) e$loo$data_keys_md5 %||% NA_character_ else NA,
    elpd_loo = if (!is.null(e$loo)) e$loo$elpd_loo else NA,
    elpd_loo_se = if (!is.null(e$loo)) e$loo$elpd_loo_se else NA,
    pareto_k_gt_0.7 = if (!is.null(e$loo)) e$loo$n_pareto_k_gt_0.7 else NA)))

  # elpd is comparable ONLY across identical LOO unit sets. Row filters change the set,
  # so their elpd is mechanically higher (fewer terms summed) and must not be read
  # against the reference. Blank it rather than print a number that invites the
  # comparison -- the header caveat alone was not enough.
  ref <- scen$data_keys_md5[match(reference, scen$scenario)]
  if (length(ref) && !is.na(ref)) {
    scen$loo_comparable <- !is.na(scen$data_keys_md5) & scen$data_keys_md5 == ref
    scen$elpd_loo[!scen$loo_comparable]    <- NA
    scen$elpd_loo_se[!scen$loo_comparable] <- NA
  } else {
    scen$loo_comparable <- NA
  }

  readr::write_csv(long, file.path(out_dir, "comparison.csv"))
  readr::write_csv(scen, file.path(out_dir, "comparison_scenarios.csv"))

  # comparison.md
  L <- c("# Scenario comparison", "",
         "## Sampler health + LOO", "",
         knitr_table(scen), "",
         paste("`elpd_loo` is blank where a scenario's LOO unit set differs from",
               sprintf("`%s`", reference), "-- a row filter changes the units, so its",
               "elpd is mechanically higher and is NOT comparable. Keyed on the",
               "OBSERVATION set (`data_keys_md5`), not the unit set: dropping a row that",
               "shares a cohort leaves the unit set unchanged while still shrinking the",
               "sum."), "",
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
  setwd(calib_dir())
  source("tier_specs.R")
  args <- commandArgs(trailingOnly = TRUE)
  tier <- if (length(args)) args[1] else DEFAULT_SCENARIO_TIER
  mod <- cmdstan_model("typhoid_dose_response.stan")
  priors0 <- load_priors()
  for (spec in SCENARIOS) {
    spec$tier <- spec$tier %||% tier
    cat(sprintf("\n=== scenario: %s / %s ===\n", spec$tier, spec$label))
    run_scenario(spec, mod, priors0 = priors0)
  }
  summarize_scenarios(vapply(SCENARIOS, function(s) s$label, character(1)), tier)
}
