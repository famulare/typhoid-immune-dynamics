#' Reusable Stan-fit diagnostics for the dose-response calibration.
#'
#' Buffalo-style workflow, R-native (bayesplot + posterior + loo + priorsense).
#' `diagnose_fit()` takes any cmdstanr fit and emits a standard battery:
#'   - sampler diagnostics table (rhat, ess_bulk/tail, mcse) + recovery vs truth
#'   - divergence / treedepth / E-BFMI counts
#'   - plot battery: trace, density, pairs (divergence-highlighted), rank, energy
#'   - prior-vs-posterior overlay (generic over prior family, from priors.yaml)
#'   - rate-space grouped-binomial PPC (from Stan p_pred -- no R likelihood mirror)
#'   - priorsense power-scaling sensitivity (needs lprior + log_lik in the fit)
#'   - summary.md (human) + results.json (machine) + fit.rds
#'
#' Sourced by fit_dose_response.R, simulate_recovery.R, run_scenarios.R.

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(bayesplot)
  library(ggplot2); library(dplyr); library(tidyr)
})

# Mike-native theming: theme_bw() everywhere (ggplot + bayesplot), small titles.
.dr_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                 plot.subtitle = ggplot2::element_text(size = 8),
                 strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey50"))
ggplot2::theme_set(.dr_theme)
bayesplot::bayesplot_theme_set(.dr_theme)
bayesplot::color_scheme_set("blue")

# ---- small helpers ------------------------------------------------------------

.save_gg <- function(p, path, w = 9, h = 6) {
  tryCatch(ggplot2::ggsave(path, p, width = w, height = h, dpi = 150, bg = "white"),
           error = function(e) message("  [skip] ", basename(path), ": ", conditionMessage(e)))
}

# bayesplot grids (e.g. mcmc_pairs) are gtables, not ggplots -> device print.
.save_grid <- function(grob, path, w = 10, h = 10) {
  tryCatch({
    grDevices::png(path, width = w, height = h, units = "in", res = 150, bg = "white")
    on.exit(grDevices::dev.off())
    if (inherits(grob, "ggplot")) print(grob) else gridExtra::grid.arrange(grob)
  }, error = function(e) message("  [skip] ", basename(path), ": ", conditionMessage(e)))
}

#' Sampler diagnostic summary table (+ recovery columns when truth supplied).
#' @param true_params named numeric of true values (recovery runs only) or NULL.
diagnostic_table <- function(fit, pars, true_params = NULL) {
  draws <- fit$draws(pars)
  # String function names (not bare mean/median/sd) so a caller variable named
  # `sd`/`mean` cannot shadow them.
  tab <- posterior::summarise_draws(
    draws, "mean", "median", "sd",
    ~posterior::quantile2(.x, probs = c(0.05, 0.95)),
    "rhat", "ess_bulk", "ess_tail", "mcse_mean"
  )
  if (!is.null(true_params)) {
    tv <- true_params[tab$variable]
    tab$true <- as.numeric(tv)
    tab$recovered_90 <- !is.na(tab$true) & tab$q5 <= tab$true & tab$true <= tab$q95
  }
  tab
}

#' Divergences / max-treedepth / E-BFMI.
sampler_health <- function(fit) {
  ds <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
  sd_arr <- fit$sampler_diagnostics()
  ndraws <- prod(dim(sd_arr)[1:2])
  div <- if ("divergent__" %in% dimnames(sd_arr)$variable) sum(sd_arr[, , "divergent__"]) else NA
  list(
    n_div = as.integer(div),
    n_draws = as.integer(ndraws),
    div_rate = if (is.na(div)) NA else div / ndraws,
    n_max_td = if (!is.null(ds)) sum(ds$num_max_treedepth) else NA_integer_,
    ebfmi = if (!is.null(ds)) ds$ebfmi else NA
  )
}

# ---- plots --------------------------------------------------------------------

plot_battery <- function(fit, out_dir, pars, true_params = NULL) {
  draws <- fit$draws(pars)
  np <- tryCatch(bayesplot::nuts_params(fit), error = function(e) NULL)

  .save_gg(bayesplot::mcmc_trace(draws), file.path(out_dir, "trace.png"),
           w = 10, h = 2 + 1.1 * length(pars))
  .save_gg(bayesplot::mcmc_dens_overlay(draws), file.path(out_dir, "density.png"),
           w = 10, h = 2 + 1.1 * length(pars))
  .save_gg(bayesplot::mcmc_rank_overlay(draws), file.path(out_dir, "rank.png"),
           w = 10, h = 2 + 1.1 * length(pars))
  if (!is.null(np)) {
    .save_grid(bayesplot::mcmc_pairs(draws, np = np, off_diag_args = list(size = 0.6, alpha = 0.3)),
               file.path(out_dir, "pairs.png"),
               w = 2 + 1.6 * length(pars), h = 2 + 1.6 * length(pars))
    .save_gg(bayesplot::mcmc_nuts_energy(np), file.path(out_dir, "energy.png"), w = 8, h = 4)
  }
  # posterior::rhat() expects one variable; on a multi-var draws_array it collapses
  # to a single scalar. Use summarise_draws for per-parameter R-hat (matches summary.md).
  rh <- tryCatch({
    rt <- posterior::summarise_draws(draws, "rhat")
    stats::setNames(rt$rhat, rt$variable)
  }, error = function(e) NULL)
  if (!is.null(rh)) .save_gg(bayesplot::mcmc_rhat(rh) + bayesplot::yaxis_text(),
                             file.path(out_dir, "rhat.png"), w = 7, h = 4)
  ne <- tryCatch(bayesplot::neff_ratio(fit, pars = pars), error = function(e) NULL)
  if (!is.null(ne)) .save_gg(bayesplot::mcmc_neff(ne) + bayesplot::yaxis_text(),
                             file.path(out_dir, "neff.png"), w = 7, h = 4)
}

#' Prior-vs-posterior overlay (generic over family, from priors.yaml).
plot_prior_posterior <- function(fit, out_dir, priors, true_params = NULL) {
  pp <- intersect(names(priors), fit$metadata()$stan_variables)
  if (!length(pp)) return(invisible())
  post <- posterior::as_draws_df(fit$draws(pp))
  rows <- lapply(pp, function(par) {
    pv <- post[[par]]
    rng <- quantile(pv, c(0.001, 0.999))
    grid <- seq(rng[1], rng[2], length.out = 256)
    dplyr::bind_rows(
      data.frame(param = par, value = pv, src = "posterior", dens = NA_real_),
      data.frame(param = par, value = grid, src = "prior",
                 dens = prior_density(priors, par, grid))
    )
  })
  df <- dplyr::bind_rows(rows)
  p <- ggplot() +
    geom_density(data = dplyr::filter(df, src == "posterior"),
                 aes(value, after_stat(density)), fill = "#4878cf", alpha = 0.45, colour = NA) +
    geom_area(data = dplyr::filter(df, src == "prior"),
              aes(value, dens), fill = "grey60", alpha = 0.25) +
    geom_line(data = dplyr::filter(df, src == "prior"), aes(value, dens),
              colour = "grey40", linetype = 2) +
    facet_wrap(~param, scales = "free") +
    labs(title = "Prior (grey) vs posterior (blue)", x = NULL, y = "density")
  if (!is.null(true_params)) {
    tdf <- data.frame(param = names(true_params), value = as.numeric(true_params)) |>
      dplyr::filter(param %in% pp)
    if (nrow(tdf)) p <- p + geom_vline(data = tdf, aes(xintercept = value),
                                       colour = "red", linewidth = 0.7)
  }
  .save_gg(p, file.path(out_dir, "prior_posterior.png"),
           w = 11, h = 2.5 * ceiling(length(pp) / 4))
}

#' Rate-space grouped-binomial PPC, fed by Stan p_pred (replaces the R mirror).
plot_ppc <- function(fit, out_dir, obs, true_params = NULL) {
  if (is.null(obs) || !"p_pred" %in% fit$metadata()$stan_variables) return(invisible())
  pp <- posterior::as_draws_matrix(fit$draws("p_pred"))   # ndraws x N_obs
  q <- t(apply(pp, 2, quantile, probs = c(0.05, 0.5, 0.95)))
  df <- dplyr::bind_cols(obs, fit_lo = q[, 1], fit_med = q[, 2], fit_hi = q[, 3])
  p <- ggplot(df, aes(obs_rate, fit_med, colour = likelihood_group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_linerange(aes(ymin = fit_lo, ymax = fit_hi), alpha = 0.5) +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(aes(label = obs_id), size = 2.5, max.overlaps = 20,
                             show.legend = FALSE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Observed attack rate (y/n)",
         y = "Posterior fitted probability (median, 90% CI)",
         title = "Posterior predictive check (rate space)", colour = "group")
  .save_gg(p, file.path(out_dir, "ppc.png"), w = 8, h = 6)
  invisible(df)
}

#' priorsense power-scaling sensitivity (needs lprior + log_lik in the fit).
#' Returns a data.frame (variable, prior, likelihood, diagnosis), filtered to `pars`.
powerscale_table <- function(fit, pars = NULL) {
  if (!requireNamespace("priorsense", quietly = TRUE)) return(NULL)
  tryCatch({
    df <- as.data.frame(priorsense::powerscale_sensitivity(fit))
    if (!is.null(pars)) df <- df[df$variable %in% pars, , drop = FALSE]
    df
  }, error = function(e) { message("  [skip] priorsense: ", conditionMessage(e)); NULL })
}

# ---- writers ------------------------------------------------------------------

.fmt <- function(x, d = 4) formatC(x, format = "g", digits = d)

write_summary_md <- function(out_dir, model_name, health, tab, corr, ps_tab) {
  L <- c(sprintf("# Fit summary: %s", model_name), "",
         sprintf("**Date:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
         "## Sampler diagnostics", "",
         "| Metric | Value |", "|---|---|",
         sprintf("| Divergent transitions | %s / %s |", health$n_div, health$n_draws),
         sprintf("| Divergence rate | %s%% |", .fmt(100 * health$div_rate, 3)),
         sprintf("| Max-treedepth hits | %s |", health$n_max_td),
         sprintf("| E-BFMI (min) | %s |", .fmt(suppressWarnings(min(health$ebfmi)), 3)), "",
         "## Parameters", "")
  hdr <- if ("true" %in% names(tab))
    "| param | true | mean | median | sd | 5% | 95% | ess_bulk | ess_tail | rhat | recovered |"
  else
    "| param | mean | median | sd | 5% | 95% | ess_bulk | ess_tail | rhat |"
  sep <- paste0("|", paste(rep("---", lengths(gregexpr("\\|", hdr)) - 1), collapse = "|"), "|")
  L <- c(L, hdr, sep)
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    if ("true" %in% names(tab)) {
      L <- c(L, sprintf("| `%s` | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                        r$variable, .fmt(r$true), .fmt(r$mean), .fmt(r$median), .fmt(r$sd),
                        .fmt(r$q5), .fmt(r$q95), .fmt(r$ess_bulk, 0), .fmt(r$ess_tail, 0),
                        .fmt(r$rhat, 4), ifelse(isTRUE(r$recovered_90), "OK", "**MISS**")))
    } else {
      L <- c(L, sprintf("| `%s` | %s | %s | %s | %s | %s | %s | %s | %s |",
                        r$variable, .fmt(r$mean), .fmt(r$median), .fmt(r$sd),
                        .fmt(r$q5), .fmt(r$q95), .fmt(r$ess_bulk, 0), .fmt(r$ess_tail, 0),
                        .fmt(r$rhat, 4)))
    }
  }
  if ("recovered_90" %in% names(tab)) {
    ok <- sum(tab$recovered_90, na.rm = TRUE); tot <- sum(!is.na(tab$true))
    L <- c(L, "", sprintf("**Recovery:** %d/%d within 90%% CI", ok, tot))
  }
  if (!is.null(corr) && nrow(corr)) {
    L <- c(L, "", "## Strong pairwise correlations (|r| > 0.7)", "",
           "| pair | r |", "|---|---|",
           sprintf("| `%s` - `%s` | %s |", corr$a, corr$b, .fmt(corr$r, 3)))
  }
  if (!is.null(ps_tab) && nrow(ps_tab)) {
    L <- c(L, "", "## priorsense power-scaling sensitivity", "",
           "| variable | prior | likelihood | diagnosis |", "|---|---|---|---|",
           sprintf("| `%s` | %s | %s | %s |", ps_tab$variable,
                   .fmt(ps_tab$prior, 3), .fmt(ps_tab$likelihood, 3), ps_tab$diagnosis))
  }
  L <- c(L, "", "## Figures", "",
         vapply(sort(basename(Sys.glob(file.path(out_dir, "*.png")))),
                function(f) sprintf("![%s](%s)\n", tools::file_path_sans_ext(f), f), character(1)))
  writeLines(L, file.path(out_dir, "summary.md"))
}

write_results_json <- function(out_dir, model_name, health, tab, corr, elapsed_s) {
  params <- setNames(lapply(seq_len(nrow(tab)), function(i) {
    r <- tab[i, ]
    base <- list(mean = r$mean, median = r$median, sd = r$sd, ci90 = c(r$q5, r$q95),
                 ess_bulk = r$ess_bulk, ess_tail = r$ess_tail, rhat = r$rhat)
    if ("true" %in% names(tab)) { base$true <- r$true; base$recovered_90 <- r$recovered_90 }
    base
  }), tab$variable)
  res <- list(model = model_name, wall_time_s = elapsed_s,
              sampler = list(n_divergent = health$n_div, n_draws = health$n_draws,
                             div_rate = health$div_rate, n_max_treedepth = health$n_max_td,
                             ebfmi_min = suppressWarnings(min(health$ebfmi))),
              parameters = params,
              correlations = if (!is.null(corr) && nrow(corr))
                setNames(as.list(corr$r), paste(corr$a, corr$b, sep = ",")) else list())
  if ("recovered_90" %in% names(tab))
    res$recovery <- list(n_recovered_90 = sum(tab$recovered_90, na.rm = TRUE),
                         n_total = sum(!is.na(tab$true)))
  jsonlite::write_json(res, file.path(out_dir, "results.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 8, na = "null")
}

#' Strong pairwise correlations among `pars`.
strong_correlations <- function(fit, pars, threshold = 0.7) {
  m <- posterior::as_draws_matrix(fit$draws(pars))
  cm <- suppressWarnings(cor(m))
  out <- data.frame(a = character(), b = character(), r = numeric())
  for (i in seq_along(pars)) for (j in seq_along(pars)) if (j > i) {
    r <- cm[pars[i], pars[j]]
    if (is.finite(r) && abs(r) > threshold) out <- rbind(out, data.frame(a = pars[i], b = pars[j], r = r))
  }
  out[order(-abs(out$r)), ]
}

# ---- main entry point ---------------------------------------------------------

#' @param fit cmdstanr CmdStanMCMC.
#' @param out_dir output directory (created).
#' @param pars character vector of interpretable parameters to diagnose/plot.
#' @param true_params named numeric of truths (recovery runs) or NULL (real data).
#' @param obs data frame from attr(stan_data,"obs") for the PPC, or NULL.
#' @param priors parsed priors (load_priors()) for the overlay, or NULL.
#' @param model_name label for the summary.
#' @param elapsed_s wall time, if known.
diagnose_fit <- function(fit, out_dir, pars, true_params = NULL, obs = NULL,
                         priors = NULL, model_name = "fit", elapsed_s = NA_real_,
                         corr_exclude = c("N50_inf", "N50_fevginf", "delta")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pars <- intersect(pars, fit$metadata()$stan_variables)

  health <- sampler_health(fit)
  tab <- diagnostic_table(fit, pars, true_params)
  corr <- strong_correlations(fit, setdiff(pars, corr_exclude))  # drop redundant log<->linear transforms
  ps_tab <- powerscale_table(fit, pars)

  plot_battery(fit, out_dir, pars, true_params)
  if (!is.null(priors)) plot_prior_posterior(fit, out_dir, priors, true_params)
  ppc_df <- plot_ppc(fit, out_dir, obs, true_params)

  write_summary_md(out_dir, model_name, health, tab, corr, ps_tab)
  write_results_json(out_dir, model_name, health, tab, corr, elapsed_s)
  tryCatch(fit$save_object(file.path(out_dir, "fit.rds")), error = function(e) NULL)
  if (!is.null(ppc_df)) readr::write_csv(ppc_df, file.path(out_dir, "ppc.csv"))
  readr::write_csv(tab, file.path(out_dir, "summary.csv"))

  cat(sprintf("\nResults saved to: %s/\n", normalizePath(out_dir)))
  cat("  summary.md / summary.csv / results.json — diagnostics + parameter tables\n")
  cat("  trace.png pairs.png rank.png energy.png density.png rhat.png neff.png\n")
  cat("  prior_posterior.png ppc.png — model checks\n")
  cat("  fit.rds — saved fit object\n")
  cat(sprintf("  divergences: %d/%d (%.2f%%)  min E-BFMI: %.2f\n",
              health$n_div, health$n_draws, 100 * health$div_rate,
              suppressWarnings(min(health$ebfmi))))
  invisible(list(table = tab, health = health, corr = corr, powerscale = ps_tab))
}
