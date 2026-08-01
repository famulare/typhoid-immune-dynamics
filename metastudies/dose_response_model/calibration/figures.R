#' Orchestrator for the whole model figure suite.
#'
#' Sourcing this file gives you everything: the math (model_math.R), the grouping
#' spec (curve_specs.R), and all three figure suites. One call, make_model_figures(),
#' writes the complete set into a run directory.
#'
#' Wired into diagnose_fit(extra_plots=) so every run directory gets the suite --
#' results/tier1, results/tier1_prior, results/scenarios/<label>, and the recovery
#' harness -- rather than only the one directory the old callsite named. That is why
#' results/tier1_prior/ previously had no dose-response figure at all.
#'
#' Standalone:
#'   Rscript figures.R                       # every run dir with a fit.rds
#'   Rscript figures.R results/tier1         # one run dir

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists("calib_dir")) source("utils.R")
source("model_math.R")
source("curve_specs.R")
source("figures_common.R")
source("dose_response_curves.R")
source("figures_grid.R")
source("figures_phi.R")
source("figures_components.R")
source("figures_targets.R")

#' Write the whole suite for one fit into one directory.
#'
#' @param verify assert the figure math reproduces this fit's own p_pred before
#'   drawing anything. Leave TRUE: it is the runtime half of the parity gate.
make_model_figures <- function(fit, stan_data, out_dir, label = basename(out_dir),
                               specs = curve_specs(), ngrid = 60, n_spaghetti = 100,
                               ndraw_ribbon = 1000, seed = 1, verify = TRUE,
                               write_columns = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (verify) {
    d <- curve_parity_check(fit, stan_data)
    jsonlite::write_json(list(max_abs_diff_vs_p_pred = d, checked = Sys.time()),
                         file.path(out_dir, "curve_parity.json"), auto_unbox = TRUE)
    message(sprintf("  curve parity vs p_pred: %.3e", d))
  }
  # saved so figures_from_dir() can reproduce a scenario's exact row selection
  saveRDS(stan_data, file.path(out_dir, "stan_data.rds"))

  plot_calibration_targets(fit, stan_data, file.path(out_dir, "calibration_targets.png"),
                           label = label)
  opts <- dose_response_plot_options(out_dir)
  plot_dose_response_fit(fit, stan_data, file.path(out_dir, "dose_response_fit.png"),
                         show_points = opts$show_points,
                         show_errorbars = opts$show_errorbars,
                         fixed_cop = opts$fixed_cop,
                         include_phi_panel = opts$include_phi_panel)
  plot_titre_protection(fit, stan_data, file.path(out_dir, "titre_protection.png"))
  plot_grouping_grid(fit, stan_data, out_dir, specs, ngrid, n_spaghetti,
                     ndraw_ribbon, seed, label, write_columns)
  plot_phi_severity(fit, stan_data, file.path(out_dir, "phi_severity.png"),
                    n_spaghetti, ndraw_ribbon, seed, label)
  plot_maryland_mixture(fit, stan_data, file.path(out_dir, "maryland_mixture.png"),
                        ndraw_ribbon = ndraw_ribbon, seed = seed, label = label)
  plot_delta_bridge(fit, stan_data, file.path(out_dir, "delta_bridge.png"),
                    ndraw_ribbon = ndraw_ribbon, seed = seed, label = label)
  plot_dose_cop_surface(fit, stan_data, file.path(out_dir, "dose_cop_surface.png"),
                        ndraw_ribbon = min(500, ndraw_ribbon), seed = seed, label = label)
  plot_eta_detection(fit, stan_data, file.path(out_dir, "eta_detection.png"),
                     ndraw_ribbon = ndraw_ribbon, seed = seed, label = label)
  invisible(out_dir)
}

#' The closure to hand to diagnose_fit(extra_plots = ...).
model_figures_hook <- function(stan_data, ...) {
  force(stan_data)
  function(fit, out_dir) make_model_figures(fit, stan_data, out_dir, ...)
}

#' Rebuild every figure for a saved run directory, without refitting.
figures_from_dir <- function(run_dir, data_csv = "dose_response_data.csv", ...) {
  fit <- readRDS(file.path(run_dir, "fit.rds"))
  # resolve_run_stan_data() prefers this run's stan_data.rds and otherwise rebuilds
  # from the tier recorded in its manifest. The old fallback rebuilt tier1_active
  # unconditionally and then hard-errored on the row count -- it guessed, and the
  # guess was wrong for every scenario and recovery dir.
  sd <- resolve_run_stan_data(run_dir, data_csv)
  n_pred <- ncol(posterior::as_draws_matrix(fit$draws("p_pred")))
  if (nrow(attr(sd, "obs")) != n_pred)
    stop("figures_from_dir: stan_data does not match this fit (", run_dir, "): ",
         nrow(attr(sd, "obs")), " obs vs ", n_pred, " p_pred columns.", call. = FALSE)
  make_model_figures(fit, sd, run_dir, ...)
}

if (sys.nframe() == 0) {
  setwd(calib_dir())
  source("priors.R"); source("data_prep.R"); source("tier_specs.R")
  suppressPackageStartupMessages(library(cmdstanr))
  args <- commandArgs(trailingOnly = TRUE)
  # Default: every REGENERABLE run dir, discovered from the audit rather than from a
  # hardcoded list -- which is how results/recovery/* used to be silently skipped.
  dirs <- if (length(args)) args else {
    a <- audit_run_dirs(quiet = TRUE)
    if (is.null(a)) character(0) else a$dir[a$regenerable]
  }
  for (d in dirs) {
    message("=== ", d, " ===")
    tryCatch(figures_from_dir(d, label = basename(d)),
             error = function(e) message("  [FAILED] ", conditionMessage(e)))
  }
}
