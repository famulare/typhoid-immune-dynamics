#' Shared, dependency-free helpers for the dose-response calibration.
#'
#' Deliberately imports nothing: every other file in this directory can source it
#' without pulling in cmdstanr/ggplot2/yaml. Everything here previously existed in
#' two to five copies; see TIER_LOCK.md for why that mattered.
#'
#' Source it with the repo's lazy idiom so double-sourcing is free:
#'   if (!exists("calib_dir")) source("utils.R")

#' Directory containing this calibration code, however the caller was invoked.
#'
#' Three cases, in order: Rscript (`--file=`), sourced from another file
#' (`sys.frames()[[1]]$ofile`), interactive (cwd). Consolidates five near-identical
#' copies -- get_script_dir() in the driver, .priors_default_path() in priors.R,
#' load_reference_points() in diagnostics.R, and inline copies in figures.R and
#' dose_response_curves.R, two of which handled only the Rscript case and crashed
#' when sourced interactively.
calib_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg))
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  getwd()
}

#' Path to a file beside the calibration code.
calib_path <- function(...) file.path(calib_dir(), ...)

#' NULL-coalescing. Was defined identically in figures_grid.R, figures_phi.R and
#' run_scenarios.R while being CONSUMED by figures_components.R and parts of
#' figures_grid.R that define it nowhere -- so correctness depended on figures.R's
#' source order. One definition removes that coupling.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Wilson score interval for a binomial proportion.
#'
#' Lived in dose_response_curves.R while being used by figures_grid.R,
#' figures_components.R, figures_phi.R and figures_targets.R -- the last of which
#' sourced dose_response_curves.R solely to borrow it.
.wilson <- function(y, n, z = 1.96) {
  p <- y / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = pmax(0, ctr - hw), hi = pmin(1, ctr + hw))
}

#' Minimal markdown table (avoids a knitr dependency).
#' Moved out of run_scenarios.R so audit_run_dirs() and `fit_tier.R --list` can use
#' it without sourcing the scenario harness (and therefore cmdstanr).
knitr_table <- function(df) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1, function(r) paste0("| ", paste(format(r, trim = TRUE), collapse = " | "), " |"))
  c(hdr, sep, body)
}

#' Carry the Stan-data attributes across a modifyList()/element-swap.
#'
#' modifyList() drops attributes, so flipping `prior_only` or swapping in synthetic
#' `y` silently loses attr "obs" (and now attr "tier") -- the row order that IS the
#' p_pred column order, and the tier identity the provenance manifest reads. Every
#' such site re-attached "obs" by hand, or forgot to.
copy_stan_attrs <- function(new, old, which = c("obs", "tier")) {
  for (a in which) if (!is.null(attr(old, a))) attr(new, a) <- attr(old, a)
  new
}
