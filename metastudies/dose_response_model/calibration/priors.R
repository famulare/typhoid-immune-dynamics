#' Prior single-source utilities for the dose-response calibration.
#'
#' Reads calibration/priors.yaml (the single source of truth) and provides:
#'   load_priors()          - parse the yaml into a named list
#'   priors_to_stan_data()  - build the Stan `pr_*` hyperparameter data list
#'   sample_prior()         - draw from a parameter's prior (generic over family)
#'   prior_density()        - prior density at x (generic over family)
#'   prior_median()         - prior median (handy as a recovery-harness truth)
#'
#' Families: normal (mu, sd) | half_normal (mu, sd, lower=0) |
#'           lognormal (meanlog, sdlog) | beta (a, b) | exponential (rate, lower=0).

suppressPackageStartupMessages(library(yaml))
if (!exists("calib_dir")) source("utils.R")   # calib_path()

#' @param path Path to priors.yaml. Defaults to the file beside this script.
load_priors <- function(path = calib_path("priors.yaml")) {
  p <- yaml::read_yaml(path)
  for (nm in names(p)) {
    fam <- p[[nm]]$family
    if (is.null(fam) || !fam %in% c("normal", "half_normal", "lognormal", "beta", "exponential"))
      stop(sprintf("priors.yaml: parameter '%s' has unknown/missing family '%s'", nm, fam))
  }
  p
}

#' Map the parsed priors to the exact `pr_*` scalar names the Stan data block expects.
priors_to_stan_data <- function(priors) {
  out <- list()
  for (nm in names(priors)) {
    spec <- priors[[nm]]
    fam <- spec$family
    if (fam %in% c("normal", "half_normal")) {
      out[[paste0("pr_", nm, "_mu")]] <- as.numeric(spec$mu)
      out[[paste0("pr_", nm, "_sd")]] <- as.numeric(spec$sd)
    } else if (fam == "lognormal") {
      out[[paste0("pr_", nm, "_mu")]] <- as.numeric(spec$meanlog)
      out[[paste0("pr_", nm, "_sd")]] <- as.numeric(spec$sdlog)
    } else if (fam == "beta") {
      out[[paste0("pr_", nm, "_a")]] <- as.numeric(spec$a)
      out[[paste0("pr_", nm, "_b")]] <- as.numeric(spec$b)
    } else if (fam == "exponential") {
      out[[paste0("pr_", nm, "_rate")]] <- as.numeric(spec$rate)
    }
  }
  out
}

#' Draw `n` samples from parameter `param`'s prior.
sample_prior <- function(priors, param, n = 10000L) {
  spec <- priors[[param]]
  if (is.null(spec)) stop(sprintf("no prior for '%s'", param))
  switch(spec$family,
    normal      = rnorm(n, spec$mu, spec$sd),
    half_normal = qnorm(runif(n, pnorm(0, spec$mu, spec$sd), 1), spec$mu, spec$sd),  # lower-truncated at 0
    lognormal   = rlnorm(n, spec$meanlog, spec$sdlog),
    beta        = rbeta(n, spec$a, spec$b),
    exponential = rexp(n, spec$rate)
  )
}

#' Prior density of `param` at `x`.
prior_density <- function(priors, param, x) {
  spec <- priors[[param]]
  switch(spec$family,
    normal      = dnorm(x, spec$mu, spec$sd),
    half_normal = ifelse(x < 0, 0, dnorm(x, spec$mu, spec$sd) / (1 - pnorm(0, spec$mu, spec$sd))),
    lognormal   = dlnorm(x, spec$meanlog, spec$sdlog),
    beta        = dbeta(x, spec$a, spec$b),
    exponential = dexp(x, spec$rate)
  )
}

#' Prior median (a convenient default "truth" for the recovery harness).
prior_median <- function(priors, param) {
  spec <- priors[[param]]
  switch(spec$family,
    normal      = spec$mu,
    half_normal = qnorm(0.75, spec$mu, spec$sd),       # median of the lower-truncated-at-0 normal (mu=0)
    lognormal   = exp(spec$meanlog),
    beta        = qbeta(0.5, spec$a, spec$b),
    exponential = qexp(0.5, spec$rate)
  )
}
