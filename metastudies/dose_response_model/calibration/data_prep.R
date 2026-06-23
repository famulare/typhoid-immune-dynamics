#' Flat Stan data assembly for the dose-response calibration.
#'
#' Builds the unified per-observation layout the refactored typhoid_dose_response.stan
#' expects (one row per observation + a `group` code), merged with the prior
#' hyperparameter data list from priors.yaml. Shared by fit_dose_response.R,
#' simulate_recovery.R, and run_scenarios.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# likelihood_group -> Stan group code (must match obs_prob() in the .stan)
.GROUP_CODE <- c(ox_fev = 1L, ox_inf = 2L, md_fev = 3L, md_inf = 4L, hornick_cond = 5L)

#' Apply nested prior overrides (for sensitivity scenarios), e.g.
#'   list(log10_delta = list(mu = 3.0))  ->  fixes the delta prior mean.
apply_prior_overrides <- function(priors, overrides = list()) {
  for (param in names(overrides)) {
    if (is.null(priors[[param]])) stop(sprintf("override for unknown prior '%s'", param))
    for (h in names(overrides[[param]])) priors[[param]][[h]] <- overrides[[param]][[h]]
  }
  priors
}

#' @param data_csv Path to dose_response_data.csv.
#' @param priors Parsed priors (from load_priors()), optionally override-mutated.
#' @param tier_col Which activation column selects rows ("tier1_active" or "tier2_active").
#' @param prior_only 1 to skip the likelihood (prior predictive).
#' @param drop_obs Character vector of obs_id to exclude (data-filter sensitivities).
#' @param keep_obs If non-NULL, restrict to these obs_id (e.g. exclude-Oxford/Maryland).
#' @return list: the flat Stan data + the prior `pr_*` scalars, plus an attribute
#'   "obs" carrying the selected data frame (obs_id, study, group, etc.) for plotting.
build_stan_data <- function(data_csv, priors,
                            tier_col = "tier1_active",
                            prior_only = 0L,
                            drop_obs = character(),
                            keep_obs = NULL) {
  d <- readr::read_csv(data_csv, show_col_types = FALSE)
  dat <- d %>% filter(.data[[tier_col]] == 1)
  if (!is.null(keep_obs)) dat <- dat %>% filter(obs_id %in% keep_obs)
  if (length(drop_obs))   dat <- dat %>% filter(!obs_id %in% drop_obs)
  dat <- dat %>% arrange(match(likelihood_group, names(.GROUP_CODE)), obs_id)

  grp <- unname(.GROUP_CODE[dat$likelihood_group])
  if (anyNA(grp)) stop("unmapped likelihood_group: ",
                       paste(unique(dat$likelihood_group[is.na(grp)]), collapse = ", "))

  stan_data <- list(
    N_obs   = nrow(dat),
    group   = as.integer(grp),
    n       = as.integer(dat$n),
    y       = as.integer(dat$y),
    dose    = as.numeric(dat$dose_cfu),
    CoP     = ifelse(is.na(dat$CoP), 1.0, as.numeric(dat$CoP)),       # used by ox groups only
    phi     = ifelse(is.na(dat$phi), 1.0, as.numeric(dat$phi)),       # used by md_fev/hornick only
    stratum = ifelse(is.na(dat$gilman_stratum), 0L, as.integer(dat$gilman_stratum)),
    prior_only = as.integer(prior_only)
  )
  stan_data <- c(stan_data, priors_to_stan_data(priors))

  stopifnot(!anyNA(unlist(stan_data[c("dose", "CoP", "phi", "n", "y", "group", "stratum")])))
  attr(stan_data, "obs") <- dat %>%
    transmute(obs_id, study, likelihood_group, group = grp,
              dose_cfu, n, y, obs_rate = y / n, CoP, phi, gilman_stratum)
  stan_data
}
