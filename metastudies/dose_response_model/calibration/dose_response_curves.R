#' Bespoke model diagnostic: posterior dose-response curves vs observed data.
#'
#' This is the model-specific PPC (in addition to the standard bayesplot/posterior
#' battery in diagnostics.R). It plots, per likelihood group, the posterior
#' dose-response curve (median + 90% ribbon) over a dose grid against the observed
#' attack rates (with Wilson 95% CIs) and the per-observation fitted `p_pred`.
#'
#' The curves are computed in R from posterior parameter draws, mirroring the
#' `obs_prob()` math in the .stan. This is a *visualization* mirror only — the
#' inference-side PPC (`p_pred`, plotted here as points) still comes from Stan, and
#' the likelihood's single-source guarantee (test_obs_prob_parity.R) is untouched.

suppressPackageStartupMessages({
  library(posterior); library(dplyr); library(tidyr); library(ggplot2)
})

# Beta-Poisson and Maryland mixture, vectorized over posterior draws (D_eff scalar).
.bp <- function(D, N50, alpha, CoP, gamma) {
  scale <- (2^(1 / alpha) - 1) / N50
  1 - (1 + D * scale)^(-alpha / CoP^gamma)
}
.wilson <- function(y, n, z = 1.96) {  # Wilson score interval
  p <- y / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = pmax(0, ctr - hw), hi = pmin(1, ctr + hw))
}

#' @param fit cmdstanr fit; @param stan_data list with attr "obs"; @param outfile png path
plot_dose_response_fit <- function(fit, stan_data, outfile) {
  obs <- attr(stan_data, "obs")
  dr  <- as_draws_df(fit$draws(c("N50_inf", "N50_fevginf", "alpha_inf", "alpha_fevginf",
                                 "gamma_inf", "gamma_fevginf", "delta",
                                 "pi_susc", "CoP_imm", "CoP_susc", "phi_md")))

  # population dose-response over a grid, per panel (median + 90% ribbon across draws)
  grid_curve <- function(doses, fn) {
    q <- vapply(doses, function(d) quantile(fn(d), c(.05, .5, .95)), numeric(3))
    tibble(dose_cfu = doses, lo = q[1, ], med = q[2, ], hi = q[3, ])
  }
  with(dr, {
    ox  <- 10^seq(2.3, 4.7, length.out = 60)
    md  <- 10^seq(2.7, 9.7, length.out = 80)
    # naive Oxford fever (CoP=1, delta=1)
    cur_ox <- grid_curve(ox, function(d)
      .bp(d, N50_inf, alpha_inf, 1, gamma_inf) * .bp(d, N50_fevginf, alpha_fevginf, 1, gamma_fevginf))
    # Maryland fever = phi * mixture of (P_inf * P_fev|inf), milk frame
    cur_mf <- grid_curve(md, function(d) { De <- d / delta
      pf <- function(C) .bp(De, N50_inf, alpha_inf, C, gamma_inf) * .bp(De, N50_fevginf, alpha_fevginf, C, gamma_fevginf)
      phi_md * (pi_susc * pf(CoP_susc) + (1 - pi_susc) * pf(CoP_imm)) })
    # Maryland infection = mixture of P_inf, milk frame
    cur_mi <- grid_curve(md, function(d) { De <- d / delta
      pi_susc * .bp(De, N50_inf, alpha_inf, CoP_susc, gamma_inf) +
        (1 - pi_susc) * .bp(De, N50_inf, alpha_inf, CoP_imm, gamma_inf) })
    curves <<- bind_rows(
      cur_ox %>% mutate(panel = "Oxford fever (bicarb, naive)"),
      cur_mf %>% mutate(panel = "Maryland fever (milk, mixture x phi)"),
      cur_mi %>% mutate(panel = "Maryland infection (milk, mixture)"))
    phi_hat <<- median(phi_md)
  })

  panel_of <- c(ox_fev = "Oxford fever (bicarb, naive)",
                md_fev = "Maryland fever (milk, mixture x phi)",
                md_inf = "Maryland infection (milk, mixture)",
                hornick_cond = "Maryland fever (milk, mixture x phi)",
                ox_inf = "Oxford fever (bicarb, naive)")
  pp <- as_draws_matrix(fit$draws("p_pred"))
  ci <- .wilson(obs$y, obs$n)
  pts <- obs %>% mutate(panel = unname(panel_of[likelihood_group]),
                        obs_rate = y / n, lo = ci$lo, hi = ci$hi,
                        fitted = apply(pp, 2, median)) %>%
    filter(likelihood_group != "hornick_cond")  # conditional isn't on this dose-response axis

  phi_line <- tibble(panel = "Maryland fever (milk, mixture x phi)", phi_hat = phi_hat)

  p <- ggplot(curves, aes(dose_cfu)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
    geom_line(aes(y = med), color = "steelblue", linewidth = 0.7) +
    geom_hline(data = phi_line, aes(yintercept = phi_hat), linetype = 3, color = "grey40") +
    geom_errorbar(data = pts, aes(ymin = lo, ymax = hi), width = 0.08, color = "grey50") +
    geom_point(data = pts, aes(y = obs_rate, color = study), size = 2.4) +
    geom_point(data = pts, aes(y = fitted), shape = 4, size = 2, stroke = 0.8) +  # x = Stan fitted
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    scale_x_log10(breaks = 10^(2:10),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "challenge dose (CFU)", y = "probability",
         color = "study",
         title = "Tier 1 posterior dose-response vs data",
         subtitle = paste0("line+ribbon: posterior median & 90% (population curve);  point: observed (Wilson 95% CI);  x: Stan p_pred;  dotted: phi_hat=",
                           sprintf("%.2f", phi_hat))) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")

  ggsave(outfile, p, width = 8.5, height = 10, dpi = 150)
  message("dose-response figure: ", outfile)
  invisible(p)
}

# Standalone: Rscript dose_response_curves.R  (regenerate from the saved tier1 fit)
if (sys.nframe() == 0) {
  here <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE))))
  setwd(here); source("priors.R"); source("data_prep.R")
  suppressPackageStartupMessages(library(cmdstanr))
  fit <- readRDS("results/tier1/fit.rds")
  sd  <- build_stan_data("dose_response_data.csv", load_priors("priors.yaml"),
                         tier_col = "tier1_active", prior_only = 0L)
  plot_dose_response_fit(fit, sd, "results/tier1/dose_response_fit.png")
}
