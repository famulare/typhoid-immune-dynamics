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
                                 "pi_susc", "CoP_imm", "CoP_susc",
                                 "phi0_a", "phi0_b")))
  T_ref_val <- if (!is.null(stan_data$T_ref)) stan_data$T_ref else 38.0
  T_curve   <- 39.4   # draw the Maryland fever curve at the Hornick threshold (spans dose)

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
    # Maryland fever = phi(T,D) * mixture of (P_inf * P_fev|inf), milk frame.
    # phi(T,D) drawn at Hornick threshold 39.4: phi0(T) + (1-phi0)*P_fev_naive(De)^beta.
    cur_mf <- grid_curve(md, function(d) { De <- d / delta
      pf <- function(C) .bp(De, N50_inf, alpha_inf, C, gamma_inf) * .bp(De, N50_fevginf, alpha_fevginf, C, gamma_fevginf)
      phi0 <- plogis(phi0_a - phi0_b * (T_curve - T_ref_val))
      p_fev_naive <- .bp(De, N50_inf, alpha_inf, 1, gamma_inf) * .bp(De, N50_fevginf, alpha_fevginf, 1, gamma_fevginf)
      phi <- phi0 + (1 - phi0) * p_fev_naive   # beta_phi pinned = 1
      phi * (pi_susc * pf(CoP_susc) + (1 - pi_susc) * pf(CoP_imm)) })
    # Maryland infection = mixture of P_inf, milk frame
    cur_mi <- grid_curve(md, function(d) { De <- d / delta
      pi_susc * .bp(De, N50_inf, alpha_inf, CoP_susc, gamma_inf) +
        (1 - pi_susc) * .bp(De, N50_inf, alpha_inf, CoP_imm, gamma_inf) })
    curves <<- bind_rows(
      cur_ox %>% mutate(panel = "Oxford fever (bicarb, naive)"),
      cur_mf %>% mutate(panel = "Maryland fever (milk, mixture x phi(39.4,D))"),
      cur_mi %>% mutate(panel = "Maryland infection (milk, mixture)"))
  })

  panel_of <- c(ox_fev = "Oxford fever (bicarb, naive)",
                md_fev = "Maryland fever (milk, mixture x phi(39.4,D))",
                md_inf = "Maryland infection (milk, mixture)",
                hornick_cond = "Maryland fever (milk, mixture x phi(39.4,D))",
                ox_inf = "Oxford fever (bicarb, naive)")
  pp <- as_draws_matrix(fit$draws("p_pred"))
  ci <- .wilson(obs$y, obs$n)
  pts <- obs %>% mutate(panel = unname(panel_of[likelihood_group]),
                        obs_rate = y / n, lo = ci$lo, hi = ci$hi,
                        fitted = apply(pp, 2, median)) %>%
    # conditional + individual single-dose endpoints aren't on the dose-response axis
    # (the Darton individuals belong on the titre panel; see tier1.5 plots task)
    filter(!likelihood_group %in% c("hornick_cond", "ox_inf_indiv", "ox_fevginf_indiv"))

  p <- ggplot(curves, aes(dose_cfu)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
    geom_line(aes(y = med), color = "steelblue", linewidth = 0.7) +
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
         subtitle = "line+ribbon: posterior median & 90% (population curve);  point: observed (Wilson 95% CI);  x: Stan p_pred") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")

  ggsave(outfile, p, width = 8.5, height = 10, dpi = 150)
  message("dose-response figure: ", outfile)
  invisible(p)
}

#' Titre -> protection (CoP-axis) figure: the view that shows the immunity slope,
#' which the dose-axis panels hide (they stack the individual Darton subjects at one
#' dose). Three facets at a fixed Oxford dose (~2e4): P(infection), P(fever|infection),
#' and composite P(fever), each the posterior CoP^gamma curve (median + 90%) vs anti-Vi
#' EU/mL. Overlays the Darton individual endpoints (jittered 0/1: infection on the
#' P(inf) facet, fever|inf on the P(fev|inf) facet) and the Jin vaccine-group fever
#' points (Wilson 95%) on the composite facet. Makes the thin Darton titre range and
#' the Jin high-titre anchors legible — the setup for +Jin-digitize.
#' @param D_ref Oxford challenge dose to evaluate the curves at (Darton 18200 ~ Jin 2e4).
plot_titre_protection <- function(fit, stan_data, outfile, D_ref = 2e4,
                                  naive_ref = 3.7) {
  obs <- attr(stan_data, "obs")
  dr  <- as_draws_df(fit$draws(c("N50_inf", "N50_fevginf", "alpha_inf",
                                 "alpha_fevginf", "gamma_inf", "gamma_fevginf")))
  cop <- 10^seq(log10(0.9), log10(200), length.out = 80)
  gc <- function(fn) {
    q <- vapply(cop, function(c) quantile(fn(c), c(.05, .5, .95)), numeric(3))
    tibble(eu = cop * naive_ref, lo = q[1, ], med = q[2, ], hi = q[3, ])
  }
  with(dr, {
    ci <- gc(function(C) .bp(D_ref, N50_inf, alpha_inf, C, gamma_inf))
    cg <- gc(function(C) .bp(D_ref, N50_fevginf, alpha_fevginf, C, gamma_fevginf))
    cf <- gc(function(C) .bp(D_ref, N50_inf, alpha_inf, C, gamma_inf) *
                         .bp(D_ref, N50_fevginf, alpha_fevginf, C, gamma_fevginf))
    curves <<- bind_rows(ci %>% mutate(panel = "P(infection)"),
                         cg %>% mutate(panel = "P(fever | infection)"),
                         cf %>% mutate(panel = "P(fever) composite"))
  })
  panel_lv <- c("P(infection)", "P(fever | infection)", "P(fever) composite")
  curves$panel <- factor(curves$panel, panel_lv)

  ind <- obs %>%
    filter(likelihood_group %in% c("ox_inf_indiv", "ox_fevginf_indiv")) %>%
    mutate(eu = CoP * naive_ref,
           panel = ifelse(likelihood_group == "ox_inf_indiv",
                          "P(infection)", "P(fever | infection)"),
           yj = y + runif(n(), -0.03, 0.03))
  ind$panel <- factor(ind$panel, panel_lv)

  jin <- obs %>% filter(study == "Jin", likelihood_group == "ox_fev")
  jci <- .wilson(jin$y, jin$n)
  jin <- jin %>% mutate(eu = CoP * naive_ref, rate = y / n,
                        lo = jci$lo, hi = jci$hi,
                        panel = factor("P(fever) composite", panel_lv))

  p <- ggplot(curves, aes(eu)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick", alpha = 0.18) +
    geom_line(aes(y = med), color = "firebrick", linewidth = 0.7) +
    geom_point(data = ind, aes(y = yj), alpha = 0.4, size = 1.3,
               color = "grey30") +
    geom_errorbar(data = jin, aes(ymin = lo, ymax = hi), width = 0.05,
                  color = "steelblue") +
    geom_point(data = jin, aes(y = rate), color = "steelblue", size = 2.6) +
    facet_wrap(~panel, ncol = 1) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "anti-Vi IgG (VaccZyme EU/mL, log scale)", y = "probability",
         title = "Titre -> protection: immunity slope (CoP^gamma) at Oxford dose",
         subtitle = paste0("red: posterior median & 90% at D=", format(D_ref, scientific = TRUE),
                           ";  grey: Darton individuals (jittered 0/1);  blue: Jin vaccine groups (Wilson 95%)")) +
    theme_minimal(base_size = 11)

  ggsave(outfile, p, width = 7.5, height = 9, dpi = 150)
  message("titre-protection figure: ", outfile)
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
  plot_titre_protection(fit, sd, "results/tier1/titre_protection.png")
}
