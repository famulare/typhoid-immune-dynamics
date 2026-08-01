#' Bespoke model diagnostic: posterior dose-response curves and CoP counterfactuals.
#'
#' This is the model-specific curve view (in addition to the standard
#' bayesplot/posterior battery in diagnostics.R). It plots the population
#' dose-response curve (median + 90% ribbon) and fixed-CoP counterfactuals over a
#' dose grid. Observed overlays are intentionally omitted: individual endpoints
#' stack at one dose and are shown in the separate grouping/titre figures.
#'
#' The curves are computed in R from posterior parameter draws via model_math.R,
#' which is gated against Stan by test_obs_prob_parity.R. The inference-side PPC
#' remains in diagnostics.R.
#'
#' 2026-07-31: the local `.bp()` beta-Poisson mirror was DELETED. It was a third,
#' ungated transcription of the model math and every panel re-derived the Maryland
#' fever product inline. All model algebra now lives in model_math.R.

suppressPackageStartupMessages({
  library(posterior); library(dplyr); library(ggplot2)
})
if (!exists("mm_pars")) source("model_math.R")
if (!exists(".wilson")) source("utils.R")
if (!exists(".fig_save")) source("figures_common.R")

#' @param fit cmdstanr fit; @param stan_data list with attr "obs"; @param outfile png path
plot_dose_response_fit <- function(fit, stan_data, outfile) {
  T_ref_val <- if (!is.null(stan_data$T_ref)) stan_data$T_ref else 38.0
  T_curve   <- 39.4   # draw the Maryland fever curve at the Hornick threshold (spans dose)
  pars <- mm_draws(fit, T_ref = T_ref_val)

  gc_ <- function(mat, doses) {
    q <- mm_quantiles(mat, doses)
    tibble(dose_cfu = q$x, lo = q$lo, med = q$med, hi = q$hi)
  }
  ox <- 10^seq(2.3, 4.7, length.out = 60)
  md <- 10^seq(2.7, 9.7, length.out = 80)
  De <- mm_by_grid(md, pars$.ndraws, length(md)) / mm_by_draw(pars$delta, pars$.ndraws, length(md))
  fixed_cop <- c(100, 1000)

  curves <- bind_rows(
    # naive Oxford fever (CoP=1, delta=1)
    gc_(mm_p_fev(ox, 1, pars), ox) %>% mutate(panel = "Oxford fever (bicarb, naive)"),
    # Maryland fever = phi(T,D) * mixture of (P_inf * P_fev|inf), milk frame
    gc_(mm_phi_td(T_curve, De, pars) * mm_md_mix(De, pars, "fev"), md) %>%
      mutate(panel = "Maryland fever (milk, mixture x phi(39.4,D))"),
    # Maryland infection = mixture of P_inf, milk frame
    gc_(mm_md_mix(De, pars, "inf"), md) %>%
      mutate(panel = "Maryland infection (milk, mixture)"),
    # phi is defined from the naive (CoP=1) fever curve and has no CoP covariate.
    gc_(mm_phi_td(T_curve, De, pars), md) %>%
      mutate(panel = "Fever-threshold sensitivity (phi(39.4,D), CoP=1)")
  )

  # These are fixed-CoP counterfactuals, not additional fitted population curves.
  # phi(T,D) remains the model's naive-CoP definition-sensitivity map; immunity
  # enters here through the fixed-CoP infection/fever kernels upstream of phi.
  cop_curves <- bind_rows(lapply(fixed_cop, function(cop) {
    series <- sprintf("fixed CoP = %g", cop)
    bind_rows(
      gc_(mm_p_fev(ox, cop, pars), ox) %>%
        mutate(panel = "Oxford fever (bicarb, naive)", series = series),
      gc_(mm_phi_td(T_curve, De, pars) * mm_p_fev(De, cop, pars), md) %>%
        mutate(panel = "Maryland fever (milk, mixture x phi(39.4,D))", series = series),
      gc_(mm_p_inf(De, cop, pars), md) %>%
        mutate(panel = "Maryland infection (milk, mixture)", series = series)
    )
  }))
  cop_colors <- setNames(scales::hue_pal()(length(fixed_cop)),
                         sprintf("fixed CoP = %g", fixed_cop))
  panel_levels <- c("Maryland fever (milk, mixture x phi(39.4,D))",
                    "Maryland infection (milk, mixture)",
                    "Oxford fever (bicarb, naive)",
                    "Fever-threshold sensitivity (phi(39.4,D), CoP=1)")
  curves$panel <- factor(curves$panel, levels = panel_levels)
  cop_curves$panel <- factor(cop_curves$panel, levels = panel_levels)

  p <- ggplot(curves, aes(dose_cfu)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
    geom_ribbon(data = cop_curves, aes(ymin = lo, ymax = hi, fill = series),
                alpha = 0.12, color = NA) +
    geom_line(aes(y = med), color = "steelblue", linewidth = 0.7) +
    geom_line(data = cop_curves, aes(y = med, color = series), linewidth = 0.7) +
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    scale_x_log10(breaks = 10^(2:10),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "challenge dose (CFU)", y = "probability",
         color = "counterfactual",
         title = "Posterior dose-response curves",
         subtitle = "line+ribbon: posterior median & 90% (population curve); colored lines+ribbons: fixed CoP counterfactuals") +
    scale_color_manual(values = cop_colors) +
    scale_fill_manual(values = cop_colors, guide = "none") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")

  ggsave(outfile, p, width = 8.5, height = 13, dpi = 150)
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
  obs  <- attr(stan_data, "obs")
  pars <- mm_draws(fit, T_ref = if (!is.null(stan_data$T_ref)) stan_data$T_ref else 38.0)
  cop  <- 10^seq(log10(0.9), log10(200), length.out = 80)
  gc_ <- function(mat) {
    q <- mm_quantiles(mat, cop * naive_ref)
    tibble(eu = q$x, lo = q$lo, med = q$med, hi = q$hi)
  }
  curves <- bind_rows(
    gc_(mm_p_inf(D_ref, cop, pars))     %>% mutate(panel = "P(infection)"),
    gc_(mm_p_fevginf(D_ref, cop, pars)) %>% mutate(panel = "P(fever | infection)"),
    gc_(mm_p_fev(D_ref, cop, pars))     %>% mutate(panel = "P(fever) composite"))
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

#' CoP -> protection across fixed raw milk doses.
#'
#' This is the CoP-axis analogue of the dose-response figures. Each dose is fixed
#' in the milk frame, so each posterior draw uses its own bicarb-equivalent dose
#' `dose_milk / delta`. Each panel overlays the latent cascade quantities from
#' model_math.R: P_inf, P_fev|inf, and P_fev = P_inf x P_fev|inf. No Maryland
#' fever-threshold definition map (phi) is applied here.
#' @param dose_milk raw challenge doses in the milk frame (CFU)
#' @param cop_range positive CoP range to display
plot_cop_response_milk <- function(fit, stan_data, outfile,
                                   dose_milk = c(1e2, 1e4, 1e7),
                                   cop_range = c(0.25, 400), ngrid = 80,
                                   label = "Tier 1") {
  if (!length(dose_milk) || any(!is.finite(dose_milk)) || any(dose_milk <= 0))
    stop("dose_milk must contain positive finite values", call. = FALSE)
  if (length(cop_range) != 2L || any(!is.finite(cop_range)) || any(cop_range <= 0) ||
      cop_range[1] >= cop_range[2])
    stop("cop_range must be two increasing positive finite values", call. = FALSE)
  if (length(ngrid) != 1L || !is.finite(ngrid) || ngrid < 2 || ngrid != as.integer(ngrid))
    stop("ngrid must be one integer >= 2", call. = FALSE)

  T_ref_val <- if (!is.null(stan_data$T_ref)) stan_data$T_ref else 38.0
  pars <- mm_draws(fit, T_ref = T_ref_val)
  cop <- 10^seq(log10(cop_range[1]), log10(cop_range[2]), length.out = ngrid)
  dose_label <- function(d) {
    e <- log10(d)
    if (abs(e - round(e)) < 1e-8) sprintf("10^%d CFU milk", round(e))
    else paste0(format(d, scientific = TRUE, trim = TRUE), " CFU milk")
  }
  dose_labels <- vapply(dose_milk, dose_label, character(1))
  curves <- bind_rows(lapply(seq_along(dose_milk), function(i) {
    D_eff <- mm_by_grid(dose_milk[i], pars$.ndraws, length(cop)) /
      mm_by_draw(pars$delta, pars$.ndraws, length(cop))
    bind_rows(
      mm_quantiles(mm_p_inf(D_eff, cop, pars), cop) %>%
        mutate(endpoint = "P_inf"),
      mm_quantiles(mm_p_fevginf(D_eff, cop, pars), cop) %>%
        mutate(endpoint = "P_fev|inf"),
      mm_quantiles(mm_p_fev(D_eff, cop, pars), cop) %>%
        mutate(endpoint = "P_fev")
    ) %>% mutate(dose = dose_labels[i])
  }))
  endpoint_levels <- c("P_inf", "P_fev|inf", "P_fev")
  endpoint_colors <- c(P_inf = "#1b7837", `P_fev|inf` = "#762a83", P_fev = "#08519c")
  curves$endpoint <- factor(curves$endpoint, levels = endpoint_levels)
  curves$dose <- factor(curves$dose, levels = dose_labels)
  delta_q <- stats::quantile(pars$delta, c(0.05, 0.5, 0.95))

  p <- ggplot(curves, aes(x, med, colour = endpoint, fill = endpoint, group = endpoint)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.75) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    facet_wrap(~dose, nrow = 1) +
    scale_x_log10(breaks = c(0.25, 1, 10, 100, 400),
                  labels = c("0.25", "1", "10", "100", "400")) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = 0.03)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_colour_manual(values = endpoint_colors,
                        labels = c(P_inf = "P_inf", `P_fev|inf` = "P_fev|inf", P_fev = "P_fev"),
                        name = NULL) +
    scale_fill_manual(values = endpoint_colors, guide = "none") +
    labs(x = "CoP (anti-Vi titre / naive)", y = "probability",
         title = sprintf("%s: CoP response across milk doses", label),
         subtitle = paste0(
           "line+ribbon: posterior median & 90%; dashed line: naive CoP = 1.\n",
           "Each panel fixes raw milk dose; each draw uses D_eff = dose / delta. ",
           sprintf("delta median %.0fx (90%% %.0f-%.0f).\n", delta_q[2], delta_q[1], delta_q[3]),
           "P_fev is the latent cascade P_inf x P_fev|inf; no Maryland fever-threshold map phi is applied.")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", strip.text = element_text(size = 10),
          plot.subtitle = element_text(size = 8.5))

  .fig_save(p, outfile, w = 12, h = 5.5)
  invisible(p)
}

# Standalone: Rscript dose_response_curves.R  (regenerate from the saved tier1 fit)
if (sys.nframe() == 0) {
  setwd(calib_dir())
  source("priors.R"); source("data_prep.R"); source("tier_specs.R"); source("provenance.R")
  suppressPackageStartupMessages(library(cmdstanr))
  args <- commandArgs(trailingOnly = TRUE)
  run_dir <- if (length(args)) args[1] else tier_out_dir(tier_spec("t1-indiv"))
  fit <- readRDS(file.path(run_dir, "fit.rds"))
  sd  <- resolve_run_stan_data(run_dir)
  plot_dose_response_fit(fit, sd, file.path(run_dir, "dose_response_fit.png"))
  plot_titre_protection(fit, sd, file.path(run_dir, "titre_protection.png"))
  plot_cop_response_milk(fit, sd, file.path(run_dir, "cop_response_milk_doses.png"),
                        label = basename(run_dir))
}
