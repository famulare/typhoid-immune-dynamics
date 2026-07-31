#' Suite 3 -- model components with biological meaning and no prior figure.
#'
#'   maryland_mixture.png  the latent susceptible/immune mixture, decomposed
#'   delta_bridge.png      the milk -> bicarb dose bridge and its confound
#'   dose_cop_surface.png  the immunity mechanism as a 2D object, not a 1D slice
#'
#' NOT drawn, deliberately:
#'   eta(D)      -- eta_lo/kappa are inert at Tier 1 (prior only; every ox_inf row
#'                  has tier1_active = 0). Emitted only when group 2 is active.
#'   VE          -- the generated-quantities VE block is stale and wrong twice
#'                  over: it uses the retired dimensionless CoP_ViTT = 5.0 /
#'                  CoP_ViPS = 2.0 while the likelihood uses 152.16 / 38.11 on the
#'                  VaccZyme EU/mL scale, AND compares a control at dose 1e4
#'                  against vaccinated arms at 2e4. Fix the .stan before plotting.
#'   sigma_study -- DELETED 2026-07-31. It was declared with a prior and used in
#'     ZERO likelihood terms at every tier (no study index existed in the Stan data);
#'     the cohort random effect it was for is LOCKED-WITHDRAWN.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists("mm_pars"))   source("model_math.R")
if (!exists("mm_curve"))  source("curve_specs.R")
if (!exists(".fig_save")) source("figures_common.R")

#' The Maryland latent-immunity mixture, decomposed into its components.
#' Only the blended curve is ever shown elsewhere; pi_susc, CoP_susc and CoP_imm
#' are otherwise invisible, and CoP_imm in particular is prior-carried rather than
#' identified (priorsense prior 0.386 vs likelihood 0.071 -- no 1960s serology
#' exists). The figure says so.
plot_maryland_mixture <- function(fit, stan_data, outfile, priors = NULL,
                                  ndraw_ribbon = 1000, seed = 1, label = "Tier 1") {
  T_ref <- stan_data$T_ref %||% 38.0
  p  <- mm_thin(mm_draws(fit, T_ref = T_ref), ndraw_ribbon, seed)
  obs <- attr(stan_data, "obs")
  D  <- 10^seq(-1, 8, length.out = 90)

  comp <- function(layer) bind_rows(
    mm_quantiles(mm_md_component(D, p, layer, "susc"), D) %>%
      mutate(k = sprintf("susceptible component (CoP_susc, median %.2f)", stats::median(p$CoP_susc))),
    mm_quantiles(mm_md_component(D, p, layer, "imm"), D) %>%
      mutate(k = sprintf("immune component (CoP_imm, median %.1f) -- PRIOR-CARRIED", stats::median(p$CoP_imm))),
    mm_quantiles(mm_md_mix(D, p, layer), D) %>%
      mutate(k = sprintf("pi-weighted mixture (pi_susc median %.2f)", stats::median(p$pi_susc))))
  cur <- bind_rows(comp("inf") %>% mutate(panel = "P(infection)"),
                   comp("fev") %>% mutate(panel = "P(fever), before phi"))
  cur$panel <- factor(cur$panel, c("P(infection)", "P(fever), before phi"))

  # Gilman strata belong on their OWN component, not on the blend. Their raw dose
  # is milk, so convert to the bicarb frame with the posterior median delta.
  dmed <- stats::median(p$delta)
  gil <- obs %>% filter(gilman_stratum %in% c(1L, 2L))
  gpts <- NULL
  if (nrow(gil)) {
    gci <- .wilson(gil$y, gil$n)
    gpts <- gil %>% mutate(x = dose_cfu / dmed, rate = y / n, lo = gci$lo, hi = gci$hi,
                           panel = factor("P(fever), before phi",
                                          c("P(infection)", "P(fever), before phi")),
                           k = ifelse(gilman_stratum == 1L, "Gilman H-Ab <1:20", "Gilman H-Ab >=1:20"))
  }

  p1 <- ggplot(cur, aes(x, med, colour = k, fill = k)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) +
    {if (!is.null(gpts)) geom_errorbar(data = gpts, aes(x, ymin = lo, ymax = hi),
                                       width = 0.1, colour = "grey35", inherit.aes = FALSE)} +
    {if (!is.null(gpts)) geom_point(data = gpts, aes(x, rate, shape = k), size = 2.8,
                                    colour = "grey10", inherit.aes = FALSE)} +
    facet_wrap(~panel, ncol = 2) +
    .fig_scale_dose() + coord_cartesian(ylim = c(0, 1)) +
    scale_colour_manual(values = c("#1b7837", "#762a83", "#08519c"), name = NULL,
                        aesthetics = c("colour", "fill")) +
    scale_shape_manual(values = c(17, 15), name = "observed (bicarb frame, D/delta)") +
    labs(x = "bicarb-equivalent dose D/delta (CFU)", y = "probability",
         title = sprintf("%s: the Maryland latent-immunity mixture, decomposed", label),
         subtitle = paste0("Only the blend is shown in the other figures. The Gilman H-antibody strata are placed on THEIR OWN component, not on the blend.\n",
                           "CoP_imm is prior-carried, not identified (priorsense prior 0.386 vs likelihood 0.071): no 1960s serology exists to pin it. ",
                           "Read the immune component as an assumption.")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.box = "vertical")
  .fig_save(p1, outfile, 12, 7)
  invisible(p1)
}

#' The delta bridge: both eras on one bicarb-equivalent dose axis.
#' delta is strongly confounded with both N50s, so the point positions carry real
#' horizontal uncertainty; the inset scatter shows the confound directly.
plot_delta_bridge <- function(fit, stan_data, outfile, ndraw_ribbon = 1000,
                              seed = 1, label = "Tier 1") {
  T_ref <- stan_data$T_ref %||% 38.0
  p   <- mm_thin(mm_draws(fit, T_ref = T_ref), ndraw_ribbon, seed)
  obs <- attr(stan_data, "obs")
  D   <- 10^seq(-1, 8, length.out = 90)

  cur <- bind_rows(
    mm_quantiles(mm_p_fev(D, 1, p), D)      %>% mutate(k = "naive fever P_inf x P_fev|inf (CoP=1)"),
    mm_quantiles(mm_md_mix(D, p, "fev"), D) %>% mutate(k = "Maryland mixture fever (before phi)"),
    mm_quantiles(mm_phi_td(39.4, D, p) * mm_md_mix(D, p, "fev"), D) %>%
      mutate(k = "Maryland observed-scale fever (x phi at 39.4C)"))

  # Maryland points get a per-draw horizontal interval from delta; Oxford sits at
  # its raw dose (delta = 1 there by construction).
  fev <- obs %>% filter(likelihood_group %in% c("ox_fev", "md_fev"))
  ci  <- .wilson(fev$y, fev$n)
  dq  <- stats::quantile(p$delta, c(0.05, 0.5, 0.95))
  pts <- fev %>% mutate(
    era  = ifelse(likelihood_group == "md_fev", "Maryland (milk)", "Oxford (bicarb)"),
    x    = ifelse(likelihood_group == "md_fev", dose_cfu / dq[2], dose_cfu),
    xlo  = ifelse(likelihood_group == "md_fev", dose_cfu / dq[3], dose_cfu),
    xhi  = ifelse(likelihood_group == "md_fev", dose_cfu / dq[1], dose_cfu),
    rate = y / n, lo = ci$lo, hi = ci$hi)

  p1 <- ggplot() +
    geom_ribbon(data = cur, aes(x, ymin = lo, ymax = hi, fill = k), alpha = 0.12) +
    geom_line(data = cur, aes(x, med, colour = k), linewidth = 0.8) +
    geom_errorbar(data = pts, aes(y = rate, xmin = xlo, xmax = xhi),
                  orientation = "y", width = 0.02, colour = "grey55", linewidth = 0.4) +
    geom_errorbar(data = pts, aes(x = x, ymin = lo, ymax = hi), width = 0.08,
                  colour = "grey45", linewidth = 0.35) +
    geom_point(data = pts, aes(x, rate, shape = era), size = 2.6, colour = "grey10") +
    .fig_scale_dose() + coord_cartesian(ylim = c(0, 1)) +
    scale_colour_brewer(palette = "Dark2", name = NULL, aesthetics = c("colour", "fill")) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    labs(x = "bicarb-equivalent dose D/delta (CFU)", y = "P(fever)",
         title = sprintf("%s: the delta bridge -- both eras on one dose axis", label),
         subtitle = sprintf(paste0("Maryland points are divided by delta (median %.0fx, 90%% %.0f-%.0f); the horizontal bars ARE that uncertainty. Oxford sits at its raw dose.\n",
                                   "How much of the ~10^3 CFU gap between the eras is absorbed by delta versus by phi. delta is strongly confounded with both N50s\n",
                                   "(log10_N50_fevginf x log10_delta r = -0.78, log10_N50_inf x log10_delta r = -0.72), so its absolute value is not separately identified."),
                            dq[2], dq[1], dq[3])) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.box = "vertical")
  .fig_save(p1, outfile, 11, 7.5)
  invisible(p1)
}

#' Dose x CoP surface: the immunity mechanism (alpha / CoP^gamma) as a 2D object.
#' titre_protection.png shows only a 1D slice at D = 2e4.
plot_dose_cop_surface <- function(fit, stan_data, outfile, naive_ref = 3.7,
                                  ndraw_ribbon = 500, seed = 1, label = "Tier 1") {
  T_ref <- stan_data$T_ref %||% 38.0
  p   <- mm_thin(mm_draws(fit, T_ref = T_ref), ndraw_ribbon, seed)
  obs <- attr(stan_data, "obs")
  D   <- 10^seq(0, 8, length.out = 70)
  C   <- 10^seq(log10(0.25), log10(400), length.out = 60)

  surf <- bind_rows(lapply(c("P(infection)", "P(fever)"), function(lab) {
    f <- if (lab == "P(infection)") mm_p_inf else mm_p_fev
    bind_rows(lapply(C, function(cc)
      data.frame(dose = D, cop = cc, z = matrixStats::colMedians(f(D, cc, p)))))  %>%
      mutate(panel = lab)
  }))
  surf$eu <- surf$cop * naive_ref
  surf$panel <- factor(surf$panel, c("P(infection)", "P(fever)"))   # cascade order

  # where data actually constrain the titre axis
  ox <- obs %>% filter(likelihood_group %in% c("ox_fev", "ox_inf", "ox_inf_indiv", "ox_fevginf_indiv"))
  cop_rng <- range(ox$CoP, na.rm = TRUE)
  opts <- ox %>% transmute(dose = dose_cfu, eu = CoP * naive_ref,
                           kind = ifelse(n == 1L, "Darton individual", "study group"))

  p1 <- ggplot(surf, aes(dose, eu)) +
    geom_raster(aes(fill = z), interpolate = TRUE) +
    geom_contour(aes(z = z), breaks = seq(0.1, 0.9, by = 0.2), colour = "white",
                 linewidth = 0.35, alpha = 0.8) +
    geom_hline(yintercept = cop_rng * naive_ref, linetype = "dashed", colour = "grey15",
               linewidth = 0.4) +
    geom_point(data = opts, aes(dose, eu, shape = kind), colour = "black",
               fill = "white", size = 1.9, stroke = 0.5, inherit.aes = FALSE) +
    facet_wrap(~panel, ncol = 2) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), name = "probability") +
    scale_shape_manual(values = c(21, 24), name = "observed") +
    .fig_scale_dose(expand = c(0, 0)) +
    scale_y_log10(expand = c(0, 0)) +
    labs(x = "bicarb-equivalent dose (CFU)", y = "anti-Vi IgG (VaccZyme EU/mL)",
         title = sprintf("%s: dose x titre surface -- the immunity mechanism as a 2D object", label),
         subtitle = paste0("Immunity enters as alpha / CoP^gamma, i.e. it rescales the beta-Poisson SHAPE parameter, not the dose. titre_protection.png shows only the D=2e4 slice of this.\n",
                           "Dashed lines bracket the observed anti-Vi range; everything outside them is extrapolation, and gamma_inf vs gamma_fevginf do not separate because\n",
                           "Darton's titre spread is thin (only 12/30 subjects above the detection limit).")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.box = "vertical")
  .fig_save(p1, outfile, 12.5, 6.5)
  invisible(p1)
}

#' eta(D) shedding detection -- only meaningful once group 2 is in the likelihood.
plot_eta_detection <- function(fit, stan_data, outfile, ndraw_ribbon = 1000,
                               seed = 1, label = "Tier 1") {
  if (!any(stan_data$group == 2L)) {
    message("  [skip] eta_detection.png: group 2 (ox_inf) not active in this run; ",
            "eta_lo/kappa are prior-only and the figure would show the prior, not a fit.")
    return(invisible(NULL))
  }
  p <- mm_thin(mm_draws(fit, T_ref = stan_data$T_ref %||% 38.0), ndraw_ribbon, seed)
  D <- 10^seq(1, 6, length.out = 80)
  cur <- bind_rows(mm_quantiles(mm_eta(D, p), D) %>% mutate(k = "eta(D) detection probability"),
                   mm_quantiles(mm_eta(D, p) * mm_p_inf(D, 1, p), D) %>%
                     mutate(k = "eta(D) x P_inf(D) = observed shedding"),
                   mm_quantiles(mm_p_inf(D, 1, p), D) %>% mutate(k = "P_inf(D) true infection"))
  p1 <- ggplot(cur, aes(x, med, colour = k, fill = k)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) + .fig_scale_dose() + coord_cartesian(ylim = c(0, 1)) +
    scale_colour_brewer(palette = "Set1", name = NULL, aesthetics = c("colour", "fill")) +
    labs(x = "challenge dose (CFU)", y = "probability",
         title = sprintf("%s: eta(D) shedding-detection correction", label),
         subtitle = "Shedding is missed only if it would have started after antibiotic treatment, so detection falls as dose rises.") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  .fig_save(p1, outfile, 9, 6)
  invisible(p1)
}
