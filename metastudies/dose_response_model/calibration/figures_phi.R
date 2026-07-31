#' Suite 2 -- the fever-severity / definition map phi(T, D).
#'
#' This is an entire fitted sub-likelihood with, before this file, ZERO visual
#' representation anywhere in the repo: phi0_a / phi0_b are identified solely by
#' the Darton placebo temperature ladder (16/10/8 of 20 TD+ subjects crossing
#' >=38.0 / 38.5 / 39.0 C), and phi(T,D) then multiplies every Maryland fever
#' observation. It is also the component whose mis-specification caused the
#' model's largest historical failure -- the phi-cap catastrophe, where a
#' dose-CONSTANT phi capped Maryland fitted fever at 0.25-0.30 against Hornick
#' high-dose data of 0.89-0.95 (tier1_lab_notebook.md:24-56, 257-278).
#'
#' Biologically, phi0(T) is the survival function of PEAK FEVER SEVERITY among
#' clinically-diagnosed cases: different studies observe different thresholds on
#' the same underlying temperature distribution.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists("mm_pars"))   source("model_math.R")
if (!exists("mm_curve"))  source("curve_specs.R")
if (!exists(".fig_save")) source("figures_common.R")

#' @param fit cmdstanr fit; @param stan_data list carrying the Darton ladder
plot_phi_severity <- function(fit, stan_data, outfile, n_spaghetti = 100,
                              ndraw_ribbon = 1000, seed = 1, label = "Tier 1") {
  T_ref <- if (!is.null(stan_data$T_ref)) stan_data$T_ref else 38.0
  p_r <- mm_thin(mm_draws(fit, T_ref = T_ref), ndraw_ribbon, seed)
  p_s <- mm_thin(p_r, n_spaghetti, seed)

  Tg  <- seq(37.0, 41.0, length.out = 120)
  # Ladder maximum: everything above it is extrapolation, and Hornick's 39.4 sits
  # outside it. The ladder is the ONLY thing identifying phi0.
  T_max_data <- max(stan_data$ladder_T)
  studies <- data.frame(T = c(T_ref, 38.3, 39.4),
                        lab = c(sprintf("T_ref %.1f (Oxford composite)", T_ref),
                                "38.3 (Levine, Gilman)", "39.4 (Hornick)"))

  # ---- A: phi0(T) with the Darton ladder ------------------------------------
  qa <- mm_quantiles(mm_phi0(Tg, p_r), Tg)
  sa <- mm_spaghetti(mm_phi0(Tg, p_s), Tg, p_s)
  lad_ci <- .wilson(stan_data$ladder_count, rep(stan_data$ladder_N, stan_data$N_ladder))
  lad <- data.frame(T = stan_data$ladder_T,
                    rate = stan_data$ladder_count / stan_data$ladder_N,
                    lo = lad_ci$lo, hi = lad_ci$hi)

  pa <- ggplot() +
    annotate("rect", xmin = T_max_data, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "grey70", alpha = 0.22) +
    geom_ribbon(data = qa, aes(x, ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.10) +
    geom_line(data = sa, aes(x, y, group = .draw), colour = "steelblue",
              linewidth = 0.18, alpha = .fig_spag_alpha(p_s$.ndraws)) +
    geom_line(data = qa, aes(x, med), colour = "steelblue4", linewidth = 0.9) +
    geom_vline(data = studies, aes(xintercept = T), linetype = "dashed",
               colour = "grey30", linewidth = 0.35) +
    geom_text(data = studies, aes(x = T, y = 0.03, label = lab), angle = 90,
              hjust = 0, vjust = -0.35, size = 2.5, colour = "grey25") +
    geom_errorbar(data = lad, aes(T, ymin = lo, ymax = hi), width = 0.05, colour = "grey35") +
    geom_point(data = lad, aes(T, rate), colour = "firebrick", size = 2.6) +
    annotate("text", x = T_max_data + 0.05, y = 0.95, hjust = 0, size = 2.6,
             colour = "grey25", label = "extrapolation:\nno ladder data") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "strict fever threshold T (deg C)", y = "phi0(T)",
         title = "A. phi0(T): low-dose definition sensitivity",
         subtitle = sprintf("red: Darton ladder %s of %d TD+ (Wilson 95%%)\nthe only data identifying phi0",
                            paste(stan_data$ladder_count, collapse = "/"), stan_data$ladder_N)) +
    theme_minimal(base_size = 10)

  # ---- B: phi(T,D) vs T, one curve per dose ---------------------------------
  doses <- 10^c(0, 2, 4, 6, 8)
  qb <- bind_rows(lapply(doses, function(d)
    mm_quantiles(mm_phi_td(Tg, rep(d, length(Tg)), p_r), Tg) %>%
      mutate(dose = sprintf("10^%d", round(log10(d))))))
  pb <- ggplot(qb, aes(x, med, colour = dose)) +
    annotate("rect", xmin = T_max_data, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "grey70", alpha = 0.22) +
    geom_line(linewidth = 0.8) +
    geom_vline(data = studies, aes(xintercept = T), linetype = "dashed",
               colour = "grey30", linewidth = 0.35, inherit.aes = FALSE) +
    scale_colour_viridis_d(option = "C", end = 0.85,
                           name = "bicarb-equivalent\ndose (CFU)") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "strict fever threshold T (deg C)", y = "phi(T, D)",
         title = "B. The dose lift flattens phi in T",
         subtitle = "phi(T,D) = phi0(T) + (1-phi0(T))*P_fev_naive(D)\nbeta_phi PINNED to 1: the lift shape is inherited, not fitted") +
    theme_minimal(base_size = 10)

  # ---- C: phi(T,D) vs dose at the study thresholds --------------------------
  Dg <- 10^seq(-1, 8, length.out = 80)
  qc <- bind_rows(lapply(seq_len(nrow(studies)), function(i)
    mm_quantiles(mm_phi_td(studies$T[i], Dg, p_r), Dg) %>%
      mutate(thr = sprintf("%.1f C", studies$T[i]))))
  sc <- mm_spaghetti(mm_phi_td(39.4, Dg, p_s), Dg, p_s)
  pc <- ggplot() +
    geom_line(data = sc, aes(x, y, group = .draw), colour = "grey45",
              linewidth = 0.18, alpha = .fig_spag_alpha(p_s$.ndraws)) +
    geom_ribbon(data = qc, aes(x, ymin = lo, ymax = hi, fill = thr), alpha = 0.15) +
    geom_line(data = qc, aes(x, med, colour = thr), linewidth = 0.8) +
    .fig_scale_dose() + coord_cartesian(ylim = c(0, 1)) +
    scale_colour_brewer(palette = "Dark2", name = "threshold") +
    scale_fill_brewer(palette = "Dark2", name = "threshold") +
    labs(x = "bicarb-equivalent dose D/delta (CFU)", y = "phi(T, D)",
         title = "C. What multiplies the Maryland fever likelihood",
         subtitle = "grey spaghetti: posterior draws at the Hornick 39.4 C threshold") +
    theme_minimal(base_size = 10)

  # ---- D: implied peak-fever severity density -------------------------------
  # phi(T,D) is a survival function in T, so -d(phi)/dT is the implied density of
  # PEAK TEMPERATURE among diagnosed cases. At low dose this is the analytic
  # logistic density phi0_b*phi0*(1-phi0).
  #
  # FINDING (2026-07-31): the high-dose density COLLAPSES rather than shifting up.
  # Because the lift term (1-phi0(T))*P_fev_naive(D) carries no T dependence,
  # phi -> 1 UNIFORMLY in T as dose grows, so the model cannot express "more dose
  # means hotter fevers" -- it expresses "more dose means threshold choice stops
  # mattering". At Hornick's top dose (1e9 milk ~ 8.5e6 bicarb-eq) that implies
  # phi(41C) = 0.95 and phi(42C) = 0.95, i.e. ~95% of diagnosed cases exceeding
  # 41C. That is a consequence of pinning beta_phi = 1, not a fitted claim.
  qd_lo <- mm_quantiles(mm_phi0_density(Tg, p_r), Tg) %>% mutate(k = "low dose (analytic phi0 density)")
  hi_d <- 1e6
  ph_hi <- mm_phi_td(Tg, rep(hi_d, length(Tg)), p_r)
  dens_hi <- -t(apply(ph_hi, 1, function(v) c(diff(v), NA) / c(diff(Tg), NA)))
  qd_hi <- mm_quantiles(dens_hi[, -length(Tg), drop = FALSE], Tg[-length(Tg)]) %>%
    mutate(k = sprintf("D = 10^%d bicarb-eq", round(log10(hi_d))))
  qd <- bind_rows(qd_lo, qd_hi)
  phi41 <- stats::median(mm_phi_td(41, hi_d, p_r))
  pd <- ggplot(qd, aes(x, med, colour = k, fill = k)) +
    annotate("rect", xmin = T_max_data, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "grey70", alpha = 0.22) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(data = studies, aes(xintercept = T), linetype = "dashed",
               colour = "grey30", linewidth = 0.35, inherit.aes = FALSE) +
    annotate("text", x = 39.6, y = Inf, hjust = 0, vjust = 1.3, size = 2.7, colour = "#b2182b",
             label = sprintf("high-dose density COLLAPSES, it does not shift up:\nthe lift has no T dependence, so phi -> 1 uniformly.\nAt D=10^6 the model implies phi(41C) = %.2f.\nArtifact of beta_phi pinned to 1.", phi41)) +
    scale_colour_manual(values = c("#2166ac", "#b2182b"), name = NULL, aesthetics = c("colour", "fill")) +
    labs(x = "peak temperature T (deg C)", y = "density  -d(phi)/dT",
         title = "D. Implied peak-fever severity distribution",
         subtitle = "phi is a survival function in T, so -d(phi)/dT is the implied\npeak-temperature density among diagnosed cases") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")

  # ---- assemble as one faceted-looking column of four ------------------------
  p <- .stack4(pa, pb, pc, pd,
    title = sprintf("%s: fever severity and the study-definition map phi(T, D)", label))
  .fig_save(p, outfile, 14, 13)
  invisible(p)
}

#' Stack four independent plots without adding a patchwork dependency.
#' gridExtra is already a dependency via diagnostics.R::.save_grid().
.stack4 <- function(a, b, c, d, title = NULL) {
  gridExtra::arrangeGrob(a, b, c, d, ncol = 2,
    top = grid::textGrob(title %||% "", gp = grid::gpar(fontsize = 13), hjust = 0.5))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
