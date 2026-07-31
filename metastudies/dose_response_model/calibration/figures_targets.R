#' Calibration targets: observed vs posterior, one dot-and-whisker pair per target.
#'
#' Complements ppc.png rather than duplicating it. ppc.png is an observed-vs-fitted
#' SCATTER against the 1:1 line, with a posterior interval on the y-axis and NO
#' interval on the observed rate. This figure puts both uncertainties on a common
#' probability axis, one row per calibration target, so you can read directly:
#'   - which targets the model misses, and in which direction
#'   - whether a miss is large relative to the DATA's own uncertainty (small n) or
#'     relative to the POSTERIOR's (tight parameter constraint)
#'   - which targets carry real weight at all (n is printed)
#'
#' Arm-level pooling of the individual (n=1) rows is delegated to
#' diagnostics.R::.ppc_rows(), so the rows here are exactly the rows ppc.png plots.
#' Wilson intervals come from dose_response_curves.R::.wilson().

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists(".fig_save"))  source("figures_common.R")
if (!exists(".wilson"))    source("dose_response_curves.R")
if (!exists(".ppc_rows"))  source("diagnostics.R")

#' @param p_star optional ndraws x N_obs matrix of the latent per-observation
#'   probability p* recovered from the beta-binomial (Step 2). When supplied it is
#'   drawn as a third series, which is what makes the overdispersion auditable per
#'   row instead of a single opaque rho. Columns must align with `obs`.
plot_calibration_targets <- function(fit, stan_data, outfile, p_star = NULL,
                                     label = "Tier 1") {
  obs <- attr(stan_data, "obs")
  pp  <- posterior::as_draws_matrix(fit$draws("p_pred"))
  if (ncol(pp) != nrow(obs))
    stop("plot_calibration_targets: obs (", nrow(obs), ") does not match p_pred (",
         ncol(pp), ") -- wrong stan_data for this fit?", call. = FALSE)

  df <- .ppc_rows(obs, pp) %>% filter(level == "arm")
  ci <- .wilson(df$y, df$n)
  df <- df %>% mutate(obs_lo = ci$lo, obs_hi = ci$hi,
                      lab = sprintf("%s  (%d/%d)", obs_id, y, n))

  series <- bind_rows(
    df %>% transmute(lab, likelihood_group, dose_cfu,
                     mid = obs_rate, lo = obs_lo, hi = obs_hi,
                     what = "observed (Wilson 95%)"),
    df %>% transmute(lab, likelihood_group, dose_cfu,
                     mid = fit_med, lo = fit_lo, hi = fit_hi,
                     what = "posterior p_pred (median, 90%)"))

  if (!is.null(p_star)) {
    if (ncol(p_star) != nrow(obs))
      stop("plot_calibration_targets: p_star does not align with obs", call. = FALSE)
    q <- t(apply(p_star, 2, stats::quantile, probs = c(0.05, 0.5, 0.95)))
    ps <- obs %>% mutate(lo = q[, 1], mid = q[, 2], hi = q[, 3]) %>%
      filter(n > 1) %>%
      transmute(lab = sprintf("%s  (%d/%d)", obs_id, y, n), likelihood_group,
                dose_cfu, mid, lo, hi, what = "latent p* (beta-binomial)")
    series <- bind_rows(series, ps %>% filter(lab %in% df$lab))
  }

  # order targets within each group by dose, then by label, and keep that order
  # top-to-bottom on the y axis
  ord <- df %>% arrange(likelihood_group, dose_cfu, lab) %>% pull(lab)
  series <- series %>% mutate(lab = factor(lab, levels = rev(unique(ord))))

  lv <- c("observed (Wilson 95%)", "posterior p_pred (median, 90%)",
          "latent p* (beta-binomial)")
  series$what <- factor(series$what, levels = intersect(lv, unique(series$what)))

  n_row <- nlevels(series$lab)
  # Non-overlap of the 90% posterior and the 95% Wilson interval is a crude but
  # useful headline: it counts targets the model cannot reach even allowing for the
  # data's own sampling uncertainty. Under a well-calibrated model with these n it
  # should be near zero.
  miss <- df %>% filter(fit_lo > obs_hi | fit_hi < obs_lo)
  sub_txt <- paste0(
    "One dot-and-whisker pair per target on a common axis. Complements ppc.png (same information as a\n",
    "scatter against the 1:1 line, but with no interval on the observed rate).  A miss is only meaningful\n",
    "relative to BOTH whiskers: a wide black interval means small n; a narrow blue interval means the\n",
    "parameters are pinned by other targets.  Individual (n=1) rows are pooled to arm level by the same\n",
    ".ppc_rows() code ppc.png uses.\n\n",
    sprintf("Non-overlapping (90%% posterior vs 95%% Wilson): %d of %d targets%s",
            nrow(miss), nrow(df),
            if (nrow(miss)) paste0(" -- ", paste(miss$obs_id, collapse = ", ")) else ""))
  p <- ggplot(series, aes(mid, lab, colour = what)) +
    geom_linerange(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = 0.6), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.6), size = 1.9) +
    facet_grid(likelihood_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = 0.01)) +
    scale_colour_manual(values = c("observed (Wilson 95%)" = "grey15",
                                   "posterior p_pred (median, 90%)" = "steelblue3",
                                   "latent p* (beta-binomial)" = "#d95f02"),
                        name = NULL) +
    labs(x = "probability", y = NULL,
         title = sprintf("%s: every calibration target, observed vs posterior", label),
         subtitle = sub_txt) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_line(colour = "grey92"),
          panel.grid.minor.x = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 8),
          strip.placement = "outside",
          axis.text.y = element_text(size = 7.5),
          plot.subtitle = element_text(size = 7.5))

  .fig_save(p, outfile, w = 10, h = max(5, 1.9 + 0.30 * n_row))
  invisible(p)
}
