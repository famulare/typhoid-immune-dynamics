#' Suite 1 -- the permutation grid.
#'
#' Five rows x N comparable-dataset-sub-groupings. Where dose_response_fit.png
#' shows three hand-picked panels at three of the model's covariate settings, this
#' shows EVERY grouping the data define, with posterior-draw spaghetti so the joint
#' parameter uncertainty and the curve-shape correlation are visible rather than
#' summarized away by a pointwise ribbon.
#'
#' Rows: P(inf|D), P(fev|D), P(fev|inf,D), phi(T,D), and the cascade factorization
#' (rows 1-4 are the factors; row 5 overlays them down to the observed quantity).
#'
#' All curve math comes from model_math.R via curve_specs.R::mm_curve(); this file
#' arranges, it does not derive. curve_parity_check() re-asserts at figure time that
#' the math agrees with THIS fit's own p_pred.

suppressPackageStartupMessages({
  library(posterior); library(dplyr); library(ggplot2)
})
if (!exists("mm_pars"))   source("model_math.R")
if (!exists("mm_curve"))  source("curve_specs.R")
if (!exists(".fig_save")) source("figures_common.R")

#' Assert the figure math reproduces this fit's p_pred at its own posterior draws.
#'
#' Strictly stronger than the offline gate in one respect: test_obs_prob_parity.R
#' checks a few hand-chosen parameter vectors, this checks the actual posterior,
#' on every figure build, forever. Costs microseconds. If it fails, every curve in
#' the suite is suspect.
#'
#' On the tolerance (measured 2026-07-31, results/tier1, do not tighten blindly):
#' cmdstan writes draws to ~8 significant figures, so recomputing p_pred from the
#' STORED (rounded) parameters is not bit-reproducible. Measured on this fit:
#'   Stan recomputed from its own stored draws vs stored p_pred : 3.0e-8
#'   model_math          vs Stan recomputed from same draws     : 5.0e-9
#'   model_math          vs stored p_pred                       : 3.4e-8
#' i.e. Stan disagrees with ITSELF (3.0e-8) by more than it disagrees with
#' model_math (5.0e-9); 3.4e-8 is CSV round-trip, not a math difference. 1e-6 sits
#' well above that floor and four orders below any real logic error, which shows up
#' as O(0.01-1). test_obs_prob_parity.R holds the tight 1e-7 bound because it feeds
#' exact parameter values in and only p_pred is rounded.
curve_parity_check <- function(fit, stan_data, tol = 1e-6) {
  obs <- attr(stan_data, "obs")
  p   <- mm_draws(fit, T_ref = stan_data$T_ref %||% 38.0)
  p_mm   <- mm_obs_prob(obs$group, obs$dose_cfu, obs$CoP, obs$T_thresh,
                        obs$gilman_stratum, p)
  p_stan <- posterior::as_draws_matrix(fit$draws("p_pred"))
  if (!identical(dim(p_mm), dim(p_stan)))
    stop("curve parity: stan_data does not match this fit (", ncol(p_mm), " obs vs ",
         ncol(p_stan), " p_pred columns)", call. = FALSE)
  d <- max(abs(p_mm - p_stan))
  if (d > tol)
    stop(sprintf("curve parity FAILED: figure math disagrees with p_pred (max %.3e > %.1e)",
                 d, tol), call. = FALSE)
  invisible(d)
}

# ---- data assembly -----------------------------------------------------------

#' Compute every curve, spaghetti, and data point for a set of groupings.
#' Computed ONCE and reused by the combined grid, the per-era grids, and the
#' per-grouping files.
curve_grid_data <- function(fit, stan_data, specs = curve_specs(),
                            ngrid = 60, n_spaghetti = 100, ndraw_ribbon = 1000,
                            seed = 1, x_scale = c("fixed", "free"),
                            drop_empty_columns = TRUE) {
  x_scale <- match.arg(x_scale)
  obs   <- attr(stan_data, "obs")
  T_ref <- stan_data$T_ref %||% 38.0
  chk   <- validate_curve_specs(specs, obs)

  # Validate against the FULL spec set above -- so a mis-specified regex still
  # hard-errors -- then draw only the columns this tier has data for. The grouped and
  # per-subject Darton columns are mutually exclusive by construction, so exactly one
  # of them is empty in every tier; a column with a curve and no data beneath it
  # invites reading the model as evidence.
  if (isTRUE(drop_empty_columns)) {
    keep <- specs$col_key %in% unique(chk$col_key)
    if (!any(keep)) stop("curve_grid_data: no spec column matches any observation",
                         call. = FALSE)
    specs <- specs[keep, , drop = FALSE]
  }

  p_full <- mm_draws(fit, T_ref = T_ref)
  p_rib  <- mm_thin(p_full, ndraw_ribbon, seed)
  # spaghetti draws are a SUBSET of the ribbon draws, and the same thin is used in
  # every panel -- so .draw = 37 is the same posterior sample everywhere and a
  # reader can trace one draw across the whole grid.
  p_spag <- mm_thin(p_rib, n_spaghetti, seed)

  rng <- if (x_scale == "fixed") c(min(specs$dose_lo), max(specs$dose_hi)) else NULL

  cur <- list(); spag <- list(); casc <- list(); copr <- list()
  for (i in seq_len(nrow(specs))) {
    sp  <- specs[i, ]
    lim <- rng %||% c(sp$dose_lo, sp$dose_hi)
    dose <- 10^seq(log10(lim[1]), log10(lim[2]), length.out = ngrid)
    cop_med <- if (sp$cop_mode == "individual")
      stats::median(chk$CoP[chk$col_key == sp$col_key], na.rm = TRUE) else NULL

    for (rk in c("p_inf", "p_fev", "p_fevginf", "phi")) {
      m_r <- mm_curve(rk, dose, sp, p_rib,  cop_override = cop_med)
      m_s <- mm_curve(rk, dose, sp, p_spag, cop_override = cop_med)
      cur[[length(cur) + 1L]]   <- mm_quantiles(m_r, dose) %>%
        mutate(col_key = sp$col_key, row_key = rk)
      spag[[length(spag) + 1L]] <- mm_spaghetti(m_s, dose, p_spag) %>%
        mutate(col_key = sp$col_key, row_key = rk)

      # For a per-subject-CoP column the single curve at the cohort median hides
      # the real spread in immunity. Add a second band spanning the 10th-90th
      # percentile of the OBSERVED subject titres (population heterogeneity), which
      # is a different object from the posterior band (parameter uncertainty).
      if (sp$cop_mode == "individual") {
        q <- stats::quantile(chk$CoP[chk$col_key == sp$col_key], c(0.1, 0.9), na.rm = TRUE)
        a <- matrixStats::colMedians(mm_curve(rk, dose, sp, p_rib, cop_override = q[1]))
        b <- matrixStats::colMedians(mm_curve(rk, dose, sp, p_rib, cop_override = q[2]))
        copr[[length(copr) + 1L]] <- data.frame(x = dose, lo = pmin(a, b), hi = pmax(a, b),
                                                col_key = sp$col_key, row_key = rk)
      }
    }
    fac <- mm_cascade(dose, sp, p_rib, cop_override = cop_med)
    casc[[length(casc) + 1L]] <- bind_rows(lapply(names(fac), function(nm)
      data.frame(x = dose, med = matrixStats::colMedians(fac[[nm]]), factor_key = nm))) %>%
      mutate(col_key = sp$col_key, row_key = "cascade")
  }

  # ---- observed data points --------------------------------------------------
  pp  <- posterior::as_draws_matrix(fit$draws("p_pred"))
  ci  <- .wilson(obs$y, obs$n)
  pts <- chk %>% mutate(.idx = seq_len(n()), lo = ci$lo, hi = ci$hi,
                        fitted = apply(pp, 2, stats::median),
                        in_likelihood = TRUE)
  # n=1 individual endpoints stack at one dose with y in {0,1}; collapse them to a
  # single Wilson-CI point per (grouping, row). Their per-subject CoP resolution
  # belongs on the titre axis (titre_protection.png), not here.
  singles <- pts %>% filter(n == 1L)
  if (nrow(singles)) {
    agg <- singles %>% group_by(col_key, row_key, study, dose_cfu) %>%
      summarise(y = sum(y), n = dplyr::n(), idx = list(.idx), .groups = "drop") %>%
      rowwise() %>%
      mutate(obs_rate = y / n,
             fitted   = stats::median(rowMeans(pp[, unlist(idx), drop = FALSE]))) %>%
      ungroup() %>% select(-idx)
    a_ci <- .wilson(agg$y, agg$n)
    agg  <- agg %>% mutate(lo = a_ci$lo, hi = a_ci$hi, in_likelihood = TRUE,
                           aggregated = TRUE)
    pts <- bind_rows(pts %>% filter(n > 1L) %>% mutate(aggregated = FALSE), agg)
  } else pts$aggregated <- FALSE

  list(curves    = bind_rows(cur),
       spaghetti = bind_rows(spag),
       cascade   = bind_rows(casc),
       cop_range = if (length(copr)) bind_rows(copr) else NULL,
       points    = pts,
       specs     = specs,
       n_spaghetti = p_spag$.ndraws, ndraw_ribbon = p_rib$.ndraws,
       x_scale = x_scale, T_ref = T_ref)
}

# ---- plotting ----------------------------------------------------------------

.ROW_LEVELS <- CURVE_ROWS$row_key
.row_labeller <- function(x) unname(setNames(CURVE_ROWS$row_label, CURVE_ROWS$row_key)[x])

#' @param gd output of curve_grid_data()
#' @param keys which col_keys to draw (default all)
plot_curve_grid <- function(gd, keys = gd$specs$col_key, title = NULL, subtitle = NULL) {
  sp <- gd$specs %>% filter(col_key %in% keys)
  lab <- setNames(sp$col_label, sp$col_key)
  xlim <- c(min(sp$dose_lo), max(sp$dose_hi))
  # Pre-filter to the file's dose range rather than letting scale_x_log10(limits=)
  # drop rows: the scale-drop path silently discards data and warns once per layer
  # per panel, and truncates lines at the stat stage instead of at draw time.
  f <- function(d, xcol = "x") d %>% filter(col_key %in% keys,
                                            .data[[xcol]] >= xlim[1] * (1 - 1e-9),
                                            .data[[xcol]] <= xlim[2] * (1 + 1e-9)) %>%
    mutate(row_key = factor(row_key, .ROW_LEVELS),
           col_key = factor(col_key, sp$col_key))
  cur <- f(gd$curves); spg <- f(gd$spaghetti); cas <- f(gd$cascade)
  pts <- f(gd$points, "dose_cfu")
  cpr <- if (!is.null(gd$cop_range)) f(gd$cop_range) else NULL
  n_col <- length(keys)

  # Layer order matters: the ribbon must sit UNDER the spaghetti, otherwise it
  # paints over the draw fan that is the whole purpose of this figure.
  p <- ggplot() +
    {if (!is.null(cpr) && nrow(cpr))
      geom_ribbon(data = cpr, aes(x, ymin = lo, ymax = hi), fill = "grey35", alpha = 0.16)} +
    geom_ribbon(data = cur, aes(x, ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.10) +
    geom_line(data = spg, aes(x, y, group = .draw), colour = "steelblue",
              linewidth = 0.18, alpha = .fig_spag_alpha(gd$n_spaghetti), lineend = "round") +
    geom_line(data = cur, aes(x, med), colour = "steelblue4", linewidth = 0.7) +
    geom_line(data = cas, aes(x, med, colour = factor_key), linewidth = 0.65) +
    geom_errorbar(data = pts, aes(x = dose_cfu, ymin = lo, ymax = hi),
                  width = 0.10, colour = "grey45", linewidth = 0.35) +
    geom_point(data = pts, aes(x = dose_cfu, y = obs_rate, fill = study),
               shape = 21, size = 2.1, colour = "grey15", stroke = 0.3) +
    geom_point(data = pts, aes(x = dose_cfu, y = fitted), shape = 4, size = 1.7,
               stroke = 0.6, colour = "grey10") +
    facet_grid(row_key ~ col_key, switch = "y",
               labeller = labeller(row_key = .row_labeller, col_key = lab)) +
    .fig_scale_dose(limits = xlim) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = 0.03)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_colour_manual(values = c("P_inf" = "#1b7837", "P_fev|inf" = "#762a83",
                                   "P_inf x P_fev|inf" = "#7f7f7f", "phi(T,D)" = "#d95f02",
                                   "observed-scale P(fever)" = "#111111"),
                        breaks = CASCADE_FACTORS, name = "cascade factor") +
    scale_fill_brewer(palette = "Set2", name = "study (observed)") +
    labs(x = "challenge dose (CFU)", y = NULL, title = title, subtitle = subtitle) +
    .fig_theme_grid(n_col) +
    guides(colour = guide_legend(order = 2, override.aes = list(linewidth = 1.2)),
           fill = guide_legend(order = 1, override.aes = list(size = 2.6)))
  p
}

.GRID_SUBTITLE <- paste0(
  "Rows are the model's factors; the last row overlays them. blue spaghetti: individual posterior draws (THE SAME draws in every panel, so one draw can be traced across the grid);  ",
  "blue band+line: 90%% and median;  grey band: 10-90%% of observed per-subject titres (population spread, not parameter uncertainty);  ",
  "filled point: observed (Wilson 95%%);  x: Stan p_pred\n",
  "phi(T,D) is the fitted fever-DEFINITION map. It multiplies the Maryland fever likelihood; for Oxford, phi==1 by construction ",
  "(the composite TD endpoint at T_ref=%.1fC IS the reference definition), so the Oxford phi row shows what fraction of TD+ subjects would cross a strict %.1fC threshold -- not a likelihood factor.\n",
  "TWO THINGS THE GRID EXPOSES: (1) the model has no definition map for INFECTION, so P(infection|D) is identical across the Hornick / Levine / Gilman columns although their infection endpoints genuinely differ ",
  "(stool+blood culture vs any-time stool vs late shedding).\n",
  "(2) Oxford infection rows are Tier-2 only, so in a Tier-1 fit the P(infection|D) panels carry a model curve with no data beneath it except for the Darton per-subject cohort.")

#' Build and write the whole Suite-1 output set for one run directory.
plot_grouping_grid <- function(fit, stan_data, out_dir, specs = curve_specs(),
                               ngrid = 60, n_spaghetti = 100, ndraw_ribbon = 1000,
                               seed = 1, label = "Tier 1", write_columns = TRUE) {
  gd <- curve_grid_data(fit, stan_data, specs, ngrid, n_spaghetti, ndraw_ribbon, seed)
  specs <- gd$specs   # curve_grid_data() drops columns with no data in this tier
  sub <- sprintf(.GRID_SUBTITLE, gd$T_ref, gd$T_ref)

  # combined grid: one common dose range so columns are directly comparable
  n <- nrow(specs); sz <- .fig_size_grid(n, 5)
  .fig_save(plot_curve_grid(gd, specs$col_key,
              sprintf("%s: model prediction for every comparable dataset sub-grouping", label), sub),
            file.path(out_dir, "dose_response_grid_all.png"), sz$w, sz$h)

  # per-era grids: the readable view (one vehicle frame, one dose range)
  for (er in unique(specs$era)) {
    k <- specs$col_key[specs$era == er]
    gde <- gd; gde$specs <- specs %>% filter(era == er)
    sz <- .fig_size_grid(length(k), 5)
    .fig_save(plot_curve_grid(gde, k,
                sprintf("%s: %s groupings (%s vehicle)", label, er,
                        unique(specs$vehicle[specs$era == er])), sub),
              file.path(out_dir, sprintf("dose_response_grid_%s.png", er)), sz$w, sz$h)
  }

  # one file per grouping (short subtitle -- the full one clips at this width)
  if (write_columns) for (k in specs$col_key) {
    gdc <- gd; gdc$specs <- specs %>% filter(col_key == k)
    sz <- .fig_size_grid(1, 5, unit_w = 5.2, pad_w = 1.6)
    .fig_save(plot_curve_grid(gdc, k, sprintf("%s: %s", label, k),
                paste0("spaghetti: posterior draws;  band: 90%;  grey band: 10-90% of observed ",
                       "subject titres;  point: observed (Wilson 95%);  x: Stan p_pred")),
              file.path(out_dir, "grid", paste0(k, ".png")), sz$w, sz$h)
  }
  invisible(gd)
}
