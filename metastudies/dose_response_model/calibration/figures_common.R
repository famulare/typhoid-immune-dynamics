#' Shared save / size / theme helpers for the model figure suite.

suppressPackageStartupMessages({library(ggplot2)})
if (!exists("calib_dir")) source("utils.R")   # %||%, .wilson(), calib_path()

#' Save a ggplot. Uses ragg when available (2-3x faster than the default device on
#' the line-heavy spaghetti layers, and smaller output), always sets a white
#' background (the pre-2026-07-31 dose_response_curves.R did not), and never lets
#' a figure error take down a fit.
.fig_save <- function(p, path, w, h, dpi = 150) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
  tryCatch({
    ggplot2::ggsave(path, p, width = w, height = h, dpi = dpi, bg = "white",
                    device = dev, limitsize = FALSE)
    message("  figure: ", path)
  }, error = function(e) message("  [skip] ", basename(path), ": ", conditionMessage(e)))
  invisible(path)
}

#' Grid sizing: wide enough to read on screen, scaling with the column count.
#' Deliberately not print-shaped -- these are meant to be scrolled, not printed.
.fig_size_grid <- function(n_col, n_row, unit_w = 2.6, unit_h = 2.2,
                           pad_w = 1.5, pad_h = 1.8, w_max = 30) {
  list(w = min(pad_w + unit_w * n_col, w_max), h = pad_h + unit_h * n_row)
}

#' Theme for the grid; strip text shrinks as columns multiply.
.fig_theme_grid <- function(n_col) {
  theme_minimal(base_size = 11) +
    theme(strip.text.x   = element_text(size = if (n_col > 6) 6.5 else 8, lineheight = 1.05),
          strip.text.y   = element_text(size = 8.5, angle = 0),
          axis.text.x    = element_text(size = 6.5),
          axis.text.y    = element_text(size = 7),
          panel.spacing.x = unit(3, "pt"),
          panel.spacing.y = unit(5, "pt"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.box = "vertical",
          legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7.5),
          plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8))
}

#' log10 dose axis, shared by every dose-axis figure in the suite.
.fig_scale_dose <- function(...)
  scale_x_log10(breaks = 10^seq(0, 12, by = 2),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                guide = guide_axis(check.overlap = TRUE), ...)

#' Alpha for spaghetti. Tuned so the draw fan is legible as the point of the
#' figure rather than washing out into the ribbon: at n=100 this gives 0.08.
#' Lower values (3/n) rendered as a faint fuzz indistinguishable from the band.
.fig_spag_alpha <- function(n_draw) max(0.03, min(0.25, 8 / n_draw))
