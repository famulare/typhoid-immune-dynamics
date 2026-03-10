#' ---
#' title: "E1 Bridge Model: Vi HA ↔ Barrett ELISA (1983)"
#' author: "Claude YOLO Implementation"
#' date: "2026-02-04"
#' output: md_document
#' ---

#' # E1 Bridge: Hemagglutination to Barrett ELISA
#'
#' This script fits a regression model to Barrett 1983 Figure 2 data,
#' bridging Vi hemagglutination (HA) titers to ELISA titers.
#'
#' ## Data Source
#' Barrett TJ et al. J Clin Microbiol 1983;17(4):625-627, Figure 2.
#' 77 sera from acute typhoid patients, El Salvador.

# Setup
library(tidyverse)

# Load extracted data
data_path <- here::here("metastudies/anti_Vi_metrology_ladder/claude_yolo/extraction/barrett_1983_fig2_data.csv")
output_path <- here::here("metastudies/anti_Vi_metrology_ladder/claude_yolo/intermediates/")

barrett_data <- read_csv(data_path, comment = "#", show_col_types = FALSE)

#' ## Data Preparation
#'
#' Both titers are on 2-fold dilution scales. We transform to log2 scale
#' for linear regression. Values <20 are coded as 10 (midpoint).

barrett_data <- barrett_data %>%
  mutate(
    log2_ELISA = log2(ELISA_titer),
    log2_HA = log2(HA_titer)
  )

# Expand aggregated data to individual observations for regression
barrett_expanded <- barrett_data %>%
  uncount(count)

cat("Data summary:\n")
cat("  Total observations:", nrow(barrett_expanded), "\n")
cat("  Unique ELISA levels:", length(unique(barrett_expanded$ELISA_titer)), "\n")
cat("  Unique HA levels:", length(unique(barrett_expanded$HA_titer)), "\n")

#' ## Model Fitting
#'
#' Model: log2(ELISA_titer) ~ Normal(α + β × log2(HA_titer), σ²)
#'
#' This treats ELISA as the "outcome" and HA as the "predictor",
#' which aligns with the direction of the bridge (HA → ELISA → modern).

# Fit linear regression
e1_model <- lm(log2_ELISA ~ log2_HA, data = barrett_expanded)

# Model summary
cat("\n=== E1 Model Summary ===\n")
print(summary(e1_model))

# Extract parameters
e1_params <- list(
  intercept = coef(e1_model)[1],
  slope = coef(e1_model)[2],
  sigma = sigma(e1_model),
  r_squared = summary(e1_model)$r.squared,
  n = nrow(barrett_expanded)
)

cat("\n=== E1 Parameters ===\n")
cat("Intercept (α):", round(e1_params$intercept, 4), "\n")
cat("Slope (β):", round(e1_params$slope, 4), "\n")
cat("Residual SD (σ):", round(e1_params$sigma, 4), "log2 units\n")
cat("R²:", round(e1_params$r_squared, 4), "\n")
cat("n:", e1_params$n, "\n")

#' ## Confidence Intervals
#'
#' We use the standard errors from the regression to quantify
#' parameter uncertainty.

e1_params$intercept_se <- summary(e1_model)$coefficients[1, 2]
e1_params$slope_se <- summary(e1_model)$coefficients[2, 2]

cat("\n=== Parameter Uncertainty ===\n")
cat("Intercept SE:", round(e1_params$intercept_se, 4), "\n")
cat("Slope SE:", round(e1_params$slope_se, 4), "\n")

# 95% CIs
e1_params$intercept_ci <- confint(e1_model)[1, ]
e1_params$slope_ci <- confint(e1_model)[2, ]

cat("Intercept 95% CI: [", round(e1_params$intercept_ci[1], 4), ",",
    round(e1_params$intercept_ci[2], 4), "]\n")
cat("Slope 95% CI: [", round(e1_params$slope_ci[1], 4), ",",
    round(e1_params$slope_ci[2], 4), "]\n")

#' ## Bridge Function
#'
#' The E1 bridge converts HA titer to predicted ELISA titer (log2 scale).

bridge_e1 <- function(log2_HA, params = e1_params, add_noise = FALSE) {
    #' Convert log2(HA titer) to predicted log2(ELISA titer)
    #'
    #' @param log2_HA Log2-transformed HA titer (can be vector)
    #' @param params List containing intercept, slope, sigma
    #' @param add_noise If TRUE, add residual noise (for simulation)
    #' @return Predicted log2(ELISA titer)

  pred <- params$intercept + params$slope * log2_HA

  if (add_noise) {
    pred <- pred + rnorm(length(pred), 0, params$sigma)
  }

  return(pred)
}

#' ## Validation: Diagnostic Plots

# Create plot data
plot_data <- barrett_data %>%
  mutate(
    predicted = bridge_e1(log2_HA),
    observed_mean = log2_ELISA
  )

# Scatter plot with regression line
p1 <- ggplot(barrett_expanded, aes(x = log2_HA, y = log2_ELISA)) +
  geom_count(alpha = 0.6) +
  geom_abline(intercept = e1_params$intercept, slope = e1_params$slope,
              color = "red", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "E1 Bridge: Vi HA → Barrett ELISA (1983)",
    subtitle = sprintf("log₂(ELISA) = %.2f + %.2f × log₂(HA), R² = %.3f",
                       e1_params$intercept, e1_params$slope, e1_params$r_squared),
    x = "log₂(HA titer)",
    y = "log₂(ELISA titer)",
    size = "Count"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = log2(c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120)),
    labels = c("<20", "20", "40", "80", "160", "320", "640", "1280", "2560", "5120")
  ) +
  scale_y_continuous(
    breaks = log2(c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120)),
    labels = c("<20", "20", "40", "80", "160", "320", "640", "1280", "2560", "5120")
  )

ggsave(file.path(output_path, "e1_diagnostic_scatter.png"), p1,
       width = 8, height = 6, dpi = 150)

# Residual plot
residuals_df <- data.frame(
  fitted = fitted(e1_model),
  residuals = residuals(e1_model)
)

p2 <- ggplot(residuals_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-2, 2) * e1_params$sigma, linetype = "dotted", color = "gray50") +
  labs(
    title = "E1 Model Residuals",
    x = "Fitted values (log₂ ELISA)",
    y = "Residuals (log₂ units)"
  ) +
  theme_minimal()

ggsave(file.path(output_path, "e1_residual_plot.png"), p2,
       width = 8, height = 5, dpi = 150)

#' ## Save Model Object

saveRDS(list(
  model = e1_model,
  params = e1_params,
  bridge_function = bridge_e1,
  data = barrett_expanded
), file.path(output_path, "e1_model_fit.rds"))

cat("\n=== Model saved to", file.path(output_path, "e1_model_fit.rds"), "===\n")

#' ## Interpretation
#'
#' The E1 bridge shows that HA and ELISA titers are positively correlated,
#' but with substantial scatter (σ ≈ 1.5-2 log2 units typical for dilution assays).
#'
#' Key findings:
#' - Slope near 1 suggests roughly linear relationship on log scale
#' - Large residual variance reflects assay imprecision and biological variability
#' - Double-negative samples (37 of 77) dominate the lower left corner
#'
#' Limitations:
#' - Left-censoring at <20 may bias intercept estimate
#' - Population is acute typhoid, not healthy vaccinees
#' - Discrete titer values create heaping in the data

#' ## Parameter Export for Cascade

e1_export <- data.frame(
  parameter = c("intercept", "slope", "sigma", "r_squared", "n",
                "intercept_se", "slope_se"),
  value = c(e1_params$intercept, e1_params$slope, e1_params$sigma,
            e1_params$r_squared, e1_params$n,
            e1_params$intercept_se, e1_params$slope_se),
  unit = c("log2", "dimensionless", "log2", "proportion", "count",
           "log2", "dimensionless"),
  description = c("Regression intercept",
                  "Regression slope",
                  "Residual standard deviation",
                  "Coefficient of determination",
                  "Sample size",
                  "Standard error of intercept",
                  "Standard error of slope")
)

write_csv(e1_export, file.path(output_path, "e1_parameters.csv"))
cat("Parameters exported to", file.path(output_path, "e1_parameters.csv"), "\n")
