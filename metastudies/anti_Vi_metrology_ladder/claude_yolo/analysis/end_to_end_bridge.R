#' ---
#' title: "End-to-End Bridge: Vi HA → VaccZyme ELISA"
#' author: "Claude YOLO Implementation"
#' date: "2026-02-04"
#' output: md_document
#' ---

#' # End-to-End Metrology Bridge
#'
#' This script implements the full cascade from historical Vi HA titers
#' to modern VaccZyme ELISA values, incorporating:
#' - E1: HA → Barrett ELISA (data-driven, Barrett 1983)
#' - E2: Barrett ELISA → Modern In-house ELISA (STRUCTURAL ASSUMPTION)
#' - E3: In-house ELISA → VaccZyme (data-driven, Lee 2020)
#'
#' ## CRITICAL CAVEAT
#' E2 has NO direct data support. The bridge assumes that Barrett 1983 ELISA
#' and modern in-house ELISAs are calibratable to the same weight-based units.
#' A structural uncertainty parameter (τ) accounts for unknown systematic bias.

# Setup
library(tidyverse)
set.seed(20260204)  # Reproducibility

# Use relative paths from project root
base_path <- "metastudies/anti_Vi_metrology_ladder/claude_yolo/"
output_path <- file.path(base_path, "intermediates/")
final_path <- file.path(base_path, "outputs/")

# Create outputs directory if needed
dir.create(final_path, showWarnings = FALSE, recursive = TRUE)

#' ## Load E1 and E3 Parameters

# E1: If already fitted, load; otherwise use extracted data
e1_data_path <- file.path(base_path, "extraction/barrett_1983_fig2_data.csv")
barrett_data <- read_csv(e1_data_path, comment = "#", show_col_types = FALSE)

# Fit E1 model
barrett_expanded <- barrett_data %>%
  mutate(
    log2_ELISA = log2(ELISA_titer),
    log2_HA = log2(HA_titer)
  ) %>%
  uncount(count)

e1_model <- lm(log2_ELISA ~ log2_HA, data = barrett_expanded)

e1_params <- list(
  intercept = coef(e1_model)[1],
  slope = coef(e1_model)[2],
  sigma = sigma(e1_model),
  intercept_se = summary(e1_model)$coefficients[1, 2],
  slope_se = summary(e1_model)$coefficients[2, 2],
  vcov = vcov(e1_model)
)

cat("=== E1 Parameters ===\n")
cat("Intercept:", round(e1_params$intercept, 3), "\n")
cat("Slope:", round(e1_params$slope, 3), "\n")
cat("Sigma:", round(e1_params$sigma, 3), "log2 units\n")

# E3: From Lee 2020
e3_params <- list(
  intercept = 3.749,          # log scale (natural log)
  slope = 0.946,
  sigma = 0.30               # conservative estimate
)

cat("\n=== E3 Parameters ===\n")
cat("Intercept:", e3_params$intercept, "\n")
cat("Slope:", e3_params$slope, "\n")
cat("Sigma:", e3_params$sigma, "log units\n")

#' ## E2: Structural Assumption
#'
#' E2 bridges Barrett 1983 ELISA to modern in-house ELISA.
#' Since both can report in μg/ml (via precipitin or Vi-IgGR1,2011 calibration),
#' we assume they measure the same quantity with an unknown offset and spread.
#'
#' Model: log(modern_μg) = log(Barrett_μg) + Normal(μ_offset, τ²)
#'
#' Where:
#' - μ_offset: Unknown systematic calibration difference (assume 0)
#' - τ: Structural uncertainty (tunable parameter)

e2_params <- list(
  mu_offset = 0,              # Assume no systematic bias (hopeful!)
  tau = 1.0                   # 1 log unit (~2.7× uncertainty factor)
)

cat("\n=== E2 Parameters (STRUCTURAL) ===\n")
cat("Offset (μ):", e2_params$mu_offset, "(assumed)\n")
cat("Uncertainty (τ):", e2_params$tau, "log units\n")
cat("WARNING: E2 has NO direct data support!\n")

#' ## Unit Conversion Notes
#'
#' - E1 operates on log2 scale (2-fold dilution titers)
#' - E2 and E3 operate on natural log scale (μg/ml, EU/ml)
#'
#' Conversion: log(x) = log2(x) × ln(2) ≈ log2(x) × 0.693
#'
#' Barrett ELISA titers don't directly report μg/ml, but the titer
#' can be interpreted as relative concentration. For bridging, we
#' assume the titer scale maps to concentration via calibration.

#' ## Monte Carlo Cascade Function

cascade_ha_to_vacczyme <- function(
    log2_HA,
    n_samples = 1000,
    e1 = e1_params,
    e2 = e2_params,
    e3 = e3_params,
    include_parameter_uncertainty = TRUE
    ) {
    #' Propagate HA titer through the full bridge to VaccZyme
    #'
    #' @param log2_HA Log2-transformed HA titer (scalar or vector)
    #' @param n_samples Number of Monte Carlo samples per input
    #' @param e1 E1 model parameters
    #' @param e2 E2 structural parameters
    #' @param e3 E3 model parameters
    #' @param include_parameter_uncertainty If TRUE, sample from parameter distributions
    #' @return Matrix of samples (rows = samples, cols = input values)

  n_inputs <- length(log2_HA)
  results <- matrix(NA, nrow = n_samples, ncol = n_inputs)

  for (i in seq_len(n_samples)) {
    # Step 1: E1 - HA to Barrett ELISA (log2 scale)
    if (include_parameter_uncertainty) {
      # Sample parameters from multivariate normal
      param_sample <- MASS::mvrnorm(1, c(e1$intercept, e1$slope), e1$vcov)
      e1_int <- param_sample[1]
      e1_slp <- param_sample[2]
    } else {
      e1_int <- e1$intercept
      e1_slp <- e1$slope
    }

    log2_barrett <- e1_int + e1_slp * log2_HA + rnorm(n_inputs, 0, e1$sigma)

    # Step 2: Convert log2(Barrett titer) to log(Barrett μg/ml)
    # Barrett titer of 320 corresponds roughly to high responder
    # Vi-IgGR1,2011 at 33 μg/ml has titer ~160-320 based on context
    # Crude calibration: titer 160 ≈ 33 μg/ml → log2(160) = 7.32 → 33 μg/ml
    # Therefore: μg/ml ≈ 33 × 2^(log2_titer - 7.32) / something
    #
    # Actually, this is problematic because Barrett used visual endpoint,
    # not standardized reference. We'll use a rough conversion factor.
    #
    # Simplification: Assume Barrett titer 160 ≈ 10 μg/ml (midrange responder)
    # log(μg/ml) = log(10) + (log2_titer - log2(160)) × ln(2)
    #            = 2.303 + (log2_titer - 7.32) × 0.693

    log2_ref_titer <- 7.32   # log2(160)
    log_ref_ug <- log(10)    # Assume titer 160 ≈ 10 μg/ml

    log_barrett_ug <- log_ref_ug + (log2_barrett - log2_ref_titer) * log(2)

    # Step 3: E2 - Barrett to Modern (STRUCTURAL)
    log_modern_ug <- log_barrett_ug + e2$mu_offset + rnorm(n_inputs, 0, e2$tau)

    # Step 4: E3 - Modern to VaccZyme
    log_vacczyme <- e3$intercept + e3$slope * log_modern_ug + rnorm(n_inputs, 0, e3$sigma)

    results[i, ] <- log_vacczyme
  }

  return(results)
}

#' ## Example: Bridge Historical HA Titers to VaccZyme
#'
#' Let's bridge a range of HA titers commonly seen in carriers/patients.

test_HA_titers <- c(20, 40, 80, 160, 320, 640, 1280)
log2_test_HA <- log2(test_HA_titers)

cat("\n=== Running Monte Carlo Cascade ===\n")
samples <- cascade_ha_to_vacczyme(log2_test_HA, n_samples = 5000)

# Summarize results
results_summary <- data.frame(
  HA_titer = test_HA_titers,
  log2_HA = log2_test_HA,
  VaccZyme_median = apply(samples, 2, function(x) exp(median(x))),
  VaccZyme_mean = apply(samples, 2, function(x) exp(mean(x))),
  VaccZyme_lo95 = apply(samples, 2, function(x) exp(quantile(x, 0.025))),
  VaccZyme_hi95 = apply(samples, 2, function(x) exp(quantile(x, 0.975))),
  VaccZyme_lo50 = apply(samples, 2, function(x) exp(quantile(x, 0.25))),
  VaccZyme_hi50 = apply(samples, 2, function(x) exp(quantile(x, 0.75)))
)

cat("\n=== Bridge Results: HA Titer → VaccZyme EU/ml ===\n")
print(results_summary %>%
        mutate(across(starts_with("VaccZyme"), ~round(., 1))) %>%
        select(HA_titer, VaccZyme_median, VaccZyme_lo95, VaccZyme_hi95))

#' ## Visualization

# Long format for plotting
samples_long <- as.data.frame(samples) %>%
  setNames(paste0("HA_", test_HA_titers)) %>%
  pivot_longer(everything(), names_to = "HA_titer", values_to = "log_VaccZyme") %>%
  mutate(
    HA_titer = as.numeric(gsub("HA_", "", HA_titer)),
    VaccZyme_EU = exp(log_VaccZyme)
  )

p1 <- ggplot(samples_long, aes(x = factor(HA_titer), y = VaccZyme_EU)) +
  geom_violin(fill = "steelblue", alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    labels = scales::comma
  ) +
  labs(
    title = "End-to-End Bridge: Vi HA Titer → VaccZyme EU/ml",
    subtitle = "E2 (Barrett→Modern) uses structural assumption (τ = 1.0 log unit)",
    x = "Historical HA Titer (reciprocal)",
    y = "Predicted VaccZyme (EU/ml)",
    caption = "Violins show Monte Carlo posterior predictive distribution"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(final_path, "end_to_end_bridge.png"), p1,
       width = 10, height = 6, dpi = 150)

#' ## Sensitivity Analysis: E2 Uncertainty (τ)
#'
#' How does the choice of τ affect the bridge uncertainty?

tau_values <- c(0.5, 1.0, 1.5, 2.0)
tau_sensitivity <- list()

for (tau in tau_values) {
  e2_test <- list(mu_offset = 0, tau = tau)
  samples_tau <- cascade_ha_to_vacczyme(log2(160), n_samples = 2000,
                                        e2 = e2_test)
  tau_sensitivity[[as.character(tau)]] <- data.frame(
    tau = tau,
    log_VaccZyme = as.vector(samples_tau)
  )
}

tau_df <- bind_rows(tau_sensitivity)

p2 <- ggplot(tau_df, aes(x = exp(log_VaccZyme), fill = factor(tau))) +
  geom_density(alpha = 0.4) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    labels = scales::comma
  ) +
  labs(
    title = "Sensitivity to E2 Structural Uncertainty (τ)",
    subtitle = "Input: HA titer = 160",
    x = "Predicted VaccZyme (EU/ml)",
    y = "Density",
    fill = "τ (log units)"
  ) +
  theme_minimal()

ggsave(file.path(final_path, "tau_sensitivity.png"), p2,
       width = 8, height = 5, dpi = 150)

#' ## Uncertainty Decomposition
#'
#' How much does each edge contribute to total uncertainty?

# Run with only residual uncertainty (no parameter uncertainty)
samples_residual_only <- cascade_ha_to_vacczyme(
  log2(160),
  n_samples = 2000,
  include_parameter_uncertainty = FALSE
)

# Run with E2 tau = 0 (impossible ideal)
e2_perfect <- list(mu_offset = 0, tau = 0.001)
samples_no_e2 <- cascade_ha_to_vacczyme(
  log2(160),
  n_samples = 2000,
  e2 = e2_perfect
)

uncertainty_decomp <- data.frame(
  source = c("Full cascade (τ=1.0)",
             "Without E2 uncertainty (τ≈0)",
             "Without parameter uncertainty"),
  SD_log = c(sd(samples[,4]),  # HA=160 is column 4
             sd(samples_no_e2),
             sd(samples_residual_only)),
  fold_95 = c(exp(2 * sd(samples[,4])),
              exp(2 * sd(samples_no_e2)),
              exp(2 * sd(samples_residual_only)))
)

cat("\n=== Uncertainty Decomposition (HA=160) ===\n")
print(uncertainty_decomp %>% mutate(across(where(is.numeric), ~round(., 2))))

#' ## Save Results

# Save samples
saveRDS(list(
  samples = samples,
  HA_titers = test_HA_titers,
  summary = results_summary,
  e1_params = e1_params,
  e2_params = e2_params,
  e3_params = e3_params,
  n_samples = 5000,
  calibration_assumption = "Barrett titer 160 ≈ 10 μg/ml"
), file.path(output_path, "cascade_samples.rds"))

# Export summary to CSV
write_csv(results_summary, file.path(final_path, "parameter_estimates.csv"))

cat("\n=== Results saved ===\n")
cat("Samples:", file.path(output_path, "cascade_samples.rds"), "\n")
cat("Summary:", file.path(final_path, "parameter_estimates.csv"), "\n")
cat("Plots:", file.path(final_path, "end_to_end_bridge.png"), "\n")
cat("       ", file.path(final_path, "tau_sensitivity.png"), "\n")

#' ## Key Findings
#'
#' 1. **Massive uncertainty**: 95% prediction intervals span 2-3 orders of magnitude
#' 2. **E2 dominates**: The structural assumption (τ=1.0) is the largest uncertainty source
#' 3. **E1 contributes**: Barrett data has substantial scatter (σ ≈ 1.5 log2 units)
#' 4. **E3 is precise**: Lee 2020 r=0.991 means E3 adds little additional uncertainty
#'
#' ## Critical Caveats
#'
#' 1. **E2 is made up**: No data supports the Barrett→Modern bridge
#' 2. **Calibration guess**: Barrett titer 160 ≈ 10 μg/ml is approximate
#' 3. **Population mismatch**: E1 from carriers, E3 from vaccinees
#' 4. **No validation**: Cannot assess accuracy without true bridging data
#'
#' ## Recommendation
#'
#' Use this bridge ONLY for:
#' - Rough order-of-magnitude comparisons
#' - Historical context (qualitative interpretation)
#' - Identifying where more data is needed
#'
#' DO NOT use for:
#' - Regulatory submissions
#' - Precise quantitative comparisons
#' - Individual-level predictions
