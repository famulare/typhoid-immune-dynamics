#' M4 parity gate (LOAD-BEARING) — see plan please-plan-a-b-smooth-ritchie.md.
#'
#' The recovery harness (simulate_recovery.R) shares obs_prob() between the data
#' simulator and the fitted model, so it is BLIND to an obs_prob() bug. This test
#' is the only guard that the unified obs_prob() faithfully reproduces the original
#' five-group likelihood. It compares the new Stan generated-quantities p_pred and
#' log_lik against an independent R transcription of the ORIGINAL per-group model
#' block, at several fixed parameter vectors, for every active Tier-1 row.
#'
#' Run:  Rscript test_obs_prob_parity.R   (exits non-zero on failure)

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})
source("priors.R"); source("data_prep.R")

# Tolerances at the floating-point / implementation-difference level (Stan pow &
# exact-lgamma binomial_lpmf vs R ^ & saddlepoint dbinom). A real logic error in
# obs_prob() would produce O(0.01-1) differences, far above these.
TOL_P  <- 1e-7   # probabilities (O(1) values)
TOL_LL <- 1e-3   # pointwise log-likelihood

# ---- Independent R reference (transcribed from the model block) ----------------
# DO NOT refactor the functions in this section to call model_math.R. They are a
# SEPARATE, hand-transcribed implementation, and that independence is the entire
# value of this gate: model_math.R is the code the figures actually run, and it is
# checked against BOTH Stan and this transcription below. DRYing these together
# would silently collapse three implementations into one and delete the check.
bp <- function(D, N50, alpha, CoP, gamma) {
  scale <- (2^(1 / alpha) - 1) / N50
  1 - (1 + D * scale)^(-alpha / CoP^gamma)
}
md_mix <- function(D, N50, alpha, gamma, pi, CoPs, CoPi)
  pi * bp(D, N50, alpha, CoPs, gamma) + (1 - pi) * bp(D, N50, alpha, CoPi, gamma)

# Dose-dependent definition sensitivity phi(T,D) = phi0(T) + (1-phi0)*P_fev_naive
# (beta_phi pinned to 1).
phi_TD_R <- function(T, D, N50i, N50f, p) {
  phi0 <- plogis(p$phi0_a - p$phi0_b * (T - T_REF))
  p_fev_naive <- bp(D, N50i, p$alpha_inf, 1, p$gamma_inf) *
                 bp(D, N50f, p$alpha_fevginf, 1, p$gamma_fevginf)
  phi0 + (1 - phi0) * p_fev_naive
}

obs_prob_R <- function(row, p) {
  N50i <- 10^p$log10_N50_inf
  N50f <- 10^(p$log10_N50_inf + p$d_fev)         # reparam: log10_N50_fevginf = inf + d_fev
  delta <- 10^p$log10_delta
  g <- row$group
  if (g == 1L) {                                  # ox_fev
    D <- row$dose_cfu
    bp(D, N50i, p$alpha_inf, row$CoP, p$gamma_inf) *
      bp(D, N50f, p$alpha_fevginf, row$CoP, p$gamma_fevginf)
  } else if (g == 3L) {                           # md_fev
    D <- row$dose_cfu / delta; st <- row$gilman_stratum
    phi <- phi_TD_R(row$T_thresh, D, N50i, N50f, p)
    pf <- function(C) bp(D, N50i, p$alpha_inf, C, p$gamma_inf) *
                      bp(D, N50f, p$alpha_fevginf, C, p$gamma_fevginf)
    if (!is.na(st) && st == 1L) phi * pf(p$CoP_susc)
    else if (!is.na(st) && st == 2L) phi * pf(p$CoP_imm)
    else phi * (p$pi_susc * pf(p$CoP_susc) + (1 - p$pi_susc) * pf(p$CoP_imm))
  } else if (g == 4L) {                           # md_inf
    md_mix(row$dose_cfu / delta, N50i, p$alpha_inf, p$gamma_inf, p$pi_susc, p$CoP_susc, p$CoP_imm)
  } else if (g == 5L) {                           # hornick_cond
    D <- row$dose_cfu / delta
    phi <- phi_TD_R(row$T_thresh, D, N50i, N50f, p)
    pf <- function(C) bp(D, N50i, p$alpha_inf, C, p$gamma_inf) *
                      bp(D, N50f, p$alpha_fevginf, C, p$gamma_fevginf)
    p_inf <- md_mix(D, N50i, p$alpha_inf, p$gamma_inf, p$pi_susc, p$CoP_susc, p$CoP_imm)
    pc <- phi * (p$pi_susc * pf(p$CoP_susc) + (1 - p$pi_susc) * pf(p$CoP_imm)) / p_inf
    min(max(pc, 1e-12), 1 - 1e-12)
  } else if (g == 6L) {                           # ox_inf_indiv: individual Oxford infection
    bp(row$dose_cfu, N50i, p$alpha_inf, row$CoP, p$gamma_inf)
  } else if (g == 7L) {                           # ox_fevginf_indiv: individual P(fever | infected)
    bp(row$dose_cfu, N50f, p$alpha_fevginf, row$CoP, p$gamma_fevginf)
  } else if (g == 2L) {                           # ox_inf: eta-corrected Oxford shedding
    # Raw dose (Oxford is the bicarbonate frame, so no /delta). eta is the
    # treatment-truncation / detection factor: 1 at zero dose, decaying to eta_lo as
    # dose rises, on the SAME N50_inf scale as the infection curve -- which is why a
    # misfit here can be absorbed by N50_inf instead of failing visibly, and why this
    # branch is the one that most needs an independent check.
    D <- row$dose_cfu
    eta <- p$eta_lo + (1 - p$eta_lo) * exp(-p$kappa * D / N50i)
    eta * bp(D, N50i, p$alpha_inf, row$CoP, p$gamma_inf)
  } else stop("obs_prob_R: unhandled likelihood group ", g)
}

# ---- Build data + compile -----------------------------------------------------
priors <- load_priors()
sd <- build_stan_data("dose_response_data.csv", priors, tier_col = "tier1_active")
obs <- attr(sd, "obs")
mod <- cmdstan_model("typhoid_dose_response.stan")

# ---- Parameter vectors to test (constrained scale; must cover the model's params) ----
PARAM_NAMES <- c("log10_N50_inf","d_fev","alpha_inf","alpha_fevginf","gamma_inf",
                 "gamma_fevginf","log10_delta","pi_susc","CoP_imm","CoP_susc",
                 "phi0_a","phi0_b","eta_lo","kappa")
# NAMED, then indexed by PARAM_NAMES: these were positional over the name vector,
# so a reordering of parameters{} would have silently permuted the gate's inputs.
vecs <- lapply(list(                   # phi0_a,phi0_b replace phi_md (beta_phi pinned=1)
  c(log10_N50_inf=2.5, d_fev=0.3, alpha_inf=0.30, alpha_fevginf=0.35, gamma_inf=0.60,
    gamma_fevginf=0.90, log10_delta=3.5, pi_susc=0.65, CoP_imm=3.0, CoP_susc=1.0,
    phi0_a=1.4, phi0_b=1.8, eta_lo=0.5, kappa=1.0),
  c(log10_N50_inf=2.0, d_fev=0.0, alpha_inf=0.15, alpha_fevginf=0.50, gamma_inf=0.20,
    gamma_fevginf=1.50, log10_delta=2.0, pi_susc=0.40, CoP_imm=5.0, CoP_susc=1.1,
    phi0_a=0.5, phi0_b=1.0, eta_lo=0.4, kappa=0.7),   # d_fev=0 edge
  c(log10_N50_inf=3.1, d_fev=1.2, alpha_inf=0.50, alpha_fevginf=0.20, gamma_inf=1.00,
    gamma_fevginf=0.30, log10_delta=4.5, pi_susc=0.80, CoP_imm=2.0, CoP_susc=0.9,
    phi0_a=2.0, phi0_b=0.5, eta_lo=0.6, kappa=1.5)
), function(v) v[PARAM_NAMES])
truth <- posterior::as_draws_matrix(do.call(rbind, lapply(vecs, function(v) setNames(v, PARAM_NAMES))))

# ---- New-model p_pred / log_lik via generate_quantities -----------------------
assert_fitted_params_match(mod, PARAM_NAMES)
gq <- mod$generate_quantities(fitted_params = truth, data = sd, seed = 1)
p_stan  <- posterior::as_draws_matrix(gq$draws("p_pred"))    # ndraws x N_obs
ll_stan <- posterior::as_draws_matrix(gq$draws("log_lik"))

# ---- Compare ------------------------------------------------------------------
# log_lik is [N_obs dose-response rows, then N_ladder phi0 threshold binomials]. The
# ladder terms were folded in so that target == lprior + sum(log_lik) (priorsense
# power-scales exactly those two objects, and previously omitted the three binomials
# that identify phi0_a/phi0_b). Check the two blocks separately.
stopifnot(ncol(ll_stan) == nrow(obs) + sd$N_ladder)
n_obs <- nrow(obs)
max_dp <- 0; max_dl <- 0; max_dlad <- 0
for (di in seq_len(nrow(p_stan))) {
  p <- as.list(vecs[[di]])
  p_ref  <- vapply(seq_len(n_obs), function(i) obs_prob_R(obs[i, ], p), numeric(1))
  ll_ref <- dbinom(obs$y, obs$n, p_ref, log = TRUE)
  max_dp <- max(max_dp, max(abs(p_stan[di, ] - p_ref)))
  max_dl <- max(max_dl, max(abs(ll_stan[di, seq_len(n_obs)] - ll_ref)))
  # independent transcription of the ladder sub-likelihood: phi0(T) is logit-linear,
  # centered at T_ref, and each threshold is a binomial over the TD+ subjects.
  phi0_ref <- plogis(p$phi0_a - p$phi0_b * (sd$ladder_T - T_REF))
  lad_ref  <- dbinom(sd$ladder_count, sd$ladder_N, phi0_ref, log = TRUE)
  max_dlad <- max(max_dlad, max(abs(ll_stan[di, n_obs + seq_len(sd$N_ladder)] - lad_ref)))
}
cat(sprintf("Parity over %d param vectors x %d rows (+%d ladder terms):\n",
            nrow(p_stan), n_obs, sd$N_ladder))
cat(sprintf("  max |p_pred_stan - p_ref|        = %.3e\n", max_dp))
cat(sprintf("  max |log_lik_stan - ll_ref|      = %.3e\n", max_dl))
cat(sprintf("  max |log_lik ladder - lad_ref|   = %.3e\n", max_dlad))
failures <- character()
if (!(max_dp < TOL_P)) failures <- c(failures, "A: Stan vs independent R on real rows (p)")
if (!(max_dl < TOL_LL)) failures <- c(failures, "A: Stan vs independent R on real rows (log_lik)")
if (!(max_dlad < TOL_LL)) failures <- c(failures, "A: Stan vs independent R on the phi0 ladder terms")

# ==============================================================================
# ARM B — model_math.R, the implementation the FIGURES run, gated against Stan.
#
# Rationale: obs_prob_R() above guards the .stan against its own refactors, but
# nothing guarded the plotting-side math (the deleted `.bp()` in
# dose_response_curves.R was an ungated third transcription). model_math.R is now
# the single plotting-side implementation, so it gets its own gate, three ways:
#   B1  Stan <-> mm_obs_prob() on real Tier-1 AND Tier-2 rows (Tier 2 is the only
#       thing that exercises group 2 / eta through the likelihood path).
#   B2  Stan <-> mm_obs_prob() on a SYNTHETIC covariate grid ~700 rows wide,
#       covering dose decades, titres beyond any observed value, thresholds
#       outside STUDY_FEVER_THRESHOLD_C, and all three Gilman strata.
#   B3  independent R <-> model_math at ~machine epsilon (two transcriptions of
#       the same arithmetic in the same engine have no excuse to differ).
#   B4  the sub-kernels obs_prob() does NOT expose -- phi0(T), phi(T,D), eta(D)
#       -- gated against the generated-quantities scalars Stan already emits.
#       This is why the .stan needs no new phi_pred output.
# ==============================================================================
source("model_math.R")
TOL_RR <- 1e-11   # R vs R; a loose tolerance here would mask a real formula difference
TOL_GQ <- 1e-7    # R vs Stan generated quantities

stopifnot(identical(MM_RAW_PARS, PARAM_NAMES))   # model_math must track parameters{}

par_bundles <- lapply(vecs, function(v) mm_pars(as.list(setNames(v, PARAM_NAMES)), T_ref = T_REF))

#' Synthetic covariate grid: structured per group rather than a full cross, so
#' every covariate a group actually READS is varied and none are wasted.
parity_grid <- function() {
  doses <- 10^seq(2, 9.7, length.out = 12)
  cops  <- c(1, 2.16, 38.11, 152.16, 500)          # incl. beyond any observed titre
  Ts    <- c(37.5, 38.0, 38.3, 39.0, 39.4, 40.5)   # incl. outside the study thresholds
  cols  <- c("likelihood_group", "dose_cfu", "CoP", "T_thresh", "gilman_stratum")
  ox <- expand.grid(likelihood_group = c("ox_fev", "ox_inf", "ox_inf_indiv", "ox_fevginf_indiv"),
                    dose_cfu = doses, CoP = cops, stringsAsFactors = FALSE)
  ox$T_thresh <- T_REF; ox$gilman_stratum <- 0L
  md <- expand.grid(likelihood_group = c("md_fev", "hornick_cond"), dose_cfu = doses,
                    T_thresh = Ts, gilman_stratum = 0:2, stringsAsFactors = FALSE)
  md$CoP <- 1
  mi <- data.frame(likelihood_group = "md_inf", dose_cfu = doses, CoP = 1,
                   T_thresh = T_REF, gilman_stratum = 0L, stringsAsFactors = FALSE)
  g <- rbind(ox[cols], md[cols], mi[cols])
  g$n <- 1L; g$y <- 0L
  g$obs_id <- sprintf("grid-%04d", seq_len(nrow(g)))
  g
}

ladder <- darton_phi0_ladder("dose_response_data.csv")
cases <- list(
  tier1 = obs,
  tier2 = attr(build_stan_data("dose_response_data.csv", priors, tier_col = "tier2_active"), "obs"),
  grid  = parity_grid()
)

for (nm in names(cases)) {
  rows <- cases[[nm]]
  sd_c <- stan_data_from_rows(rows, priors, ladder, T_ref = T_REF)
  gq_c <- mod$generate_quantities(fitted_params = truth, data = sd_c, seed = 1)
  p_s  <- posterior::as_draws_matrix(gq_c$draws("p_pred"))
  grp  <- as.integer(.GROUP_CODE[rows$likelihood_group])

  d_mm <- 0; d_rr <- 0
  for (di in seq_len(nrow(p_s))) {
    p_mm <- as.vector(mm_obs_prob(grp, rows$dose_cfu, rows$CoP, rows$T_thresh,
                                  rows$gilman_stratum, par_bundles[[di]]))
    d_mm <- max(d_mm, max(abs(p_s[di, ] - p_mm)))
    # obs_prob_R() now covers all 7 groups including ox_inf, so eta gets the same
    # THREE-implementation check (Stan / model_math / independent R) as every other
    # branch. It previously had only two, on exactly the branch Tier 2 turns on.
    pl <- as.list(vecs[[di]])
    r  <- rows; r$group <- grp
    p_rr <- vapply(seq_len(nrow(r)), function(i) obs_prob_R(r[i, ], pl), numeric(1))
    d_rr <- max(d_rr, max(abs(p_rr - p_mm)))
  }
  cat(sprintf("model_math parity [%s: %d rows x %d param vectors]\n", nm, nrow(rows), nrow(p_s)))
  cat(sprintf("  max |p_pred_stan - p_model_math| = %.3e\n", d_mm))
  cat(sprintf("  max |p_indep_R   - p_model_math| = %.3e\n", d_rr))
  if (!(d_mm < TOL_P))  failures <- c(failures, sprintf("B: Stan vs model_math (%s)", nm))
  if (!(d_rr < TOL_RR)) failures <- c(failures, sprintf("B: independent R vs model_math (%s)", nm))
}

# ---- B4: sub-kernels obs_prob() does not expose, vs Stan generated quantities --
gq_s  <- mod$generate_quantities(fitted_params = truth, data = sd, seed = 1)
gqv   <- function(v) as.vector(posterior::as_draws_matrix(gq_s$draws(v)))
p_all <- mm_pars(as.data.frame(do.call(rbind, lapply(vecs, function(v) setNames(v, PARAM_NAMES)))),
                 T_ref = T_REF)
delta <- p_all$delta
kern <- list(
  phi0_38_3        = list(gq = gqv("phi0_38_3"),        mm = as.vector(mm_phi0(38.3, p_all))),
  phi0_39_4        = list(gq = gqv("phi0_39_4"),        mm = as.vector(mm_phi0(39.4, p_all))),
  phi_hornick_1e3  = list(gq = gqv("phi_hornick_1e3"),  mm = as.vector(mm_phi_td(39.4, cbind(1e3 / delta), p_all))),
  phi_hornick_1e5  = list(gq = gqv("phi_hornick_1e5"),  mm = as.vector(mm_phi_td(39.4, cbind(1e5 / delta), p_all))),
  phi_hornick_1e9  = list(gq = gqv("phi_hornick_1e9"),  mm = as.vector(mm_phi_td(39.4, cbind(1e9 / delta), p_all))),
  eta_1e3          = list(gq = gqv("eta_1e3"),          mm = as.vector(mm_eta(1e3, p_all))),
  eta_1e4          = list(gq = gqv("eta_1e4"),          mm = as.vector(mm_eta(1e4, p_all))),
  p_inf_1e3_naive  = list(gq = gqv("p_inf_1e3_naive"),  mm = as.vector(mm_p_inf(1e3, 1, p_all))),
  p_inf_1e4_naive  = list(gq = gqv("p_inf_1e4_naive"),  mm = as.vector(mm_p_inf(1e4, 1, p_all))),
  p_fev_1e3_naive  = list(gq = gqv("p_fev_1e3_naive"),  mm = as.vector(mm_p_fev(1e3, 1, p_all))),
  p_fev_1e4_naive  = list(gq = gqv("p_fev_1e4_naive"),  mm = as.vector(mm_p_fev(1e4, 1, p_all))),
  delta_fold       = list(gq = gqv("delta_fold"),       mm = delta)
)
# RELATIVE difference: most of these are probabilities (O(1)) but delta_fold is
# O(1e3), where an absolute 1e-7 tolerance is below cmdstan's CSV output precision.
d_k <- vapply(kern, function(z) max(abs(z$gq - z$mm) / pmax(1, abs(z$gq))), numeric(1))
cat("sub-kernel parity vs Stan generated quantities (phi0, phi(T,D), eta, reference doses):\n")
cat(sprintf("  max relative over %d quantities = %.3e   (worst: %s)\n",
            length(d_k), max(d_k), names(which.max(d_k))))
if (!(max(d_k) < TOL_GQ)) failures <- c(failures, "B4: sub-kernels vs Stan generated quantities")

# ---- B5: the one quantity that exists nowhere in Stan -------------------------
# The implied peak-fever severity density -d(phi0)/dT used by figures_phi.R.
# Gate the analytic form against a central finite difference of mm_phi0().
Tg   <- seq(37, 41, by = 0.25)
d_fd <- max(abs(mm_phi0_density(Tg, p_all) -
                -(mm_phi0(Tg + 1e-5, p_all) - mm_phi0(Tg - 1e-5, p_all)) / 2e-5))
cat(sprintf("severity density -d(phi0)/dT vs finite difference: %.3e\n", d_fd))
if (!(d_fd < 1e-6)) failures <- c(failures, "B5: severity density vs finite difference")

# ---- Verdict -----------------------------------------------------------------
if (!length(failures)) {
  cat("PARITY PASS\n")
} else {
  cat("PARITY FAIL:\n"); for (f in failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
