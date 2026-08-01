#' Draws-vectorized R transcription of the Stan model's `functions{}` block.
#'
#' THE ONE R implementation of the model math used by every figure. It exists so
#' that plotting code contains no model algebra: before this file, `.bp()` in
#' dose_response_curves.R was an UNGATED third transcription (alongside the .stan
#' and test_obs_prob_parity.R's obs_prob_R()), and each new panel re-derived the
#' Maryland fever product inline.
#'
#' Gating (test_obs_prob_parity.R):
#'   - mm_obs_prob() is compared against Stan p_pred over a synthetic covariate
#'     grid spanning all 7 groups, plus the real Tier-1 and Tier-2 rows.
#'   - obs_prob_R() in that test remains an INDEPENDENT transcription and must
#'     never be refactored to call these functions -- that independence is the
#'     whole point of the gate.
#'   - figures.R additionally asserts, at figure time, that mm_obs_prob() equals
#'     THIS fit's p_pred at its own posterior draws.
#'
#' Source of truth: typhoid_dose_response.stan lines 48-159. Every function below
#' is written to stay eyeball-comparable to its Stan counterpart.
#'
#' Shape contract: parameters are draws-indexed (length ndraws), covariates are
#' grid-indexed (length 1 or ngrid), every exported mm_* returns an
#' `ndraws x ngrid` matrix. To pass a PARAMETER where a covariate is expected
#' (e.g. CoP = CoP_susc), conform it explicitly with mm_by_draw().

# The model's parameters{} block, in declaration order. Must match
# names(mod$variables()$parameters); assert_fitted_params_match() enforces it.
MM_RAW_PARS <- c("log10_N50_inf", "d_fev", "alpha_inf", "alpha_fevginf",
                 "gamma_inf", "gamma_fevginf", "log10_delta", "pi_susc",
                 "CoP_imm", "CoP_susc", "phi0_a", "phi0_b", "eta_lo", "kappa",
                 "grand_overdispersion_rho", "log_V_M01ZH09", "log_V_Ty21a")

# ---- shape helpers -----------------------------------------------------------

#' Conform a draws-indexed quantity to ndraws x ngrid (constant along the grid).
mm_by_draw <- function(x, nd, ng) {
  if (is.matrix(x)) return(x)
  matrix(x, nd, ng)                      # length-nd recycles DOWN columns: [i,j] = x[i]
}

#' Conform a grid-indexed quantity to ndraws x ngrid (constant across draws).
mm_by_grid <- function(x, nd, ng) {
  if (is.matrix(x)) return(x)
  matrix(x, nd, ng, byrow = TRUE)        # length-ng fills ACROSS rows: [i,j] = x[j]
}

.mm_ng <- function(nd, ...) {
  lens <- vapply(list(...), function(x) if (is.matrix(x)) ncol(x) else length(x), integer(1))
  lens <- lens[lens > 1L]
  if (!length(lens)) return(1L)
  if (length(unique(lens)) != 1L)
    stop("model_math: covariates have inconsistent grid lengths: ",
         paste(unique(lens), collapse = ", "), call. = FALSE)
  unique(lens)
}

# ---- parameter bundle --------------------------------------------------------

#' Build a par bundle from anything coercible to a data.frame of the constrained
#' parameters. Derives N50_inf / N50_fevginf / delta in R so the Stan
#' reparameterization (log10_N50_fevginf = log10_N50_inf + d_fev) is mirrored
#' HERE and therefore gated, rather than silently assumed by callers.
#' @param raw data.frame / draws_df / named list of MM_RAW_PARS
#' @param T_ref phi0(T) logit centre; must match data_prep.R T_REF
mm_pars <- function(raw, T_ref = 38.0) {
  raw <- as.data.frame(raw)
  missing <- setdiff(MM_RAW_PARS, names(raw))
  if (length(missing))
    stop("mm_pars: missing parameters: ", paste(missing, collapse = ", "), call. = FALSE)
  p <- lapply(MM_RAW_PARS, function(nm) as.numeric(raw[[nm]]))
  names(p) <- MM_RAW_PARS
  p$log10_N50_fevginf <- p$log10_N50_inf + p$d_fev        # transformed parameters{}
  p$N50_inf           <- 10^p$log10_N50_inf
  p$N50_fevginf       <- 10^p$log10_N50_fevginf
  p$delta             <- 10^p$log10_delta
  # beta-binomial concentration (Step 2); mirrors the .stan's transformed parameters{}.
  p$grand_concentration_k <- (1 - p$grand_overdispersion_rho) / p$grand_overdispersion_rho
  # Per-vaccine non-anti-Vi protection factors, natural scale (+vaccine-terms).
  p$V_M01ZH09 <- exp(p$log_V_M01ZH09)
  p$V_Ty21a   <- exp(p$log_V_Ty21a)
  p$.ndraws <- length(p$log10_N50_inf)
  p$.draw   <- seq_len(p$.ndraws)
  p$T_ref   <- T_ref
  p
}

#' Par bundle straight from a CmdStanMCMC fit.
mm_draws <- function(fit, T_ref = 38.0) {
  d <- posterior::as_draws_df(fit$draws(MM_RAW_PARS))
  mm_pars(d, T_ref = T_ref)
}

#' Systematic (chain-balanced) thin to `ndraw` draws, preserving .draw ids.
#'
#' Every panel in a figure MUST use the same thin, so that a given .draw is the
#' same posterior sample everywhere and a reader can trace one draw across the
#' grid. Systematic rather than random because cmdstanr draws are stored
#' chain-major: seq() spreads the sample evenly over all chains.
mm_thin <- function(p, ndraw, seed = 1) {
  n <- p$.ndraws
  if (is.null(ndraw) || n <= ndraw) return(p)
  set.seed(seed)
  off <- sample.int(max(1L, n %/% ndraw), 1L) - 1L      # jitter the phase, not the spacing
  idx <- unique(pmin(n, round(seq(1, n, length.out = ndraw)) + off))
  q <- p
  for (nm in c(MM_RAW_PARS, "log10_N50_fevginf", "N50_inf", "N50_fevginf", "delta",
              "grand_concentration_k", "V_M01ZH09", "V_Ty21a", ".draw"))
    q[[nm]] <- p[[nm]][idx]
  q$.ndraws <- length(idx)
  q
}

# ---- kernels (mirror of typhoid_dose_response.stan functions{}) --------------

#' beta_poisson_vax(): P = 1 - (1+D_eff*(2^(1/alpha)-1)/N50)^(-alpha/(CoP^gamma * V))
#' (+vaccine-terms). V=1 recovers mm_bp() exactly. Elementwise and shape-free.
mm_bp_vax <- function(D_eff, N50, alpha, CoP, gamma, V) {
  scale <- (2^(1 / alpha) - 1) / N50
  1 - (1 + D_eff * scale)^(-alpha / (CoP^gamma * V))
}

#' beta_poisson(): P = 1 - (1 + D_eff*(2^(1/alpha)-1)/N50)^(-alpha/CoP^gamma)
#' A 1-line wrapper at V=1, mirroring the .stan's own beta_poisson()/beta_poisson_vax()
#' split -- every pre-existing caller of mm_bp() is untouched by +vaccine-terms.
mm_bp <- function(D_eff, N50, alpha, CoP, gamma) mm_bp_vax(D_eff, N50, alpha, CoP, gamma, 1)

#' phi0(T) = inv_logit(phi0_a - phi0_b*(T - T_ref)); low-dose definition sensitivity.
mm_phi0 <- function(T, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, T)
  plogis(mm_by_draw(p$phi0_a, nd, ng) -
         mm_by_draw(p$phi0_b, nd, ng) * (mm_by_grid(T, nd, ng) - p$T_ref))
}

#' Analytic severity density -d(phi0)/dT = phi0_b * phi0 * (1 - phi0).
#' The implied peak-temperature density among low-dose TD+ cases. Not in Stan;
#' gated against a central finite difference of mm_phi0().
mm_phi0_density <- function(T, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, T)
  ph <- mm_phi0(T, p)
  mm_by_draw(p$phi0_b, nd, ng) * ph * (1 - ph)
}

#' P_inf: beta-Poisson infection kernel.
mm_p_inf <- function(D_eff, CoP = 1, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, D_eff, CoP)
  mm_bp(mm_by_grid(D_eff, nd, ng), mm_by_draw(p$N50_inf, nd, ng),
        mm_by_draw(p$alpha_inf, nd, ng), mm_by_grid(CoP, nd, ng),
        mm_by_draw(p$gamma_inf, nd, ng))
}

#' P_fev|inf: beta-Poisson fever-given-infection kernel.
mm_p_fevginf <- function(D_eff, CoP = 1, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, D_eff, CoP)
  mm_bp(mm_by_grid(D_eff, nd, ng), mm_by_draw(p$N50_fevginf, nd, ng),
        mm_by_draw(p$alpha_fevginf, nd, ng), mm_by_grid(CoP, nd, ng),
        mm_by_draw(p$gamma_fevginf, nd, ng))
}

#' The cascade product P_inf x P_fev|inf (= obs_prob group 1, ox_fev).
mm_p_fev <- function(D_eff, CoP = 1, p) mm_p_inf(D_eff, CoP, p) * mm_p_fevginf(D_eff, CoP, p)

#' Per-observation additional protection factor V (+vaccine-terms): 0=none/Placebo -> 1,
#' 1=M01ZH09, 2=Ty21a. `vaccine_id` is grid-indexed (constant across draws) while
#' V_M01ZH09/V_Ty21a are draw-indexed, so this returns the full ndraws x length(vaccine_id)
#' selection directly (mm_by_grid/mm_by_draw's matrix passthrough then treats it as
#' already-conformed).
mm_vaccine_V <- function(vaccine_id, p) {
  nd <- p$.ndraws; ng <- length(vaccine_id)
  V <- matrix(1, nd, ng)
  j <- which(vaccine_id == 1L); if (length(j)) V[, j] <- p$V_M01ZH09
  j <- which(vaccine_id == 2L); if (length(j)) V[, j] <- p$V_Ty21a
  V
}

#' P_inf with an additional vaccine protection factor V (+vaccine-terms). V=1 recovers
#' mm_p_inf() exactly. Used only by mm_obs_prob()'s group 6 (ox_inf_indiv).
mm_p_inf_vax <- function(D_eff, CoP, V, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, D_eff, CoP, V)
  mm_bp_vax(mm_by_grid(D_eff, nd, ng), mm_by_draw(p$N50_inf, nd, ng),
           mm_by_draw(p$alpha_inf, nd, ng), mm_by_grid(CoP, nd, ng),
           mm_by_draw(p$gamma_inf, nd, ng), mm_by_grid(V, nd, ng))
}

#' P_fev|inf with an additional vaccine protection factor V (+vaccine-terms). V=1
#' recovers mm_p_fevginf() exactly. Used only by mm_obs_prob()'s group 7 (ox_fevginf_indiv).
mm_p_fevginf_vax <- function(D_eff, CoP, V, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, D_eff, CoP, V)
  mm_bp_vax(mm_by_grid(D_eff, nd, ng), mm_by_draw(p$N50_fevginf, nd, ng),
           mm_by_draw(p$alpha_fevginf, nd, ng), mm_by_grid(CoP, nd, ng),
           mm_by_draw(p$gamma_fevginf, nd, ng), mm_by_grid(V, nd, ng))
}

#' One Maryland mixture component evaluated at CoP_susc or CoP_imm.
#' layer: "inf" = P_inf, "fevginf" = P_fev|inf, "fev" = their product.
mm_md_component <- function(D_eff, p, layer = c("inf", "fevginf", "fev"),
                            which = c("susc", "imm")) {
  layer <- match.arg(layer); which <- match.arg(which)
  nd <- p$.ndraws; ng <- .mm_ng(nd, D_eff)
  CoP <- mm_by_draw(if (which == "susc") p$CoP_susc else p$CoP_imm, nd, ng)
  switch(layer, inf = mm_p_inf(D_eff, CoP, p),
                fevginf = mm_p_fevginf(D_eff, CoP, p),
                fev = mm_p_fev(D_eff, CoP, p))
}

#' Maryland latent-immunity mixture: pi_susc*P(CoP_susc) + (1-pi_susc)*P(CoP_imm).
#'
#' NOTE the "fev" layer mixes the PRODUCTS (pi*p_fev_susc + (1-pi)*p_fev_imm, as
#' written inline in obs_prob group 3), which is NOT the product of the mixtures.
#' Getting this wrong is the easiest way to silently mis-draw every Maryland panel.
mm_md_mix <- function(D_eff, p, layer = c("inf", "fevginf", "fev")) {
  layer <- match.arg(layer)
  nd <- p$.ndraws; ng <- .mm_ng(nd, D_eff)
  pi_s <- mm_by_draw(p$pi_susc, nd, ng)
  pi_s        * mm_md_component(D_eff, p, layer, "susc") +
  (1 - pi_s)  * mm_md_component(D_eff, p, layer, "imm")
}

#' phi(T, D_eff) = phi0(T) + (1 - phi0(T)) * P_fev_naive(D_eff).
#' P_fev_naive is the CoP=1 cascade so immunity is not double-counted; the
#' dose-lift exponent beta_phi is PINNED to 1 in the .stan (2026-07-15).
mm_phi_td <- function(T, D_eff, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, T, D_eff)
  p_fev_naive <- mm_p_fev(mm_by_grid(D_eff, nd, ng), 1, p)
  phi0 <- mm_phi0(mm_by_grid(T, nd, ng), p)
  phi0 + (1 - phi0) * p_fev_naive
}

#' eta(D_eff) shedding-detection probability. INERT at Tier 1 (eta_lo/kappa are
#' prior-only; every ox_inf row has tier1_active = 0).
mm_eta <- function(D_eff, p) {
  nd <- p$.ndraws; ng <- .mm_ng(nd, D_eff)
  eta_lo <- mm_by_draw(p$eta_lo, nd, ng)
  eta_lo + (1 - eta_lo) * exp(-mm_by_draw(p$kappa, nd, ng) *
                              mm_by_grid(D_eff, nd, ng) / mm_by_draw(p$N50_inf, nd, ng))
}

# ---- the mirror: obs_prob() dispatch -----------------------------------------

#' Per-observation binomial success probability, mirroring obs_prob() exactly.
#' Covariates are per-observation (grid-indexed, length 1 or N_obs).
#' @param vaccine_id 0=none/Placebo, 1=M01ZH09, 2=Ty21a (+vaccine-terms; groups 6/7 only).
#'   Kept LAST (after `p`, with a default) so every pre-existing positional call site
#'   -- `mm_obs_prob(group, dose, CoP, T_thresh, stratum, p)` -- is unaffected.
#' @return ndraws x N_obs matrix, column order = the order of `group`
mm_obs_prob <- function(group, dose, CoP, T_thresh, stratum, p, vaccine_id = 0L) {
  nd <- p$.ndraws
  ng <- max(length(group), length(dose), length(CoP), length(T_thresh), length(stratum),
           length(vaccine_id))
  rep_to <- function(x) if (length(x) == 1L) rep(x, ng) else x
  group <- as.integer(rep_to(group)); dose <- rep_to(dose); CoP <- rep_to(CoP)
  T_thresh <- rep_to(T_thresh); stratum <- as.integer(rep_to(stratum))
  vaccine_id <- as.integer(rep_to(vaccine_id))
  stopifnot(length(group) == ng, length(dose) == ng, length(CoP) == ng,
            length(T_thresh) == ng, length(stratum) == ng, length(vaccine_id) == ng)

  out <- matrix(NA_real_, nd, ng)
  D_ox <- mm_by_grid(dose, nd, ng)                                # delta = 1
  D_md <- D_ox / mm_by_draw(p$delta, nd, ng)                      # dose / delta
  CoPm <- mm_by_grid(CoP, nd, ng)
  sel  <- function(m, j) m[, j, drop = FALSE]

  g <- group
  j <- which(g == 1L)                                             # ox_fev
  if (length(j)) out[, j] <- mm_p_fev(sel(D_ox, j), sel(CoPm, j), p)
  j <- which(g == 2L)                                             # ox_inf (eta x P_inf)
  if (length(j)) out[, j] <- mm_eta(sel(D_ox, j), p) * mm_p_inf(sel(D_ox, j), sel(CoPm, j), p)
  j <- which(g == 6L)                                             # ox_inf_indiv (+vaccine V)
  if (length(j)) out[, j] <- mm_p_inf_vax(sel(D_ox, j), sel(CoPm, j),
                                          mm_vaccine_V(vaccine_id[j], p), p)
  j <- which(g == 7L)                                             # ox_fevginf_indiv (+vaccine V)
  if (length(j)) out[, j] <- mm_p_fevginf_vax(sel(D_ox, j), sel(CoPm, j),
                                              mm_vaccine_V(vaccine_id[j], p), p)
  j <- which(g == 4L)                                             # md_inf
  if (length(j)) out[, j] <- mm_md_mix(sel(D_md, j), p, "inf")

  j <- which(g == 3L)                                             # md_fev
  if (length(j)) {
    D <- sel(D_md, j)
    phi <- mm_phi_td(T_thresh[j], D, p)
    st <- stratum[j]
    pf <- matrix(NA_real_, nd, length(j))
    k <- which(st == 1L); if (length(k)) pf[, k] <- mm_md_component(D[, k, drop = FALSE], p, "fev", "susc")
    k <- which(st == 2L); if (length(k)) pf[, k] <- mm_md_component(D[, k, drop = FALSE], p, "fev", "imm")
    k <- which(!st %in% c(1L, 2L)); if (length(k)) pf[, k] <- mm_md_mix(D[, k, drop = FALSE], p, "fev")
    out[, j] <- phi * pf
  }

  j <- which(g == 5L)                                             # hornick_cond
  if (length(j)) {
    D <- sel(D_md, j)
    phi <- mm_phi_td(T_thresh[j], D, p)
    pc <- phi * mm_md_mix(D, p, "fev") / mm_md_mix(D, p, "inf")
    out[, j] <- pmin(pmax(pc, 1e-12), 1 - 1e-12)                  # guard the division
  }

  if (anyNA(out)) stop("mm_obs_prob: unhandled group code(s): ",
                       paste(sort(unique(group[apply(is.na(out), 2, any)])), collapse = ", "),
                       call. = FALSE)
  out
}

# ---- summarizing draws -------------------------------------------------------

#' Pointwise quantiles of an ndraws x ngrid matrix.
mm_quantiles <- function(mat, x, probs = c(0.05, 0.5, 0.95)) {
  q <- matrixStats::colQuantiles(mat, probs = probs, useNames = FALSE)
  data.frame(x = x, lo = q[, 1], med = q[, 2], hi = q[, 3])
}

#' Long-format draw curves for spaghetti layers. `p` supplies the .draw ids so
#' the same draw is identifiable across every panel of a figure.
mm_spaghetti <- function(mat, x, p) {
  data.frame(.draw = rep(p$.draw, times = length(x)),
             x     = rep(x, each = nrow(mat)),
             y     = as.vector(mat))
}
