#' Figure-math gate: exercise mm_curve()/mm_cascade() over every (row x grouping)
#' cell at prior draws, with NO fit and NO MCMC.
#'
#' Why this exists separately from test_obs_prob_parity.R: that gate checks the
#' quantities the LIKELIHOOD uses (obs_prob / p_pred). The figure suite also draws
#' quantities the likelihood never evaluates -- latent P_inf on an Oxford column, the
#' observed-scale eta(D) x P_inf row, the cascade factors -- and those had no gate at
#' all. A `t2-*` grid is currently unfittable (both rungs are blocked), so without
#' this the eta curve branch would be entirely unexercised until someone unblocks a
#' tier and discovers it at figure time.
#'
#' Run:  Rscript test_curve_math.R   (exits non-zero on failure)

suppressPackageStartupMessages({library(dplyr)})
if (!exists("calib_dir")) source("utils.R")
setwd(calib_dir())
source("priors.R"); source("data_prep.R"); source("model_math.R"); source("curve_specs.R")

NDRAW <- 50L
SEED  <- 7L
failures <- character()
chk <- function(cond, msg) if (!isTRUE(cond)) failures <<- c(failures, msg)

priors <- load_priors()
set.seed(SEED)
raw <- vapply(MM_RAW_PARS, function(p) sample_prior(priors, p, NDRAW), numeric(NDRAW))
p <- mm_pars(as.data.frame(raw), T_ref = T_REF)

specs <- curve_specs()
dose  <- 10^seq(2, 9, length.out = 25)
rows  <- setdiff(CURVE_ROWS$row_key, "cascade")

for (i in seq_len(nrow(specs))) {
  sp  <- specs[i, ]
  cop <- if (sp$cop_mode == "individual") 1.98 else NULL   # cohort median stand-in
  m <- list()
  for (rk in rows) {
    x <- mm_curve(rk, dose, sp, p, cop_override = cop)
    m[[rk]] <- x
    chk(identical(dim(x), c(NDRAW, length(dose))),
        sprintf("%s/%s: shape %s, expected %dx%d", sp$col_key, rk,
                paste(dim(x), collapse = "x"), NDRAW, length(dose)))
    chk(all(is.finite(x)), sprintf("%s/%s: non-finite values", sp$col_key, rk))
    chk(all(x >= -1e-12 & x <= 1 + 1e-12),
        sprintf("%s/%s: outside [0,1] (min %.3g max %.3g)", sp$col_key, rk,
                min(x), max(x)))
  }
  # The invariant the eta row exists to express: eta <= 1, so the observed-scale
  # shedding curve can never sit above latent P(infection|D). If this inverts, the
  # figure would show data above a curve it is an observation of.
  chk(all(m$p_inf_obs <= m$p_inf + 1e-12),
      sprintf("%s: p_inf_obs > p_inf (eta > 1?)", sp$col_key))
  # The cascade identity, which must hold in EVERY cop_mode: the marginal is the
  # product of the infection curve and the conditional. For the Maryland mixture it
  # holds because p_fevginf is defined as phi*mix_fev/mix_inf, so the mixture cancels.
  # This is the check that would catch a wrong CoP or a dropped phi in one branch.
  chk(max(abs(m$p_fev - m$p_inf * m$p_fevginf)) < 1e-9,
      sprintf("%s: p_fev != p_inf * p_fevginf (max dev %.3g)", sp$col_key,
              max(abs(m$p_fev - m$p_inf * m$p_fevginf))))
  # Fever is nested inside infection.
  chk(all(m$p_fev <= m$p_inf + 1e-9),
      sprintf("%s: p_fev > p_inf (fever is nested in infection)", sp$col_key))
  # Note phi is NOT asserted to be 1 on Oxford columns: phi_obs (the likelihood
  # factor) is 1 there, but the drawn phi row is the diagnostic "what fraction of
  # TD+ would cross strict T_ref", which is a genuine dose-dependent quantity.

  fac <- mm_cascade(dose, sp, p, cop_override = cop)
  for (nm in names(fac)) {
    chk(all(is.finite(fac[[nm]])), sprintf("%s/cascade/%s: non-finite", sp$col_key, nm))
    chk(all(fac[[nm]] >= -1e-12 & fac[[nm]] <= 1 + 1e-12),
        sprintf("%s/cascade/%s: outside [0,1]", sp$col_key, nm))
  }
}

cat(sprintf("curve math over %d groupings x %d rows x %d prior draws\n",
            nrow(specs), length(rows), NDRAW))
if (length(failures)) {
  cat("CURVE MATH FAIL:\n"); cat(paste0("  - ", failures, "\n"), sep = "")
  quit(status = 1)
}
cat("CURVE MATH PASS\n")
