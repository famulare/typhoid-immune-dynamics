#' Declarative definition of the comparable dataset sub-groupings, and the
#' assembly of model curves for each.
#'
#' One grouping = one column of the permutation grid = a set of observations that
#' share an era/vehicle, a fever-threshold definition, an infection definition,
#' and an immunity setting -- REGARDLESS of study. Waddington and Gibani share the
#' `ox_naive` column for exactly that reason: same composite TD definition, same
#' CoP = 1, different study (study is the point colour).
#'
#' Design constraints:
#'   - Membership is a PREDICATE over the obs tibble, never a column in
#'     dose_response_data.csv, because the 56 Darton per-subject rows do not exist
#'     in that CSV at all -- data_prep.R synthesizes them in
#'     darton_placebo_individual_rows(). A CSV column would silently miss exactly
#'     the rows the cascade decomposition is about.
#'   - T_thresh is DERIVED from data_prep.R's STUDY_FEVER_THRESHOLD_C / T_REF, and
#'     Stan group codes from .GROUP_CODE. Nothing already defined elsewhere is
#'     re-typed here, so the spec cannot drift from the data assembly.
#'   - validate_curve_specs() hard-errors if any observation matches != 1 spec, or
#'     if a spec's asserted covariates disagree with what build_stan_data() built.
#'     Adding a data row without a spec is a figure-time failure, not a silently
#'     dropped point.

suppressPackageStartupMessages({library(dplyr)})
if (!exists("mm_pars")) source("model_math.R")

# ---- rows of the grid --------------------------------------------------------
# Five quantities per grouping. Rows 1-4 are the factors of the model; row 5 is
# the full multiplicative factorization overlaid on one axis.
CURVE_ROWS <- tibble::tribble(
  ~row_key,    ~row_label,
  "p_inf",     "P(infection | D)",
  "p_inf_obs", "eta(D) x P(inf | D)",
  "p_fev",     "P(fever | D)",
  "p_fevginf", "P(fever | inf, D)",
  "phi",       "phi(T, D)",
  "cascade",   "cascade"
)

# Which ROW each observation's data points belong on. Single source for the join,
# so a point can never land on a curve it is not an observation of.
# ox_inf is the eta-corrected SHEDDING endpoint, so its points belong on the
# observed-scale p_inf_obs row, not on latent P(infection|D). Putting them on p_inf
# drew the data below a curve it is not an observation of.
.LG_ROW <- c(ox_fev = "p_fev", ox_inf = "p_inf_obs", md_fev = "p_fev", md_inf = "p_inf",
             hornick_cond = "p_fevginf", ox_inf_indiv = "p_inf",
             ox_fevginf_indiv = "p_fevginf")

# The five factors drawn on the cascade row, in multiplication order.
CASCADE_FACTORS <- c("P_inf", "P_fev|inf", "P_inf x P_fev|inf", "phi(T,D)", "observed-scale P(fever)")

#' The comparable dataset sub-groupings (columns of the grid).
#'
#' cop_mode drives which immunity setting the model curve is evaluated at:
#'   "fixed"      - a known group-average anti-Vi CoP (Oxford)
#'   "individual" - per-subject CoP; the dose-axis curve uses the cohort median
#'                  and a second ribbon spans the observed CoP range
#'   "mixture"    - the Maryland latent pi_susc/CoP_susc/CoP_imm mixture
#'   "stratum"    - a single Maryland mixture component (Gilman H-antibody strata)
curve_specs <- function(era_dose_range = list(oxford = c(1e2, 1e6),
                                              maryland = c(1e2, 1e10))) {
  s <- tibble::tribble(
    ~col_key,         ~study_set,              ~obs_match,               ~arm_label,                  ~cop_mode,    ~CoP,   ~stratum, ~thr_study, ~inf_def,
    "ox_naive",       "Waddington + Gibani",   "^(W-[FI]-|G20-)",        "naive",                     "fixed",       1.00,  0L,       NA,         "stool shedding",
    "ox_darton_plac", "Darton",                "^D-(I|FgI)-plac-",       "placebo, per-subject anti-Vi", "individual", NA,   0L,       NA,         "bacteraemia or stool",
    "ox_darton_plac_grp", "Darton",            "^D-[FI]-plac$",          "placebo, cohort GMT anti-Vi", "fixed",      1.98,  0L,       NA,         "stool shedding (grouped)",
    "ox_jin_ctrl",    "Jin",                   "^J-[FI]-ctrl$",          "control",                   "fixed",       2.16,  0L,       NA,         "stool shedding",
    "ox_jin_vips",    "Jin",                   "^J-[FI]-ViPS$",          "Vi-PS vaccinated",          "fixed",      38.11,  0L,       NA,         "stool shedding",
    "ox_jin_vitt",    "Jin",                   "^J-[FI]-ViTT$",          "Vi-TT vaccinated",          "fixed",     152.16,  0L,       NA,         "stool shedding",
    "md_hornick",     "Hornick",               "^H-(F|I|FgI)-[0-9]",     "latent mixture",            "mixture",       NA,  0L,       "Hornick",  "stool/blood culture (Table 2)",
    "md_levine",      "Levine",                "^Lev-",                  "latent mixture",            "mixture",       NA,  0L,       "Levine",   "any-time stool positive",
    "md_gilman_mix",  "Gilman",                "^Gil-(F-rest|I-ctrl)$",  "no H-Ab data (mixture)",    "mixture",       NA,  0L,       "Gilman",   "late shedding (4-30 d)",
    "md_gilman_susc", "Gilman",                "^Gil-F-Hlo$",            "H-Ab <1:20 (susceptible)",  "stratum",       NA,  1L,       "Gilman",   "fever + culture",
    "md_gilman_imm",  "Gilman",                "^Gil-F-Hhi$",            "H-Ab >=1:20 (immune)",      "stratum",       NA,  2L,       "Gilman",   "fever + culture"
  )
  # thr_study (fever-threshold lookup; a column may hold several studies) is a tribble
  # column, not a parallel named vector -- a new col_key can no longer be added to one
  # and silently forgotten in the other.
  s %>% mutate(
    era      = if_else(grepl("^md_", col_key), "maryland", "oxford"),
    vehicle  = if_else(era == "maryland", "milk", "bicarbonate"),
    # DERIVED from data_prep.R constants -- never re-typed, so it cannot drift
    T_thresh = unname(if_else(era == "maryland",
                              STUDY_FEVER_THRESHOLD_C[thr_study], T_REF)),
    dose_lo  = vapply(era, function(e) era_dose_range[[e]][1], numeric(1)),
    dose_hi  = vapply(era, function(e) era_dose_range[[e]][2], numeric(1)),
    col_label = sprintf("%s (%s)\n%s\nfever T>=%.1fC | inf: %s",
                        study_set, vehicle, arm_label, T_thresh, inf_def))
}

#' Assign each observation to exactly one grouping, and verify the spec's asserted
#' covariates against what build_stan_data() actually produced. This is the sync
#' mechanism: it fails loudly rather than dropping points.
validate_curve_specs <- function(specs, obs) {
  hit <- lapply(obs$obs_id, function(id)
    specs$col_key[vapply(specs$obs_match, grepl, logical(1), x = id)])
  n <- lengths(hit)
  if (any(n != 1L))
    stop("curve_specs: ", sum(n != 1L), " observation(s) match != 1 grouping:\n  ",
         paste(sprintf("%s (%d matches)", obs$obs_id[n != 1L], n[n != 1L]), collapse = "\n  "),
         "\n  -> add or fix a spec in curve_specs.R", call. = FALSE)
  chk <- obs %>% mutate(col_key = unlist(hit)) %>% left_join(specs, by = "col_key")
  bad <- chk %>% filter(
    (likelihood_group %in% c("md_fev", "hornick_cond") & abs(T_thresh.x - T_thresh.y) > 1e-9) |
    (cop_mode == "fixed"   & abs(CoP.x - CoP.y) > 1e-6) |
    (cop_mode == "stratum" & gilman_stratum != stratum))
  if (nrow(bad))
    stop("curve_specs disagrees with stan_data for: ",
         paste(bad$obs_id, collapse = ", "), call. = FALSE)
  chk %>% mutate(row_key = unname(.LG_ROW[likelihood_group]),
                 T_thresh = T_thresh.x, CoP = CoP.x) %>%
    select(-T_thresh.x, -T_thresh.y, -CoP.x, -CoP.y)
}

#' Model curve for one (quantity, grouping) cell.
#'
#' Each branch below is the corresponding obs_prob() branch: the "mixture" p_fev
#' is group 3 stratum 0, "stratum" p_fev is group 3 stratum 1/2, "fixed" p_fev is
#' group 1, and the "mixture" p_fevginf is group 5 (hornick_cond) verbatim --
#' including the fact that it does NOT factor, because the mixture doesn't.
#'
#' @param dose raw challenge dose (CFU) in the column's own vehicle frame
#' @return ndraws x length(dose) matrix
mm_curve <- function(quantity, dose, spec, p, cop_override = NULL) {
  nd <- p$.ndraws; ng <- length(dose)
  D  <- mm_by_grid(dose, nd, ng)
  De <- if (spec$era == "maryland") D / mm_by_draw(p$delta, nd, ng) else D
  phi <- if (spec$era == "maryland") mm_phi_td(spec$T_thresh, De, p) else
         mm_phi_td(p$T_ref, De, p)          # Oxford: computed but NOT applied (see below)
  is_mix <- spec$cop_mode == "mixture"
  CoP <- if (!is.null(cop_override)) cop_override else switch(
    spec$cop_mode,
    fixed      = spec$CoP,
    individual = spec$CoP,                   # supplied by the caller (cohort median)
    stratum    = mm_by_draw(if (spec$stratum == 1L) p$CoP_susc else p$CoP_imm, nd, ng),
    mixture    = NULL)

  # phi multiplies the observed-scale Maryland fever quantities only. Oxford's
  # composite TD endpoint IS the reference definition, so phi == 1 there by
  # construction of T_ref -- the phi row for an Oxford column is drawn as
  # "what fraction of TD+ would cross strict 38.0", not as a likelihood factor.
  phi_obs <- if (spec$era == "maryland") phi else matrix(1, nd, ng)

  # eta is the Oxford shedding detection/truncation factor and multiplies P_inf for
  # group 2 only. It is a SEPARATE row rather than a modification of p_inf: p_inf means
  # latent P(infection|D) and is correct as-is for the Darton/Waddington columns.
  switch(quantity,
    phi       = phi,
    p_inf     = if (is_mix) mm_md_mix(De, p, "inf") else mm_p_inf(De, CoP, p),
    p_inf_obs = mm_eta(De, p) * (if (is_mix) mm_md_mix(De, p, "inf")
                                 else mm_p_inf(De, CoP, p)),
    p_fevginf = if (is_mix) phi_obs * mm_md_mix(De, p, "fev") / mm_md_mix(De, p, "inf")
                else phi_obs * mm_p_fevginf(De, CoP, p),
    p_fev     = if (is_mix) phi_obs * mm_md_mix(De, p, "fev")
                else phi_obs * mm_p_fev(De, CoP, p),
    stop("mm_curve: unknown quantity '", quantity, "'"))
}

#' The five cascade factors for one grouping, as a named list of matrices.
mm_cascade <- function(dose, spec, p, cop_override = NULL) {
  nd <- p$.ndraws; ng <- length(dose)
  D  <- mm_by_grid(dose, nd, ng)
  De <- if (spec$era == "maryland") D / mm_by_draw(p$delta, nd, ng) else D
  is_mix <- spec$cop_mode == "mixture"
  CoP <- if (!is.null(cop_override)) cop_override else switch(
    spec$cop_mode,
    fixed = spec$CoP, individual = spec$CoP,
    stratum = mm_by_draw(if (spec$stratum == 1L) p$CoP_susc else p$CoP_imm, nd, ng),
    mixture = NULL)
  p_i <- if (is_mix) mm_md_mix(De, p, "inf")     else mm_p_inf(De, CoP, p)
  p_c <- if (is_mix) mm_md_mix(De, p, "fevginf") else mm_p_fevginf(De, CoP, p)
  p_f <- if (is_mix) mm_md_mix(De, p, "fev")     else mm_p_fev(De, CoP, p)
  phi <- if (spec$era == "maryland") mm_phi_td(spec$T_thresh, De, p) else matrix(1, nd, ng)
  out <- list(p_i, p_c, p_f, phi, phi * p_f)
  names(out) <- CASCADE_FACTORS
  # For a mixture column, P_inf x P_fev|inf is NOT the mixture fever probability
  # (mixing does not commute with multiplication); drop the misleading product line.
  if (is_mix) out[["P_inf x P_fev|inf"]] <- NULL
  # Oxford: phi == 1 in the likelihood by construction of T_ref, so the phi line is
  # a flat 1 and the observed-scale curve is identical to the product. Drawing both
  # is redundant AND hides the product line under the observed-scale one. The
  # definition map itself is still shown, on its own row.
  if (spec$era != "maryland") {
    out[["phi(T,D)"]] <- NULL
    out[["observed-scale P(fever)"]] <- NULL
  }
  out
}
