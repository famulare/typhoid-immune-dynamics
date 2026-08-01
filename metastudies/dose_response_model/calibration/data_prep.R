#' Flat Stan data assembly for the dose-response calibration.
#'
#' Builds the unified per-observation layout the refactored typhoid_dose_response.stan
#' expects (one row per observation + a `group` code), merged with the prior
#' hyperparameter data list from priors.yaml. Shared by fit_dose_response.R,
#' simulate_recovery.R, and run_scenarios.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
if (!exists("calib_dir")) source("utils.R")

# likelihood_group -> Stan group code (must match obs_prob() in the .stan)
.GROUP_CODE <- c(ox_fev = 1L, ox_inf = 2L, md_fev = 3L, md_inf = 4L, hornick_cond = 5L,
                 ox_inf_indiv = 6L, ox_fevginf_indiv = 7L)

#' Fail loudly when supplied parameter names don't match the model's parameters{}.
#' generate_quantities(fitted_params=) requires *exactly* the model's parameter
#' block; a stale name list otherwise surfaces as a cryptic cmdstan
#' "Mismatch between model and fitted_parameters csv file" error. Call this before
#' any generate_quantities() that feeds a fitted-params matrix. Reads names from
#' the .stan via mod$variables(), so it works without a compiled binary.
#' @param mod a cmdstanr model object
#' @param provided character vector of parameter names being supplied
#' @return invisibly, the model's parameter names (so callers can reuse them)
assert_fitted_params_match <- function(mod, provided) {
  model_pars <- names(mod$variables()$parameters)
  missing <- setdiff(model_pars, provided)
  extra   <- setdiff(provided, model_pars)
  if (length(missing) || length(extra)) {
    stop(sprintf(paste0(
      "fitted_params do not match typhoid_dose_response.stan parameters{}:\n",
      "  missing (in model, not provided): %s\n",
      "  extra   (provided, not in model): %s\n",
      "  -> update the PARAM_NAMES / truth vector to match the model's parameter block."),
      if (length(missing)) paste(missing, collapse = ", ") else "<none>",
      if (length(extra))   paste(extra,   collapse = ", ") else "<none>"),
      call. = FALSE)
  }
  invisible(model_pars)
}

# Naive anti-Vi reference (VaccZyme EU/mL, <LLD imputed). CoP = anti-Vi / ref, so
# CoP=1 at naive. See tier1.5_plan.md / tier1_lab_notebook.md D1.
NAIVE_VI_REF <- 3.7

# +vaccine-terms (2026-07-31): which non-anti-Vi protection channel a subject sits in.
# 0 = none (Placebo, and every non-Darton-individual row via the harmless default).
# Both M01ZH09 and Ty21a protect via mechanisms OTHER than anti-Vi (Ty21a is
# Vi-negative; M01ZH09 did not raise anti-Vi IgG) -- see tier1.5_plan.md
# "+vaccine-terms". This does NOT mean their anti-Vi titre is unusable: both arms
# have real per-subject vi_igg_prechallenge measurements (M01ZH09 31/31, Ty21a
# 29/30) that feed the SAME shared CoP^gamma channel as Placebo; V_v is the
# ADDITIONAL residual protection their own titre does not explain.
VACCINE_ID <- c(Placebo = 0L, M01ZH09 = 1L, Ty21a = 2L)

# ---- cohort_id: PROVENANCE ONLY (added 2026-07-31) ---------------------------
# `cohort_id` records which observations come from the SAME GROUP OF VOLUNTEERS.
# It is NOT passed to Stan and NO likelihood term reads it. Its two uses are both
# post-fit / bookkeeping:
#   1. provenance -- shared-subject structure is recoverable from the data;
#   2. LOO unit grouping in run_scenarios.R::compute_loo_units(), so model comparison
#      does not count the same volunteers twice (Tier 1: 80 rows -> 46 units).
# It is also the precondition for any future hierarchical term
# (see cohort_random_effects_design.md).
#
# Why it cannot be derived from existing columns:
#   - `study` is the PAPER. Levine's four trials (1970-73) are one study but four
#     separate challenge cohorts, and their 25-55% fever spread is WITHIN Levine.
#   - `study x year` fails the other way: Hornick's rows are all year 1965 but are
#     five different volunteer groups, distinguished only by challenge dose.
#
# What it encodes:
#   - Levine Lev-F-k and Lev-I-k are the SAME men (e.g. 13/26 fever and 19/26
#     infection in trial 1); same for the Oxford fever/infection pairs.
#   - H-I-7 (30) and H-FgI-7 (28) are NESTED, not disjoint; the cascade
#     factorization already handles that correctly.
#
# Two documented imperfections -- do not treat cohort_id as exact:
#   - MD-HOR-1965-D5 (H-F-5, n=116) is itself a POOL of many challenges across
#     years. It is not a single cohort even in principle; the id is a placeholder.
#   - MD-GIL-1975-CTRL covers both the three H-strata (which partition all 64
#     controls) and Gil-I-ctrl (the trials-1&3 subset, 43 of 64). Those PARTIALLY
#     overlap rather than being disjoint or nested.

# ---- phi(T,D) definition-map inputs (C3) ------------------------------------
# T_ref = Oxford composite fever threshold (>=38 degC). phi0(T) is logit-centered here.
T_REF <- 38.0
# Per-study strict fever threshold (degC) for the Maryland fever obs. The dose-
# response is threshold-free; only the phi(T,D) definition map reads T. Hornick
# >=103F/24-36h ~ 39.4; Levine >=101F ~ 38.3; Gilman fever+culture ~ 38.3.
#
# PROVENANCE [2026-07-31]: Hornick and Levine are EXTRACTED; Gilman is an ASSUMPTION.
#   Levine 1976 p.426 defines the endpoint numerically: "acute illness with oral
#     temperature of >=101 F accompanied by isolation of S. typhi from blood or stool"
#     -> 38.3 C, extracted.
#   Gilman 1977 p.719 NEVER numerically defines fever in its case definition. It says
#     "For purposes of therapy, typhoid fever was defined as the presence of fever and
#     a blood or stool culture positive for S. typhi", then gives 39.4/38.3/37.8 C as
#     CHLORAMPHENICOL TREATMENT TRIGGERS (>103F 1 day, >101F 3 days, >100F 5 days),
#     not as the case definition. The 38.3 below is imported from Levine on the
#     assumption that the two Maryland programs used the same operational threshold.
#     It propagates into phi(T,D) for Gil-F-Hlo/Hhi/rest.
#   -> If that assumption matters, the sensitivity to run is Gilman in {37.8, 38.3, 39.4}.
STUDY_FEVER_THRESHOLD_C <- c(Hornick = 39.4, Levine = 38.3, Gilman = 38.3)
# Darton placebo temperature-ladder thresholds fed to the phi0(T) sub-likelihood.
LADDER_THRESHOLDS_C <- c(38.0, 38.5, 39.0)

# ---- psi-correction (Tier 2, tier2_plan.md; joint_inference_plan.md Sec 2.8) -------
# psi_def per md_inf row, by study: 0 = broad reference (Hornick stool-or-blood,
# psi=1), 1 = Levine any-time-stool (psi_stool), 2 = Gilman late shedding (psi_late =
# psi_stool*frac_late). Applies ONLY to group 4 (md_inf); every other group's psi_def
# is 0 (harmless default -- obs_prob()/mm_obs_prob() never read psi_def outside group 4).
PSI_DEF_BY_STUDY <- c(Hornick = 0L, Levine = 1L, Gilman = 2L)

#' Darton placebo stool-vs-broad cross-tab: identifies psi_stool*eta(18200) jointly
#' (a decoupled binomial sub-likelihood, same pattern as darton_phi0_ladder()).
#'
#' CORRECTION 2026-07-31 [Mike, caught mid-implementation; see tier2_plan.md]:
#' joint_inference_plan.md Sec 2.8 states "19 ~ binomial(26, psi_stool)" using the
#' MARGINAL totals (19 total stool_positive, 26 total bact_or_stool+ among the 30
#' Placebo subjects). Those do NOT nest at the subject level -- 4 subjects have
#' stool_positive=1 with bact_or_stool=0 AND bacteremia=0 (bact_or_stool is not a
#' simple OR of its visible columns in this extract). Uses the TRUE intersection
#' (stool_positive AND bact_or_stool), computed from the extract rather than
#' hardcoded, so a future correction to the extract is reflected automatically.
darton_psi_crosstab <- function(data_csv, dose_cfu = 18200) {
  s1 <- file.path(dirname(data_csv), "..", "analysis_data", "darton_individual_endpoints.csv")
  d  <- readr::read_csv(s1, show_col_types = FALSE) %>% filter(group == "Placebo")
  n_bos <- sum(d$bact_or_stool == 1L)
  y_stool_and_bos <- sum(d$stool_positive == 1L & d$bact_or_stool == 1L)
  list(psi_crosstab_y = as.integer(y_stool_and_bos),
       psi_crosstab_n = as.integer(n_bos),
       psi_crosstab_dose = as.numeric(dose_cfu))
}

#' Darton placebo temperature ladder: among the TD+ placebo subjects, how many
#' crossed each strict threshold. Pins phi0(T) as a decoupled binomial sub-model.
#' Single source of truth: the individual-endpoints extract (same as C1).
darton_phi0_ladder <- function(data_csv, thresholds = LADDER_THRESHOLDS_C) {
  s1 <- file.path(dirname(data_csv), "..", "analysis_data", "darton_individual_endpoints.csv")
  d  <- readr::read_csv(s1, show_col_types = FALSE) %>%
    filter(group == "Placebo", fever_td == 1)
  col <- function(t) paste0("fever_", sub("\\.", "_", format(t, trim = TRUE, nsmall = 1)))
  cols <- vapply(thresholds, col, character(1))          # 38.0 -> fever_38_0 ...
  # the extract uses fever_38 / fever_38_5 / fever_39 (no trailing _0); normalize
  cols <- sub("_0$", "", cols)
  missing <- setdiff(cols, names(d))
  if (length(missing)) stop("darton ladder: missing threshold columns: ",
                            paste(missing, collapse = ", "))
  list(
    N_ladder     = length(thresholds),
    ladder_T     = as.numeric(thresholds),
    ladder_count = as.integer(vapply(cols, function(c) sum(d[[c]]), numeric(1))),
    ladder_N     = nrow(d)
  )
}

# The grouped Darton placebo rows, which the per-subject cascade REPLACES. Both are
# the same 30 volunteers under different markers: D-F-plac is composite fever (TD),
# D-I-plac is stool shedding (tier2_active only). Under psi (joint_inference_plan
# Sec 2.8) the stool-vs-broad contrast is carried by a decoupled sub-likelihood on the
# S1 cross-tab, so neither grouped row is needed in the main likelihood.
DARTON_PLACEBO_GROUPED_OBS <- c("D-F-plac", "D-I-plac")

#' obs_id of the per-subject Darton rows, so a scenario's drop_obs/keep_obs can
#' address them (they exist only after individualization, so they are not in the CSV).
darton_indiv_obs_ids <- function(data_csv) {
  darton_placebo_individual_rows(data_csv)$obs_id
}

#' Tier 1.5 +cascade (issue #15): a Darton arm's subjects as individual n=1 rows,
#' decomposed into the proper CASCADE — an infection endpoint (bact_or_stool) for ALL
#' subjects (group ox_inf_indiv, P_inf), and a fever|infection endpoint (fever_td) for
#' the INFECTED subjects only (group ox_fevginf_indiv, P_fev|inf). This replaces the C1
#' composite-fever rows (which conflated the two layers); the product recovers the
#' composite while separating gamma_inf from gamma_fevginf via the per-subject titre
#' spread. Infection = bact_or_stool (broadest marker, includes bacteremia; treated
#' as true infection — no eta treatment-truncation correction at Tier 1). Single
#' source of truth: analysis_data/darton_individual_endpoints.csv.
#'
#' +vaccine-terms (2026-07-31): generalized from Placebo-only to any of the trial's
#' three arms. Subjects with a missing pre-challenge titre are DROPPED (cannot compute
#' CoP) rather than imputed -- affects exactly one Ty21a subject. Each arm's own
#' anti-Vi titre feeds the SAME CoP^gamma channel; `vaccine_id` (0/1/2) is read by
#' obs_prob() to apply the arm's ADDITIONAL non-anti-Vi protection factor V_v (0 for
#' Placebo, i.e. no additional effect -- see VACCINE_ID above).
#' @param arm one of names(VACCINE_ID): "Placebo", "M01ZH09", "Ty21a".
darton_arm_individual_rows <- function(data_csv, arm) {
  if (!arm %in% names(VACCINE_ID)) stop("darton_arm_individual_rows: unknown arm '", arm, "'")
  s1 <- file.path(dirname(data_csv), "..", "analysis_data", "darton_individual_endpoints.csv")
  d <- readr::read_csv(s1, show_col_types = FALSE) %>%
    filter(group == arm, !is.na(vi_igg_prechallenge)) %>%
    mutate(CoP = vi_igg_prechallenge / NAIVE_VI_REF)   # VaccZyme EU/mL, ref naive
  # "Placebo" keeps the PRE-EXISTING "plac" tag (obs_id "D-I-plac-...", cohort_id
  # "OX-DAR-2013-PLAC") that issue #15's original individualization already committed
  # to -- curve_specs.R's regex and every doc/plot referencing it depend on this exact
  # string. Deriving the tag generically from the arm name (as for the two NEW arms)
  # would have silently renamed it to "placebo" and broken all of that.
  tag <- if (arm == "Placebo") "plac" else tolower(gsub("[^A-Za-z0-9]", "", arm))
  cohort  <- sprintf("OX-DAR-2013-%s", toupper(tag))
  vid     <- unname(VACCINE_ID[arm])
  infection_rows <- d %>% transmute(
    obs_id = sprintf("D-I-%s-%s", tag, subject_id), study = "Darton",
    cohort_id = cohort,
    likelihood_group = "ox_inf_indiv", dose_cfu = 18200, n = 1L,
    y = as.integer(bact_or_stool), CoP = CoP, phi = 1.0, gilman_stratum = 0L,
    vaccine_id = vid)
  fevginf_rows <- d %>% filter(bact_or_stool == 1) %>% transmute(
    obs_id = sprintf("D-FgI-%s-%s", tag, subject_id), study = "Darton",
    cohort_id = cohort,   # all subjects in one arm are ONE challenge cohort
    likelihood_group = "ox_fevginf_indiv", dose_cfu = 18200, n = 1L,
    y = as.integer(fever_td), CoP = CoP, phi = 1.0, gilman_stratum = 0L,
    vaccine_id = vid)
  bind_rows(infection_rows, fevginf_rows)
}

#' Backward-compatible name: the Placebo arm only (issue #15's original scope).
darton_placebo_individual_rows <- function(data_csv) darton_arm_individual_rows(data_csv, "Placebo")

#' Apply nested prior overrides (for sensitivity scenarios), e.g.
#'   list(log10_delta = list(mu = 3.0))  ->  fixes the delta prior mean.
apply_prior_overrides <- function(priors, overrides = list()) {
  for (param in names(overrides)) {
    if (is.null(priors[[param]])) stop(sprintf("override for unknown prior '%s'", param))
    for (h in names(overrides[[param]])) priors[[param]][[h]] <- overrides[[param]][[h]]
  }
  priors
}

#' Assemble the flat Stan data list from an arbitrary covariate frame.
#'
#' Factored out of build_stan_data() so that test_obs_prob_parity.R can gate the
#' model over a SYNTHETIC covariate grid through the same production assembly
#' path the real fit uses -- a change to how stratum/CoP defaults are coerced is
#' then automatically reflected in the gate instead of drifting away from it.
#'
#' @param rows data frame with likelihood_group, dose_cfu, n, y, CoP,
#'   gilman_stratum, T_thresh (already resolved; see STUDY_FEVER_THRESHOLD_C).
#' @param priors Parsed priors (from load_priors()).
#' @param ladder Darton phi0 ladder list (from darton_phi0_ladder()).
#' @param psi_crosstab Darton psi cross-tab list (from darton_psi_crosstab()).
#' @param psi_active Tier 2 psi-correction gate (tier2_plan.md decision B): 0 forces
#'   psi=1 for every row regardless of psi_def and skips the crosstab likelihood term
#'   entirely -- the bit-identical-to-before behavior for t1-* tiers.
#' @return the flat Stan data list (no "obs" attribute; callers attach it).
stan_data_from_rows <- function(rows, priors, ladder, psi_crosstab = NULL,
                                T_ref = T_REF, prior_only = 0L, psi_active = 0L) {
  grp <- unname(.GROUP_CODE[rows$likelihood_group])
  if (anyNA(grp)) stop("unmapped likelihood_group: ",
                       paste(unique(rows$likelihood_group[is.na(grp)]), collapse = ", "))
  # vaccine_id (+vaccine-terms): 0 = no additional non-anti-Vi effect. Harmless default
  # for every row that isn't a Darton M01ZH09/Ty21a individual row, and for any caller
  # (e.g. the parity gate's synthetic grid) that predates this column entirely.
  vacc <- if ("vaccine_id" %in% names(rows)) rows$vaccine_id else 0L
  # psi_def (Tier 2, tier2_plan.md): 0 = broad/none. Harmless default for every row
  # that predates this column (e.g. the parity gate's synthetic grid).
  psid <- if ("psi_def" %in% names(rows)) rows$psi_def else rep(0L, nrow(rows))
  if (is.null(psi_crosstab))
    psi_crosstab <- list(psi_crosstab_y = 0L, psi_crosstab_n = 0L, psi_crosstab_dose = 18200)
  stan_data <- list(
    N_obs    = nrow(rows),
    group    = as.integer(grp),
    n        = as.integer(rows$n),
    y        = as.integer(rows$y),
    dose     = as.numeric(rows$dose_cfu),
    CoP      = ifelse(is.na(rows$CoP), 1.0, as.numeric(rows$CoP)),   # used by ox groups only
    stratum  = ifelse(is.na(rows$gilman_stratum), 0L, as.integer(rows$gilman_stratum)),
    T_thresh = as.numeric(rows$T_thresh),
    T_ref    = T_ref,
    vaccine_id = ifelse(is.na(vacc), 0L, as.integer(vacc)),
    psi_def  = ifelse(is.na(psid), 0L, as.integer(psid)),
    psi_active = as.integer(psi_active),
    prior_only = as.integer(prior_only)
  )
  stan_data <- c(stan_data, ladder, psi_crosstab, priors_to_stan_data(priors))
  stopifnot(!anyNA(unlist(stan_data[c("dose", "CoP", "n", "y", "group", "stratum",
                                      "T_thresh", "vaccine_id", "psi_def", "ladder_count")])))
  stan_data
}

#' @param data_csv Path to dose_response_data.csv.
#' @param priors Parsed priors (from load_priors()), optionally override-mutated.
#' @param tier_col Which activation column selects rows ("tier1_active" or "tier2_active").
#' @param prior_only 1 to skip the likelihood (prior predictive).
#' @param drop_obs Character vector of obs_id to exclude (data-filter sensitivities).
#' @param keep_obs If non-NULL, restrict to these obs_id (e.g. exclude-Oxford/Maryland).
#' @param include_vaccine_arms +vaccine-terms: also individualize Darton's M01ZH09 and
#'   Ty21a arms (same cascade shape as Placebo), each carrying its own vaccine_id so
#'   obs_prob() can apply its ADDITIONAL non-anti-Vi protection factor V_v. Requires
#'   individualize_darton = TRUE (there is no grouped representation of these arms in
#'   the active tiers -- the CSV's D-F-Ty21a/D-F-M01 rows are validation_only and
#'   inert regardless of this flag).
#' @param psi_active Tier 2 (tier2_plan.md): apply the psi infection-definition
#'   correction to group 4 (md_inf) rows (Levine/Gilman vs the Hornick broad
#'   reference) and add the Darton psi cross-tab decoupled likelihood term. FALSE
#'   (default) keeps t1-* tiers bit-identical to before this increment: psi_def is
#'   0 for every row and the crosstab term is skipped.
#' @return list: the flat Stan data + the prior `pr_*` scalars, plus an attribute
#'   "obs" carrying the selected data frame (obs_id, study, group, etc.) for plotting.
build_stan_data <- function(data_csv, priors,
                            tier_col = "tier1_active",
                            prior_only = 0L,
                            drop_obs = character(),
                            keep_obs = NULL,
                            individualize_darton = TRUE,
                            include_vaccine_arms = FALSE,
                            psi_active = FALSE) {
  if (include_vaccine_arms && !individualize_darton)
    stop("build_stan_data: include_vaccine_arms requires individualize_darton = TRUE",
        call. = FALSE)
  d <- readr::read_csv(data_csv, show_col_types = FALSE)
  dat <- d %>% filter(.data[[tier_col]] == 1)
  # The per-subject rows REPLACE **both** grouped Darton placebo rows: same 30
  # volunteers, nested markers. Dropping only D-F-plac left the grouped D-I-plac
  # shedding row (tier2_active, n=30 y=19, eta-corrected) in the likelihood beside the
  # 30 ox_inf_indiv rows (bact_or_stool, 26/30) for those same men -- the same
  # infection events entering twice as independent observations.
  #
  # Individualization runs BEFORE drop_obs/keep_obs so that (a) a scenario can address
  # the per-subject rows at all, and (b) drop_obs = "D-F-plac" no longer silently
  # removes all 56 Darton rows by making the guard below false.
  if (individualize_darton && any(DARTON_PLACEBO_GROUPED_OBS %in% dat$obs_id)) {
    dat <- dat %>% filter(!obs_id %in% DARTON_PLACEBO_GROUPED_OBS) %>%
      bind_rows(darton_placebo_individual_rows(data_csv))
  }
  # +vaccine-terms: M01ZH09/Ty21a have NO grouped representation in any active tier
  # (their CSV rows are validation_only, tier1_active = tier2_active = 0), so there is
  # nothing to filter out here -- just bind the individual rows.
  if (include_vaccine_arms) {
    dat <- dat %>% bind_rows(darton_arm_individual_rows(data_csv, "M01ZH09"),
                             darton_arm_individual_rows(data_csv, "Ty21a"))
  }
  if (!is.null(keep_obs)) dat <- dat %>% filter(obs_id %in% keep_obs)
  if (length(drop_obs))   dat <- dat %>% filter(!obs_id %in% drop_obs)
  dat <- dat %>% arrange(match(likelihood_group, names(.GROUP_CODE)), obs_id)
  if (!"vaccine_id" %in% names(dat)) dat$vaccine_id <- 0L
  dat$vaccine_id[is.na(dat$vaccine_id)] <- 0L

  # phi(T,D) reads a strict fever threshold only for Maryland fever obs; elsewhere
  # T_thresh is unused, set to T_REF as a harmless default.
  is_md_fever <- dat$likelihood_group %in% c("md_fev", "hornick_cond")
  dat$T_thresh <- ifelse(is_md_fever,
                         unname(STUDY_FEVER_THRESHOLD_C[dat$study]), T_REF)
  if (anyNA(dat$T_thresh)) stop("no fever threshold mapped for Maryland study: ",
                                paste(unique(dat$study[is_md_fever & is.na(dat$T_thresh)]), collapse = ", "))

  # psi_def (Tier 2, tier2_plan.md): only meaningful for group 4 (md_inf); 0
  # (broad/none) elsewhere, and 0 everywhere when psi_active is FALSE (so t1-*
  # tiers are bit-identical regardless of what this column would otherwise say).
  is_md_inf <- dat$likelihood_group == "md_inf"
  dat$psi_def <- 0L
  if (isTRUE(psi_active) && any(is_md_inf)) {
    dat$psi_def[is_md_inf] <- unname(PSI_DEF_BY_STUDY[dat$study[is_md_inf]])
    if (anyNA(dat$psi_def)) stop("no psi_def mapped for Maryland infection study: ",
                                 paste(unique(dat$study[is_md_inf & is.na(dat$psi_def)]), collapse = ", "))
  }

  stan_data <- stan_data_from_rows(dat, priors, darton_phi0_ladder(data_csv),
                                   psi_crosstab = darton_psi_crosstab(data_csv),
                                   T_ref = T_REF, prior_only = prior_only,
                                   psi_active = as.integer(isTRUE(psi_active)))
  attr(stan_data, "obs") <- dat %>%
    mutate(T_thresh = as.numeric(T_thresh),
           group = as.integer(unname(.GROUP_CODE[likelihood_group]))) %>%
    transmute(obs_id, study, cohort_id, likelihood_group, group,
              dose_cfu, n, y, obs_rate = y / n, CoP, phi, T_thresh, gilman_stratum,
              vaccine_id, psi_def)
  stan_data
}
