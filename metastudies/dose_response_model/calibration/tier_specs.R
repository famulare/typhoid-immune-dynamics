#' The tier ladder, declaratively. THE single source of truth for what a "tier" is.
#'
#' Why this file exists: `build_stan_data()` defaults to `individualize_darton = TRUE`
#' and no driver ever overrode it, so every fit since +cascade was the 80-observation
#' configuration while being labelled "Tier 1", reported as 25 observations, and
#' written to results/tier1/. One undocumented default was the whole ambiguity.
#'
#' The ladder is TWO-DIMENSIONAL, not a line:
#'   row set                 tier_col in {tier1_active, tier2_active}
#'   Darton representation   individualize_darton in {FALSE grouped, TRUE per-subject}
#' Naming both axes in the key is why results/<key>__<stage>/ cannot misdescribe
#' itself. A third axis -- the parameter set -- moves at every increment, so it is
#' the run dir's STAGE token, derived from the .stan and asserted, never typed.
#'
#' See TIER_LOCK.md (normative) and TIER_LADDER.md (generated counts).

if (!exists("calib_dir")) source("utils.R")
if (!exists("MM_RAW_PARS")) source("model_math.R")

# ---- model stage tokens -------------------------------------------------------

#' Parameter increments, in the order they entered the model. A stage token is the
#' hyphen-joined list of increments ACTIVE for a given tier: an increment counts
#' only when the .stan declares all of its parameters, and (for eta) only when the
#' tier's data actually contain the rows that make it enter the likelihood.
#'
#' eta is the reason `needs_group2` exists: eta_lo/kappa have been declared in the
#' .stan since before Tier 1, but they appear in exactly one likelihood branch
#' (group 2, ox_inf). Declared-but-unreached is not a model stage.
STAGE_INCREMENTS <- list(
  list(token = "phi", pars = c("phi0_a", "phi0_b"),          needs_group2 = FALSE),
  list(token = "rho", pars = "grand_overdispersion_rho",     needs_group2 = FALSE),
  list(token = "eta", pars = c("eta_lo", "kappa"),           needs_group2 = TRUE),
  list(token = "psi", pars = "psi_stool",                    needs_group2 = FALSE)
)

#' Derive the stage token from the model's parameter block and the tier's data.
#' @param model_pars character vector, names(mod$variables()$parameters).
#' @param has_group2 does this tier's data contain any ox_inf (group 2) row?
derive_stage_token <- function(model_pars, has_group2) {
  toks <- vapply(STAGE_INCREMENTS, function(inc) {
    active <- all(inc$pars %in% model_pars) && (!inc$needs_group2 || has_group2)
    if (active) inc$token else NA_character_
  }, character(1))
  toks <- toks[!is.na(toks)]
  if (!length(toks)) stop("no model stage increments active -- is this the right .stan?",
                          call. = FALSE)
  paste(toks, collapse = "-")
}

# ---- the registry -------------------------------------------------------------

#' Keys carry hyphens, so access is TIER_SPECS[["t1-indiv"]].
#'
#' `expect` is a hard assertion, not documentation: validate_tier_spec() compares the
#' realized N_obs AND the per-group row counts against it. A data edit that silently
#' changes a tier's composition then costs seconds instead of a refit. Deliberate
#' changes are a registry edit -- which is the point, because that edit shows up in
#' `git diff` where the old default never did.
TIER_SPECS <- list(
  `t1-grouped` = list(
    key = "t1-grouped",
    label = "Tier 1 rows, Darton placebo grouped",
    doc_ref = "CALIBRATION_WORKFLOW.md Step 1",
    tier_col = "tier1_active", individualize_darton = FALSE,
    drop_obs = character(), keep_obs = NULL,
    stage = "phi",
    expect = list(N_obs = 25L,
                  groups = c(ox_fev = 7L, md_fev = 11L, md_inf = 6L,
                             hornick_cond = 1L)),
    inert_pars = c("eta_lo", "kappa"),
    requires_params = character(),
    status_declared = "runnable",
    blocked_reason = NULL,
    note = paste("The diagnostic contrast that shows what individualizing Darton buys.",
                 "Keeps D-F-plac as one grouped binomial with the cohort GMT CoP.")
  ),
  `t1-indiv` = list(
    key = "t1-indiv",
    label = "Tier 1 rows, Darton placebo individualized",
    doc_ref = "CALIBRATION_WORKFLOW.md Step 1.5",
    tier_col = "tier1_active", individualize_darton = TRUE,
    drop_obs = character(), keep_obs = NULL,
    stage = "phi",
    expect = list(N_obs = 80L,
                  groups = c(ox_fev = 6L, md_fev = 11L, md_inf = 6L,
                             hornick_cond = 1L, ox_inf_indiv = 30L,
                             ox_fevginf_indiv = 26L)),
    inert_pars = c("eta_lo", "kappa"),
    requires_params = character(),
    status_declared = "runnable",
    blocked_reason = NULL,
    note = "THE FIT. Was labelled 'Tier 1' and written to results/tier1/."
  ),
  `t2-grouped` = list(
    key = "t2-grouped",
    label = "Tier 2 rows (+Oxford shedding), Darton placebo grouped",
    doc_ref = "CALIBRATION_WORKFLOW.md Step 3",
    tier_col = "tier2_active", individualize_darton = FALSE,
    drop_obs = character(), keep_obs = NULL,
    stage = "phi-eta",
    expect = list(N_obs = 31L,
                  groups = c(ox_fev = 7L, ox_inf = 6L, md_fev = 11L, md_inf = 6L,
                             hornick_cond = 1L)),
    inert_pars = character(),
    requires_params = c("eta_lo", "kappa"),
    status_declared = "blocked",
    blocked_reason = paste(
      "SCIENTIFIC DECISION, not a code defect. (a) eta Option A (parametric",
      "eta_lo/kappa, what the .stan implements) vs Option C (fixed eta_fixed_optC,",
      "a CSV column no code reads) is undecided; eta_detection() is monotone",
      "DECREASING in dose while eta_fixed_optC is non-monotone (1.00@1e3, 0.62@1e4,",
      "0.94@1.82e4, 0.92@2e4). Because eta multiplies P_inf and shares N50_inf in its",
      "exponent, a misfit moves the BIOLOGICAL parameters instead of failing visibly.",
      "(b) ../joint_inference_plan.md Sec 2.6 EXCLUDES Oxford shedding on",
      "treatment-truncation grounds while the Tier 2 design restores it with eta --",
      "unresolved tension. (c) psi (Sec 2.8, adopted 34aac76) is unimplemented, and",
      "psi_stool is confounded with eta at the single Darton dose.",
      "Unblock deliberately with allow_blocked = TRUE.")
  ),
  `t2-indiv` = list(
    key = "t2-indiv",
    label = "Tier 2 rows (+Oxford shedding), Darton placebo individualized",
    doc_ref = "unnamed in the docs before TIER_LOCK.md",
    tier_col = "tier2_active", individualize_darton = TRUE,
    drop_obs = character(), keep_obs = NULL,
    stage = "phi-eta",
    expect = list(N_obs = 85L,
                  groups = c(ox_fev = 6L, ox_inf = 5L, md_fev = 11L, md_inf = 6L,
                             hornick_cond = 1L, ox_inf_indiv = 30L,
                             ox_fevginf_indiv = 26L)),
    inert_pars = character(),
    requires_params = c("eta_lo", "kappa"),
    status_declared = "blocked",
    blocked_reason = paste(
      "Everything blocking t2-grouped, plus: Darton contributes no group-2 row here",
      "(the grouped D-I-plac is dropped as a double count of the 30 ox_inf_indiv rows",
      "for the same volunteers), so eta is identified by 5 rows -- W-I-3/4 and the",
      "three Jin arms. N_obs 85 = 86 - 1 for that drop.")
  )
)

# ---- accessors ----------------------------------------------------------------

#' Derived parameter lists. `MM_RAW_PARS` (model_math.R) stays the ONE list of the
#' model's parameters; everything the harness reports is computed from it minus what
#' a tier declares inert. This replaced 12 hardcoded name vectors that encoded four
#' distinct lists -- so eta_lo/kappa now migrate from inert to reported by themselves
#' when a t2-* tier is fit, which was the whole point.
DERIVED_REPORT_PARS <- "log10_N50_fevginf"                    # transformed, log scale
DERIVED_LINEAR_PARS <- c("N50_inf", "N50_fevginf", "delta")   # redundant log<->linear

tier_spec <- function(key) {
  if (!is.character(key) || length(key) != 1L || !key %in% names(TIER_SPECS))
    stop(sprintf("unknown tier '%s'. Valid: %s",
                 paste(key, collapse = ", "), paste(names(TIER_SPECS), collapse = ", ")),
         call. = FALSE)
  TIER_SPECS[[key]]
}

#' @param status optional filter, e.g. "runnable" (uses the DECLARED status; call
#'   tier_status() for the effective one, which needs the model).
tier_keys <- function(status = NULL) {
  ks <- names(TIER_SPECS)
  if (is.null(status)) return(ks)
  ks[vapply(ks, function(k) TIER_SPECS[[k]]$status_declared %in% status, logical(1))]
}

tier_report_pars <- function(spec)
  setdiff(c(MM_RAW_PARS, DERIVED_REPORT_PARS), spec$inert_pars)

tier_plot_pars <- function(spec)
  c(tier_report_pars(spec), DERIVED_LINEAR_PARS)

tier_key_pars <- function(spec)
  setdiff(tier_report_pars(spec), c(DERIVED_REPORT_PARS, "phi0_a", "phi0_b"))

#' Run directory: <data config>__<model stage>, the two things that determine what a
#' fit IS. A stage is never reused, so a new stage is a new directory and the
#' pre-increment fit is not clobbered.
tier_out_dir <- function(spec, root = "results", prior = FALSE, stage = spec$stage)
  file.path(root, sprintf("%s__%s%s", spec$key, stage, if (prior) "-prior" else ""))

#' Effective status: declared, auto-downgraded when the .stan cannot support the tier.
#' A blocked rung stays DECLARED -- "not started" was never true of Tier 2, and
#' deleting it from the registry would lose the reason.
tier_status <- function(spec, mod = NULL) {
  if (is.null(mod)) return(spec$status_declared)
  missing <- setdiff(spec$requires_params, names(mod$variables()$parameters))
  if (length(missing))
    return("blocked")
  spec$status_declared
}

# ---- validation ---------------------------------------------------------------

.group_table <- function(obs) {
  tb <- table(obs$likelihood_group)
  stats::setNames(as.integer(tb), names(tb))
}

#' Assert the realized data match the declaration, and the stage token match the .stan.
#' @param mod optional cmdstan_model(compile = FALSE) -- stanc only, no compilation.
validate_tier_spec <- function(spec, stan_data, mod = NULL, check_curve_specs = TRUE) {
  obs <- attr(stan_data, "obs")
  if (is.null(obs)) stop("stan_data has no \"obs\" attribute", call. = FALSE)

  if (!identical(as.integer(stan_data$N_obs), as.integer(spec$expect$N_obs)))
    stop(sprintf("tier '%s': N_obs is %d, registry declares %d. If this change is
  intended, edit TIER_SPECS and regenerate TIER_LADDER.md.",
                 spec$key, stan_data$N_obs, spec$expect$N_obs), call. = FALSE)

  got <- .group_table(obs); want <- spec$expect$groups
  if (!identical(got[order(names(got))], want[order(names(want))])) {
    fmt <- function(x) paste(sprintf("%s=%d", names(x), x), collapse = ", ")
    stop(sprintf("tier '%s': group composition changed.\n  realized: %s\n  declared: %s",
                 spec$key, fmt(got), fmt(want)), call. = FALSE)
  }

  if (!is.null(mod)) {
    model_pars <- names(mod$variables()$parameters)
    unknown_inert <- setdiff(spec$inert_pars, model_pars)
    if (length(unknown_inert))
      stop(sprintf("tier '%s': inert_pars not in the .stan: %s", spec$key,
                   paste(unknown_inert, collapse = ", ")), call. = FALSE)
    got_stage <- derive_stage_token(model_pars, any(stan_data$group == 2L))
    if (!identical(got_stage, spec$stage))
      stop(sprintf("tier '%s': stage token is '%s' but the registry declares '%s'.
  The .stan's parameter set changed. Update TIER_SPECS$%s$stage deliberately -- that
  rename is what keeps a run dir from misdescribing its own model.",
                   spec$key, got_stage, spec$stage, spec$key), call. = FALSE)
  }

  if (check_curve_specs && exists("validate_curve_specs") && exists("curve_specs"))
    validate_curve_specs(curve_specs(), obs)

  invisible(TRUE)
}

# ---- the one data entry point -------------------------------------------------

#' Build Stan data for a tier, validated, with the tier spec attached.
#'
#' Every entry point goes through here instead of calling build_stan_data() with a
#' literal tier_col. `attr(, "tier")` travels with the data so figures, the
#' provenance manifest and stan_data.rds can recover the tier without being told.
#'
#' Two-phase on purpose: the UNFILTERED tier is asserted against `expect` first, so a
#' scenario's drop_obs cannot mask a data-composition change; then the filter is
#' applied and the curve specs re-checked on what remains.
build_tier_data <- function(key,
                            data_csv = "dose_response_data.csv",
                            priors = load_priors(),
                            prior_only = 0L,
                            drop_obs = character(), keep_obs = NULL,
                            prior_overrides = list(),
                            allow_blocked = FALSE,
                            validate = TRUE,
                            mod = NULL) {
  if (!exists("build_stan_data")) source("data_prep.R")
  spec <- tier_spec(key)

  status <- tier_status(spec, mod)
  if (identical(status, "blocked") && !isTRUE(allow_blocked))
    stop(sprintf("tier '%s' is BLOCKED -- refusing to build.\n%s\nPass allow_blocked = TRUE to override.",
                 spec$key, spec$blocked_reason %||% "(no reason recorded)"), call. = FALSE)

  if (length(prior_overrides)) priors <- apply_prior_overrides(priors, prior_overrides)

  base <- build_stan_data(data_csv, priors,
                          tier_col = spec$tier_col,
                          prior_only = prior_only,
                          drop_obs = spec$drop_obs, keep_obs = spec$keep_obs,
                          individualize_darton = spec$individualize_darton)
  if (validate) validate_tier_spec(spec, base, mod, check_curve_specs = TRUE)

  drop_all <- union(spec$drop_obs, drop_obs)
  if (length(drop_all) || !is.null(keep_obs)) {
    sd <- build_stan_data(data_csv, priors,
                          tier_col = spec$tier_col, prior_only = prior_only,
                          drop_obs = drop_all, keep_obs = keep_obs,
                          individualize_darton = spec$individualize_darton)
    if (validate && exists("validate_curve_specs") && exists("curve_specs"))
      validate_curve_specs(curve_specs(), attr(sd, "obs"))
  } else {
    sd <- base
  }

  attr(sd, "tier") <- c(spec, list(status = status,
                                   n_obs_tier = as.integer(base$N_obs),
                                   dropped = setdiff(attr(base, "obs")$obs_id,
                                                     attr(sd, "obs")$obs_id)))
  sd
}

#' Registry-wide pre-flight. Zero Stan sampling; stanc only when stanc = TRUE.
#' This is the check that catches a CSV edit, and it costs seconds.
validate_all_tiers <- function(data_csv = "dose_response_data.csv", stanc = TRUE,
                               stan_file = "typhoid_dose_response.stan",
                               quiet = FALSE) {
  if (!exists("load_priors")) source("priors.R")
  mod <- NULL
  if (isTRUE(stanc)) {
    mod <- tryCatch(cmdstanr::cmdstan_model(stan_file, compile = FALSE),
                    error = function(e) { message("  [skip] stanc: ", conditionMessage(e)); NULL })
  }
  priors <- load_priors()
  rows <- lapply(tier_keys(), function(k) {
    spec <- tier_spec(k)
    sd <- build_tier_data(k, data_csv, priors, allow_blocked = TRUE, validate = TRUE,
                          mod = mod)
    obs <- attr(sd, "obs")
    data.frame(key = k, tier_col = spec$tier_col,
               individualize_darton = spec$individualize_darton,
               N_obs = sd$N_obs,
               grouped = sum(obs$n > 1), individual = sum(obs$n == 1),
               groups = paste(sprintf("%s=%d", names(.group_table(obs)),
                                      .group_table(obs)), collapse = " "),
               stage = spec$stage,
               status = tier_status(spec, mod),
               run_dir = tier_out_dir(spec),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (!quiet) {
    cat("\nTier registry pre-flight — all declarations verified against the data\n\n")
    cat(paste(knitr_table(out[, c("key", "N_obs", "grouped", "individual", "stage",
                                  "status")]), collapse = "\n"), "\n\n")
  }
  invisible(out)
}

# ---- generated documentation --------------------------------------------------

#' Write TIER_LADDER.md from the registry, with every count COMPUTED.
#'
#' Tracked in git precisely because results/ is not: a change to the data or to
#' build_stan_data()'s defaults then appears as a diff in this file, which is the
#' alarm that was missing while "Tier 1" silently meant 80 observations. No other
#' document in the repo states a tier observation count.
write_tier_ladder_md <- function(path = "TIER_LADDER.md",
                                 data_csv = "dose_response_data.csv") {
  tab <- validate_all_tiers(data_csv, stanc = TRUE, quiet = TRUE)
  sha <- tryCatch(sub("\\s.*$", "", system2("git", c("rev-parse", "--short", "HEAD"),
                                            stdout = TRUE, stderr = FALSE)),
                  error = function(e) "unknown")
  L <- c(
    "<!-- GENERATED FILE — do not hand-edit. Regenerate:",
    "     Rscript -e 'source(\"tier_specs.R\"); write_tier_ladder_md()' -->",
    "# Tier ladder — generated from tier_specs.R",
    "",
    sprintf("Generated %s from `tier_specs.R` + `%s` +", format(Sys.Date()), data_csv),
    sprintf("`../analysis_data/darton_individual_endpoints.csv` at git `%s`.", sha),
    "",
    "Counts are **computed** by calling `build_stan_data()` for each spec and counting",
    "the rows it returns; the generator hard-errors when a computed count disagrees",
    "with the spec's declared `expect`. If a number here is wrong, the registry or the",
    "data is wrong — fix it there and regenerate. **No other document in this repo",
    "states a tier observation count.**",
    "",
    paste(knitr_table(tab[, c("key", "tier_col", "individualize_darton", "N_obs",
                              "grouped", "individual", "stage", "status")]),
          collapse = "\n"),
    "",
    "## Likelihood groups per tier", "",
    paste(knitr_table(tab[, c("key", "groups")]), collapse = "\n"),
    "",
    "## Run directories", "",
    paste(knitr_table(tab[, c("key", "run_dir")]), collapse = "\n"),
    "",
    "Run dirs are `<data config>__<model stage>`. The stage token is derived from the",
    "`.stan`'s `parameters{}` and asserted against the registry, so a directory name",
    "cannot drift from the model that produced it. Prior-predictive companions take a",
    "`-prior` suffix. See `TIER_LOCK.md` for the naming rule and the blocked rungs.",
    "")
  writeLines(L, path)
  message("wrote ", path)
  invisible(tab)
}

if (sys.nframe() == 0) {
  setwd(calib_dir())
  source("priors.R"); source("data_prep.R")
  suppressPackageStartupMessages(library(cmdstanr))
  validate_all_tiers()
}
