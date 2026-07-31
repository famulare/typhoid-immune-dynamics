#' Per-run provenance: what code, what data, what settings produced a run directory.
#'
#' Why: results/ is gitignored and results.json recorded only what came OUT (model
#' name, wall time, sampler health, parameter table). Nothing tied a run dir to the
#' code that made it -- no tier, no input hashes, no git SHA, no seeds, no tool
#' versions. The practical cost was that 9 of 11 existing run dirs turned out to be
#' un-regenerable, with nothing in them saying so.
#'
#' A SEPARATE run_manifest.json rather than more fields in results.json, because:
#'   1. a manifest is written BEFORE sampling and updated after, so a crashed or
#'      killed run is still self-describing -- results.json can only exist post-fit;
#'   2. results.json is consumed by summarize_scenarios(), and should stay "what came
#'      out" rather than also carrying "what went in";
#'   3. ABSENCE of the manifest is then an unambiguous legacy marker. Inside
#'      results.json, "old file without the new keys" and "new file where a field
#'      failed to resolve" would look identical.
#' results.json gains exactly two scalars (tier, n_obs) so the scenario table can
#' show them without learning this schema.

if (!exists("calib_dir")) source("utils.R")

MANIFEST_FILE <- "run_manifest.json"
MANIFEST_SCHEMA_VERSION <- 1L

# ---- primitives ---------------------------------------------------------------

#' Git state of the working tree that produced a run. Every call is tryCatch'd to NA:
#' provenance recording must never be the thing that kills a fit.
.git_state <- function(dir = calib_dir()) {
  g <- function(...) tryCatch(
    sub("\\s+$", "", paste(system2("git", c("-C", dir, ...), stdout = TRUE,
                                   stderr = FALSE), collapse = "\n")),
    error = function(e) NA_character_, warning = function(w) NA_character_)
  sha <- g("rev-parse", "HEAD")
  dirty <- g("status", "--porcelain", "--", dir)
  list(git_sha = sha,
       git_short = if (is.na(sha)) NA_character_ else substr(sha, 1, 7),
       git_branch = g("rev-parse", "--abbrev-ref", "HEAD"),
       git_dirty = !is.na(dirty) && nzchar(dirty),
       git_dirty_files = if (!is.na(dirty) && nzchar(dirty))
         strsplit(dirty, "\n")[[1]] else character(0))
}

#' md5 + size + mtime for each input file. tools::md5sum() is base R -- no `digest`
#' dependency for something this peripheral.
file_fingerprints <- function(paths) {
  paths <- paths[!is.na(paths)]
  lapply(paths, function(p) {
    if (!file.exists(p)) return(list(path = p, exists = FALSE))
    fi <- file.info(p)
    list(path = p, exists = TRUE,
         md5 = unname(tools::md5sum(p)),
         bytes = as.numeric(fi$size),
         mtime_utc = format(as.POSIXct(fi$mtime, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"))
  })
}

#' The input files whose content determines a fit.
manifest_input_paths <- function(data_csv = "dose_response_data.csv",
                                 stan_file = "typhoid_dose_response.stan") {
  c(stan_file, data_csv, "priors.yaml", "reference_points.csv",
    file.path(dirname(data_csv), "..", "analysis_data",
              "darton_individual_endpoints.csv"))
}

.env_state <- function() {
  pkg <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                              error = function(e) NA_character_)
  list(r_version = R.version.string,
       platform = R.version$platform,
       cmdstan_version = tryCatch(cmdstanr::cmdstan_version(),
                                  error = function(e) NA_character_),
       packages = list(cmdstanr = pkg("cmdstanr"), posterior = pkg("posterior"),
                       bayesplot = pkg("bayesplot"), loo = pkg("loo"),
                       priorsense = pkg("priorsense"), ggplot2 = pkg("ggplot2"),
                       dplyr = pkg("dplyr"), matrixStats = pkg("matrixStats"),
                       jsonlite = pkg("jsonlite"), yaml = pkg("yaml"),
                       ragg = pkg("ragg")))
}

# ---- build / write ------------------------------------------------------------

#' Assemble a manifest. Pure: no I/O beyond hashing the declared input files.
#'
#' @param kind "posterior" | "prior_predictive" | "recovery_point" |
#'   "recovery_prior_draw" | "scenario" | "declared_only" | "unknown".
#' @param tier resolved tier spec, normally attr(stan_data, "tier").
#' @param sampler named list of the settings actually used (chains, iter_*, seed, ...).
run_manifest <- function(kind, tier = NULL, stan_data = NULL, priors = NULL,
                         model_name = NULL,
                         data_csv = "dose_response_data.csv",
                         stan_file = "typhoid_dose_response.stan",
                         sampler = list(), extra = list(), out_dir = NULL) {
  if (is.null(tier) && !is.null(stan_data)) tier <- attr(stan_data, "tier")
  obs <- if (!is.null(stan_data)) attr(stan_data, "obs") else NULL

  dat <- if (!is.null(obs)) {
    gt <- table(obs$likelihood_group)
    list(N_obs = nrow(obs),
         # obs_id order IS the p_pred / log_lik column order. Recorded so a stored
         # fit can be re-associated with its rows even if the CSV later changes.
         obs_id = obs$obs_id,
         groups = as.list(stats::setNames(as.integer(gt), names(gt))),
         n_total = sum(obs$n), y_total = sum(obs$y),
         n_obs_tier = tier$n_obs_tier %||% NA_integer_,
         dropped = tier$dropped %||% character(0),
         ladder = if (!is.null(stan_data$N_ladder))
           stan_data[c("N_ladder", "ladder_T", "ladder_count", "ladder_N")] else NULL)
  } else NULL

  list(
    schema_version = MANIFEST_SCHEMA_VERSION,
    run = list(kind = kind, model_name = model_name, out_dir = out_dir,
               started_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")),
    tier = if (!is.null(tier))
      tier[intersect(c("key", "label", "doc_ref", "tier_col", "individualize_darton",
                       "drop_obs", "keep_obs", "stage", "expect", "inert_pars",
                       "requires_params", "status", "status_declared",
                       "blocked_reason"), names(tier))] else NULL,
    data = dat,
    inputs = list(files = file_fingerprints(manifest_input_paths(data_csv, stan_file))),
    # The resolved pr_* scalars ARE the prior-configuration identity: an override is
    # then diffable from an edited yaml without hashing anything, and `normal` vs
    # `lognormal` (which collapse to the same pr_*_mu/_sd names) stay distinguishable
    # because the family list is carried alongside.
    priors = if (!is.null(priors)) list(
      families = lapply(priors, function(s) s$family),
      stan_values = tryCatch(priors_to_stan_data(priors), error = function(e) NULL),
      yaml_md5 = unname(tools::md5sum(calib_path("priors.yaml")))) else NULL,
    sampler = sampler,
    code = .git_state(),
    env = .env_state(),
    extra = extra
  )
}

#' Write (or update) run_manifest.json, folding in what is now known post-fit.
manifest_write <- function(manifest, out_dir, fit = NULL, health = NULL, tab = NULL,
                           outcome = list(), pars_dropped = character()) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  manifest$run$out_dir <- out_dir
  manifest$run$finished_utc <- format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")

  if (!is.null(fit)) {
    md <- tryCatch(fit$metadata(), error = function(e) NULL)
    if (!is.null(md))
      for (k in c("iter_sampling", "iter_warmup", "adapt_delta", "seed", "thin"))
        if (is.null(manifest$sampler[[k]]) && !is.null(md[[k]]))
          manifest$sampler[[k]] <- md[[k]]
  }
  if (!is.null(health))
    manifest$fit <- list(n_divergent = health$n_div, n_draws = health$n_draws,
                         div_rate = health$div_rate,
                         n_max_treedepth = health$n_max_td,
                         ebfmi_min = suppressWarnings(min(health$ebfmi)))
  if (!is.null(tab) && "rhat" %in% names(tab)) {
    manifest$fit$max_rhat <- suppressWarnings(max(tab$rhat, na.rm = TRUE))
    manifest$fit$min_ess_bulk <- suppressWarnings(min(tab$ess_bulk, na.rm = TRUE))
    manifest$fit$pars_diagnosed <- tab$variable
  }
  if (length(pars_dropped)) manifest$outputs$pars_dropped <- pars_dropped

  pngs <- basename(Sys.glob(file.path(out_dir, "*.png")))
  sub <- basename(list.dirs(out_dir, recursive = FALSE))
  manifest$outputs <- utils::modifyList(manifest$outputs %||% list(), list(
    figures = pngs,
    subdir_figures = as.list(stats::setNames(
      lapply(sub, function(s) basename(Sys.glob(file.path(out_dir, s, "*.png")))), sub)),
    has_fit_rds = file.exists(file.path(out_dir, "fit.rds")),
    stan_data_rds = file.exists(file.path(out_dir, "stan_data.rds")),
    curve_parity_max_abs_diff = if (file.exists(file.path(out_dir, "curve_parity.json")))
      tryCatch(jsonlite::read_json(file.path(out_dir, "curve_parity.json"))$max_abs_diff_vs_p_pred,
               error = function(e) NULL) else NULL))
  manifest$outputs <- utils::modifyList(manifest$outputs, outcome)
  manifest$regenerable <- isTRUE(manifest$outputs$has_fit_rds) &&
                          isTRUE(manifest$outputs$stan_data_rds)

  jsonlite::write_json(manifest, file.path(out_dir, MANIFEST_FILE),
                       auto_unbox = TRUE, pretty = TRUE, digits = 10, na = "null")
  invisible(manifest)
}

#' The ONE writer of stan_data.rds. Shared by diagnose_fit() and make_model_figures()
#' so a run without figures is still regenerable -- previously only the figure path
#' wrote it, which is why recover_once() and run_scenario(figures = FALSE) produced
#' dead directories.
save_stan_data <- function(stan_data, out_dir) {
  if (is.null(stan_data)) return(invisible(NULL))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(stan_data, file.path(out_dir, "stan_data.rds"))
  invisible(file.path(out_dir, "stan_data.rds"))
}

# ---- reading back -------------------------------------------------------------

read_manifest <- function(run_dir) {
  p <- file.path(run_dir, MANIFEST_FILE)
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
}

#' Recover the stan_data for a saved run: prefer stan_data.rds, else rebuild from the
#' manifest's tier. Replaces the old figures_from_dir() fallback, which rebuilt
#' tier1_active unconditionally and then hard-errored on the row-count mismatch --
#' i.e. it guessed, and the guess was wrong for every scenario and recovery dir.
resolve_run_stan_data <- function(run_dir, data_csv = "dose_response_data.csv") {
  p <- file.path(run_dir, "stan_data.rds")
  if (file.exists(p)) return(readRDS(p))
  m <- read_manifest(run_dir)
  if (is.null(m) || is.null(m$tier$key))
    stop("cannot resolve stan_data for ", run_dir,
         ": no stan_data.rds and no tier in run_manifest.json. This run predates
  provenance recording; refit it.", call. = FALSE)
  if (isTRUE(m$extra$y_synthetic))
    stop("refusing to rebuild stan_data for ", run_dir,
         ": this run was fit to SYNTHETIC y, which is not reconstructable from the
  data files. Its stan_data.rds is required and missing.", call. = FALSE)
  if (!exists("build_tier_data")) source("tier_specs.R")
  build_tier_data(m$tier$key, data_csv,
                  drop_obs = m$tier$drop_obs %||% character(),
                  keep_obs = m$tier$keep_obs,
                  allow_blocked = TRUE)
}

#' What is in results/, whether it can be regenerated, and what code state made it.
#' Cheap by construction: reads json and file.exists, never loads a fit.
#'
#' @param stamp_legacy write an HONEST stub manifest into pre-provenance dirs, so the
#'   invariant "every run dir has a manifest" holds without inventing a SHA or a tier
#'   the audit cannot know.
audit_run_dirs <- function(root = "results", stamp_legacy = FALSE, out_md = NULL,
                           data_csv = "dose_response_data.csv",
                           stan_file = "typhoid_dose_response.stan",
                           quiet = FALSE) {
  # Depth-independent: Sys.glob patterns are fixed-depth, so a hand-written set of
  # them silently misses anything deeper -- e.g. results/scenarios/<tier>/<label>/,
  # which is three levels down and was invisible to the earlier two-level version.
  cand <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  dirs <- sort(cand[vapply(cand, function(d)
    file.exists(file.path(d, "fit.rds")) || file.exists(file.path(d, "results.json")),
    logical(1))])
  if (!length(dirs)) {
    if (!quiet) message("no run directories under ", root)
    return(invisible(NULL))
  }

  cur <- stats::setNames(vapply(manifest_input_paths(data_csv, stan_file),
                                function(p) if (file.exists(p)) unname(tools::md5sum(p))
                                            else NA_character_, character(1)),
                         basename(manifest_input_paths(data_csv, stan_file)))
  rows <- lapply(dirs, function(d) {
    m <- read_manifest(d)
    has_fit <- file.exists(file.path(d, "fit.rds"))
    has_sd  <- file.exists(file.path(d, "stan_data.rds"))
    md5_of <- function(base) {
      if (is.null(m$inputs$files)) return(NA_character_)
      f <- m$inputs$files
      hit <- which(basename(unlist(f$path)) == base)
      if (!length(hit)) NA_character_ else unlist(f$md5)[hit[1]]
    }
    same <- function(base) {
      h <- md5_of(base)
      if (is.na(h) || is.na(cur[[base]])) "unknown" else if (h == cur[[base]]) "Y" else "N"
    }
    if (stamp_legacy && is.null(m)) {
      jsonlite::write_json(list(
        schema_version = MANIFEST_SCHEMA_VERSION, legacy = TRUE,
        run = list(kind = "unknown", out_dir = d),
        code_state = "unknown", regenerable = has_fit && has_sd,
        note = paste("Stamped retroactively by audit_run_dirs(); no provenance was",
                     "recorded when this run was produced. Tier, git SHA, seeds and",
                     "input hashes are NOT recoverable -- refit to obtain them.")),
        file.path(d, MANIFEST_FILE), auto_unbox = TRUE, pretty = TRUE)
    }
    data.frame(
      dir = d,
      kind = m$run$kind %||% (if (is.null(m)) "legacy" else "unknown"),
      tier = m$tier$key %||% NA_character_,
      stage = m$tier$stage %||% NA_character_,
      n_obs = m$data$N_obs %||% NA_integer_,
      fit = has_fit, stan_data = has_sd,
      manifest = !is.null(m),
      regenerable = has_fit && has_sd,
      git = m$code$git_short %||% NA_character_,
      dirty = m$code$git_dirty %||% NA,
      stan_cur = same("typhoid_dose_response.stan"),
      data_cur = same("dose_response_data.csv"),
      n_png = length(Sys.glob(file.path(d, "*.png"))),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (!quiet) cat("\nRun directory audit —", nrow(out), "dirs,", sum(out$regenerable),
      "regenerable,", sum(!out$manifest), "without a manifest\n\n")
  if (!quiet)
    cat(paste(knitr_table(out[, c("dir", "tier", "stage", "n_obs", "regenerable",
                                  "manifest", "git", "stan_cur")]), collapse = "\n"), "\n")
  if (!is.null(out_md)) {
    dir.create(dirname(out_md), showWarnings = FALSE, recursive = TRUE)
    writeLines(c("# Run directory audit", "",
                 sprintf("Generated %s. `stan_cur`/`data_cur`: does the run's recorded",
                         format(Sys.time(), "%Y-%m-%d %H:%M")),
                 "input hash still match the file on disk?", "",
                 knitr_table(out), ""), out_md)
    message("wrote ", out_md)
  }
  invisible(out)
}

if (sys.nframe() == 0) {
  setwd(calib_dir())
  suppressPackageStartupMessages(library(jsonlite))
  invisible(audit_run_dirs())
}
