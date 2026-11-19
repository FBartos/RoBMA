#!/usr/bin/env Rscript
# Verification of the coarse-corrected selection sampler as the default
# (RoBMA.options() 'selection.sampler' and 'selection.coarse_min_acceptance').
# Not part of any test profile. Every run records its build (package versions,
# DLL md5 sums, first library path).
#
#   Rscript tools/verify-coarse-default.R --part=wider --out=<dir>
#       [--workers=6] [--seeds=10] [--models=<regex>]
#   Rscript tools/verify-coarse-default.R --part=rerun --out=<dir>
#       [--workers=6] [--seeds=10]
#   Rscript tools/verify-coarse-default.R --part=cost --out=<dir> --other-library=<dir>
#       [--seeds=3]
#   Rscript tools/verify-coarse-default.R --part=identity --out=<dir> --other-library=<dir>
#   Rscript tools/verify-coarse-default.R --part=selftest
#
# wider: arms 'coarse' (selection.sampler = "coarse_corrected", the default
# threshold) and 'default' (selection.sampler = "default"), --seeds runs per arm
# of 3 chains x 4,000 draws (package defaults otherwise, serial chains in every
# fit), arms interleaved by seed, the runs spread over --workers worker
# processes of this build. Per run: the mean of every monitored quantity
# (categorical indicators as the frequency of each level), the corrected
# sampler's record (selection_sampler_info(fit = )), the warnings, the JAGS
# seconds and the effective draws per second of the global quantities and,
# separately, of the model-inclusion indicators. Decision per model, on the
# Welch statistics of the run means of every non-degenerate quantity: pass when
# no |t| > 3.5 and the count of |t| > 2 is at most twice its expectation under
# the Welch t distribution; a failing model is rerun once ('--part=rerun') with
# fresh seeds, judged on its failing quantities only. The count treats the
# quantities as independent, though some are one (complementary weights and
# their prior parameters, strongly correlated run means), so a single deviation
# can count several times. Every planned run is
# counted (planned, completed, failed, missing); a failed or missing run, an arm
# with fewer than two completed runs or no comparable quantity makes the model
# 'incomplete', never a pass. Models (dat.assink2016 unless named, V from
# metafor::vcalc()): bselmodel.mv and RoBMA.mv with the full V; RoBMA.mv on the
# sampling variances only; 'dense', blocks of three Markov-correlated rows
# (twelve studies) fitted in both arms with relative_tolerance = .1, 16
# scrambles and up to 32,768 points per scramble (under the default selection
# control the dense normalizer is rejected by its QMC diagnostics in both arms);
# 'rank_one' (constant correlation .6); 'factor' (the three studies whose full-V
# blocks route to the factor rules); 'best_ensemble' (RoBMA.mv whose bias
# branches all use the best-weight rule, six studies); 'meta_regression'
# (bselmodel.mv, full V, ~ deltype); dat.knapp2017 and dat.berkey1998
# (bselmodel.mv). A run stopped by the selection-normalizer diagnostic is
# replaced by the next unused seed of its arm (after the planned seeds), at most
# three per model and arm; replaced runs are listed and written to
# '<part>-replaced.csv'. A quantity whose mean is missing from a completed run
# (other than an indicator level the run never visited, frequency zero) has no
# statistic and makes the model 'incomplete'; such cases are listed. The
# efficiency table compares the arms, with the attempts BayesTools stopped and
# restarted, the runs the diagnostic stopped and the wall hours; for the
# model-averaged models the 'default' arm declares the allocation block and the
# 'coarse' arm declares it too unless the corrected sampler takes the
# allocation nodes (the block methods and samplers of each run show which).
#
# cost: Assink et al. bselmodel.mv (full V) with sampler.blocks = FALSE and the
# corrected sampler; arm 'this' is this build with the fallback off
# (selection.coarse_min_acceptance = 0), arm 'other' the build of
# --other-library with the options it has; --seeds pairs, one R process per arm
# and seed, run in the order this, other, this, other. The draws of a pair must
# be identical (any difference is reported) and the seconds per 1,000
# iterations of the arms are compared.
#
# identity: with selection.sampler = "default" in this build and the defaults of
# the build of --other-library, one R process each: the sampler assignment and
# the draws of a univariate selection model, a .mv selection model, a .mv
# ensemble and a model without selection must be equal; and, in this build, a
# univariate selection model, a univariate ensemble, a model without selection
# and a .mv ensemble on sampling variances (no coarse-capable selection child;
# the ensemble keeps its allocation block) must be equal under both values of
# selection.sampler.
#
# --part=selftest runs the drivers on synthetic run files in a temporary
# directory (no fits) and stops on a wrong decision, shape or count, or on any
# error.

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}
part      <- opt("part", "wider")
out_dir   <- opt("out", if (part %in% c("selftest", "job")) tempfile("verify-coarse-") else NULL)
workers   <- as.integer(opt("workers", "6"))
seeds     <- as.integer(opt("seeds", if (identical(part, "cost")) "3" else "10"))
models_re <- opt("models", ".")
other_lib <- opt("other-library")
if (is.null(out_dir)) stop("Supply --out=<directory>.", call. = FALSE)
if (!part %in% c("wider", "rerun", "cost", "identity", "selftest", "job")) {
  stop("--part must be wider, rerun, cost, identity or selftest.", call. = FALSE)
}
if (part %in% c("cost", "identity") && is.null(other_lib)) {
  stop("--part=", part, " needs --other-library=<library of the build to compare with>.", call. = FALSE)
}
dir.create(file.path(out_dir, "runs"), recursive = TRUE, showWarnings = FALSE)
options(warn = 1, width = 200)

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]),
                             winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(RoBMA))

# The build a process runs: package versions, DLL md5 sums and the library
# paths it loads from.
build_info <- function() {

  list(
    RoBMA          = as.character(utils::packageVersion("RoBMA")),
    BayesTools     = as.character(utils::packageVersion("BayesTools")),
    RoBMA_dll      = unname(tools::md5sum(getLoadedDLLs()[["RoBMA"]][["path"]])),
    BayesTools_dll = unname(tools::md5sum(getLoadedDLLs()[["BayesTools"]][["path"]])),
    RoBMA_path     = getNamespaceInfo(asNamespace("RoBMA"), "path"),
    library        = .libPaths()[1L]
  )
}

build_text <- function(build) {

  paste0("RoBMA ", build[["RoBMA"]], " (dll ", build[["RoBMA_dll"]], "), BayesTools ", build[["BayesTools"]],
         " (dll ", build[["BayesTools_dll"]], "), library ", build[["library"]])
}

if (!identical(part, "job")) cat("build:", build_text(build_info()), "\n")


# Models -------------------------------------------------------------------------

# Markov correlations rho^|i - j| within each cluster: blocks of three or more
# rows that no factor structure reproduces, so they keep the dense route.
markov_covariance <- function(vi, cluster, rho = .8) {

  out <- diag(vi, length(vi))
  for (group in unique(cluster)) {
    rows <- which(cluster == group)
    if (length(rows) > 1L) {
      out[rows, rows] <- outer(seq_along(rows), seq_along(rows), function(i, j) rho^abs(i - j)) *
        tcrossprod(sqrt(vi[rows]))
    }
  }

  return(out)
}

model_specs <- function() {

  data("dat.assink2016", package = "metadat", envir = environment())
  assink <- dat.assink2016
  V_full <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid, rho = c(.7, .5), data = assink)
  V_one  <- metafor::vcalc(vi, cluster = study, obs = esid, rho = .6, data = assink)
  # the three studies whose full-V blocks route to the factor rules
  factor_studies <- assink[assink[["study"]] %in% c(11, 12, 16), ]
  V_factor <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid, rho = c(.7, .5),
                             data = factor_studies)
  # the first three estimates of every study with at least three: blocks of
  # three rows with Markov correlations (.8, .64, .8), which route to the dense
  # rules; the default selection control rejects their normalizer in both arms
  # (and relative_tolerance .05 with 16 scrambles still did, at a relative MCSE
  # of .023-.028), so both arms fit them with a looser tolerance and a larger
  # QMC budget
  position <- stats::ave(assink[["esid"]], assink[["study"]], FUN = seq_along)
  count    <- stats::ave(assink[["esid"]], assink[["study"]], FUN = length)
  dense    <- assink[position <= 3L & count >= 3L, ]
  dense_control <- RoBMA::set_selection_likelihood_control(
    relative_tolerance = .1, scrambles = 16L, max_points_per_scramble = 32768L)
  best  <- RoBMA::selection_model(weight_rule = "best", group = "study")
  best_studies <- assink[assink[["study"]] %in% c(1, 2, 4, 6, 7, 10), ]
  V_best <- metafor::vcalc(vi, cluster = study, obs = esid, rho = .6, data = best_studies)

  data("dat.knapp2017", package = "metadat", envir = environment())
  knapp <- dat.knapp2017
  # estimates of a study (several per task) correlated .5
  knapp[["id"]] <- seq_len(nrow(knapp))
  V_knapp <- metafor::vcalc(vi, cluster = study, obs = id, rho = .5, data = knapp)
  data("dat.berkey1998", package = "metadat", envir = environment())
  berkey <- dat.berkey1998
  V_berkey <- metafor::bldiag(lapply(split(berkey[, c("v1i", "v2i")], berkey[["trial"]]), as.matrix))

  mv <- function(constructor, V, data, ...) {
    list(constructor = constructor, args = c(list(yi = quote(yi), V = V, data = data, measure = "SMD"),
                                            list(...)))
  }
  list(
    bselmodel.mv         = mv("bselmodel.mv", V_full, assink, random = ~ 1 | study / esid),
    RoBMA.mv             = mv("RoBMA.mv", V_full, assink, random = ~ 1 | study / esid),
    RoBMA.mv_vi          = mv("RoBMA.mv", quote(vi), assink, random = ~ 1 | study / esid),
    dense                = mv("bselmodel.mv", markov_covariance(dense[["vi"]], dense[["study"]]), dense,
                              random = ~ 1 | study / esid, selection_control = dense_control),
    rank_one             = mv("bselmodel.mv", V_one, assink, random = ~ 1 | study / esid),
    factor               = mv("bselmodel.mv", V_factor, factor_studies, random = ~ 1 | study / esid),
    # six studies with rank-one blocks: a best-weight ensemble on the full data
    # takes about 0.1 s per iteration, on the full V far longer
    best_ensemble        = mv("RoBMA.mv", V_best, best_studies, random = ~ 1 | study / esid, prior_bias = list(
                             BayesTools::prior_weightfunction("one-sided", steps = .025, model = best),
                             BayesTools::prior_weightfunction("one-sided", steps = c(.025, .05), model = best))),
    meta_regression      = mv("bselmodel.mv", V_full, assink, random = ~ 1 | study / esid, mods = ~ deltype),
    knapp2017            = list(constructor = "bselmodel.mv", args = list(
                             yi = quote(yi), V = V_knapp, data = knapp, measure = "SMD",
                             random = ~ 1 | study / id)),
    berkey1998           = list(constructor = "bselmodel.mv", args = list(
                             yi = quote(yi), V = V_berkey, ni = quote(ni), data = berkey, measure = "GEN",
                             random = ~ 1 | trial / outcome))
  )
}

arm_options <- list(
  coarse  = list(selection.sampler = "coarse_corrected"),
  default = list(selection.sampler = "default")
)


# One run of the wider testing: fit with the arm's options and summarize the
# draws. Runs in a worker; returns a list and never throws.
run_job <- function(job, specs, arm_options, sample = 4000L) {

  spec <- specs[[job[["model"]]]]
  do.call(RoBMA::RoBMA.options, arm_options[[job[["arm"]]]])
  warnings <- character()
  classes  <- character()
  started  <- Sys.time()
  fit <- tryCatch(withCallingHandlers(
    do.call(getExportedValue("RoBMA", spec[["constructor"]]), c(spec[["args"]], list(
      chains = 3L, sample = sample, seed = job[["seed"]], parallel = FALSE, silent = TRUE
    ))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      classes  <<- c(classes, class(w)[[1L]])
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  wall <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  control <- spec[["args"]][["selection_control"]]
  result <- list(job = job, build = build_info(), wall = wall, warnings = unique(warnings),
                 warning_classes = unique(classes), jags_seconds = NA_real_,
                 control = if (is.null(control)) "default" else paste0(
                   "relative_tolerance ", control[["relative_tolerance"]], ", scrambles ", control[["scrambles"]],
                   ", points ", control[["points_per_scramble"]], "-", control[["max_points_per_scramble"]]))
  if (inherits(fit, "error") || is.null(fit[["fit"]][["mcmc"]])) {
    result[["error"]] <- if (inherits(fit, "error")) conditionMessage(fit) else
      paste(c("no posterior draws", fit[["fit"]][["errors"]]), collapse = ": ")
    return(result)
  }

  draws <- fit[["fit"]][["mcmc"]]
  samplers <- fit[["fit"]][["samplers"]]
  result[["jags_seconds"]] <- as.numeric(fit[["fit"]][["timetaken"]], units = "secs")
  result[["samplers"]]     <- if (!is.null(samplers)) table(samplers[["Sampler"]]) else NULL
  result[["record"]]       <- RoBMA::selection_sampler_info(fit = fit)
  result[["blocks"]]       <- block_methods(fit)
  result[["means"]]        <- quantity_means(draws)
  result[["ess"]]          <- coda::effectiveSize(draws)
  result[["indicator_ess"]] <- indicator_ess(draws)
  result[["roles"]]        <- coordinate_roles(fit)
  result[["fit_warnings"]] <- unique(attr(fit[["fit"]], "warnings", exact = TRUE))

  return(result)
}

# The selection blocks of a fit by method and size ("dense:3 x12, singleton:1 x3").
block_methods <- function(fit) {

  plan  <- utils::getFromNamespace(".data_selection_execution_plan", "RoBMA")(fit[["data"]])
  sizes <- vapply(plan[["row_blocks"]], length, integer(1))
  tally <- table(paste0(plan[["block_methods"]], ":", sizes))

  return(paste(names(tally), as.integer(tally), sep = " x", collapse = ", "))
}

# Run means of the monitored quantities; a categorical indicator gives the
# frequency of each of its levels.
quantity_means <- function(draws) {

  matrix <- as.matrix(draws)
  means  <- numeric()
  for (name in colnames(matrix)) {
    values <- matrix[, name]
    if (grepl("indicator", name, fixed = TRUE)) {
      for (level in sort(unique(values))) {
        means[paste0(name, "==", level)] <- mean(values == level)
      }
    } else {
      means[name] <- mean(values)
    }
  }

  return(means)
}

# Effective draws of each level of each model-inclusion indicator (the series
# of draws at that level), summed over chains; levels the chains never leave or
# never visit have none and are dropped.
indicator_ess <- function(draws) {

  names <- grep("indicator", coda::varnames(draws), fixed = TRUE, value = TRUE)
  ess   <- numeric()
  for (name in names) {
    values <- as.matrix(draws[, name])
    for (level in sort(unique(as.numeric(values)))) {
      series <- coda::as.mcmc.list(lapply(draws[, name], function(chain) {
        coda::mcmc(as.numeric(chain == level))
      }))
      if (stats::var(unlist(series)) > 0) {
        ess[paste0(name, "==", level)] <- sum(coda::effectiveSize(series))
      }
    }
  }

  return(ess)
}

# The convergence role of each coordinate, declared from the priors at fit time.
coordinate_roles <- function(fit) {

  coordinates <- BayesTools::parameter_coordinates(fit[["fit"]])

  return(stats::setNames(coordinates[["convergence_role"]], coordinates[["coordinate_name"]]))
}

# Global scalar parameters: by monitored name and their role, sampled
# coordinates and those derived from them only (not structural, auxiliary or
# indicator ones), not unit-level latents or random-effect coordinates.
global_parameters <- function(names, roles) {

  keep <- !grepl("indicator", names, fixed = TRUE) &
    !grepl("_xRE_(Zx|COEFx|GROUP_Zx|UNIT_COEFx)", names) &
    !grepl("^prior_par_", names) &
    !grepl("^(theta|gamma|sampling_z|pi|phi)\\[", names) &
    names %in% names(roles)[roles %in% c("sampled", "derived")]

  return(names[keep])
}

run_file <- function(job, dir = out_dir) {

  file.path(dir, "runs", sprintf("%s-%s-%s-%03d.rds", job[["part"]], job[["model"]], job[["arm"]], job[["seed"]]))
}

run_jobs <- function(jobs, specs, n_workers) {

  todo <- Filter(function(job) !file.exists(run_file(job)), jobs)
  cat(length(jobs), "jobs,", length(jobs) - length(todo), "done before,", length(todo), "to run\n")
  if (length(todo) == 0L) return(invisible(NULL))
  cluster <- parallel::makePSOCKcluster(max(1L, min(n_workers, length(todo))), rscript_args = "--vanilla")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(paths) {
    .libPaths(paths)
    Sys.setenv(OMP_NUM_THREADS = "1")
    suppressPackageStartupMessages(library(RoBMA))
    invisible(TRUE)
  }, .libPaths())
  parallel::clusterExport(cluster, c("run_job", "quantity_means", "indicator_ess", "coordinate_roles",
                                     "block_methods", "build_info", "arm_options", "specs"), envir = environment())
  builds <- parallel::clusterCall(cluster, build_info)
  if (!all(vapply(builds, identical, logical(1), build_info()))) {
    stop("A worker loaded a different RoBMA or BayesTools build than this session.", call. = FALSE)
  }
  worker <- function(i, todo, files) {
    job    <- todo[[i]]
    result <- run_job(job, specs, arm_options)
    saveRDS(result, files[[i]])
    sprintf("%s %s %s seed %d: %.0f s%s", job[["part"]], job[["model"]], job[["arm"]], job[["seed"]],
            result[["wall"]], if (!is.null(result[["error"]])) paste(" ERROR", result[["error"]]) else "")
  }
  environment(worker) <- globalenv()
  files <- vapply(todo, run_file, character(1))
  log <- parallel::clusterApplyLB(cluster, seq_along(todo), worker, todo = todo, files = files)
  cat(unlist(log), sep = "\n")
}

# The fields every run record has; a failed or missing run has NA timings.
complete_run <- function(run) {

  defaults <- list(wall = NA_real_, jags_seconds = NA_real_, warnings = character(),
                   warning_classes = character())
  for (field in names(defaults)) {
    if (is.null(run[[field]])) run[[field]] <- defaults[[field]]
  }

  return(run)
}

# The saved runs of the jobs; a job without a saved run is marked missing.
read_runs <- function(jobs, dir = out_dir) {

  lapply(jobs, function(job) {
    file <- run_file(job, dir)
    if (!file.exists(file)) return(complete_run(list(job = job, error = "missing run file", missing = TRUE)))
    complete_run(readRDS(file))
  })
}

# Planned, completed, failed and missing runs per model and arm.
run_counts <- function(runs) {

  if (length(runs) == 0L) {
    return(data.frame(model = character(), arm = character(), planned = integer(), completed = integer(),
                      failed = integer(), missing = integer(), stringsAsFactors = FALSE))
  }
  rows <- data.frame(
    model   = vapply(runs, function(run) run[["job"]][["model"]], character(1)),
    arm     = vapply(runs, function(run) run[["job"]][["arm"]], character(1)),
    missing = vapply(runs, function(run) isTRUE(run[["missing"]]), logical(1)),
    failed  = vapply(runs, function(run) !is.null(run[["error"]]) && !isTRUE(run[["missing"]]), logical(1)),
    stringsAsFactors = FALSE
  )
  rows[["completed"]] <- as.integer(!rows[["missing"]] & !rows[["failed"]])
  rows[["failed"]]    <- as.integer(rows[["failed"]])
  rows[["missing"]]   <- as.integer(rows[["missing"]])
  rows[["planned"]]   <- 1L
  counts <- stats::aggregate(cbind(planned, completed, failed, missing) ~ model + arm, data = rows, FUN = sum)

  return(counts[order(counts[["model"]], counts[["arm"]]), , drop = FALSE])
}

# A run-by-quantity matrix (one row per run, one column per name).
run_matrix <- function(runs, names, extract) {

  values <- matrix(NA_real_, nrow = length(runs), ncol = length(names), dimnames = list(NULL, names))
  for (i in seq_along(runs)) {
    if (length(names) > 0L) values[i, ] <- unname(extract(runs[[i]], names))
  }

  return(values)
}


# Welch comparison -----------------------------------------------------------------

welch_table <- function(runs, model, quantities = NULL) {

  runs <- Filter(function(run) run[["job"]][["model"]] == model && is.null(run[["error"]]), runs)
  arms <- vapply(runs, function(run) run[["job"]][["arm"]], character(1))
  # requested quantities (a rerun) are compared exactly, present in the runs or not
  names_all <- if (is.null(quantities)) unique(unlist(lapply(runs, function(run) names(run[["means"]])))) else
    unique(quantities)
  # quantity_means() names only the levels of an indicator that a run visited:
  # a level is zero in a run whose means name no such level but hold another
  # level of the same indicator; a stored NA, an indicator the run lacks
  # altogether and any other quantity missing from a completed run stay missing
  means <- run_matrix(runs, names_all, function(run, names) {
    stored    <- run[["means"]]
    values    <- stored[names]
    indicator <- sub("==.*$", "", names)
    visited   <- sub("==.*$", "", names(stored)[grepl("==", names(stored), fixed = TRUE) & !is.na(stored)])
    unvisited <- grepl("indicator", names, fixed = TRUE) & grepl("==", names, fixed = TRUE) &
      !names %in% names(stored) & indicator %in% visited
    values[unvisited] <- 0
    values
  })
  missing <- which(is.na(means), arr.ind = TRUE)

  return(list(
    table   = welch_rows(means, arms, model),
    missing = data.frame(
      model = rep(model, nrow(missing)), quantity = colnames(means)[missing[, "col"]],
      arm = arms[missing[, "row"]],
      seed = vapply(runs[missing[, "row"]], function(run) as.integer(run[["job"]][["seed"]]), integer(1)),
      stringsAsFactors = FALSE)
  ))
}

# Welch statistics of the per-run values of each column (runs in rows). With
# fewer than two completed runs in an arm, or a value missing from a completed
# run, there is no statistic: the rows carry NA statistics (the decision is then
# 'incomplete').
welch_rows <- function(values, arms, model) {

  row <- function(name, a, b, se, t, df) {
    data.frame(model = model, quantity = name,
               mean_coarse = if (length(a) > 0L) mean(a) else NA_real_,
               mean_default = if (length(b) > 0L) mean(b) else NA_real_,
               se = se, t = t, df = df,
               p_beyond_2 = if (is.na(df)) NA_real_ else 2 * stats::pt(-2, df), stringsAsFactors = FALSE)
  }
  testable <- sum(arms == "coarse") >= 2L && sum(arms == "default") >= 2L
  rows <- lapply(colnames(values), function(name) {
    a <- values[arms == "coarse", name]
    b <- values[arms == "default", name]
    if (!testable || anyNA(c(a, b))) return(row(name, a, b, NA_real_, NA_real_, NA_real_))
    # degenerate: the same value in every run of both arms
    if (stats::var(c(a, b)) == 0) return(NULL)
    se2 <- stats::var(a) / length(a) + stats::var(b) / length(b)
    if (se2 > 0) {
      df <- se2^2 / ((stats::var(a) / length(a))^2 / (length(a) - 1) +
                     (stats::var(b) / length(b))^2 / (length(b) - 1))
      t  <- (mean(a) - mean(b)) / sqrt(se2)
    } else {
      # constant within each arm but different between them
      df <- length(a) + length(b) - 2
      t  <- Inf
    }
    row(name, a, b, sqrt(se2), t, df)
  })

  do.call(rbind, rows)
}

# The decision of one model: 'incomplete' when any planned run failed or is
# missing, an arm has fewer than two completed runs, a quantity has no statistic
# (a mean missing from a completed run), or no quantity could be compared;
# otherwise 'pass' or 'fail' by the rule.
decide <- function(table, model, counts) {

  counts <- counts[counts[["model"]] == model, , drop = FALSE]
  arm_count <- function(arm, field) {
    value <- counts[[field]][counts[["arm"]] == arm]
    if (length(value) == 0L) 0L else as.integer(value)
  }
  t <- if (is.null(table)) numeric() else table[["t"]]
  p <- if (is.null(table)) numeric() else table[["p_beyond_2"]]
  complete <- sum(counts[["failed"]]) == 0L && sum(counts[["missing"]]) == 0L &&
    all(vapply(c("coarse", "default"), arm_count, integer(1), field = "completed") >= 2L) &&
    length(t) > 0L && !anyNA(t)
  rule <- complete && !any(abs(t) > 3.5) && sum(abs(t) > 2) <= 2 * sum(p)
  data.frame(
    model = model,
    planned_coarse = arm_count("coarse", "planned"), completed_coarse = arm_count("coarse", "completed"),
    failed_coarse = arm_count("coarse", "failed"), missing_coarse = arm_count("coarse", "missing"),
    planned_default = arm_count("default", "planned"), completed_default = arm_count("default", "completed"),
    failed_default = arm_count("default", "failed"), missing_default = arm_count("default", "missing"),
    quantities = length(t),
    max_abs_t = if (any(!is.na(t))) max(abs(t), na.rm = TRUE) else NA_real_,
    beyond_3.5 = sum(abs(t) > 3.5, na.rm = TRUE),
    beyond_2 = sum(abs(t) > 2, na.rm = TRUE),
    expected_beyond_2 = sum(p, na.rm = TRUE),
    decision = if (!complete) "incomplete" else if (rule) "pass" else "fail",
    stringsAsFactors = FALSE
  )
}

# Stopped attempts of a run: the restart messages BayesTools stores among the
# fit's warnings (older run files have them among the signalled ones). A run
# without a posterior states its earlier attempts only when its warnings or
# error carry such messages; otherwise their number is unknown (NA).
restart_pattern <- "fitting attempt [0-9]+ failed and was restarted"

stopped_attempts <- function(run) {

  text  <- unique(c(run[["warnings"]], run[["fit_warnings"]], run[["error"]]))
  count <- sum(grepl(restart_pattern, text))
  if (!is.null(run[["error"]]) && count == 0L) return(NA_integer_)

  return(count)
}

# Whether a run has no posterior because the selection-normalizer diagnostic
# stopped every attempt.
stopped_by_diagnostic <- function(run) {

  !isTRUE(run[["missing"]]) && !is.null(run[["error"]]) && grepl(normalizer_diagnostic, run[["error"]], fixed = TRUE)
}

# Per model and arm over all runs, stopped ones included (the planned runs and
# the runs a replacement took the place of): runs stopped by the diagnostic and
# wall hours.
arm_totals <- function(runs) {

  if (length(runs) == 0L) return(NULL)
  rows <- data.frame(
    model        = vapply(runs, function(run) run[["job"]][["model"]], character(1)),
    arm          = vapply(runs, function(run) run[["job"]][["arm"]], character(1)),
    stopped_runs = vapply(runs, function(run) as.integer(stopped_by_diagnostic(run)), integer(1)),
    wall_hours   = vapply(runs, function(run) as.numeric(run[["wall"]]) / 3600, numeric(1)),
    stringsAsFactors = FALSE
  )

  return(stats::aggregate(rows[c("stopped_runs", "wall_hours")], rows[c("model", "arm")], sum, na.rm = TRUE))
}

# One row per completed run: its record, warnings, stopped attempts and efficiency.
run_rows <- function(runs) {

  rows <- lapply(Filter(function(run) is.null(run[["error"]]), runs), function(run) {
    record <- run[["record"]]
    global <- global_parameters(names(run[["ess"]]), run[["roles"]])
    ess    <- run[["ess"]][global]
    ess    <- ess[is.finite(ess) & ess > 0]
    indicators <- run[["indicator_ess"]]
    data.frame(
      model = run[["job"]][["model"]], arm = run[["job"]][["arm"]], seed = run[["job"]][["seed"]],
      blocks = if (is.null(run[["blocks"]])) NA_character_ else run[["blocks"]],
      control = if (is.null(run[["control"]])) NA_character_ else run[["control"]],
      jags_seconds = as.numeric(run[["jags_seconds"]]),
      samplers = record[["samplers"]], fallback_samplers = record[["fallback_samplers"]],
      fallback_few_proposals = record[["fallback_few_proposals"]],
      fallback_surrogate_exact = record[["fallback_surrogate_exact"]],
      fallback_low_acceptance = record[["fallback_low_acceptance"]],
      acceptance_rate = record[["acceptance_rate"]], min_sampler_rate = record[["min_sampler_rate"]],
      complete = record[["complete"]],
      low_acceptance_warning = "RoBMA_selection_sampler_low_acceptance" %in% run[["warning_classes"]],
      warnings = length(run[["warnings"]]),
      stopped_attempts = stopped_attempts(run),
      wall = as.numeric(run[["wall"]]),
      min_ess = if (length(ess)) min(ess) else NA_real_,
      min_ess_quantity = if (length(ess)) names(ess)[which.min(ess)] else NA_character_,
      ess_per_second = if (length(ess)) min(ess) / as.numeric(run[["jags_seconds"]]) else NA_real_,
      indicator_min_ess = if (length(indicators)) min(indicators) else NA_real_,
      indicator_ess_per_second = if (length(indicators)) min(indicators) / as.numeric(run[["jags_seconds"]]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) return(NULL)

  return(do.call(rbind, rows))
}

# Means over the completed runs per model and arm, and the ratio of the
# 'coarse' to the 'default' arm of the JAGS seconds and the effective draws per
# second; the stopped attempts (total) and restarted runs over the completed
# runs; and from 'totals' (arm_totals() over all runs) the stopped runs, the wall
# hours and their coarse / default ratio.
efficiency_table <- function(rows, totals = NULL) {

  metrics <- c("jags_seconds", "samplers", "fallback_samplers", "fallback_few_proposals",
               "fallback_surrogate_exact", "fallback_low_acceptance", "acceptance_rate", "min_sampler_rate",
               "warnings", "low_acceptance_warning", "min_ess", "ess_per_second", "indicator_min_ess",
               "indicator_ess_per_second")
  if (is.null(totals)) {
    totals <- data.frame(model = character(), arm = character(), stopped_runs = integer(), wall_hours = numeric())
  }
  # one row per model and arm with any attempted run; NA means where none completed
  table <- unique(rbind(if (!is.null(rows)) rows[c("model", "arm")], totals[c("model", "arm")]))
  if (nrow(table) == 0L) return(NULL)
  if (!is.null(rows)) {
    restarts <- data.frame(stopped_attempts = rows[["stopped_attempts"]],
                           runs_restarted = as.integer(rows[["stopped_attempts"]] > 0L))
    table <- merge(table, stats::aggregate(rows[metrics], rows[c("model", "arm")], function(x) mean(as.numeric(x))),
                   by = c("model", "arm"), all.x = TRUE)
    table <- merge(table, stats::aggregate(restarts, rows[c("model", "arm")], sum), by = c("model", "arm"),
                   all.x = TRUE)
  } else {
    for (name in c(metrics, "stopped_attempts", "runs_restarted")) table[[name]] <- NA_real_
  }
  table <- merge(table, totals, by = c("model", "arm"), all.x = TRUE)
  wall_default <- table[["wall_hours"]][match(paste(table[["model"]], "default"),
                                              paste(table[["model"]], table[["arm"]]))]
  table[["wall_ratio"]] <- table[["wall_hours"]] / wall_default
  default <- table[table[["arm"]] == "default", c("model", "jags_seconds", "ess_per_second",
                                                  "indicator_ess_per_second")]
  names(default)[-1L] <- paste0(names(default)[-1L], "_default")
  table <- merge(table, default, by = "model", all.x = TRUE)
  table[["time_ratio"]] <- table[["jags_seconds"]] / table[["jags_seconds_default"]]
  table[["speedup"]] <- table[["ess_per_second"]] / table[["ess_per_second_default"]]
  table[["indicator_speedup"]] <- table[["indicator_ess_per_second"]] / table[["indicator_ess_per_second_default"]]
  settings <- if (!is.null(rows)) unique(rows[c("model", "blocks", "control")]) else
    data.frame(model = unique(table[["model"]]), blocks = NA_character_, control = NA_character_)
  settings <- settings[!duplicated(settings[["model"]]), , drop = FALSE]
  table <- merge(table, settings, by = "model", all.x = TRUE)
  table <- table[order(table[["model"]], table[["arm"]]), c("model", "arm", "blocks", "control", metrics,
                                                            "stopped_attempts", "runs_restarted", "stopped_runs",
                                                            "wall_hours", "time_ratio", "wall_ratio", "speedup",
                                                            "indicator_speedup")]
  rownames(table) <- NULL

  return(table)
}

# One row per run; 'replaced' (the table of replace_stopped_runs()) gives the
# seed that took the place of a stopped run ('replaced_by', NA otherwise).
run_summary <- function(runs, replaced = NULL) {

  do.call(rbind, lapply(runs, function(run) {
    run <- complete_run(run)
    replaced_by <- if (is.null(replaced)) integer() else replaced[["replacement_seed"]][
      replaced[["model"]] == run[["job"]][["model"]] & replaced[["arm"]] == run[["job"]][["arm"]] &
        replaced[["stopped_seed"]] == run[["job"]][["seed"]]]
    data.frame(
      model = run[["job"]][["model"]], arm = run[["job"]][["arm"]], seed = run[["job"]][["seed"]],
      replaced_by = if (length(replaced_by)) as.integer(replaced_by[[1L]]) else NA_integer_,
      wall = round(as.numeric(run[["wall"]]), 1), jags_seconds = round(as.numeric(run[["jags_seconds"]]), 1),
      stopped_attempts = if (isTRUE(run[["missing"]])) NA_integer_ else stopped_attempts(run),
      samplers = if (is.null(run[["samplers"]])) "" else
        paste(names(run[["samplers"]]), run[["samplers"]], sep = ":", collapse = " "),
      build = if (is.null(run[["build"]])) "" else build_text(run[["build"]]),
      error = if (is.null(run[["error"]])) "" else run[["error"]],
      warnings = paste(c(run[["warnings"]], run[["fit_warnings"]]), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
}


# Drivers ----------------------------------------------------------------------------

# The replacement rule: a run stopped by the selection-normalizer diagnostic is
# replaced by the next unused seed of its arm (base + seeds + 1, + 2, ...), at
# most three replacements per model and arm. Returns the jobs with the
# replacements in place of the stopped runs and the table of replaced runs.
normalizer_diagnostic <- "Selection normalizer was rejected by diagnostics"

replace_stopped_runs <- function(jobs, base, seeds, specs, workers, dir, execute, limit = 3L) {

  replaced <- data.frame(model = character(), arm = character(), stopped_seed = integer(),
                         replacement_seed = integer(), stringsAsFactors = FALSE)
  repeat {
    runs    <- read_runs(jobs, dir)
    stopped <- which(vapply(runs, stopped_by_diagnostic, logical(1)))
    new_jobs <- list()
    for (k in stopped) {
      job  <- jobs[[k]]
      used <- sum(replaced[["model"]] == job[["model"]] & replaced[["arm"]] == job[["arm"]])
      if (used >= limit) next
      replacement <- job
      replacement[["seed"]] <- base[[job[["arm"]]]] + seeds + used + 1L
      replaced[nrow(replaced) + 1L, ] <- list(job[["model"]], job[["arm"]], job[["seed"]], replacement[["seed"]])
      jobs[[k]] <- replacement
      new_jobs[[length(new_jobs) + 1L]] <- replacement
    }
    if (length(new_jobs) == 0L) break
    execute(new_jobs, specs, workers)
  }

  return(list(jobs = jobs, replaced = replaced))
}

empty_welch_table <- function() {

  data.frame(model = character(), quantity = character(), mean_coarse = numeric(), mean_default = numeric(),
             se = numeric(), t = numeric(), df = numeric(), p_beyond_2 = numeric(), stringsAsFactors = FALSE)
}

# The wider testing (--part=wider or --part=rerun): runs the missing jobs with
# 'execute', reads every planned run from 'dir', and writes the run table, the
# Welch table, the decisions and the efficiency table there.
wider_driver <- function(part, models, seeds, workers, specs, dir = out_dir, execute = run_jobs) {

  # seeds: arm 'coarse' 1000 + i, arm 'default' 2000 + i; the rerun uses fresh
  # seeds (5000 + i, 6000 + i)
  base <- if (identical(part, "rerun")) c(coarse = 5000L, default = 6000L) else c(coarse = 1000L, default = 2000L)
  failing <- NULL
  if (identical(part, "rerun")) {
    first     <- utils::read.csv(file.path(dir, "wider.csv"), stringsAsFactors = FALSE)
    failing   <- first[(abs(first[["t"]]) > 2) %in% TRUE, c("model", "quantity")]
    decisions <- utils::read.csv(file.path(dir, "wider-decisions.csv"), stringsAsFactors = FALSE)
    models    <- intersect(models, decisions[["model"]][decisions[["decision"]] == "fail"])
    cat("rerun models:", paste(models, collapse = ", "), "\n")
    if (length(models) == 0L) {
      cat("Nothing to rerun: no model failed the first stage.\n")
      utils::write.csv(empty_welch_table(), file.path(dir, paste0(part, ".csv")), row.names = FALSE)
      utils::write.csv(data.frame(model = character(), decision = character()),
                       file.path(dir, paste0(part, "-decisions.csv")), row.names = FALSE)
      # no file of an earlier rerun in this directory may contradict that
      utils::write.csv(data.frame(model = character(), arm = character(), seed = integer(),
                                  replaced_by = integer(), wall = numeric(), jags_seconds = numeric(),
                                  stopped_attempts = integer(), samplers = character(), build = character(),
                                  error = character(), warnings = character()),
                       file.path(dir, paste0(part, "-runs.csv")), row.names = FALSE)
      utils::write.csv(data.frame(model = character(), arm = character(), stopped_seed = integer(),
                                  replacement_seed = integer()),
                       file.path(dir, paste0(part, "-replaced.csv")), row.names = FALSE)
      unlink(file.path(dir, paste0(part, c("-efficiency.csv", "-efficiency-runs.csv"))))
      return(invisible(NULL))
    }
  }
  # arms interleaved by seed
  jobs <- list()
  for (i in seq_len(seeds)) for (model in models) for (arm in c("coarse", "default")) {
    jobs[[length(jobs) + 1L]] <- list(part = part, model = model, arm = arm, seed = base[[arm]] + i)
  }
  execute(jobs, specs, workers)
  replacement <- replace_stopped_runs(jobs, base, seeds, specs, workers, dir, execute)
  jobs        <- replacement[["jobs"]]
  replaced    <- replacement[["replaced"]]
  utils::write.csv(replaced, file.path(dir, paste0(part, "-replaced.csv")), row.names = FALSE)
  if (nrow(replaced) > 0L) {
    cat("\nRuns stopped by the selection-normalizer diagnostic and their replacements:\n")
    print(replaced, row.names = FALSE)
  }
  runs    <- read_runs(jobs, dir)
  counts  <- run_counts(runs)
  summary <- run_summary(runs)
  # the stopped runs that a replacement took the place of: listed and in the
  # totals, not in the counts or the comparison
  stopped_jobs <- lapply(seq_len(nrow(replaced)), function(i) {
    list(part = part, model = replaced[["model"]][[i]], arm = replaced[["arm"]][[i]],
         seed = replaced[["stopped_seed"]][[i]])
  })
  stopped_runs <- read_runs(stopped_jobs, dir)
  utils::write.csv(run_summary(c(runs, stopped_runs), replaced), file.path(dir, paste0(part, "-runs.csv")),
                   row.names = FALSE)
  cat("\nRuns with errors:", sum(nzchar(summary[["error"]])), "\n")
  if (any(nzchar(summary[["error"]]))) print(summary[nzchar(summary[["error"]]), c("model", "arm", "seed", "error")])
  builds <- unique(summary[["build"]][nzchar(summary[["build"]])])
  cat("\nBuilds of the runs:\n")
  cat(builds, sep = "\n")
  cat("\nSamplers per model and arm:\n")
  print(unique(summary[, c("model", "arm", "samplers")]), row.names = FALSE)

  welch <- lapply(models, function(model) {
    quantities <- if (!is.null(failing)) failing[["quantity"]][failing[["model"]] == model] else NULL
    welch_table(runs, model, quantities)
  })
  tables  <- lapply(welch, `[[`, "table")
  missing <- do.call(rbind, lapply(welch, `[[`, "missing"))
  if (!is.null(missing) && nrow(missing) > 0L) {
    cat("\nQuantities missing from a completed run (no statistic; the model is incomplete):\n")
    print(missing, row.names = FALSE)
  }
  table <- do.call(rbind, tables)
  if (is.null(table)) table <- empty_welch_table()
  utils::write.csv(table, file.path(dir, paste0(part, ".csv")), row.names = FALSE)
  decisions <- do.call(rbind, lapply(seq_along(models), function(i) decide(tables[[i]], models[[i]], counts)))
  utils::write.csv(decisions, file.path(dir, paste0(part, "-decisions.csv")), row.names = FALSE)
  cat("\nRuns per model and arm (planned, completed, failed, missing):\n"); print(counts, row.names = FALSE)
  cat("\nDecision per model:\n"); print(decisions, row.names = FALSE)
  incomplete <- decisions[["model"]][decisions[["decision"]] == "incomplete"]
  if (length(incomplete) > 0L) {
    cat("INCOMPLETE: failed or missing runs, fewer than two completed runs in an arm, a quantity missing from a",
        "completed run, or no comparable quantity in model(s)", paste(incomplete, collapse = ", "),
        "; the comparison is not a pass.\n")
  }
  cat("\nQuantities with |t| > 2:\n")
  print(table[(abs(table[["t"]]) > 2) %in% TRUE, , drop = FALSE], row.names = FALSE)

  rows <- run_rows(runs)
  if (!is.null(rows)) utils::write.csv(rows, file.path(dir, paste0(part, "-efficiency-runs.csv")), row.names = FALSE)
  # all runs of an arm for the totals: the planned and replacement runs and the
  # stopped runs that a replacement took the place of
  efficiency <- efficiency_table(rows, arm_totals(c(runs, stopped_runs)))
  if (!is.null(efficiency)) {
    utils::write.csv(efficiency, file.path(dir, paste0(part, "-efficiency.csv")), row.names = FALSE)
    cat("\nRecord, warnings and efficiency (means over completed runs; stopped_attempts and runs_restarted",
        "summed over completed runs; stopped_runs and wall_hours over all runs, stopped ones included;",
        "ratios coarse / default):\n")
    print(efficiency, row.names = FALSE, digits = 4)
    ensembles <- intersect(c("RoBMA.mv", "RoBMA.mv_vi", "best_ensemble"), efficiency[["model"]])
    if (length(ensembles) > 0L) {
      cat("\nAllocation block (model-averaged models: 'default' declares the block; 'coarse' declares it",
          "unless the corrected sampler takes the allocation nodes):\n")
      print(efficiency[efficiency[["model"]] %in% ensembles,
                       c("model", "arm", "jags_seconds", "ess_per_second", "indicator_ess_per_second",
                         "time_ratio", "speedup", "indicator_speedup")], row.names = FALSE, digits = 4)
    }
  }

  invisible(decisions)
}


# Cross-build parts: one R process per arm --------------------------------------------

# Libraries of an arm: this session's, or the other build's library in place of
# the first entry.
arm_libraries <- function(arm) {

  if (identical(arm, "this")) .libPaths() else c(normalizePath(other_lib, winslash = "/", mustWork = TRUE),
                                                 .libPaths()[-1L])
}

# Runs one job in a new R process with the arm's libraries (--part=job).
run_process <- function(job, dir = out_dir) {

  job_file <- cross_file(job, dir, ".job.rds")
  saveRDS(job, job_file)
  libraries <- paste(arm_libraries(job[["arm"]]), collapse = .Platform$path.sep)
  previous  <- Sys.getenv(c("R_LIBS", "R_LIBS_USER"), unset = NA)
  Sys.setenv(R_LIBS = libraries, R_LIBS_USER = libraries)
  on.exit({
    for (name in names(previous)) {
      if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
    }
  }, add = TRUE)
  log_file <- sub("\\.job\\.rds$", ".log", job_file)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(script_path), "--part=job", paste0("--job=", shQuote(job_file)),
                      paste0("--out=", shQuote(dir))), stdout = log_file, stderr = log_file)
  cat(sprintf("%s %s %s seed %d: exit %d\n", job[["part"]], job[["model"]], job[["arm"]], job[["seed"]], status))

  invisible(status)
}

# A fit of a cross-build job in this process: the options the build has, the
# draws, the sampler assignment and the timing.
fit_job <- function(job) {

  options <- job[["options"]]
  options <- options[names(options) %in% names(RoBMA::RoBMA.options())]
  if (length(options)) do.call(RoBMA::RoBMA.options, options)
  spec <- model_specs_cross()[[job[["model"]]]]
  started <- Sys.time()
  fit <- tryCatch(suppressWarnings(do.call(getExportedValue("RoBMA", spec[["constructor"]]), c(spec[["args"]],
    job[["fit"]], list(seed = job[["seed"]], parallel = FALSE, silent = TRUE)))), error = function(e) e)
  result <- list(job = job, build = build_info(), options_set = options,
                 wall = as.numeric(difftime(Sys.time(), started, units = "secs")))
  if (inherits(fit, "error") || is.null(fit[["fit"]][["mcmc"]])) {
    result[["error"]] <- if (inherits(fit, "error")) conditionMessage(fit) else "no posterior draws"
    return(result)
  }
  result[["draws"]]    <- as.matrix(fit[["fit"]][["mcmc"]])
  result[["samplers"]] <- fit[["fit"]][["samplers"]][, c("Node", "Sampler")]
  result[["jags_seconds"]] <- as.numeric(fit[["fit"]][["timetaken"]], units = "secs")
  result[["iterations"]]   <- job[["fit"]][["chains"]] *
    (job[["fit"]][["adapt"]] + job[["fit"]][["burnin"]] + job[["fit"]][["sample"]])

  return(result)
}

model_specs_cross <- function() {

  data("dat.assink2016", package = "metadat", envir = environment())
  assink <- dat.assink2016
  V_full <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid, rho = c(.7, .5), data = assink)
  list(
    bselmodel    = list(constructor = "bselmodel", args = list(yi = quote(yi), vi = quote(vi), data = assink,
                                                                measure = "SMD")),
    bselmodel.mv = list(constructor = "bselmodel.mv", args = list(yi = quote(yi), V = V_full, data = assink,
                                                                   measure = "SMD", random = ~ 1 | study / esid)),
    RoBMA        = list(constructor = "RoBMA", args = list(yi = quote(yi), vi = quote(vi), data = assink,
                                                            measure = "SMD")),
    RoBMA.mv     = list(constructor = "RoBMA.mv", args = list(yi = quote(yi), V = V_full, data = assink,
                                                               measure = "SMD", random = ~ 1 | study / esid)),
    brma.mv      = list(constructor = "brma.mv", args = list(yi = quote(yi), V = V_full, data = assink,
                                                              measure = "SMD", random = ~ 1 | study / esid)),
    RoBMA.mv_vi  = list(constructor = "RoBMA.mv", args = list(yi = quote(yi), V = quote(vi), data = assink,
                                                               measure = "SMD", random = ~ 1 | study / esid))
  )
}

cross_file <- function(job, dir = out_dir, extension = ".rds") {

  variant <- if (is.null(job[["variant"]])) "" else paste0("-", job[["variant"]])
  file.path(dir, "runs", sprintf("%s-%s-%s%s-%03d%s", job[["part"]], job[["model"]], job[["arm"]], variant,
                                 job[["seed"]], extension))
}

read_cross <- function(job, dir = out_dir) {

  file <- cross_file(job, dir)
  if (!file.exists(file)) return(list(job = job, error = "missing run file", missing = TRUE))
  readRDS(file)
}

# The cost of the diagnostic (--part=cost).
cost_driver <- function(seeds, dir = out_dir, execute = run_process) {

  fit  <- list(chains = 3L, adapt = 500L, burnin = 2000L, sample = 4000L)
  jobs <- list()
  for (i in seq_len(seeds)) for (arm in c("this", "other")) {
    jobs[[length(jobs) + 1L]] <- list(part = "cost", model = "bselmodel.mv", arm = arm, seed = 4000L + i, fit = fit,
      options = list(sampler.blocks = FALSE, selection.sampler = "coarse_corrected",
                     selection.coarse_min_acceptance = 0))
  }
  for (job in jobs) if (!file.exists(cross_file(job, dir))) execute(job, dir)
  runs <- lapply(jobs, read_cross, dir = dir)
  rows <- do.call(rbind, lapply(seq_len(seeds), function(i) {
    this  <- runs[[2L * i - 1L]]
    other <- runs[[2L * i]]
    ok <- is.null(this[["error"]]) && is.null(other[["error"]])
    seconds <- function(run) if (is.null(run[["error"]])) 1000 * run[["jags_seconds"]] / run[["iterations"]] else NA_real_
    data.frame(seed = this[["job"]][["seed"]],
               status = if (ok) "complete" else "incomplete",
               identical_draws = ok && identical(this[["draws"]], other[["draws"]]),
               identical_samplers = ok && identical(this[["samplers"]], other[["samplers"]]),
               max_abs_draw_difference = if (ok && identical(dim(this[["draws"]]), dim(other[["draws"]])))
                 max(abs(this[["draws"]] - other[["draws"]])) else NA_real_,
               seconds_per_1000_this = seconds(this), seconds_per_1000_other = seconds(other),
               options_this = paste(names(this[["options_set"]]), collapse = " "),
               options_other = paste(names(other[["options_set"]]), collapse = " "),
               build_this = if (is.null(this[["build"]])) "" else build_text(this[["build"]]),
               build_other = if (is.null(other[["build"]])) "" else build_text(other[["build"]]),
               stringsAsFactors = FALSE)
  }))
  rows[["time_ratio"]] <- rows[["seconds_per_1000_this"]] / rows[["seconds_per_1000_other"]]
  utils::write.csv(rows, file.path(dir, "cost.csv"), row.names = FALSE)
  cat("\nCost of the diagnostic (this build, fallback off, against the other build):\n")
  print(rows, row.names = FALSE, digits = 4)
  cat("mean time ratio:", mean(rows[["time_ratio"]]), "\n")
  if (any(rows[["status"]] != "complete")) cat("INCOMPLETE: a run failed or is missing.\n")
  if (any(rows[["status"]] == "complete" & !rows[["identical_draws"]])) {
    cat("FINDING: the draws of the builds differ.\n")
  }

  invisible(rows)
}

# The default identity (--part=identity).
identity_driver <- function(dir = out_dir, execute = run_process) {

  fit  <- list(chains = 2L, adapt = 500L, burnin = 500L, sample = 500L)
  jobs <- list()
  # the other build with its defaults against this build with selection.sampler = "default"
  for (model in c("bselmodel", "bselmodel.mv", "RoBMA.mv", "brma.mv")) {
    jobs[[length(jobs) + 1L]] <- list(part = "identity", model = model, arm = "this", seed = 7001L, fit = fit,
                                      options = list(selection.sampler = "default"), pair = "builds")
    jobs[[length(jobs) + 1L]] <- list(part = "identity", model = model, arm = "other", seed = 7001L, fit = fit,
                                      options = list(), pair = "builds")
  }
  # this build alone: no coarse-capable selection child (RoBMA.mv on sampling
  # variances has single-estimate blocks only, so it also keeps its allocation
  # block), both option values
  for (model in c("bselmodel", "RoBMA", "brma.mv", "RoBMA.mv_vi")) {
    for (sampler in c("default", "coarse_corrected")) {
      jobs[[length(jobs) + 1L]] <- list(part = "identity", model = model, arm = "this", seed = 7001L, fit = fit,
                                        options = list(selection.sampler = sampler), pair = "options",
                                        variant = sampler)
    }
  }
  for (job in jobs) if (!file.exists(cross_file(job, dir))) execute(job, dir)
  runs <- lapply(jobs, read_cross, dir = dir)
  rows <- do.call(rbind, lapply(seq(1L, length(jobs), by = 2L), function(i) {
    a <- runs[[i]]
    b <- runs[[i + 1L]]
    ok <- is.null(a[["error"]]) && is.null(b[["error"]])
    data.frame(model = a[["job"]][["model"]], pair = a[["job"]][["pair"]],
               first = paste(a[["job"]][["arm"]], paste(names(a[["job"]][["options"]]),
                                                       unlist(a[["job"]][["options"]]), sep = "=", collapse = " ")),
               second = paste(b[["job"]][["arm"]], paste(names(b[["job"]][["options"]]),
                                                        unlist(b[["job"]][["options"]]), sep = "=", collapse = " ")),
               status = if (ok) "complete" else "incomplete",
               identical_samplers = ok && identical(a[["samplers"]], b[["samplers"]]),
               identical_draws = ok && identical(a[["draws"]], b[["draws"]]),
               corrected_samplers = if (ok) sum(a[["samplers"]][["Sampler"]] == "RoBMA::CoarseCorrectedSlice") else NA_integer_,
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(rows, file.path(dir, "identity.csv"), row.names = FALSE)
  cat("\nDefault identity:\n")
  print(rows, row.names = FALSE)
  if (any(rows[["status"]] != "complete")) cat("INCOMPLETE: a run failed or is missing.\n")
  if (any(rows[["status"]] == "complete" & !(rows[["identical_samplers"]] & rows[["identical_draws"]]))) {
    cat("FINDING: a pair differs.\n")
  }

  invisible(rows)
}


# Self-test ------------------------------------------------------------------------

# Synthetic run files in a temporary directory through the complete drivers;
# every case reads back the written tables and checks the decision, the counts
# and the shapes. Stops on a wrong result or on any error.
selftest <- function() {

  root <- tempfile("verify-coarse-selftest-")
  no_fits <- function(...) invisible(NULL)
  check <- function(ok, what) {
    if (!isTRUE(ok)) stop("Self-test failed: ", what, call. = FALSE)
    cat("ok:", what, "\n")
  }
  # Deterministic run means: run i of an arm is offset by v[i]; the 'default'
  # arm uses the offsets in reverse order, so both arms have the same run values
  # (Welch t = 0) unless shifted.
  offsets <- stats::qnorm(stats::ppoints(10L)) / 5
  record  <- function(samplers = 6) data.frame(sampler = "coarse_corrected", min_acceptance = .9,
    samplers = samplers, fallback_samplers = 1, fallback_few_proposals = 0, fallback_surrogate_exact = 0,
    fallback_low_acceptance = 1, sampling_proposed = 1000, sampling_accepted = 990, acceptance_rate = .99,
    min_sampler_rate = .98, complete = TRUE, stringsAsFactors = FALSE)
  wider_run <- function(job, names, shift = 0, constant = FALSE, i = job[["seed"]] %% 1000L) {
    offset <- if (job[["arm"]] == "coarse") offsets[i] else offsets[11L - i]
    means  <- vapply(seq_along(names), function(j) if (constant) j else offset + j + shift[[min(j, length(shift))]],
                     numeric(1L))
    list(job = job, build = list(RoBMA = "0", BayesTools = "0", RoBMA_dll = "x", BayesTools_dll = "y",
                                 RoBMA_path = "p", library = "l"),
         wall = 1, jags_seconds = if (job[["arm"]] == "coarse") 1 else 2, warnings = character(),
         warning_classes = character(), means = stats::setNames(means, names),
         ess = stats::setNames(rep(100, length(names)), names),
         indicator_ess = c("bias_indicator==1" = 50), roles = stats::setNames(rep("sampled", length(names)), names),
         blocks = "rank_one:3 x2", control = "default",
         record = record(if (job[["arm"]] == "coarse") 6 else 0))
  }
  failed_run <- function(job) list(job = job, wall = 1, error = "synthetic fit failure")
  write_runs <- function(dir, part, make) {
    dir.create(file.path(dir, "runs"), recursive = TRUE, showWarnings = FALSE)
    for (arm in c("coarse", "default")) for (i in seq_len(10L)) {
      seed <- switch(part, wider = c(coarse = 1000L, default = 2000L)[[arm]] + i,
                     rerun = c(coarse = 5000L, default = 6000L)[[arm]] + i)
      job  <- list(part = part, model = "m", arm = arm, seed = seed)
      run  <- make(job, i)
      if (!is.null(run)) saveRDS(run, run_file(job, dir))
    }
  }
  wider_case <- function(label, make, part = "wider", dir = file.path(root, label), execute = no_fits) {
    write_runs(dir, part, make)
    output <- utils::capture.output(wider_driver(part, "m", 10L, 1L, list(m = list()), dir = dir,
                                                 execute = function(jobs, ...) execute(jobs, dir)))
    efficiency <- file.path(dir, paste0(part, "-efficiency.csv"))
    list(decisions = utils::read.csv(file.path(dir, paste0(part, "-decisions.csv")), stringsAsFactors = FALSE),
         table = utils::read.csv(file.path(dir, paste0(part, ".csv")), stringsAsFactors = FALSE),
         replaced = utils::read.csv(file.path(dir, paste0(part, "-replaced.csv")), stringsAsFactors = FALSE),
         runs = utils::read.csv(file.path(dir, paste0(part, "-runs.csv")), stringsAsFactors = FALSE),
         efficiency = if (file.exists(efficiency)) utils::read.csv(efficiency, stringsAsFactors = FALSE) else NULL,
         output = output)
  }
  counts_are <- function(d, coarse, default) {
    all(unlist(d[1L, c("planned_coarse", "completed_coarse", "failed_coarse", "missing_coarse")]) == coarse) &&
      all(unlist(d[1L, c("planned_default", "completed_default", "failed_default", "missing_default")]) == default)
  }

  x <- wider_case("complete", function(job, i) wider_run(job, c("a", "b", "c")))
  check(x$decisions$decision == "pass" && x$decisions$quantities == 3L &&
          counts_are(x$decisions, c(10, 10, 0, 0), c(10, 10, 0, 0)) && nrow(x$efficiency) == 2L &&
          abs(x$efficiency$speedup[x$efficiency$arm == "coarse"] - 2) < 1e-12 &&
          abs(x$efficiency$time_ratio[x$efficiency$arm == "coarse"] - .5) < 1e-12 &&
          x$efficiency$samplers[x$efficiency$arm == "coarse"] == 6 &&
          all(x$efficiency$blocks == "rank_one:3 x2") && all(x$efficiency$control == "default"),
        "complete: pass, 3 quantities, 10 of 10 runs per arm, efficiency rows, ratios, blocks and control")

  # one coarse run restarted twice by BayesTools (messages among the fit's warnings)
  restart <- function(n) sprintf("JAGS fitting attempt %d failed and was restarted: synthetic stop.", n)
  x <- wider_case("restarts", function(job, i) {
    run <- wider_run(job, c("a", "b"))
    if (job$arm == "coarse" && i == 3L) run$fit_warnings <- c(restart(1L), restart(2L))
    run
  })
  coarse  <- x$efficiency[x$efficiency$arm == "coarse", ]
  default <- x$efficiency[x$efficiency$arm == "default", ]
  check(x$decisions$decision == "pass" && coarse$stopped_attempts == 2 && coarse$runs_restarted == 1 &&
          default$stopped_attempts == 0 && default$runs_restarted == 0 &&
          all(x$efficiency$stopped_runs == 0) && all(abs(x$efficiency$wall_hours - 10 / 3600) < 1e-12) &&
          abs(coarse$wall_ratio - 1) < 1e-12 && coarse$warnings == 0 &&
          sum(x$runs$stopped_attempts) == 2,
        "two restart messages in one run: 2 stopped attempts, 1 restarted run, wall hours and ratio")

  # a quantity missing from one completed run: no statistic, incomplete, listed
  x <- wider_case("missing_mean", function(job, i) {
    run <- wider_run(job, c("a", "b"))
    if (job$arm == "default" && i == 4L) run$means <- run$means["a"]
    run
  })
  check(x$decisions$decision == "incomplete" && x$decisions$quantities == 2L &&
          is.na(x$table$t[x$table$quantity == "b"]) && !is.na(x$table$t[x$table$quantity == "a"]) &&
          counts_are(x$decisions, c(10, 10, 0, 0), c(10, 10, 0, 0)) &&
          any(grepl("missing from a completed run", x$output)) && any(grepl("default +2004", x$output)),
        "a mean missing from a completed run: row with NA statistics, incomplete, case listed")

  # an indicator level a run never visited has frequency zero there, not missing
  x <- wider_case("unvisited_level", function(job, i) {
    run <- wider_run(job, c("a", "bias_indicator==1", "bias_indicator==2"))
    if (job$arm == "default" && i == 4L) run$means <- run$means[c("a", "bias_indicator==1")]
    run
  })
  check(x$decisions$decision %in% c("pass", "fail") && !anyNA(x$table$t) &&
          !any(grepl("missing from a completed run", x$output)),
        "an unvisited indicator level: frequency zero, the model stays complete")

  # a stored NA level, or an indicator a run lacks altogether, stays missing
  for (case in list(list(label = "stored_na_level", keep = NULL),
                    list(label = "absent_indicator", keep = "a"))) {
    x <- wider_case(case$label, function(job, i) {
      run <- wider_run(job, c("a", "bias_indicator==1", "bias_indicator==2"))
      if (job$arm == "default" && i == 4L) {
        if (is.null(case$keep)) run$means["bias_indicator==2"] <- NA else run$means <- run$means[case$keep]
      }
      run
    })
    check(x$decisions$decision == "incomplete" && is.na(x$table$t[x$table$quantity == "bias_indicator==2"]) &&
            any(grepl("bias_indicator==2 +default +2004", x$output)),
          paste0(case$label, ": no statistic, incomplete, listed"))
  }

  # a requested rerun quantity that the reruns lack: no statistic, incomplete, listed
  dir <- file.path(root, "rerun_absent")
  x <- wider_case("rerun_absent", function(job, i) {
    wider_run(job, c("a", "b", "c"), shift = if (job$arm == "coarse") c(1, 1, 0) else 0)
  }, dir = dir)
  check(x$decisions$decision == "fail" && setequal(x$table$quantity[abs(x$table$t) > 2], c("a", "b")),
        "stage one fails on 'a' and 'b'")
  x <- wider_case("rerun_absent", function(job, i) wider_run(job, "a"), part = "rerun", dir = dir)
  check(x$decisions$decision == "incomplete" && setequal(x$table$quantity, c("a", "b")) &&
          is.na(x$table$t[x$table$quantity == "b"]) && sum(grepl("m +b +(coarse|default) +[56]0", x$output)) == 20L,
        "a rerun quantity absent from the reruns: incomplete, listed once per run")

  # a rerun after a first stage that passed: nothing to rerun, empty tables
  dir <- file.path(root, "rerun_none")
  x <- wider_case("rerun_none", function(job, i) wider_run(job, c("a", "b")), dir = dir)
  # the files of an earlier rerun with rows in the same directory
  earlier <- list.files(file.path(root, "rerun_absent"), pattern = "^rerun", full.names = TRUE)
  file.copy(earlier, dir, overwrite = TRUE)
  utils::write.csv(data.frame(model = "m", arm = "default", stopped_seed = 6009L, replacement_seed = 6011L),
                   file.path(dir, "rerun-replaced.csv"), row.names = FALSE)
  stale <- file.path(dir, c("rerun-runs.csv", "rerun-replaced.csv", "rerun-efficiency.csv",
                            "rerun-efficiency-runs.csv"))
  check(all(file.exists(stale)) && all(vapply(stale, function(file) nrow(utils::read.csv(file)) > 0L, logical(1))),
        "the files of an earlier rerun with rows are in place")
  output <- utils::capture.output(wider_driver("rerun", "m", 10L, 1L, list(m = list()), dir = dir,
                                               execute = no_fits))
  check(x$decisions$decision == "pass" && any(grepl("Nothing to rerun", output)) &&
          nrow(utils::read.csv(file.path(dir, "rerun.csv"))) == 0L &&
          nrow(utils::read.csv(file.path(dir, "rerun-decisions.csv"))) == 0L &&
          nrow(utils::read.csv(file.path(dir, "rerun-runs.csv"))) == 0L &&
          "replaced_by" %in% names(utils::read.csv(file.path(dir, "rerun-runs.csv"))) &&
          nrow(utils::read.csv(file.path(dir, "rerun-replaced.csv"))) == 0L &&
          !any(file.exists(stale[3:4])),
        "a rerun after an all-pass first stage: nothing to rerun, empty tables, no stale rerun files")

  x <- wider_case("one_failed", function(job, i) {
    if (job$arm == "default" && i == 10L) failed_run(job) else wider_run(job, c("a", "b"))
  })
  check(x$decisions$decision == "incomplete" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 9, 1, 0)) &&
          any(grepl("INCOMPLETE", x$output)) && nrow(x$replaced) == 0L,
        "one failed run: incomplete, default 10 planned / 9 completed / 1 failed, not replaced")

  # a run stopped by the selection-normalizer diagnostic is replaced by the next
  # unused seed of its arm (2011 after the planned 2001-2010)
  stopped_run <- function(job) {
    list(job = job, wall = 1, error = paste0("RUNTIME ERROR:\n", normalizer_diagnostic, ": relative MCSE 0.12"))
  }
  # like run_jobs(): runs only the jobs without a saved run
  write_jobs <- function(make) function(jobs, dir) {
    for (job in jobs) if (!file.exists(run_file(job, dir))) saveRDS(make(job), run_file(job, dir))
  }
  x <- wider_case("replaced", function(job, i) {
    if (job$arm == "default" && i == 9L) stopped_run(job) else wider_run(job, c("a", "b"))
  }, execute = write_jobs(function(job) wider_run(job, c("a", "b"), i = 9L)))
  check(x$decisions$decision == "pass" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 10, 0, 0)) &&
          nrow(x$replaced) == 1L && x$replaced$arm == "default" && x$replaced$stopped_seed == 2009L &&
          x$replaced$replacement_seed == 2011L && any(grepl("replacements", x$output)) &&
          x$efficiency$stopped_runs[x$efficiency$arm == "default"] == 1 &&
          x$efficiency$stopped_runs[x$efficiency$arm == "coarse"] == 0 &&
          abs(x$efficiency$wall_hours[x$efficiency$arm == "default"] - 11 / 3600) < 1e-12 &&
          identical(x$runs$stopped_attempts[x$runs$seed == 2011L], 0L) &&
          identical(x$runs$replaced_by[x$runs$arm == "default" & x$runs$seed == 2009L], 2011L) &&
          all(is.na(x$runs$replaced_by[x$runs$seed != 2009L])),
        "a run stopped by the normalizer diagnostic: replaced by seed 2011, listed, complete; counted in the totals")

  # at most three replacements per model and arm; the fourth stopped run stays
  x <- wider_case("replaced_cap", function(job, i) {
    if (job$arm == "default" && i >= 7L) stopped_run(job) else wider_run(job, c("a", "b"))
  }, execute = write_jobs(stopped_run))
  check(x$decisions$decision == "incomplete" && nrow(x$replaced) == 3L &&
          identical(x$replaced$replacement_seed, 2011:2013) &&
          counts_are(x$decisions, c(10, 10, 0, 0), c(10, 6, 4, 0)),
        "replacements capped at three per arm: the stopped runs remain, incomplete")

  # an arm whose runs and replacements all stopped keeps its efficiency row;
  # the replaced seeds stay in the run listing
  x <- wider_case("arm_stopped", function(job, i) {
    if (job$arm == "default") stopped_run(job) else wider_run(job, c("a", "b"))
  }, execute = write_jobs(stopped_run))
  default <- x$efficiency[x$efficiency$arm == "default", ]
  coarse  <- x$efficiency[x$efficiency$arm == "coarse", ]
  check(x$decisions$decision == "incomplete" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 0, 10, 0)) &&
          nrow(default) == 1L && is.na(default$jags_seconds) && default$stopped_runs == 13 &&
          abs(default$wall_hours - 13 / 3600) < 1e-12 && is.finite(coarse$wall_ratio) &&
          abs(coarse$wall_ratio - 10 / 13) < 1e-12 &&
          identical(x$runs$replaced_by[x$runs$arm == "default" & x$runs$seed == 2001L], 2011L) &&
          sum(x$runs$arm == "default") == 13L,
        "an arm without completed runs: its row with wall hours, a finite wall ratio, replaced seeds listed")

  x <- wider_case("one_missing", function(job, i) {
    if (job$arm == "default" && i == 10L) NULL else wider_run(job, c("a", "b"))
  })
  check(x$decisions$decision == "incomplete" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 9, 0, 1)),
        "one missing run file: incomplete, default 10 planned / 9 completed / 1 missing")

  x <- wider_case("no_success", function(job, i) if (job$arm == "default") failed_run(job) else wider_run(job, "a"))
  check(x$decisions$decision == "incomplete" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 0, 10, 0)) &&
          all(is.na(x$table$t)), "no successful run in an arm: incomplete, rows with NA statistics")

  x <- wider_case("no_quantity", function(job, i) wider_run(job, c("a", "b"), constant = TRUE))
  check(x$decisions$decision == "incomplete" && x$decisions$quantities == 0L && nrow(x$table) == 0L,
        "no comparable quantity: incomplete")

  dir <- file.path(root, "rerun")
  x <- wider_case("rerun", function(job, i) {
    wider_run(job, c("a", "b", "c"), shift = if (job$arm == "coarse") c(0, 1, 0) else 0)
  }, dir = dir)
  check(x$decisions$decision == "fail" && identical(x$table$quantity[abs(x$table$t) > 2], "b"),
        "a shifted quantity: fail on that quantity only")
  x <- wider_case("rerun", function(job, i) wider_run(job, c("a", "b", "c")), part = "rerun", dir = dir)
  check(x$decisions$decision == "pass" && x$decisions$quantities == 1L && all(x$table$quantity == "b"),
        "a single-quantity rerun: pass on 'b' alone")

  # cross-build drivers on synthetic run files
  cross_run <- function(job, draws, seconds = 10) {
    list(job = job, build = list(RoBMA = job[["arm"]], BayesTools = "0", RoBMA_dll = "x", BayesTools_dll = "y",
                                 RoBMA_path = "p", library = job[["arm"]]),
         options_set = job[["options"]], wall = 1, draws = draws,
         samplers = data.frame(Node = c("mu", "tau"), Sampler = c("RoBMA::CoarseCorrectedSlice", "base::RealSlicer"),
                               stringsAsFactors = FALSE),
         jags_seconds = seconds, iterations = 1000)
  }
  draws <- matrix(stats::qnorm(stats::ppoints(20L)), 10L, 2L)
  cost_case <- function(label, change) {
    dir <- file.path(root, label)
    dir.create(file.path(dir, "runs"), recursive = TRUE, showWarnings = FALSE)
    write <- function(job, dir) {
      changed <- change(job)
      saveRDS(cross_run(job, if (changed) draws + 1e-9 else draws, if (changed) 11 else 10), cross_file(job, dir))
    }
    output <- utils::capture.output(cost_driver(3L, dir = dir, execute = write))
    list(rows = utils::read.csv(file.path(dir, "cost.csv"), stringsAsFactors = FALSE), output = output)
  }
  y <- cost_case("cost_same", function(job) FALSE)
  check(nrow(y$rows) == 3L && all(y$rows$identical_draws) && all(abs(y$rows$time_ratio - 1) < 1e-12) &&
          !any(grepl("FINDING", y$output)), "cost: identical draws, time ratio 1")
  y <- cost_case("cost_differs", function(job) job$arm == "other" && job$seed == 4002L)
  check(sum(!y$rows$identical_draws) == 1L && any(grepl("FINDING", y$output)) &&
          abs(y$rows$time_ratio[2L] - 10 / 11) < 1e-12, "cost: a differing pair is reported as a finding")

  identity_case <- function(label, change) {
    dir <- file.path(root, label)
    dir.create(file.path(dir, "runs"), recursive = TRUE, showWarnings = FALSE)
    write <- function(job, dir) saveRDS(cross_run(job, if (change(job)) draws * 2 else draws), cross_file(job, dir))
    output <- utils::capture.output(identity_driver(dir = dir, execute = write))
    list(rows = utils::read.csv(file.path(dir, "identity.csv"), stringsAsFactors = FALSE), output = output)
  }
  z <- identity_case("identity_same", function(job) FALSE)
  check(nrow(z$rows) == 8L && all(z$rows$status == "complete") && all(z$rows$identical_draws) &&
          sum(z$rows$pair == "builds") == 4L && sum(z$rows$pair == "options") == 4L,
        "identity: 4 build pairs and 4 option pairs, all identical")
  z <- identity_case("identity_differs", function(job) job$model == "RoBMA.mv" && job$arm == "other")
  check(sum(!z$rows$identical_draws) == 1L && any(grepl("FINDING", z$output)),
        "identity: a differing pair is reported as a finding")

  unlink(root, recursive = TRUE)
  cat("Self-test passed.\n")
}


# Parts ----------------------------------------------------------------------------

if (identical(part, "selftest")) {
  selftest()
  quit(save = "no", status = 0L)
}

if (identical(part, "job")) {
  job    <- readRDS(opt("job"))
  result <- fit_job(job)
  saveRDS(result, cross_file(job, out_dir))
  cat(build_text(result[["build"]]), "\n")
  if (!is.null(result[["error"]])) cat("ERROR", result[["error"]], "\n")
  quit(save = "no", status = 0L)
}

cat("part:", part, " seeds:", seeds, " workers:", workers, "\n")

if (part %in% c("wider", "rerun")) {
  specs  <- model_specs()
  models <- grep(models_re, names(specs), value = TRUE, perl = TRUE)
  cat("models:", paste(models, collapse = ", "), "\n")
  wider_driver(part, models, seeds, workers, specs)
}

if (identical(part, "cost")) cost_driver(seeds)

if (identical(part, "identity")) identity_driver()
