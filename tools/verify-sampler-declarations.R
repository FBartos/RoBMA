#!/usr/bin/env Rscript
# Verification of the sampler declarations RoBMA passes to BayesTools
# (RoBMA.options() 'sampler.gates' and 'sampler.blocks'): replicated
# comparison against the stock samplers, and a timing comparison. Not part of
# any test profile; uses the installed RoBMA in the parent and the workers.
#
#   Rscript tools/verify-sampler-declarations.R --part=compare --out=<dir>
#       [--workers=6] [--seeds=10] [--models=<regex>]
#   Rscript tools/verify-sampler-declarations.R --part=rerun --out=<dir>
#       [--workers=6] [--seeds=10]
#   Rscript tools/verify-sampler-declarations.R --part=timing --out=<dir>
#       [--seeds=3] [--models=<regex>]
#   Rscript tools/verify-sampler-declarations.R --part=selftest
#
# Every planned run is counted per model and arm (planned, completed, failed,
# missing); a comparison or timing with a failed or missing run, an arm with
# fewer than two completed runs (comparison) or none (timing), or no comparable
# or assessable quantity is reported as 'incomplete', never as a pass, and the
# script still writes its tables. --part=selftest runs the complete comparison
# and timing drivers on synthetic run files in a temporary directory (no fits)
# and stops on a wrong decision, shape or count, or on any error.
#
# Models (dat.assink2016, V from metafor::vcalc() as in the assink2016
# scenario): 'bselmodel' (univariate, heterogeneity), 'bselmodel.mv' (full V,
# ~ 1 | study / esid), 'RoBMA' (univariate bias ensemble), 'RoBMA.mv' (full V,
# allocation), and the stress model 'bselmodel_small' (univariate, studies 1,
# 2 and 4 only: heterogeneity weakly identified under its default prior).
# Package defaults for priors, adaptation and burn-in; 3 chains of 8,000 draws
# (4,000 on the full-V models), serial chains in every fit.
#
# compare: arms 'on' (both options TRUE) and 'off' (both FALSE), --seeds runs
# per arm. Each run gives the mean of every monitored quantity over its draws
# (categorical indicators as the frequency of each level); quantities constant
# in every run are degenerate and dropped. Welch statistic of the run means of
# the two arms per quantity. Decision rule (plan D4): pass when no |t| > 3.5
# and the count of |t| > 2 is at most twice its expectation under the Welch
# t distribution; a failing model is rerun once ('--part=rerun') with fresh
# seeds, judged on its failing quantities only. A fit that stops with an error
# (e.g. a failure to calculate a density) is reported as such. Secondary rows
# (marked 'secondary', judged by the same rule but not part of the decision):
# Welch statistics of the run SDs and of the run frequencies above the pooled
# 90% quantile of the 'off' arm, for every non-indicator quantity; the tail
# frequency is computed from every 10th draw, which each run saves. Runs saved
# before these were recorded give the primary rows only.
#
# timing: arms 'gates' (gates only), 'blocks' (blocks only), 'both', 'neither',
# --seeds runs per arm, arms interleaved within each seed, in this session.
# Seconds per 1,000 iterations (JAGS time over chains x (adapt + burnin +
# sample)) and effective sample size per second, the effective sample size
# being the minimum over the global scalar parameters (effect, moderators,
# heterogeneity and allocation parameters, weights, PET, PEESE; not unit-level
# latents or random-effect coordinates, not indicators, not auxiliaries) of
# coda::effectiveSize() summed over chains. Structurally fixed quantities are
# excluded by their convergence role, which BayesTools declares from the priors
# when the model is fitted (BayesTools::parameter_coordinates(): 'structural'
# for point priors, fixed weights and the like; also 'auxiliary' and
# 'indicator'; 'sampled' and 'derived' coordinates are kept). A
# free parameter whose draws are constant in a run (e.g. a model-averaged
# weight whose branches the chain never visits) has no effective sample size
# there and is reported as not assessable, by name, in the per-run table and in
# the summary; the minimum is taken over the assessable parameters. Runs saved
# without their roles take them from one short fit of the model (roles do not
# depend on the draws), cached as roles-<model>.rds in --out.

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}
part    <- opt("part", "compare")
out_dir <- opt("out", if (identical(part, "selftest")) tempfile("verify-selftest-") else NULL)
workers <- as.integer(opt("workers", "6"))
seeds   <- as.integer(opt("seeds", if (identical(part, "timing")) "3" else "10"))
models_re <- opt("models", ".")
if (is.null(out_dir)) stop("Supply --out=<directory>.", call. = FALSE)
if (!part %in% c("compare", "rerun", "timing", "selftest")) {
  stop("--part must be compare, rerun, timing or selftest.", call. = FALSE)
}
dir.create(file.path(out_dir, "runs"), recursive = TRUE, showWarnings = FALSE)
options(warn = 1, width = 200)

suppressPackageStartupMessages(library(RoBMA))
cat("RoBMA:", getNamespaceInfo(asNamespace("RoBMA"), "path"), as.character(utils::packageVersion("RoBMA")), "\n")
cat("BayesTools:", getNamespaceInfo(asNamespace("BayesTools"), "path"), as.character(utils::packageVersion("BayesTools")), "\n")
for (package in c("RoBMA", "BayesTools")) {
  cat(package, "DLL md5:", unname(tools::md5sum(getLoadedDLLs()[[package]][["path"]])), "\n")
}


# Models -------------------------------------------------------------------------

model_specs <- function() {

  data("dat.assink2016", package = "metadat", envir = environment())
  data  <- dat.assink2016
  V     <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid,
                          rho = c(0.7, 0.5), data = data)
  small <- data[data[["study"]] %in% c(1, 2, 4), ]

  list(
    bselmodel       = list(constructor = "bselmodel", sample = 8000L, full_v = FALSE,
                           args = list(yi = quote(yi), vi = quote(vi), measure = "SMD", data = data)),
    bselmodel.mv    = list(constructor = "bselmodel.mv", sample = 4000L, full_v = TRUE,
                           args = list(yi = quote(yi), V = V, random = ~ 1 | study / esid,
                                       measure = "SMD", data = data)),
    RoBMA           = list(constructor = "RoBMA", sample = 8000L, full_v = FALSE,
                           args = list(yi = quote(yi), vi = quote(vi), measure = "SMD", data = data)),
    RoBMA.mv        = list(constructor = "RoBMA.mv", sample = 4000L, full_v = TRUE,
                           args = list(yi = quote(yi), V = V, random = ~ 1 | study / esid,
                                       measure = "SMD", data = data)),
    bselmodel_small = list(constructor = "bselmodel", sample = 8000L, full_v = FALSE,
                           args = list(yi = quote(yi), vi = quote(vi), measure = "SMD", data = small))
  )
}

arm_options <- list(
  on      = c(gates = TRUE,  blocks = TRUE),
  off     = c(gates = FALSE, blocks = FALSE),
  gates   = c(gates = TRUE,  blocks = FALSE),
  blocks  = c(gates = FALSE, blocks = TRUE),
  both    = c(gates = TRUE,  blocks = TRUE),
  neither = c(gates = FALSE, blocks = FALSE)
)


# One run: fit with the arm's options and summarize the draws. Runs in a worker
# or in this session; returns a list and never throws.
run_job <- function(job, specs, arm_options) {

  spec    <- specs[[job[["model"]]]]
  options <- arm_options[[job[["arm"]]]]
  RoBMA::RoBMA.options(sampler.gates = options[["gates"]], sampler.blocks = options[["blocks"]])
  warnings <- character()
  started  <- Sys.time()
  fit <- tryCatch(withCallingHandlers(
    do.call(getExportedValue("RoBMA", spec[["constructor"]]), c(spec[["args"]], list(
      chains = 3L, sample = spec[["sample"]], seed = job[["seed"]], parallel = FALSE, silent = TRUE
    ))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  wall <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  result <- list(job = job, wall = wall, warnings = unique(warnings), jags_seconds = NA_real_,
                 iterations = NA_real_)
  if (inherits(fit, "error") || is.null(fit[["fit"]][["mcmc"]])) {
    result[["error"]] <- if (inherits(fit, "error")) conditionMessage(fit) else
      paste(c("no posterior draws", fit[["fit"]][["errors"]]), collapse = ": ")
    return(result)
  }

  draws <- fit[["fit"]][["mcmc"]]
  declarations <- attr(fit[["fit"]], "sampler_declarations", exact = TRUE)
  samplers <- fit[["fit"]][["samplers"]]
  result[["jags_seconds"]] <- as.numeric(fit[["fit"]][["timetaken"]], units = "secs")
  result[["iterations"]]   <- 3L * (500L + 2000L + spec[["sample"]])
  result[["samplers"]]     <- if (!is.null(samplers)) table(samplers[["Sampler"]]) else NULL
  result[["declared"]]     <- c(
    gates  = if (is.null(declarations)) 0L else nrow(declarations[["gates"]]),
    blocks = if (is.null(declarations)) 0L else nrow(declarations[["blocks"]])
  )
  result[["means"]] <- quantity_means(draws)
  result[["sds"]]   <- quantity_sds(draws)
  result[["thinned"]] <- thin_draws(draws)
  result[["ess"]]   <- coda::effectiveSize(draws)
  result[["roles"]] <- coordinate_roles(fit)
  result[["fit_warnings"]] <- unique(c(fit[["fit"]][["warnings"]], attr(fit[["fit"]], "warnings", exact = TRUE)))

  return(result)
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

# Run SDs of the monitored non-indicator quantities.
quantity_sds <- function(draws) {

  matrix <- as.matrix(draws)
  matrix <- matrix[, !grepl("indicator", colnames(matrix), fixed = TRUE), drop = FALSE]

  return(apply(matrix, 2L, stats::sd))
}

# Every 10th draw of the monitored non-indicator quantities (chains stacked),
# for the tail frequencies of the secondary comparison.
thin_draws <- function(draws) {

  matrix <- as.matrix(draws)
  matrix <- matrix[, !grepl("indicator", colnames(matrix), fixed = TRUE), drop = FALSE]

  return(matrix[seq(1L, nrow(matrix), by = 10L), , drop = FALSE])
}

# The convergence role of each coordinate, declared from the priors at fit time.
coordinate_roles <- function(fit) {

  coordinates <- BayesTools::parameter_coordinates(fit[["fit"]])

  return(stats::setNames(coordinates[["convergence_role"]], coordinates[["coordinate_name"]]))
}

# The roles of a model from one short fit (they do not depend on the draws).
model_roles <- function(model, specs, dir = out_dir) {

  file <- file.path(dir, paste0("roles-", model, ".rds"))
  if (file.exists(file)) return(readRDS(file))
  spec <- specs[[model]]
  RoBMA::RoBMA.options(sampler.gates = FALSE, sampler.blocks = FALSE)
  fit  <- suppressWarnings(do.call(getExportedValue("RoBMA", spec[["constructor"]]), c(spec[["args"]], list(
    chains = 1L, adapt = 50L, burnin = 50L, sample = 100L, seed = 1L, parallel = FALSE, silent = TRUE
  ))))
  roles <- coordinate_roles(fit)
  saveRDS(roles, file)

  return(roles)
}

# Global scalar parameters for the timing ESS: by monitored name, and, with the
# roles, only sampled coordinates and those derived from them (e.g. the
# component SDs of an allocation), not structural, auxiliary or indicator ones.
global_parameters <- function(names, roles = NULL) {

  keep <- !grepl("indicator", names, fixed = TRUE) &
    !grepl("_xRE_(Zx|COEFx|GROUP_Zx|UNIT_COEFx)", names) &
    !grepl("^prior_par_", names) &
    !grepl("^(theta|gamma|sampling_z|pi|phi)\\[", names)
  if (!is.null(roles)) {
    keep <- keep & names %in% names(roles)[roles %in% c("sampled", "derived")]
  }

  return(names[keep])
}

run_file <- function(job, dir = out_dir) {

  file.path(dir, "runs", sprintf("%s-%s-%s-%03d.rds", job[["part"]], job[["model"]], job[["arm"]], job[["seed"]]))
}

run_jobs <- function(jobs, specs, n_workers) {

  todo <- Filter(function(job) !file.exists(run_file(job)), jobs)
  cat(length(jobs), "jobs,", length(jobs) - length(todo), "done before,", length(todo), "to run\n")
  if (length(todo) == 0L) return(invisible(NULL))
  if (n_workers > 1L) {
    cluster <- parallel::makePSOCKcluster(min(n_workers, length(todo)), rscript_args = "--vanilla")
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(paths) {
      .libPaths(paths)
      Sys.setenv(OMP_NUM_THREADS = "1")
      suppressPackageStartupMessages(library(RoBMA))
      invisible(TRUE)
    }, .libPaths())
    builds <- parallel::clusterCall(cluster, function() {
      vapply(c("RoBMA", "BayesTools"), function(package) {
        unname(tools::md5sum(getLoadedDLLs()[[package]][["path"]]))
      }, character(1))
    })
    parent <- vapply(c("RoBMA", "BayesTools"), function(package) {
      unname(tools::md5sum(getLoadedDLLs()[[package]][["path"]]))
    }, character(1))
    if (!all(vapply(builds, identical, logical(1), parent))) {
      stop("A worker loaded a different RoBMA or BayesTools DLL than this session.", call. = FALSE)
    }
    parallel::clusterExport(cluster, c("run_job", "quantity_means", "quantity_sds", "thin_draws",
                                       "coordinate_roles", "arm_options", "specs"),
                            envir = environment())
    worker <- function(job, file) {
      result <- run_job(job, specs, arm_options)
      saveRDS(result, file)
      sprintf("%s %s %s seed %d: %.0f s%s", job[["part"]], job[["model"]], job[["arm"]], job[["seed"]],
              result[["wall"]], if (!is.null(result[["error"]])) paste(" ERROR", result[["error"]]) else "")
    }
    environment(worker) <- globalenv()
    parallel::clusterExport(cluster, "worker", envir = environment())
    files <- vapply(todo, run_file, character(1))
    log <- parallel::clusterApplyLB(cluster, seq_along(todo), function(i, todo, files) {
      worker(todo[[i]], files[[i]])
    }, todo = todo, files = files)
    cat(unlist(log), sep = "\n")
  } else {
    for (job in todo) {
      result <- run_job(job, specs, arm_options)
      saveRDS(result, run_file(job))
      cat(sprintf("%s %s %s seed %d: %.0f s%s\n", job[["part"]], job[["model"]], job[["arm"]], job[["seed"]],
                  result[["wall"]], if (!is.null(result[["error"]])) paste(" ERROR", result[["error"]]) else ""))
      flush.console()
    }
  }
}

# The fields every run record has; a failed or missing run has NA timings.
complete_run <- function(run) {

  defaults <- list(wall = NA_real_, jags_seconds = NA_real_, iterations = NA_real_, warnings = character())
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

# A run-by-quantity matrix (one row per run, one column per name) for zero, one
# or several names; 'extract' returns a run's named values.
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
  names_all <- unique(unlist(lapply(runs, function(run) names(run[["means"]]))))
  if (!is.null(quantities)) {
    names_all <- intersect(names_all, quantities)
  }
  # a level of an indicator absent from a run has frequency zero there
  means <- run_matrix(runs, names_all, function(run, names) {
    values <- run[["means"]][names]
    values[is.na(values) & grepl("==", names, fixed = TRUE)] <- 0
    values
  })
  table <- welch_rows(means, arms, model, "mean", "primary")

  # secondary: run SDs and tail frequencies, when the runs recorded them and
  # each arm has at least two completed runs
  if (sum(arms == "on") >= 2L && sum(arms == "off") >= 2L &&
      all(vapply(runs, function(run) !is.null(run[["sds"]]) && !is.null(run[["thinned"]]), logical(1)))) {
    names_sd <- Reduce(intersect, lapply(runs, function(run) names(run[["sds"]])))
    if (!is.null(quantities)) names_sd <- intersect(names_sd, quantities)
    if (length(names_sd) > 0L) {
      sds <- run_matrix(runs, names_sd, function(run, names) run[["sds"]][names])
      # the tail: above the pooled 90% quantile of the 'off' arm
      off_draws <- do.call(rbind, lapply(runs[arms == "off"], function(run) run[["thinned"]][, names_sd, drop = FALSE]))
      cut <- apply(off_draws, 2L, stats::quantile, probs = .9, names = FALSE)
      tails <- run_matrix(runs, names_sd, function(run, names) {
        colMeans(sweep(run[["thinned"]][, names, drop = FALSE], 2L, cut, ">"))
      })
      table <- rbind(table, welch_rows(sds, arms, model, "sd", "secondary"),
                     welch_rows(tails, arms, model, "tail90", "secondary"))
    }
  }

  return(table)
}

# Welch statistics of the per-run values of each column (runs in rows). With
# fewer than two completed runs in an arm there is no statistic: the rows carry
# NA statistics (the decision is then 'incomplete').
welch_rows <- function(values, arms, model, statistic, kind) {

  row <- function(name, on, off, se, t, df) {
    data.frame(model = model, kind = kind, statistic = statistic, quantity = name,
               mean_on = if (length(on) > 0L) mean(on) else NA_real_,
               mean_off = if (length(off) > 0L) mean(off) else NA_real_,
               se = se, t = t, df = df,
               p_beyond_2 = if (is.na(df)) NA_real_ else 2 * stats::pt(-2, df), stringsAsFactors = FALSE)
  }
  testable <- sum(arms == "on") >= 2L && sum(arms == "off") >= 2L
  rows <- lapply(colnames(values), function(name) {
    on  <- values[arms == "on", name]
    off <- values[arms == "off", name]
    if (!testable) return(row(name, on, off, NA_real_, NA_real_, NA_real_))
    if (anyNA(c(on, off))) return(NULL)
    # degenerate: the same value in every run of both arms
    if (stats::var(c(on, off)) == 0) return(NULL)
    se2 <- stats::var(on) / length(on) + stats::var(off) / length(off)
    if (se2 > 0) {
      df <- se2^2 / ((stats::var(on) / length(on))^2 / (length(on) - 1) +
                     (stats::var(off) / length(off))^2 / (length(off) - 1))
      t  <- (mean(on) - mean(off)) / sqrt(se2)
    } else {
      # constant within each arm but different between them
      df <- length(on) + length(off) - 2
      t  <- Inf
    }
    row(name, on, off, sqrt(se2), t, df)
  })

  do.call(rbind, rows)
}

# The decision of one model and kind: 'incomplete' when any planned run of the
# model failed or is missing, an arm has fewer than two completed runs, or no
# quantity could be compared; otherwise 'pass' or 'fail' by the D4 rule.
decide <- function(table, model, kind, counts) {

  counts <- counts[counts[["model"]] == model, , drop = FALSE]
  arm_count <- function(arm, field) {
    value <- counts[[field]][counts[["arm"]] == arm]
    if (length(value) == 0L) 0L else as.integer(value)
  }
  t <- if (is.null(table)) numeric() else table[["t"]]
  p <- if (is.null(table)) numeric() else table[["p_beyond_2"]]
  complete <- sum(counts[["failed"]]) == 0L && sum(counts[["missing"]]) == 0L &&
    all(vapply(c("on", "off"), arm_count, integer(1), field = "completed") >= 2L) &&
    length(t) > 0L && !anyNA(t)
  rule <- complete && !any(abs(t) > 3.5) && sum(abs(t) > 2) <= 2 * sum(p)
  data.frame(
    model     = model,
    kind      = kind,
    planned_on = arm_count("on", "planned"), completed_on = arm_count("on", "completed"),
    failed_on = arm_count("on", "failed"), missing_on = arm_count("on", "missing"),
    planned_off = arm_count("off", "planned"), completed_off = arm_count("off", "completed"),
    failed_off = arm_count("off", "failed"), missing_off = arm_count("off", "missing"),
    quantities = length(t),
    max_abs_t = if (any(!is.na(t))) max(abs(t), na.rm = TRUE) else NA_real_,
    beyond_3.5 = sum(abs(t) > 3.5, na.rm = TRUE),
    beyond_2  = sum(abs(t) > 2, na.rm = TRUE),
    expected_beyond_2 = sum(p, na.rm = TRUE),
    pass = rule,
    decision = if (!complete) "incomplete" else if (rule) "pass" else "fail",
    stringsAsFactors = FALSE
  )
}

# The decisions of every model for one kind of rows ('primary', 'secondary').
decide_all <- function(tables, models, kind, counts) {

  do.call(rbind, lapply(seq_along(models), function(i) {
    table <- tables[[i]]
    if (!is.null(table)) table <- table[table[["kind"]] == kind, , drop = FALSE]
    if (identical(kind, "secondary") && (is.null(table) || nrow(table) == 0L)) return(NULL)
    decide(table, models[[i]], kind, counts)
  }))
}

run_summary <- function(runs) {

  do.call(rbind, lapply(runs, function(run) {
    run <- complete_run(run)
    data.frame(
      model = run[["job"]][["model"]], arm = run[["job"]][["arm"]], seed = run[["job"]][["seed"]],
      wall = round(as.numeric(run[["wall"]]), 1),
      jags_seconds = round(as.numeric(run[["jags_seconds"]]), 1),
      gates = if (is.null(run[["declared"]])) NA_integer_ else run[["declared"]][["gates"]],
      blocks = if (is.null(run[["declared"]])) NA_integer_ else run[["declared"]][["blocks"]],
      samplers = if (is.null(run[["samplers"]])) "" else paste(names(run[["samplers"]]), run[["samplers"]], sep = ":", collapse = " "),
      error = if (is.null(run[["error"]])) "" else run[["error"]],
      warnings = paste(c(run[["warnings"]], run[["fit_warnings"]]), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
}


# Drivers ----------------------------------------------------------------------------

# The comparison (--part=compare or --part=rerun): runs the missing jobs with
# 'execute', reads every planned run from 'dir', and writes the run table, the
# Welch table and the decisions there.
compare_driver <- function(part, models, seeds, workers, specs, dir = out_dir, execute = run_jobs) {

  # seeds: arm 'on' 1000 + i, arm 'off' 2000 + i; the rerun uses fresh seeds
  # (5000 + i, 6000 + i)
  base <- if (identical(part, "rerun")) c(on = 5000L, off = 6000L) else c(on = 1000L, off = 2000L)
  if (identical(part, "rerun")) {
    first <- utils::read.csv(file.path(dir, "compare.csv"), stringsAsFactors = FALSE)
    kind    <- if (is.null(first[["kind"]])) rep("primary", nrow(first)) else first[["kind"]]
    failing <- first[(abs(first[["t"]]) > 2) %in% TRUE & kind == "primary", c("model", "quantity")]
    decisions <- utils::read.csv(file.path(dir, "compare-decisions.csv"), stringsAsFactors = FALSE)
    failed  <- if (is.null(decisions[["decision"]])) !decisions[["pass"]] else decisions[["decision"]] == "fail"
    models <- intersect(models, decisions[["model"]][failed %in% TRUE])
    cat("rerun models:", paste(models, collapse = ", "), "\n")
  }
  jobs <- list()
  for (model in models) for (arm in c("on", "off")) for (i in seq_len(seeds)) {
    jobs[[length(jobs) + 1L]] <- list(part = part, model = model, arm = arm, seed = base[[arm]] + i)
  }
  execute(jobs, specs, workers)
  runs <- read_runs(jobs, dir)
  counts <- run_counts(runs)
  summary <- run_summary(runs)
  utils::write.csv(summary, file.path(dir, paste0(part, "-runs.csv")), row.names = FALSE)
  cat("\nRuns with errors:", sum(nzchar(summary[["error"]])), "\n")
  if (any(nzchar(summary[["error"]]))) print(summary[nzchar(summary[["error"]]), c("model", "arm", "seed", "error")])
  if (any(!is.na(summary[["wall"]]))) {
    cat("\nRun times (wall s, median per model and arm):\n")
    print(stats::aggregate(wall ~ model + arm, summary, stats::median))
  }
  cat("\nSamplers per arm (first run):\n")
  print(unique(summary[, c("model", "arm", "gates", "blocks", "samplers")]), row.names = FALSE)

  tables <- lapply(models, function(model) {
    quantities <- if (identical(part, "rerun")) failing[["quantity"]][failing[["model"]] == model] else NULL
    welch_table(runs, model, quantities)
  })
  table <- do.call(rbind, tables)
  if (is.null(table)) {
    table <- data.frame(model = character(), kind = character(), statistic = character(), quantity = character(),
                        mean_on = numeric(), mean_off = numeric(), se = numeric(), t = numeric(), df = numeric(),
                        p_beyond_2 = numeric(), stringsAsFactors = FALSE)
  }
  utils::write.csv(table, file.path(dir, paste0(part, ".csv")), row.names = FALSE)
  # a rerun is judged by the same rule on the first run's failing quantities
  # (|t| > 2 in a failing model) only; the secondary rows are reported apart
  decisions <- decide_all(tables, models, "primary", counts)
  utils::write.csv(decisions, file.path(dir, paste0(part, "-decisions.csv")), row.names = FALSE)
  cat("\nRuns per model and arm (planned, completed, failed, missing):\n"); print(counts, row.names = FALSE)
  cat("\nDecision per model (primary: run means):\n"); print(decisions, row.names = FALSE)
  incomplete <- decisions[["model"]][decisions[["decision"]] == "incomplete"]
  if (length(incomplete) > 0L) {
    cat("INCOMPLETE: failed or missing runs, fewer than two completed runs in an arm, or no comparable quantity",
        "in model(s)", paste(incomplete, collapse = ", "), "; the comparison is not a pass.\n")
  }
  secondary <- decide_all(tables, models, "secondary", counts)
  if (!is.null(secondary)) {
    utils::write.csv(secondary, file.path(dir, paste0(part, "-decisions-secondary.csv")), row.names = FALSE)
    cat("\nSecondary (run SDs and tail frequencies, same rule, not part of the decision):\n")
    print(secondary, row.names = FALSE)
  } else {
    cat("\nSecondary rows unavailable: the runs did not record SDs and thinned draws, or an arm has fewer than two completed runs.\n")
  }
  cat("\nQuantities with |t| > 2:\n")
  print(table[(abs(table[["t"]]) > 2) %in% TRUE, , drop = FALSE], row.names = FALSE)

  invisible(decisions)
}

# One row per completed timing run; a run without an assessable quantity has
# NA effective sample sizes.
timing_rows <- function(runs, specs, dir = out_dir) {

  rows <- lapply(runs, function(run) {
    roles  <- run[["roles"]]
    if (is.null(roles)) roles <- model_roles(run[["job"]][["model"]], specs, dir)
    global <- global_parameters(names(run[["ess"]]), roles)
    ess    <- run[["ess"]][global]
    assessable <- is.finite(ess) & ess > 0
    data.frame(model = run[["job"]][["model"]], arm = run[["job"]][["arm"]], seed = run[["job"]][["seed"]],
               jags_seconds = as.numeric(run[["jags_seconds"]]),
               seconds_per_1000 = 1000 * as.numeric(run[["jags_seconds"]]) / as.numeric(run[["iterations"]]),
               assessable = sum(assessable),
               min_ess = if (any(assessable)) min(ess[assessable]) else NA_real_,
               min_ess_quantity = if (any(assessable)) names(ess)[assessable][which.min(ess[assessable])] else NA_character_,
               ess_per_second = if (any(assessable)) min(ess[assessable]) / as.numeric(run[["jags_seconds"]]) else NA_real_,
               not_assessable = paste(names(ess)[!assessable], collapse = " "),
               stringsAsFactors = FALSE)
  })
  if (length(rows) == 0L) {
    return(data.frame(model = character(), arm = character(), seed = integer(), jags_seconds = numeric(),
                      seconds_per_1000 = numeric(), assessable = integer(), min_ess = numeric(),
                      min_ess_quantity = character(), ess_per_second = numeric(), not_assessable = character(),
                      stringsAsFactors = FALSE))
  }

  return(do.call(rbind, rows))
}

# The timing summary: one row per planned model and arm with the run counts,
# means over the completed runs and ratios to 'neither'. A model is
# 'incomplete' when a run failed or is missing, an arm has no completed run,
# or a run has no assessable quantity.
timing_summary <- function(timing, counts, arms) {

  metrics <- c("seconds_per_1000", "assessable", "min_ess", "ess_per_second")
  summary <- do.call(rbind, lapply(seq_len(nrow(counts)), function(i) {
    x <- timing[timing[["model"]] == counts[["model"]][i] & timing[["arm"]] == counts[["arm"]][i], , drop = FALSE]
    values <- vapply(metrics, function(metric) if (nrow(x) > 0L) mean(x[[metric]]) else NA_real_, numeric(1))
    names  <- unlist(strsplit(x[["not_assessable"]], " ", fixed = TRUE))
    tally  <- table(names)
    cbind(counts[i, , drop = FALSE], as.data.frame(as.list(values)),
          not_assessable = if (length(tally) == 0L) "none" else
            paste0(names(tally), " (", as.integer(tally), "/", nrow(x), ")", collapse = ", "),
          no_assessable_run = any(x[["assessable"]] == 0L),
          stringsAsFactors = FALSE)
  }))
  if (is.null(summary) || nrow(summary) == 0L) {
    return(data.frame(model = character(), arm = character(), planned = integer(), completed = integer(),
                      failed = integer(), missing = integer(), status = character(),
                      seconds_per_1000 = numeric(), assessable = numeric(), min_ess = numeric(),
                      ess_per_second = numeric(), time_ratio = numeric(), speedup = numeric(),
                      not_assessable = character(), stringsAsFactors = FALSE))
  }
  neither <- summary[summary[["arm"]] == "neither", c("model", "seconds_per_1000", "ess_per_second"), drop = FALSE]
  names(neither)[2:3] <- c("seconds_neither", "ess_per_second_neither")
  summary <- merge(summary, neither, by = "model", all.x = TRUE)
  summary[["time_ratio"]] <- summary[["seconds_per_1000"]] / summary[["seconds_neither"]]
  summary[["speedup"]]    <- summary[["ess_per_second"]] / summary[["ess_per_second_neither"]]
  incomplete <- unique(summary[["model"]][summary[["failed"]] > 0L | summary[["missing"]] > 0L |
                                            summary[["completed"]] == 0L | summary[["no_assessable_run"]] %in% TRUE])
  summary[["status"]] <- ifelse(summary[["model"]] %in% incomplete, "incomplete", "complete")
  summary <- summary[order(summary[["model"]], match(summary[["arm"]], arms)),
                     c("model", "arm", "planned", "completed", "failed", "missing", "status",
                       "seconds_per_1000", "assessable", "min_ess", "ess_per_second", "time_ratio",
                       "speedup", "not_assessable")]
  rownames(summary) <- NULL

  return(summary)
}

# The timing comparison (--part=timing): runs the missing jobs with 'execute'
# (serially), reads every planned run from 'dir', and writes the per-run and
# summary tables there.
timing_driver <- function(models, seeds, specs, dir = out_dir, execute = run_jobs) {

  arms <- c("gates", "blocks", "both", "neither")
  jobs <- list()
  for (i in seq_len(seeds)) {
    # rotate the arm order by seed so that no arm always runs first
    order <- arms[(seq_along(arms) + i - 2L) %% length(arms) + 1L]
    for (model in models) for (arm in order) {
      jobs[[length(jobs) + 1L]] <- list(part = "timing", model = model, arm = arm, seed = 3000L + i)
    }
  }
  execute(jobs, specs, 1L)
  all_runs <- read_runs(jobs, dir)
  counts   <- run_counts(all_runs)
  runs     <- Filter(function(run) is.null(run[["error"]]), all_runs)
  timing   <- timing_rows(runs, specs, dir)
  utils::write.csv(timing, file.path(dir, "timing-runs.csv"), row.names = FALSE)
  summary  <- timing_summary(timing, counts, arms)
  utils::write.csv(summary, file.path(dir, "timing.csv"), row.names = FALSE)
  cat("\nTiming (means over seeds; ratios to 'neither'; minimum ESS over the assessable global parameters):\n")
  print(summary, row.names = FALSE, digits = 4)
  incomplete <- unique(summary[["model"]][summary[["status"]] == "incomplete"])
  if (length(incomplete) > 0L) {
    cat("INCOMPLETE: failed or missing runs, an arm without a completed run, or a run without an assessable",
        "quantity in model(s)", paste(incomplete, collapse = ", "), "; their timing is incomplete.\n")
  }
  cat("\nWorst-mixing assessable global quantity and the not-assessable ones per run:\n")
  print(timing[, c("model", "arm", "seed", "assessable", "min_ess_quantity", "min_ess", "jags_seconds",
                   "not_assessable")], row.names = FALSE)

  invisible(summary)
}


# Self-test ------------------------------------------------------------------------

# Synthetic run files in a temporary directory through the complete drivers
# (compare_driver(), timing_driver()); every case reads back the written
# decision or timing file and checks the decision and the run counts. Stops on
# a wrong decision, shape or count, or on any error.
selftest <- function() {

  root <- tempfile("verify-selftest-")
  no_fits <- function(...) invisible(NULL)
  # Deterministic draws: run i of an arm is a normal grid moved by offset v[i]
  # and scaled by 1 + v[i]; the 'off' arm uses the offsets in reverse order, so
  # both arms have the same set of run values (Welch t = 0) unless shifted.
  offsets <- stats::qnorm(stats::ppoints(10L)) / 5
  compare_run <- function(job, names, shift = 0, constant = FALSE) {
    i      <- job[["seed"]] %% 1000L
    offset <- if (job[["arm"]] == "on") offsets[i] else offsets[11L - i]
    draws  <- vapply(seq_along(names), function(j) {
      if (constant) rep(j, 400L) else
        stats::qnorm(stats::ppoints(400L)) * (1 + offset) + offset + shift[[min(j, length(shift))]] + j
    }, numeric(400L))
    draws  <- matrix(draws, ncol = length(names), dimnames = list(NULL, names))
    list(job = job, wall = 1, jags_seconds = 1, iterations = 100,
         means = colMeans(draws), sds = apply(draws, 2L, stats::sd),
         thinned = draws[seq(1L, nrow(draws), by = 10L), , drop = FALSE])
  }
  timing_run <- function(job, ess = 100) {
    list(job = job, wall = 1, jags_seconds = 2, iterations = 1000,
         ess = c(a = ess, b = 2 * ess), roles = c(a = "sampled", b = "sampled"))
  }
  failed_run <- function(job) list(job = job, wall = 1, error = "synthetic fit failure")
  # write the run files of one case; 'make' returns a run record or NULL (missing file)
  write_runs <- function(dir, part, arms, seeds, make) {
    dir.create(file.path(dir, "runs"), recursive = TRUE, showWarnings = FALSE)
    for (arm in arms) for (i in seq_len(seeds)) {
      seed <- switch(part, compare = c(on = 1000L, off = 2000L)[[arm]] + i,
                     rerun = c(on = 5000L, off = 6000L)[[arm]] + i, timing = 3000L + i)
      job  <- list(part = part, model = "m", arm = arm, seed = seed)
      run  <- make(job, i)
      if (!is.null(run)) saveRDS(run, run_file(job, dir))
    }
  }
  check <- function(ok, what) {
    if (!isTRUE(ok)) stop("Self-test failed: ", what, call. = FALSE)
    cat("ok:", what, "\n")
  }
  compare_case <- function(label, make, part = "compare", dir = file.path(root, label)) {
    write_runs(dir, part, c("on", "off"), 10L, make)
    output <- utils::capture.output(compare_driver(part, "m", 10L, 1L, list(m = list()), dir = dir, execute = no_fits))
    list(decisions = utils::read.csv(file.path(dir, paste0(part, "-decisions.csv")), stringsAsFactors = FALSE),
         table = utils::read.csv(file.path(dir, paste0(part, ".csv")), stringsAsFactors = FALSE),
         output = output, dir = dir)
  }
  timing_case <- function(label, make) {
    dir <- file.path(root, label)
    write_runs(dir, "timing", c("gates", "blocks", "both", "neither"), 3L, make)
    output <- utils::capture.output(timing_driver("m", 3L, list(m = list()), dir = dir, execute = no_fits))
    list(summary = utils::read.csv(file.path(dir, "timing.csv"), stringsAsFactors = FALSE), output = output)
  }
  counts_are <- function(d, on, off) {
    all(unlist(d[1L, c("planned_on", "completed_on", "failed_on", "missing_on")]) == on) &&
      all(unlist(d[1L, c("planned_off", "completed_off", "failed_off", "missing_off")]) == off)
  }

  # complete: three quantities, no difference
  x <- compare_case("complete", function(job, i) compare_run(job, c("a", "b", "c")))
  check(x$decisions$decision == "pass" && x$decisions$quantities == 3L &&
          counts_are(x$decisions, c(10, 10, 0, 0), c(10, 10, 0, 0)) &&
          sum(x$table$kind == "primary") == 3L && sum(x$table$kind == "secondary") == 6L,
        "complete: pass, 3 quantities, 10 of 10 runs per arm, 3 primary and 6 secondary rows")

  # one failed run
  x <- compare_case("one_failed", function(job, i) {
    if (job$arm == "off" && i == 10L) failed_run(job) else compare_run(job, c("a", "b"))
  })
  check(x$decisions$decision == "incomplete" && !x$decisions$pass &&
          counts_are(x$decisions, c(10, 10, 0, 0), c(10, 9, 1, 0)) && any(grepl("INCOMPLETE", x$output)),
        "one failed run: incomplete, off 10 planned / 9 completed / 1 failed")

  # one missing run file
  x <- compare_case("one_missing", function(job, i) {
    if (job$arm == "off" && i == 10L) NULL else compare_run(job, c("a", "b"))
  })
  check(x$decisions$decision == "incomplete" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 9, 0, 1)),
        "one missing run file: incomplete, off 10 planned / 9 completed / 1 missing")

  # one successful run in an arm
  x <- compare_case("one_success", function(job, i) {
    if (job$arm == "off" && i > 1L) failed_run(job) else compare_run(job, c("a", "b"))
  })
  check(x$decisions$decision == "incomplete" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 1, 9, 0)) &&
          nrow(x$table) == 2L && all(is.na(x$table$t)),
        "one successful run in an arm: incomplete, rows with NA statistics")

  # no successful run in an arm
  x <- compare_case("no_success", function(job, i) {
    if (job$arm == "off") failed_run(job) else compare_run(job, c("a", "b"))
  })
  check(x$decisions$decision == "incomplete" && counts_are(x$decisions, c(10, 10, 0, 0), c(10, 0, 10, 0)) &&
          all(is.na(x$table$t)),
        "no successful run in an arm: incomplete, rows with NA statistics")

  # no assessable (comparable) quantity: constant in every run
  x <- compare_case("no_quantity", function(job, i) compare_run(job, c("a", "b"), constant = TRUE))
  check(x$decisions$decision == "incomplete" && x$decisions$quantities == 0L && nrow(x$table) == 0L &&
          counts_are(x$decisions, c(10, 10, 0, 0), c(10, 10, 0, 0)),
        "no comparable quantity: incomplete with complete run counts")

  # a single quantity
  x <- compare_case("single", function(job, i) compare_run(job, "a"))
  check(x$decisions$decision == "pass" && x$decisions$quantities == 1L &&
          sum(x$table$kind == "primary") == 1L && sum(x$table$kind == "secondary") == 2L,
        "a single quantity: pass, 1 primary and 2 secondary rows")

  # a single-quantity rerun: the first comparison fails on 'b' only
  dir <- file.path(root, "rerun")
  x <- compare_case("rerun", function(job, i) {
    compare_run(job, c("a", "b", "c"), shift = if (job$arm == "on") c(0, 1, 0) else 0)
  }, dir = dir)
  check(x$decisions$decision == "fail" &&
          identical(x$table$quantity[x$table$kind == "primary" & abs(x$table$t) > 2], "b"),
        "a shifted quantity: fail on that quantity only")
  x <- compare_case("rerun", function(job, i) compare_run(job, c("a", "b", "c")), part = "rerun", dir = dir)
  check(x$decisions$decision == "pass" && x$decisions$quantities == 1L &&
          all(x$table$quantity == "b") && sum(x$table$kind == "secondary") == 2L &&
          counts_are(x$decisions, c(10, 10, 0, 0), c(10, 10, 0, 0)),
        "a single-quantity rerun: pass on 'b' alone")

  # timing: complete
  y <- timing_case("timing_complete", function(job, i) timing_run(job))
  check(nrow(y$summary) == 4L && all(y$summary$status == "complete") && all(y$summary$completed == 3L) &&
          all(abs(y$summary$speedup - 1) < 1e-12), "timing complete: 4 arms, 3 of 3 runs, complete")

  # timing: every run failed
  y <- timing_case("timing_all_failed", function(job, i) failed_run(job))
  check(nrow(y$summary) == 4L && all(y$summary$status == "incomplete") && all(y$summary$failed == 3L) &&
          all(y$summary$completed == 0L) && all(is.na(y$summary$ess_per_second)) &&
          any(grepl("INCOMPLETE", y$output)),
        "every timing run failed: incomplete, 3 failed per arm, NA metrics")

  # timing: no assessable quantity in any run
  y <- timing_case("timing_no_assessable", function(job, i) timing_run(job, ess = 0))
  check(nrow(y$summary) == 4L && all(y$summary$status == "incomplete") && all(y$summary$completed == 3L) &&
          all(is.na(y$summary$min_ess)) && all(y$summary$assessable == 0),
        "no assessable timing quantity: incomplete, NA effective sample sizes")

  # timing: the 'neither' arm failed, so the ratios are unavailable
  y <- timing_case("timing_neither_failed", function(job, i) {
    if (job$arm == "neither") failed_run(job) else timing_run(job)
  })
  check(nrow(y$summary) == 4L && all(y$summary$status == "incomplete") && all(is.na(y$summary$speedup)),
        "failed 'neither' arm: incomplete, no ratios")

  unlink(root, recursive = TRUE)
  cat("Self-test passed.\n")
}

if (identical(part, "selftest")) {
  selftest()
  quit(save = "no", status = 0L)
}


# Parts ----------------------------------------------------------------------------

specs  <- model_specs()
models <- grep(models_re, names(specs), value = TRUE)
cat("part:", part, " models:", paste(models, collapse = ", "), " seeds:", seeds, " workers:", workers, "\n")

if (part %in% c("compare", "rerun")) compare_driver(part, models, seeds, workers, specs)

if (identical(part, "timing")) timing_driver(models, seeds, specs)
