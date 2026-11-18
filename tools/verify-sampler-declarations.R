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
# missing); a comparison or timing with a failed or missing run is reported as
# 'incomplete', not as a pass. --part=selftest checks the comparison code on
# synthetic run records (no fits) and stops on a wrong shape or decision.
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
  result <- list(job = job, wall = wall, warnings = unique(warnings))
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
model_roles <- function(model, specs) {

  file <- file.path(out_dir, paste0("roles-", model, ".rds"))
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

run_file <- function(job) {

  file.path(out_dir, "runs", sprintf("%s-%s-%s-%03d.rds", job[["part"]], job[["model"]], job[["arm"]], job[["seed"]]))
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

# The saved runs of the jobs; a job without a saved run is marked missing.
read_runs <- function(jobs) {

  lapply(jobs, function(job) {
    file <- run_file(job)
    if (!file.exists(file)) return(list(job = job, error = "missing run file", missing = TRUE))
    readRDS(file)
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

  # secondary: run SDs and tail frequencies, when the runs recorded them
  if (length(runs) > 0L &&
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

# Welch statistics of the per-run values of each column (runs in rows).
welch_rows <- function(values, arms, model, statistic, kind) {

  rows <- lapply(colnames(values), function(name) {
    on  <- values[arms == "on", name]
    off <- values[arms == "off", name]
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
    data.frame(model = model, kind = kind, statistic = statistic, quantity = name,
               mean_on = mean(on), mean_off = mean(off),
               se = sqrt(se2), t = t, df = df,
               p_beyond_2 = 2 * stats::pt(-2, df), stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}

# The decision of one model and kind: 'incomplete' when any planned run of the
# model failed or is missing, otherwise 'pass' or 'fail' by the D4 rule.
decide <- function(table, model, kind, counts) {

  counts <- counts[counts[["model"]] == model, , drop = FALSE]
  arm_count <- function(arm, field) {
    value <- counts[[field]][counts[["arm"]] == arm]
    if (length(value) == 0L) 0L else as.integer(value)
  }
  t <- if (is.null(table)) numeric() else table[["t"]]
  p <- if (is.null(table)) numeric() else table[["p_beyond_2"]]
  complete <- sum(counts[["failed"]]) == 0L && sum(counts[["missing"]]) == 0L &&
    all(vapply(c("on", "off"), arm_count, integer(1), field = "planned") > 0L)
  rule <- length(t) > 0L && !any(abs(t) > 3.5) && sum(abs(t) > 2) <= 2 * sum(p)
  data.frame(
    model     = model,
    kind      = kind,
    planned_on = arm_count("on", "planned"), completed_on = arm_count("on", "completed"),
    failed_on = arm_count("on", "failed"), missing_on = arm_count("on", "missing"),
    planned_off = arm_count("off", "planned"), completed_off = arm_count("off", "completed"),
    failed_off = arm_count("off", "failed"), missing_off = arm_count("off", "missing"),
    quantities = length(t),
    max_abs_t = if (length(t) > 0L) max(abs(t)) else NA_real_,
    beyond_3.5 = sum(abs(t) > 3.5),
    beyond_2  = sum(abs(t) > 2),
    expected_beyond_2 = sum(p),
    pass = complete && rule,
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
    data.frame(
      model = run[["job"]][["model"]], arm = run[["job"]][["arm"]], seed = run[["job"]][["seed"]],
      wall = round(run[["wall"]], 1),
      jags_seconds = if (is.null(run[["jags_seconds"]])) NA_real_ else round(run[["jags_seconds"]], 1),
      gates = if (is.null(run[["declared"]])) NA_integer_ else run[["declared"]][["gates"]],
      blocks = if (is.null(run[["declared"]])) NA_integer_ else run[["declared"]][["blocks"]],
      samplers = if (is.null(run[["samplers"]])) "" else paste(names(run[["samplers"]]), run[["samplers"]], sep = ":", collapse = " "),
      error = if (is.null(run[["error"]])) "" else run[["error"]],
      warnings = paste(c(run[["warnings"]], run[["fit_warnings"]]), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
}


# Self-test ------------------------------------------------------------------------

# Synthetic run records through the comparison code: one quantity, several
# quantities with secondary rows, a rerun restricted to one quantity, a failed
# run in an arm, a clean case and a shifted case. Stops on a wrong shape or
# decision.
selftest <- function() {

  # Deterministic draws: run i of an arm is a normal grid moved by offset v[i]
  # and scaled by 1 + v[i]; the 'off' arm uses the offsets in reverse order, so
  # both arms have the same set of run values (Welch t = 0) unless shifted.
  offsets <- stats::qnorm(stats::ppoints(10L)) / 5
  make_run <- function(arm, seed, names, shift = 0, secondary = TRUE, error = NULL) {
    job <- list(part = "selftest", model = "m", arm = arm, seed = seed)
    if (!is.null(error)) return(list(job = job, error = error))
    i      <- seed %% 1000L
    offset <- if (arm == "on") offsets[i] else offsets[11L - i]
    draws  <- vapply(seq_along(names), function(j) {
      stats::qnorm(stats::ppoints(400L)) * (1 + offset) + offset + shift + j
    }, numeric(400L))
    draws  <- matrix(draws, ncol = length(names), dimnames = list(NULL, names))
    run <- list(job = job, means = colMeans(draws))
    if (secondary) {
      run[["sds"]] <- apply(draws, 2L, stats::sd)
      run[["thinned"]] <- draws[seq(1L, nrow(draws), by = 10L), , drop = FALSE]
    }
    run
  }
  make_runs <- function(names, shift = 0, secondary = TRUE, fail_off = FALSE) {
    runs <- c(lapply(1:10, function(i) make_run("on", 1000L + i, names, shift, secondary)),
              lapply(1:10, function(i) make_run("off", 2000L + i, names, 0, secondary)))
    if (fail_off) runs[[20L]] <- make_run("off", 2010L, names, error = "synthetic failure")
    runs
  }
  check <- function(ok, what) {
    if (!isTRUE(ok)) stop("Self-test failed: ", what, call. = FALSE)
    cat("ok:", what, "\n")
  }
  primary <- function(table) table[table[["kind"]] == "primary", , drop = FALSE]
  count_rows <- function(table, kind) if (is.null(table)) 0L else sum(table[["kind"]] == kind)

  # one quantity, primary and secondary
  runs  <- make_runs("a")
  table <- welch_table(runs, "m")
  check(count_rows(table, "primary") == 1L && count_rows(table, "secondary") == 2L,
        "one quantity gives one primary and two secondary rows")
  check(decide(primary(table), "m", "primary", run_counts(runs))[["decision"]] == "pass",
        "one quantity, no difference: pass")

  # several quantities
  runs  <- make_runs(c("a", "b", "c"))
  table <- welch_table(runs, "m")
  check(count_rows(table, "primary") == 3L && count_rows(table, "secondary") == 6L,
        "three quantities give three primary and six secondary rows")
  decision <- decide(primary(table), "m", "primary", run_counts(runs))
  check(decision[["decision"]] == "pass" && decision[["completed_on"]] == 10L && decision[["completed_off"]] == 10L,
        "clean case: pass with 10 of 10 runs per arm")

  # a rerun restricted to one quantity of several
  table <- welch_table(runs, "m", quantities = "b")
  check(count_rows(table, "primary") == 1L && count_rows(table, "secondary") == 2L &&
          all(table[["quantity"]] == "b"), "a rerun of one quantity keeps the matrix shape")
  table <- welch_table(runs, "m", quantities = "absent")
  check(count_rows(table, "primary") == 0L && count_rows(table, "secondary") == 0L,
        "a rerun of no available quantity gives no rows")
  check(decide(NULL, "m", "primary", run_counts(runs))[["decision"]] == "fail",
        "no quantities is not a pass")

  # runs without secondary records
  runs  <- make_runs(c("a", "b"), secondary = FALSE)
  table <- welch_table(runs, "m")
  check(count_rows(table, "primary") == 2L && count_rows(table, "secondary") == 0L,
        "runs without SDs and thinned draws give primary rows only")

  # a failed run in the 'off' arm
  runs  <- make_runs(c("a", "b"), fail_off = TRUE)
  table <- welch_table(runs, "m")
  decision <- decide(primary(table), "m", "primary", run_counts(runs))
  check(decision[["decision"]] == "incomplete" && !decision[["pass"]] &&
          decision[["failed_off"]] == 1L && decision[["completed_off"]] == 9L,
        "one failed run: incomplete with the counts")
  missing <- runs
  missing[[20L]] <- list(job = missing[[20L]][["job"]], error = "missing run file", missing = TRUE)
  decision <- decide(primary(table), "m", "primary", run_counts(missing))
  check(decision[["decision"]] == "incomplete" && decision[["missing_off"]] == 1L,
        "one missing run: incomplete with the counts")

  # a shifted 'on' arm fails
  runs  <- make_runs(c("a", "b"), shift = 1)
  table <- welch_table(runs, "m")
  check(decide(primary(table), "m", "primary", run_counts(runs))[["decision"]] == "fail",
        "a shifted arm: fail")

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

if (part %in% c("compare", "rerun")) {

  # seeds: arm 'on' 1000 + i, arm 'off' 2000 + i; the rerun uses fresh seeds
  # (5000 + i, 6000 + i)
  base <- if (identical(part, "rerun")) c(on = 5000L, off = 6000L) else c(on = 1000L, off = 2000L)
  if (identical(part, "rerun")) {
    first <- utils::read.csv(file.path(out_dir, "compare.csv"), stringsAsFactors = FALSE)
    kind    <- if (is.null(first[["kind"]])) rep("primary", nrow(first)) else first[["kind"]]
    failing <- first[abs(first[["t"]]) > 2 & kind == "primary", c("model", "quantity")]
    decisions <- utils::read.csv(file.path(out_dir, "compare-decisions.csv"), stringsAsFactors = FALSE)
    failed  <- if (is.null(decisions[["decision"]])) !decisions[["pass"]] else decisions[["decision"]] == "fail"
    models <- intersect(models, decisions[["model"]][failed])
    cat("rerun models:", paste(models, collapse = ", "), "\n")
  }
  jobs <- list()
  for (model in models) for (arm in c("on", "off")) for (i in seq_len(seeds)) {
    jobs[[length(jobs) + 1L]] <- list(part = part, model = model, arm = arm, seed = base[[arm]] + i)
  }
  run_jobs(jobs, specs, workers)
  runs <- read_runs(jobs)
  counts <- run_counts(runs)
  summary <- run_summary(runs)
  utils::write.csv(summary, file.path(out_dir, paste0(part, "-runs.csv")), row.names = FALSE)
  cat("\nRuns with errors:", sum(nzchar(summary[["error"]])), "\n")
  if (any(nzchar(summary[["error"]]))) print(summary[nzchar(summary[["error"]]), c("model", "arm", "seed", "error")])
  cat("\nRun times (wall s, median per model and arm):\n")
  print(stats::aggregate(wall ~ model + arm, summary, stats::median))
  cat("\nSamplers per arm (first run):\n")
  print(unique(summary[, c("model", "arm", "gates", "blocks", "samplers")]), row.names = FALSE)

  tables <- lapply(models, function(model) {
    quantities <- if (identical(part, "rerun")) failing[["quantity"]][failing[["model"]] == model] else NULL
    welch_table(runs, model, quantities)
  })
  table <- do.call(rbind, tables)
  utils::write.csv(table, file.path(out_dir, paste0(part, ".csv")), row.names = FALSE)
  # a rerun is judged by the same rule on the first run's failing quantities
  # (|t| > 2 in a failing model) only; the secondary rows are reported apart
  decisions <- decide_all(tables, models, "primary", counts)
  utils::write.csv(decisions, file.path(out_dir, paste0(part, "-decisions.csv")), row.names = FALSE)
  cat("\nRuns per model and arm (planned, completed, failed, missing):\n"); print(counts, row.names = FALSE)
  cat("\nDecision per model (primary: run means):\n"); print(decisions, row.names = FALSE)
  if (any(decisions[["decision"]] == "incomplete")) {
    cat("INCOMPLETE: failed or missing runs in model(s)",
        paste(decisions[["model"]][decisions[["decision"]] == "incomplete"], collapse = ", "),
        "; the comparison is not a pass.\n")
  }
  secondary <- decide_all(tables, models, "secondary", counts)
  if (!is.null(secondary)) {
    utils::write.csv(secondary, file.path(out_dir, paste0(part, "-decisions-secondary.csv")), row.names = FALSE)
    cat("\nSecondary (run SDs and tail frequencies, same rule, not part of the decision):\n")
    print(secondary, row.names = FALSE)
  } else {
    cat("\nSecondary rows unavailable: the runs did not record SDs and thinned draws.\n")
  }
  cat("\nQuantities with |t| > 2:\n")
  print(table[abs(table[["t"]]) > 2, ], row.names = FALSE)
}

if (identical(part, "timing")) {

  arms <- c("gates", "blocks", "both", "neither")
  jobs <- list()
  for (i in seq_len(seeds)) {
    # rotate the arm order by seed so that no arm always runs first
    order <- arms[(seq_along(arms) + i - 2L) %% length(arms) + 1L]
    for (model in models) for (arm in order) {
      jobs[[length(jobs) + 1L]] <- list(part = part, model = model, arm = arm, seed = 3000L + i)
    }
  }
  run_jobs(jobs, specs, 1L)
  all_runs <- read_runs(jobs)
  counts   <- run_counts(all_runs)
  runs     <- Filter(function(run) is.null(run[["error"]]), all_runs)
  timing <- do.call(rbind, lapply(runs, function(run) {
    roles  <- run[["roles"]]
    if (is.null(roles)) roles <- model_roles(run[["job"]][["model"]], specs)
    global <- global_parameters(names(run[["ess"]]), roles)
    ess    <- run[["ess"]][global]
    assessable <- is.finite(ess) & ess > 0
    data.frame(model = run[["job"]][["model"]], arm = run[["job"]][["arm"]], seed = run[["job"]][["seed"]],
               jags_seconds = run[["jags_seconds"]],
               seconds_per_1000 = 1000 * run[["jags_seconds"]] / run[["iterations"]],
               assessable = sum(assessable),
               min_ess = min(ess[assessable]), min_ess_quantity = names(ess)[assessable][which.min(ess[assessable])],
               ess_per_second = min(ess[assessable]) / run[["jags_seconds"]],
               not_assessable = paste(names(ess)[!assessable], collapse = " "),
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(timing, file.path(out_dir, "timing-runs.csv"), row.names = FALSE)
  # the free parameters without an effective sample size, with the share of the
  # model's runs in that arm where it happened
  not_assessable <- do.call(rbind, lapply(split(timing, list(timing[["model"]], timing[["arm"]]), drop = TRUE), function(x) {
    names <- unlist(strsplit(x[["not_assessable"]], " ", fixed = TRUE))
    counts <- table(names)
    data.frame(model = x[["model"]][1L], arm = x[["arm"]][1L],
               not_assessable = if (length(counts) == 0L) "none" else
                 paste0(names(counts), " (", as.integer(counts), "/", nrow(x), ")", collapse = ", "),
               stringsAsFactors = FALSE)
  }))
  summary <- stats::aggregate(cbind(seconds_per_1000, assessable, min_ess, ess_per_second) ~ model + arm, timing, mean)
  neither <- summary[summary[["arm"]] == "neither", c("model", "seconds_per_1000", "ess_per_second")]
  names(neither)[2:3] <- c("seconds_neither", "ess_per_second_neither")
  summary <- merge(summary, neither, by = "model")
  summary[["time_ratio"]] <- summary[["seconds_per_1000"]] / summary[["seconds_neither"]]
  summary[["speedup"]]    <- summary[["ess_per_second"]] / summary[["ess_per_second_neither"]]
  summary <- merge(summary, not_assessable, by = c("model", "arm"))
  # every planned run of a model and arm, and whether the model's timing is
  # complete (its ratios use the 'neither' arm)
  summary <- merge(counts, summary, by = c("model", "arm"), all.x = TRUE)
  incomplete_models <- unique(summary[["model"]][summary[["failed"]] > 0L | summary[["missing"]] > 0L])
  summary[["status"]] <- ifelse(summary[["model"]] %in% incomplete_models, "incomplete", "complete")
  summary <- summary[order(summary[["model"]], match(summary[["arm"]], arms)),
                     c("model", "arm", "planned", "completed", "failed", "missing", "status",
                       "seconds_per_1000", "assessable", "min_ess", "ess_per_second", "time_ratio",
                       "speedup", "not_assessable")]
  utils::write.csv(summary, file.path(out_dir, "timing.csv"), row.names = FALSE)
  cat("\nTiming (means over seeds; ratios to 'neither'; minimum ESS over the assessable global parameters):\n")
  print(summary, row.names = FALSE, digits = 4)
  if (length(incomplete_models) > 0L) {
    cat("INCOMPLETE: failed or missing runs in model(s)", paste(incomplete_models, collapse = ", "),
        "; their timing is incomplete.\n")
  }
  cat("\nWorst-mixing assessable global quantity and the not-assessable ones per run:\n")
  print(timing[, c("model", "arm", "seed", "assessable", "min_ess_quantity", "min_ess", "jags_seconds",
                   "not_assessable")], row.names = FALSE)
}
