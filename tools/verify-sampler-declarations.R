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
# (e.g. a failure to calculate a density) is reported as such.
#
# timing: arms 'gates' (gates only), 'blocks' (blocks only), 'both', 'neither',
# --seeds runs per arm, arms interleaved within each seed, in this session.
# Seconds per 1,000 iterations (JAGS time over chains x (adapt + burnin +
# sample)) and effective sample size per second, the effective sample size
# being the minimum over the global scalar parameters (effect, moderators,
# heterogeneity and allocation parameters, weights, PET, PEESE; not unit-level
# latents or random-effect coordinates, not indicators, not auxiliaries) of
# coda::effectiveSize() summed over chains.

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}
part    <- opt("part", "compare")
out_dir <- opt("out")
workers <- as.integer(opt("workers", "6"))
seeds   <- as.integer(opt("seeds", if (identical(part, "timing")) "3" else "10"))
models_re <- opt("models", ".")
if (is.null(out_dir)) stop("Supply --out=<directory>.", call. = FALSE)
if (!part %in% c("compare", "rerun", "timing")) stop("--part must be compare, rerun or timing.", call. = FALSE)
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
  result[["ess"]]   <- coda::effectiveSize(draws)
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

# Global scalar parameters for the timing ESS (by monitored name).
global_parameters <- function(names) {

  keep <- !grepl("indicator", names, fixed = TRUE) &
    !grepl("_xRE_(Zx|COEFx|GROUP_Zx|UNIT_COEFx)", names) &
    !grepl("^prior_par_", names) &
    !grepl("^(theta|gamma|sampling_z|pi|phi)\\[", names)

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
    parallel::clusterExport(cluster, c("run_job", "quantity_means", "arm_options", "specs"),
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

read_runs <- function(jobs) {

  lapply(jobs, function(job) readRDS(run_file(job)))
}


# Welch comparison -----------------------------------------------------------------

welch_table <- function(runs, model, quantities = NULL) {

  runs <- Filter(function(run) run[["job"]][["model"]] == model && is.null(run[["error"]]), runs)
  arms <- vapply(runs, function(run) run[["job"]][["arm"]], character(1))
  names_all <- unique(unlist(lapply(runs, function(run) names(run[["means"]]))))
  # a level of an indicator absent from a run has frequency zero there
  means <- t(vapply(runs, function(run) {
    values <- run[["means"]][names_all]
    values[is.na(values) & grepl("==", names_all, fixed = TRUE)] <- 0
    values
  }, numeric(length(names_all))))
  colnames(means) <- names_all
  if (!is.null(quantities)) {
    means <- means[, intersect(quantities, names_all), drop = FALSE]
  }
  rows <- lapply(colnames(means), function(name) {
    on  <- means[arms == "on", name]
    off <- means[arms == "off", name]
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
    data.frame(model = model, quantity = name, mean_on = mean(on), mean_off = mean(off),
               se = sqrt(se2), t = t, df = df,
               p_beyond_2 = 2 * stats::pt(-2, df), stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}

decide <- function(table) {

  data.frame(
    model     = table[["model"]][1L],
    quantities = nrow(table),
    max_abs_t = max(abs(table[["t"]])),
    beyond_3.5 = sum(abs(table[["t"]]) > 3.5),
    beyond_2  = sum(abs(table[["t"]]) > 2),
    expected_beyond_2 = sum(table[["p_beyond_2"]]),
    pass = !any(abs(table[["t"]]) > 3.5) && sum(abs(table[["t"]]) > 2) <= 2 * sum(table[["p_beyond_2"]]),
    stringsAsFactors = FALSE
  )
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
    failing <- first[abs(first[["t"]]) > 2, c("model", "quantity")]
    decisions <- utils::read.csv(file.path(out_dir, "compare-decisions.csv"), stringsAsFactors = FALSE)
    models <- intersect(models, decisions[["model"]][!decisions[["pass"]]])
    cat("rerun models:", paste(models, collapse = ", "), "\n")
  }
  jobs <- list()
  for (model in models) for (arm in c("on", "off")) for (i in seq_len(seeds)) {
    jobs[[length(jobs) + 1L]] <- list(part = part, model = model, arm = arm, seed = base[[arm]] + i)
  }
  run_jobs(jobs, specs, workers)
  runs <- read_runs(jobs)
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
  # (|t| > 2 in a failing model) only
  decisions <- do.call(rbind, lapply(Filter(Negate(is.null), tables), decide))
  utils::write.csv(decisions, file.path(out_dir, paste0(part, "-decisions.csv")), row.names = FALSE)
  cat("\nDecision per model:\n"); print(decisions, row.names = FALSE)
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
  runs <- Filter(function(run) is.null(run[["error"]]), read_runs(jobs))
  timing <- do.call(rbind, lapply(runs, function(run) {
    global <- global_parameters(names(run[["ess"]]))
    ess    <- run[["ess"]][global]
    ess    <- ess[is.finite(ess) & ess > 0]
    data.frame(model = run[["job"]][["model"]], arm = run[["job"]][["arm"]], seed = run[["job"]][["seed"]],
               jags_seconds = run[["jags_seconds"]],
               seconds_per_1000 = 1000 * run[["jags_seconds"]] / run[["iterations"]],
               min_ess = min(ess), min_ess_quantity = names(ess)[which.min(ess)],
               ess_per_second = min(ess) / run[["jags_seconds"]], stringsAsFactors = FALSE)
  }))
  utils::write.csv(timing, file.path(out_dir, "timing-runs.csv"), row.names = FALSE)
  summary <- stats::aggregate(cbind(seconds_per_1000, min_ess, ess_per_second) ~ model + arm, timing, mean)
  neither <- summary[summary[["arm"]] == "neither", c("model", "seconds_per_1000", "ess_per_second")]
  names(neither)[2:3] <- c("seconds_neither", "ess_per_second_neither")
  summary <- merge(summary, neither, by = "model")
  summary[["time_ratio"]] <- summary[["seconds_per_1000"]] / summary[["seconds_neither"]]
  summary[["speedup"]]    <- summary[["ess_per_second"]] / summary[["ess_per_second_neither"]]
  summary <- summary[order(summary[["model"]], match(summary[["arm"]], arms)),
                     c("model", "arm", "seconds_per_1000", "min_ess", "ess_per_second", "time_ratio", "speedup")]
  utils::write.csv(summary, file.path(out_dir, "timing.csv"), row.names = FALSE)
  cat("\nTiming (means over seeds; ratios to 'neither'):\n")
  print(summary, row.names = FALSE, digits = 4)
  cat("\nWorst-mixing global quantity per run:\n")
  print(timing[, c("model", "arm", "seed", "min_ess_quantity", "min_ess", "jags_seconds")], row.names = FALSE)
}
