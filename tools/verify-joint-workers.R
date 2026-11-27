#!/usr/bin/env Rscript
# Exercise public PSOCK fitting and serialized continuation with the joint
# bias sampler and a finite tau readout, plus the automatic Gaussian row
# budget. Install both packages into the worker library before starting.

suppressPackageStartupMessages(library(RoBMA))

packages <- c("RoBMA", "BayesTools", "rjags")
versions <- vapply(packages, function(package) {
  path <- getNamespaceInfo(asNamespace(package), "path")
  if (!file.exists(file.path(path, "Meta", "package.rds"))) {
    stop("Install ", package, " into the worker library before running this verification.", call. = FALSE)
  }
  as.character(utils::packageVersion(package, lib.loc = dirname(path)))
}, character(1))
stopifnot(rjags::jags.version() >= "5", utils::packageVersion("rjags") >= "5.5")
cat("Installed package versions:", paste(paste(packages, versions), collapse = ", "), "\n")

verify_joint_workers <- function() {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  data("Bem2011", package = "RoBMA", envir = environment())

  requested_tau <- " tau_indicator [ 1 ] "
  RoBMA.options(max_cores = 2L, sampler.joint = "auto", sampler.nested.weights = FALSE,
                sampler.probabilities = requested_tau, selection.sampler = "slice.coarse")
  fitted <- RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD",
    chains = 3L, parallel = TRUE, adapt = 250L, burnin = 50L,
    sample = 100L, thin = 2L, autofit = FALSE, silent = TRUE, seed = 806L,
    convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))
  backend <- fitted[["fit"]]
  stopifnot(isTRUE(backend[["has_posterior"]]), is.null(backend[["errors"]]))

  prior_list <- attr(backend, "prior_list", exact = TRUE)
  declarations <- attr(backend, "sampler_declarations", exact = TRUE)
  stopifnot(identical(prior_list, fitted[["priors"]][["outcome"]]))
  bias <- BayesTools::JAGS_mixture_sources(prior_list, "bias")
  mu   <- BayesTools::JAGS_mixture_sources(prior_list, "mu")
  tau_gate <- BayesTools::JAGS_mixture_sources(prior_list, "tau")[["gate"]]
  active_mu <- Filter(function(branch) isTRUE(branch[["stochastic"]]), mu[["branches"]])
  expected <- c(bias[["gate"]], mu[["gate"]], active_mu[[1L]][["source"]])

  check_backend <- function(backend, samples, workers) {

    declarations <- attr(backend, "sampler_declarations", exact = TRUE)
    recorded_tau <- declarations[["probabilities"]][["node"]]
    audit <- declarations[["robma"]][["probabilities"]]
    stopifnot(identical(recorded_tau, paste0(tau_gate, "[1]")),
              identical(audit[["node"]], recorded_tau), identical(audit[["requested_node"]], requested_tau),
              identical(audit[["sampler"]], "BayesTools::Finite"), identical(audit[["status"]], "active"))
    reports <- attr(backend, "runtime_report_values", exact = TRUE)
    stopifnot(length(reports) == 1L, length(reports[[1L]]) == workers,
              all(!vapply(reports[[1L]], is.null, logical(1))))
    assignment <- backend[["samplers"]]
    owners <- assignment[["Node"]][assignment[["Sampler"]] == "RoBMA::IndependentBiasJoint"]
    stopifnot(length(owners) == 3L, setequal(owners, expected))

    draws <- backend[["mcmc"]]
    trace <- backend[["sampler_probabilities"]][[recorded_tau]]
    stopifnot(length(draws) == 3L, length(trace) == 3L)
    for (chain in seq_along(trace)) {
      probability <- trace[[chain]][, c("p[1]", "p[2]"), drop = FALSE]
      stopifnot(nrow(probability) == samples, all(is.finite(probability)),
                all(probability >= 0 & probability <= 1),
                max(abs(rowSums(probability) - 1)) < 1e-12,
                identical(stats::time(trace[[chain]]), stats::time(draws[[chain]])),
                identical(as.numeric(trace[[chain]][, "selected"]),
                          as.numeric(draws[[chain]][, tau_gate])))
    }
  }

  check_backend(backend, 100L, workers = 3L)
  saveRDS(fitted, path)
  RoBMA.options(sampler.joint = "off", sampler.probabilities = character(),
                sampler.nested.weights = TRUE, selection.sampler = "jags.default")
  extended <- update(readRDS(path), sample_extend = 50L, adapt = 250L,
                     parallel = TRUE, cores = 2L, silent = TRUE)
  rebuilt <- extended[["fit"]]
  check_backend(rebuilt, 150L, workers = 2L)
  stopifnot(identical(attr(rebuilt, "sampler_declarations", exact = TRUE), declarations),
            identical(attr(rebuilt, "prior_list", exact = TRUE), prior_list))
  recorded_tau <- declarations[["probabilities"]][["node"]]
  for (chain in seq_len(3L)) {
    stopifnot(identical(as.matrix(rebuilt[["mcmc"]][[chain]])[seq_len(100L), , drop = FALSE],
                        as.matrix(backend[["mcmc"]][[chain]])),
              identical(as.matrix(rebuilt[["sampler_probabilities"]][[recorded_tau]][[chain]])[seq_len(100L), , drop = FALSE],
                        as.matrix(backend[["sampler_probabilities"]][[recorded_tau]][[chain]])))
  }

  for (fit in list(fitted, extended)) {
    estimates  <- summary(fit, include_mcmc_diagnostics = FALSE)[["estimates"]]
    predictions <- as.matrix(predict(fit, type = "terms"))
    stopifnot(nrow(estimates) > 0L, all(is.finite(estimates[, "Mean"])),
              length(predictions) > 0L, all(is.finite(predictions)))
  }
  cat("Three PSOCK chains used three workers at fit and two after serialized continuation.\n")
  cat("Joint roots:", paste(expected, collapse = ", "), "\n")
  cat("Finite tau readouts align with all chains. Prior lists, posterior prefixes,",
      "sampler declarations, summaries, and predictions passed.\n")
}

verify_joint_workers()

verify_gaussian_workers <- function() {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  RoBMA.options(max_cores = 2L, sampler.probabilities = character(), selection.sampler = "slice.coarse")

  rows <- seq_len(801L)
  dat <- data.frame(y = .14 + .05 * sin(rows) + .001 * rows, se = .2, x = rows)
  fit <- function(policy) {
    RoBMA.options(sampler.joint = policy)
    # JAGS_fit checks the package builds on the two actual workers before fitting.
    model <- BMA(yi = y, sei = se, mods = ~ x, data = dat, measure = "SMD",
      chains = 2L, parallel = TRUE, adapt = 500L, burnin = 50L, sample = 100L, seed = 165L,
      convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))
    stopifnot(model[["fit_control"]][["cores"]] == 2L)
    model[["fit"]]
  }

  bypassed <- fit("auto")
  ordinary <- fit("off")
  declarations <- attr(bypassed, "sampler_declarations", exact = TRUE)
  report <- declarations[["robma"]][["joint"]]
  stopifnot(identical(report[["status"]], "budget_bypassed"), is.null(declarations[["joint"]]),
            length(bypassed[["mcmc"]]) == 2L, length(bypassed[["end_state"]]) == 2L,
            nrow(bypassed[["samplers"]]) > 0L,
            is.null(bypassed[["samplers_by_worker"]]), is.null(ordinary[["samplers_by_worker"]]))
  for (field in c("model", "data", "samplers", "mcmc", "end_state")) {
    stopifnot(identical(bypassed[[field]], ordinary[[field]]))
  }
  cat("At 801 estimates two PSOCK workers compiled the same model and data with auto and off.",
      "Sampler assignments, draws, and final chain states were identical.\n")
}

verify_gaussian_workers()
