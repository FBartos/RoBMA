#!/usr/bin/env Rscript
# Exercise public PSOCK fitting and serialized continuation with the joint
# bias sampler and a finite tau readout. Use the configured JAGS 5 runtime
# in the main R process and in its PSOCK workers.

suppressPackageStartupMessages(library(RoBMA))

verify_joint_workers <- function() {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  data("Bem2011", package = "RoBMA", envir = environment())

  RoBMA.options(max_cores = 2L, sampler.joint = "auto", sampler.nested.weights = FALSE,
                sampler.probabilities = "tau_indicator", selection.sampler = "slice.coarse")
  fitted <- RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD",
    chains = 3L, parallel = TRUE, adapt = 500L, burnin = 100L,
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

    reports <- attr(backend, "runtime_report_values", exact = TRUE)
    stopifnot(length(reports) == 1L, length(reports[[1L]]) == workers,
              all(!vapply(reports[[1L]], is.null, logical(1))))
    assignment <- backend[["samplers"]]
    owners <- assignment[["Node"]][assignment[["Sampler"]] == "RoBMA::IndependentBiasJoint"]
    stopifnot(length(owners) == 3L, setequal(owners, expected),
              tau_gate %in% assignment[["Node"]][assignment[["Sampler"]] == "BayesTools::Finite"])

    draws <- backend[["mcmc"]]
    trace <- backend[["sampler_probabilities"]][[tau_gate]]
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
  extended <- update(readRDS(path), sample_extend = 50L, adapt = 500L,
                     parallel = TRUE, cores = 2L, silent = TRUE)
  rebuilt <- extended[["fit"]]
  check_backend(rebuilt, 150L, workers = 2L)
  stopifnot(identical(attr(rebuilt, "sampler_declarations", exact = TRUE), declarations),
            identical(attr(rebuilt, "prior_list", exact = TRUE), prior_list))
  for (chain in seq_len(3L)) {
    stopifnot(identical(as.matrix(rebuilt[["mcmc"]][[chain]])[seq_len(100L), , drop = FALSE],
                        as.matrix(backend[["mcmc"]][[chain]])),
              identical(as.matrix(rebuilt[["sampler_probabilities"]][[tau_gate]][[chain]])[seq_len(100L), , drop = FALSE],
                        as.matrix(backend[["sampler_probabilities"]][[tau_gate]][[chain]])))
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
