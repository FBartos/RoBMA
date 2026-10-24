context("Exactness of the corrected sampler on cluster and factor selection blocks")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The corrected sampler draws a slice transition on the surrogate (the exact
# numerator minus the log normalizer of rule 3 of the block's own schedule) and
# accepts it with the exact-to-surrogate ratio, so its invariant law is the exact
# posterior whatever the quality of the surrogate. The models are a location
# mu ~ t(4) with y = 0 observed on the rows of one selection block; their exact
# posterior moments follow from QUADPACK (helper-selection-coarse-exactness.R).
# The schedules are crude, so that rule 3 has three nodes (per axis) while the
# exact kernel accepts a later rung; an uncorrected chain on the surrogate would
# sample another posterior, which each test asserts to be many Monte Carlo
# standard errors away from the exact one.

.expect_corrected_chain_is_exact <- function(kind, draws = 10000L) {

  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  setup <- .coarse_exactness_setup(kind)

  # The exact kernel accepts a rule later than the surrogate's third at states
  # across the posterior bulk, at a visibly different normalizer.
  rungs <- .coarse_exactness_rungs(setup, c(-2, -1.3, -.6, .3))
  expect_false(anyNA(rungs[["accepted_rule"]]))
  expect_true(all(rungs[["accepted_rule"]] > 3L))
  expect_true(all(rungs[["gap"]] > 1e-2))
  expect_true(all(rungs[["surrogate_vs_rule"]] < 1e-12))

  exact_tracker <- new.env()
  exact_tracker$relative_error <- 0
  exact <- .coarse_exactness_moments(setup, function(mu) {
    .coarse_exactness_quadpack_log_normalizer(setup, mu, exact_tracker)
  }, exact_tracker)
  surrogate_tracker <- new.env()
  surrogate_tracker$relative_error <- 0
  surrogate_only <- .coarse_exactness_moments(setup, function(mu) {
    .coarse_exactness_rule_log_normalizer(setup, mu, setup$orders[3L])
  }, surrogate_tracker)

  chain <- .coarse_exactness_chain(setup, draws)
  expect_identical(chain$samplers, "RoBMA::CoarseCorrectedSlice")
  counters <- chain$counters
  expect_identical(counters[["proposed"]], as.numeric(draws))
  expect_identical(counters[["accepted"]] + counters[["correction_rejected"]], as.numeric(draws))
  expect_gte(counters[["correction_rejected"]], 0.05 * draws)

  values      <- chain$draws
  observables <- cbind(mean = values, variance = (values - mean(values))^2)
  ess         <- as.numeric(coda::effectiveSize(coda::mcmc(observables)))
  mcse        <- sqrt(apply(observables, 2L, stats::var) / ess)
  estimate    <- c(mean(values), stats::var(values))
  reference   <- c(exact$mean, exact$variance)
  reference_error <- c(exact$mean_error, exact$variance_error)
  uncorrected <- c(surrogate_only$mean, surrogate_only$variance)
  expect_true(all(is.finite(mcse)) && all(mcse > 0))
  expect_true(all(ess >= 200))
  # the surrogate-only target is another distribution: an uncorrected chain fails
  expect_true(all(abs(reference - uncorrected) >= 8 * mcse))
  # the corrected chain is the exact one
  expect_true(all(abs(estimate - reference) <= 6 * mcse + 2 * reference_error))
  invisible(list(estimate = estimate, reference = reference, uncorrected = uncorrected,
                 mcse = mcse, ess = ess, counters = counters, rungs = rungs))
}


test_that("corrected chains on a cluster and a rank-two factor block run and correct (smoke)", {

  skip_if_not_installed("rjags")
  for (kind in c("cluster", "factor")) {
    setup  <- .coarse_exactness_setup(kind)
    chain  <- .coarse_exactness_chain(setup, 500L)
    counts <- chain$counters
    expect_identical(chain$samplers, "RoBMA::CoarseCorrectedSlice", info = kind)
    expect_identical(counts[["proposed"]], 500, info = kind)
    expect_identical(counts[["accepted"]] + counts[["correction_rejected"]], 500, info = kind)
    # the crude third rule is wrong, so the correction has something to reject
    expect_gt(counts[["correction_rejected"]], 0, label = kind)
    expect_true(all(is.finite(chain$draws)), info = kind)
  }
})


test_that("the corrected chain on a rank-one cluster block targets the exact posterior", {

  result <- .expect_corrected_chain_is_exact("cluster")
})


test_that("the corrected chain on a rank-two factor block targets the exact posterior", {

  skip_if_not_certification("The nested QUADPACK reference runs in the numerical-kernels case.")
  result <- .expect_corrected_chain_is_exact("factor")
})
