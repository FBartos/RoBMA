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


# The residual SD and the loadings of a block can be sampled (models that
# integrate random effects): the parents that change during sampling are then
# the residual SD of a cluster, or one loading column of a factor block, and the
# sampler must read them at every update. The location is fixed and the sampled
# parameter has the posterior
#   p(s | y) ~ dunif(s) N(y | mu, S(s)) / A(s),
# whose moments follow from QUADPACK over s, with A by QUADPACK (cluster) or by
# a Gauss-Hermite rule of 40 nodes per axis (factor), whose distance from the rule
# of 64 nodes bounds its error.
.expect_sampled_parameter_chain_is_exact <- function(kind, draws = 8000L) {

  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  setup    <- .coarse_exactness_setup(kind)
  location <- -2
  y        <- location + c(-.3, .2, -.8, .1)
  pattern  <- c(0, 0, 1, 1)
  if (kind == "cluster") {
    range <- c(.3, 1.5)
    varied <- function(value) { changed <- setup; changed$residual_sd <- rep(value, 4L); changed }
    numerator <- function(value) {
      .coarse_gaussian_log_density(y, location, rep(value, 4L), setup$loading)
    }
    hand <- .coarse_hand_model("cluster", 4L, 1L, setup$rules,
      data = list(y = y, means = rep(location, 4L), loading = as.numeric(setup$loading)),
      definitions = "  sdr ~ dunif(0.3, 1.5)\n  for (j in 1:4) { residual_sd[j] <- sdr }",
      inits = list(sdr = .8))
    monitor <- "sdr"
  } else {
    range <- c(.05, 1)
    varied <- function(value) {
      changed <- setup
      changed$loading <- cbind(setup$loading[, 1L], value * pattern)
      changed
    }
    numerator <- function(value) {
      .coarse_gaussian_log_density(y, location, setup$residual_sd, varied(value)$loading)
    }
    hand <- .coarse_hand_model("factor", 4L, 2L, setup$rules,
      data = list(y = y, means = rep(location, 4L), residual_sd = setup$residual_sd,
                  load1 = setup$loading[, 1L], pattern = pattern),
      definitions = paste0("  lam ~ dunif(0.05, 1)\n  for (j in 1:4) {\n",
                           "    loading[j, 1] <- load1[j]\n    loading[j, 2] <- lam * pattern[j]\n  }"),
      inits = list(lam = .5))
    monitor <- "lam"
  }
  # the sampled parameter's node gets the corrected sampler
  expect_identical(hand$samplers, "RoBMA::CoarseCorrectedSlice", info = kind)

  exact_tracker <- new.env()
  exact_tracker$relative_error <- 0
  rule_log_a <- function(value, order) {
    .coarse_exactness_rule_log_normalizer(varied(value), location, order)
  }
  exact_log_a <- if (kind == "cluster") {
    function(value) .coarse_exactness_quadpack_log_normalizer(varied(value), location, exact_tracker)
  } else {
    function(value) rule_log_a(value, 40L)
  }
  if (kind == "factor") {
    # the largest distance of the 40-node rule from the 64-node rule, across the range, is
    # the error allowed for the reference (far below the Monte Carlo error)
    distance <- vapply(seq(.05, 1, length.out = 12L), function(value) {
      abs(rule_log_a(value, 40L) - rule_log_a(value, 64L))
    }, numeric(1L))
    expect_lte(max(distance), 1e-5)
    exact_tracker$relative_error <- 2 * max(distance)
  }
  exact <- .coarse_moments_1d(function(value) numerator(value) - exact_log_a(value),
                              range[1L], range[2L], exact_tracker)
  surrogate_tracker <- new.env()
  surrogate_tracker$relative_error <- 0
  surrogate_only <- .coarse_moments_1d(function(value) numerator(value) - rule_log_a(value, setup$orders[3L]),
                                       range[1L], range[2L], surrogate_tracker)

  chain <- .coarse_hand_chain(hand, monitor, draws)
  counters <- chain$counters
  expect_identical(counters[["proposed"]], as.numeric(draws), info = kind)
  expect_identical(counters[["accepted"]] + counters[["correction_rejected"]], as.numeric(draws), info = kind)
  expect_gt(counters[["correction_rejected"]], 0, label = kind)

  values      <- chain$draws
  observables <- cbind(mean = values, variance = (values - mean(values))^2)
  ess         <- as.numeric(coda::effectiveSize(coda::mcmc(observables)))
  mcse        <- sqrt(apply(observables, 2L, stats::var) / ess)
  estimate    <- c(mean(values), stats::var(values))
  reference   <- c(exact$mean, exact$variance)
  uncorrected <- c(surrogate_only$mean, surrogate_only$variance)
  expect_true(all(ess >= 200), info = kind)
  # the surrogate-only target is another distribution, at least in the mean
  expect_gte(abs(reference[1L] - uncorrected[1L]), 8 * mcse[1L])
  expect_true(all(abs(estimate - reference) <= 6 * mcse + 2 * c(exact$mean_error, exact$variance_error)),
              info = kind)
  invisible(list(estimate = estimate, reference = reference, surrogate_only = uncorrected,
                 mcse = mcse, ess = ess, counters = counters))
}


test_that("the corrected chain is exact when the residual SD of a cluster block is sampled", {

  result <- .expect_sampled_parameter_chain_is_exact("cluster")
})


test_that("the corrected chain is exact when a loading column of a factor block is sampled", {

  result <- .expect_sampled_parameter_chain_is_exact("factor")
})
