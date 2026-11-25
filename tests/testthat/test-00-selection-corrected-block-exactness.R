probe_state("test-00-selection-corrected-block-exactness.R")
context("Exactness of the corrected block sampler")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The corrected block sampler makes one factor slice sweep on the surrogate in a
# random direction order (forward or backward, each with probability 1/2) and
# accepts it with the ratio of the full conditional to the surrogate; its
# invariant law is the exact posterior whatever the quality of the surrogate.
# The fixtures (helper-selection-corrected-block.R) have exact references by
# QUADPACK: a power fixture (four members, two of them bounded, one the parent
# of another, strongly correlated, a crude nonseparable surrogate, no
# adaptation, the corrected kernel throughout), a stationary start from exact
# draws, a dense selection child; and the surrogate's history independence.
# Certification case 'numerical-kernels'.


test_that("the corrected block chain targets the exact posterior and a fixed-order sweep does not", {

  skip_if_not_certification("The power fixture (about 45 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  setup <- .block_power_setup()

  # the surrogate (rule 3 of the crude schedule) is far from the exact normalizer
  # across the bulk of mu = a + b (about -2.7, SD .5): log gaps of 1 to 5 per child
  rungs <- .coarse_exactness_rungs(setup$cluster, c(-3.5, -3, -2.7, -2.4, -2))
  expect_true(all(rungs[["gap"]] > 1))
  tracker <- new.env()
  tracker$relative_error <- 0
  exact <- .block_power_reference(setup, function(mu) {
    .coarse_exactness_quadpack_log_normalizer(setup$cluster, mu, tracker)
  }, tracker)
  expect_true(all(exact$allowance < 1e-4))
  expect_lt(exact$correlation_ab, -.7)

  draws <- 80000L
  members <- c("a", "b", "s", "r")
  corrected <- .block_chain(.block_power_model(setup), members, draws, adapt = 0L, burnin = 1000L)
  expect_identical(unname(corrected$samplers[members]), rep("RoBMA::CoarseCorrectedBlock", 4L))
  expect_identical(corrected$counters[["proposed"]], as.numeric(draws))
  expect_identical(corrected$counters[["accepted"]] + corrected$counters[["correction_rejected"]], as.numeric(draws))
  expect_gte(corrected$counters[["correction_rejected"]], .1 * draws)
  comparison <- .block_compare(.block_power_observables(corrected$draws), exact)
  expect_true(all(comparison$ess >= 400), info = paste(round(comparison$ess), collapse = " "))
  expect_true(all(abs(comparison$estimate - comparison$reference) <= 4.5 * comparison$mcse + 2 * comparison$allowance),
              info = paste(round(comparison$z, 2), collapse = " "))

  # power: the sweep in its fixed order with the same correction is not exact
  mutant <- .block_chain(.block_power_model(setup, fixed_order = TRUE), members, draws, adapt = 0L, burnin = 1000L)
  biased <- .block_compare(.block_power_observables(mutant$draws), exact)
  expect_gt(max(abs(biased$z)), 6)
  invisible(list(corrected = comparison, mutant = biased))
})


test_that("one update from exact draws of independent blocks leaves their law unchanged", {

  skip_if_not_certification("The stationary start (about 20 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .block_stationary_setup()
  n     <- 4000L
  start <- withr::with_seed(20261119L, .block_stationary_draws(setup, n))
  tracker <- new.env()
  tracker$relative_error <- 0
  reference <- .block_stationary_reference(setup, tracker)
  expect_true(all(reference$allowance < 1e-4))
  # iid moments: standardized distances from the reference
  iid_z <- function(a, b) {
    observables <- .block_ab_observables(cbind(a = a, b = b))
    estimate <- colMeans(observables)
    estimate[c("var_a", "var_b")] <- estimate[c("var_a", "var_b")] * n / (n - 1)
    (estimate - reference$estimate[colnames(observables)]) / (apply(observables, 2L, stats::sd) / sqrt(n))
  }
  # the exact draws are draws of the reference law
  expect_true(all(abs(iid_z(start$a, start$b)) <= 4.5))

  hand  <- .block_stationary_model(setup, start)
  expect_true(all(hand$samplers[paste0("a[", seq_len(n), "]")] == "RoBMA::CoarseCorrectedBlock"))
  chain <- .block_chain(hand, c("a", "b"), 1L, adapt = 0L, burnin = 0L)
  expect_identical(chain$counters[["proposed"]], as.numeric(n))
  expect_gte(chain$counters[["correction_rejected"]], .2 * n)
  a <- unname(chain$draws[1L, paste0("a[", seq_len(n), "]")])
  b <- unname(chain$draws[1L, paste0("b[", seq_len(n), "]")])
  expect_true(all(abs(iid_z(a, b)) <= 4.5))
  # paired: the mean change of each observable is zero
  change <- .block_ab_observables(cbind(a = a, b = b)) - .block_ab_observables(cbind(a = start$a, b = start$b))
  paired_z <- colMeans(change) / (apply(change, 2L, stats::sd) / sqrt(n))
  expect_true(all(abs(paired_z) <= 4.5), info = paste(round(paired_z, 2), collapse = " "))
})


test_that("the corrected block chain is exact on a dense selection child with a crude anchor grid", {

  skip_if_not_certification("The dense fixture (about 15 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  setup   <- .block_dense_setup()
  tracker <- new.env()
  tracker$relative_error <- 0
  exact <- .block_dense_reference(setup, tracker)
  expect_true(all(exact$allowance < 1e-4))
  chain <- .block_chain(.block_dense_model(setup), c("a", "b"), 30000L, adapt = 0L, burnin = 500L)
  expect_identical(unname(chain$samplers[c("a", "b")]), rep("RoBMA::CoarseCorrectedBlock", 2L))
  expect_gt(chain$counters[["correction_rejected"]], 0)
  comparison <- .block_compare(.block_ab_observables(chain$draws), exact)
  expect_true(all(comparison$ess >= 400), info = paste(round(comparison$ess), collapse = " "))
  expect_true(all(abs(comparison$estimate - comparison$reference) <= 4.5 * comparison$mcse + 2 * comparison$allowance),
              info = paste(round(comparison$z, 2), collapse = " "))
})


# Every surrogate evaluation of a traced chain at the same member values has the
# same value: revisits after other states, after full evaluations and after
# rejections (the next update starts where the rejected one did), and in a new
# model with cleared selection caches (cold memos) started at a visited state.
# `build(inits)` compiles the fixture (the trace starts after compilation).
.block_expect_history_independent <- function(build, members, label) {

  on.exit(.block_test_hook(trace = FALSE), add = TRUE)
  hand <- build(NULL)
  .block_test_hook(trace = TRUE)
  chain <- .block_chain(hand, members, 300L, adapt = 0L, burnin = 0L)
  .block_test_hook(trace = FALSE)
  trace <- .block_surrogate_trace(chain$trace)
  expect_gt(chain$counters[["correction_rejected"]], 0, label = label)
  revisits <- tapply(trace$value, trace$key, length)
  expect_gt(sum(revisits > 1L), 10L, label = label)
  spread <- tapply(trace$value, trace$key, function(values) length(unique(values)))
  expect_true(all(spread == 1L), info = label)

  # a cold start at the state after draw 100
  RoBMA::selection_cache_info(clear = TRUE)
  cold <- build(as.list(chain$draws[100L, members]))
  .block_test_hook(trace = TRUE)
  rjags::adapt(cold$model, 0L, end.adaptation = TRUE)
  stats::update(cold$model, n.iter = 1L, progress.bar = "none")
  first <- .block_surrogate_trace(.block_test_hook(trace = FALSE)$records)[1L, ]
  expect_true(first$key %in% trace$key, info = label)
  expect_identical(first$value, trace$value[match(first$key, trace$key)], info = label)
}


test_that("the block surrogate does not depend on the evaluation history or on warm memos", {

  skip_if_not_certification("These traced chains (about 5 s) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  power <- .block_power_setup()
  .block_expect_history_independent(function(inits) {
    if (is.null(inits)) .block_power_model(power) else .block_power_model(power, inits = inits)
  }, c("a", "b", "s", "r"), "cluster children")
  dense <- .block_dense_setup()
  .block_expect_history_independent(function(inits) {
    if (is.null(inits)) .block_dense_model(dense) else .block_dense_model(dense, inits = inits)
  }, c("a", "b"), "dense child")
})
