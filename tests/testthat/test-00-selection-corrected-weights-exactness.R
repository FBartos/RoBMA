context("Exactness of the corrected samplers of selection weights and indicators")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The gated weight block, the screened indicator and the gamma-sum coordinates
# leave the exact posterior invariant for any surrogate that is one function
# of the evaluated state (the crude schedule's three-node rule, the row cache
# and its perturbations). The references integrate the package's exact
# likelihood of the children (the distribution class the sampler evaluates)
# over the weights and the allocation by fixed Gauss-Legendre rules
# (helper-selection-corrected-weights.R).
#
# Rule, fixed before the runs: each observable's z-score is its distance from
# the reference over the standard error of 50 batch means (and the reference's
# allowance); an exact kernel passes when every required observable is
# compared and every |z| is finite and at most 4 (a frozen chain, with no
# variance, fails; only the predeclared structural constant of a fixture with
# one label, its indicator, is compared by value), and a mutant of the test
# hook (the gamma groups without the Jacobian of the proportions, the
# screened indicator without its correction, a fixed sweep order) is detected
# when some |z| >= 6. Certification case 'numerical-kernels'.

.weights_strong_x <- rbind(c(2.5, 2.3, 2.8, 2.1, 3.0, 1.2, .5, -.4),
                           c(2.2, 2.6, 3.1, 2.0, 2.4, .9, 1.4, -.9),
                           c(2.9, 2.1, 2.4, 2.7, 1.5, .3, -.2, -1.1))

# z-scores of a weights-fixture chain against the reference.
.weights_z <- function(setup, reference, draws, flags = list(), fixed_order = FALSE) {

  monitor <- c("b", if (2L %in% setup$branches) "r", if (3L %in% setup$branches) "eta")
  chain <- .with_corrected_flags(flags, .block_chain(.weights_model(setup, fixed_order = fixed_order), monitor,
                                                     draws, adapt = 500L, burnin = 500L))
  observables <- .weights_observables(chain$draws, setup)
  targets <- c(reference$estimate, if (3L %in% setup$branches) .weights_total_reference(setup))
  allowance <- c(reference$allowance, if (3L %in% setup$branches) c(T = 0, T2 = 0))
  # a single label's indicator is a structural constant
  .weights_compare(observables, targets, allowance, constants = .weights_constants(setup))
}

.weights_expect_exact <- function(comparison, label) {

  expect_true(.weights_is_exact(comparison),
              info = paste(label, paste(sprintf("%s=%.2f", comparison$quantity, comparison$z), collapse = " ")))
}

.weights_expect_detected <- function(comparison, label) {

  expect_gte(max(abs(comparison$z), na.rm = TRUE), 6, label = paste(label, "largest |z|"))
}


test_that("the comparison fails a frozen chain and needs every observable", {

  # a gamma total frozen at 7 against its reference 3 (no variance, no
  # allowance) and a mixing indicator
  observables <- cbind(T = rep(7, 100), I1 = rep(c(0, 1, 1, 1), 25))
  frozen <- .weights_compare(observables, c(T = 3, I1 = .75))
  expect_identical(frozen$z[frozen$quantity == "T"], Inf)
  expect_false(.weights_is_exact(frozen))
  # frozen at its reference value: undefined, not exact unless declared constant
  observables[, "T"] <- 3
  expect_false(.weights_is_exact(.weights_compare(observables, c(T = 3, I1 = .75))))
  expect_true(.weights_is_exact(.weights_compare(observables, c(T = 3, I1 = .75), constants = "T")))
  # a declared constant is checked against its reference value
  constant <- .weights_compare(observables, c(T = 2, I1 = .75), constants = "T")
  expect_identical(constant$z[constant$quantity == "T"], Inf)
  expect_false(.weights_is_exact(constant))
  # every required observable is compared, and every observable has a reference
  expect_error(.weights_compare(observables, c(T = 3, I1 = .75, I2 = .1)), "needs a reference")
  expect_error(.weights_compare(observables, c(T = 3)), "needs a reference")
  expect_false(.weights_is_exact(frozen[0L, ]))
})


test_that("the three-branch ensemble is exact and its mutants are not", {

  skip_if_not_certification("The three-branch ensemble (about 60 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .weights_setup(1:3, x = .weights_strong_x, correlation = .9)
  reference <- .weights_reference(setup)
  expect_true(all(reference$allowance < 1e-3))
  # every label is visited often
  expect_true(all(reference$estimate[c("I1", "I2", "I3")] > .02))

  draws <- 40000L
  # the default screened proposal mixes Liu's step and the plain surrogate
  # draw; each alone stays exact
  .weights_expect_exact(.weights_z(setup, reference, draws), "corrected, default mixture")
  .weights_expect_exact(.weights_z(setup, reference, draws, list(finite_proposal = 2)), "Liu's step")
  .weights_expect_exact(.weights_z(setup, reference, draws, list(finite_proposal = 0)), "plain surrogate draw")
  .weights_expect_detected(.weights_z(setup, reference, draws, list(finite_proposal = 1)), "uncorrected indicator")
  .weights_expect_detected(.weights_z(setup, reference, draws, list(drop_sum_jacobian = TRUE)), "no Jacobian")
})


test_that("the two-branch ensemble is exact", {

  skip_if_not_certification("The two-branch ensemble (about 20 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .weights_setup(c(1L, 3L), prior = c(.5, .5))
  reference <- .weights_reference(setup)
  expect_true(all(reference$estimate[c("I1", "I3")] > .1))
  .weights_expect_exact(.weights_z(setup, reference, 40000L), "corrected")
})


test_that("the gated block is exact in its random order and not in a fixed one", {

  skip_if_not_certification("The open gated block (about 70 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  # one open branch of eight-row blocks under strong selection: the block
  # accepts about half of its sweeps, and a fixed direction order with one
  # correction is biased by about 4 standard errors per 40,000 draws
  setup <- .weights_setup(3L, x = .weights_strong_x, correlation = .9)
  reference <- .weights_reference(setup)
  draws <- 200000L
  .weights_expect_exact(.weights_z(setup, reference, draws), "random order")
  .weights_expect_detected(.weights_z(setup, reference, draws, fixed_order = TRUE), "fixed order")
})


# Singleton step-selection rows (the row cache of helper-step-row-cache.R): the
# rows of the weights fixtures as independent dselnorm_step_switch rows at SD
# and selection SE one. The test hook perturbs the cache (each bin probability
# scaled by 1 + c (b + 1) / 3) so that the corrections reject a part of the
# proposals; the references integrate the rows' likelihood written in R.
.rows_mild_x <- rbind(c(2.4, 1.3, .6, -.3), c(2.1, 1.7, .9, -.8), c(2.9, 1.1, .2, -.5), c(2.2, 1.5, .4, -.2))

test_that("the three-branch ensemble of step-selection rows is exact with a perturbed cache and its mutants are not", {

  skip_if_not_certification("The step-row ensemble (about 10 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .rows_weights_setup(1:3, x = .rows_mild_x)
  reference <- .weights_reference(setup, children_log = .rows_children_log)
  expect_true(all(reference$allowance < 1e-3))
  expect_true(all(reference$estimate[c("I1", "I2", "I3")] > .02))

  draws <- 40000L
  perturbed <- list(perturbation = .25)
  corrected <- .rows_weights_z(setup, reference, draws, row_flags = perturbed)
  .weights_expect_exact(corrected, "corrected")
  # the perturbed cache leaves the corrections work to do
  kinds <- attr(corrected, "kinds")
  expect_lt(kinds["finite", "accepted"] / kinds["finite", "proposed"], .9)
  expect_lt(kinds["gated", "accepted"] / kinds["gated", "proposed"], .99)
  .weights_expect_exact(.rows_weights_z(setup, reference, draws, list(finite_proposal = 2), perturbed), "Liu proposal")
  .weights_expect_detected(.rows_weights_z(setup, reference, draws, list(finite_proposal = 1), perturbed),
                           "uncorrected indicator")
  .weights_expect_detected(.rows_weights_z(setup, reference, draws, list(drop_sum_jacobian = TRUE), perturbed),
                           "no Jacobian")
})


test_that("the nested weight layout with a PET label is exact with a perturbed cache and its mutant is not", {

  skip_if_not_certification("The nested PET fixture (about 10 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  # helper-step-row-cache.R: PET moves the rows' means, so a selection label's
  # surrogate must be one function of its state from every entry level
  setup <- .rows_nested_setup()
  reference <- .rows_nested_reference(setup)
  expect_true(all(reference$estimate[c("I1", "I2", "I3", "I4")] > .2))
  expect_true(all(reference$allowance < .01))
  draws <- 60000L
  perturbed <- list(perturbation = .25)
  corrected <- .rows_nested_z(setup, reference, draws, row_flags = perturbed)
  expect_identical(unname(attr(corrected, "samplers")[c("b", "a", "v", "u", "pet")]),
                   c("RoBMA::ScreenedFinite", rep("RoBMA::GatedCorrectedBlock", 3L), "BayesTools::GatedSlice"))
  .weights_expect_exact(corrected, "corrected, perturbed cache")
  # the perturbed cache leaves the corrections work to do
  kinds <- attr(corrected, "kinds")
  expect_lt(kinds["finite", "accepted"] / kinds["finite", "proposed"], .9)
  .weights_expect_detected(.rows_nested_z(setup, reference, draws, list(finite_proposal = 1), perturbed),
                           "uncorrected indicator")
})


test_that("the gated block of step-selection rows is exact in its random order and not in a fixed one", {

  skip_if_not_certification("The open step-row block (about 35 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  # one open branch of 24 rows under strong selection and a perturbed cache:
  # the block accepts about two thirds of its sweeps. A fixed direction order
  # with one correction biases the weights along its adapted directions, which
  # vary with the RNG stream (the one-label indicator's screened proposal
  # included); the log proportions show it in every chain (max |z| of the
  # weight moments 2.3-11.3, of the log proportions 13.7-30.5 over seeds 1-8 at
  # 400,000 draws; random order at most 2.7 over 18 quantities)
  setup <- .rows_weights_setup(3L, x = .weights_strong_x)
  setup$log_moments <- TRUE
  reference <- .weights_reference(setup, children_log = .rows_children_log)
  draws <- 400000L
  perturbed <- list(perturbation = 1)
  .weights_expect_exact(.rows_weights_z(setup, reference, draws, row_flags = perturbed), "random order")
  .weights_expect_detected(.rows_weights_z(setup, reference, draws, row_flags = perturbed, fixed_order = TRUE),
                           "fixed order")
})


test_that("the heterogeneity gate and the allocation block are exact", {

  skip_if_not_certification("The heterogeneity fixture (about 30 s) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .heterogeneity_setup()
  reference <- .heterogeneity_reference(setup)
  expect_true(reference$estimate[["h"]] > .1 && reference$estimate[["h"]] < .9)
  run <- function(flags = list()) {
    chain <- .with_corrected_flags(flags, .block_chain(.heterogeneity_model(setup), c("h", "s", "e"), 40000L,
                                                       adapt = 500L, burnin = 500L))
    .weights_compare(.heterogeneity_observables(chain$draws), c(reference$estimate, T = 2, T2 = 6),
                     c(reference$allowance, T = 0, T2 = 0))
  }
  .weights_expect_exact(run(), "corrected")
  .weights_expect_detected(run(list(drop_sum_jacobian = TRUE)), "no Jacobian")
})
