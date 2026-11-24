context("Exactness of the location block sampler")
source(testthat::test_path("common-functions.R"))
skip_on_cran()

# The location block sampler (RoBMA::LocationBlock) composes tilted flips,
# refreshes, swaps and slice updates of the study effects on the corrected
# sampler's surrogate and accepts each sweep with one exact correction; its
# invariant law is the exact posterior. The fixture (helper-location-block.R)
# is a small selection meta-regression with three correlated clusters, a crude
# surrogate (rule 3 has three nodes), a gated intercept with a nonzero slab
# mean, a gated between-cluster moderator and a gated moderator collinear with
# the intercept; its reference integrates every gate pattern by Gauss-Hermite
# quadrature with the package's exact cluster likelihood. The rule is fixed
# beforehand: every |z| < 4 and at most one above 3 (batch means, 50 batches).
# Certification case 'numerical-kernels'.


test_that("the location block chain targets the exact posterior and broken kernels do not", {

  skip_if_not_certification("The location block fixture (about 3 min) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .location_selection_setup()
  model <- .location_selection_model(setup)
  reference <- .location_selection_reference(setup, model, coefficient_nodes = 20L, latent_nodes = 40L)
  coarser   <- .location_selection_reference(setup, model, coefficient_nodes = 16L, latent_nodes = 32L)
  # the quadrature error is far below the Monte Carlo error of the chains
  expect_lt(max(abs(reference - coarser)), 1e-4)

  draws <- 100000L
  run <- function(mutant) {
    chain <- .location_hand_chain(model, c("g", "b"), draws = draws, burnin = 1000L, adapt = 500L,
                                  mutant = mutant)
    list(chain = chain, comparison = .location_compare(.location_selection_observables(chain$draws), reference))
  }
  exact <- run(0L)
  expect_identical(exact$chain$decision, "claimed")
  expect_true(all(exact$chain$samplers[c(paste0("g[", 1:3, "]"), paste0("slab[", 1:3, "]"),
                                         paste0("Z[", 1:3, "]"))] == "RoBMA::LocationBlock"))
  expect_true(attr(exact$comparison, "passes"), info = paste(round(exact$comparison$z, 2), collapse = " "))
  counters <- exact$chain$counters
  # the surrogate is crude: the correction rejects a part of the sweeps
  expect_equal(counters[["sweeps_selection"]], draws)
  expect_lt(counters[["corrections_selection_accepted"]], 0.95 * draws)
  expect_gt(counters[["flip_on_accepted"]] + counters[["flip_off_accepted"]], 0.5 * draws)

  # power: an asymmetric shear, flips without their reverse proposal density
  # and the sweep without its exact correction are not exact
  for (mutant in 1:3) {
    broken <- run(mutant)
    expect_false(attr(broken$comparison, "passes"),
                 info = paste("mutant", mutant, ":", paste(round(broken$comparison$z, 2), collapse = " ")))
  }
})
