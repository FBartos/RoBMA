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
  run <- function(mutant, config = NULL) {
    chain <- .location_hand_chain(model, c("g", "b"), draws = draws, burnin = 1000L, adapt = 500L,
                                  mutant = mutant, config = config)
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

  # the elliptical slice step E as the only refresh of the active coefficients
  # is exact; an ellipse around the slab prior means instead of the Gaussian
  # working mean is not
  ellipse_only <- c(refresh = 2L, fit_rule = -1L, timing = -1L)
  ellipse <- run(0L, ellipse_only)
  expect_true(attr(ellipse$comparison, "passes"), info = paste(round(ellipse$comparison$z, 2), collapse = " "))
  expect_identical(ellipse$chain$counters[["refresh_proposed"]], 0)
  expect_gt(ellipse$chain$counters[["ellipse_proposed"]], 0.5 * draws)
  expect_gte(ellipse$chain$counters[["ellipse_accepted"]], 0.99 * ellipse$chain$counters[["ellipse_proposed"]])
  broken <- run(4L, ellipse_only)
  expect_false(attr(broken$comparison, "passes"),
               info = paste("ellipse mutant:", paste(round(broken$comparison$z, 2), collapse = " ")))
})


test_that("the location block's bias-label moves target the exact posterior and broken swaps do not", {

  skip_if_not_certification("The bias-label fixtures (about 6 min) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  # The fixture of helper-location-bias.R: labels none, step selection, PET and
  # PEESE (and a second selection label), a gated intercept with a nonzero
  # slab mean, a gated moderator, three clusters, a crude surrogate; once in a
  # nested layout (a beta step weight whose parameters depend on the label).
  # Reference: Gaussian labels in closed form with the slope integrated by
  # adaptive quadrature, selection labels by Gauss-Hermite in the cluster
  # effects and a Laplace-adapted grid in the coefficients (and Gauss-Legendre
  # in the nested weight). Rule fixed beforehand: every |z| < 4, at most one > 3.
  draws <- 100000L
  run <- function(setup, mutant = 0L, n = draws) {
    model <- .bias_model(setup)
    # L2 in every sweep (the test switch; the default schedule is one quarter)
    chain <- .location_hand_chain(model, c("bias_indicator", "g", "b", "PET", "PEESE"), draws = n,
                                  burnin = 1000L, adapt = 500L, mutant = mutant,
                                  config = c(refresh = -1L, fit_rule = -1L, timing = -1L, swap = 1L, trace = -1L))
    labels <- if (setup$second_selection) 5L else 4L
    list(chain = chain, model = model, labels = labels)
  }
  compare <- function(result, reference) {
    .location_compare(.bias_observables(result$chain$draws, result$labels), reference)
  }
  for (variant in list(list(), list(nested = TRUE), list(second_selection = TRUE))) {
    setup <- do.call(.bias_setup, variant)
    model <- .bias_model(setup)
    reference <- .bias_reference(setup, model, coefficient_nodes = 20L, latent_nodes = 40L, weight_nodes = 24L)
    coarser   <- .bias_reference(setup, model, coefficient_nodes = 16L, latent_nodes = 32L, weight_nodes = 16L)
    label <- paste(c("bias fixture", names(variant)), collapse = " ")
    expect_lt(max(abs(reference - coarser)), 1e-4, label = label)
    exact <- run(setup)
    expect_identical(exact$chain$decision, "claimed with the bias label")
    comparison <- compare(exact, reference)
    expect_true(attr(comparison, "passes"), info = paste(label, ":", paste(round(comparison$z, 2), collapse = " ")))
    counters <- exact$chain$counters
    expect_gt(counters[["bias_label_swap_accepted"]], 0.2 * draws, label = label)
    expect_gt(counters[["bias_gaussian_family_accepted"]], 0.05 * draws, label = label)
    expect_identical(counters[["bias_slope_accepted"]], counters[["bias_slope_proposed"]])
    # power: the swap without its label Hastings factor, without the gate
    # share, without the prior density of the slope that becomes dormant, and
    # (nested layout) without the nested node's density is not exact
    mutants <- if (isTRUE(variant$nested)) 9L else if (length(variant) == 0L) c(6L, 7L, 8L) else integer()
    for (mutant in mutants) {
      broken <- compare(run(setup, mutant), reference)
      expect_false(attr(broken, "passes"),
                   info = paste(label, "mutant", mutant, ":", paste(round(broken$z, 2), collapse = " ")))
    }
    if (length(variant) == 0L) {
      # the within-label contrast instead of the absolute evidence keeps the
      # swap reversible but wastes it: few swaps pass stage 2
      wasted <- run(setup, 10L, n = 20000L)$chain$counters
      expect_lt(wasted[["bias_label_swap_accepted"]] / 20000, 0.2 * counters[["bias_label_swap_accepted"]] / draws)
    }
  }
})
