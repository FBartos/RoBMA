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
    chain <- .location_hand_chain(model, c("g", "b", "Z"), draws = draws, burnin = 1000L, adapt = 500L,
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
  # an ellipse around the slab prior means, and one without the latent shear
  # (the joint moment of the intercept and the first cluster effect), are not
  for (mutant in 4:5) {
    broken <- run(mutant, ellipse_only)
    expect_false(attr(broken$comparison, "passes"),
                 info = paste("ellipse mutant", mutant, ":", paste(round(broken$comparison$z, 2), collapse = " ")))
  }
})


test_that("the location block over step-selection rows targets the exact posterior in both effect directions", {

  skip_if_not_certification("The step-switch fixtures (about 1 min) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  # Independent dselnorm_step_switch rows (the children of RoBMA()
  # meta-regressions, which the automatic policy gives the block): the
  # closed-form tilt, the screen's step branch and the block's own row density
  # in the inner target and the exact correction, with row observation weights
  # and a bias-exempt row; the reference writes the row density in R.
  draws <- 100000L
  for (direction in c(1, -1)) {
    setup <- .location_step_setup(direction)
    model <- .location_step_model(setup)
    reference <- .location_step_reference(setup, 30L)
    expect_lt(max(abs(reference - .location_step_reference(setup, 20L))), 1e-6)
    run <- function(mutant) {
      chain <- .location_hand_chain(model, c("g", "b"), draws = draws, burnin = 1000L, adapt = 500L, mutant = mutant)
      list(chain = chain, comparison = .location_compare(.location_selection_observables(chain$draws), reference))
    }
    exact <- run(0L)
    label <- paste("direction", direction)
    expect_identical(exact$chain$decision, "claimed without latents", label = label)
    expect_true(attr(exact$comparison, "passes"), info = paste(label, ":", paste(round(exact$comparison$z, 2), collapse = " ")))
    counters <- exact$chain$counters
    expect_gt(counters[["swap_accepted"]], 0.05 * draws, label = label)
    expect_gt(counters[["flip_on_accepted"]] + counters[["flip_off_accepted"]], 0.2 * draws, label = label)
    if (direction == 1) {
      # power: flips without their reverse proposal density
      broken <- run(2L)
      expect_false(attr(broken$comparison, "passes"),
                   info = paste("mutant 2:", paste(round(broken$comparison$z, 2), collapse = " ")))
    }
  }
})


test_that("the location block over dense multivariate selection blocks targets the exact posterior", {

  skip_if_not_certification("The dense fixture (about 1 min) runs in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  # dselnorm_mnorm_step children (known-V fits), whose density the reference
  # writes in R with a one-dimensional normalizer per block
  setup <- .location_dense_setup()
  model <- .location_dense_model(setup)
  reference <- .location_dense_reference(setup, 20L)
  expect_lt(max(abs(reference - .location_dense_reference(setup, 14L))), 1e-6)
  run <- function(mutant) {
    chain <- .location_hand_chain(model, c("g", "b"), draws = 100000L, burnin = 1000L, adapt = 500L, mutant = mutant)
    list(chain = chain, comparison = .location_compare(.location_dense_observables(chain$draws), reference))
  }
  exact <- run(0L)
  expect_identical(exact$chain$decision, "claimed without latents")
  expect_true(attr(exact$comparison, "passes"), info = paste(round(exact$comparison$z, 2), collapse = " "))
  expect_gt(exact$chain$counters[["surrogate_evaluations"]], 0)
  broken <- run(2L)
  expect_false(attr(broken$comparison, "passes"),
               info = paste("mutant 2:", paste(round(broken$comparison$z, 2), collapse = " ")))
})


test_that("a two-coordinate slab, the forced slice refresh and the exact mode after the fallback are exact", {

  skip_if_not_certification("The slab, slice and fallback fixtures (about 3 min) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  draws <- 100000L
  # A gated group of two coordinates (correlated dmnorm slab) in selection
  # labels: the multivariate Newton fits, the floored precision's eigen path,
  # the multivariate proposal densities and swaps between groups of different
  # sizes; with refresh variant 5 every group's refresh is the slice update,
  # a composite of two coordinate updates that reverses in reversed sweeps.
  setup <- .location_selection_setup()
  model <- .location_factor_model(setup)
  reference <- .location_factor_reference(setup, model, coefficient_nodes = 20L, latent_nodes = 40L)
  expect_lt(max(abs(reference - .location_factor_reference(setup, model))), 1e-4)
  run <- function(config = NULL, mutant = 0L) {
    chain <- .location_hand_chain(model, c("g1", "g2", "b1", "b2"), draws = draws, burnin = 1000L, adapt = 500L,
                                  mutant = mutant, config = config)
    list(chain = chain, comparison = .location_compare(.location_factor_observables(chain$draws), reference))
  }
  for (variant in list(NULL, c(refresh = 5L, fit_rule = -1L, timing = -1L))) {
    exact <- run(variant)
    label <- if (is.null(variant)) "default refresh" else "forced slice refresh"
    expect_identical(exact$chain$decision, "claimed", label = label)
    expect_true(attr(exact$comparison, "passes"), info = paste(label, ":", paste(round(exact$comparison$z, 2), collapse = " ")))
    expect_gt(exact$chain$counters[["swap_accepted"]], 0.05 * draws, label = label)
  }

  # The selection fixture after the fallback to exact mode (the corrected
  # samplers' test flag): no surrogate is evaluated and every correction is
  # one up to rounding.
  model <- .location_selection_model(setup)
  reference <- .location_selection_reference(setup, model, coefficient_nodes = 20L, latent_nodes = 40L)
  chain <- .with_corrected_flags(list(force_fallback = TRUE),
    .location_hand_chain(model, c("g", "b", "Z"), draws = 50000L, burnin = 1000L, adapt = 500L))
  comparison <- .location_compare(.location_selection_observables(chain$draws), reference)
  expect_true(attr(comparison, "passes"), info = paste(round(comparison$z, 2), collapse = " "))
  expect_identical(chain$counters[["surrogate_evaluations"]], 0)
  expect_identical(chain$counters[["sweeps_evidence"]], 0)
  expect_lt(chain$counters[["max_correction"]], 1e-10)
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
    # (nested layout) without the nested node's density in its acceptance
    # ratio is not exact (a term left out of both the ratio and the recorded
    # inner-target change is restored by the sweep's exact correction)
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


# The heterogeneity scale move (src/samplers/LocationScale.cc, fixtures in
# helper-location-scale.R): the location sweep, the sources' embedded ordinary
# updates and the scale move together target the exact posterior. A Gaussian
# fixture (four clusters, a sampled and an integrated level, a half-normal
# total SD, an inclusion indicator, gamma weight sources) has a two-dimensional
# quadrature reference; the selection fixture above with a random cluster SD
# (L = 1) composes the move with the corrected location sweep and a corrected
# source update. Same rule as above.
test_that("the scale move keeps the exact posterior of the heterogeneity and broken moves do not", {

  skip_if_not_certification("The scale move fixtures (about 2 min) run in the numerical-kernels case.")
  skip_if_not_installed("rjags")
  setup <- .scale_allocation_setup(shapes = c(1, 1))
  model <- .scale_allocation_model(setup)
  reference <- .scale_allocation_reference(setup)
  expect_lt(max(abs(reference - .scale_allocation_reference(setup, a_nodes = 80L, w_nodes = 80L))), 1e-6)
  gaussian <- function(mutant) {
    compiled <- .scale_compile(model, list(I = 1), mutant = mutant, seed = 9157L)
    run <- .scale_run(compiled, c("g", "I", "a", "w", "sigma"), draws = 200000L, burnin = 1000L, adapt = 500L)
    list(decision = compiled$decisions, report = run$report,
         comparison = .location_compare(.scale_allocation_observables(run$draws), reference))
  }
  exact <- gaussian(0L)
  expect_identical(exact$decision$status, "owned")
  expect_true(attr(exact$comparison, "passes"), info = paste(round(exact$comparison$z, 2), collapse = " "))
  # the move runs while the heterogeneity is included and is the identity otherwise
  expect_gt(exact$report$counters[["moves"]], 0.5 * 200000)
  expect_gt(exact$report$counters[["identity"]], 0)
  # the exponent c^(1 - K), the factor a'^-(2L - 1) omitted, the gamma total
  # not preserved
  for (mutant in c(1L, 2L, 4L)) {
    broken <- gaussian(mutant)
    expect_false(attr(broken$comparison, "passes"),
                 info = paste("mutant", mutant, ":", paste(round(broken$comparison$z, 2), collapse = " ")))
  }
  # one latent left out of the transformation changes the children: the
  # move's check undoes it with the classed error
  expect_error(gaussian(3L), "RoBMA_location_scale_invariance", fixed = TRUE)

  options <- RoBMA.options()[c("selection.sampler", "selection.coarse_max_rules", "selection.coarse_min_acceptance")]
  on.exit(do.call(RoBMA.options, options), add = TRUE)
  RoBMA.options(selection.sampler = "slice.coarse", selection.coarse_max_rules = 3L,
                selection.coarse_min_acceptance = 0)
  setup <- .location_selection_setup()
  model <- .scale_selection_model(setup)
  reference <- .scale_selection_reference(setup, model)
  expect_lt(max(abs(reference - .scale_selection_reference(setup, model, 30L, 16L, 32L))), 1e-4)
  selection <- function(mutant) {
    compiled <- .scale_compile(model, list(I = 1), mutant = mutant, seed = 4127L)
    run <- .scale_run(compiled, c("g", "I", "b", "a", "sigma"), draws = 100000L, burnin = 1000L, adapt = 500L)
    list(decision = compiled$decisions,
         comparison = .location_compare(.scale_selection_observables(run$draws), reference))
  }
  exact <- selection(0L)
  expect_identical(exact$decision$status, "owned")
  expect_identical(exact$decision$embedded, "RoBMA::CoarseCorrectedSlice")
  expect_true(attr(exact$comparison, "passes"), info = paste(round(exact$comparison$z, 2), collapse = " "))
  broken <- selection(1L)
  expect_false(attr(broken$comparison, "passes"),
               info = paste("mutant 1:", paste(round(broken$comparison$z, 2), collapse = " ")))
})
