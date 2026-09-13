test_that("Gaussian zplot reference matches dense parts for all source choices", {

  data <- data.frame(yi = c(.1, -.2), paper = c("a", "a"))
  V <- matrix(c(.16, .06, .06, .49), 2L)
  samples <- matrix(0, 2L, 1L)
  mean <- matrix(c(-.4, .7, .2, -.1), 2L)
  # No random formula is supplied: this full-V model has no extra variance.
  within <- matrix(0, 2L, 2L)
  predictive <- list(mu = mean, mu_extrapolated = mean + .25,
    tau_within = within, sei = sqrt(diag(V)))
  expected_variance <- sweep(within^2, 2L, diag(V), "+")
  for (estimate in c("integrate", "condition")) {
    for (sampling in c("integrate", "condition")) {
      object <- bselmodel.mv(yi = yi, V = V, data = data, measure = "GEN",
        prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE,
        selection = selection_model(estimate_random_effects = estimate,
          known_sampling_variance = sampling, group = paper))
      dense <- .predict_joint_selection_gaussian_parts(
        object, object$data, samples, mean, within, matrix(0, 2L, 2L),
        draw_context = FALSE)
      dense_diagonal <- vapply(seq_len(2L), function(row) {
        dense$covariance[, row, row] + dense$context_covariance[, row, row]
      }, numeric(2L))
      reference <- .zplot_gaussian_marginal_reference(object, samples, predictive)
      expect_identical(reference$mu, predictive$mu_extrapolated)
      expect_identical(reference$sei, predictive$sei)
      expect_equal(unname(reference$variance), expected_variance, tolerance = 1e-14)
      expect_equal(unname(reference$variance), dense_diagonal, tolerance = 1e-14)
      z <- c(-3, -.2, 1, 4)
      expected_density <- vapply(z, function(value) {
        rowMeans(stats::dnorm(value,
          sweep(predictive$mu_extrapolated, 2L, predictive$sei, "/"),
          sweep(sqrt(expected_variance), 2L, predictive$sei, "/")))
      }, numeric(2L))
      expect_equal(.zplot_normal_density_matrix(z, reference$mu,
        sqrt(reference$variance), reference$sei), expected_density, tolerance = 1e-14)
    }
  }
})

test_that("Gaussian reference includes fresh estimate and study variation", {

  data <- data.frame(yi = c(.1, -.2, .3), sei = c(.3, .4, .7),
    paper = c("a", "a", "b"))
  samples <- matrix(0, 2L, 1L)
  mean <- matrix(c(-.4, .7, .2, -.1, 0, .3), 2L)
  within <- matrix(c(0, .3, .5, .8, .2, .4), 2L)
  between <- matrix(c(.2, .6), 2L, 3L)
  extra_variance <- within^2 + between^2
  predictive <- list(mu_extrapolated = mean, tau_within = sqrt(extra_variance),
    sei = data$sei)
  expected <- sweep(extra_variance, 2L, data$sei^2, "+")
  for (estimate in c("integrate", "condition")) {
    for (study in c("integrate", "condition")) {
      for (sampling in c("integrate", "condition")) {
        object <- bselmodel(yi = yi, sei = sei, cluster = paper,
          data = data, measure = "GEN", prior_unit_information_sd = 1,
          only_priors = TRUE, silent = TRUE,
          selection = selection_model(estimate_random_effects = estimate,
            other_random_effects = study, known_sampling_variance = sampling,
            group = paper))
        dense <- .predict_joint_selection_gaussian_parts(
          object, object$data, samples, mean, within, between, draw_context = FALSE)
        dense_diagonal <- vapply(seq_len(3L), function(row) {
          dense$covariance[, row, row] + dense$context_covariance[, row, row]
        }, numeric(2L))
        reference <- .zplot_gaussian_marginal_reference(object, samples, predictive)
        expect_equal(reference$variance, expected, tolerance = 1e-14)
        expect_equal(reference$variance, dense_diagonal, tolerance = 1e-14)
      }
    }
  }
})

test_that("Gaussian references reuse prepared diagonal heterogeneity without dense reconstruction", {

  samples <- matrix(0, 2L, 1L)
  object <- bselmodel.mv(yi = c(0, 0), V = diag(c(.3, .7)^2),
    data = data.frame(paper = c("a", "a")), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE,
    selection = selection_model(group = paper))
  predictive <- list(mu_extrapolated = matrix(c(1, 2), 2L, 2L),
    tau_within = matrix(c(0, .5), 2L, 1L), sei = c(.3, .7))
  helper <- .zplot_gaussian_marginal_reference
  helper_calls <- 0L
  testthat::local_mocked_bindings(
    .zplot_gaussian_marginal_reference = function(...) {
      helper_calls <<- helper_calls + 1L
      helper(...)
    },
    .predict_joint_selection_gaussian_parts = function(...) stop("Dense parts must not be requested."),
    .brma_mv_random_effects_marginal_vcov = function(...) stop("Prepared diagonals must be reused."),
    .zplot_predictive_heterogeneity = function(...) stop("Prepared heterogeneity must be reused."),
    .package = "RoBMA"
  )
  reference <- .zplot_gaussian_marginal_reference(object, samples, predictive)
  expect_equal(reference$variance, matrix(c(.09, .34, .49, .74), 2L), tolerance = 1e-14)
  expect_identical(reference$mu, predictive$mu_extrapolated)
  control <- set_selection_likelihood_control()
  z <- c(-2, 0, 2)
  density <- .zplot_joint_marginal(object, samples, predictive, NULL,
    z, FALSE, control, extrapolate_only = TRUE)
  expect_null(density$fitted)
  expect_null(density$EDR)
  expect_identical(density$weights, c(1, 1))
  expect_equal(density$extrapolated, .zplot_normal_density_matrix(
    z, reference$mu, sqrt(reference$variance), reference$sei), tolerance = 1e-14)
  tail <- .zplot_joint_marginal(object, samples, predictive, NULL,
    2, TRUE, control, extrapolate_only = TRUE)
  threshold <- matrix(2 * reference$sei, 2L, 2L, byrow = TRUE)
  expected_tail <- rowMeans(stats::pnorm(threshold, reference$mu,
    sqrt(reference$variance), lower.tail = FALSE) +
    stats::pnorm(-threshold, reference$mu, sqrt(reference$variance)))
  expect_null(tail$fitted)
  expect_equal(tail$EDR, expected_tail, tolerance = 1e-14)
  expect_equal(as.numeric(tail$extrapolated), expected_tail, tolerance = 1e-14)
  expect_identical(helper_calls, 3L)
})

test_that("selected full-V zplot chunks retain density tails and posterior row order", {

  # Same bivariate rank-one full-V construction as the independent selected
  # law test in test-00-zplot-compute.R; no fitted model or mocked numerics.
  S <- 5L
  sei <- c(.7, 1.3)
  rho <- .7
  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)),
    model = selection_model(other_random_effects = "integrate",
      known_sampling_variance = "integrate", group = "paper"))
  object <- bselmodel.mv(yi = c(0, 0),
    data = data.frame(paper = c("p1", "p1")),
    V = known_v_factor((1 - rho) * sei^2, matrix(sqrt(rho) * sei, ncol = 1)),
    prior_bias = prior, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)
  selection <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), sei,
    effect_direction = "positive", signed_data = FALSE)
  selection$omega <- matrix(rep(c(1, .2), S), S, 2L, byrow = TRUE)
  selection$alpha <- numeric(S)
  selection$phack_kind <- integer(S)
  selection$kernel_mode <- rep(SELKERNEL_STEP, S)
  selection$vector_rule <- integer(S)
  selection$use_normal <- rep(FALSE, S)
  samples <- matrix(0, S, 1L)
  mean <- matrix(c(-.7, .3, 1.1, -.2, .8, .2, -.4, .9, -.8, .1), S)
  predictive <- list(mu = mean, mu_extrapolated = mean,
    tau_within = matrix(0, S, 2L), sei = sei)
  control <- set_selection_likelihood_control(points_per_scramble = 256L,
    max_points_per_scramble = 16384L, relative_tolerance = .002)
  # Four coexisting arrays, with the existing peak estimate admitting two
  # rows per recursive call: chunks are 1:2, 3:4, and 5.
  chunk_budget <- 4 * .known_v_covariance_peak_bytes(2L, 2L)
  withr::local_options(list(RoBMA.known_v_covariance_max_bytes = Inf))
  permutation <- c(5L, 2L, 4L, 1L, 3L)
  permuted_predictive <- predictive
  for (name in c("mu", "mu_extrapolated", "tau_within")) {
    permuted_predictive[[name]] <- predictive[[name]][permutation, , drop = FALSE]
  }
  for (probability in c(FALSE, TRUE)) {
    z <- if (probability) stats::qnorm(.975) else c(-2, -.3, .8, 3)
    options(RoBMA.known_v_covariance_max_bytes = Inf)
    whole <- .zplot_joint_marginal(object, samples, predictive, selection,
      z, probability, control)
    options(RoBMA.known_v_covariance_max_bytes = chunk_budget)
    expect_identical(.known_v_covariance_chunk_indices(S, 2L,
      max_bytes = chunk_budget / 4), list(1:2, 3:4, 5L))
    chunked <- .zplot_joint_marginal(object, samples, predictive, selection,
      z, probability, control)
    expect_equal(chunked, whole, tolerance = 1e-12)
    expect_gt(length(unique(whole$fitted[, 1L])), 1L)
    reordered <- .zplot_joint_marginal(object, samples[permutation, , drop = FALSE],
      permuted_predictive, BayesTools::selection_context_subset_rows(selection, permutation),
      z, probability, control)
    expect_equal(reordered$fitted, whole$fitted[permutation, , drop = FALSE], tolerance = 1e-12)
    expect_equal(reordered$extrapolated, whole$extrapolated[permutation, , drop = FALSE], tolerance = 1e-12)
    expect_equal(reordered$weights, whole$weights[permutation], tolerance = 1e-12)
    if (probability) expect_equal(reordered$EDR, whole$EDR[permutation], tolerance = 1e-12) else
      expect_null(reordered$EDR)
  }
})
