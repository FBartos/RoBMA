.zplot_density_reuse_context <- function(use_normal, sei) {

  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)))
  context <- .selection_spec(list(outcome = list(bias = prior)), rep(0, length(sei)), sei,
    effect_direction = "positive", signed_data = FALSE)
  S <- length(use_normal)
  context$omega <- matrix(rep(c(1, .2), S), S, byrow = TRUE)
  context$alpha <- numeric(S)
  context$phack_kind <- integer(S)
  context$kernel_mode <- ifelse(use_normal, SELKERNEL_NORMAL, SELKERNEL_STEP)
  context$vector_rule <- integer(S)
  context$use_normal <- use_normal
  context
}


test_that("latent mixtures retain full-matrix results while evaluating only normal rows", {

  sei <- c(.5, .8)
  mean <- matrix(c(.1, -.2, .3, 0, .2, -.1, .2, .1, -.2, .3), 5L,
    dimnames = list(paste0("draw", 1:5), c("a", "b")))
  sd <- matrix(seq(.6, 1, length.out = 10L), 5L)
  latent <- matrix(seq(0, .15, length.out = 10L), 5L)
  extrapolated_mean <- mean + .125
  total_sd <- .root_sum_squares(sd, latent)
  control <- set_selection_likelihood_control()
  original_normal <- .zplot_normal_density_matrix

  for (normal_rows in list(1:5, c(2L, 4L), integer())) {
    context <- .zplot_density_reuse_context(seq_len(5L) %in% normal_rows, sei)
    active <- which(!context$use_normal)
    for (probability in c(FALSE, TRUE)) {
      z <- if (probability) 2 else c(-3, -.5, 0, 2, 4)
      normal <- function(location) {
        if (!probability) return(original_normal(z, location, total_sd, sei))
        q <- matrix(z * sei, 5L, length(sei), byrow = TRUE)
        matrix(rowMeans(stats::pnorm(q, location, total_sd, lower.tail = FALSE) +
          stats::pnorm(-q, location, total_sd)), 5L, 1L)
      }
      for (fitted_only in c(FALSE, TRUE)) {
        # The original route evaluated both complete Gaussian matrices, then
        # overwrote selection rows. Retain that arithmetic as the reference.
        reference <- list(fitted = normal(mean), extrapolated = normal(extrapolated_mean),
          weights = rep(1, 5L))
        if (length(active)) {
          selected <- .zplot_latent_mixture(z, mean[active, , drop = FALSE],
            sd[active, , drop = FALSE], latent[active, , drop = FALSE], sei,
            BayesTools::selection_context_subset_rows(context, active), probability,
            control, extrapolated_mean[active, , drop = FALSE], fitted_only)
          reference$fitted[active, ] <- selected$fitted
          reference$extrapolated[active, ] <- selected$extrapolated
          reference$weights[active] <- selected$weights
        }
        reference$EDR <- if (probability && !fitted_only) {
          reference$extrapolated[, 1L] / reference$weights
        } else NULL
        calls <- list()
        actual <- testthat::with_mocked_bindings(
          .zplot_latent_mixture(z, mean, sd, latent, sei, context, probability,
            control, extrapolated_mean, fitted_only),
          .zplot_normal_density_matrix = function(z_sequence, mean, sd, sei) {
            calls[[length(calls) + 1L]] <<- mean
            original_normal(z_sequence, mean, sd, sei)
          }, .package = "RoBMA")
        expect_identical(actual$fitted, reference$fitted)
        expect_identical(actual$weights, reference$weights)
        expect_identical(actual$EDR, reference$EDR)
        if (!fitted_only) expect_identical(actual$extrapolated, reference$extrapolated)
        expected_calls <- if (probability || !length(normal_rows)) 0L else if (fitted_only) 1L else 2L
        expect_length(calls, expected_calls)
        if (length(calls)) {
          expect_identical(calls[[1L]], mean[normal_rows, , drop = FALSE])
          if (!fitted_only) expect_identical(calls[[2L]], extrapolated_mean[normal_rows, , drop = FALSE])
        }
      }
    }
  }
})


test_that("fitted-only mixtures skip the unused population precheck and retain density gates", {

  context <- .zplot_density_reuse_context(c(FALSE, TRUE, FALSE, TRUE, FALSE), .5)
  mean <- matrix(.1, 5L)
  sd <- matrix(.6, 5L)
  latent <- matrix(.1, 5L)
  control <- set_selection_likelihood_control()
  evaluate <- function(probability = FALSE) {
    .zplot_latent_mixture(if (probability) 0 else c(-1, 0, 1), mean, sd, latent,
      .5, context, probability, control, fitted_only = TRUE)
  }
  original_nodes <- .gauss_hermite_nodes_only
  original_rule <- .gauss_hermite_nodes
  requested <- integer()
  testthat::with_mocked_bindings({
    expect_true(all(is.finite(evaluate()$fitted)))
    expect_equal(as.numeric(evaluate(TRUE)$fitted), rep(1, 5L), tolerance = 1e-14)
  },
    .gauss_hermite_nodes_only = function(n) {
      requested <<- c(requested, n)
      original_nodes(n)
    },
    .gauss_hermite_nodes = function(n) {
      requested <<- c(requested, n)
      original_rule(n)
    }, .package = "RoBMA")
  expect_false(utils::tail(SELNORM_CLUSTER_QUADRATURE_ORDERS, 1L) %in% requested)

  cases <- list(
    list(error = NULL, resolved = TRUE, message = "Zplot density integration diagnostics are unavailable."),
    list(error = 0, resolved = NULL, message = "Zplot density selection-resolution diagnostics are unavailable."),
    list(error = .02, resolved = TRUE, message = "relative integration error was 0.02"),
    list(error = 0, resolved = FALSE, message = "selection transitions remain unresolved")
  )
  for (case in cases) {
    testthat::with_mocked_bindings(
      expect_error(evaluate(), case$message, fixed = TRUE),
      .zplot_selnorm_density_matrix = function(z_sequence, mean, ...) {
        result <- matrix(.5, nrow(mean), length(z_sequence))
        if (!is.null(case$error)) attr(result, "relative_integration_error") <- rep(case$error, nrow(mean))
        if (!is.null(case$resolved)) attr(result, "selection_transition_resolved") <- rep(case$resolved, nrow(mean))
        result
      }, .package = "RoBMA")
  }
})


test_that("joint scalar return routes do not calculate a discarded Gaussian reference", {

  S <- 5L
  K <- 2L
  predictive <- list(mu = matrix(.1, S, K), mu_extrapolated = matrix(.2, S, K),
    sei = c(.5, .8), tau_within = matrix(.1, S, K))
  context <- .zplot_density_reuse_context(rep(FALSE, S), predictive$sei)
  samples <- matrix(0, S, 1L)
  scalar <- list(fitted = matrix(.25, S, 3L), extrapolated = matrix(.5, S, 3L),
    weights = rep(2, S), EDR = NULL)
  testthat::local_mocked_bindings(
    .zplot_vector_selection_target = function(...) FALSE,
    .is_data_random = function(...) FALSE,
    .zplot_tau_samples = function(...) list(tau_within = matrix(.1, S, K),
      tau_between = matrix(.2, S, K)),
    .zplot_normal_density_matrix = function(...) stop("Discarded Gaussian reference requested."),
    .zplot_latent_mixture = function(...) scalar,
    .zplot_integrated_marginal = function(...) scalar,
    .package = "RoBMA")
  for (diagonal in c(TRUE, FALSE)) {
    sigma <- if (diagonal) diag(c(.25, .64)) else matrix(c(.25, .1, .1, .64), K)
    parts <- list(covariance = .as_block_covariance(array(rep(sigma, each = S), c(S, K, K))),
      context_covariance = if (diagonal) .as_block_covariance(array(
        rep(diag(.04, K), each = S), c(S, K, K))) else .block_covariance_zero(S, K))
    actual <- testthat::with_mocked_bindings(
      .zplot_joint_marginal(list(data = list()), samples, predictive, context,
        c(-1, 0, 1), FALSE, set_selection_likelihood_control()),
      .predict_joint_selection_gaussian_parts = function(...) parts,
      .package = "RoBMA")
    expect_identical(actual, scalar)
  }
})
