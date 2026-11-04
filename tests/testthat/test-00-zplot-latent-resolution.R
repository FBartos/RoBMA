.zplot_latent_resolution_reference <- function(nodes, mean, sd, latent_sd, sei, selection) {

  vapply(seq_len(nrow(mean)), function(row) {
    .selection_transition_rule_resolves(
      nodes, mean[row, ], sd[row, ], matrix(latent_sd[row, ], ncol = 1L), sei,
      list(kernel_mode = selection$kernel_mode[row],
        omega = selection$omega[row, , drop = FALSE],
        z_lower = selection$z_lower, sign = selection$sign)
    )
  }, logical(1L))
}


test_that("latent transition checks preserve every row's original multi-observation decision", {

  withr::local_seed(183L)
  S <- 40L
  K <- 3L
  sei <- c(.5, .25, 1)
  mean <- matrix(stats::rnorm(S * K, 0, 2), S, K)
  sd <- matrix(exp(stats::runif(S * K, -2, 1)), S, K)
  latent <- matrix(stats::runif(S * K, .1, 2), S, K)
  selection <- list(
    kernel_mode = rep(SELKERNEL_STEP, S),
    omega = cbind(1, matrix(sample(c(.2, .5, 1), S * 2L, replace = TRUE), S)),
    z_lower = c(2, .5, -Inf), sign = 1L
  )
  selection$omega[1:6, ] <- rep(c(1, .2, .2), each = 6L)
  # Exact node ties use the wider adjacent gap; exterior centers are left to
  # the existing tail diagnostics. Normal and flat-weight rows bypass nodes.
  latent[1:6, ] <- 1
  sd[1:6, ] <- .75
  mean[1, ] <- 2 * sei - c(-1, 0, 0)
  mean[2, ] <- 2 * sei - c(-2, 2, 3)
  mean[3, ] <- mean[1, ]
  sd[3, ] <- 1
  selection$kernel_mode[4] <- SELKERNEL_NORMAL
  selection$omega[5, ] <- 1
  mean[4:6, ] <- Inf
  sd[4:6, ] <- NA_real_
  latent[6, ] <- 0
  latent[7, 1] <- Inf
  sd[8, 2] <- 0
  mean[9, 3] <- NA_real_
  grid <- c(-2, -1, -.5, 0, .5, 1, 2)
  node_sets <- list(grid, rev(grid), .gauss_hermite_nodes(15L)$nodes,
    c(grid, NA_real_), c(grid, grid[1L]), c(grid, Inf), c(-Inf, grid), 0, numeric())
  seen <- logical()
  for (direction in c(1L, -1L)) {
    selection$sign <- direction
    for (nodes in node_sets) {
      reference <- .zplot_latent_resolution_reference(
        nodes, direction * mean, sd, latent, sei, selection
      )
      actual <- .zplot_latent_transition_resolves(
        nodes, direction * mean, sd, latent, sei, selection
      )
      expect_identical(actual, reference)
      if (identical(nodes, grid)) {
        expect_identical(actual[1:6], c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
      }
      seen <- c(seen, reference)
    }
  }
  expect_true(any(seen) && any(!seen))
})


test_that("latent summaries retain the original resolution gate and reuse identical means", {

  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)))
  sei <- c(.2, .3)
  context <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), sei,
    effect_direction = "positive", signed_data = FALSE)
  context$omega <- matrix(rep(c(1, .2), 3L), 3L, byrow = TRUE)
  context$alpha <- numeric(3L)
  context$phack_kind <- integer(3L)
  context$kernel_mode <- rep(SELKERNEL_STEP, 3L)
  context$vector_rule <- integer(3L)
  context$use_normal <- rep(FALSE, 3L)
  mean <- matrix(c(.1, -.1, .3, -.2, .3, .2), 3L)
  sd <- matrix(c(.4, .5, .3, .7, .9, .6), 3L)
  latent <- matrix(c(.8, 0, .2, .3, .1, .5), 3L)
  control <- set_selection_likelihood_control(relative_tolerance = 1e-7)
  evaluate <- function(extrapolated) {
    .zplot_latent_mixture(2, mean, sd, latent, sei, context, TRUE, control,
      mean_extrapolated = extrapolated)
  }
  original <- .zplot_latent_transition_resolves
  for (extrapolated in list(mean, mean + .125)) {
    calls <- list()
    actual <- testthat::with_mocked_bindings(evaluate(extrapolated),
      .zplot_latent_transition_resolves = function(nodes, mean, ...) {
        calls[[length(calls) + 1L]] <<- list(nodes = nodes, mean = mean)
        original(nodes, mean, ...)
      }, .package = "RoBMA")
    reference <- testthat::with_mocked_bindings(evaluate(extrapolated),
      .zplot_latent_transition_resolves = .zplot_latent_resolution_reference,
      .package = "RoBMA")
    expect_identical(actual, reference)
    same_mean <- identical(mean, extrapolated)
    # The finest precheck occurs before posterior rows can be retired.
    expect_identical(calls[[1L]]$mean, mean)
    if (same_mean) {
      expect_false(identical(calls[[1L]]$nodes, calls[[2L]]$nodes))
    } else {
      expect_identical(calls[[2L]]$mean, extrapolated)
      expect_identical(calls[[1L]]$nodes, calls[[2L]]$nodes)
    }
  }
})
