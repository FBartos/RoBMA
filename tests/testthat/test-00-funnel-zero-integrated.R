context("(00) funnel contours of cells with no integrated variance")

# A conditioned source keeps its population law outside the selection
# normalization. A (posterior row, standard error) cell with no integrated
# variance has nothing left to select among: the selected normal at the
# candidate location m is a point mass at m with weight w(m) / w(m), so under a
# step weight function with every weight positive the observed law is the
# retained Normal itself, at the conditioned spread. The limit is continuous in
# the integrated SD.


# A funnel setup of `S` posterior rows under the four-bin step weight function
# of .test_step_spec(), with the variance buckets of the selection sources.
# `omega` replaces the rows' weights.
.funnel_zero_setup <- function(mu, integrated, conditioned, retains_sampling,
                               omega = NULL, kernel_mode = SELKERNEL_STEP) {

  S    <- length(mu)
  spec <- .test_step_spec(c(.2, .1), c(.2, .2))
  selection <- spec
  if (is.null(omega)) {
    omega <- spec[["fixed_omega"]]
  }
  selection[["omega"]]       <- matrix(omega, nrow = S, ncol = length(omega), byrow = TRUE)
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- rep(as.integer(kernel_mode), S)
  selection[["use_normal"]]  <- rep(FALSE, S)
  selection[["vector_rule"]] <- rep(0L, S)

  list(
    mu = mu, tau = rep(0, S), PET = 0, PEESE = 0,
    is_weightfunction = rep(TRUE, S), selection = selection,
    weights = rep(1, S),
    sources = list(integrated = integrated, conditioned = conditioned,
                   retains_sampling = retains_sampling, common = TRUE)
  )
}


# The mixture over the rows of their observed laws: a cancelled row is the
# Normal at its conditioned spread, any other row the mixture over the
# conditioned offset of the selected normals at its integrated SD, each with its
# own normalizer (301-node Gauss-Hermite rule). Inverted by root finding on the
# mixture CDF.
.funnel_zero_reference <- function(setup, se, cancelled, probs = c(.025, .975)) {

  static  <- BayesTools::selection_native_static_args(setup[["selection"]])
  omega   <- as.numeric(setup[["selection"]][["omega"]][1L, ])
  rule    <- .gauss_hermite_nodes(301L)
  nodes   <- rule[["nodes"]]
  node_weights <- rule[["weights"]] / sum(rule[["weights"]])
  S <- length(setup[["mu"]])
  integrated  <- sqrt(setup[["sources"]][["integrated"]])
  conditioned <- sqrt(setup[["sources"]][["conditioned"]] + se^2)
  selected_cdf <- function(q, centre, sd_kernel) {
    lower <- static[["z_lower"]] * se
    upper <- static[["z_upper"]] * se
    inside <- stats::pnorm((pmin(pmax(q, lower), upper) - centre) / sd_kernel) -
      stats::pnorm((lower - centre) / sd_kernel)
    total <- stats::pnorm((upper - centre) / sd_kernel) -
      stats::pnorm((lower - centre) / sd_kernel)
    sum(omega * inside) / sum(omega * total)
  }
  mixture <- function(q) {
    mean(vapply(seq_len(S), function(row) {
      if (cancelled[[row]]) {
        return(stats::pnorm(q, setup[["mu"]][[row]], conditioned[[row]]))
      }
      sum(node_weights * vapply(nodes, function(node) {
        selected_cdf(q, setup[["mu"]][[row]] + node * conditioned[[row]],
                     integrated[[row]])
      }, numeric(1)))
    }, numeric(1)))
  }
  vapply(probs, function(p) {
    stats::uniroot(function(q) mixture(q) - p, interval = c(-5, 5),
                   tol = 1e-12)[["root"]]
  }, numeric(1))
}


test_that("a cell with no integrated variance is the retained Normal law", {

  skip_if_not(.has_native_selnorm_kernel())

  # The seed: one row, all variance in the retained sampling error. The band is
  # the Normal at the sampling error, +-qnorm(.975) standard errors, not an atom
  # of a 401-node quadrature rule at 1.88.
  setup <- .funnel_zero_setup(0, integrated = 0, conditioned = 0,
                              retains_sampling = TRUE)
  se    <- c(.05, .1, .2)
  band  <- .funnel_model_averaged_quantiles_native(se, setup, "positive")
  expect_equal(band[["upper"]] / se, rep(stats::qnorm(.975), 3L), tolerance = 1e-9)
  expect_equal(band[["lower"]] / se, rep(stats::qnorm(.025), 3L), tolerance = 1e-9)
  expect_equal(2 * stats::pnorm(band[["upper"]] / se) - 1, rep(.95, 3L),
               tolerance = 1e-9)
  expect_equal(band[["mid"]], rep(0, 3L), tolerance = 1e-12)
})


test_that("rows with and without integrated variance mix in one contour", {

  skip_if_not(.has_native_selnorm_kernel())

  # Half the rows have no integrated variance (the retained sampling error is
  # their only spread) and half have integrated variance beside it.
  mu         <- c(0, .1, .2, -.1)
  integrated <- c(0, .01, 0, .04)
  setup <- .funnel_zero_setup(mu, integrated = integrated, conditioned = rep(0, 4L),
                              retains_sampling = TRUE)
  cancelled <- integrated == 0
  for (se in c(.1, .25)) {
    band <- .funnel_model_averaged_quantiles_native(se, setup, "positive")
    expected <- .funnel_zero_reference(setup, se, cancelled)
    expect_equal(c(band[["lower"]], band[["upper"]]), expected, tolerance = 1e-5,
                 info = paste("se", se))
  }

  # The quadrature rule is chosen from the rows that still need it: the cancelled
  # row's zero integrated SD no longer sets the finest rule.
  integrated_sd  <- cbind(c(1, 0), c(1, 0))
  conditioned_sd <- cbind(c(.1, 1), c(.1, 1))
  two      <- .funnel_zero_setup(c(0, 0), integrated = c(1, 0), conditioned = c(0, 0),
                                 retains_sampling = TRUE)
  location <- matrix(0, 2L, 2L)
  reduced  <- .funnel_conditioned_expansion(location, integrated_sd, conditioned_sd,
                                            two)
  # The selection-cancellation certificate is resolved before budgeting.
  expect_identical(nrow(reduced[["mean"]]), .funnel_conditioned_order(
    integrated_sd[1L, , drop = FALSE], conditioned_sd[1L, , drop = FALSE]) + 1L)
  expect_identical(reduced[["selected_rows"]],
                   c(rep(TRUE, nrow(reduced[["mean"]]) - 1L), FALSE))
  expect_equal(sum(reduced[["weights"]]), 1)
})

test_that("constant and ordinary funnel rows retain exact Gaussian components", {

  setup <- .funnel_zero_setup(c(0, .1), integrated = c(1e-12, .04),
    conditioned = c(0, 0), retains_sampling = TRUE, omega = c(1, 1, 1, 1))
  setup$is_weightfunction <- c(TRUE, FALSE)
  setup$weights <- c(.3, .7)
  for (v in c(0, 1e-12, 1e-10, 1e-8, 1e-6)) {
    setup$sources$integrated[1L] <- v
    band <- .funnel_model_averaged_quantiles_native(.1, setup, "positive")
    sd <- sqrt(c(v + .01, .05))
    expected <- vapply(c(.025, .975, .5), function(p) stats::uniroot(
      function(q) sum(c(.3, .7) * stats::pnorm(q, c(0, .1), sd)) - p,
      c(-3, 3), tol = 1e-12)$root, numeric(1))
    expect_lt(max(abs(unlist(band) - expected)), 1e-11)
    expansion <- .funnel_conditioned_expansion(matrix(c(0, .1)),
      matrix(sqrt(c(v, .04))), matrix(.1, 2L, 1L), setup)
    expect_identical(nrow(expansion$mean), 2L)
    expect_identical(expansion$selected_rows, c(FALSE, FALSE))
    expect_identical(expansion$weights, c(.3, .7))
  }
})

test_that("funnel quadrature refuses an unresolved required order before expansion", {

  expect_identical(vapply(c(1, 6, 8), function(ratio) {
    .funnel_conditioned_order(1, ratio)
  }, integer(1)), c(49L, 289L, 385L))
  for (ratio in c(9, 1e5, Inf)) {
    expect_error(.funnel_conditioned_order(1, ratio),
      "Selected funnel contours were rejected by diagnostics: the conditioned-to-integrated standard deviation ratio",
      fixed = TRUE)
  }
  condition <- tryCatch(.funnel_conditioned_order(1, 9), error = identity)
  expect_s3_class(condition, "RoBMA_funnel_quadrature_rejected")
  expect_identical(class(condition), c("RoBMA_funnel_quadrature_rejected",
    "RoBMA_numerical_unavailable", "error", "condition"))
  expect_identical(conditionMessage(condition), paste0(
    "Selected funnel contours were rejected by diagnostics: the conditioned-to-integrated ",
    "standard deviation ratio of 9 requires a quadrature order of 433 that the funnel ",
    "quadrature cannot resolve. Use 'sampling_bias = FALSE' for an unselected reference."
  ))
  expect_null(conditionCall(condition))
  expect_identical(condition$ratio, 9)
  expect_identical(condition$required_order, 433)
  expect_identical(condition$maximum, 401L)
  expect_identical(.funnel_conditioned_order(0, 0), 21L)
  setup <- .funnel_zero_setup(0, integrated = 1, conditioned = 81,
    retains_sampling = FALSE)
  testthat::local_mocked_bindings(
    .gauss_hermite_nodes = function(...) stop("Node expansion must not start."),
    .package = "RoBMA")
  expect_error(.funnel_conditioned_expansion(matrix(0), matrix(1), matrix(9), setup),
    "ratio of 9 requires a quadrature order of 433", fixed = TRUE)
  setup$sources$integrated <- c(1, 2)
  expect_error(.funnel_model_averaged_quantiles_native(.1, setup, "positive"),
    "Funnel source variances must contain one finite non-negative value per posterior row.", fixed = TRUE)
})


test_that("a cell whose selection does not cancel keeps the selected route", {

  skip_if_not(.has_native_selnorm_kernel())

  base <- .funnel_zero_setup(c(0, .1), integrated = c(0, 0), conditioned = c(0, 0),
                             retains_sampling = TRUE)
  integrated <- matrix(0, 2L, 1L)
  expect_identical(.funnel_selection_cancelled_cells(integrated, base),
                   matrix(TRUE, 2L, 1L))

  # A zero bin: a candidate in it is never observed.
  zero_bin <- .funnel_zero_setup(c(0, .1), integrated = c(0, 0),
                                 conditioned = c(0, 0), retains_sampling = TRUE,
                                 omega = c(1, .7, 0, .2))
  expect_identical(.funnel_selection_cancelled_cells(integrated, zero_bin),
                   matrix(FALSE, 2L, 1L))
  # A p-hacking kernel.
  phack <- .funnel_zero_setup(c(0, .1), integrated = c(0, 0), conditioned = c(0, 0),
                              retains_sampling = TRUE,
                              kernel_mode = SELKERNEL_STEP_PHACK_POWER)
  expect_identical(.funnel_selection_cancelled_cells(integrated, phack),
                   matrix(FALSE, 2L, 1L))
  # A positive integrated SD, and rows on a non-selection branch.
  expect_identical(.funnel_selection_cancelled_cells(matrix(c(0, .1), 2L, 1L), base),
                   matrix(c(TRUE, FALSE), 2L, 1L))
  plain <- base
  plain[["is_weightfunction"]] <- c(FALSE, TRUE)
  expect_identical(.funnel_selection_cancelled_cells(integrated, plain),
                   matrix(c(FALSE, TRUE), 2L, 1L))

  # The zero-bin route has no certified degenerate selected law.
  location <- matrix(0, 2L, 1L)
  error <- tryCatch(.funnel_conditioned_expansion(location, integrated,
    matrix(.1, 2L, 1L), zero_bin), error = function(e) e)
  expect_identical(conditionMessage(error), paste0(
    "Selected funnel contours are unavailable for exact-zero integrated ",
    "cells without a valid selection-cancellation support certificate. ",
    "Use 'sampling_bias = FALSE' for an unselected reference."))
  expect_null(conditionCall(error))
})


test_that("the columns of a contour are inverted in groups of the same cancelling rows", {

  skip_if_not(.has_native_selnorm_kernel())

  # Without retained sampling the integrated SD of a cell is the root of the
  # integrated variance and se^2: it is zero at se = 0 alone. That column
  # cancels and the others do not, and each is the same contour as the column
  # inverted by itself.
  setup <- .funnel_zero_setup(c(0, .1), integrated = c(0, 0), conditioned = c(.01, .02),
                              retains_sampling = FALSE)
  se <- c(0, .1, .2)
  together <- .funnel_model_averaged_quantiles_native(se, setup, "positive")
  for (column in seq_along(se)) {
    alone <- .funnel_model_averaged_quantiles_native(se[[column]], setup, "positive")
    expect_equal(c(together[["lower"]][[column]], together[["upper"]][[column]]),
                 c(alone[["lower"]], alone[["upper"]]), tolerance = 1e-9)
  }
  # At se = 0 the band is the mixture of the two retained Normals.
  mixture <- function(q) mean(c(stats::pnorm(q, 0, .1), stats::pnorm(q, .1, sqrt(.02))))
  expected <- vapply(c(.025, .975), function(p) {
    stats::uniroot(function(q) mixture(q) - p, c(-2, 2), tol = 1e-12)$root
  }, numeric(1))
  expect_equal(c(together[["lower"]][[1L]], together[["upper"]][[1L]]), expected,
               tolerance = 1e-8)
})


test_that("each expanded draw gets the quadrature order of its own ratio", {

  # Conditioned-to-integrated ratios of .5, 1.5 and 2.2 (the largest over the two
  # columns of a draw) need 49, 73 and 107 nodes. The rows hold the nodes and
  # weights of each draw's own rule, not the rule of the most extreme draw.
  setup <- .funnel_zero_setup(c(0, .1, -.1), integrated = rep(1, 3), conditioned = rep(0, 3),
    retains_sampling = FALSE)
  integrated_sd  <- matrix(1, 3L, 2L)
  conditioned_sd <- cbind(c(.5, 1.5, 2.2), c(.2, 1, 2))
  location       <- cbind(c(0, .1, -.1), c(.05, .15, -.05))
  orders <- vapply(1:3, function(draw) .funnel_conditioned_order(
    integrated_sd[draw, , drop = FALSE], conditioned_sd[draw, , drop = FALSE]), integer(1L))
  expect_identical(orders, c(49L, 73L, 107L))

  expansion <- .funnel_conditioned_expansion(location, integrated_sd, conditioned_sd, setup)
  expect_identical(nrow(expansion[["mean"]]), sum(orders))
  expect_identical(nrow(expansion[["sd"]]), sum(orders))
  expect_equal(sum(expansion[["weights"]]), 1)
  expect_identical(expansion[["selected_rows"]], rep(TRUE, sum(orders)))

  rows <- do.call(rbind, lapply(1:3, function(draw) {
    rule <- .gauss_hermite_nodes(orders[[draw]])
    cbind(
      mean   = outer(rule[["nodes"]], conditioned_sd[draw, ]) +
        matrix(location[draw, ], orders[[draw]], 2L, byrow = TRUE),
      sd     = matrix(integrated_sd[draw, ], orders[[draw]], 2L, byrow = TRUE),
      weight = rule[["weights"]] / sum(rule[["weights"]]) / 3
    )
  }))
  observed <- cbind(expansion[["mean"]], expansion[["sd"]], expansion[["weights"]])
  sort_rows <- function(x) unname(x[do.call(order, as.data.frame(x)), ])
  expect_equal(sort_rows(observed), sort_rows(rows), tolerance = 1e-14)

  # Draws that share an order share one rule: a single group is the plain rule.
  same <- .funnel_conditioned_expansion(location[1L, , drop = FALSE],
    integrated_sd[1L, , drop = FALSE], conditioned_sd[1L, , drop = FALSE],
    .funnel_zero_setup(0, integrated = 1, conditioned = 0, retains_sampling = FALSE))
  expect_identical(nrow(same[["mean"]]), 49L)
})


test_that("per-draw quadrature orders keep the band within the order-401 reference", {

  skip_if_not(.has_native_selnorm_kernel())

  # Six draws with ratios from .3 to 2.5 need orders from 49 to 121. The
  # reference gives every draw the largest order the rule allows (401).
  se <- c(.1, .2)
  ratio <- c(.3, .6, 1, 1.5, 2, 2.5)
  integrated  <- rep(.03, 6L)
  conditioned <- (ratio * sqrt(integrated + se[[1L]]^2))^2
  setup <- .funnel_zero_setup(seq(-.2, .3, length.out = 6L), integrated = integrated,
    conditioned = conditioned, retains_sampling = FALSE)

  calls <- new.env()
  calls$rows <- integer()
  original <- .plot_mixture_quantiles_native
  testthat::local_mocked_bindings(
    .plot_mixture_quantiles_native = function(mean_samples, ...) {
      calls$rows <- c(calls$rows, nrow(mean_samples))
      original(mean_samples = mean_samples, ...)
    },
    .package = "RoBMA"
  )
  band <- .funnel_model_averaged_quantiles_native(se[[1L]], setup, "positive")
  sd_integrated  <- rep(sqrt(.03 + se[[1L]]^2), 6L)
  sd_conditioned <- ratio * sd_integrated
  expect_identical(calls$rows, sum(vapply(seq_len(6L), function(draw) {
    .funnel_conditioned_order(sd_integrated[draw], sd_conditioned[draw])
  }, integer(1L))))
  expect_identical(calls$rows, 49L * 3L + 73L + 97L + 121L)

  reference <- testthat::with_mocked_bindings(
    .funnel_model_averaged_quantiles_native(se[[1L]], setup, "positive"),
    .funnel_conditioned_order = function(...) 401L,
    .package = "RoBMA"
  )
  expect_lt(max(abs(unlist(band) - unlist(reference))), 4e-8)
})


test_that("a draw that no quadrature order resolves refuses the whole expansion", {

  # One draw of ratio 9 beside draws of ratio .5: the refusal names the largest
  # ratio and no node is built, as when every draw is expanded with one order.
  setup <- .funnel_zero_setup(c(0, .1), integrated = c(1, 1), conditioned = c(0, 0),
    retains_sampling = FALSE)
  testthat::local_mocked_bindings(
    .gauss_hermite_nodes = function(...) stop("Node expansion must not start."),
    .package = "RoBMA")
  condition <- tryCatch(.funnel_conditioned_expansion(matrix(c(0, .1)),
    matrix(1, 2L, 1L), matrix(c(.5, 9), 2L, 1L), setup), error = identity)
  expect_s3_class(condition, "RoBMA_funnel_quadrature_rejected")
  expect_identical(condition$ratio, 9)
  expect_identical(condition$required_order, 433)
})
