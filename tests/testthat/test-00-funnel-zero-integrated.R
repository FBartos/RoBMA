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
  finest   <- .funnel_conditioned_expansion(location, integrated_sd, conditioned_sd,
                                            two)
  reduced  <- .funnel_conditioned_expansion(location, integrated_sd, conditioned_sd,
                                            two, cancelled = c(FALSE, TRUE))
  # The existing certificate is resolved before budgeting even when the
  # caller omits the optional precomputed cancellation flags.
  expect_identical(nrow(finest[["mean"]]), nrow(reduced[["mean"]]))
  expect_identical(nrow(reduced[["mean"]]), .funnel_conditioned_order(
    integrated_sd[1L, , drop = FALSE], conditioned_sd[1L, , drop = FALSE]) + 1L)
  expect_identical(reduced[["selected_rows"]],
                   c(rep(TRUE, nrow(reduced[["mean"]]) - 1L), FALSE))
  expect_equal(sum(reduced[["weights"]]), sum(finest[["weights"]]))
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
      "Selected funnel contours are numerically unavailable: required quadrature order", fixed = TRUE)
  }
  expect_identical(.funnel_conditioned_order(0, 0), 21L)
  setup <- .funnel_zero_setup(0, integrated = 1, conditioned = 81,
    retains_sampling = FALSE)
  testthat::local_mocked_bindings(
    .gauss_hermite_nodes = function(...) stop("Node expansion must not start."),
    .package = "RoBMA")
  expect_error(.funnel_conditioned_expansion(matrix(0), matrix(1), matrix(9), setup),
    "required quadrature order 433 exceeds the current cap 401", fixed = TRUE)
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
    matrix(.1, 2L, 1L), zero_bin,
    cancelled = .funnel_selection_cancelled_cells(integrated, zero_bin)[, 1L]),
    error = function(e) e)
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
