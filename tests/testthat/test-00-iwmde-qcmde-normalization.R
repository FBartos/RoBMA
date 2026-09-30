context("qCMDE normalization pass: call structure and grid requests")

# The joint density is replaced by analytic rows so that the calls the pass
# makes, and what it assembles from them, can be compared exactly.

.qcmde_pass_rows <- function(values) {

  cbind(stats::dnorm(values, log = TRUE),
        stats::dnorm(values, .5, 1.2, log = TRUE))
}


# Laplace rows are estimate rows whose tails decay slowly, so a normalization
# range of [-1, 1] has to be extended by many steps.
.qcmde_pass_laplace_rows <- function(values) {

  cbind(-abs(values - .2) / .5, -abs(values + .1) / .8)
}


# One normalization pass over `rows` with the joint density and the call-owned
# grid builder mocked. `grid` makes the builder return a grid; `fail_in_grid`
# names a predicate on the values of a request of that grid that makes it fail.
.qcmde_pass_run <- function(rows, range, n_points, display_grid = c(-.33, 1.97),
                            estimate_only = TRUE, grid = FALSE,
                            fail_in_grid = function(values) FALSE,
                            .env = parent.frame()) {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  values    <- seq(range[[1L]], range[[2L]], length.out = n_points)
  calls     <- list()
  requested <- new.env()
  requested$queries <- NULL
  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      in_grid <- !is.null(context[["normalizer_grid"]])
      calls[[length(calls) + 1L]] <<- list(values = values, grid = in_grid)
      if (in_grid && fail_in_grid(values)) {
        stop("the grid cannot answer this request")
      }
      out <- rows(values)
      attr(out, "max_quadrature_relative_change") <- .001
      out
    },
    .iwmde_qcmde_estimate_only = function(context, row_states, rows) {
      estimate_only
    },
    .selection_normalizer_grid = function(context, queries, rows) {
      if (!grid) {
        return(NULL)
      }
      requested$queries <- sort(unique(queries[is.finite(queries)]))
      new.env()
    },
    .selection_normalizer_grid_diagnostics = function(shared) NULL,
    .package = "RoBMA",
    .env     = .env
  )
  pass <- .iwmde_qcmde_normalization_pass(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = list(x = values, z = values,
                              log_jacobian = rep(0, length(values))),
    transform          = transform,
    normalization_prob = .999,
    row_states         = list(list(), list()),
    replacement        = list(type = "scalar"),
    estimator_rows     = 1:2
  )

  list(pass = pass, calls = calls, queries = requested$queries)
}


# The pass without the attributes of the display evaluation, for exact
# comparison of everything it assembled.
.qcmde_pass_plain <- function(pass) {

  display <- pass[["log_q_display"]]
  pass[["log_q_display"]] <- matrix(as.numeric(display), nrow(display),
                                    ncol(display))

  pass
}


test_that("an estimate-only batch evaluates the display values, nodes and midpoints in one call", {

  merged     <- .qcmde_pass_run(.qcmde_pass_rows, c(-6, 6), 25L)
  sequential <- .qcmde_pass_run(.qcmde_pass_rows, c(-6, 6), 25L,
                                estimate_only = FALSE)

  # One call with the display values, the 25 nodes and the 24 midpoints between
  # them; the sequential path makes a call for each.
  expect_length(sequential[["calls"]], 3L)
  expect_identical(vapply(sequential[["calls"]], function(call) {
    length(call[["values"]])
  }, integer(1)), c(2L, 25L, 24L))
  expect_length(merged[["calls"]], 1L)
  expect_identical(
    merged[["calls"]][[1L]][["values"]],
    unlist(lapply(sequential[["calls"]], `[[`, "values"), use.names = FALSE)
  )

  # The grids, normalizers and diagnostics assembled are those of the separate
  # calls, value for value.
  expect_identical(.qcmde_pass_plain(merged[["pass"]]),
                   .qcmde_pass_plain(sequential[["pass"]]))
  expect_identical(merged[["pass"]][["quadrature_change"]], .001)
})


test_that("the extension of an estimate-only batch follows the joint call and assembles the same grids", {

  merged     <- .qcmde_pass_run(.qcmde_pass_laplace_rows, c(-1, 1), 21L)
  sequential <- .qcmde_pass_run(.qcmde_pass_laplace_rows, c(-1, 1), 21L,
                                estimate_only = FALSE)
  steps <- merged[["pass"]][["extension_steps"]]
  expect_gt(steps[["lower"]], 0L)
  expect_gt(steps[["upper"]], 0L)

  expect_identical(.qcmde_pass_plain(merged[["pass"]]),
                   .qcmde_pass_plain(sequential[["pass"]]))
  # The first call holds the display values, the nodes and the initial
  # midpoints; the extension nodes and midpoints are evaluated after it, and
  # no value twice.
  expect_length(merged[["calls"]][[1L]][["values"]], 2L + 21L + 20L)
  values <- unlist(lapply(merged[["calls"]], `[[`, "values"), use.names = FALSE)
  expect_length(unique(values), length(values))
  expect_true(all(merged[["pass"]][["nested"]][["x"]] %in% values))
  expect_lt(length(merged[["calls"]]), length(sequential[["calls"]]))
})


test_that("call-owned grids keep the display values alone in their first request and reserve the extension steps", {

  reserve <- .iwmde_qcmde_extension_reserve(41L)
  expect_identical(reserve, 5L)
  display      <- c(-.33, 1.97)
  with_grid    <- .qcmde_pass_run(.qcmde_pass_laplace_rows, c(-1, 1), 41L,
                                  display_grid = display, grid = TRUE)
  without_grid <- .qcmde_pass_run(.qcmde_pass_laplace_rows, c(-1, 1), 41L,
                                  display_grid = display)

  # The queries hold the display values, the nodes and midpoints of the range,
  # and the reserve of extension nodes and midpoints on each side.
  lattice  <- .iwmde_qcmde_lattice(c(-1, 1), 41L,
                                   .iwmde_parameter_transform(c(-Inf, Inf)))
  reserved <- .iwmde_qcmde_lattice_points(
    lattice, seq(-reserve, 40L + reserve, by = .5)
  )[["x"]]
  expect_identical(with_grid[["queries"]], sort(unique(c(display, reserved))))

  # The display values are alone in the first request. The nodes, midpoints
  # and the extension values of the reserve are requests of the grid; the
  # extension beyond the reserve is evaluated without it.
  calls <- with_grid[["calls"]]
  expect_identical(calls[[1L]][["values"]], display)
  expect_true(calls[[1L]][["grid"]])
  by_grid <- vapply(calls, `[[`, logical(1), "grid")
  in_grid <- unlist(lapply(calls[by_grid], `[[`, "values"), use.names = FALSE)
  outside <- unlist(lapply(calls[!by_grid], `[[`, "values"), use.names = FALSE)
  expect_true(all(in_grid %in% c(display, reserved)))
  expect_false(any(outside %in% c(display, reserved)))
  expect_gt(length(outside), 0L)
  extension <- in_grid[in_grid < -1 - 1e-9 | in_grid > 1 + 1e-9]
  extension <- setdiff(extension, display)
  expect_gt(length(extension), 0L)
  expect_lte(sum(extension < -1), 2L * reserve)
  expect_lte(sum(extension > 1), 2L * reserve)

  # The grids assembled from the requests are those of the run without grids.
  expect_identical(.qcmde_pass_plain(with_grid[["pass"]]),
                   .qcmde_pass_plain(without_grid[["pass"]]))
})


test_that("an extension value the grid cannot answer is evaluated without it", {

  inside   <- function(values) any(values < -1.001 | values > 1.001)
  display  <- c(-.3, .4)
  reference <- .qcmde_pass_run(.qcmde_pass_laplace_rows, c(-1, 1), 41L,
                               display_grid = display)
  failing   <- .qcmde_pass_run(.qcmde_pass_laplace_rows, c(-1, 1), 41L,
                               display_grid = display, grid = TRUE,
                               fail_in_grid = inside)

  # The requests of the grid for values beyond the range failed, and the pass
  # evaluated those values directly; what it assembled is unchanged.
  calls <- failing[["calls"]]
  by_grid <- vapply(calls, `[[`, logical(1), "grid")
  failed  <- vapply(calls[by_grid], function(call) inside(call[["values"]]),
                    logical(1))
  expect_true(any(failed))
  expect_identical(.qcmde_pass_plain(failing[["pass"]]),
                   .qcmde_pass_plain(reference[["pass"]]))
})


test_that("only batches whose rows cannot report a Gaussian kernel are estimate-only", {

  state <- function(weightfunction) {
    list(active_setup = list(is_weightfunction = weightfunction))
  }
  data <- function(type) structure(list(), outcome_type = type)
  rows <- 1:2

  expect_true(.iwmde_qcmde_estimate_only(
    list(data = data("bin")), list(state(NULL), state(NULL)), rows))
  expect_true(.iwmde_qcmde_estimate_only(
    list(data = data("pois")), list(state(NULL), state(NULL)), rows))
  expect_true(.iwmde_qcmde_estimate_only(
    list(data = data("norm")), list(state(TRUE), state(TRUE)), rows))
  # A row without a weight function may report the kernel of its
  # normal-location route.
  expect_false(.iwmde_qcmde_estimate_only(
    list(data = data("norm")), list(state(TRUE), state(FALSE)), rows))
  expect_false(.iwmde_qcmde_estimate_only(
    list(data = data("norm")), list(state(TRUE), list()), rows))
  # Unknown models and empty batches report nothing certain.
  expect_false(.iwmde_qcmde_estimate_only(list(), list(state(TRUE)), 1L))
  expect_false(.iwmde_qcmde_estimate_only(
    list(data = data("norm")), list(state(TRUE)), integer()))
  expect_false(.iwmde_qcmde_estimate_only(
    list(data = data("other")), list(state(TRUE)), 1L))
})
