.density_cell_cache_inputs <- function() {

  # Unequal multiplicities and nonadjacent repeats exercise the original K
  # denominator. Equal SEs alone do not identify a cell: each parameter varies.
  columns <- c(1, 6, 1, 6, 2, 2, 1, 3, 3, 4, 4, 5, 5, 1)
  sei <- c(.7, .7, .7, .7, .7, 1.1)[columns]
  location <- c(.3, -.2, .3, .3, .3, -.4)[columns]
  integrated <- c(.8, .8, .9, .8, .15, .6)[columns]
  retained <- c(.1, .1, .1, .6, .6, 0)[columns]
  prior <- BayesTools::prior_weightfunction("one-sided", .5,
    BayesTools::wf_fixed(c(1, .2)))
  selection <- .selection_spec(list(outcome = list(bias = prior)),
    rep(0, length(sei)), sei, effect_direction = "positive", signed_data = FALSE)
  selection$omega <- rbind(c(1, .2), c(1, 0), c(1, .2))
  selection$alpha <- numeric(3L)
  selection$phack_kind <- integer(3L)
  selection$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_STEP, SELKERNEL_NORMAL)
  selection$vector_rule <- integer(3L)
  selection$use_normal <- c(FALSE, FALSE, TRUE)
  list(z_sequence = c(-3, -.3, 0, .1, .8, 2, 4),
    mean = outer(c(0, .1, -.2), location, "+"),
    sd = matrix(rep(integrated, each = 3L), 3L), sei = sei,
    selection_context = selection, extrapolate = c(FALSE, TRUE),
    latent_sd = matrix(rep(retained, each = 3L), 3L),
    quadrature = .gauss_hermite_nodes(15L))
}


.density_cell_uncached <- function(args) {

  # A one-column call cannot hit the repeated-cell memo. For this small grid
  # interpolation is disabled, so a cell has no full-row reference dependence.
  # Sum the returned cells in the original observation order, without grouping.
  cells <- lapply(seq_along(args$sei), function(column) {
    one <- args
    for (name in c("mean", "sd", "latent_sd")) {
      if (!is.null(one[[name]])) one[[name]] <- one[[name]][, column, drop = FALSE]
    }
    one$sei <- one$sei[column]
    do.call(.zplot_selnorm_density_matrix, one)
  })
  lapply(c("fitted", "extrapolated"), function(kind) {
    value <- cells[[1L]][[kind]]
    for (cell in seq_along(value)) {
      value[cell] <- sum(vapply(cells, function(x) x[[kind]][cell], numeric(1L))) /
        length(cells)
    }
    attr(value, "selection_transition_resolved") <- Reduce(`&`,
      lapply(cells, function(x) attr(x[[kind]], "selection_transition_resolved")))
    value
  }) |> stats::setNames(c("fitted", "extrapolated"))
}


test_that("density cell reuse preserves original-order uncached contributions", {

  args <- .density_cell_cache_inputs()
  # Native cells are accumulated in long double before the final conversion;
  # the oracle converts each cell first. Bound that extra rounding explicitly.
  rounding <- 64 * .Machine$double.eps
  for (direction in c(1L, -1L)) {
    args$selection_context$sign <- direction
    args$mean <- direction * abs(args$mean)
    reference <- .density_cell_uncached(args)
    actual <- do.call(.zplot_selnorm_density_matrix, args)
    for (kind in names(actual)) {
      expect_lte(max(abs(actual[[kind]] - reference[[kind]])),
        rounding * max(1, abs(reference[[kind]])))
      expect_identical(attr(actual[[kind]], "selection_transition_resolved"),
        attr(reference[[kind]], "selection_transition_resolved"))
      one <- args
      one$extrapolate <- identical(kind, "extrapolated")
      expect_identical(do.call(.zplot_selnorm_density_matrix, one), actual[[kind]])
    }
    coarse <- args
    coarse$quadrature <- .gauss_hermite_nodes(7L)
    coarse_reference <- .density_cell_uncached(coarse)
    attr(args$quadrature, "density_control") <- list(.001, coarse$quadrature)
    checked <- do.call(.zplot_selnorm_density_matrix, args)
    for (kind in names(checked)) {
      error <- apply(abs(reference[[kind]] - coarse_reference[[kind]]), 1L, max)
      peak <- apply(reference[[kind]], 1L, max)
      # The reported error is e/(peak-e). Compare its bounded equivalent
      # e/peak: zero-weight examples can have peak-e at rounding precision,
      # where converting individual cells first cannot recover that ratio.
      relative <- attr(checked[[kind]], "relative_integration_error")
      error_fraction <- ifelse(is.infinite(relative), 1, relative / (1 + relative))
      expect_equal(error_fraction, pmin(1, error / peak), tolerance = rounding)
      expect_identical(as.numeric(checked[[kind]]), as.numeric(actual[[kind]]))
      expect_identical(attr(checked[[kind]], "selection_transition_resolved"),
        attr(actual[[kind]], "selection_transition_resolved"))
    }
    attr(args$quadrature, "density_control") <- NULL
  }
})


test_that("density cell reuse keeps interpolation errors and row isolation", {

  args <- .density_cell_cache_inputs()
  args$z_sequence <- seq(-4, 4, length.out = 161L)
  args$quadrature <- .gauss_hermite_nodes(31L)
  direct <- do.call(.zplot_selnorm_density_matrix, args)
  attr(args$quadrature, "density_control") <- list(.001, .gauss_hermite_nodes(15L))
  interpolated <- do.call(.zplot_selnorm_density_matrix, args)
  for (kind in names(interpolated)) {
    difference <- apply(abs(interpolated[[kind]] - direct[[kind]]), 1L, max)
    peak <- apply(interpolated[[kind]], 1L, max)
    error <- attr(interpolated[[kind]], "relative_integration_error")
    expect_true(all(is.finite(error[c(1L, 3L)])))
    expect_true(all(difference <= peak * error + 64 * .Machine$double.eps))
  }
  for (row in seq_len(nrow(args$mean))) {
    one <- args
    for (name in c("mean", "sd", "latent_sd")) one[[name]] <- one[[name]][row, , drop = FALSE]
    one$selection_context <- BayesTools::selection_context_subset_rows(args$selection_context, row)
    actual <- do.call(.zplot_selnorm_density_matrix, one)
    for (kind in names(actual)) {
      expect_identical(as.numeric(actual[[kind]]), as.numeric(interpolated[[kind]][row, ]))
      for (attribute in c("relative_integration_error", "selection_transition_resolved")) {
        expect_identical(attr(actual[[kind]], attribute), attr(interpolated[[kind]], attribute)[row])
      }
    }
  }
})


test_that("density cell reuse preserves unresolved and invalid-input behavior", {

  args <- .density_cell_cache_inputs()
  args$mean[,] <- 0
  args$sd[,] <- .001
  args$latent_sd[,] <- .6
  args$z_sequence <- c(-1e-4, -5e-5)
  args$quadrature <- .gauss_hermite_nodes(3L)
  attr(args$quadrature, "density_control") <- list(.001, .gauss_hermite_nodes(1L))
  actual <- do.call(.zplot_selnorm_density_matrix, args)
  expect_identical(attr(actual$fitted, "selection_transition_resolved"), c(FALSE, TRUE, TRUE))
  expect_identical(attr(actual$extrapolated, "selection_transition_resolved"), c(FALSE, FALSE, TRUE))
  expect_identical(as.numeric(actual$fitted[2L, ]), c(0, 0))
  expect_identical(attr(actual$fitted, "relative_integration_error")[2L], 0)

  invalid <- args
  invalid$latent_sd[1L, 3L] <- -1
  expect_error(do.call(.zplot_selnorm_density_matrix, invalid),
    "positive integrated SDs and nonnegative retained SDs", fixed = TRUE)
  invalid <- args
  invalid$sd[1L, 3L] <- NaN
  expect_error(do.call(.zplot_selnorm_density_matrix, invalid),
    "positive integrated SDs and nonnegative retained SDs", fixed = TRUE)
  invalid <- args
  invalid$quadrature$nodes[1L] <- Inf
  expect_error(do.call(.zplot_selnorm_density_matrix, invalid),
    "Quadrature nodes and log weights must be finite", fixed = TRUE)
  # A normal row with an infinite location retains its unavailable density.
  # Nonfinite inputs must not become reusable cache entries.
  args$mean[3L, ] <- Inf
  actual <- do.call(.zplot_selnorm_density_matrix, args)
  expect_true(all(is.nan(actual$fitted[3L, ])))
  expect_identical(attr(actual$fitted, "relative_integration_error")[3L], Inf)
})


test_that("density cell results are unchanged by one or two native threads", {

  previous <- .native_threads_configure(1L)
  on.exit(.native_threads_configure(previous), add = TRUE)
  args <- .density_cell_cache_inputs()
  args$z_sequence <- seq(-4, 4, length.out = 161L)
  attr(args$quadrature, "density_control") <- list(.001, .gauss_hermite_nodes(7L))
  latent <- do.call(.zplot_selnorm_density_matrix, args)
  .native_threads_configure(2L)
  expect_identical(do.call(.zplot_selnorm_density_matrix, args), latent)

  # Nonlatent calls retain their original path and enough rows to schedule
  # two workers. They include normal, step and zero-weight posterior rows.
  args$latent_sd <- args$quadrature <- NULL
  rows <- rep(seq_len(3L), 8L)
  args$mean <- args$mean[rows, , drop = FALSE]
  args$sd <- args$sd[rows, , drop = FALSE]
  args$selection_context <- BayesTools::selection_context_subset_rows(args$selection_context, rows)
  for (extrapolate in list(FALSE, TRUE, c(FALSE, TRUE))) {
    args$extrapolate <- extrapolate
    .native_threads_configure(1L)
    reference <- do.call(.zplot_selnorm_density_matrix, args)
    .native_threads_configure(2L)
    expect_identical(do.call(.zplot_selnorm_density_matrix, args), reference)
  }
})
