# ============================================================================ #
# test-00-selection-normalizer-grid-vectorized.R
# ============================================================================ #

context("Selection normalizer grid vectorized helpers")
skip_on_cran()


# The per-draw and per-packed-column loops these helpers replaced. Both are
# transcribed from the implementations they superseded, so the comparison is a
# value-identity check and not a second derivation of the same quantity.
.reference_grid_gaussian_lpdf <- function(yi, means, covariance_lower, block_size) {

  S        <- nrow(means)
  out      <- numeric(S)
  lower    <- lower.tri(matrix(0, block_size, block_size), diag = TRUE)
  upper    <- upper.tri(matrix(0, block_size, block_size))
  previous <- NULL
  root     <- NULL
  log_det  <- NA_real_
  for (draw in seq_len(S)) {
    packed <- covariance_lower[draw, ]
    if (is.null(previous) || any(packed != previous)) {
      covariance        <- matrix(0, block_size, block_size)
      covariance[lower] <- packed
      covariance[upper] <- t(covariance)[upper]
      root <- tryCatch(chol(covariance), error = function(e) NULL)
      if (is.null(root)) {
        return(NULL)
      }
      log_det  <- 2 * sum(log(diag(root)))
      previous <- packed
    }
    whitened    <- backsolve(root, yi - means[draw, ], transpose = TRUE)
    out[[draw]] <- -0.5 *
      (block_size * log(2 * pi) + log_det + sum(whitened^2))
  }

  return(out)
}


.reference_factor_covariance_lower <- function(residual_sd, loading, block_size) {

  S     <- nrow(residual_sd)
  rank  <- if (block_size == 0L) 0L else ncol(loading) %/% block_size
  index <- which(lower.tri(matrix(0, block_size, block_size), diag = TRUE),
                 arr.ind = TRUE)
  out   <- matrix(0, nrow = S, ncol = nrow(index))
  for (column in seq_len(nrow(index))) {
    i     <- index[[column, 1L]]
    j     <- index[[column, 2L]]
    value <- if (i == j) residual_sd[, i]^2 else numeric(S)
    for (factor in seq_len(rank)) {
      value <- value + loading[, (factor - 1L) * block_size + i] *
        loading[, (factor - 1L) * block_size + j]
    }
    out[, column] <- value
  }

  return(out)
}


# Packed lower triangles of positive-definite blocks, laid out so that the run
# structure the vectorized whitening exploits is exercised deliberately:
# a run of several equal draws, single-draw runs, and a returning covariance.
.grid_packed_covariance <- function(block_size, seed) {

  withr::local_seed(seed)
  root  <- matrix(stats::rnorm(block_size * block_size), block_size)
  value <- crossprod(root) + diag(block_size) * (block_size + 1)

  return(as.numeric(value[lower.tri(value, diag = TRUE)]))
}


test_that("the vectorized normalizer-grid Gaussian density is unchanged", {

  for (block_size in c(1L, 2L, 3L, 5L, 8L)) {
    packed <- lapply(seq_len(4L), function(index) {
      .grid_packed_covariance(block_size, 20260918L + index)
    })
    # Run lengths 3, 1, 1, 2, 4 with the first covariance returning at the end.
    pattern <- c(1L, 1L, 1L, 2L, 3L, 4L, 4L, 1L, 1L, 1L, 1L)
    covariance_lower <- do.call(rbind, packed[pattern])
    S <- nrow(covariance_lower)

    withr::local_seed(900L + block_size)
    yi    <- stats::rnorm(block_size)
    means <- matrix(stats::rnorm(S * block_size, sd = 0.7), nrow = S)

    expect_identical(
      .selection_grid_gaussian_lpdf(yi, means, covariance_lower, block_size),
      .reference_grid_gaussian_lpdf(yi, means, covariance_lower, block_size)
    )

    # A single draw is its own run of length one.
    expect_identical(
      .selection_grid_gaussian_lpdf(yi, means[1L, , drop = FALSE],
        covariance_lower[1L, , drop = FALSE], block_size),
      .reference_grid_gaussian_lpdf(yi, means[1L, , drop = FALSE],
        covariance_lower[1L, , drop = FALSE], block_size)
    )

    # Every draw its own covariance: only runs of length one.
    distinct <- do.call(rbind, lapply(seq_len(S), function(draw) {
      .grid_packed_covariance(block_size, 5000L + draw)
    }))
    expect_identical(
      .selection_grid_gaussian_lpdf(yi, means, distinct, block_size),
      .reference_grid_gaussian_lpdf(yi, means, distinct, block_size)
    )
  }
})


test_that("a row-compacted covariance gives the same Gaussian density", {

  for (block_size in c(1L, 3L, 6L)) {
    withr::local_seed(700L + block_size)
    states      <- rep(seq_len(7L), each = 5L)
    S           <- length(states)
    residual_sd <- matrix(stats::runif(7L * block_size, 0.2, 1.1),
                          nrow = 7L)[states, , drop = FALSE]
    loading     <- matrix(stats::rnorm(7L * block_size * 2L),
                          nrow = 7L)[states, , drop = FALSE]
    yi          <- stats::rnorm(block_size)
    means       <- matrix(stats::rnorm(S * block_size), nrow = S)

    full    <- .selection_factor_covariance_lower(residual_sd, loading, block_size)
    compact <- .selection_factor_covariance_rows(residual_sd, loading, block_size, states)
    expect_false(is.null(compact[["rows"]]))
    expect_identical(nrow(compact[["values"]]), 7L)
    expect_identical(compact[["values"]][compact[["rows"]], , drop = FALSE], full)

    expect_identical(
      .selection_grid_gaussian_lpdf(yi, means, compact[["values"]], block_size,
        covariance_rows = compact[["rows"]]),
      .reference_grid_gaussian_lpdf(yi, means, full, block_size)
    )

    # A subset of rows keeps the mapping.
    rows <- c(2L, 3L, 9L, 10L, 11L, 30L)
    expect_identical(
      .selection_grid_gaussian_lpdf(yi, means[rows, , drop = FALSE],
        compact[["values"]], block_size,
        covariance_rows = compact[["rows"]][rows]),
      .reference_grid_gaussian_lpdf(yi, means[rows, , drop = FALSE],
        full[rows, , drop = FALSE], block_size)
    )

    # Factor inputs that differ inside a state keep the full matrix.
    varied <- residual_sd
    varied[2L, 1L] <- varied[2L, 1L] + 0.25
    plain <- .selection_factor_covariance_rows(varied, loading, block_size, states)
    expect_null(plain[["rows"]])
    expect_identical(plain[["values"]],
      .selection_factor_covariance_lower(varied, loading, block_size))

    # No state information at all also keeps the full matrix.
    expect_null(.selection_factor_covariance_rows(residual_sd, loading,
      block_size, NULL)[["rows"]])
  }
})


test_that("a non-positive-definite normalizer-grid block still declines", {

  block_size <- 3L
  good       <- .grid_packed_covariance(block_size, 77L)
  bad        <- good
  bad[[1L]]  <- -1

  withr::local_seed(31L)
  yi    <- stats::rnorm(block_size)
  means <- matrix(stats::rnorm(4L * block_size), nrow = 4L)

  # The failing block is reached on the first draw and on a later draw.
  for (pattern in list(c(1L, 1L, 2L, 2L), c(2L, 2L, 1L, 1L), c(2L, 1L, 2L, 1L))) {
    covariance_lower <- do.call(rbind, list(bad, good)[pattern])
    expect_null(.selection_grid_gaussian_lpdf(yi, means, covariance_lower, block_size))
    expect_null(.reference_grid_gaussian_lpdf(yi, means, covariance_lower, block_size))
  }

  # A zero-row batch keeps its empty result.
  expect_identical(
    .selection_grid_gaussian_lpdf(yi, means[0L, , drop = FALSE],
      matrix(good, nrow = 0L, ncol = length(good)), block_size),
    numeric(0)
  )
})


test_that("the vectorized packed factor covariance is unchanged", {

  for (block_size in c(1L, 2L, 4L, 7L)) {
    for (rank in c(0L, 1L, 2L, 3L)) {
      withr::local_seed(4000L + 100L * block_size + rank)
      S           <- 11L
      residual_sd <- matrix(stats::runif(S * block_size, 0.05, 1.5), nrow = S)
      loading     <- matrix(stats::rnorm(S * block_size * rank),
                            nrow = S, ncol = block_size * rank)

      expect_identical(
        .selection_factor_covariance_lower(residual_sd, loading, block_size),
        .reference_factor_covariance_lower(residual_sd, loading, block_size)
      )
    }
  }

  # One draw, and a zero-size block, keep their shapes.
  residual_sd <- matrix(c(0.4, 1.1, 0.7), nrow = 1L)
  loading     <- matrix(c(0.2, -0.5, 0.9, 1.3, -0.1, 0.3), nrow = 1L)
  expect_identical(
    .selection_factor_covariance_lower(residual_sd, loading, 3L),
    .reference_factor_covariance_lower(residual_sd, loading, 3L)
  )
  expect_identical(
    .selection_factor_covariance_lower(matrix(numeric(0), nrow = 2L, ncol = 0L),
      matrix(numeric(0), nrow = 2L, ncol = 0L), 0L),
    .reference_factor_covariance_lower(matrix(numeric(0), nrow = 2L, ncol = 0L),
      matrix(numeric(0), nrow = 2L, ncol = 0L), 0L)
  )
})
