# ============================================================================ #
# Block-diagonal posterior covariances
# ============================================================================ #

# The joint selection prediction parts carry their covariances one block at a
# time instead of as dense S x K x K cubes. Every accessor has to answer what
# the dense cube answered, element for element, on a partition with several
# block sizes at once.

.block_covariance_toy <- function(S = 4L, seed = 4181L) {

  set.seed(seed)
  blocks <- list(1:3, 4:5, 6L)
  K      <- 6L
  values <- lapply(blocks, function(rows) {
    n     <- length(rows)
    value <- array(0, dim = c(S, n, n))
    for (draw in seq_len(S)) {
      basis <- matrix(stats::rnorm(n * n), n, n)
      value[draw, , ] <- crossprod(basis) + diag(n) * n
    }
    value
  })
  dense <- lapply(seq_len(S), function(draw) {
    out <- matrix(0, K, K)
    for (index in seq_along(blocks)) {
      rows <- blocks[[index]]
      out[rows, rows] <- matrix(values[[index]][draw, , ],
                                length(rows), length(rows))
    }
    out
  })
  list(object = .block_covariance(blocks, values, S, K), blocks = blocks,
       values = values, dense = dense, S = S, K = K)
}


test_that("a block covariance answers what its dense cube answered", {

  toy <- .block_covariance_toy()
  x   <- toy[["object"]]

  expect_true(.is_block_covariance(x))
  expect_identical(.block_covariance_dim(x), c(toy[["S"]], toy[["K"]], toy[["K"]]))

  # Dense assembly, element for element.
  for (draw in seq_len(toy[["S"]])) {
    expect_identical(.block_covariance_dense(x, draw), toy[["dense"]][[draw]])
  }

  # The per-row variances `cube[, k, k]` gave.
  expect_identical(
    .block_covariance_diag_matrix(x),
    matrix(vapply(seq_len(toy[["K"]]), function(k) {
      vapply(toy[["dense"]], function(matrix) matrix[k, k], numeric(1L))
    }, numeric(toy[["S"]])), toy[["S"]], toy[["K"]])
  )

  # `cube[draw, rows, rows]` for one block, a union of blocks, a subset of a
  # block, and an unsorted selection that straddles two blocks.
  for (rows in list(1:3, 4:5, 6L, c(1:3, 6L), c(2L, 3L, 5L), c(5L, 2L, 6L, 1L))) {
    for (draw in seq_len(toy[["S"]])) {
      expect_identical(.block_covariance_sub(x, draw, rows),
                       toy[["dense"]][[draw]][rows, rows, drop = FALSE])
    }
  }

  # `cube[draw, , ] %*% v` is bit-identical block-wise: what the block form
  # leaves out is exactly zero.
  vector <- stats::rnorm(toy[["K"]])
  for (draw in seq_len(toy[["S"]])) {
    expect_identical(.block_covariance_matvec(x, draw, vector),
                     as.vector(toy[["dense"]][[draw]] %*% vector))
  }

  # The block-wise solve against the dense factorization. LAPACK factorizes an
  # assembled matrix and its blocks with different blockings, so these agree to
  # rounding, not bit for bit; callers that need bit-identical results take the
  # assembled matrix from `.block_covariance_dense()` instead.
  rhs <- stats::rnorm(toy[["K"]])
  for (draw in seq_len(toy[["S"]])) {
    factor <- chol(toy[["dense"]][[draw]])
    expect_equal(
      .block_covariance_chol_solve(x, draw, rhs),
      backsolve(factor, forwardsolve(t(factor), rhs)),
      tolerance = 1e-12
    )
  }

  # A subset of the draws.
  subset <- .block_covariance_draws(x, c(3L, 1L))
  expect_identical(.block_covariance_dim(subset), c(2L, toy[["K"]], toy[["K"]]))
  expect_identical(.block_covariance_dense(subset, 1L), toy[["dense"]][[3L]])
  expect_identical(.block_covariance_dense(subset, 2L), toy[["dense"]][[1L]])

  expect_true(.block_covariance_all_finite(x))
  expect_false(.block_covariance_is_zero(x))
  expect_identical(.block_covariance_zero_draws(x), rep(FALSE, toy[["S"]]))
  expect_identical(.block_covariance_diagonal_draws(x), rep(FALSE, toy[["S"]]))
})


test_that("zero and diagonal block covariances answer structurally", {

  S <- 3L
  K <- 4L
  empty <- .block_covariance_zero(S, K)
  expect_true(.block_covariance_is_zero(empty))
  expect_identical(.block_covariance_zero_draws(empty), rep(TRUE, S))
  expect_identical(.block_covariance_diagonal_draws(empty), rep(TRUE, S))
  expect_identical(.block_covariance_dense(empty, 2L), matrix(0, K, K))
  expect_identical(.block_covariance_diag_matrix(empty), matrix(0, S, K))

  diagonal <- rbind(c(.4, .7, 0, .2), c(0, 0, 0, 0), c(.3, .8, .1, .5))
  x <- .block_covariance_from_diagonal(diagonal)
  expect_identical(.block_covariance_diag_matrix(x), diagonal)
  expect_identical(.block_covariance_zero_draws(x), c(FALSE, TRUE, FALSE))
  # Singleton blocks have no off-diagonal element to scan.
  expect_identical(.block_covariance_diagonal_draws(x), rep(TRUE, S))
  expect_identical(.block_covariance_dense(x, 1L), diag(diagonal[1L, ]))
  expect_false(.block_covariance_is_zero(x))

  # A block that is diagonal in one draw and not in another.
  varying <- array(0, c(3L, 2L, 2L))
  varying[, 1L, 1L] <- 1
  varying[, 2L, 2L] <- 1
  varying[2L, 1L, 2L] <- varying[2L, 2L, 1L] <- 2
  mixed <- .block_covariance(
    list(1:2, 3:4), list(varying, array(0, c(1L, 2L, 2L))), S, K
  )
  expect_identical(.block_covariance_diagonal_draws(mixed), c(TRUE, FALSE, TRUE))
  expect_identical(.block_covariance_zero_draws(mixed), rep(FALSE, S))
})


test_that("a constant block covariance is stored once and adds over partitions", {

  S <- 5L
  K <- 4L
  sampling <- matrix(0, K, K)
  sampling[1:2, 1:2] <- matrix(c(.09, .02, .02, .12), 2L)
  sampling[3L, 3L]   <- .1
  sampling[4L, 4L]   <- .2
  constant <- .block_covariance_from_matrix(sampling, S)
  expect_identical(lapply(constant[["blocks"]], as.integer), list(1:2, 3L, 4L))
  # One value per block, not one per draw.
  expect_identical(
    vapply(constant[["values"]], function(value) dim(value)[[1L]], integer(1L)),
    rep(1L, 3L)
  )
  for (draw in seq_len(S)) {
    expect_identical(.block_covariance_dense(constant, draw), sampling)
  }

  # Adding a finer partition assembles on the coarser one and keeps the values.
  diagonal <- matrix(seq_len(S * K) / 10, S, K)
  total <- .block_covariance_add(
    .block_covariance_from_diagonal(diagonal), constant
  )
  expect_identical(lapply(total[["blocks"]], as.integer), list(1:2, 3L, 4L))
  for (draw in seq_len(S)) {
    expect_identical(.block_covariance_dense(total, draw),
                     diag(diagonal[draw, ]) + sampling)
  }

  # Two partitions where neither refines the other assemble on their join.
  left  <- .block_covariance(list(1:2, 3L, 4L),
    list(array(1, c(1L, 2L, 2L)), array(1, c(1L, 1L, 1L)), array(1, c(1L, 1L, 1L))),
    S, K)
  right <- .block_covariance(list(1L, 2:3, 4L),
    list(array(1, c(1L, 1L, 1L)), array(1, c(1L, 2L, 2L)), array(1, c(1L, 1L, 1L))),
    S, K)
  joined <- .block_covariance_add(left, right)
  expect_identical(lapply(joined[["blocks"]], as.integer), list(1:3, 4L))

  # The diagonal update promotes a constant block to a per-draw one.
  updated <- .block_covariance_add_diagonal(constant, diagonal)
  for (draw in seq_len(S)) {
    expect_identical(.block_covariance_dense(updated, draw),
                     sampling + diag(diagonal[draw, ]))
  }
})


test_that("a block covariance refuses values outside its declared blocks", {

  S <- 2L
  K <- 3L
  samples <- array(0, c(S, K, K))
  samples[, 1L, 1L] <- 1
  samples[, 2L, 2L] <- 1
  samples[, 3L, 3L] <- 1
  blocks <- list(1:2, 3L)
  expect_identical(
    .block_covariance_dense(.block_covariance_from_array(samples, blocks), 1L),
    diag(3L)
  )
  samples[1L, 1L, 3L] <- .5
  error <- tryCatch(.block_covariance_from_array(samples, blocks),
                    error = identity)
  expect_identical(
    conditionMessage(error),
    paste0("Covariance samples are nonzero outside their declared dependency ",
           "blocks.")
  )
  expect_error(.block_covariance(list(1:2), list(array(0, c(S, 2L, 2L))), S, K),
               "Known-V block metadata must partition the fitted rows.")
  expect_error(.block_covariance(list(1:2, 3L),
                                 list(array(0, c(S, 3L, 3L)), array(0, c(S, 1L, 1L))),
                                 S, K),
               "Block covariance values must be draw x row x row arrays.")
})


test_that("the selected response kernel draws the same from blocks as from a cube", {

  S <- 400L
  K <- 4L
  block_matrix <- matrix(c(1, .3, .3, 1), 2L)
  cube <- array(0, c(S, K, K))
  cube[, 1:2, 1:2] <- array(rep(block_matrix, each = S), c(S, 2L, 2L))
  cube[, 3L, 3L]   <- 1
  cube[, 4L, 4L]   <- 1
  blocked <- .block_covariance_from_array(cube, list(1:2, 3L, 4L))
  arguments <- list(
    matrix(0, S, K), NULL, rep(1, K),
    matrix(rep(c(1, .3), each = S), S), c(0, -Inf), c(Inf, 0),
    1L, 1L, list(1:2, 3:4), 1000L, 0L
  )
  dense_call <- function(covariance) {
    arguments[[2L]] <- covariance
    do.call(.Call, c(list("RoBMA_selnorm_mnorm_step_rng_batch"), arguments,
                     list(PACKAGE = "RoBMA")))
  }
  set.seed(4182)
  dense <- dense_call(cube)
  set.seed(4182)
  blocks <- dense_call(.block_covariance_native_parts(blocked))
  expect_identical(dense[["failure_code"]], 0L)
  expect_identical(blocks, dense)

  # A storage partition that straddles two dependency blocks is refused, and a
  # cube that is nonzero across them still is.
  straddling <- .block_covariance_from_array(cube, list(1:2, 3L, 4L))
  arguments[[9L]] <- list(1L, 2L, 3L, 4L)
  expect_error(
    dense_call(.block_covariance_native_parts(straddling)),
    "'dependency_blocks' must include all observations with nonzero covariance."
  )
  expect_error(
    dense_call(cube),
    "'dependency_blocks' must include all observations with nonzero covariance."
  )
  arguments[[9L]] <- list(1:2, 3:4)
  expect_error(
    dense_call(list(list(array(0, c(S, 2L, 2L))), list(1:2))),
    "'covariance' blocks must partition the observations."
  )
  expect_error(
    dense_call(list(list(array(0, c(S, 3L, 3L))), list(1:2))),
    "'covariance' block values must be draw x row x row arrays."
  )
})
