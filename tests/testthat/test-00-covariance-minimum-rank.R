# Exact minimum-rank recovery. Covariance from overlapping samples or shared
# arms is `diag(d) + U U'` with few columns without being block constant, and a
# representation is accepted only when it reproduces every supplied entry.

.minimum_rank_reconstruction <- function(factor) {

  out <- tcrossprod(factor[["loading"]])
  diag(out) <- diag(out) + factor[["diagonal"]]
  out
}


.minimum_rank_residual <- function(factor, covariance) {

  max(abs(.minimum_rank_reconstruction(factor) - covariance)) /
    max(abs(covariance))
}


.minimum_rank_context <- function(sei, draws, weights = c(1, .7, .35)) {

  prior <- BayesTools::prior_weightfunction("one-sided", steps = c(.025, .05),
    weights = BayesTools::wf_fixed(weights))
  spec <- .selection_spec(priors = list(outcome = list(bias = prior)),
    yi = rep(0, length(sei)), sei = sei, effect_direction = "positive",
    signed_data = FALSE)
  spec[["omega"]]        <- matrix(weights, draws, length(weights), byrow = TRUE)
  spec[["kernel_mode"]]  <- rep(SELKERNEL_STEP, draws)
  spec[["vector_rule"]]  <- rep(0L, draws)
  spec[["native_cache"]] <- new.env(parent = emptyenv())
  spec
}


.minimum_rank_plan <- function(k, method, rank) {

  .selection_joint_execution_plan(
    row_blocks = list(seq_len(k)), block_methods = method,
    factor_ranks = rank, selection_control = set_selection_likelihood_control(),
    sampling = NULL, sampling_factor_blocks = NULL, random_covariance = NULL)
}


.minimum_rank_fixture <- function(size, rank, seed) {

  set.seed(seed)
  loading  <- matrix(stats::rnorm(size * rank, sd = .3), size, rank)
  diagonal <- stats::runif(size, .15, .3)
  covariance <- tcrossprod(loading)
  diag(covariance) <- diag(covariance) + diagonal
  .known_v_exact_symmetrize(covariance)
}


test_that("exactly representable blocks recover exactly, at or below their rank", {

  # The minimum rank can be below the rank a block was built from: a small
  # block has more free loadings than off-diagonal entries, so a lower rank can
  # reproduce it exactly. Only the certificate is a contract.
  for (rank in 1:6) {
    for (size in c(max(2L * rank + 2L, 5L), 14L)) {
      covariance <- .minimum_rank_fixture(size, rank, 100L * rank + size)
      factor <- .covariance_minimum_rank_factor(covariance)
      label <- sprintf("size %d rank %d", size, rank)
      expect_false(is.null(factor), info = label)
      expect_lte(factor[["rank"]], rank, label = label)
      expect_identical(factor[["method"]], "minimum_rank", info = label)
      expect_lte(.minimum_rank_residual(factor, covariance),
                 8 * size * .Machine$double.eps)
      expect_true(all(factor[["diagonal"]] > 0), info = label)
      # The reported residual is the certificate the recovery accepted on.
      expect_equal(factor[["residual"]],
                   .minimum_rank_residual(factor, covariance),
                   tolerance = 0)
    }
  }

  # Exactly representable blocks whose rank the kernels cannot accept keep
  # their supplied entries and their general route.
  expect_null(.covariance_minimum_rank_factor(
    .minimum_rank_fixture(24L, 10L, 991L)
  ))
})


test_that("a signed two-row block recovers where block-constant declines", {

  covariance <- matrix(c(.04, -.012, -.012, .05), 2L)
  expect_null(.covariance_block_constant_factor(covariance))

  factor <- .covariance_minimum_rank_factor(covariance)
  expect_identical(factor[["rank"]], 1L)
  expect_lt(prod(factor[["loading"]]), 0)
  expect_equal(.minimum_rank_reconstruction(factor), covariance,
               tolerance = 1e-15)
  expect_true(all(factor[["diagonal"]] > 0))
})


test_that("an overlapping shared-arm block recovers at rank two", {

  # Two overlapping groups of rows: comparisons one to three share one arm and
  # comparisons three to five share another. The groups are not nested, so no
  # correlation level partitions the block and block-constant recovery stops.
  loading <- matrix(0, 5L, 2L)
  loading[1:3, 1L] <- .18
  loading[3:5, 2L] <- .15
  covariance <- tcrossprod(loading)
  diag(covariance) <- diag(covariance) + seq(.02, .06, length.out = 5L)
  covariance <- .known_v_exact_symmetrize(covariance)

  expect_null(.covariance_block_constant_factor(covariance))
  factor <- .covariance_minimum_rank_factor(covariance)
  expect_identical(factor[["rank"]], 2L)
  expect_lte(.minimum_rank_residual(factor, covariance),
             8 * nrow(covariance) * .Machine$double.eps)
})


test_that("recovery declines what it cannot represent exactly", {

  # The cap is a choice about what the kernels can integrate, not about what
  # exists: the same block is representable once the search may reach its rank.
  expect_null(.covariance_minimum_rank_factor(
    .minimum_rank_fixture(24L, 10L, 4242L)
  ))
  deep <- .covariance_minimum_rank_factor(
    .minimum_rank_fixture(24L, 10L, 4242L), max_rank = 10L
  )
  expect_lte(deep[["rank"]], 10L)

  # A singular block has an exact representation only with a vanishing
  # residual variance, which is not a certifiably positive one.
  set.seed(515L)
  loading <- matrix(stats::rnorm(8L, sd = .3), 4L, 2L)
  expect_null(.covariance_minimum_rank_factor(
    .known_v_exact_symmetrize(tcrossprod(loading))
  ))

  # Invalid or trivial inputs are declined rather than repaired.
  expect_null(.covariance_minimum_rank_factor(matrix(.04, 1L, 1L)))
  expect_null(.covariance_minimum_rank_factor(diag(c(.04, .05))))
  expect_null(.covariance_minimum_rank_factor(
    matrix(c(.04, .01, .02, .05), 2L)
  ))
})


test_that("recovery does not depend on the row order", {

  covariance <- .minimum_rank_fixture(9L, 3L, 7171L)
  order      <- c(5L, 1L, 9L, 3L, 7L, 2L, 8L, 4L, 6L)
  permuted   <- covariance[order, order, drop = FALSE]

  reference <- .covariance_minimum_rank_factor(covariance)
  observed  <- .covariance_minimum_rank_factor(permuted)
  expect_identical(observed[["rank"]], reference[["rank"]])
  expect_equal(observed[["diagonal"]], reference[["diagonal"]][order],
               tolerance = 1e-12)
  expect_lte(.minimum_rank_residual(observed, permuted),
             8 * nrow(permuted) * .Machine$double.eps)
})


test_that("a minimum-rank block reaches the certified factor metadata", {

  loading <- matrix(0, 5L, 2L)
  loading[1:3, 1L] <- .18
  loading[3:5, 2L] <- .15
  covariance <- tcrossprod(loading)
  diag(covariance) <- diag(covariance) + seq(.02, .06, length.out = 5L)
  covariance <- .known_v_exact_symmetrize(covariance)

  known_V <- .known_v_canonicalize(covariance)
  expect_identical(.known_v_certified_factor_status(known_V), "recovered")
  blocks <- .known_v_certified_factor_blocks(known_V)
  expect_length(blocks, 1L)
  expect_identical(ncol(blocks[[1L]][["loading"]]), 2L)
  expect_identical(known_V[["certified_factor"]][["blocks"]][[1L]][["method"]],
                   "minimum_rank")
  # The validator re-checks the identity, so a drifted representation fails.
  expect_silent(.known_v_validate_certified_factor(known_V))
  damaged <- known_V
  damaged[["certified_factor"]][["blocks"]][[1L]][["loading"]][1L, 1L] <- .3
  expect_error(.known_v_validate_certified_factor(damaged),
               "no longer reproduces", fixed = TRUE)
})


test_that("a minimum-rank block evaluates the same selected law as the dense route", {

  loading <- matrix(0, 5L, 2L)
  loading[1:3, 1L] <- .18
  loading[3:5, 2L] <- .15
  covariance <- tcrossprod(loading)
  diag(covariance) <- diag(covariance) + seq(.02, .06, length.out = 5L)
  covariance <- .known_v_exact_symmetrize(covariance)

  known_V <- .known_v_canonicalize(covariance)
  block   <- .known_v_certified_factor_blocks(known_V)[[1L]]
  index   <- block[["index"]]
  k       <- length(index)
  rank    <- ncol(block[["loading"]])
  yi      <- c(-.30, .12, .41, -.08, .22)
  sei     <- sqrt(diag(covariance))
  tolerance <- set_selection_likelihood_control()[["relative_tolerance"]]

  set.seed(3350)
  draws   <- 60L
  context <- .minimum_rank_context(sei, draws)
  means   <- matrix(stats::rnorm(draws * k, .15, .10), draws, k)
  extra   <- stats::runif(draws, .01, .06)

  residual_sd <- sqrt(matrix(block[["diagonal"]], draws, k, byrow = TRUE) +
                        matrix(extra, draws, k))
  factor_loading <- matrix(as.numeric(block[["loading"]]), draws, k * rank,
                           byrow = TRUE)
  observed <- .selection_joint_factor_loglik_block(
    yi = yi, means = means, residual_sd = residual_sd, loading = factor_loading,
    sei = sei, selection_context = context,
    execution_plan = .minimum_rank_plan(k, "factor", rank),
    block_index = 1L, return_normalizer = TRUE)

  pairs <- .selection_joint_lower_pairs(
    .minimum_rank_plan(k, "dense", NA_integer_), seq_len(k))
  packed <- matrix(covariance[cbind(pairs[["row_1"]], pairs[["row_2"]])],
                   draws, length(pairs[["row_1"]]), byrow = TRUE)
  diagonal <- pairs[["row_1"]] == pairs[["row_2"]]
  packed[, diagonal] <- packed[, diagonal] + extra
  reference <- .selection_joint_dense_loglik_block(
    yi = yi, means = means, covariance_lower = packed, sei = sei,
    selection_context = context,
    execution_plan = .minimum_rank_plan(k, "dense", NA_integer_),
    block_size = k, return_normalizer = TRUE)

  expect_true(all(is.finite(observed[["log_density"]])))
  expect_identical(observed[["relative_mcse"]], rep(0, draws))
  budget <- tolerance + max(observed[["relative_change"]],
                            reference[["relative_mcse"]])
  expect_lt(max(abs(observed[["log_density"]] - reference[["log_density"]])),
            budget)
})


test_that("a block with no low-rank structure declines without refinement", {

  # An exact `D + U U'` makes every off-diagonal submatrix `V[I, J]` equal
  # `U[I, ] U[J, ]'`, so its rank cannot exceed the factor rank. That bound
  # settles an unstructured or slowly decaying block before any Gauss-Newton
  # runs, which is what keeps model construction with an arbitrary dense 'V'
  # from paying a minute at input.
  set.seed(4242L)
  unstructured <- function(size) {
    entries <- matrix(stats::rnorm(size * size), size, size)
    .known_v_exact_symmetrize(crossprod(entries) / size + diag(size))
  }
  autoregressive <- function(size, phi = .8) {
    correlation <- outer(seq_len(size), seq_len(size),
                         function(i, j) phi^abs(i - j))
    .known_v_exact_symmetrize(
      correlation * tcrossprod(sqrt(seq(.02, .06, length.out = size)))
    )
  }

  for (block in list(unstructured(80L), autoregressive(80L))) {
    elapsed <- system.time(
      factor <- .covariance_minimum_rank_factor(block)
    )[["elapsed"]]
    expect_null(factor)
    expect_lt(elapsed, 1)
  }

  # The same for the whole input path a user reaches with a dense 'V'.
  dense   <- unstructured(150L)
  elapsed <- system.time(known_V <- .known_v_canonicalize(dense))[["elapsed"]]
  expect_false(.known_v_has_certified_factor(known_V))
  expect_lt(elapsed, 2)

  # The bound is necessary, never sufficient: exact structure at the same sizes
  # is still recovered.
  for (rank in 1:3) {
    exact     <- .minimum_rank_fixture(80L, rank, 909L + rank)
    recovered <- .covariance_minimum_rank_factor(exact)
    expect_false(is.null(recovered))
    expect_lte(recovered[["rank"]], rank)
    expect_lte(.minimum_rank_residual(recovered, exact),
               8 * 80 * .Machine$double.eps)
  }
})
