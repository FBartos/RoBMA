.guarded_prediction_fixture <- function(variable = FALSE, single = FALSE) {

  blocks <- list(c(1L, 4L), c(2L, 3L, 5L))
  G <- matrix(0, 5L, 5L)
  G[blocks[[1L]], blocks[[1L]]] <- tcrossprod(c(.5, 1))
  G[blocks[[2L]], blocks[[2L]]] <- tcrossprod(c(.25, .5, 1))
  N <- matrix(0, 5L, 5L)
  N[blocks[[1L]], blocks[[1L]]] <- matrix(c(.25, .0625, .0625, .5), 2L)
  N[blocks[[2L]], blocks[[2L]]] <- matrix(c(.5, .125, 0, .125, .75, .0625, 0, .0625, .25), 3L)
  if (single) blocks <- list(1:5)
  covariance <- function(matrix) {
    if (!variable) return(.block_covariance_from_matrix(matrix, 3L, blocks))
    values <- lapply(blocks, function(rows) {
      value <- array(0, c(3L, length(rows), length(rows)))
      for (s in 1:3) value[s, , ] <- matrix[rows, rows, drop = FALSE] * c(1, 4, .25)[s]
      value
    })
    .block_covariance(blocks, values, 3L, 5L)
  }
  list(G = covariance(G), N = covariance(N), blocks = blocks)
}


test_that("guarded block sampling preserves odd draw and row normal assignments", {

  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  for (kind in c("Inversion", "Box-Muller")) for (variable in c(FALSE, TRUE)) {
    for (single in c(FALSE, TRUE)) for (source in c("G", "N")) {
      fixture <- .guarded_prediction_fixture(variable, single)
      x <- fixture[[source]]
      do.call(RNGkind, list("Mersenne-Twister", kind))
      set.seed(11621)
      expected <- matrix(0, 3L, 5L)
      for (s in 1:3) {
        if (single) {
          factor <- .covariance_sampling_factor(.covariance_factorization(.block_covariance_dense(x, s)))
          expected[s, ] <- as.vector(stats::rnorm(5L) %*% factor)
        } else {
          normals <- stats::rnorm(5L)
          for (rows in fixture$blocks) {
            C <- .block_covariance_sub(x, s, rows)
            factor <- if (source == "N") chol(C) else {
              result <- matrix(0, length(rows), length(rows))
              result[1L, ] <- sqrt(diag(C))
              result
            }
            expected[s, rows] <- as.vector(normals[rows] %*% factor)
          }
        }
      }
      expected_seed <- .Random.seed
      expected_next <- stats::rnorm(5L)
      do.call(RNGkind, list("Mersenne-Twister", kind))
      set.seed(11621)
      actual <- .outcome_rng.norm_known_v_covariance(matrix(0, 3L, 5L), x)
      expect_equal(actual, expected, tolerance = 0)
      expect_identical(.Random.seed, expected_seed)
      expect_identical(stats::rnorm(5L), expected_next)
    }
  }
})


test_that("guard eligibility preserves the whole-matrix boundary classifier", {

  K <- 20L
  blocks <- split(seq_len(K), rep(1:10, each = 2L))
  cases <- list(rho_half = .5, gap_1e12 = 1 - 1e-12, gap_1e14 = 1 - 1e-14,
    rho_one = 1, gap_negative = 1 + 1e-14, zero = NA_real_)
  for (name in names(cases)) {
    rho <- cases[[name]]
    C <- matrix(0, K, K)
    if (!is.na(rho)) for (rows in blocks) C[rows, rows] <- matrix(c(1, rho, rho, 1), 2L)
    x <- .block_covariance_from_matrix(C, 3L, blocks)
    classification <- .covariance_factorization(C)
    expected_rank <- if (name == "zero") 0L else
      if (name %in% c("rho_half", "gap_1e12")) K else 10L
    expect_identical(nrow(classification$sampling_factor), expected_rank)
    expect_identical(classification$status, if (expected_rank == K)
      "positive_definite" else "positive_semidefinite")
    plan <- .block_covariance_guard_plan(x, 1L, .block_covariance_guard_cache(x))
    if (name %in% c("gap_1e12", "gap_1e14", "gap_negative")) {
      expect_false(plan$fast)
      set.seed(13051)
      expected <- matrix(stats::rnorm(3L * K), 3L, K, byrow = TRUE) %*%
        .covariance_sampling_factor(classification)
      set.seed(13051)
      expect_equal(.outcome_rng.norm_known_v_covariance(matrix(0, 3L, K), x), expected, tolerance = 0)
    } else expect_true(plan$fast)
  }
  ordinary <- matrix(c(1, .25, .25, 1), 2L)
  expect_identical(.block_covariance_guard_block(ordinary, 1e308)$reason, "invalid_gamma")
  expect_true(.block_covariance_guard_block(ordinary * 1e200, 20L)$eligible)
  expect_true(.block_covariance_guard_block(ordinary * 1e-200, 20L)$eligible)
  embedded <- matrix(0, 3L, 3L)
  embedded[c(1L, 3L), c(1L, 3L)] <- tcrossprod(c(.5, 1))
  plan <- .block_covariance_guard_block(embedded, 20L)
  expect_identical(plan$reason, "exact_rank_one")
  expect_equal(crossprod(plan$factor), embedded, tolerance = 0)
  expect_false(.block_covariance_guard_block(embedded, 20L, pd_only = TRUE)$eligible)
})


test_that("unexpected eligible block failures retain the original full fallback", {

  guard <- .block_covariance_guard_block
  fault <- new.env(parent = environment(guard))
  fault$.covariance_sampling_factor <- function(...) stop("injected factor failure")
  environment(guard) <- fault
  expect_identical(guard(matrix(c(1, .25, .25, 1), 2L), 5L)$reason,
    "eligible_factor_unavailable")
  fixture <- .guarded_prediction_fixture()
  x <- .block_covariance_add(fixture$G, fixture$N)
  solve <- .block_covariance_guarded_chol_solve
  fault <- new.env(parent = environment(solve))
  fault$.block_covariance_chol_solve <- function(...) stop("injected block solve failure")
  environment(solve) <- fault
  rhs <- c(.2, -.1, .4, .8, -.2)
  cache <- .block_covariance_guard_cache(x, factors = FALSE, pd_only = TRUE)
  C <- .block_covariance_dense(x, 1L)
  factor <- chol(C)
  expected <- backsolve(factor, forwardsolve(t(factor), rhs))
  expect_equal(solve(x, 1L, rhs, cache), expected, tolerance = 0)
})


test_that("invalid draw failure retains its original normal consumption point", {

  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  fixture <- .guarded_prediction_fixture(variable = TRUE)
  invalid <- fixture$N
  invalid$values[[1L]][2L, 1L, 2L] <- invalid$values[[1L]][2L, 1L, 2L] + .01
  for (kind in c("Inversion", "Box-Muller")) {
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(11621)
    stats::rnorm(5L)
    expected_seed <- .Random.seed
    expected_next <- stats::rnorm(5L)
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(11621)
    expect_error(.outcome_rng.norm_known_v_covariance(matrix(0, 3L, 5L), invalid),
      "Covariance must be symmetric.", fixed = TRUE)
    expect_identical(.Random.seed, expected_seed)
    expect_identical(stats::rnorm(5L), expected_next)
    parts <- list(means = matrix(0, 3L, 5L), latent_means = matrix(0, 3L, 5L),
      random_covariance = fixture$G, sampling_covariance = invalid,
      covariance = .block_covariance_add(fixture$G, invalid), dependency_blocks = fixture$blocks)
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(11621)
    stats::rnorm(20L)
    expected_seed <- .Random.seed
    expected_next <- stats::rnorm(5L)
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(11621)
    expect_error(.predict_joint_selection_source_posterior(parts, rep(1, 5L)),
      "Covariance must be symmetric.", fixed = TRUE)
    expect_identical(.Random.seed, expected_seed)
    expect_identical(stats::rnorm(5L), expected_next)
  }
})
