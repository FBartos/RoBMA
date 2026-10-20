.guarded_prediction_fixture <- function(variable = FALSE, single = FALSE, copies = 1L) {

  blocks <- list(c(1L, 4L), c(2L, 3L, 5L))
  G <- matrix(0, 5L, 5L)
  G[blocks[[1L]], blocks[[1L]]] <- tcrossprod(c(.5, 1))
  G[blocks[[2L]], blocks[[2L]]] <- tcrossprod(c(.25, .5, 1))
  N <- matrix(0, 5L, 5L)
  N[blocks[[1L]], blocks[[1L]]] <- matrix(c(.25, .0625, .0625, .5), 2L)
  N[blocks[[2L]], blocks[[2L]]] <- matrix(c(.5, .125, 0, .125, .75, .0625, 0, .0625, .25), 3L)
  if (copies > 1L) {
    G <- kronecker(diag(copies), G)
    N <- kronecker(diag(copies), N)
    blocks <- unlist(lapply(seq_len(copies), function(copy) {
      lapply(blocks, function(rows) rows + 5L * (copy - 1L))
    }), recursive = FALSE)
  }
  K <- nrow(G)
  if (single) blocks <- list(seq_len(K))
  covariance <- function(matrix) {
    if (!variable) return(.block_covariance_from_matrix(matrix, 3L, blocks))
    values <- lapply(blocks, function(rows) {
      value <- array(0, c(3L, length(rows), length(rows)))
      for (s in 1:3) value[s, , ] <- matrix[rows, rows, drop = FALSE] * c(1, 4, .25)[s]
      value
    })
    .block_covariance(blocks, values, 3L, K)
  }
  list(G = covariance(G), N = covariance(N), blocks = blocks)
}


test_that("large guarded blocks preserve odd draw and row normal assignments", {

  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  # A covariance that varies by draw is sampled blockwise only from the sampling
  # threshold on (30 copies give K = 150); below it the assembled route draws.
  for (kind in c("Inversion", "Box-Muller")) for (variable in c(FALSE, TRUE)) for (copies in c(7L, 30L)) {
    for (single in c(FALSE, TRUE)) for (source in c("G", "N")) {
      fixture <- .guarded_prediction_fixture(variable, single, copies = copies)
      K <- 5L * copies
      assembled <- single || (variable && K < BLOCK_COVARIANCE_GUARD_MIN_K_SAMPLING_VARYING)
      x <- fixture[[source]]
      do.call(RNGkind, list("Mersenne-Twister", kind))
      set.seed(11621)
      expected <- matrix(0, 3L, K)
      for (s in 1:3) {
        if (assembled) {
          factor <- .covariance_sampling_factor(.covariance_factorization(.block_covariance_dense(x, s)))
          expected[s, ] <- as.vector(stats::rnorm(K) %*% factor)
        } else {
          normals <- stats::rnorm(K)
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
      actual <- .outcome_rng.norm_known_v_covariance(matrix(0, 3L, K), x)
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
  # The block solve is reached from the solve threshold on (80 copies give K = 400).
  fixture <- .guarded_prediction_fixture(copies = 80L)
  x <- .block_covariance_add(fixture$G, fixture$N)
  solve <- .block_covariance_guarded_chol_solve
  fault <- new.env(parent = environment(solve))
  reached <- FALSE
  fault$.block_covariance_chol_solve <- function(...) {
    reached <<- TRUE
    stop("injected block solve failure")
  }
  environment(solve) <- fault
  rhs <- rep(c(.2, -.1, .4, .8, -.2), 80L)
  cache <- .block_covariance_guard_cache(x, factors = FALSE, pd_only = TRUE)
  C <- .block_covariance_dense(x, 1L)
  factor <- chol(C)
  expected <- backsolve(factor, forwardsolve(t(factor), rhs))
  expect_equal(solve(x, 1L, rhs, cache), expected, tolerance = 0)
  expect_true(reached)
})


test_that("declared singleton sampling batches row-ordered normals including zeros", {

  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  order <- c(5L, 2L, 4L, 1L, 3L)
  for (kind in c("Inversion", "Box-Muller")) for (variable in c(FALSE, TRUE)) {
    for (diagonal in list(c(.25, 0, 1, .5, 0), rep(0, 5L))) {
      variances <- matrix(diagonal, 3L, 5L, byrow = TRUE)
      if (variable) variances <- variances * c(1, 4, .25)
      x <- .block_covariance(lapply(order, as.integer), lapply(order, function(row) {
        array(if (variable) variances[, row] else variances[1L, row],
          c(if (variable) 3L else 1L, 1L, 1L))
      }), 3L, 5L)
      mu <- matrix(seq_len(15L) / 10, 3L, 5L)
      do.call(RNGkind, list("Mersenne-Twister", kind))
      set.seed(13061)
      expected <- mu + matrix(stats::rnorm(15L), 3L, 5L, byrow = TRUE) * sqrt(variances)
      expected_seed <- .Random.seed
      expected_next <- stats::rnorm(5L)
      do.call(RNGkind, list("Mersenne-Twister", kind))
      set.seed(13061)
      expect_equal(.outcome_rng.norm_known_v_covariance(mu, x), expected, tolerance = 0)
      expect_identical(.Random.seed, expected_seed)
      expect_identical(stats::rnorm(5L), expected_next)
    }
  }
})


test_that("bad later singleton draws retain the original classifier failure phase", {

  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  for (kind in c("Inversion", "Box-Muller")) for (bad in c(-1, Inf, NaN)) {
    x <- .block_covariance(lapply(1:5, as.integer),
      lapply(1:5, function(row) array(1, c(3L, 1L, 1L))), 3L, 5L)
    # Inject an invalid later draw after the constructor's finite-value check.
    x$values[[3L]][2L, 1L, 1L] <- bad
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(13062)
    stats::rnorm(5L)
    expected_seed <- .Random.seed
    expected_next <- stats::rnorm(5L)
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(13062)
    expect_error(.outcome_rng.norm_known_v_covariance(matrix(0, 3L, 5L), x))
    expect_identical(.Random.seed, expected_seed)
    expect_identical(stats::rnorm(5L), expected_next)
  }
})


test_that("small whole matrices bypass block guards for sampling and solving", {

  fixture <- .guarded_prediction_fixture(variable = TRUE)
  parts <- list(means = matrix(.1, 3L, 5L), latent_means = matrix(.2, 3L, 5L),
    random_covariance = fixture$G, sampling_covariance = fixture$N,
    covariance = .block_covariance_add(fixture$G, fixture$N), dependency_blocks = fixture$blocks)
  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  for (kind in c("Inversion", "Box-Muller")) {
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(13063)
    random <- sampling <- matrix(0, 3L, 5L)
    for (s in 1:3) random[s, ] <- as.vector(stats::rnorm(5L) %*%
      .covariance_sampling_factor(.covariance_factorization(.block_covariance_dense(fixture$G, s))))
    for (s in 1:3) sampling[s, ] <- as.vector(stats::rnorm(5L) %*%
      .covariance_sampling_factor(.covariance_factorization(.block_covariance_dense(fixture$N, s))))
    for (s in 1:3) {
      factor <- chol(.block_covariance_dense(parts$covariance, s))
      residual <- rep(.7, 5L) - .1 - random[s, ] - sampling[s, ]
      correction <- backsolve(factor, forwardsolve(t(factor), residual))
      random[s, ] <- random[s, ] + as.vector(.block_covariance_dense(fixture$G, s) %*% correction)
    }
    expected <- parts$latent_means + random
    expected_seed <- .Random.seed
    expected_next <- stats::rnorm(5L)
    posterior <- .predict_joint_selection_source_posterior
    fault <- new.env(parent = environment(posterior))
    fault$.block_covariance_guard_cache <- fault$.block_covariance_guard_plan <- function(...) {
      stop("small matrix reached block guard")
    }
    for (name in c(".outcome_rng.norm_known_v_covariance", ".block_covariance_guarded_chol_solve")) {
      fun <- get(name, envir = environment(posterior))
      environment(fun) <- fault
      assign(name, fun, envir = fault)
    }
    environment(posterior) <- fault
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(13063)
    expect_equal(posterior(parts, rep(.7, 5L)), expected, tolerance = 0)
    expect_identical(.Random.seed, expected_seed)
    expect_identical(stats::rnorm(5L), expected_next)
  }
  expect_false(.block_covariance_use_guarded_execution(.block_covariance_zero(3L, 32L)))
  expect_true(.block_covariance_use_guarded_execution(.block_covariance_zero(3L, 33L)))
})


test_that("guarded execution follows the draw dependence of the blocks and the consumer", {

  rule <- .block_covariance_use_guarded_execution
  covariance <- function(K, varying) {
    blocks <- lapply(seq(1L, K, by = 2L), function(first) as.integer(first:min(first + 1L, K)))
    values <- lapply(seq_along(blocks), function(index) {
      n <- length(blocks[[index]])
      draws <- if (varying[[index %% length(varying) + 1L]]) 3L else 1L
      array(rep(diag(n), each = draws), c(draws, n, n))
    })
    .block_covariance(blocks, values, 3L, K)
  }
  sampling_varying <- BLOCK_COVARIANCE_GUARD_MIN_K_SAMPLING_VARYING
  solve_invariant <- BLOCK_COVARIANCE_GUARD_MIN_K_SOLVE_INVARIANT
  solve_varying <- BLOCK_COVARIANCE_GUARD_MIN_K_SOLVE_VARYING
  # Sampling: draw-invariant blocks from the whole-matrix boundary, varying ones
  # from the sampling threshold (one varying block is enough to pay the guard).
  expect_false(rule(covariance(32L, FALSE), "sampling"))
  expect_true(rule(covariance(33L, FALSE), "sampling"))
  expect_false(rule(covariance(33L, TRUE), "sampling"))
  expect_false(rule(covariance(33L, c(FALSE, TRUE)), "sampling"))
  expect_false(rule(covariance(sampling_varying - 1L, TRUE), "sampling"))
  expect_true(rule(covariance(sampling_varying, TRUE), "sampling"))
  expect_true(rule(covariance(sampling_varying, c(FALSE, TRUE)), "sampling"))
  # Solving: invariant blocks from their threshold, varying ones from a larger one.
  expect_false(rule(covariance(33L, FALSE), "solve"))
  expect_false(rule(covariance(solve_invariant - 1L, FALSE), "solve"))
  expect_true(rule(covariance(solve_invariant, FALSE), "solve"))
  expect_false(rule(covariance(solve_varying - 1L, TRUE), "solve"))
  expect_true(rule(covariance(solve_varying, TRUE), "solve"))
  expect_false(rule(covariance(solve_varying - 1L, c(FALSE, TRUE)), "solve"))
  expect_true(rule(covariance(solve_varying, c(FALSE, TRUE)), "solve"))
  expect_identical(rule(covariance(sampling_varying, TRUE)),
    rule(covariance(sampling_varying, TRUE), "sampling"))
})


test_that("large eligible sampling factor faults actually reach the assembled fallback", {

  fixture <- .guarded_prediction_fixture(copies = 7L)
  sampling <- .outcome_rng.norm_known_v_covariance
  guards <- new.env(parent = environment(sampling))
  factor_fault <- new.env(parent = environment(.block_covariance_guard_block))
  reached <- 0L
  factor_fault$.covariance_sampling_factor <- function(...) {
    reached <<- reached + 1L
    stop("injected eligible factor failure")
  }
  guard <- .block_covariance_guard_block
  environment(guard) <- factor_fault
  guards$.block_covariance_guard_block <- guard
  for (name in c(".block_covariance_guard_cache", ".block_covariance_guard_plan")) {
    fun <- get(name, envir = environment(sampling))
    environment(fun) <- guards
    assign(name, fun, envir = guards)
  }
  environment(sampling) <- guards
  set.seed(13064)
  expected <- matrix(0, 3L, 35L)
  for (s in 1:3) expected[s, ] <- as.vector(stats::rnorm(35L) %*%
    .covariance_sampling_factor(.covariance_factorization(.block_covariance_dense(fixture$N, s))))
  expected_seed <- .Random.seed
  set.seed(13064)
  expect_equal(sampling(matrix(0, 3L, 35L), fixture$N), expected, tolerance = 0)
  expect_gt(reached, 0L)
  expect_identical(.Random.seed, expected_seed)
})


test_that("large source phases and the small execution boundary retain normal budgets", {

  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_kind)))
  fixture <- .guarded_prediction_fixture(variable = TRUE, copies = 7L)
  for (kind in c("Inversion", "Box-Muller")) for (invalid in c(FALSE, TRUE)) {
    noise <- fixture$N
    if (invalid) noise$values[[1L]][2L, 1L, 2L] <- noise$values[[1L]][2L, 1L, 2L] + .01
    parts <- list(means = matrix(0, 3L, 35L), latent_means = matrix(0, 3L, 35L),
      random_covariance = fixture$G, sampling_covariance = noise,
      covariance = .block_covariance_add(fixture$G, noise), dependency_blocks = fixture$blocks)
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(13065)
    stats::rnorm(if (invalid) 140L else 210L)
    expected_seed <- .Random.seed
    expected_next <- stats::rnorm(5L)
    do.call(RNGkind, list("Mersenne-Twister", kind))
    set.seed(13065)
    if (invalid) expect_error(.predict_joint_selection_source_posterior(parts, rep(1, 35L)),
      "Covariance must be symmetric.", fixed = TRUE)
    else expect_true(all(is.finite(.predict_joint_selection_source_posterior(parts, rep(1, 35L)))))
    expect_identical(.Random.seed, expected_seed)
    expect_identical(stats::rnorm(5L), expected_next)
  }
  # Draw-invariant blocks: sampling is guarded above K = 32 (one plan per draw),
  # solving above the solve threshold (one plan per solve).
  for (K in c(32L, 33L, BLOCK_COVARIANCE_GUARD_MIN_K_SOLVE_INVARIANT - 1L,
              BLOCK_COVARIANCE_GUARD_MIN_K_SOLVE_INVARIANT)) {
    covariance <- diag(K)
    covariance[1:2, 1:2] <- covariance[1:2, 1:2] + .1
    covariance[3:K, 3:K] <- covariance[3:K, 3:K] + .1
    x <- .block_covariance_from_matrix(covariance, 3L, list(1:2, 3:K))
    calls <- 0L
    fault <- new.env(parent = environment(.outcome_rng.norm_known_v_covariance))
    fault$.block_covariance_guard_plan <- function(...) {
      calls <<- calls + 1L
      .block_covariance_guard_plan(...)
    }
    sampling <- .outcome_rng.norm_known_v_covariance
    solve <- .block_covariance_guarded_chol_solve
    environment(sampling) <- environment(solve) <- fault
    sampling(matrix(0, 3L, K), x)
    expect_identical(calls, if (K > 32L) 3L else 0L)
    calls <- 0L
    cache <- if (K >= BLOCK_COVARIANCE_GUARD_MIN_K_SOLVE_INVARIANT) {
      .block_covariance_guard_cache(x, factors = FALSE, pd_only = TRUE)
    } else NULL
    solve(x, 1L, rep(1, K), cache)
    expect_identical(calls, if (K >= BLOCK_COVARIANCE_GUARD_MIN_K_SOLVE_INVARIANT) 1L else 0L)
  }
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
