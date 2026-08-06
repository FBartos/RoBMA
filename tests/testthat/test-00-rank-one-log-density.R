# test-00-rank-one-log-density.R

.aligned_rank_one_oracle <- function(diagonal, rank_one, residual_scale) {

  size        <- length(rank_one)
  rank_energy <- sum(rank_one^2)
  log_det     <- (size - 1L) * log(diagonal) +
    log(diagonal + rank_energy)
  quadratic   <- residual_scale^2 * rank_energy /
    (diagonal + rank_energy)

  list(
    log_det     = log_det,
    quadratic   = quadratic,
    log_density = -0.5 * (size * log(2 * pi) + log_det + quadratic)
  )
}


test_that("diagonal rank-one log density is stable for aligned residuals", {

  rank_one <- c(1, -2, 3)
  cases    <- expand.grid(
    diagonal       = c(1e-15, 1e-18),
    residual_scale = c(1e-4, 1e8)
  )

  for (i in seq_len(nrow(cases))) {
    diagonal       <- cases[["diagonal"]][[i]]
    residual_scale <- cases[["residual_scale"]][[i]]
    oracle          <- .aligned_rank_one_oracle(
      diagonal       = diagonal,
      rank_one       = rank_one,
      residual_scale = residual_scale
    )
    log_density     <- .log_dmvnorm_diag_rank_one(
      x        = residual_scale * rank_one,
      mean     = rep(0, length(rank_one)),
      diagonal = rep(diagonal, length(rank_one)),
      rank_one = rank_one
    )
    quadratic       <- -2 * log_density -
      length(rank_one) * log(2 * pi) - oracle[["log_det"]]

    expect_true(is.finite(log_density))
    expect_gte(quadratic, 0)
    expect_lt(
      abs(quadratic - oracle[["quadratic"]]),
      max(1e-12, 1e-11 * oracle[["quadratic"]])
    )
    expect_equal(log_density, oracle[["log_density"]], tolerance = 1e-11)
  }
})


test_that("vectorized cluster log density is stable for aligned residuals", {

  rank_one       <- c(1, -2, 3)
  diagonal       <- c(1e-15, 1e-18, 1e-15, 1e-18)
  residual_scale <- c(1e-4, 1e-4, 1e8, 1e8)
  sample_size    <- length(diagonal)
  estimate_size  <- length(rank_one)
  setup <- list(
    S           = sample_size,
    mu          = -residual_scale *
      matrix(rank_one, nrow = sample_size, ncol = estimate_size, byrow = TRUE),
    tau_within  = matrix(
      rep(sqrt(diagonal), each = estimate_size),
      nrow = sample_size,
      byrow = TRUE
    ),
    tau_between = matrix(
      rank_one,
      nrow = sample_size,
      ncol = estimate_size,
      byrow = TRUE
    ),
    cluster     = list(seq_len(estimate_size))
  )

  log_density <- as.vector(.log_lik_cluster_norm_analytic(
    setup = setup,
    yi    = rep(0, estimate_size),
    vi    = rep(0, estimate_size)
  ))
  oracle <- lapply(seq_len(sample_size), function(i) {
    .aligned_rank_one_oracle(
      diagonal       = diagonal[[i]],
      rank_one       = rank_one,
      residual_scale = residual_scale[[i]]
    )
  })
  oracle_log_density <- vapply(
    oracle,
    function(x) x[["log_density"]],
    numeric(1)
  )

  expect_true(all(is.finite(log_density)))
  expect_equal(log_density, oracle_log_density, tolerance = 1e-11)

  quadratic <- vapply(seq_len(sample_size), function(i) {
    -2 * log_density[[i]] - estimate_size * log(2 * pi) -
      oracle[[i]][["log_det"]]
  }, numeric(1))
  oracle_quadratic <- vapply(
    oracle,
    function(x) x[["quadratic"]],
    numeric(1)
  )
  tolerance <- pmax(1e-12, 1e-11 * oracle_quadratic)

  expect_true(all(quadratic >= 0))
  expect_true(all(abs(quadratic - oracle_quadratic) < tolerance))
})
