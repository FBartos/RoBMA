test_that("residual standard deviations retain small variances", {

  expect_identical(
    .hat_variance_sd(c(0, .Machine$double.xmin), "Test variance"),
    sqrt(c(0, .Machine$double.xmin))
  )
  expect_error(
    .hat_variance_sd(-.Machine$double.xmin, "Test variance"),
    "finite and non-negative"
  )
})


test_that("unit-leverage residual rows are identified structurally", {

  X <- cbind(1, c(1, 0, 0))

  expect_identical(.hat_zero_residual_rows(X), 1L)

  nearly_dependent <- rbind(
    c(0, 1),
    c(1, 1),
    c(2, 2 + 1e-10)
  )
  expect_identical(.hat_zero_residual_rows(nearly_dependent), integer())

  alloc <- factor(c(
    rep("random", 4), "alternate", "alternate", rep("random", 3),
    rep("systematic", 4)
  ))
  year <- factor(c(
    rep("before", 3), "after", "after", "before", "after", "after",
    "before", "before", rep("after", 3)
  ))
  factorial <- stats::model.matrix(~ alloc * year)

  expect_identical(.hat_zero_residual_rows(factorial), c(5L, 6L, 10L))
})


test_that("transformed residual variances use a covariance factor", {

  X          <- cbind(1, c(-1, 0, 1))
  diagonal   <- c(0.5, 1, 2)
  rank_one   <- c(0.2, 0.3, 0.4)
  block      <- list(seq_len(3L))
  covariance <- diag(diagonal) + tcrossprod(rank_one)
  precision  <- solve(covariance)
  WX         <- precision %*% X
  XB         <- X %*% solve(crossprod(X, WX))
  A          <- diag(3L) - XB %*% t(WX)

  variance <- .hat_transformed_covariance_diag(
    transform     = A,
    diagonal      = diagonal,
    rank_one      = rank_one,
    block_indices = block
  )

  expect_true(all(variance >= 0))
  expect_equal(variance, diag(A %*% covariance %*% t(A)), tolerance = 1e-15)
})


test_that("estimate and cluster residual variances avoid PSD subtraction", {

  X <- matrix(c(
    1, 1, 1, 1,
    6.072075e-7, -3.901161e5, 4.130463, 104.9131,
    -5.899230e-6, -1.276173e5, 4.659023e-6, 1270.221
  ), nrow = 4L)
  diagonal <- c(4.240359, 2.437391e-6, 1.313608, 6.024016)
  rank_one <- c(-67.93884, 844.9597, .02247927, -2779.021)
  blocks   <- list(seq_len(4L))
  W        <- .hat_precision_matrix(diagonal, rank_one, blocks)
  WX       <- W %*% X
  B        <- .hat_solve_crossprod(crossprod(X, WX))
  A        <- diag(4L) - X %*% B %*% t(WX)

  old_estimate <- diag(W) - rowSums((WX %*% B) * WX)
  estimate_transform <- W - WX %*% B %*% t(WX)
  estimate_variance <- .hat_transformed_covariance_diag(
    transform     = estimate_transform,
    diagonal      = diagonal,
    rank_one      = rank_one,
    block_indices = blocks
  )

  between_covariance <- tcrossprod(rank_one)
  cluster_adjustment <- diag(4L) - between_covariance %*% W
  cluster_variance <- .hat_transformed_covariance_diag(
    transform     = cluster_adjustment %*% A,
    diagonal      = diagonal,
    rank_one      = rank_one,
    block_indices = blocks
  )

  expect_true(any(old_estimate < 0))
  expect_true(all(estimate_variance >= 0))
  expect_true(all(cluster_variance >= 0))
})


test_that("known-V estimate residual variance avoids PSD subtraction", {

  X <- matrix(c(
    1, 1, 1, 1,
    6.072075e-7, -3.901161e5, 4.130463, 104.9131,
    -5.899230e-6, -1.276173e5, 4.659023e-6, 1270.221
  ), nrow = 4L)
  sampling_variance <- c(
    151.324100309942,
    5463391.08469307,
    1.54136201531135e-8,
    1.8986576975829e-4
  )
  extra_variance <- c(
    .0931437304737963,
    186841.573373015,
    861201.17955025,
    1.0401676712931e-4
  )
  known_V <- .known_v_prepare(
    V                         = sampling_variance,
    keep_rows                 = rep(TRUE, 4L),
    known_v_parameterization  = "auto",
    known_v_residual_fraction = NULL,
    warn_singular             = FALSE
  )

  projection <- .known_v_gls_projection_blocks(
    X              = X,
    y              = c(-.2, .1, .3, -.1),
    known_V        = known_V,
    extra_variance = extra_variance
  )

  marginal_variance <- sampling_variance + extra_variance
  W                 <- diag(1 / marginal_variance)
  WX                <- W %*% X
  B                 <- .hat_solve_crossprod(crossprod(X, WX))
  V_W_X             <- sampling_variance * WX
  subtractive <- sampling_variance^2 / marginal_variance -
    rowSums((V_W_X %*% B) * V_W_X)

  expect_true(any(subtractive < 0))
  expect_true(all(projection[["sampling_residual_variance"]] >= 0))
})


test_that("known-V rank-one residuals allow a regularized zero diagonal", {

  V <- tcrossprod(c(.2, .3, .4))
  known_V <- .known_v_canonicalize(V, warn_singular = FALSE)

  projection <- .known_v_gls_projection_blocks(
    X              = cbind(1, c(-1, 0, 1)),
    y              = c(-.1, .2, .3),
    known_V        = known_V,
    extra_variance = c(0, .1, .2)
  )

  expect_true(all(is.finite(unlist(projection))))
  expect_true(all(projection[["sampling_residual_variance"]] >= 0))
})


test_that("hat-matrix diagnostics retain small nonzero residuals", {

  tiny <- 2^-100
  object <- brma.norm(
    yi        = c(0, tiny, 0),
    sei       = rep(1, 3),
    only_data = TRUE
  )
  object[["priors"]] <- list(
    outcome = list(tau = NULL, bias = NULL),
    scale   = NULL
  )
  object[["fit"]] <- coda::mcmc.list(coda::mcmc(matrix(
    0,
    nrow     = 1,
    ncol     = 1,
    dimnames = list(NULL, "tau")
  )))
  class(object) <- c("brma.norm", "brma")

  diagnostics <- .compute_hat_matrix_samples(
    object       = object,
    return_se    = TRUE,
    return_resid = TRUE
  )

  expect_true(all(diagnostics[["resid"]] != 0))
  expect_equal(
    as.vector(diagnostics[["resid"]]),
    c(-tiny / 3, 2 * tiny / 3, -tiny / 3),
    tolerance = 1e-15
  )
  expect_equal(
    as.vector(diagnostics[["se"]]),
    rep(sqrt(2 / 3), 3),
    tolerance = 1e-15
  )
})


test_that("LOO-PIT transformation preserves extreme normal quantiles", {

  z <- c(-100, -40, 0, 40, 100)
  log_lower <- matrix(stats::pnorm(z, log.p = TRUE), nrow = 1L)
  log_upper <- matrix(stats::pnorm(
    z,
    lower.tail = FALSE,
    log.p      = TRUE
  ), nrow = 1L)

  expect_equal(
    .loo_pit_z_from_log_tails(
      log_lower    = log_lower,
      log_upper    = log_upper,
      psis_weights = matrix(1, nrow = 1L, ncol = length(z))
    ),
    z,
    tolerance = 1e-12
  )
})


test_that("LOO-PIT tail mixtures are invariant to weight scale", {

  z <- c(-8, 8)
  log_lower <- matrix(
    rep(stats::pnorm(z, log.p = TRUE), each = 2L),
    nrow = 2L
  )
  log_upper <- matrix(
    rep(stats::pnorm(z, lower.tail = FALSE, log.p = TRUE), each = 2L),
    nrow = 2L
  )

  expect_equal(
    .loo_pit_z_from_log_tails(
      log_lower    = log_lower,
      log_upper    = log_upper,
      psis_weights = matrix(c(2, 3, 20, 30), nrow = 2L)
    ),
    z,
    tolerance = 1e-12
  )
})


test_that("normal LOO-PIT tails remain finite after probability underflow", {

  setup <- list(
    outcome_type = "norm",
    data         = list(),
    yi           = c(-100, 100),
    sei          = c(1, 1),
    S            = 2L,
    K            = 2L,
    mu           = matrix(0, nrow = 2L, ncol = 2L),
    tau_within   = matrix(0, nrow = 2L, ncol = 2L)
  )
  tails <- .loo_predictive_log_tails_estimate(setup)

  expect_true(all(is.finite(tails[["log_lower"]])))
  expect_true(all(is.finite(tails[["log_upper"]])))
  expect_equal(
    .loo_pit_z_from_log_tails(
      tails[["log_lower"]],
      tails[["log_upper"]],
      matrix(0.5, nrow = 2L, ncol = 2L)
    ),
    c(-100, 100),
    tolerance = 1e-12
  )
})


test_that("GLMM LOO-PIT is rejected before PSIS work", {

  expect_error(
    .check_residual_type_availability("LOO-PIT", "bin", FALSE),
    "discrete PIT convention"
  )
  expect_error(
    .loo_predictive_log_tails_estimate(list(outcome_type = "pois")),
    "discrete PIT convention"
  )
})


test_that("predictive variance uses centered mixture moments", {

  mean <- matrix(1e12 + c(-1, 1), ncol = 1L)
  moments <- .weighted_predictive_moments(
    mean     = mean,
    variance = matrix(4, nrow = 2L),
    weights  = matrix(c(1, 1), ncol = 1L)
  )

  expect_equal(moments[["mean"]], 1e12)
  expect_equal(moments[["variance"]], 5)
  expect_equal(moments[["se"]], sqrt(5))
  expect_error(
    .weighted_predictive_moments(
      mean     = mean,
      variance = matrix(c(1, -1), ncol = 1L),
      weights  = matrix(c(1, 1), ncol = 1L)
    ),
    "invalid"
  )
})
