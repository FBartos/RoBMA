

test_that("inactive multivariate selection ignores a zero weight vector", {

  covariance <- matrix(c(1, .4, .4, 1), 2L)
  omega <- rbind(c(0, 0), c(1, .2), c(0, 0))
  reference_omega <- omega
  reference_omega[c(1L, 3L), ] <- 1
  generate <- function(weights, mode = c(0L, 1L, 0L)) {

    .Call("RoBMA_selnorm_mnorm_step_rng_batch",
      matrix(c(-.2, .3, .1, -.1, .2, .4), 3L, 2L),
      array(rep(covariance, each = 3L), c(3L, 2L, 2L)), c(1, 1), weights,
      c(0, -Inf), c(Inf, 0), 1L, mode, list(1:2), 1000L, 0L,
      PACKAGE = "RoBMA")
  }
  withr::local_seed(73019)
  initial_seed <- .Random.seed
  reference <- generate(reference_omega)
  reference_seed <- .Random.seed
  assign(".Random.seed", initial_seed, envir = .GlobalEnv)
  result <- generate(omega)
  expect_identical(result, reference)
  expect_identical(.Random.seed, reference_seed)
  expect_identical(result$failure_code, 0L)
  expect_true(all(is.finite(result$draws)))
  # The active zero-acceptance model remains rejected.
  rejected <- generate(omega, rep(1L, 3L))
  expect_identical(rejected$failure_code, 3L)
})

test_that("selected response blocks preserve every supplied covariance", {

  covariance <- matrix(c(1, .8, .8, 1), 2L)
  generate <- function(blocks, covariance_value = covariance) {

    .Call("RoBMA_selnorm_mnorm_step_rng_batch", matrix(0, 1L, 2L),
      array(covariance_value, c(1L, 2L, 2L)), c(1, 1), matrix(1, 1L, 1L),
      -Inf, Inf, 1L, 0L, blocks, 1L, 0L, PACKAGE = "RoBMA")
  }
  withr::local_seed(73021)
  initial_seed <- .Random.seed
  expect_error(generate(list(1L, 2L)),
    "'dependency_blocks' must include all observations with nonzero covariance.", fixed = TRUE)
  expect_identical(.Random.seed, initial_seed)
  expect_error(generate(list(c(NA_integer_, 2L))),
    "'dependency_blocks' must partition the observations.", fixed = TRUE)
  expect_identical(.Random.seed, initial_seed)
  expect_identical(generate(list(1:2))$failure_code, 0L)
  expect_identical(generate(list(1L, 2L), diag(2))$failure_code, 0L)
})

.native_srs_draws <- function(covariance, samples = 8L) {

  K <- nrow(covariance)
  .Call("RoBMA_selnorm_mnorm_step_rng_batch", matrix(0, samples, K),
    array(rep(covariance, each = samples), c(samples, K, K)), rep(1, K),
    matrix(1, samples, 1L), -Inf, Inf, 1L, 0L, list(seq_len(K)), 1L, 0L,
    PACKAGE = "RoBMA")
}

test_that("native recognized rank-one sampling retains its original first-normal path", {

  loading <- c(.7, .2)
  covariance <- tcrossprod(loading)
  # Rounded products can leave a positive Cholesky pivot despite the established
  # factor-recognition equality. This case distinguishes the branch priority.
  withr::local_seed(73027)
  result <- .native_srs_draws(covariance)
  observed_seed <- .Random.seed
  set.seed(73027)
  standard <- matrix(stats::rnorm(8L * 2L), ncol = 2L, byrow = TRUE)
  expected <- outer(standard[, 1L], loading)
  expect_identical(result$failure_code, 0L)
  expect_identical(unname(result$draws), unname(expected))
  expect_identical(.Random.seed, observed_seed)
})

test_that("native spectral fallback preserves marginal scales and structural zeros", {

  withr::local_seed(73029)
  small <- .native_srs_draws(diag(c(1, 1e-20, 0)))
  small_seed <- .Random.seed
  set.seed(73029)
  unit <- .native_srs_draws(diag(c(1, 1, 0)))
  expect_identical(small$failure_code, 0L)
  expect_identical(unit$failure_code, 0L)
  expect_identical(.Random.seed, small_seed)
  expect_identical(small$draws[, 1L], unit$draws[, 1L])
  expect_equal(small$draws[, 2L] / sqrt(1e-20), unit$draws[, 2L],
    tolerance = 8 * .Machine$double.eps)
  expect_true(all(small$draws[, 2L] != 0))
  expect_identical(as.numeric(small$draws[, 3L]), rep(0, 8L))
})

test_that("native SRS fallback keeps rank-two Gaussian support and draw count", {

  # Exact dyadic Gram matrix of columns (1,1,1,1) and (.5,-.5,.5,-.5).
  # Its exact eigenvalues are4,1,0,0; rows1/3 and2/4 are equal as random variables.
  covariance <- matrix(c(1.25, .75, 1.25, .75,
    .75, 1.25, .75, 1.25, 1.25, .75, 1.25, .75,
    .75, 1.25, .75, 1.25), 4L)
  withr::local_seed(73031)
  result <- .native_srs_draws(covariance)
  observed_seed <- .Random.seed
  expect_identical(result$failure_code, 0L)
  bound <- 64 * .Machine$double.eps * max(abs(result$draws))
  expect_lte(max(abs(result$draws[, 1L] - result$draws[, 3L])), bound)
  expect_lte(max(abs(result$draws[, 2L] - result$draws[, 4L])), bound)
  expect_equal(qr(result$draws)$rank, 2L)
  set.seed(73031)
  invisible(stats::rnorm(8L * 4L))
  expect_identical(.Random.seed, observed_seed)
})

test_that("nonfinite standardized covariance fails before spectral evaluation", {

  withr::local_seed(73033)
  initial_seed <- .Random.seed
  invalid <- matrix(c(1e-308, 1e308, 1e308, 1e-308), 2L)
  result <- .native_srs_draws(invalid, samples = 1L)
  expect_identical(result$failure_code, 1L)
  expect_identical(.Random.seed, initial_seed)
})

test_that("JAGS selection controls reject fractional modes and bins", {
  skip_if_not_installed("rjags")


  base <- list(y = c(0, 0), mu = c(0, 0), covariance = c(1, 0, 1),
    sd = c(1, 1), loading = c(.1, .1), sei = c(1, 1), omega = c(1, .2),
    lower = c(0, -1e300), upper = c(1e300, 0), bins = c(1, 1),
    sign = 1, telescope = 1, mode = 0, points = 4, maximum = 4,
    scrambles = 2, tolerance = .005, rule = 0,
    nodes = rep(0, 9), log_weights = rep(0, 9), orders = c(1, 3, 5))
  cases <- list(
    dense = list(distribution = "dselnorm_mnorm_step", arguments = c(
      "mu", "covariance", "sei", "omega", "lower", "upper", "bins",
      "sign", "telescope", "mode", "qmc", "points", "scrambles", "tolerance",
      "rule", "nodes", "log_weights", "orders"),
      extra = list(qmc = rep(.5, 32)), bad = c("mode", "sign", "telescope", "bins", "points")),
    cluster = list(distribution = "dselnorm_cluster_step", arguments = c(
      "mu", "sd", "loading", "sei", "omega", "lower", "upper", "bins",
      "sign", "telescope", "mode", "nodes", "log_weights", "orders", "qmc",
      "points", "maximum", "scrambles", "tolerance", "rule"),
      extra = list(qmc = rep(.5, 16)), bad = c("mode", "sign", "telescope", "bins")),
    factor = list(distribution = "dselnorm_factor_step", arguments = c(
      "mu", "sd", "factor_loading", "sei", "omega", "lower", "upper", "bins",
      "sign", "telescope", "mode", "nodes", "log_weights", "orders", "rule_count", "qmc",
      "points", "maximum", "scrambles", "tolerance", "rule"),
      extra = list(factor_loading = rep(.1, 4), rule_count = 3, qmc = rep(.5, 32)),
      bad = c("mode", "bins", "points", "rule_count")),
    conditioned = list(distribution = "dselnorm_sampling_conditioned", arguments = c(
      "mu", "covariance", "residual", "loading", "rank", "mu", "mu", "sei", "omega",
      "lower", "upper", "bins", "sign", "mode", "telescope", "rule", "groups", "qmc",
      "points", "maximum", "scrambles", "tolerance", "nodes", "log_weights", "orders",
      "nodes", "log_weights", "orders", "rule_count"),
      extra = list(residual = c(0, 0), rank = 1, groups = c(1, 2), qmc = rep(.5, 32),
        rule_count = rep(3, 3)), bad = c("rank", "bins", "groups", "rule", "scrambles"))
  )
  compile <- function(specification, data) {

    used <- unique(c("y", specification$arguments))
    data <- data[used]
    syntax <- paste0("model { y[1:2] ~ ", specification$distribution, "(",
      paste(specification$arguments, collapse = ", "), ") }")
    connection <- textConnection(syntax)
    on.exit(close(connection), add = TRUE)
    model <- rjags::jags.model(connection, data = data,
      n.chains = 1L, n.adapt = 0L, quiet = TRUE)
    invisible(TRUE)
  }
  for (name in names(cases)) {
    specification <- cases[[name]]
    data <- utils::modifyList(base, specification$extra)
    expect_true(isTRUE(compile(specification, data)))
    for (field in intersect(c("mode", "bins"), specification$bad)) {
      changed <- data
      # mode0.5 would formerly have been silently interpreted as NORMAL0.
      changed[[field]][1L] <- changed[[field]][1L] + .5
      error <- tryCatch(compile(specification, changed), error = identity)
      expect_s3_class(error, "error")
      expect_match(conditionMessage(error), "Invalid parent values", fixed = TRUE)
    }
  }
  invisible(TRUE)
})

test_that("native covariance ingress preserves explicit structural failures", {

  withr::local_seed(73037)
  asymmetric <- matrix(c(1, .4, .4, 1), 2L)
  asymmetric[1L, 2L] <- asymmetric[1L, 2L] + .Machine$double.eps
  result <- .native_srs_draws(asymmetric, samples = 1L)
  expect_identical(result$failure_code, 4L)
  expect_identical(result$failure_size, 2L)
  for (covariance in list(diag(c(1, -1e-20, 0)),
      matrix(c(0, 1e-20, 1e-20, 1), 2L))) {
    expect_identical(.native_srs_draws(covariance, samples = 1L)$failure_code, 1L)
  }
})

# Append to test-00-selection-native-inputs.R, reusing its .native_srs_draws.
test_that("native Gaussian sampling preserves exact integer covariance support", {

  # Exact Gram of rows(1,1),(1,0),(0,1): rank2 and Y1-Y2-Y3=0.
  covariance <- matrix(c(2, 1, 1, 1, 1, 0, 1, 0, 1), 3L)
  null <- c(1, -1, -1)
  permutations <- list(1:3, c(1L, 3L, 2L), c(2L, 1L, 3L),
    c(2L, 3L, 1L), c(3L, 1L, 2L), 3:1)
  withr::local_seed(73037)
  for (index in permutations) {
    result <- .native_srs_draws(covariance[index, index, drop = FALSE])
    expect_identical(result$failure_code, 0L)
    residual <- as.numeric(result$draws %*% null[index])
    rounding_bound <- 64 * .Machine$double.eps *
      max(abs(result$draws)) * sum(abs(null))
    expect_lte(mean(residual^2), rounding_bound^2)
    expect_equal(qr(result$draws)$rank, 2L)
  }
})
