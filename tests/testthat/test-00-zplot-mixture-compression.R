test_that("compressed mixed-sign moments retain the finite-rule selected density", {

  sei <- c(.7, 1.3)
  covariance <- diag(.6 * sei^2) + .4 * tcrossprod(sei)
  residual_sd <- sqrt(.6) * sei
  loading <- sqrt(.4) * sei
  means <- cbind(c(-1.2, -.4, 0, .3, 1.1), c(.8, -.7, 0, .5, -1.3))
  context_weight <- c(1, .8, .3, 1e-4, 1e-12)
  context_weight <- context_weight / sum(context_weight)
  rule <- .gauss_hermite_nodes(31L)
  attr(rule$nodes, "bounded_cdf") <- FALSE
  z <- seq(-8, 10, length.out = 321L)
  expected <- numeric(length(z))
  normalizers <- numeric(nrow(means))
  for (context in seq_len(nrow(means))) {
    location <- outer(rule$nodes, loading) +
      matrix(means[context, ], nrow = length(rule$nodes), ncol = 2L, byrow = TRUE)
    acceptance <- .2 + .8 * stats::pnorm(sweep(location, 2L, residual_sd, "/"))
    normalizers[context] <- sum(rule$weights * acceptance[, 1L] * acceptance[, 2L])
    for (row in seq_along(sei)) {
      contribution <- vapply(z, function(point) {
        sum(rule$weights * acceptance[, 3L - row] *
          stats::dnorm(point * sei[row], location[, row], residual_sd[row])) * sei[row]
      }, 0)
      expected <- expected + context_weight[context] * contribution /
        (2 * normalizers[context])
    }
  }
  expected <- expected * ifelse(z >= 0, 1, .2)
  actual <- .Call("RoBMA_selnorm_context_star_zplot", means,
    as.double(covariance[lower.tri(covariance, diag = TRUE)]), sei,
    c(1, .2), c(0, -Inf), c(Inf, 0), 1L, FALSE, log(normalizers),
    log(context_weight), rule$nodes, rule$log_weights, z, 1e-6, PACKAGE = "RoBMA")
  expect_true(actual$available)
  expect_gt(as.double(actual$compression_error), 0)
  expect_lte(as.double(actual$compression_error), 1e-6)
  expect_lte(max(abs(actual$density - expected)),
    as.double(actual$compression_error) + 1e-12)
  expect_equal(actual$mass, 1, tolerance = 1e-12)
  expect_equal(as.double(actual$compact_log_normalizers), log(normalizers), tolerance = 1e-12)
  expect_identical(actual$mixture_components,
    as.double(nrow(means) * ncol(means) * length(rule$nodes)))
})


test_that("singleton Gaussians retain their scalar density after recurrence refusal", {

  # Repeated central nodes form one constant group; the two distinct endpoints
  # remain substantial singleton components. This is a finite-mixture oracle,
  # not a claim that this artificial rule integrates the continuous model.
  nodes <- c(rep(0, 29L), -8, 8)
  weights <- c(rep(.5 / 29, 29L), .25, .25)
  attr(nodes, "bounded_cdf") <- FALSE
  location <- sqrt(.3) * nodes
  sd <- sqrt(.7)
  acceptance <- .2 + .8 * stats::pnorm(location / sd)
  normalizer <- sum(weights * acceptance^2)
  regular <- seq(-10, 10, length.out = 321L)
  duplicate <- regular
  duplicate[20L] <- duplicate[19L]
  irregular <- regular + seq_along(regular)^2 * 1e-6
  for (z in list(regular, rev(regular), duplicate, irregular,
                 seq(-160, 160, length.out = 321L))) {
    expected <- vapply(z, function(point) {
      sum(weights * acceptance * stats::dnorm(point, location, sd)) *
        ifelse(point >= 0, 1, .2) / normalizer
    }, 0)
    actual <- .Call("RoBMA_selnorm_context_star_zplot", matrix(0, 1L, 2L),
      c(1, .3, 1), c(1, 1), c(1, .2), c(0, -Inf), c(Inf, 0),
      1L, FALSE, log(normalizer), 0, nodes, log(weights), z, 1e-6,
      PACKAGE = "RoBMA")
    expect_true(actual$available)
    expect_lte(max(abs(actual$density - expected)),
      as.double(actual$compression_error) + 1e-12)
    expect_equal(actual$mass, 1, tolerance = 1e-12)
    expect_equal(as.double(actual$compact_log_normalizers), log(normalizer), tolerance = 1e-12)
    expect_identical(actual$mixture_components, 62)
  }
})
