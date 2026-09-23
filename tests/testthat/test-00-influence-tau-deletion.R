test_that("heterogeneity deletion preserves small remaining scale rows", {

  tau <- rbind(c(1, 1e-10, 2e-10), c(2, 3e-10, 4e-10))
  weights <- cbind(c(.25, .75), c(.6, .4), c(.8, .2))
  expected <- vapply(seq_len(ncol(tau)), function(i) {
    remaining <- tau[, -i, drop = FALSE]
    sum(weights[, i] * sqrt(rowMeans(remaining^2)))
  }, numeric(1))

  observed <- .influence_tau_del_from_samples(tau, weights)
  expect_equal(observed / expected, rep(1, ncol(tau)), tolerance = 1e-14)

  # The deletion target is independent of observation ordering.
  order <- c(3L, 1L, 2L)
  permuted <- .influence_tau_del_from_samples(tau[, order], weights[, order])
  expect_equal(permuted / expected[order], rep(1, ncol(tau)), tolerance = 1e-14)
})

test_that("heterogeneity deletion is stable across finite scale units", {

  tau <- rbind(c(1, 2, 3), c(2, 3, 4), c(0, 0, 0))
  weights <- cbind(c(.2, .3, .5), c(.3, .5, .2), c(.5, .2, .3))
  expected <- vapply(seq_len(ncol(tau)), function(i) {
    sum(weights[, i] * sqrt(rowMeans(tau[, -i, drop = FALSE]^2)))
  }, numeric(1))

  for (scale in c(1e-200, 1e200, 1e307)) {
    observed <- .influence_tau_del_from_samples(tau * scale, weights)
    expect_equal(observed / scale, expected, tolerance = 1e-14)
  }

  expect_equal(
    .influence_tau_del_from_samples(tau[, 1L, drop = FALSE], weights[, 1L, drop = FALSE]),
    sum(tau[, 1L] * weights[, 1L])
  )
  expect_equal(
    .influence_tau_del_from_samples(tau[, 1:2], weights[, 1:2]),
    c(sum(tau[, 2L] * weights[, 1L]), sum(tau[, 1L] * weights[, 2L]))
  )
})
