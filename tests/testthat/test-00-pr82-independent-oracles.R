test_that("known-V BLUP matches direct Gaussian conditioning with heterogeneous scales", {

  V <- matrix(c(.09, .018, -.006, .018, .16, .012, -.006, .012, .04), 3L)
  yi <- c(.2, -.1, .4)
  mu <- rbind(c(.1, .1, .1), c(-.2, .05, .15), c(0, 0, 0))
  tau <- rbind(c(.3, .1, .2), c(0, .2, .4), c(0, 0, 0))
  bias <- rbind(c(.02, -.03, .01), c(.01, .01, .01), c(0, 0, 0))
  expected <- mu
  for (draw in seq_len(nrow(mu))) {
    Q <- diag(tau[draw, ]^2)
    expected[draw, ] <- mu[draw, ] + as.vector(
      Q %*% solve(Q + V, yi - mu[draw, ] - bias[draw, ])
    )
  }
  actual <- .evaluate.brma.known_v_blup.norm(
    mu_samples = mu, tau_within = tau, yi = yi,
    known_V = .known_v_canonicalize(V), bias_offset = bias
  )
  expect_lt(max(abs(actual - expected)), 1e-12)
  expect_identical(unname(actual[3L, ]), rep(0, 3L))
})


test_that("working VIF covariance matches independent weighted least-squares solves", {

  X <- cbind(intercept = 1, x = c(-2, -1, 0, 1, 2), z = c(0, 1, -1, 1, 0))
  vi <- c(.1, .2, .15, .08, .12)
  weights <- c(1, .5, 2, 1.2, .8)
  tau <- rbind(rep(0, 5L), rep(.3, 5L), c(.1, .2, .3, .4, .5))
  expected <- Reduce(`+`, lapply(seq_len(nrow(tau)), function(draw) {
    W <- diag(weights / (vi + tau[draw, ]^2))
    solve(t(X) %*% W %*% X)
  })) / nrow(tau)
  actual <- .vif_vcov_from_tau_samples(X, vi, weights, tau)
  expect_lt(max(abs(actual - expected)) / max(abs(expected)), 1e-12)
})
