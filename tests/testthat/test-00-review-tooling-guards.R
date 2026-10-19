test_that("review metadata helpers preserve named fields and reject discarded fields", {

  valid <- with_draw_metadata(1:2, linear_weights = c(mu = 1))
  expect_identical(BayesTools::posterior_metadata(valid, "linear_weights"), c(mu = 1))
  expect_identical(with_draw_metadata(1:2), 1:2)
  malformed <- list(list(c(mu = 1)),
    stats::setNames(list(c(mu = 1), c(mu = 2)), c("linear_weights", "")),
    stats::setNames(list(c(mu = 1), c(mu = 2)), rep("linear_weights", 2L)))
  for (fields in malformed) {
    expect_error(do.call(with_draw_metadata, c(list(x = 1:2, class = NULL), fields)),
      "Draw metadata fields must have nonempty unique names.", fixed = TRUE)
  }
})


test_that("review covariance helpers refuse malformed geometry before draw extraction", {

  expect_error(.as_block_covariance(array(0, c(1L, 2L, 3L))),
    "A fixture covariance must be a draw x row x row array.", fixed = TRUE)
  invalid <- list(list(rho = -.6), list(rho = NA_real_), list(rho = Inf), list(rho = 0),
    list(variance = c(0, .1, .2)), list(variance = c(-.1, .1, .2)),
    list(variance = c(Inf, .1, .2)))
  for (args in invalid) {
    expect_error(do.call(.dense_route_negative_covariance, args),
      "A negative dense-route block needs positive finite variances, at least three rows, and -1/(size - 1) <= rho < 0.",
      fixed = TRUE)
  }
  expect_equal(diag(.dense_route_negative_covariance()), c(.04, .09, .05),
    tolerance = 8 * .Machine$double.eps)
  for (V in list(diag(2L), matrix(.1, 3L, 3L))) {
    expect_error(known_v_pair_conditional_log_lik(NULL, rep(0, nrow(V)), V),
      "A pair-conditional reference requires exactly one partner per row.", fixed = TRUE)
  }
  draws <- cbind(mu_intercept = c(.1, .2), mu__xREx__estimate_intercept = c(.3, .4))
  reference <- list(fit = coda::mcmc.list(coda::mcmc(draws)))
  V <- matrix(c(.04, .02, .02, .09), 2L)
  yi <- c(-.1, .3)
  actual <- known_v_pair_conditional_log_lik(reference, yi, V)
  expected <- vapply(1:2, function(i) {
    j <- 3L - i
    vapply(1:2, function(s) {
      covariance <- V + diag(draws[s, 2L]^2, 2L)
      beta <- covariance[i, j] / covariance[j, j]
      stats::dnorm(yi[[i]], draws[s, 1L] + beta * (yi[[j]] - draws[s, 1L]),
        sqrt(covariance[i, i] - beta * covariance[i, j]), log = TRUE)
    }, numeric(1L))
  }, numeric(2L))
  expect_equal(unname(actual), unname(expected), tolerance = 1e-14)
})
