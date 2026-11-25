test_that("factor proposals retain the complete Gaussian product-selection law", {

  S <- 600L
  K <- 8L
  mu <- -.5
  rho <- .1
  lower_weight <- .1
  residual_sd <- sqrt(1 - rho)
  loading <- sqrt(rho)
  # Conditional independence given one normal factor gives a separate
  # one-dimensional reference for all eight selection events jointly.
  integrand <- function(factor, moment) {
    location <- mu + loading * factor
    score <- location / residual_sd
    weight <- lower_weight + (1 - lower_weight) * stats::pnorm(score)
    first <- location * weight +
      (1 - lower_weight) * residual_sd * stats::dnorm(score)
    value <- switch(moment,
      mass = weight^K,
      first = first * weight^(K - 1L),
      second = ((location^2 + residual_sd^2) * weight +
        (1 - lower_weight) * residual_sd * location * stats::dnorm(score)) * weight^(K - 1L),
      cross = first^2 * weight^(K - 2L)
    )
    stats::dnorm(factor) * value
  }
  integral <- function(moment) stats::integrate(integrand, -Inf, Inf,
    moment = moment, rel.tol = 1e-10, abs.tol = 1e-14)$value
  mass <- integral("mass")
  expected <- integral("first") / mass
  variance <- integral("second") / mass - expected^2
  covariance <- integral("cross") / mass - expected^2
  mean_se <- sqrt((variance + (K - 1L) * covariance) / (S * K))

  covariance_matrix <- matrix(rho, K, K)
  diag(covariance_matrix) <- 1
  selection_context <- list(
    omega = matrix(c(1, lower_weight), S, 2L, byrow = TRUE),
    kernel_mode = rep(SELKERNEL_STEP, S), vector_rule = 0L,
    use_normal = rep(FALSE, S), p_cuts = c(0, .5, 1), sign = 1L
  )
  withr::local_seed(365)
  draws <- .outcome_rng.selnorm_mvn(
    matrix(mu, S, K), array(rep(covariance_matrix, each = S), c(S, K, K)),
    rep(1, K), selection_context, list(seq_len(K))
  )
  expect_true(all(is.finite(draws)))
  expect_lt(abs(mean(draws) - expected), 6 * mean_se)
})
