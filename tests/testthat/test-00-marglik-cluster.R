test_that("cluster bridge likelihood equals independent Gaussian integration", {

  yi          <- c(0.15, -0.10, 0.40, 0.05, 0.30)
  mu          <- c(0.05, 0.02, 0.25, 0.10, 0.20)
  tau_within  <- c(0.20, 0.25, 0.15, 0.30, 0.18)
  tau_between <- c(0.35, 0.30, 0.25, 0.40, 0.28)
  sei         <- c(0.10, 0.12, 0.08, 0.15, 0.11)
  cluster     <- c(1L, 1L, 2L, 2L, 2L)
  weights     <- c(0.75, 1.25, 0.50, 1.50, 1.00)

  observed <- .marglik_cluster_norm_log_lik(
    yi, mu, tau_within, tau_between, sei, cluster, weights
  )
  reference <- sum(vapply(split(seq_along(yi), cluster), function(index) {
    integrand <- function(gamma) {
      vapply(gamma, function(value) {
        conditional <- stats::dnorm(
          yi[index],
          mean = mu[index] + value * tau_between[index],
          sd   = sqrt(sei[index]^2 + tau_within[index]^2),
          log  = TRUE
        )
        exp(sum(weights[index] * conditional)) * stats::dnorm(value)
      }, numeric(1))
    }
    log(stats::integrate(
      integrand,
      lower     = -Inf,
      upper     = Inf,
      rel.tol   = 1e-12,
      stop.on.error = TRUE
    )[["value"]])
  }, numeric(1)))

  expect_equal(observed, reference, tolerance = 2e-11)
})


test_that("cluster bridge likelihood preserves allocation boundaries", {

  yi      <- c(-0.20, 0.10, 0.35, 0.50)
  mu      <- rep(0.15, 4L)
  sei     <- c(0.10, 0.20, 0.15, 0.12)
  cluster <- c(1L, 1L, 2L, 2L)
  tau     <- 0.30

  rho_zero <- .marglik_cluster_norm_log_lik(
    yi, mu,
    tau_within  = rep(tau, 4L),
    tau_between = rep(0, 4L),
    sei, cluster
  )
  independent <- sum(stats::dnorm(
    yi,
    mean = mu,
    sd   = sqrt(sei^2 + tau^2),
    log  = TRUE
  ))
  expect_equal(rho_zero, independent, tolerance = 1e-14)

  rho_one <- .marglik_cluster_norm_log_lik(
    yi, mu,
    tau_within  = rep(0, 4L),
    tau_between = rep(tau, 4L),
    sei, cluster
  )
  covariance <- matrix(0, nrow = 4L, ncol = 4L)
  covariance[1:2, 1:2] <- tau^2
  covariance[3:4, 3:4] <- tau^2
  diag(covariance) <- diag(covariance) + sei^2
  reference <- mvtnorm::dmvnorm(yi, mean = mu, sigma = covariance, log = TRUE)
  expect_equal(rho_one, unname(reference), tolerance = 1e-12)
})


test_that("cluster bridge likelihood subtracts no explained quadratic", {

  # Two rows with the sampling variance d share one unit latent effect and the
  # residual (1, 1) lies along its loading. In the eigenbasis (1, 1) / sqrt(2)
  # and (1, -1) / sqrt(2) the variances are 2 + d and d, so
  #   log L = -log(2 pi) - (log(d) + log(2 + d)) / 2 - 1 / (2 + d).
  # Subtracting the explained quadratic from the total lost 0.5 at d = 1e-16.
  for (d in c(1e-8, 1e-12, 1e-14, 1e-16, 1e-18)) {
    closed_form <- -log(2 * pi) - (log(d) + log(2 + d)) / 2 - 1 / (2 + d)
    observed <- .marglik_cluster_norm_log_lik(
      yi = c(1, 1), mu = c(0, 0), tau_within = c(0, 0), tau_between = c(1, 1),
      sei = rep(sqrt(d), 2L), cluster = c(1L, 1L)
    )
    expect_equal(observed, closed_form, tolerance = 1e-14)
  }

  # A residual of 10 against the sd 2^-30 makes the total quadratic 2.3e20
  # while the explained part leaves 100: the subtraction returned 18.61.
  d <- 2^-60
  closed_form <- -log(2 * pi) - (log(d) + log(2 + d)) / 2 - 100 / (2 + d)
  expect_equal(closed_form, -31.39003523989096, tolerance = 1e-15)
  observed <- .marglik_cluster_norm_log_lik(
    yi = c(10, 10), mu = c(0, 0), tau_within = c(0, 0), tau_between = c(1, 1),
    sei = rep(2^-30, 2L), cluster = c(1L, 1L)
  )
  expect_equal(observed, closed_form, tolerance = 1e-14)

  # The R reference without the native kernel agrees.
  testthat::local_mocked_bindings(
    .has_native_norm_cluster_analytic = function(...) FALSE,
    .package = "RoBMA"
  )
  expect_equal(.marglik_cluster_norm_log_lik(
    yi = c(10, 10), mu = c(0, 0), tau_within = c(0, 0), tau_between = c(1, 1),
    sei = rep(2^-30, 2L), cluster = c(1L, 1L)
  ), closed_form, tolerance = 1e-14)
})


test_that("cluster bridge likelihood honors zero powers and row-specific loadings", {

  # One cluster of four rows with row-specific loadings (one of them zero) and
  # unequal powers; 80-digit mpmath value of the powered Gaussian integral.
  yi          <- c(.3, -.2, .5, .1)
  sei         <- c(.1, .2, .05, .3)
  tau_between <- c(.7, .35, 0, 1.4)
  weights     <- c(.5, 1, 1.5, 2)
  expect_equal(
    .marglik_cluster_norm_log_lik(yi, rep(0, 4L), rep(0, 4L), tau_between, sei,
      rep(1L, 4L), weights),
    -73.95729500926897, tolerance = 1e-13
  )

  # Rows with power 0 are omitted and a cluster without a positive power adds
  # nothing.
  extended <- .marglik_cluster_norm_log_lik(
    c(yi, 7, 8), rep(0, 6L), rep(0, 6L), c(tau_between, 3, 2), c(sei, 2, 1),
    c(rep(1L, 4L), 1L, 2L), c(weights, 0, 0)
  )
  expect_equal(extended, -73.95729500926897, tolerance = 1e-13)

  # A powered single-row cluster against numerical integration of the latent.
  integrated <- stats::integrate(function(g) vapply(g, function(value) {
    exp(2 * stats::dnorm(.1, mean = 1.4 * value, sd = .3, log = TRUE)) *
      stats::dnorm(value)
  }, numeric(1)), -Inf, Inf, rel.tol = 1e-12)[["value"]]
  expect_equal(
    .marglik_cluster_norm_log_lik(.1, 0, 0, 1.4, .3, 1L, 2),
    log(integrated), tolerance = 1e-11
  )
})


test_that("the bridge callback marginalizes the cluster effect without cancellation", {

  point  <- function(value) BayesTools::prior("point", list(value))
  object <- brma(yi = c(1, 1), sei = c(1e-8, 1e-8), cluster = c("a", "a"),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_effect = point(0), prior_heterogeneity = point(1),
    prior_heterogeneity_allocation = point(1))
  callback <- .log_posterior(
    parameters = list(mu = 0, tau = 1, rho = 1),
    data = .create_fit_data(object$data, object$priors),
    is_scale = FALSE, is_multilevel = TRUE, is_weights = FALSE,
    is_known_v = FALSE, is_PET = FALSE, is_PEESE = FALSE,
    is_weightfunction = FALSE, effect_direction = "positive",
    outcome_type = "norm", cluster_effects_marginalized = TRUE
  )
  d <- 1e-16
  expect_equal(callback,
    -log(2 * pi) - (log(d) + log(2 + d)) / 2 - 1 / (2 + d), tolerance = 1e-14)
})
