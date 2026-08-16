context("Prediction conditioning depth")

test_that("marginal prediction is independent of implicit versus explicit design", {

  dat <- data.frame(
    yi = c(0.2, 0.5),
    vi = c(0.04, 0.09)
  )
  object <- brma(
    yi                        = yi,
    vi                        = vi,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- matrix(
    c(
      0.10, 0.30,
      0.20, 0.40,
      0.15, 0.35
    ),
    nrow     = 3L,
    byrow    = TRUE,
    dimnames = list(NULL, c("mu", "tau"))
  )

  set.seed(804)
  implicit <- predict(
    object,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  set.seed(804)
  explicit <- predict(
    object,
    newdata            = dat,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_equal(implicit, explicit, tolerance = 0)
  expect_error(
    predict(
      object,
      newdata            = dat,
      type               = "estimate",
      conditioning_depth = "estimate",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "fitted observation identities"
  )
  expect_error(
    predict(
      object,
      type               = "terms",
      conditioning_depth = "estimate",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "only available for type = 'estimate' or type = 'response'",
    fixed = TRUE
  )
})


test_that("estimate depth includes conditional latent uncertainty", {

  S   <- 30000L
  mu  <- 0.2
  tau <- 0.4
  yi  <- 0.8
  sei <- 0.3
  tau2 <- tau^2
  vi   <- sei^2
  expected_mean <- mu + tau2 / (tau2 + vi) * (yi - mu)
  expected_var  <- tau2 * vi / (tau2 + vi)

  set.seed(805)
  observed <- .evaluate.brma.true_effects_posterior.norm(
    mu_samples = matrix(mu, nrow = S, ncol = 1L),
    tau_within = matrix(tau, nrow = S, ncol = 1L),
    yi         = yi,
    sei        = sei
  )

  expect_equal(mean(observed), expected_mean, tolerance = 0.004)
  expect_equal(as.numeric(stats::var(observed)), expected_var,
               tolerance = 0.002)
})


test_that("known-V estimate depth preserves the joint conditional posterior", {

  S   <- 30000L
  mu  <- c(0.1, -0.2)
  tau <- c(0.4, 0.3)
  yi  <- c(0.5, 0.1)
  V   <- matrix(c(0.09, 0.03, 0.03, 0.16), nrow = 2L)
  D   <- diag(tau^2)
  conditional_mean <- mu + D %*% solve(D + V, yi - mu)
  conditional_v    <- D - D %*% solve(D + V, D)
  known_V <- .known_v_newdata_prepare(V, k = 2L)

  set.seed(806)
  observed <- .evaluate.brma.known_v_posterior.norm(
    mu_samples = matrix(mu, nrow = S, ncol = 2L, byrow = TRUE),
    tau_within = matrix(tau, nrow = S, ncol = 2L, byrow = TRUE),
    yi         = yi,
    known_V    = known_V
  )

  expect_equal(colMeans(observed), as.vector(conditional_mean), tolerance = 0.004)
  expect_equal(stats::cov(observed), conditional_v, tolerance = 0.003)
})


test_that("legacy multilevel estimate depth preserves joint latent uncertainty", {

  S       <- 15000L
  mu      <- c(0.1, -0.2)
  tau_b   <- 0.3
  tau_w   <- c(0.4, 0.25)
  yi      <- c(0.5, 0.1)
  vi      <- c(0.09, 0.16)
  latent_v <- tau_b^2 * matrix(1, nrow = 2L, ncol = 2L) +
    diag(tau_w^2)
  marginal_v     <- latent_v + diag(vi)
  conditional_mean <- mu + latent_v %*% solve(marginal_v, yi - mu)
  conditional_v    <- latent_v -
    latent_v %*% solve(marginal_v, latent_v)

  set.seed(807)
  observed <- matrix(mu, nrow = S, ncol = 2L, byrow = TRUE) +
    .evaluate.brma.multilevel_posterior.norm(
      mu_samples = matrix(mu, nrow = S, ncol = 2L, byrow = TRUE),
      tau_within = matrix(tau_w, nrow = S, ncol = 2L, byrow = TRUE),
      tau_between = matrix(tau_b, nrow = S, ncol = 2L),
      yi          = yi,
      vi          = vi,
      cluster     = c(1L, 1L)
    )

  expect_equal(colMeans(observed), as.vector(conditional_mean), tolerance = 0.006)
  expect_equal(stats::cov(observed), conditional_v, tolerance = 0.005)
})
