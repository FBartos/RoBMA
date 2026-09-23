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
