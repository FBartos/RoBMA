probe_state("test-00-selection-dense-triple.R")
test_that("the rejected Assink triple is recovered within its reported radius", {

  # The first rejected JAGS state from issue 115. Preserve the full supplied
  # covariance; the independent oracle integrates Gaussian orthant probabilities.
  se <- c(.3920459156782532, .48311489316724648, .40914545090957566)
  V <- matrix(c(.1979656538431033, .15152257653564372, .10265843388635931,
                .15152257653564372, .27766565384310327, .15813140864483566,
                .10265843388635931, .15813140864483566, .21166565384310329), 3L)
  mu <- rep(-.18060343668680501, 3L)
  omega <- c(1, .20760940515684104)
  cut <- stats::qnorm(.025, lower.tail = FALSE)
  rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  evaluate <- function(order = 1:3, quadrature = rules, tolerance = .005,
                       points = 512L, scrambles = 8L, covariance = V,
                       mu_vec = mu, sei = se) {
    covariance <- covariance[order, order]
    packed <- covariance[lower.tri(covariance, diag = TRUE)]
    design <- BayesTools::selection_qmc_design(6L, points, scrambles, seed = 1L)
    .Call("RoBMA_selnorm_mnorm_step_loglik_batch", numeric(3), matrix(mu_vec[order], 2, 3, byrow = TRUE),
          matrix(packed, 2, 6, byrow = TRUE), sei[order], matrix(omega, 2, 2, byrow = TRUE),
          c(cut, -Inf), c(Inf, cut), rep(2L, 3), 1L, TRUE, 1L, as.double(design),
          points, scrambles, tolerance, TRUE, 0L, quadrature, PACKAGE = "RoBMA")
  }

  # w(y) = omega[2] + (1 - omega[2]) I(y >= cut*se) for each row.
  reference_mass <- function(V, mu_vec = mu, sei = se) {
    mass <- omega[2]^3
    for (size in 1:3) {
      orthants <- utils::combn(3, size, function(rows) {
        if (size == 1L) return(stats::pnorm(cut*sei[rows], mu_vec[rows], sqrt(V[rows, rows]),
                                            lower.tail = FALSE))
        as.numeric(mvtnorm::pmvnorm(lower = cut*sei[rows], mean = mu_vec[rows],
          sigma = V[rows, rows, drop = FALSE], algorithm = mvtnorm::Miwa(steps = 512L)))
      })
      mass <- mass + omega[2]^(3-size) * (1-omega[2])^size * sum(orthants)
    }
    mass
  }
  mass <- reference_mass(V)
  qmc <- evaluate(quadrature = NULL)
  expect_gt(.selection_qmc_z() * qmc$relative_mcse[1], .005)
  recovered <- evaluate()
  radius <- attr(recovered$integration_diagnostics, "rank_one_relative_error")
  expect_length(radius, 2)
  expect_identical(radius[1], radius[2])
  expect_equal(recovered$relative_mcse, c(0, 0))
  expect_lte(radius[1], .005)
  expect_lte(abs(expm1(recovered$log_normalizer[1] - log(mass))), radius[1] + 1e-8)
  expect_equal(recovered$log_density[1], mvtnorm::dmvnorm(numeric(3), mu, V, log = TRUE) +
    3*log(omega[2]) - recovered$log_normalizer[1], tolerance = 1e-12)
  permuted <- evaluate(c(2, 1, 3))
  expect_lte(abs(expm1(permuted$log_normalizer[1] - log(mass))), .005)

  event <- function(upper, quadrature) {
    design <- BayesTools::selection_qmc_design(6L, 512L, 8L, seed = 1L)
    .Call("RoBMA_selnorm_gaussian_event_mass_batch", matrix(mu, 1), matrix(V[lower.tri(V, diag = TRUE)], 1),
      se, matrix(omega, 1L), c(cut, -Inf), c(Inf, cut), 1L, 1L, 0L, matrix(-Inf, 1, 3),
      matrix(upper, 1), as.double(design), 512L, 8L, NULL, quadrature, .005, NULL, PACKAGE = "RoBMA")
  }
  full <- event(rep(Inf, 3), rules)
  expect_equal(full$log_mass, recovered$log_normalizer[1], tolerance = 1e-12)
  expect_equal(full$relative_quadrature_error, radius[1], tolerance = 1e-12)
  expect_equal(full$relative_mcse, 0)
  expect_equal(event(c(0, Inf, Inf), rules), event(c(0, Inf, Inf), NULL))

  negative <- V * tcrossprod(c(1, -1, 1))
  unchanged <- evaluate(covariance = negative)
  baseline <- evaluate(covariance = negative, quadrature = NULL)
  expect_identical(unchanged$relative_mcse, baseline$relative_mcse)
  expect_identical(unchanged$log_normalizer, baseline$log_normalizer)

  # The fixed-effect AR(1) model is positive definite; its auxiliary factor
  # has zero uniqueness on the middle row. Check the supplied V against an
  # independent integral, also just inside that boundary.
  sampling <- tcrossprod(se) * .8^abs(outer(1:3, 1:3, "-"))
  for (tau2 in c(0, 1e-12)) {
    cov <- sampling + diag(tau2, 3)
    result <- evaluate(covariance = cov)
    bound <- attr(result$integration_diagnostics, "rank_one_relative_error")
    expect_length(bound, 2L)
    expect_equal(result$relative_mcse, c(0, 0))
    expect_lte(bound[1L], .005)
    expect_lte(abs(expm1(result$log_normalizer[1] - log(reference_mass(cov)))), bound[1] + 1e-8)
  }
  heywood <- matrix(c(1,.8,.4, .8,1,.8, .4,.8,1), 3L)
  refused <- evaluate(covariance = heywood)
  expect_null(attr(refused$integration_diagnostics, "rank_one_relative_error"))
  expect_gt(.selection_qmc_z() * refused$relative_mcse[1L], .005)

  # An unequal-SE AR(.95) block needs the full covariance eigenvalue bound.
  alt_se <- c(.15, .3, .6)
  alt_V <- tcrossprod(alt_se) * .95^abs(outer(1:3, 1:3, "-"))
  alt <- evaluate(covariance = alt_V, mu_vec = -alt_se, sei = alt_se)
  alt_radius <- attr(alt$integration_diagnostics, "rank_one_relative_error")
  expect_length(alt_radius, 2L)
  expect_lte(alt_radius[1L], .005)
  expect_lte(abs(expm1(alt$log_normalizer[1L] - log(reference_mass(alt_V, -alt_se, alt_se)))),
             alt_radius[1L] + 1e-8)

  # Deliberately insufficient work must never count as a bounded success.
  limited <- evaluate(points = 8L, scrambles = 2L, tolerance = 1e-9)
  expect_null(attr(limited$integration_diagnostics, "rank_one_relative_error"))
  expect_gt(.selection_qmc_z() * limited$relative_mcse[1], 1e-9)
})
