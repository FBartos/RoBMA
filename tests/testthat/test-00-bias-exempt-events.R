context("Publication-bias exemptions in Gaussian selection events")

.bias_exempt_half_space <- function(rule, draws = 1L) {

  list(omega = matrix(rep(c(1, 0), each = draws), draws), sign = 1L,
       vector_rule = rep.int(rule, draws), kernel_mode = rep.int(1L, draws),
       use_normal = rep(FALSE, draws), p_cuts = c(0, .5, 1),
       z_lower = c(0, -Inf), z_upper = c(Inf, 0),
       bias_exempt = c(FALSE, TRUE))
}

test_that("mixed partial events preserve the exempt coordinate's covariance", {

  rho <- .6
  gaussian <- matrix(c(1, rho, 1), 1L)
  mu <- matrix(c(0, 0), 1L)
  se <- c(1, 1)
  plan <- list(points_per_scramble = 4096L, scrambles = 8L, seed = 551L)
  for (rule in 0:1) {
    selection <- .bias_exempt_half_space(rule)
    numerator <- function(rows, y) .selection_joint_event_numerator(
      mu, gaussian, se, selection, plan, rows, y
    )[["log_numerator"]]
    # Given Y2 = y, Y1 is normal with mean rho*y and variance 1-rho^2.
    expect_equal(numerator(2L, -1), dnorm(-1, log = TRUE) +
      pnorm(-rho / sqrt(1 - rho^2), log.p = TRUE), tolerance = 1e-12)
    expect_equal(numerator(1L, 1), dnorm(1, log = TRUE), tolerance = 1e-12)
    expect_identical(numerator(1L, -1), -Inf)

    mass <- .selection_gaussian_event_mass(mu, gaussian, se, selection, plan,
      upper = c(Inf, 0))
    # P(Y1 > 0, Y2 < 0) for a standard bivariate Gaussian.
    expect_equal(exp(mass$log_mass), .25 - asin(rho) / (2 * pi), tolerance = 2e-4)
    line <- .selection_gaussian_event_mass(mu, selection_se = se,
      selection_context = selection, execution_plan = plan,
      rank_one_loading = matrix(c(1, -1), 1L), upper = c(Inf, 0))
    expect_equal(exp(line$log_mass), .5, tolerance = 1e-14)
  }
})

test_that("selected predictions preserve the Gaussian conditional of exempt rows", {

  S <- 250L
  withr::local_seed(4401)
  for (rule in 0:1) {
    # A known rank-one Gaussian has Y2 = -Y1 exactly. Exempting Y2 still
    # requires this identity even though its sign is opposite to Y1's.
    covariance <- .block_covariance_from_matrix(matrix(c(1, -1, -1, 1), 2L),
      S, list(1:2))
    draws <- .outcome_rng.selnorm_mvn(matrix(0, S, 2L), covariance, c(1, 1),
      .bias_exempt_half_space(rule, S), list(1:2))
    expect_true(all(draws[, 1L] > 0))
    expect_equal(draws[, 2L], -draws[, 1L], tolerance = 1e-14)
  }
})
