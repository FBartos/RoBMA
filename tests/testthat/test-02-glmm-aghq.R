context("Adaptive GLMM quadrature")


# Independently integrate the binomial GLMM over theta and the base rate.
.binom_aghq_oracle <- function(a, c, n1, n2, mu, tau, alpha, beta,
                               weight = 1) {

  log_coefficient <- lchoose(n1, a) + lchoose(n2, c)
  integrate_theta <- function(theta) {

    vapply(theta, function(theta_value) {
      effect <- mu + tau * theta_value
      inner <- stats::integrate(
        function(pi0) {
          p1 <- stats::plogis(stats::qlogis(pi0) + 0.5 * effect)
          p2 <- stats::plogis(stats::qlogis(pi0) - 0.5 * effect)
          exp(weight * (
            log_coefficient +
              a * log(p1) + (n1 - a) * log1p(-p1) +
              c * log(p2) + (n2 - c) * log1p(-p2)
          )) * stats::dbeta(pi0, alpha, beta)
        },
        lower         = 0,
        upper         = 1,
        subdivisions  = 1000L,
        rel.tol       = 1e-10,
        stop.on.error = TRUE
      )[["value"]]

      inner * stats::dnorm(theta_value)
    }, numeric(1))
  }

  value <- stats::integrate(
    integrate_theta,
    lower         = -9,
    upper         = 9,
    subdivisions  = 1000L,
    rel.tol       = 1e-9,
    stop.on.error = TRUE
  )[["value"]]

  return(log(value))
}


# Independently integrate the Poisson GLMM over theta and the log base rate.
.pois_aghq_oracle <- function(x1, x2, t1, t2, mu, tau, mean, sd,
                              weight = 1) {

  integrate_theta <- function(theta) {

    vapply(theta, function(theta_value) {
      effect <- mu + tau * theta_value
      inner <- stats::integrate(
        function(phi) {
          lambda1 <- t1 * exp(phi + 0.5 * effect)
          lambda2 <- t2 * exp(phi - 0.5 * effect)
          exp(weight * (
            stats::dpois(x1, lambda1, log = TRUE) +
              stats::dpois(x2, lambda2, log = TRUE)
          )) * stats::dnorm(phi, mean, sd)
        },
        lower         = mean - 10 * sd,
        upper         = mean + 10 * sd,
        subdivisions  = 1000L,
        rel.tol       = 1e-10,
        stop.on.error = TRUE
      )[["value"]]

      inner * stats::dnorm(theta_value)
    }, numeric(1))
  }

  value <- stats::integrate(
    integrate_theta,
    lower         = -9,
    upper         = 9,
    subdivisions  = 1000L,
    rel.tol       = 1e-9,
    stop.on.error = TRUE
  )[["value"]]

  return(log(value))
}


# Independently refine the broad-prior Poisson integral when tau is zero.
.pois_zero_tau_oracle <- function(x1, x2, t1, t2, mu, mean, sd) {

  log_kernel <- function(phi) {
    stats::dpois(x1, t1 * exp(phi + 0.5 * mu), log = TRUE) +
      stats::dpois(x2, t2 * exp(phi - 0.5 * mu), log = TRUE) +
      stats::dnorm(phi, mean, sd, log = TRUE)
  }
  mode <- stats::optimize(
    log_kernel,
    interval = c(-30, 5),
    maximum  = TRUE,
    tol      = 1e-13
  )
  scaled <- stats::integrate(
    function(phi) exp(log_kernel(phi) - mode[["objective"]]),
    lower         = -Inf,
    upper         = Inf,
    subdivisions  = 2000L,
    rel.tol       = 1e-12,
    stop.on.error = TRUE
  )[["value"]]

  return(mode[["objective"]] + log(scaled))
}


test_that("binomial AGHQ matches closed form and streams row sums", {

  ai        <- c(505L, 5L)
  ci        <- c(499L, 3L)
  n1i       <- c(88391L, 2498L)
  n2i       <- c(88391L, 2341L)
  mu        <- matrix(c(0, 0, 0, 0), nrow = 2)
  tau       <- matrix(0, nrow = 2, ncol = 2)
  weights   <- c(1, 0.5)
  alpha     <- 1.5
  beta      <- 2.5
  prior     <- c(alpha = alpha, beta = beta)

  out     <- .glmm_binom_aghq(
    ai, ci, n1i, n2i, mu, tau, weights, prior
  )
  out_sum <- .glmm_binom_aghq(
    ai, ci, n1i, n2i, mu, tau, weights, prior,
    row_sum = TRUE
  )

  expected <- weights * (lchoose(n1i, ai) + lchoose(n2i, ci)) +
    lbeta(alpha + weights * (ai + ci),
          beta + weights * (n1i + n2i - ai - ci)) -
    lbeta(alpha, beta)

  expect_equal(out[["value"]][1, ], expected, tolerance = 1e-7,
               info = "closed-form powered beta-binomial ordinate")
  expect_equal(out[["value"]][2, ], expected, tolerance = 1e-7,
               info = "posterior-row invariance at zero heterogeneity")
  expect_equal(out_sum[["value"]], rowSums(out[["value"]]),
               tolerance = 1e-10,
               info = "streaming binomial row sums")
  expect_equal(out[["max_change"]], 0,
               info = "exact beta-binomial shortcut")
  expect_equal(out[["max_order"]], 0L,
               info = "exact rows do not enter adaptive quadrature")
  expect_equal(sum(out[["order_counts"]]), 0L,
               info = "exact rows do not increment quadrature counts")
  expect_equal(out[["exact_count"]], length(mu),
               info = "one exact diagnostic count per cell")
})


test_that("binomial AGHQ initializes extreme beta priors on the log-odds scale", {

  events <- .glmm_binom_aghq(
    1L, 1L, 2L, 2L,
    matrix(.1), matrix(.2), NULL,
    c(alpha = 1e16, beta = 1)
  )
  non_events <- .glmm_binom_aghq(
    1L, 1L, 2L, 2L,
    matrix(-.1), matrix(.2), NULL,
    c(alpha = 1, beta = 1e16)
  )

  expect_true(is.finite(events[["value"]][1L, 1L]))
  expect_equal(events[["value"]], non_events[["value"]], tolerance = 1e-10)
  expect_identical(events[["max_mode_iterations"]], 1L)
  expect_identical(non_events[["max_mode_iterations"]], 1L)
})


test_that("exact beta-binomial shortcut rejects nonfinite powered values", {

  expect_error(
    .glmm_binom_aghq(
      1L, 1L, 2L, 2L,
      matrix(0), matrix(0), 1e308,
      c(alpha = 1, beta = 1)
    ),
    "Exact powered beta-binomial value is nonfinite at sample 1, observation 1"
  )
})


test_that("AGHQ rejects matrices beyond its integer diagnostic capacity", {

  skip_if(.Machine[["sizeof.pointer"]] < 8L,
          "long-vector ALTREP requires a 64-bit R build")

  oversized <- seq_len(65536 * 32768)
  dim(oversized) <- c(65536L, 32768L)

  expect_error(
    .glmm_binom_aghq(
      0L, 0L, 1L, 1L,
      oversized, oversized, NULL,
      c(alpha = 1, beta = 1)
    ),
    "AGHQ supports at most 2147483647 matrix cells",
    info = "cell cap is checked before compact ALTREP matrices materialize"
  )
  expect_error(
    .glmm_pois_aghq(
      0L, 0L, 1, 1,
      oversized, oversized, NULL,
      c(mean = 0, sd = 1)
    ),
    "AGHQ supports at most 2147483647 matrix cells",
    info = "Poisson dispatch shares the overflow-safe extent validation"
  )
})


test_that("AGHQ certifies difficult rows from the cached GLMM fits", {

  skip_on_cran()

  ai  <- c(505L, 186L)
  ci  <- c(499L, 141L)
  n1i <- c(88391L, 50634L)
  n2i <- c(88391L, 27338L)
  bin <- .glmm_binom_aghq(
    ai, ci, n1i, n2i,
    matrix(0, nrow = 1, ncol = 2),
    matrix(0, nrow = 1, ncol = 2),
    NULL,
    c(alpha = 1, beta = 1)
  )
  expected_bin <- lchoose(n1i, ai) + lchoose(n2i, ci) +
    lbeta(ai + ci + 1, n1i + n2i - ai - ci + 1)

  expect_equal(bin[["value"]][1, ], expected_bin, tolerance = 1e-7,
               info = "BCG rows 8 and 11 beta-binomial ordinates")
  expect_equal(
    unname(bin[["value"]][1, ]),
    c(-15.779757485139271, -18.919776158465083),
    tolerance = 1e-7,
    info = "fixed BCG regression references"
  )

  mean <- -6.01874892013504
  sd   <- 20.2747133148659
  pois <- .glmm_pois_aghq(
    c(7L, 17L), c(11L, 15L),
    c(1344, 6760), c(1988, 6840),
    matrix(-0.45, nrow = 1, ncol = 2),
    matrix(0, nrow = 1, ncol = 2),
    NULL,
    c(mean = mean, sd = sd)
  )
  expected_pois <- c(
    .pois_zero_tau_oracle(7L, 11L, 1344, 1988, -0.45, mean, sd),
    .pois_zero_tau_oracle(17L, 15L, 6760, 6840, -0.45, mean, sd)
  )

  expect_equal(pois[["value"]][1, ], expected_pois, tolerance = 1e-6,
               info = "Nielweise rows 1 and 7 refined broad-prior ordinates")
  expect_equal(pois[["exact_count"]], 0L)
  expect_equal(sum(pois[["order_counts"]]), length(pois[["value"]]))
  expect_equal(
    unname(pois[["value"]][1, ]),
    c(-8.7933665774352523, -10.728389896154845),
    tolerance = 1e-6,
    info = "fixed Nielweise regression references"
  )
})


test_that("joint binomial AGHQ matches independent integration", {

  skip_on_cran()

  expected <- .binom_aghq_oracle(
    a      = 7L,
    c      = 3L,
    n1     = 24L,
    n2     = 21L,
    mu     = 0.35,
    tau    = 0.7,
    alpha  = 1.5,
    beta   = 2.5,
    weight = 0.8
  )
  out <- .glmm_binom_aghq(
    7L, 3L, 24L, 21L,
    matrix(0.35), matrix(0.7), 0.8,
    c(alpha = 1.5, beta = 2.5)
  )

  expect_equal(out[["value"]][1, 1], expected, tolerance = 2e-6,
               info = "independent joint binomial integral")
})


test_that("joint Poisson AGHQ matches independent integration and row sums", {

  skip_on_cran()

  x1i     <- c(7L, 1L)
  x2i     <- c(11L, 2L)
  t1i     <- c(15, 370)
  t2i     <- c(18, 483)
  mu      <- matrix(c(0.35, -0.2, 0.35, -0.2), nrow = 2)
  tau     <- matrix(c(0.6, 0.3, 0.6, 0.3), nrow = 2)
  weights <- c(0.8, 1)
  prior   <- c(mean = -1, sd = 1.3)

  out <- .glmm_pois_aghq(
    x1i, x2i, t1i, t2i, mu, tau, weights, prior
  )
  out_sum <- .glmm_pois_aghq(
    x1i, x2i, t1i, t2i, mu, tau, weights, prior,
    row_sum = TRUE
  )
  expected <- .pois_aghq_oracle(
    x1 = 7L, x2 = 11L, t1 = 15, t2 = 18,
    mu = 0.35, tau = 0.6, mean = -1, sd = 1.3,
    weight = 0.8
  )

  expect_equal(out[["value"]][1, 1], expected, tolerance = 2e-6,
               info = "independent joint Poisson integral")
  expect_equal(out_sum[["value"]], rowSums(out[["value"]]),
               tolerance = 1e-10,
               info = "streaming Poisson row sums")
  expect_true(out[["max_change"]] <= 1e-6,
              info = "reported Poisson convergence")
  expect_equal(sum(out[["order_counts"]]), length(mu),
               info = "one Poisson diagnostic count per cell")
})


test_that("ordinary wrappers require AGHQ-supported default priors", {

  prior_pi  <- BayesTools::prior("beta", list(1.5, 2.5))
  prior_phi <- BayesTools::prior("normal", list(-1, 1.3))

  bin_raw  <- .glmm_binom_aghq(
    7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), NULL,
    c(alpha = 1.5, beta = 2.5)
  )
  pois_raw <- .glmm_pois_aghq(
    7L, 11L, 15, 18, matrix(0.35), matrix(0.6), NULL,
    c(mean = -1, sd = 1.3)
  )
  expected_bin  <- .glmm_aghq_value(bin_raw)
  expected_pois <- .glmm_aghq_value(pois_raw)

  expect_equal(
    as.numeric(.outcome_pdf.binom(
      7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), prior_pi
    )),
    as.numeric(expected_bin),
    tolerance = 1e-12,
    info = "default binomial wrapper AGHQ dispatch"
  )
  expect_equal(
    as.numeric(.outcome_pdf.pois(
      7L, 11L, 15, 18, matrix(0.35), matrix(0.6), prior_phi
    )),
    as.numeric(expected_pois),
    tolerance = 1e-12,
    info = "default Poisson wrapper AGHQ dispatch"
  )
  bin_sum  <- .outcome_pdf_sum.binom(
    7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), prior_pi
  )
  pois_sum <- .outcome_pdf_sum.pois(
    7L, 11L, 15, 18, matrix(0.35), matrix(0.6), prior_phi
  )
  expect_equal(
    as.numeric(bin_sum),
    rowSums(bin_raw[["value"]]),
    tolerance = 1e-12,
    info = "default binomial streaming AGHQ dispatch"
  )
  expect_equal(
    as.numeric(pois_sum),
    rowSums(pois_raw[["value"]]),
    tolerance = 1e-12,
    info = "default Poisson streaming AGHQ dispatch"
  )
  expect_named(
    attr(bin_sum, "glmm_aghq_diagnostics"),
    c(
      "max_order", "max_change", "max_mode_iterations", "order_counts",
      "exact_count", "grid_columns", "grid_max_theta_order",
      "grid_max_nuisance_order", "grid_max_change"
    ),
    info = "binomial wrapper diagnostics"
  )
  expect_named(
    attr(pois_sum, "glmm_aghq_diagnostics"),
    c(
      "max_order", "max_change", "max_mode_iterations", "order_counts",
      "exact_count", "grid_columns", "grid_max_theta_order",
      "grid_max_nuisance_order", "grid_max_change"
    ),
    info = "Poisson wrapper diagnostics"
  )

  legacy <- .outcome_pdf.binom(
    7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), prior_pi,
    n_theta = 7L, n_pi = 9L
  )
  expect_true(is.matrix(legacy) && is.finite(legacy[1, 1]),
              info = "explicit quadrature arguments preserve legacy path")

  truncated_phi <- BayesTools::prior(
    "normal",
    list(-1, 1.3),
    truncation = list(-3, 2)
  )
  truncated_pi <- BayesTools::prior(
    "beta",
    list(1.5, 2.5),
    truncation = list(.01, .99)
  )
  expect_null(.glmm_aghq_prior_spec(truncated_phi, "pois"),
              info = "truncated log-rate prior is outside certified AGHQ")
  expect_null(.glmm_aghq_prior_spec(truncated_pi, "bin"),
              info = "truncated baserate prior is outside certified AGHQ")
  expect_error(
    .outcome_pdf.pois(
      7L, 11L, 15, 18, matrix(0.35), matrix(0.6), truncated_phi
    ),
    "requires an untruncated normal 'prior_phi'",
    info = "unsupported default Poisson quadrature fails explicitly"
  )
  expect_error(
    .outcome_pdf.binom(
      7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), truncated_pi
    ),
    "requires an untruncated beta 'prior_pi'",
    info = "unsupported default binomial quadrature fails explicitly"
  )
  expect_true(is.finite(.outcome_pdf.pois(
    7L, 11L, 15, 18, matrix(0.35), matrix(0.6), truncated_phi,
    n_theta = 7L, n_phi = 9L
  )[1, 1]), info = "explicit Poisson legacy grid remains opt-in")
  expect_true(is.finite(.outcome_pdf.binom(
    7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), truncated_pi,
    n_theta = 7L, n_pi = 9L
  )[1, 1]), info = "explicit binomial legacy grid remains opt-in")
})


test_that("ordinary AGHQ recovers only failed observations", {

  prior <- BayesTools::prior("beta", list(1.5, 2.5))
  control <- .glmm_aghq_control(
    orders      = c(5L, 9L, 13L),
    tolerance   = 1e-15,
    consecutive = 2L
  )
  args <- list(
    ai          = c(0L, 7L),
    ci          = c(0L, 3L),
    n1i         = c(1L, 24L),
    n2i         = c(1L, 21L),
    mu_samples  = matrix(c(0, 0.35), nrow = 1),
    tau_within  = matrix(c(0, 0.7), nrow = 1),
    weights     = NULL,
    prior_pi    = prior,
    prior_spec  = .glmm_aghq_prior_spec(prior, "bin"),
    control     = control
  )

  pointwise <- do.call(.glmm_binom_marginal, args)
  summed <- do.call(
    .glmm_binom_marginal,
    c(args, list(row_sum = TRUE))
  )
  expected_second <- .binom_aghq_oracle(
    a = 7L, c = 3L, n1 = 24L, n2 = 21L,
    mu = 0.35, tau = 0.7, alpha = 1.5, beta = 2.5
  )
  diagnostics <- attr(pointwise, "glmm_aghq_diagnostics")

  expect_equal(pointwise[1L, 2L], expected_second, tolerance = 2e-6,
               info = "failed observation uses independently checked grid")
  expect_identical(diagnostics[["grid_columns"]], 2L,
                   info = "the exact first observation does not fall back")
  expect_equal(as.numeric(summed), rowSums(pointwise), tolerance = 1e-12,
               info = "pointwise and sum routes share recovered values")
})


test_that("AGHQ recovery does not catch unrelated errors", {

  expect_error(
    .glmm_aghq_dispatch(
      S             = 1L,
      K             = 1L,
      outcome_type  = "pois",
      evaluate_aghq = function(k) stop("unrelated programming error"),
      evaluate_grid = function(k) stop("fallback must not run"),
      row_sum       = FALSE
    ),
    "unrelated programming error"
  )
})


test_that("Poisson AGHQ handles extreme predictors without split exponentiation", {

  skip_on_cran()

  expected <- log(stats::integrate(
    function(z) {
      exp(
        stats::dpois(1L, exp(z), log = TRUE) +
          stats::dpois(0L, exp(z - 1500), log = TRUE)
      ) * stats::dnorm(z)
    },
    lower         = -10,
    upper         = 10,
    subdivisions  = 1000L,
    rel.tol       = 1e-12,
    stop.on.error = TRUE
  )[["value"]])
  out <- .glmm_pois_aghq(
    1L, 0L, 1, 1,
    matrix(1500), matrix(0), NULL,
    c(mean = -750, sd = 1)
  )

  expect_true(is.finite(out[["value"]][1, 1]),
              info = "finite result when separate exponentials would overflow")
  expect_equal(out[["value"]][1, 1], expected, tolerance = 1e-6,
               info = "stable extreme-predictor reference integral")
})


test_that("fixed Poisson grids preserve finite combined log rates", {

  prior <- BayesTools::prior("point", list(-750))
  args <- list(
    x1i        = 1L,
    x2i        = 0L,
    t1i        = 1,
    t2i        = 1,
    mu_samples = matrix(1500),
    tau_within = matrix(0),
    prior_phi  = prior,
    n_theta    = 15L,
    n_phi      = 1L
  )

  expect_equal(do.call(.outcome_pdf.pois, args)[1L, 1L], -1,
               tolerance = 1e-12, info = "native combined-log-rate kernel")
  expect_equal(do.call(.outcome_pdf.pois_r, args)[1L, 1L], -1,
               tolerance = 1e-12, info = "R combined-log-rate reference")
})


test_that("unstable high-order Hermite rules fail structurally", {

  rules <- .glmm_aghq_control()[["rules"]][["nodes"]]
  expect_identical(length(rules[[length(rules)]]), 49L)
  expect_error(
    .gauss_hermite_nodes(65L),
    "Gauss-Hermite weights are numerically unstable at order 65"
  )
})


test_that("AGHQ reports failure when refinement is capped", {

  control <- .glmm_aghq_control(
    orders      = c(5L, 9L, 13L),
    tolerance   = 1e-15,
    consecutive = 2L
  )

  expect_error(
    .glmm_pois_aghq(
      1L, 1L, 370, 483,
      matrix(-0.45), matrix(0.5), NULL,
      c(mean = -6.0187, sd = 20.27),
      control = control
    ),
    regexp = "failed to converge.*sample 1, observation 1",
    info   = "hard failure includes cell and terminal refinement"
  )
})
