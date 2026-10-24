source(testthat::test_path("common-functions.R"))

context("Adaptive GLMM quadrature")


# Independently integrate the binomial GLMM over theta and the base rate.
.binom_aghq_oracle <- function(a, c, n1, n2, mu, tau, alpha, beta,
                               weight = 1, lower = 0, upper = 1) {

  log_coefficient <- lchoose(n1, a) + lchoose(n2, c)
  log_mass <- if (lower == 0) {
    stats::pbeta(upper, alpha, beta, log.p = TRUE)
  } else if (upper == 1) {
    stats::pbeta(lower, alpha, beta, lower.tail = FALSE, log.p = TRUE)
  } else {
    use_upper_tail <- stats::pbeta(lower, alpha, beta) > .5
    large <- stats::pbeta(if (use_upper_tail) lower else upper,
      alpha, beta, lower.tail = !use_upper_tail, log.p = TRUE)
    small <- stats::pbeta(if (use_upper_tail) upper else lower,
      alpha, beta, lower.tail = !use_upper_tail, log.p = TRUE)
    large + log(-expm1(small - large))
  }
  inner_results <- list()
  count_log <- function(count, log_probability) {
    if (count == 0L) rep(0, length(log_probability)) else count * log_probability
  }
  log_integrate_theta <- function(theta) {

    vapply(theta, function(theta_value) {
      effect <- mu + tau * theta_value
      inner <- stats::integrate(
        function(pi0) {
          eta1 <- log(pi0) - log1p(-pi0) + 0.5 * effect
          eta2 <- log(pi0) - log1p(-pi0) - 0.5 * effect
          log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
          exp(weight * (
            log_coefficient +
              count_log(a, -log1pexp(-eta1)) + count_log(n1 - a, -log1pexp(eta1)) +
              count_log(c, -log1pexp(-eta2)) + count_log(n2 - c, -log1pexp(eta2))
          ) + stats::dbeta(pi0, alpha, beta, log = TRUE) - log_mass)
        },
        lower         = lower,
        upper         = upper,
        subdivisions  = 1000L,
        rel.tol       = 1e-10,
        stop.on.error = TRUE
      )
      inner_results[[length(inner_results) + 1L]] <<- list(theta = theta_value, integral = inner)

      log(inner[["value"]]) + stats::dnorm(theta_value, log = TRUE)
    }, numeric(1))
  }

  mode <- stats::optimize(
    log_integrate_theta,
    interval = c(-12, 12),
    maximum  = TRUE,
    tol      = 1e-12
  )
  step      <- 1e-4
  curvature <- -(
    log_integrate_theta(mode[["maximum"]] + step) -
      2 * mode[["objective"]] +
      log_integrate_theta(mode[["maximum"]] - step)
  ) / step^2
  value <- stats::integrate(
    function(z) {
      exp(log_integrate_theta(
        mode[["maximum"]] + z / sqrt(curvature)
      ) - mode[["objective"]]) / sqrt(curvature)
    },
    lower         = -Inf,
    upper         = Inf,
    subdivisions  = 1000L,
    rel.tol       = 1e-9,
    stop.on.error = TRUE
  )

  # QUADPACK errors and messages are retained; these are estimates, not a
  # rounding enclosure or a total nested-integration error certificate.
  return(list(value = mode[["objective"]] + log(value[["value"]]),
    inner = inner_results, outer = value, mode = mode, log_prior_mass = log_mass))
}


# Independently integrate the binomial GLMM over the base rate and theta for a
# concentrated beta prior. The base rate is integrated in the coordinate of the
# smaller probability p (the mirror image of the prior when alpha > beta), so
# that no exp() of a logit enters: dbeta(p, ..., log = TRUE) is evaluated at a
# representable p and the likelihood at log(p) and log1p(-p). R's dbeta shares
# the node-independent offset of up to about half of eps^2 * (alpha + beta) at
# large shapes (see the noise bound in src/glmm-aghq.cc); at the shape sums of
# the tests that use this reference (up to 1e20, 2.5e-12) it is below the
# changes the results are compared with. The range is the prior mean +/-
# `width` prior SDs, the integrand is scaled by its value at the mean, and the
# theta integral (tau > 0) is nested. Valid while the prior SD is resolved by
# the doubles around the mean of p, that is for skewed priors (the symmetric
# prior at 0.5 is not).
.binom_concentrated_reference <- function(a, c, n1, n2, mu, tau, alpha, beta,
                                          weight = 1, width = 12, theta_range = 10) {

  swapped <- alpha > beta
  shapes  <- if (swapped) c(beta, alpha) else c(alpha, beta)
  total   <- sum(shapes)
  centre  <- shapes[[1L]] / total
  spread  <- sqrt(shapes[[1L]] * shapes[[2L]] / (total^2 * (total + 1)))
  lower   <- max(centre - width * spread, 0)
  upper   <- min(centre + width * spread, 1)
  logit   <- function(p) if (swapped) log1p(-p) - log(p) else log(p) - log1p(-p)
  log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  coefficient <- lchoose(n1, a) + lchoose(n2, c)
  log_lik <- function(x, theta) {
    effect <- mu + tau * theta
    eta1   <- x + .5 * effect
    eta2   <- x - .5 * effect
    weight * (coefficient - a * log1pexp(-eta1) - (n1 - a) * log1pexp(eta1) -
                c * log1pexp(-eta2) - (n2 - c) * log1pexp(eta2))
  }
  log_prior_centre <- stats::dbeta(centre, shapes[[1L]], shapes[[2L]], log = TRUE)
  integrate_checked <- function(f, lower, upper) {
    value <- stats::integrate(f, lower, upper, rel.tol = 1e-13, abs.tol = 0,
      subdivisions = 1000L, stop.on.error = FALSE)
    if (!identical(value[["message"]], "OK")) {
      stop("The reference integral reported: ", value[["message"]], call. = FALSE)
    }
    value[["value"]]
  }
  inner <- function(theta) {
    scale <- log_lik(logit(centre), theta)
    integral <- integrate_checked(function(p) {
      exp(stats::dbeta(p, shapes[[1L]], shapes[[2L]], log = TRUE) - log_prior_centre +
            log_lik(logit(p), theta) - scale)
    }, lower, upper)
    log(integral) + log_prior_centre + scale
  }
  if (tau == 0) {
    return(inner(0))
  }
  log_theta <- function(theta) {
    vapply(theta, function(t) inner(t) + stats::dnorm(t, log = TRUE), numeric(1L))
  }
  log_max <- max(log_theta(seq(-theta_range, theta_range, length.out = 401L)))
  log(integrate_checked(function(t) exp(log_theta(t) - log_max),
    -theta_range, theta_range)) + log_max
}


# Independently integrate the Poisson GLMM over theta and the log base rate.
# For an untruncated target, the default +/-10 SD nuisance interval is
# renormalized; omitted prior mass is 2*pnorm(-10), about 1.524e-23. Positive
# powered probabilities are at most one, so the absolute mass change is bounded
# by that tail mass (before converting to a log-integral comparison).
.pois_aghq_oracle <- function(x1, x2, t1, t2, mu, tau, mean, sd,
                              weight = 1, lower = mean - 10 * sd,
                              upper = mean + 10 * sd) {

  prior_mass <- stats::pnorm(upper, mean, sd) -
    stats::pnorm(lower, mean, sd)
  log_integrate_theta <- function(theta) {

    vapply(theta, function(theta_value) {
      effect <- mu + tau * theta_value
      inner <- stats::integrate(
        function(phi) {
          lambda1 <- t1 * exp(phi + 0.5 * effect)
          lambda2 <- t2 * exp(phi - 0.5 * effect)
          exp(weight * (
            stats::dpois(x1, lambda1, log = TRUE) +
              stats::dpois(x2, lambda2, log = TRUE)
          )) * stats::dnorm(phi, mean, sd) / prior_mass
        },
        lower         = lower,
        upper         = upper,
        subdivisions  = 1000L,
        rel.tol       = 1e-10,
        stop.on.error = TRUE
      )[["value"]]

      log(inner) + stats::dnorm(theta_value, log = TRUE)
    }, numeric(1))
  }

  mode <- stats::optimize(
    log_integrate_theta,
    interval = c(-12, 12),
    maximum  = TRUE,
    tol      = 1e-12
  )
  step      <- 1e-4
  curvature <- -(
    log_integrate_theta(mode[["maximum"]] + step) -
      2 * mode[["objective"]] +
      log_integrate_theta(mode[["maximum"]] - step)
  ) / step^2
  value <- stats::integrate(
    function(z) {
      exp(log_integrate_theta(
        mode[["maximum"]] + z / sqrt(curvature)
      ) - mode[["objective"]]) / sqrt(curvature)
    },
    lower         = -Inf,
    upper         = Inf,
    subdivisions  = 1000L,
    rel.tol       = 1e-9,
    stop.on.error = TRUE
  )[["value"]]

  return(mode[["objective"]] + log(value))
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
  # Independently checked beta-one moment bounds at mu=.1, tau=.2 enclose
  # [-71.60328143412962686785, -71.60328143412962566032]. The 1e-6 absolute
  # allowance includes quadrature/rounding rather than claiming those bounds
  # enclose the floating-point computation itself.
  expect_lt(abs(events[["value"]][1L, 1L] + 71.603281434129627), 1e-6)
  expect_identical(events[["exact_count"]], 0L)
  expect_identical(sum(events[["order_counts"]]), 1L)
  expect_identical(events[["max_mode_iterations"]], 1L)
  expect_identical(non_events[["max_mode_iterations"]], 1L)
})


test_that("binomial AGHQ is continuous across the closed-form beta density bound", {

  # Shapes with alpha + beta <= 1e3 use the closed-form log density and larger
  # ones the density at the smaller probability. Scaling the shape sum by
  # 1 -/+ 1e-9 around the bound moves the true log likelihood by far less than
  # the 1e-8 relative allowance, which covers the closed form's rounding error
  # (about 1e-13 absolute) and the quadrature of each side.
  ai  <- c(2L, 5L)
  ci  <- c(4L, 9L)
  n1i <- c(30L, 60L)
  n2i <- c(30L, 55L)
  mu  <- matrix(c(.2, -.1), 1L, 2L)
  tau <- matrix(.3, 1L, 2L)
  for (fraction in c(.5, .2, .93, 1e-4)) {
    value <- vapply(1e3 * c(1 - 1e-9, 1 + 1e-9), function(shape_sum) {
      .glmm_binom_aghq(ai, ci, n1i, n2i, mu, tau, NULL,
        c(alpha = fraction * shape_sum, beta = (1 - fraction) * shape_sum))[["value"]]
    }, numeric(2L))
    expect_true(all(is.finite(value)))
    expect_lt(max(abs(value[, 1L] - value[, 2L]) / abs(value[, 1L])), 1e-8,
      label = paste("relative change across the bound at prior mean", fraction))
  }
})


test_that("binomial AGHQ mode search converges for concentrated beta priors over many draws", {

  # The closed-form log density has rounding noise of about
  # eps * (alpha * |lp| + beta * |lq|); with shape sums of 1e4 and more it
  # reaches mode_tolerance (1e-12) and stalls the Newton steps of the mode
  # search at some draws (RoBMA_glmm_aghq_nonconvergence). Above the bound the
  # density is evaluated at the smaller probability. Reference values: the
  # 4.1.0.420 build, which used that route for every non-uniform prior, at
  # draws 1, 44, 140 and 300 of observations 1, 4 and 6.
  ai  <- c(2L, 5L, 12L, 20L, 0L, 3L)
  ci  <- c(4L, 9L, 10L, 25L, 1L, 0L)
  n1i <- c(30L, 60L, 40L, 80L, 15L, 22L)
  n2i <- c(30L, 55L, 41L, 75L, 14L, 25L)
  draws <- withr::with_seed(21, {
    list(mu = matrix(rnorm(300L * 6L, 0, .4), 300L, 6L),
         tau = matrix(abs(rnorm(300L * 6L, .2, .1)), 300L, 6L))
  }, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  reference <- list(
    "25000_25000" = list(
      sum = -43119.633805268,
      cells = matrix(c(
        -25.77374187554875, -25.32112309540952, -25.22699212109874, -26.25838901753625,
        -19.22587533341839, -19.27964685074012, -18.70743191483005, -18.97089339229578,
        -25.16711273003426, -24.91001705588358, -26.00400693510275, -25.61741944804021),
        4L)),
    "10000_20000" = list(
      sum = -19995.851651314,
      cells = matrix(c(
        -12.65969280828293, -12.2106049738906, -12.09034974667973, -13.11930364628969,
        -5.897572385329831, -5.922391967190909, -5.411420863735282, -5.657378313055217,
        -13.74614753827331, -13.51290323260307, -14.49199384278415, -14.16031300434924),
        4L)))
  for (shapes in list(c(25000, 25000), c(10000, 20000))) {
    key <- paste(shapes, collapse = "_")
    out <- .glmm_binom_aghq(ai, ci, n1i, n2i, draws$mu, draws$tau, NULL,
      c(alpha = shapes[[1L]], beta = shapes[[2L]]))
    expect_true(all(is.finite(out[["value"]])), info = key)
    expect_lt(abs(sum(out[["value"]]) / reference[[key]]$sum - 1), 1e-12, label = key)
    expect_lt(max(abs(out[["value"]][c(1L, 44L, 140L, 300L), c(1L, 4L, 6L)] /
                        reference[[key]]$cells - 1)), 1e-12, label = key)
  }
})


test_that("binomial AGHQ does not report a change below the rounding noise of concentrated beta priors", {

  # A zero-event study at mu = .1 and tau = 0 integrates to its likelihood at
  # the prior center, -log1p(exp(.05)) - log1p(exp(-.05)), for alpha = beta.
  # Above the closed-form bound the beta log density carries a rounding noise
  # of about eps * sqrt(alpha + beta) that the change between consecutive
  # orders does not show (at alpha = beta = 1e21 the error was 2.6e-6 with a
  # reported change of 2.8e-8). The reported change is never below that noise,
  # so shape sums beyond (tolerance / eps)^2, about 2e19, are refused instead
  # of accepted with a certificate that does not hold.
  limit <- -log1p(exp(.05)) - log1p(exp(-.05))
  aghq <- function(shape) {
    .glmm_binom_aghq(0L, 0L, 1L, 1L, matrix(.1), matrix(0), NULL,
      c(alpha = shape, beta = shape))
  }

  for (shape in c(1e20, 1e21, 1e22)) {
    expect_error(aghq(shape), class = "RoBMA_glmm_aghq_nonconvergence",
      info = paste("shapes", shape))
  }
  # The refusal reports the noise bound, which the message does not call a
  # change between orders.
  condition <- tryCatch(aghq(1e21), error = identity)
  expect_gte(condition[["change"]], .Machine[["double.eps"]] * sqrt(2e21))
  expect_gt(condition[["noise"]], 0)
  expect_identical(condition[["change"]], condition[["noise"]])
  expect_match(conditionMessage(condition), paste0("the rounding-noise bound of the prior density was ",
    format(condition[["change"]], digits = 8L), " at quadrature order"), fixed = TRUE)
  expect_false(grepl("last log-likelihood change", conditionMessage(condition), fixed = TRUE))

  accepted <- aghq(1e18)
  expect_gte(accepted[["max_change"]], .Machine[["double.eps"]] * sqrt(2e18))
  expect_lte(abs(accepted[["value"]][1L, 1L] - limit), accepted[["max_change"]])

  # The public route takes the accurate grid where the quadrature is refused.
  out <- .outcome_pdf.binom(0L, 0L, 1L, 1L, matrix(.1), matrix(0),
    BayesTools::prior("beta", list(alpha = 1e21, beta = 1e21)))
  diagnostics <- attr(out, "glmm_aghq_diagnostics", exact = TRUE)
  expect_identical(diagnostics$grid_columns, 1L)
  expect_lt(abs(out[1L, 1L] - limit), 1e-9)
})


test_that("binomial AGHQ accepts skewed concentrated beta priors within its reported change", {

  # The rounding noise of the beta log density scales with the curvature
  # sqrt(alpha * beta / (alpha + beta)) of the logit-scale prior, not with the
  # shape sum: skewed priors with shape sums far above 2e19 are computed
  # accurately and their reported change bounds the error. Base AGHQ accepted
  # these; a bound that grows with the shape sum refused them and the grid then
  # failed to find finite interior nodes.
  for (shapes in list(c(1e12, 1e20), c(1e3, 1e20), c(1e20, 1e3), c(1e10, 1e20))) {
    label <- paste("shapes", shapes[[1L]], shapes[[2L]])
    reference <- .binom_concentrated_reference(3L, 7L, 20L, 20L, .1, .2, shapes[[1L]], shapes[[2L]])
    wider     <- .binom_concentrated_reference(3L, 7L, 20L, 20L, .1, .2, shapes[[1L]], shapes[[2L]],
      width = 16, theta_range = 12)
    expect_lt(abs(reference - wider), 1e-10, label = label)

    out <- .outcome_pdf.binom(3L, 7L, 20L, 20L, matrix(.1), matrix(.2),
      BayesTools::prior("beta", list(alpha = shapes[[1L]], beta = shapes[[2L]])))
    diagnostics <- attr(out, "glmm_aghq_diagnostics", exact = TRUE)
    expect_length(diagnostics$grid_columns, 0L)
    expect_gte(diagnostics$max_change, 2 * .Machine[["double.eps"]] *
      sqrt(shapes[[1L]] * shapes[[2L]] / (shapes[[1L]] + shapes[[2L]])))
    expect_lte(abs(out[1L, 1L] - reference), diagnostics$max_change, label = label)
  }
})


test_that("binomial AGHQ bounds the node-independent dbeta error of very skewed concentrated priors", {

  # R's dbeta is off by a constant, measured up to about half of
  # eps^2 * (alpha + beta) at large shapes. It shifts the log likelihood of
  # every quadrature order alike, so the change between orders does not show it,
  # and for a skewed prior the slope term of the noise bound no longer contains
  # it. The references are literals of a 60-digit mpmath integral in
  # probability coordinates (the coordinate of the smaller probability, so no
  # exp() of a logit enters), cross-checked by a 300-bit Rmpfr integral; they
  # are not computed with stats::dbeta. Data 19/27 and 6/17, mu = .738, tau = 0.
  eps <- .Machine[["double.eps"]]
  aghq <- function(alpha, beta) {
    .glmm_binom_aghq(19L, 6L, 27L, 17L, matrix(.738), matrix(0), NULL,
      c(alpha = alpha, beta = beta))
  }

  # Accepted: the bound (4.9e-7 at the shape sum 1e25) is below the tolerance
  # and covers the error of 2.3e-7, which a bound without the offset (2.7e-10)
  # did not.
  accepted <- aghq(1e9, 1e25)
  reference <- -892.2004709683037946278289
  expect_gte(accepted[["max_change"]], eps^2 * (1e9 + 1e25))
  expect_lte(abs(accepted[["value"]][1L, 1L] - reference), accepted[["max_change"]])

  # Refused: the error is 2.4e-4 under a change of 8.4e-9 without the offset;
  # the bound (4.9e-4) refuses the prior.
  condition <- tryCatch(aghq(1e12, 1e28), error = identity)
  expect_s3_class(condition, "RoBMA_glmm_aghq_nonconvergence")
  expect_gte(condition[["change"]], eps^2 * (1e12 + 1e28))
  expect_gte(condition[["noise"]], 1e-6)
  expect_match(conditionMessage(condition), "the rounding-noise bound of the prior density was", fixed = TRUE)

  # With a tolerance above the bound the same prior is accepted and the change
  # covers the error against the reference value below.
  reference_refused <- -892.2004712680037888707065
  loose <- .glmm_binom_aghq(19L, 6L, 27L, 17L, matrix(.738), matrix(0), NULL,
    c(alpha = 1e12, beta = 1e28), control = .glmm_aghq_control(tolerance = 1e-3))
  expect_lte(abs(loose[["value"]][1L, 1L] - reference_refused), loose[["max_change"]])

  # The public route then goes to the prior-CDF grid, which cannot represent a
  # prior with shapes this extreme; the error names the prior and the remedies.
  condition <- tryCatch(.outcome_pdf.binom(19L, 6L, 27L, 17L, matrix(.738), matrix(0),
    BayesTools::prior("beta", list(alpha = 1e12, beta = 1e28))), error = identity)
  expect_s3_class(condition, "error")
  expect_identical(conditionMessage(condition), paste0(
    "The Beta(1e+12, 1e+28) base-rate prior in 'prior_baserate' is unavailable: ",
    "its shapes are too extreme for the quadrature nodes and weights of the ",
    "prior-CDF grid to be represented in double precision (a very concentrated ",
    "or an extremely diffuse prior). Use a beta prior with less extreme shapes ",
    "or a fixed base rate, e.g. 'prior_baserate = BayesTools::prior(\"point\", list(0.1))'."
  ))
})


test_that("exact beta-binomial shortcut keeps event increments for large shapes", {

  # One event and one non-event: the sequence probability is
  # a / (a + b) * b / (a + b + 1). The ratio of Beta functions through
  # log-gamma values of size 1e21 lost the increments (0 at the shape 1e20).
  for (shapes in list(c(1e20, 1e20), c(1e16, 2e16), c(1e14, 1e14),
                      c(1e12, 3e12), c(1e8, 1e8), c(30, 70))) {
    out <- .glmm_binom_aghq(
      1L, 0L, 1L, 1L, matrix(0), matrix(0), NULL,
      c(alpha = shapes[[1L]], beta = shapes[[2L]])
    )
    a <- shapes[[1L]]
    b <- shapes[[2L]]
    expect_equal(out[["value"]][1L, 1L], log(a / (a + b)) + log(b / (a + b + 1)),
      tolerance = 1e-12, info = paste("shapes", a, b))
    expect_equal(out[["exact_count"]], 1L)
  }
})


test_that("exact beta-binomial shortcut keeps shapes whose sum overflows", {

  # One event and one non-event again, with alpha + beta beyond the largest
  # double (the sum overflowed to Inf and the shortcut stopped with the
  # nonfinite-value error although the ratio is representable). 900-digit
  # mpmath values of log(a / (a + b)) + log(b / (a + b + 1)).
  cases <- list(
    list(shapes = c(1e308, 1e308), expected = -1.3862943611198906),
    list(shapes = c(1.5e308, 1.5e308), expected = -1.3862943611198906),
    list(shapes = c(1.7e308, 1e307), expected = -2.9475301717361133),
    list(shapes = c(1e307, 1.7e308), expected = -2.9475301717361133)
  )
  for (case in cases) {
    out <- .glmm_binom_aghq(
      1L, 0L, 1L, 1L, matrix(0), matrix(0), NULL,
      c(alpha = case$shapes[[1L]], beta = case$shapes[[2L]])
    )
    expect_equal(out[["value"]][1L, 1L], case$expected, tolerance = 1e-12,
      info = paste("shapes", case$shapes[[1L]], case$shapes[[2L]]))
    expect_equal(out[["exact_count"]], 1L)
  }
})


test_that("exact beta-binomial shortcut equals lbeta differences for moderate shapes", {

  ai  <- c(3L, 0L, 12L)
  ci  <- c(2L, 1L, 0L)
  n1i <- c(10L, 4L, 40L)
  n2i <- c(8L, 9L, 0L)
  for (shapes in list(c(.7, 2.3), c(1e-3, .5), c(45.5, 3.25), c(2, 2))) {
    for (weights in list(NULL, c(.5, 1.5, 2.5))) {
      out <- .glmm_binom_aghq(ai, ci, n1i, n2i, matrix(0, 1L, 3L),
        matrix(0, 1L, 3L), weights, c(alpha = shapes[[1L]], beta = shapes[[2L]]))
      power <- if (is.null(weights)) rep(1, 3L) else weights
      observed <- ai + ci
      expected <- power * (lchoose(n1i, ai) + lchoose(n2i, ci)) +
        lbeta(shapes[[1L]] + power * observed,
              shapes[[2L]] + power * (n1i + n2i - observed)) -
        lbeta(shapes[[1L]], shapes[[2L]])
      expect_equal(out[["value"]][1L, ], expected, tolerance = 1e-12,
        info = paste("shapes", shapes[[1L]], shapes[[2L]]))
    }
  }
})


test_that("exact beta-binomial shortcut keeps the log-gamma increments below 16", {

  # Increments log Gamma(z + t) - log Gamma(z) with z < 16 (z = alpha, beta and
  # alpha + beta: all six shape pairs have alpha + beta < 16) are the plain
  # difference of two log-gamma values. Counts (a, c, n1, n2) and likelihood
  # weights w give t = w (a + c) from 0.13 to 5000 and w (n1 + n2 - a - c) from
  # 0.97 to 13000. 60-digit mpmath values of
  #   w (log C(n1, a) + log C(n2, c)) + D(alpha, e) + D(beta, n - e) -
  #   D(alpha + beta, n),  D(z, t) = lgamma(z + t) - lgamma(z),
  # e = w (a + c), n = w (n1 + n2); compared per cell with an absolute
  # tolerance of 1e-13 for values below 1 and a relative one above.
  ai <- c(1L, 3L, 8L, 3L, 3L, 8L, 3L, 3L)
  ci <- c(0L, 2L, 5L, 2L, 2L, 5L, 2L, 2L)
  n1i <- c(1L, 10L, 60L, 10L, 10L, 60L, 10L, 10L)
  n2i <- c(1L, 8L, 50L, 8L, 8L, 50L, 8L, 8L)
  weights <- c(1, 1, 1, .5, 2.5, .01, 50, 1000)
  cases <- list(
    list(shapes = c(0.5, 0.5), value = c(
      -2.0794415416798359283, -4.2025887170617863969, -6.3128805713994477185,
      -2.6143108468821985218, -8.4242489408065669333, -0.59572108178606644264,
      -129.40052109612378402, -2520.5890510424101252)),
    list(shapes = c(1, 1), value = c(
      -1.7917594692280550008, -3.8805323383367752382, -6.2969490466446527857,
      -2.3092266074910831562, -8.090832744028430503, -0.49942413195023106158,
      -129.05938752016674643, -2520.2475202181369027)),
    list(shapes = c(2.5, 7.3), value = c(
      -1.7577634397781997241, -3.0870022605110920958, -5.3748848887333106882,
      -1.620961314612767903, -7.2123531844171271289, -0.13207092068295665132,
      -128.11383179058469678, -2519.2982584250863632)),
    list(shapes = c(0.01, 15.5), value = c(
      -7.409781397183374506, -10.377787518691391299, -10.446505795696202998,
      -8.0629461241824517361, -15.331089520943065029, -2.6982201646095804546,
      -137.03395998320866366, -2528.2679688944615381)),
    list(shapes = c(14, 1.5), value = c(
      -2.4996779671083128661, -10.818985324643229082, -25.616185664127602963,
      -7.2221976479569756632, -17.491831211851447491, -1.9325410176264752254,
      -141.53357272390982945, -2532.9448118143548388)),
    list(shapes = c(15.9, 0.05), value = c(
      -5.829139827385117962, -17.443877587396626256, -35.071548364923925978,
      -12.99010730939690261, -25.075248598703094999, -5.3115793696306189951,
      -150.33448996348256676, -2541.8392555745124572))
  )
  for (case in cases) {
    out <- .glmm_binom_aghq(ai, ci, n1i, n2i, matrix(0, 1L, 8L),
      matrix(0, 1L, 8L), weights,
      c(alpha = case$shapes[[1L]], beta = case$shapes[[2L]]))
    expect_equal(out[["exact_count"]], 8L)
    error <- abs(out[["value"]][1L, ] - case$value) / pmax(1, abs(case$value))
    expect_lt(max(error), 1e-13,
      label = paste("error at shapes", case$shapes[[1L]], case$shapes[[2L]]))
  }
})


test_that("exact beta-binomial shortcut is a normalized distribution", {

  n <- 6L
  for (shapes in list(c(.3, .7), c(1e12, 3e12), c(1e20, 1e20))) {
    out <- .glmm_binom_aghq(0:n, rep(0L, n + 1L), rep(n, n + 1L),
      rep(0L, n + 1L), matrix(0, 1L, n + 1L), matrix(0, 1L, n + 1L), NULL,
      c(alpha = shapes[[1L]], beta = shapes[[2L]]))
    expect_equal(sum(exp(out[["value"]][1L, ])), 1, tolerance = 1e-12)
  }
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

  expect_equal(out[["value"]][1, 1], expected[["value"]], tolerance = 2e-6,
               info = "independent joint binomial integral")
  expect_identical(expected[["outer"]][["message"]], "OK")
  expect_true(all(vapply(expected[["inner"]], function(x) x$integral$message == "OK", logical(1))))
})


test_that("active binomial boundary, sharp, point and truncated routes retain physical references", {

  # Fixed independent physical-p/theta integrals use full p support, positive
  # powered coefficients once and theta[-10,10] without renormalization. Omitted
  # normal mass <=1.524e-23. All final integration messages were OK; QUADPACK
  # estimates (and shared Rmath beta/CDF evaluation) are not rounding bounds.
  # Estimated log errors: zero/all6.78e-11, power2.33e-11, sharp6.83e-12,
  # point1.16e-14, truncated1.42e-11. Sharp observed rounding exceeded its
  # quadrature estimate; the existing absolute comparison scale is2e-6.
  cases <- list(
    zero = list(a = 0L, c = 0L, n1 = 10L, n2 = 10L, mu = .35, tau = .7,
      shapes = c(.01, 2.5), power = 1, reference = -.024476521042923388, grid = TRUE),
    all = list(a = 10L, c = 10L, n1 = 10L, n2 = 10L, mu = .35, tau = .7,
      shapes = c(2.5, .01), power = 1, reference = -.024476521042923388, grid = TRUE),
    power = list(a = 7L, c = 3L, n1 = 24L, n2 = 21L, mu = .35, tau = .7,
      shapes = c(1.5, 2.5), power = .8, reference = -4.2314458086137074, grid = FALSE),
    sharp = list(a = 505L, c = 499L, n1 = 88391L, n2 = 88391L,
      mu = -.6308474, tau = .373213, shapes = c(1, 1), power = 1,
      reference = -18.992479602352965, grid = FALSE),
    point = list(a = 7L, c = 3L, n1 = 24L, n2 = 21L, mu = .35, tau = .7,
      point = .4, power = .8, reference = -5.7572219676012706, grid = FALSE),
    truncated = list(a = 7L, c = 3L, n1 = 24L, n2 = 21L, mu = .35, tau = .7,
      shapes = c(.01, 2.5), upper = .9, power = 1,
      reference = -8.8888302714960936, grid = TRUE)
  )
  for (name in names(cases)) {
    x <- cases[[name]]
    prior <- if (!is.null(x$point)) BayesTools::prior("point", list(x$point)) else
      BayesTools::prior("beta", as.list(x$shapes),
        truncation = list(lower = 0, upper = if (is.null(x$upper)) 1 else x$upper))
    out <- .outcome_pdf.binom(x$a, x$c, x$n1, x$n2,
      matrix(x$mu), matrix(x$tau), prior, weights = x$power)
    diagnostics <- attr(out, "glmm_aghq_diagnostics", exact = TRUE)
    expect_lt(abs(out[1L, 1L] - x$reference), 2e-6, label = name)
    expect_identical(diagnostics$exact_count, 0L, info = name)
    if (x$grid) {
      expect_identical(diagnostics$grid_columns, 1L, info = name)
      expect_identical(sum(diagnostics$order_counts), 0L, info = name)
      expect_gte(diagnostics$grid_max_theta_order, 29L, label = name)
      expect_gte(diagnostics$grid_max_nuisance_order, 75L, label = name)
      expect_lte(diagnostics$grid_max_change, 1e-6, label = name)
    } else {
      expect_length(diagnostics$grid_columns, 0L)
      expect_identical(sum(diagnostics$order_counts), 1L, info = name)
      expect_gt(diagnostics$max_order, 0L, label = name)
      expect_lte(diagnostics$max_change, 1e-6, label = name)
      if (name == "point") {
        # The point route is the R one-dimensional evaluator.
        direct <- .glmm_binom_point_aghq(x$a, x$c, x$n1, x$n2,
          matrix(x$mu), matrix(x$tau), x$power, x$point)
        expect_identical(diagnostics$max_mode_iterations, 55L)
      } else {
        direct <- .glmm_binom_aghq(x$a, x$c, x$n1, x$n2,
          matrix(x$mu), matrix(x$tau), x$power,
          c(alpha = x$shapes[1L], beta = x$shapes[2L]))
      }
      expect_lt(abs(direct$value[1L, 1L] - x$reference), 2e-6, label = name)
      expect_identical(direct$exact_count, 0L)
      expect_identical(sum(direct$order_counts), 1L)
    }
    if (name %in% c("zero", "all")) {
      expect_error(.glmm_binom_aghq(x$a, x$c, x$n1, x$n2,
        matrix(x$mu), matrix(x$tau), x$power,
        c(alpha = x$shapes[1L], beta = x$shapes[2L])),
        "Binomial AGHQ was rejected by diagnostics.*quadrature order 49",
        class = "RoBMA_glmm_aghq_nonconvergence")
    }
  }
})


test_that("fitted multilevel GLMM estimate scores retain cluster and scale coordinates", {

  object <- load_fit("bcg_glmm_3lvl_scale")
  posterior <- .get_posterior_samples(object$fit)
  setup <- .log_lik_posterior_setup(object$fit, posterior, object$data,
    object$priors, unit = "estimate", data_hash = .get_outcome_hash(object))
  design <- attr(object$fit, "formula_design", exact = TRUE)[["log_tau"]]
  mapping <- BayesTools::JAGS_formula_name_map(object$fit, "log_tau")
  mapping <- mapping[mapping$kind == "fixed" & mapping$role == "coefficient", ]
  coefficients <- as.numeric(posterior[1L, mapping$jags_name[match(design$model_terms, mapping$term)]])
  coefficients[design$model_terms == "intercept"] <- log(coefficients[design$model_terms == "intercept"])
  total <- exp(sum(coefficients * design$model_matrix[1L, ]))
  rho <- posterior[1L, "rho"]
  data <- .create_fit_data(object$data, object$priors)
  offset <- posterior[1L, paste0("gamma[", data$cluster[1L], "]")] * total * sqrt(rho)
  mu <- posterior[1L, "mu"] + offset
  within <- total * sqrt(1 - rho)
  expect_lt(abs(mu - setup$mu[1L, 1L]), 1e-12)
  expect_lt(abs(within - setup$tau_within[1L, 1L]), 1e-12)
  score <- log_lik(object, unit = "estimate")
  diagnostics <- attr(score, "glmm_aghq_diagnostics", exact = TRUE)
  expect_identical(dim(score), dim(setup$mu))
  expect_true(all(is.finite(score)))
  expect_identical(diagnostics$exact_count, 0L)
  expect_length(diagnostics$grid_columns, 0L)
  expect_identical(sum(diagnostics$order_counts), length(score))
  first <- object$data$outcome[1L, ]
  prior <- object$priors$outcome$pi$parameters
  direct <- .glmm_binom_aghq(first$ai, first$ci, first$n1i, first$n2i,
    matrix(mu), matrix(within), NULL, c(alpha = prior$alpha, beta = prior$beta))
  expect_lt(abs(score[1L, 1L] - direct$value[1L, 1L]), 1e-10)
  expect_error(log_lik(object, unit = "cluster"), "Cluster-unit GLMM log-likelihood is unavailable")
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


test_that("ordinary wrappers dispatch supported nuisance priors", {

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
  expected_bin  <- bin_raw[["value"]]
  expected_pois <- pois_raw[["value"]]

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
  expect_true(is.finite(
    .outcome_pdf.pois(
      7L, 11L, 15, 18, matrix(0.35), matrix(0.6), truncated_phi
    )[1L, 1L]
  ), info = "truncated Poisson prior uses refined prior-CDF quadrature")
  expect_true(is.finite(
    .outcome_pdf.binom(
      7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), truncated_pi
    )[1L, 1L]
  ), info = "truncated binomial prior uses refined prior-CDF quadrature")
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

  expect_equal(pointwise[1L, 2L], expected_second[["value"]], tolerance = 2e-6,
               info = "failed observation uses independently checked grid")
  expect_identical(diagnostics[["grid_columns"]], 2L,
                   info = "the exact first observation does not fall back")
  expect_equal(as.numeric(summed), rowSums(pointwise), tolerance = 1e-12,
               info = "pointwise and sum routes share recovered values")
})


test_that("point nuisance priors collapse their quadrature dimension", {

  prior_pi  <- BayesTools::prior("point", list(0.4))
  prior_phi <- BayesTools::prior("point", list(-1))
  bin <- .outcome_pdf.binom(
    7L, 3L, 24L, 21L, matrix(0.35), matrix(0), prior_pi
  )
  pois <- .outcome_pdf.pois(
    7L, 11L, 15, 18, matrix(0.35), matrix(0), prior_phi
  )
  expected_bin <- stats::dbinom(
    7L, 24L, stats::plogis(stats::qlogis(0.4) + 0.35 / 2), log = TRUE
  ) + stats::dbinom(
    3L, 21L, stats::plogis(stats::qlogis(0.4) - 0.35 / 2), log = TRUE
  )
  expected_pois <- stats::dpois(
    7L, 15 * exp(-1 + 0.35 / 2), log = TRUE
  ) + stats::dpois(
    11L, 18 * exp(-1 - 0.35 / 2), log = TRUE
  )

  expect_equal(bin[1L, 1L], expected_bin, tolerance = 1e-12)
  expect_equal(pois[1L, 1L], expected_pois, tolerance = 1e-12)
  expect_identical(
    attr(bin, "glmm_aghq_diagnostics")[["grid_columns"]],
    integer(0)
  )
  expect_identical(
    attr(pois, "glmm_aghq_diagnostics")[["grid_columns"]],
    integer(0)
  )
})


test_that("point nuisance AGHQ centers active heterogeneity integrals", {

  scaled_integral <- function(log_kernel, negative_hessian) {
    mode <- stats::optimize(
      log_kernel,
      interval = c(-12, 12),
      maximum  = TRUE,
      tol      = 1e-12
    )
    hessian <- negative_hessian(mode[["maximum"]])
    scaled <- stats::integrate(
      function(z) {
        exp(
          log_kernel(mode[["maximum"]] + z / sqrt(hessian)) -
            mode[["objective"]]
        ) / sqrt(hessian)
      },
      lower         = -Inf,
      upper         = Inf,
      subdivisions  = 2000L,
      rel.tol       = 1e-11,
      stop.on.error = TRUE
    )[["value"]]
    mode[["objective"]] + log(scaled)
  }

  a   <- 505L
  c   <- 499L
  n1  <- 88391L
  n2  <- 88391L
  mu  <- -0.6308474
  tau <- 0.373213
  pi  <- 0.5
  bin_log_kernel <- function(theta) {
    effect <- mu + tau * theta
    stats::dbinom(
      a, n1, stats::plogis(stats::qlogis(pi) + 0.5 * effect), log = TRUE
    ) + stats::dbinom(
      c, n2, stats::plogis(stats::qlogis(pi) - 0.5 * effect), log = TRUE
    ) + stats::dnorm(theta, log = TRUE)
  }
  bin_hessian <- function(theta) {
    effect <- mu + tau * theta
    p1     <- stats::plogis(stats::qlogis(pi) + 0.5 * effect)
    p2     <- stats::plogis(stats::qlogis(pi) - 0.5 * effect)
    1 + 0.25 * tau^2 * (n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2))
  }
  bin <- .glmm_binom_point_aghq(
    a, c, n1, n2, matrix(mu), matrix(tau), NULL, pi
  )
  expected_bin <- scaled_integral(bin_log_kernel, bin_hessian)

  x1  <- 7L
  x2  <- 11L
  t1  <- 15
  t2  <- 18
  mu  <- 0.35
  tau <- 0.6
  phi <- -1
  pois_log_kernel <- function(theta) {
    effect  <- mu + tau * theta
    lambda1 <- t1 * exp(phi + 0.5 * effect)
    lambda2 <- t2 * exp(phi - 0.5 * effect)
    stats::dpois(x1, lambda1, log = TRUE) +
      stats::dpois(x2, lambda2, log = TRUE) +
      stats::dnorm(theta, log = TRUE)
  }
  pois_hessian <- function(theta) {
    effect  <- mu + tau * theta
    lambda1 <- t1 * exp(phi + 0.5 * effect)
    lambda2 <- t2 * exp(phi - 0.5 * effect)
    1 + 0.25 * tau^2 * (lambda1 + lambda2)
  }
  pois <- .glmm_pois_point_aghq(
    x1, x2, t1, t2, matrix(mu), matrix(tau), NULL, phi
  )

  expect_equal(
    bin[["value"]][1L, 1L],
    expected_bin,
    tolerance = 1e-7,
    info = "sharp BCG posterior row"
  )
  expect_equal(
    pois[["value"]][1L, 1L],
    scaled_integral(pois_log_kernel, pois_hessian),
    tolerance = 1e-7,
    info = "Poisson active-heterogeneity row"
  )
  expect_lte(bin[["max_order"]], 13L)
  expect_lte(pois[["max_order"]], 13L)
})


test_that("point nuisance AGHQ fails when refinement cannot certify accuracy", {

  control <- .glmm_aghq_control(
    orders      = c(5L, 9L, 13L),
    tolerance   = 1e-15,
    consecutive = 2L
  )
  expect_error(
    .glmm_binom_point_aghq(
      4L, 11L, 123L, 139L,
      matrix(-0.3614195), matrix(1.45076), NULL, 0.5,
      control = control
    ),
    "Binomial AGHQ was rejected by diagnostics.*sample 1, observation 1",
    class = "RoBMA_glmm_aghq_nonconvergence",
    info = "capped refinement does not silently accept the last rule"
  )
})


test_that("point Poisson AGHQ retains zero-count finite log-rate tails", {

  for (counts in list(c(0L, 0L), c(1L, 0L), c(0L, 1L))) {
    out <- .glmm_pois_point_aghq(counts[1L], counts[2L], 1, 1,
      matrix(0), matrix(.7), .8, -800)
    # exp(-800 +/- .35*theta) is negligible at this comparison scale.
    # A single event gives the normal MGF exp((.8*.35)^2/2).
    expected <- if (sum(counts) == 0L) 0 else -.8 * 800 + (.8 * .35)^2 / 2
    expect_lt(abs(out[["value"]][1L, 1L] - expected), 1e-10)
    expect_identical(out[["exact_count"]], 0L)
    expect_identical(sum(out[["order_counts"]]), 1L)
    expect_identical(out[["max_mode_iterations"]], 55L)
  }
})


test_that("GLMM refinements reject persistent non-finite ordinates", {

  grid_control <- .glmm_grid_control(
    theta_orders    = c(3L, 5L, 7L),
    nuisance_orders = c(3L, 5L, 7L),
    consecutive     = 2L
  )
  expect_error(
    .glmm_grid_refine(
      evaluate     = function(theta_order, nuisance_order) -Inf,
      outcome_type = "bin",
      observation  = 1L,
      control      = grid_control
    ),
    "prior-CDF quadrature was rejected by diagnostics",
    class = "RoBMA_glmm_grid_nonconvergence",
    info = "equal non-finite grid ordinates are numerical failure"
  )

  point_control <- .glmm_aghq_control(
    orders      = c(3L, 5L, 7L),
    consecutive = 2L
  )
  expect_error(
    .glmm_point_aghq_refine(
      log_density      = function(theta) theta * 0 - Inf,
      gradient         = function(theta) theta * 0,
      negative_hessian = function(theta) theta * 0 + 1,
      n                = 1L,
      outcome_type     = "bin",
      observation      = 1L,
      control          = point_control
    ),
    "AGHQ was rejected by diagnostics",
    class = "RoBMA_glmm_aghq_nonconvergence",
    info = "equal non-finite point-prior ordinates are numerical failure"
  )
})


test_that("a fitted point-baserate GLMM supports log-likelihood criteria", {

  fit <- load_fit("lane_glmm_point_baserate")
  log_lik <- log_lik(fit)
  fit_loo <- suppressWarnings(add_loo(fit))
  fit_waic <- suppressWarnings(add_waic(fit))

  expect_gt(fit[["summary"]]["tau", "Mean"], 0)
  expect_identical(dim(log_lik), c(100L, 3L))
  expect_true(all(is.finite(log_lik)))
  expect_s3_class(loo(fit_loo), "loo")
  expect_s3_class(waic(fit_waic), "waic")
})


test_that("boundary-concentrated beta priors support likelihood criteria", {

  fit <- load_fit("lane_glmm_boundary_beta")
  log_lik <- log_lik(fit)
  fit_loo <- suppressWarnings(add_loo(fit))
  fit_waic <- suppressWarnings(add_waic(fit))

  expect_identical(dim(log_lik), c(100L, 3L))
  expect_true(all(is.finite(log_lik)))
  expect_s3_class(loo(fit_loo), "loo")
  expect_s3_class(waic(fit_waic), "waic")
})


test_that("boundary-matched beta quadrature integrates the normalized prior", {

  interval_mass <- function(alpha, beta, lower, upper) {

    if (lower == 0) {
      return(stats::pbeta(upper, alpha, beta))
    }
    if (upper == 1) {
      return(stats::pbeta(lower, alpha, beta, lower.tail = FALSE))
    }
    return(stats::pbeta(upper, alpha, beta) -
             stats::pbeta(lower, alpha, beta))
  }
  expected_zero_probability <- function(alpha, beta, lower, upper) {

    moment <- exp(lbeta(alpha, beta + 2) - lbeta(alpha, beta))
    return(moment *
      interval_mass(alpha, beta + 2, lower, upper) /
      interval_mass(alpha, beta, lower, upper))
  }

  cases <- list(
    list(alpha = 0.01, beta = 0.01, lower = 0, upper = 1),
    list(alpha = 0.01, beta = 2.5, lower = 0, upper = 0.9),
    list(alpha = 2.5, beta = 0.01, lower = 0.1, upper = 1)
  )
  for (case in cases) {
    prior_pi <- BayesTools::prior(
      "beta",
      list(case[["alpha"]], case[["beta"]]),
      truncation = list(case[["lower"]], case[["upper"]])
    )
    observed <- .outcome_pdf.binom(
      0L, 0L, 1L, 1L, matrix(0), matrix(0), prior_pi,
      n_theta = 1L, n_pi = 25L
    )[1L, 1L]
    expected <- log(expected_zero_probability(
      case[["alpha"]], case[["beta"]],
      case[["lower"]], case[["upper"]]
    ))

    expect_equal(observed, expected, tolerance = 1e-11)
  }
})


test_that("truncated nuisance priors match independent integration", {

  prior_pi <- BayesTools::prior(
    "beta", list(1.5, 2.5), truncation = list(0.1, 0.9)
  )
  prior_phi <- BayesTools::prior(
    "normal", list(-1, 1.3), truncation = list(-3, 2)
  )
  bin <- .outcome_pdf.binom(
    7L, 3L, 24L, 21L, matrix(0.35), matrix(0.7), prior_pi
  )
  pois <- .outcome_pdf.pois(
    7L, 11L, 15, 18, matrix(0.35), matrix(0.6), prior_phi
  )
  expected_bin <- .binom_aghq_oracle(
    a = 7L, c = 3L, n1 = 24L, n2 = 21L,
    mu = 0.35, tau = 0.7, alpha = 1.5, beta = 2.5,
    lower = 0.1, upper = 0.9
  )
  expected_pois <- .pois_aghq_oracle(
    x1 = 7L, x2 = 11L, t1 = 15, t2 = 18,
    mu = 0.35, tau = 0.6, mean = -1, sd = 1.3,
    lower = -3, upper = 2
  )

  expect_equal(bin[1L, 1L], expected_bin[["value"]], tolerance = 2e-6,
               info = "truncated beta prior")
  expect_equal(pois[1L, 1L], expected_pois, tolerance = 2e-6,
               info = "truncated normal prior")
})


test_that("truncated nuisance recovery certifies sharp theta integrals", {

  prior_pi <- BayesTools::prior(
    "beta", list(1.5, 2.5), truncation = list(0.1, 0.9)
  )
  bin <- .outcome_pdf.binom(
    6L, 29L, 306L, 303L,
    matrix(-0.4289519), matrix(1.018359), prior_pi
  )
  expected_bin <- .binom_aghq_oracle(
    a = 6L, c = 29L, n1 = 306L, n2 = 303L,
    mu = -0.4289519, tau = 1.018359,
    alpha = 1.5, beta = 2.5, lower = 0.1, upper = 0.9
  )

  prior_phi <- BayesTools::prior(
    "normal", list(-1, 1.3), truncation = list(-3, 2)
  )
  pois <- .outcome_pdf.pois(
    7L, 11L, 18, 21,
    matrix(-0.2936048), matrix(0.9191402), prior_phi
  )
  expected_pois <- .pois_aghq_oracle(
    x1 = 7L, x2 = 11L, t1 = 18, t2 = 21,
    mu = -0.2936048, tau = 0.9191402,
    mean = -1, sd = 1.3, lower = -3, upper = 2
  )

  expect_equal(bin[1L, 1L], expected_bin[["value"]], tolerance = 2e-6)
  expect_equal(pois[1L, 1L], expected_pois, tolerance = 2e-6)
  expect_gt(
    attr(bin, "glmm_aghq_diagnostics")[["grid_max_theta_order"]],
    49L
  )
  expect_gt(
    attr(pois, "glmm_aghq_diagnostics")[["grid_max_theta_order"]],
    49L
  )
})


test_that("fitted truncated-nuisance GLMMs support likelihood criteria", {

  fit_bin <- load_fit("lane_glmm_truncated_beta")
  fit_pois <- load_fit("lane_glmm_truncated_lograte")

  for (fit in list(binomial = fit_bin, Poisson = fit_pois)) {
    log_lik <- log_lik(fit)
    fit_loo <- suppressWarnings(add_loo(fit))
    fit_waic <- suppressWarnings(add_waic(fit))

    expect_identical(dim(log_lik), c(100L, 3L))
    expect_true(all(is.finite(log_lik)))
    expect_s3_class(loo(fit_loo), "loo")
    expect_s3_class(waic(fit_waic), "waic")
  }
})


test_that("quadrature recovery catches only its classed failures", {

  grid_failure <- .glmm_grid_nonconvergence_condition(
    outcome_type     = "bin",
    observation      = 7L,
    theta_changes    = c(2e-5, 3e-6),
    nuisance_changes = c(4e-7, 5e-8)
  )
  expect_true(.glmm_grid_numerical_failure(grid_failure, "bin"))
  expect_false(.glmm_grid_numerical_failure(grid_failure, "pois"))
  expect_identical(grid_failure[["observation"]], 7L)
  expect_equal(grid_failure[["theta_changes"]], c(2e-5, 3e-6))
  expect_match(
    conditionMessage(grid_failure),
    "^Binomial prior-CDF quadrature was rejected by diagnostics at observation 7:"
  )

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

  matching_plain_error <- conditionMessage(grid_failure)
  testthat::local_mocked_bindings(
    .outcome_pdf.binom = function(...) stop(matching_plain_error),
    .package = "RoBMA"
  )
  prior_pi <- BayesTools::prior(
    "beta", list(1.5, 2.5), truncation = list(0.1, 0.9)
  )
  expect_error(
    .glmm_binom_marginal(
      ai = 1L, ci = 1L, n1i = 10L, n2i = 10L,
      mu_samples = matrix(0), tau_within = matrix(0.5),
      weights = NULL, prior_pi = prior_pi, prior_spec = NULL
    ),
    "observation 7"
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


test_that("AGHQ production refinement remains capped at order 49", {

  rules <- .glmm_aghq_control()[["rules"]][["nodes"]]
  expect_identical(length(rules[[length(rules)]]), 49L)
})


test_that("AGHQ reports failure when refinement is capped", {

  control <- .glmm_aghq_control(
    orders      = c(5L, 9L, 13L),
    tolerance   = 1e-15,
    consecutive = 2L
  )

  condition <- tryCatch(.glmm_pois_aghq(
    1L, 1L, 370, 483,
    matrix(-0.45), matrix(0.5), NULL,
    c(mean = -6.0187, sd = 20.27),
    control = control
  ), error = identity)
  expect_s3_class(condition, "RoBMA_glmm_aghq_nonconvergence")
  expect_s3_class(condition, "RoBMA_glmm_numerical_error")
  expect_identical(condition[["outcome_type"]], "pois")
  expect_identical(condition[["reason"]], "quadrature_nonconvergence")
  expect_identical(condition[["sample"]], 1L)
  expect_identical(condition[["observation"]], 1L)
  expect_identical(condition[["order"]], 13L)
  expect_true(is.finite(condition[["change"]]))
  expect_identical(conditionMessage(condition), paste0(
    "Poisson AGHQ was rejected by diagnostics: last log-likelihood change ",
    format(condition[["change"]], digits = 8L),
    " at quadrature order 13 at sample 1, observation 1. ",
    "Inspect the data and 'prior_lograte' before refitting the model."))
  expect_identical(condition[["noise"]], 0)
})


test_that("AGHQ refusals name the rounding-noise bound only where it set the change", {

  wording <- function(change, noise) conditionMessage(.glmm_numerical_condition(
    "bin", "aghq", "quadrature_nonconvergence", 1L, sample = 1L, order = 13L,
    change = change, noise = noise))
  # no bound (Poisson, fixed base rate, closed-form beta route): the change
  for (noise in list(NULL, 0)) {
    expect_match(wording(3e-7, noise), "AGHQ was rejected by diagnostics: last log-likelihood change 3e-07 at quadrature order 13",
      fixed = TRUE)
  }
  # a bound above the difference of the last two orders: the bound
  expect_match(wording(3e-6, 3e-6), "diagnostics: the rounding-noise bound of the prior density was 3e-06 at quadrature order 13",
    fixed = TRUE)
  # a difference above the bound: the change
  expect_match(wording(3e-6, 3e-7), "diagnostics: last log-likelihood change 3e-06 at quadrature order 13",
    fixed = TRUE)

  # Through the native route: a beta base rate on the density route (shape sum
  # above 1e3) whose refinement cap is reached before the bound matters keeps
  # the change; the point base-rate route has no bound.
  control <- .glmm_aghq_control(orders = c(3L, 5L, 7L), tolerance = 1e-15, consecutive = 2L)
  condition <- tryCatch(.glmm_binom_aghq(1L, 1L, 2L, 2L, matrix(-.45), matrix(.5), NULL,
    c(alpha = 2e4, beta = 3e4), control = control), error = identity)
  expect_s3_class(condition, "RoBMA_glmm_aghq_nonconvergence")
  expect_gt(condition[["noise"]], 0)
  expect_gt(condition[["change"]], condition[["noise"]])
  expect_match(conditionMessage(condition), "diagnostics: last log-likelihood change ", fixed = TRUE)
  condition <- tryCatch(.glmm_binom_point_aghq(1L, 1L, 2L, 2L, matrix(-.45), matrix(.5), NULL, .3,
    control = control), error = identity)
  expect_s3_class(condition, "RoBMA_glmm_aghq_nonconvergence")
  expect_match(conditionMessage(condition), "diagnostics: last log-likelihood change ", fixed = TRUE)
})
