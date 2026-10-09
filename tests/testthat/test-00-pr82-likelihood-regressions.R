.cluster_norm_dispatch_setup <- function(yi, sei, tau_within, tau_between,
                                         weights = NULL, mu = 0,
                                         effect_direction = "positive") {

  K <- length(yi)
  list(
    outcome_type      = "norm",
    is_weightfunction = FALSE,
    effect_direction  = effect_direction,
    data              = list(),
    yi                = yi,
    sei               = sei,
    S                 = 1L,
    K                 = K,
    mu                = matrix(mu, 1L, K),
    tau_within        = matrix(tau_within, 1L, K),
    tau_between       = matrix(tau_between, 1L, K),
    cluster           = list(seq_len(K)),
    weights           = weights
  )
}

test_that("powered normal cluster likelihood integrates the latent effect in closed form", {

  # Reference values: 80-digit mpmath evaluation of
  #   log int prod_i N(yi_i; mu_i + tau_i g, sd_i)^w_i N(g; 0, 1) dg
  # through the completed square in the unit-normal latent.
  cases <- list(
    # a dominant latent effect against two precise rows; Gauss-Hermite on the
    # prior scale returned -186.4 for the powers 2 and 2, -93.77 for 1 and 1
    list(yi = c(.1, .1), sei = c(.01, .01), tau_within = 0, tau_between = 1,
      weights = c(2, 2), expected = 9.441596869738760),
    list(yi = c(.1, .1), sei = c(.01, .01), tau_within = 0, tau_between = 1,
      weights = c(1, 1), expected = 2.415694779911253),
    list(yi = c(.3, -.2, .5), sei = c(.1, .2, .05), tau_within = 0,
      tau_between = .7, weights = c(.5, 1, 1.5), expected = -5.184067167560930),
    # row-specific loadings (a scale model) and unequal powers
    list(yi = c(.3, -.2, .5, .1), sei = c(.1, .2, .05, .3), tau_within = 0,
      tau_between = c(.7, .35, 0, 1.4), weights = c(.5, 1, 1.5, 2),
      expected = -73.95729500926897),
    list(yi = c(.3, -.2, .5), sei = c(.1, .2, .05), tau_within = 0,
      tau_between = c(.7, -.35, 1.4), weights = c(1, 1, 1),
      expected = 0.5355236481055064),
    list(yi = c(1e-3, 2e-3, -4e-3, 5e-4), sei = c(1e-6, 1e-3, 1e-9, 2e-5),
      tau_within = 0, tau_between = c(3, 1e-4, 2e5, 0),
      weights = c(1, 2, .5, 1), expected = -500364.4073559806),
    # a unit information of 1e800 and a representable density
    list(yi = 0, sei = 1e-200, tau_within = 0, tau_between = 1e200,
      weights = 1, expected = -461.4359571320138),
    # the same row next to a zero-power companion, which is omitted
    list(yi = c(0, 7), sei = c(1e-200, 3), tau_within = 0,
      tau_between = c(1e200, 2), weights = c(1, 0), expected = -461.4359571320138)
  )
  for (case in cases) {
    setup <- do.call(.cluster_norm_dispatch_setup, case[names(case) != "expected"])
    expected <- matrix(case$expected, 1L, 1L)
    expect_equal(.log_lik_cluster_from_setup(setup), expected, tolerance = 1e-13)
    expect_equal(.log_lik_cluster_sum_from_setup(setup), as.vector(expected),
      tolerance = 1e-13)
    expect_equal(.log_lik_cluster_norm_analytic_r(setup, setup$yi, setup$sei),
      expected, tolerance = 1e-13)
  }

  # A cluster without a positive power contributes nothing.
  setup <- .cluster_norm_dispatch_setup(c(.1, .2), c(.1, .1), 0, 1,
    weights = c(0, 0))
  expect_identical(.log_lik_cluster_from_setup(setup), matrix(0, 1L, 1L))
  expect_identical(.log_lik_cluster_norm_analytic_r(setup, setup$yi, setup$sei),
    matrix(0, 1L, 1L))
})

test_that("powered cluster likelihood agrees with numerical integration and the allocation grid", {

  S <- 3L
  K <- 5L
  set.seed(20261009)
  setup <- list(
    S           = S,
    mu          = matrix(rnorm(S * K, 0, .25), S, K),
    tau_within  = matrix(runif(S * K, .03, .25), S, K),
    tau_between = matrix(runif(S * K, .02, .20), S, K),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L, 5L)),
    weights     = c(1, .5, 1.25, 2, .75)
  )
  yi  <- c(.02, .15, -.05, .30, .08)
  sei <- c(.10, .12, .08, .20, .15)

  integrated <- matrix(NA_real_, S, length(setup$cluster))
  for (s in seq_len(S)) {
    for (g in seq_along(setup$cluster)) {
      idx <- setup$cluster[[g]]
      integrand <- function(gamma) {
        vapply(gamma, function(value) {
          exp(sum(setup$weights[idx] * stats::dnorm(yi[idx],
            mean = setup$mu[s, idx] + value * setup$tau_between[s, idx],
            sd = sqrt(sei[idx]^2 + setup$tau_within[s, idx]^2), log = TRUE))) *
            stats::dnorm(value)
        }, numeric(1))
      }
      integrated[s, g] <- log(stats::integrate(integrand, -Inf, Inf,
        rel.tol = 1e-12, stop.on.error = TRUE)$value)
    }
  }
  expect_equal(.log_lik_cluster_norm_analytic(setup, yi, sei), integrated,
    tolerance = 1e-10)
  expect_equal(.log_lik_cluster_norm_analytic_r(setup, yi, sei), integrated,
    tolerance = 1e-10)
  expect_equal(.log_lik_cluster_norm_analytic_sum(setup, yi, sei),
    rowSums(integrated), tolerance = 1e-10)

  # The allocation grid applies the powers to the same cluster integral.
  grid <- list(S = S, mu = setup$mu, cluster = setup$cluster,
    weights = setup$weights,
    tau_total = matrix(runif(S * K, .1, .4), S, K))
  rho <- c(0, .37, 1)
  evaluated <- .log_lik_cluster_norm_analytic_rho_grid_sum(grid, yi, sei, rho)
  expected <- vapply(seq_len(S), function(s) vapply(rho, function(value) {
    allocation <- grid
    allocation$S <- 1L
    allocation$mu <- grid$mu[s, , drop = FALSE]
    allocation$tau_within <- grid$tau_total[s, , drop = FALSE] * sqrt(1 - value)
    allocation$tau_between <- grid$tau_total[s, , drop = FALSE] * sqrt(value)
    sum(.log_lik_cluster_norm_analytic_r(allocation, yi, sei))
  }, numeric(1)), numeric(length(rho)))
  expect_equal(evaluated, expected, tolerance = 1e-12)
})

test_that("the native cluster integral agrees with the R helper", {

  # Both evaluate the same two-regime form (standardized units when no loading
  # exceeds the smallest sd, units scaled by the largest loading otherwise), so
  # clusters of both regimes must agree to rounding.
  set.seed(20261010)
  S <- 4L
  G <- 150L
  sizes <- sample.int(5L, G, replace = TRUE)
  K <- sum(sizes)
  cluster <- split(seq_len(K), rep(seq_len(G), sizes))
  scale <- function(n, low, high) exp(stats::runif(n, log(low), log(high)))
  setup <- list(
    S           = S,
    mu          = matrix(stats::rnorm(S * K, 0, .3), S, K),
    tau_within  = matrix(scale(S * K, 1e-8, 1), S, K) *
      matrix(stats::rbinom(S * K, 1L, .5), S, K),
    tau_between = matrix(scale(S * K, 1e-8, 5), S, K) *
      matrix(stats::rbinom(S * K, 1L, .85), S, K),
    cluster     = cluster,
    weights     = sample(c(0, .5, 1, 2, 3), K, replace = TRUE, prob = c(.1, .2, .4, .2, .1))
  )
  yi  <- stats::rnorm(K, 0, .3)
  sei <- scale(K, 1e-8, 2)

  native <- .log_lik_cluster_norm_analytic_native(setup, yi, sei)
  log_form <- .log_lik_cluster_norm_analytic_r(setup, yi, sei)
  expect_true(all(is.finite(native)))
  expect_lte(max(abs(native - log_form) / pmax(abs(log_form), 1)), 1e-14)
  expect_lte(max(abs(
    .log_lik_cluster_norm_analytic_native(setup, yi, sei, row_sum = TRUE) -
      rowSums(log_form)) / pmax(abs(rowSums(log_form)), 1)), 1e-14)
})

test_that("cluster integrals over a wide range of loadings and sds match mpmath", {

  # Rows (.3, -.2) with sd s, powers (1, 2) and a common loading, from a
  # loading of 1e-100 times the sd to one of 1e160 times the sd (b = 1e160
  # cannot be squared in plain arithmetic). Reference values: 120-digit
  # mpmath.
  cases <- list(
    list(tau = 1e100, s = 1, expected = -233.6479643766859748),
    list(tau = 1e-100, s = 1, expected = -2.841815599614018225),
    list(tau = 1e-50, s = 1e50, expected = -348.1445795487208708),
    list(tau = 1e50, s = 1e-50, expected = -8.333333333333333333e98),
    list(tau = 1.1e100, s = 1, expected = -233.7432745564902997),
    list(tau = 9e-101, s = 1, expected = -2.841815599614018225),
    list(tau = 1e101, s = 1e50, expected = -466.1257254357512556),
    list(tau = 1e100, s = .5, expected = -232.5116700155660842),
    list(tau = 1e-50, s = 1e51, expected = -355.0523348277030079),
    list(tau = 10, s = 1e101, expected = -700.4400987768098605),
    list(tau = 1e80, s = 1e-80, expected = -8.333333333333333333e158)
  )
  for (case in cases) {
    setup <- list(S = 1L, mu = matrix(0, 1L, 2L), tau_within = matrix(0, 1L, 2L),
      tau_between = matrix(case$tau, 1L, 2L), cluster = list(1:2),
      weights = c(1, 2))
    sei <- rep(case$s, 2L)
    expect_equal(.log_lik_cluster_norm_analytic_native(setup, c(.3, -.2), sei),
      matrix(case$expected, 1L, 1L), tolerance = 1e-13,
      info = paste("tau", case$tau, "sd", case$s))
    expect_equal(.log_lik_cluster_norm_analytic_r(setup, c(.3, -.2), sei),
      matrix(case$expected, 1L, 1L), tolerance = 1e-13,
      info = paste("tau", case$tau, "sd", case$s))
  }
})

test_that("normal cluster likelihood refuses invalid rows without a substituted value", {

  setup <- .cluster_norm_dispatch_setup(c(.1, .2), c(.1, .1), 0, 1,
    weights = c(1, -1))
  expect_true(is.na(.log_lik_cluster_norm_analytic_native(
    setup, setup$yi, setup$sei)))
  expect_true(is.na(.log_lik_cluster_norm_analytic_r(setup, setup$yi, setup$sei)))
  setup <- .cluster_norm_dispatch_setup(c(.1, NA), c(.1, .1), 0, 1)
  expect_true(is.na(.log_lik_cluster_norm_analytic_r(setup, setup$yi, setup$sei)))

  # An invalid cluster does not spoil the others; a cluster without a positive
  # power contributes nothing.
  value <- .log_lik_cluster_norm_stable(residual = c(.1, NA, .3, .2),
    sd = rep(.1, 4L), loading = rep(1, 4L), weights = c(1, 1, 1, 0),
    group = c(1L, 1L, 2L, 3L))
  expect_true(is.na(value[1L]))
  expect_equal(value[2L], stats::dnorm(.3, sd = sqrt(.1^2 + 1), log = TRUE),
    tolerance = 1e-14)
  expect_identical(value[3L], 0)
})

test_that("normal cluster R reference integrates weighted rows without the native kernel", {

  setup <- list(S = 2L, cluster = list(1:2), mu = matrix(0, 2L, 2L),
    tau_within = matrix(0.2, 2L, 2L), tau_between = matrix(0, 2L, 2L),
    weights = c(2, 0.5))
  local_mocked_bindings(.has_native_norm_cluster_analytic = function(...) FALSE,
                        .package = "RoBMA")
  actual <- .log_lik_cluster_norm_analytic(setup, c(0.1, -0.3), c(0.3, 0.4))
  expected <- sum(dnorm(c(0.1, -0.3), 0, sqrt(c(0.3, 0.4)^2 + 0.2^2),
                        log = TRUE) * c(2, 0.5))
  expect_equal(actual, matrix(expected, 2L, 1L), tolerance = 1e-12)
  expect_equal(.log_lik_cluster_norm_analytic_sum(setup, c(0.1, -0.3),
    c(0.3, 0.4)), rep(expected, 2L), tolerance = 1e-12)
})

test_that("cluster-unit likelihoods of selection models never reach the normal kernels", {

  # A selection setup carries the fitted joint selection model, which the
  # dispatchers evaluate first. Without it the normal cluster kernels would
  # silently ignore the selection.
  calls <- 0L
  sentinel <- matrix(0, 1L, 1L)
  local_mocked_bindings(
    .log_lik_cluster_norm_analytic     = function(...) {
      calls <<- calls + 1L
      sentinel
    },
    .log_lik_cluster_norm_analytic_sum = function(...) {
      calls <<- calls + 1L
      sentinel
    },
    .package = "RoBMA"
  )
  setup <- function(is_weightfunction) {
    list(
      outcome_type      = "norm",
      is_weightfunction = is_weightfunction,
      effect_direction  = "positive",
      yi                = 0.1,
      sei               = 0.2,
      weights           = 1,
      data              = list()
    )
  }

  # The same weighted setup without selection reaches the closed-form kernel ...
  expect_identical(.log_lik_cluster_from_setup(setup(FALSE)), sentinel)
  expect_identical(.log_lik_cluster_sum_from_setup(setup(FALSE)), sentinel)
  expect_identical(calls, 2L)

  # ... and with selection it stops before either kernel.
  expect_error(.log_lik_cluster_from_setup(setup(TRUE)))
  expect_error(.log_lik_cluster_sum_from_setup(setup(TRUE)))
  expect_identical(calls, 2L)
})


test_that("selection deletion rejects a nonfinite full event density", {

  local_mocked_bindings(
    .selection_retains_sampling = function(...) FALSE,
    .data_selection_execution_plan = function(...) list(row_blocks = list(1:2)),
    .estimate_normal_covariance_target_location_from_setup = function(...) {
      list(means = matrix(0, 1L, 2L), y = c(0, 0))
    },
    .selection_joint_signed_context = function(...) list(),
    .selection_joint_random_factor_samples = function(...) list(),
    .selection_joint_covariance_lower = function(...) matrix(c(1, 0, 1), 1L),
    .selection_joint_event_numerator = function(...) list(log_numerator = -Inf),
    .package = "RoBMA"
  )
  local_mocked_bindings(selection_context_subset_observations = function(...) list(),
                        .package = "BayesTools")
  expect_error(.selection_joint_deletion_loglik_from_setup(
    list(S = 1L, data = list(), selection_sei = c(1, 1)), list(1L)),
    "observed outcomes have zero or invalid event density.", fixed = TRUE)
})

test_that("deletion fallback preserves coordinates resolved by the fast path", {

  local_mocked_bindings(
    .data_selection_execution_plan = function(...) list(row_blocks = list(1L, 2L),
      relative_tolerance = 1e-6, max_points_per_scramble = 128L),
    .data_selection_model = function(...) list(groups = list(row_blocks = list(1:2))),
    .selection_joint_signed_context = function(...) list(),
    .selection_conditioned_sampling_independent_targets = function(...) {
      list(mean = matrix(c(5, NA_real_), 1L))
    },
    .selection_conditioned_sampling_state = function(...) list(
      baseline_mu = matrix(0, 1L, 2L), e = matrix(0, 1L, 2L),
      sampling_covariance = diag(2), integrated_covariance = array(0, c(1L, 2L, 2L)),
      total_covariance = array(diag(c(0, 1)), c(1L, 2L, 2L))),
    .selection_deleted_gaussian_coordinates = function(...) {
      list(mean = c(0, 0), variance = c(0, 1))
    }, .package = "RoBMA"
  )
  out <- .selection_conditioned_sampling_estimate_targets(
    list(S = 1L, K = 2L, yi = c(0, 0), mu = matrix(c(2, 3), 1L)), "mean")
  expect_identical(out$mean, matrix(c(5, 3), 1L))
})

test_that("the cluster integral forms no quotient by a loading", {

  # References: 400-digit mpmath evaluation of
  #   -sum(w (log(2 pi) + 2 log(s))) / 2 - log(A) / 2 - (C - B^2 / A) / 2
  # with b = loading / s, z = residual / s, A = 1 + sum(w b^2), B = sum(w b z)
  # and C = sum(w z^2). The residual divided by the loading overflows in every
  # case below; the density does not, except in the last.
  cases <- list(
    # a standardized residual of 1e108 against a loading of +-1e-100
    list(yi = 1e308, sei = 1e200, tau_between = 1e-100,
      expected = -5.0000000000000004e+215),
    list(yi = 1e308, sei = 1e200, tau_between = -1e-100,
      expected = -5.0000000000000004e+215),
    # a latent value pinned by a loading of 1e100 next to a loading of 1e-300
    list(yi = c(.3, .2), sei = c(1, 1), tau_between = c(1e-300, 1e100),
      expected = -232.1413863658139),
    list(yi = c(.2, .3), sei = c(1, 1), tau_between = c(1e100, 1e-300),
      expected = -232.1413863658139),
    # exact value -5.0000000000000001e+615: not representable
    list(yi = 1e308, sei = 1, tau_between = 1e-100, expected = -Inf),
    list(yi = 1e308, sei = 1, tau_between = -1e-100, expected = -Inf),
    list(yi = 1e308, sei = 1, tau_between = 1e-50, expected = -Inf)
  )
  for (i in seq_along(cases)) {
    case  <- cases[[i]]
    info  <- paste("case", i)
    setup <- .cluster_norm_dispatch_setup(case$yi, case$sei, 0, case$tau_between)
    expected <- matrix(case$expected, 1L, 1L)
    agree <- function(observed, label) {
      if (is.finite(case$expected)) {
        expect_equal(observed, case$expected, tolerance = 1e-14,
          info = paste(info, label))
      } else {
        expect_identical(observed, case$expected, info = paste(info, label))
      }
    }
    agree(as.numeric(.log_lik_cluster_from_setup(setup)), "dispatch")
    agree(.log_lik_cluster_sum_from_setup(setup), "sum")
    agree(as.numeric(.log_lik_cluster_norm_analytic_r(setup, setup$yi, setup$sei)),
      "R helper")
    agree(.log_lik_cluster_norm_stable(case$yi, case$sei, case$tau_between), "stable")
    if (.has_native_norm_cluster_analytic()) {
      agree(as.numeric(.log_lik_cluster_norm_analytic_native(
        setup, setup$yi, setup$sei)), "native")
      agree(as.numeric(.log_lik_cluster_norm_analytic_native(
        setup, setup$yi, setup$sei, row_sum = TRUE)), "native row sum")
    }
  }

  # An unrepresentable cluster spoils only itself.
  setup <- .cluster_norm_dispatch_setup(c(1e308, .3, .2), rep(1, 3L), 0,
    c(1e-100, 1e-300, 1e100))
  setup$cluster <- list(1L, 2:3)
  expected <- c(-Inf, -232.1413863658139)
  expect_equal(as.numeric(.log_lik_cluster_from_setup(setup)), expected,
    tolerance = 1e-14)
  expect_equal(
    .log_lik_cluster_norm_stable(
      c(1e308, .3, .2), rep(1, 3L), c(1e-100, 1e-300, 1e100),
      group = c(1L, 2L, 2L)
    ),
    expected, tolerance = 1e-14
  )

  # Invalid inputs are still NA, not -Inf.
  expect_identical(.log_lik_cluster_norm_stable(NA_real_, 1, .5), NA_real_)
  expect_identical(.log_lik_cluster_norm_stable(.1, 0, .5), NA_real_)
})

test_that("the cluster integral keeps aligned residuals, standardized overflow and extreme powers", {

  # References: 1500-digit mpmath evaluation of
  #   -sum(w (log(2 pi) + 2 log(s))) / 2 - log(A) / 2 - (C - B^2 / A) / 2
  # with b = loading / s, z = residual / s, A = 1 + sum(w b^2), B = sum(w b z)
  # and C = sum(w z^2), on the exact double inputs. The two-regime form keeps
  # the residual after the latent projection exactly zero when the residual lies
  # along the loadings, never forms residual / sd for a loaded row whose
  # latent absorbs it, and never squares a scale before dividing by another.
  cases <- list(
    # residuals along the loadings at sd 1e-150 (the latent absorbs them exactly)
    list(yi = c(1, 3, 5), sei = rep(1e-150, 3L), tau_between = c(1, 3, 5),
      expected = 685.74103826785498),
    list(yi = c(2, 6, 10), sei = rep(1e-100, 3L), tau_between = c(1, 3, 5),
      weights = c(1, 2, 3), expected = 1141.5072679066598),
    list(yi = c(3, -9, 6), sei = rep(1e-120, 3L), tau_between = c(1, -3, 2),
      weights = c(1, 1, 2), expected = 819.30969346608967),
    # residual / sd overflows although the density is representable
    list(yi = 1e200, sei = 1e-150, tau_between = 1e200,
      expected = -461.93595713201381),
    # the information underflows
    list(yi = c(.5, .2), sei = c(1e200, 1e-10), tau_between = c(1e200, 0),
      expected = -2.0000000000000005e+18),
    list(yi = c(.5, .2, .1), sei = c(1e250, 1e-100, 1),
      tau_between = c(1e250, 0, 1e-120), weights = c(1, 2, 1),
      expected = -4.0000000000000003e+198),
    # a dominant row with a power that is not a power of two, and a residual
    # of 1e10 against an sd of 1e-300
    list(yi = 1, sei = 1e-120, tau_between = 1, weights = .1,
      expected = -248.11979135018038),
    list(yi = 1, sei = 1e-120, tau_between = 1, weights = 1.5,
      expected = 136.07396522578165),
    list(yi = 1, sei = 3e-110, tau_between = 1, weights = .3,
      expected = -176.70371871627228),
    list(yi = 1e10, sei = 1e-300, tau_between = 1, expected = -5e+19),
    # numerator terms that underflow next to a representable information
    list(yi = c(5e-60, -8e-201), sei = c(9e-60, 9e-200),
      tau_between = c(5e-56, 0), weights = c(3, 1),
      expected = 853.34187241648172),
    list(yi = c(-2e125, -3e-193, 2e-118), sei = c(8e123, 3e-193, 7e-118),
      tau_between = c(1.2e-17, -5e-162, 0), weights = c(3, 3, .1),
      expected = -514.56525521598696),
    # a loading 1e-395 times the largest one in the cluster
    list(yi = c(1.2e226, 1.1e20, -1.1e-247), sei = c(5.6e225, 4.5e19, 1.4e-249),
      tau_between = c(-4.4e230, 0, 8.3e-165), weights = c(7.3, 7.3, .5),
      expected = -4085.7133308171571),
    # powers near the largest double, with a loaded and an unloaded row
    list(yi = 0, sei = 1.1, tau_between = .5, weights = 1e308,
      expected = -1.0142487130089977e+308),
    list(yi = 0, sei = .1, tau_between = .5, weights = 1e308,
      expected = 1.3836465597893729e+308),
    list(yi = 0, sei = 1.1, tau_between = 0, weights = 1e308,
      expected = -1.0142487130089977e+308)
  )
  for (i in seq_along(cases)) {
    case  <- cases[[i]]
    info  <- paste("case", i)
    setup <- .cluster_norm_dispatch_setup(case$yi, case$sei, 0, case$tau_between,
      weights = case$weights)
    agree <- function(observed, label) {
      expect_equal(as.numeric(observed), case$expected, tolerance = 1e-14,
        info = paste(info, label))
    }
    agree(.log_lik_cluster_from_setup(setup), "dispatch")
    agree(.log_lik_cluster_sum_from_setup(setup), "sum")
    agree(.log_lik_cluster_norm_analytic_r(setup, setup$yi, setup$sei), "R helper")
    agree(.log_lik_cluster_norm_stable(case$yi, case$sei, case$tau_between,
      weights = case$weights), "stable")
    if (.has_native_norm_cluster_analytic()) {
      agree(.log_lik_cluster_norm_analytic_native(setup, setup$yi, setup$sei),
        "native")
      agree(.log_lik_cluster_norm_analytic_native(setup, setup$yi, setup$sei,
        row_sum = TRUE), "native row sum")
    }
  }

  # Clusters of different regimes in one call of the helper give what each
  # gives alone: loadings below the smallest sd, above it, none, and a
  # denominator that needs a log-scale factor.
  residual <- c(.3, -.2, .1, 1, 2, .5, .4, .5, .2, .1)
  sd       <- c(.5, .4, .3, 1e-3, 2e-3, .2, 1e200, 1e-10, .2, 1e250)
  loading  <- c(.1, -.2, 0, 1, -3, 0, 1e200, 0, 0, 1e250)
  weights  <- c(1, 2, .5, 1, 1, 3, 1, 1, 0, 1)
  group    <- c(1L, 1L, 1L, 2L, 2L, 3L, 4L, 4L, 5L, 6L)
  together <- .log_lik_cluster_norm_stable(residual, sd, loading, weights, group)
  alone    <- vapply(1:6, function(g) {
    use <- group == g
    .log_lik_cluster_norm_stable(residual[use], sd[use], loading[use],
      weights[use])
  }, numeric(1))
  expect_identical(together, alone)
  expect_identical(together[5L], 0)
})
