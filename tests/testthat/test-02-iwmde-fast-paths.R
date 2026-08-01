# ============================================================================ #
# test-02-iwmde-fast-paths.R
# ============================================================================ #

context("IWMDE fast-path equivalence")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


test_that("IWMDE predictor cache keys remain compact", {

  key <- .iwmde_predictor_cache_key(
    prefix     = "test",
    active_key = paste(rep("branch", 100), collapse = "|"),
    rows       = seq_len(5000),
    unit       = "estimate"
  )
  cache <- new.env(parent = emptyenv())

  expect_lt(nchar(key), 64)
  assign(key, TRUE, envir = cache)
  expect_true(get(key, envir = cache, inherits = FALSE))
})


test_that("IWMDE numeric cache keys preserve represented coordinates", {

  values <- c(
    0,
    .Machine$double.xmin,
    1e-20,
    2e-20,
    1,
    1 + .Machine$double.eps
  )
  keys <- .iwmde_key_number(values)

  expect_length(unique(keys), length(values))
  expect_identical(.iwmde_key_number(-0), .iwmde_key_number(0))
  expect_identical(
    .iwmde_key_number(c(0, 1, -1)),
    c("0000000000000000", "3ff0000000000000", "bff0000000000000")
  )
})


test_that("IWMDE identifies sampled random SD focal parameters", {

  dat <- data.frame(
    yi    = c(.10, .20, .30),
    x     = c(0, 1, 2),
    study = c("s1", "s2", "s3")
  )
  object <- brma.mv(
    yi                         = yi,
    V                          = diag(rep(.04, 3L)),
    random                     = ~ 1 | study,
    scale                      = ~ x,
    data                       = dat,
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    marginalize_estimate_level = FALSE,
    only_priors                = TRUE
  )
  prior <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  context <- list(
    object          = object,
    data            = object[["data"]],
    flat_prior_list = list(
      tau               = prior,
      log_tau_intercept = prior,
      log_tau_x         = prior,
      mu                = prior
    )
  )

  expect_true(.iwmde_parameter_controls_sampled_random_sd(
    context,
    "log_tau_intercept"
  ))
  expect_true(.iwmde_parameter_controls_sampled_random_sd(context, "log_tau_x"))
  expect_true(.iwmde_parameter_controls_sampled_random_sd(context, "tau[1]"))
  expect_false(.iwmde_parameter_controls_sampled_random_sd(context, "mu"))
  expect_null(.iwmde_predictor_column_basis(
    context = context,
    column  = "tau",
    state   = list(active_setup = list())
  ))
})


test_that("IWMDE disables focal prior delta for sampled random SD rows", {

  .skip_if_missing_raw_fits("brma.mv_block_mvn_random_scale")

  context   <- .iwmde_context(load_fit("brma.mv_block_mvn_random_scale", validate = FALSE))
  parameter <- "log_tau_intercept"
  if (!parameter %in% colnames(context[["posterior_samples"]])) {
    skip("brma.mv random-scale fixture does not contain log_tau_intercept.")
  }
  rows <- which(is.finite(context[["posterior_samples"]][, parameter]))
  rows <- head(rows, 3L)

  expect_gt(length(rows), 0L)
  for (row in rows) {
    state <- .iwmde_row_state(context, row, parameter)
    expect_false(state[["use_focal_prior_delta"]])
  }
})


test_that("IWMDE density aggregation matches row-wise reference", {

  log_terms <- matrix(c(
    0, -1, -Inf, NA,
    -Inf, -Inf, NA, NA,
    2, 1, 0, -1
  ), nrow = 3, byrow = TRUE)
  active_mass <- .75
  denominator <- ncol(log_terms)

  fast <- .iwmde_density_aggregate(
    log_terms   = log_terms,
    active_mass = active_mass,
    denominator = denominator
  )

  y                <- numeric(nrow(log_terms))
  finite_terms     <- integer(nrow(log_terms))
  max_log_ratio    <- rep(Inf, nrow(log_terms))
  ess              <- numeric(nrow(log_terms))
  max_weight_share <- rep(1, nrow(log_terms))
  contributions    <- matrix(0, nrow = nrow(log_terms), ncol = ncol(log_terms))

  for (g in seq_len(nrow(log_terms))) {
    finite <- is.finite(log_terms[g, ])
    finite_terms[g] <- sum(finite)
    if (any(finite)) {
      max_term            <- max(log_terms[g, finite])
      scaled_terms        <- exp(log_terms[g, finite] - max_term)
      sum_scaled_terms    <- sum(scaled_terms)
      y[g]                <- active_mass * exp(max_term) *
        sum_scaled_terms / denominator
      contributions[g, finite] <- active_mass * exp(log_terms[g, finite])
      max_log_ratio[g]    <- max_term - stats::median(log_terms[g, finite])
      ess[g]              <- sum_scaled_terms^2 / sum(scaled_terms^2)
      max_weight_share[g] <- max(scaled_terms) / sum_scaled_terms
    }
  }

  expect_equal(fast[["y"]], y)
  expect_equal(fast[["finite_terms"]], finite_terms)
  expect_equal(fast[["max_log_ratio"]], max_log_ratio)
  expect_equal(fast[["ess"]], ess)
  expect_equal(fast[["max_weight_share"]], max_weight_share)
  expect_equal(fast[["contributions"]], contributions)

  padded <- .iwmde_density_aggregate(
    log_terms   = matrix(log(2), nrow = 1L),
    active_mass = 1,
    denominator = 2L
  )
  expect_equal(padded[["y"]], 1)
  expect_equal(rowMeans(padded[["contributions"]]), padded[["y"]])
})


test_that("qCMDE density normalization is row-wise and mass preserving", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(list(), list())

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      cbind(
        rep(log(2), length(values)),
        rep(log(8), length(values))
      )
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = .6,
    replacement        = list(type = "scalar")
  )

  expect_equal(density[["estimator"]], "q_grid_cmde")
  expect_equal(density[["normalization_mass_ratio"]], 1)
  final_width <- diff(density[["normalization_range"]])
  pilot_width <- diff(density[["normalization_initial_range"]])
  expect_equal(density[["pilot_normalization_integral"]],
               .6 * pilot_width / final_width,
               tolerance = 1e-10)
  expect_equal(density[["final_normalization_integral"]], .6,
               tolerance = 1e-10)
  expect_equal(density[["normalization_relative_error"]], 0,
               tolerance = 1e-10)
  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]), .6 / final_width,
               tolerance = 1e-10)
  expect_equal(density[["y"]], rep(.6 / final_width, length(grid)),
               tolerance = 1e-10)
  expect_equal(density[["pilot_y"]], rep(.6 / pilot_width, length(grid)),
               tolerance = 1e-10)
  expect_equal(density[["n_normalized_rows"]], 2L)
})


test_that("qCMDE density matches an analytic row-normalized normal mixture", {

  display_grid <- seq(-3, 3, length.out = 101)
  norm_grid    <- seq(-6, 6, length.out = 601)
  means        <- c(-.75, 1.25)
  sds          <- c(.5, 1.1)
  active_mass  <- .8
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      cbind(
        stats::dnorm(values, means[[1L]], sds[[1L]], log = TRUE),
        stats::dnorm(values, means[[2L]], sds[[2L]], log = TRUE)
      )
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = list(list(), list()),
    active_mass        = active_mass,
    replacement        = list(type = "scalar")
  )
  normalizer <- stats::pnorm(
    density[["normalization_range"]][[2L]],
    means,
    sds
  ) - stats::pnorm(
    density[["normalization_range"]][[1L]],
    means,
    sds
  )
  expected <- active_mass / 2 * (
    stats::dnorm(display_grid, means[[1L]], sds[[1L]]) / normalizer[[1L]] +
      stats::dnorm(display_grid, means[[2L]], sds[[2L]]) / normalizer[[2L]]
  )

  expect_equal(density[["y"]], expected, tolerance = 2e-4)
  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]),
               .iwmde_trapz(display_grid, expected),
               tolerance = 1e-6)
})


test_that("qCMDE and IWMDE match conjugate normal-normal posterior oracle", {

  yi <- c(-.25, .10, .35, .60)
  sei <- c(.30, .45, .40, .50)
  prior_mean <- .15
  prior_sd   <- 1.20

  posterior_var <- 1 / (1 / prior_sd^2 + sum(1 / sei^2))
  posterior_sd  <- sqrt(posterior_var)
  posterior_mean <- posterior_var * (
    prior_mean / prior_sd^2 + sum(yi / sei^2)
  )
  log_q <- function(mu) {
    vapply(mu, function(value) {
      sum(stats::dnorm(yi, mean = value, sd = sei, log = TRUE)) +
        stats::dnorm(value, mean = prior_mean, sd = prior_sd, log = TRUE)
    }, numeric(1))
  }

  display_grid <- seq(
    posterior_mean - 3 * posterior_sd,
    posterior_mean + 3 * posterior_sd,
    length.out = 91
  )
  norm_grid <- seq(
    posterior_mean - 8 * posterior_sd,
    posterior_mean + 8 * posterior_sd,
    length.out = 501
  )
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )
  active_values <- stats::qnorm(
    stats::ppoints(80),
    mean = posterior_mean,
    sd   = posterior_sd
  )
  row_states <- lapply(active_values, function(value) {
    list(baseline_log_q = log_q(value))
  })

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        log_q(values),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .iwmde_chen_log_weight = function(context, parameter, parameter_spec,
                                      active_rows, active_values,
                                      weight_rows, weight_values, support) {

      list(
        log_weight = stats::dnorm(
          active_values,
          mean = posterior_mean,
          sd   = posterior_sd,
          log  = TRUE
        ),
        method = "oracle_posterior"
      )
    },
    .package = "RoBMA"
  )

  qcmde <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )
  qcmde_mass <- diff(stats::pnorm(
    qcmde[["normalization_range"]],
    mean = posterior_mean,
    sd   = posterior_sd
  ))
  expect_equal(
    qcmde[["y"]],
    stats::dnorm(display_grid, posterior_mean, posterior_sd) / qcmde_mass,
    tolerance = 1e-3
  )

  iwmde <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu",
    parameter_spec     = list(type = "primitive"),
    display_grid       = display_grid,
    row_states         = row_states,
    active_rows        = seq_along(active_values),
    active_values      = active_values,
    weight_rows        = seq_along(active_values),
    weight_values      = active_values,
    support            = c(-Inf, Inf),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_grid = normalization_grid
  )
  expect_equal(
    iwmde[["y"]],
    stats::dnorm(display_grid, posterior_mean, posterior_sd),
    tolerance = 1e-10
  )
  expect_equal(iwmde[["weight_method"]], "oracle_posterior")
})


test_that("qCMDE and IWMDE match correlated linear-contrast normal oracle", {

  mean_beta <- c(mu_intercept = .20, mu_ablat = -.35)
  sigma <- matrix(
    c(.50^2, .50 * .30 * .65,
      .50 * .30 * .65, .30^2),
    nrow = 2,
    dimnames = list(names(mean_beta), names(mean_beta))
  )
  weights <- c(mu_intercept = 1, mu_ablat = -.5)
  target_mean <- sum(weights * mean_beta[names(weights)])
  target_sd   <- sqrt(as.numeric(t(weights) %*% sigma[names(weights),
                                                        names(weights)] %*%
                                   weights))
  log_q <- function(value) {
    stats::dnorm(value, mean = target_mean, sd = target_sd, log = TRUE)
  }

  display_grid <- seq(target_mean - 3 * target_sd,
                      target_mean + 3 * target_sd,
                      length.out = 81)
  norm_grid <- seq(target_mean - 8 * target_sd,
                   target_mean + 8 * target_sd,
                   length.out = 501)
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )
  active_values <- stats::qnorm(
    stats::ppoints(80),
    mean = target_mean,
    sd   = target_sd
  )
  row_states <- lapply(active_values, function(value) {
    list(baseline_log_q = log_q(value))
  })
  spec <- list(type = "linear", weights = weights)

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        log_q(values),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .iwmde_chen_log_weight = function(context, parameter, parameter_spec,
                                      active_rows, active_values,
                                      weight_rows, weight_values, support) {

      list(
        log_weight = stats::dnorm(
          active_values,
          mean = target_mean,
          sd   = target_sd,
          log  = TRUE
        ),
        method = "oracle_linear_contrast"
      )
    },
    .package = "RoBMA"
  )

  qcmde <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu_intercept - .5 * mu_ablat",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = 1,
    replacement        = spec
  )
  qcmde_mass <- diff(stats::pnorm(
    qcmde[["normalization_range"]],
    mean = target_mean,
    sd   = target_sd
  ))
  expect_equal(
    qcmde[["y"]],
    stats::dnorm(display_grid, target_mean, target_sd) / qcmde_mass,
    tolerance = 1e-3
  )

  iwmde <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu_intercept - .5 * mu_ablat",
    parameter_spec     = spec,
    display_grid       = display_grid,
    row_states         = row_states,
    active_rows        = seq_along(active_values),
    active_values      = active_values,
    weight_rows        = seq_along(active_values),
    weight_values      = active_values,
    support            = c(-Inf, Inf),
    active_mass        = 1,
    replacement        = spec,
    normalization_grid = normalization_grid
  )
  expect_equal(
    iwmde[["y"]],
    stats::dnorm(display_grid, target_mean, target_sd),
    tolerance = 1e-10
  )
  expect_equal(iwmde[["weight_method"]], "oracle_linear_contrast")
})


test_that("qCMDE point ordinates are independent of unrelated display values", {

  norm_grid <- seq(-6, 6, length.out = 301)
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        stats::dnorm(values, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  base <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = 0,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 30),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )
  with_far_value <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = c(0, 100),
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 30),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )

  expect_equal(with_far_value[["y"]][[1L]], base[["y"]][[1L]],
               tolerance = 1e-12)
  expect_equal(with_far_value[["normalization_range"]],
               base[["normalization_range"]])
})


test_that("qCMDE point ordinates are stable under doubled integration budget", {

  normalizer_grid <- function(n) {
    z <- seq(-6, 6, length.out = n)
    list(
      x            = z,
      z            = z,
      log_jacobian = rep(0, length(z))
    )
  }

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        stats::dnorm(values, mean = .2, sd = 1.1, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  base <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = 0,
    normalization_grid = normalizer_grid(151),
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 40),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )
  doubled <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = 0,
    normalization_grid = normalizer_grid(301),
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 40),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )

  expect_equal(doubled[["y"]][[1L]], base[["y"]][[1L]], tolerance = 1e-4)
  expect_lt(doubled[["ordinate_relative_change"]][[1L]], .01)
})


test_that("qCMDE dropped normalizer rows keep target denominator", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(list(), list())

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      cbind(
        rep(log(2), length(values)),
        rep(-Inf, length(values))
      )
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = .6,
    replacement        = list(type = "scalar")
  )

  final_width <- diff(density[["normalization_range"]])
  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]), .3 / final_width,
               tolerance = 1e-10)
  expect_equal(density[["y"]], rep(.3 / final_width, length(grid)),
               tolerance = 1e-10)
  expect_equal(density[["n_candidate_rows"]], 2L)
  expect_equal(density[["n_normalized_rows"]], 1L)
  expect_equal(density[["n_dropped_normalizer"]], 1L)
  expect_equal(density[["row_drop_fraction"]], .5)
})


test_that("qCMDE uses refined normalizers and diagnoses pilot-grid impact", {

  display_grid <- 0
  norm_grid    <- seq(-.5, .5, length.out = 21)
  row_states   <- rep(list(list()), 20)
  normalization_grid <- list(
    x            = norm_grid,
    z            = norm_grid,
    log_jacobian = rep(0, length(norm_grid))
  )

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(
        stats::dnorm(values, log = TRUE),
        nrow = length(values),
        ncol = length(row_states)
      )
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = row_states,
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )

  expect_lt(density[["pilot_normalization_integral"]], 1)
  expect_equal(density[["final_normalization_integral"]], 1, tolerance = 1e-10)
  expect_gt(density[["ordinate_relative_change"]][[1]], .05)
  expect_gt(density[["max_normalizer_relative_change"]], .05)
  expect_match(
    .iwmde_diagnostics_bf_failure_reason(list(
      estimator              = "q_grid_cmde",
      ordinate               = density[["y"]][[1]],
      relative_mcse          = .1,
      finite_terms           = 20,
      ess                    = 20,
      max_weight_share       = .2,
      evaluation_value       = 0,
      normalization_range    = range(norm_grid),
      normalization_relative_error =
        density[["normalization_relative_error"]],
      ordinate_relative_change =
        density[["ordinate_relative_change"]][[1]],
      max_normalizer_relative_change =
        density[["max_normalizer_relative_change"]]
    )),
    "qCMDE.*ordinate"
  )
})


test_that("qCMDE BF gate rejects validation movement despite perfect mass", {

  diagnostics <- list(
    estimator              = "q_grid_cmde",
    ordinate               = 1,
    relative_mcse          = .01,
    finite_terms           = 80L,
    ess                    = 40,
    max_weight_share       = .10,
    evaluation_value       = 0,
    final_normalization_integral = 1,
    normalization_relative_error = 0,
    active_mass            = 1,
    normalization_mass_ratio = 1,
    row_drop_fraction      = 0,
    ordinate_relative_change = .06,
    max_normalizer_relative_change = 0
  )

  expect_match(
    .iwmde_diagnostics_bf_failure_reason(diagnostics),
    "qCMDE.*ordinate"
  )
})


test_that("IWMDE density reports raw support mass without scaling the curve", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(
    list(baseline_log_q = 0),
    list(baseline_log_q = 0)
  )

  testthat::local_mocked_bindings(
    .iwmde_chen_log_weight = function(context, parameter, parameter_spec,
                                      active_rows, active_values,
                                      weight_rows, weight_values, support) {

      list(log_weight = rep(log(.5), length(active_values)), method = "mock")
    },
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(0, nrow = length(values), ncol = length(row_states))
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu",
    parameter_spec     = list(type = "primitive"),
    display_grid       = grid,
    row_states         = row_states,
    active_rows        = 1:2,
    active_values      = c(.25, .75),
    weight_rows        = 1:2,
    weight_values      = c(.25, .75),
    support            = c(0, 1),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_grid = normalization_grid
  )

  expect_equal(density[["estimator"]], "iwmde")
  expect_equal(density[["weight_method"]], "mock")
  expect_equal(
    density[["support_grid_normalization_integral"]],
    .5,
    tolerance = 1e-10
  )
  expect_equal(density[["normalization_mass_ratio"]], 2, tolerance = 1e-10)
  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]), .5,
               tolerance = 1e-10)
  expect_equal(density[["y"]], rep(.5, length(grid)), tolerance = 1e-10)
})


test_that("IWMDE dropped weight rows keep target denominator", {

  grid <- seq(0, 1, length.out = 101)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(
    list(baseline_log_q = 0),
    list(baseline_log_q = 0)
  )

  testthat::local_mocked_bindings(
    .iwmde_chen_log_weight = function(context, parameter, parameter_spec,
                                      active_rows, active_values,
                                      weight_rows, weight_values, support) {

      list(log_weight = c(log(.5), -Inf), method = "mock")
    },
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(0, nrow = length(values), ncol = length(row_states))
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu",
    parameter_spec     = list(type = "primitive"),
    display_grid       = grid,
    row_states         = row_states,
    active_rows        = 1:2,
    active_values      = c(.25, .75),
    weight_rows        = 1:2,
    weight_values      = c(.25, .75),
    support            = c(0, 1),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_grid = normalization_grid
  )

  expect_equal(.iwmde_trapz(density[["x"]], density[["y"]]), .25,
               tolerance = 1e-10)
  expect_equal(density[["y"]], rep(.25, length(grid)), tolerance = 1e-10)
  expect_equal(
    density[["support_grid_normalization_integral"]],
    .25,
    tolerance = 1e-10
  )
  expect_equal(density[["normalization_mass_ratio"]], 4, tolerance = 1e-10)
  expect_equal(density[["n_candidate_rows"]], 2L)
  expect_equal(density[["n_normalized_rows"]], 1L)
  expect_equal(density[["n_dropped_weight"]], 1L)
  expect_equal(density[["row_drop_fraction"]], .5)
})


test_that("restricted Chen fallback weights are moment matched", {

  p <- stats::ppoints(500)

  distance_fit <- stats::qgamma(p, shape = 18, scale = .04)
  active       <- seq(.02, 2, length.out = 100)
  gamma_weight <- .iwmde_chen_gamma_log_weight_single(
    active_values = active,
    weight_values = distance_fit,
    support       = c(0, Inf)
  )

  shape <- mean(distance_fit)^2 / stats::var(distance_fit)
  rate  <- mean(distance_fit) / stats::var(distance_fit)

  expect_equal(gamma_weight[["method"]], "chen_gamma")
  expect_equal(
    gamma_weight[["log_weight"]],
    stats::dgamma(active, shape = shape, rate = rate, log = TRUE),
    tolerance = 1e-12
  )

  prob_fit <- stats::qbeta(p, shape1 = 6, shape2 = 3)
  active   <- seq(.01, .99, length.out = 100)
  beta_weight <- .iwmde_chen_beta_log_weight(
    active_values = active,
    weight_values = prob_fit,
    support       = c(0, 1)
  )

  location <- mean(prob_fit)
  variance <- stats::var(prob_fit)
  common   <- location * (1 - location) / variance - 1

  expect_equal(beta_weight[["method"]], "chen_beta")
  expect_equal(
    beta_weight[["log_weight"]],
    stats::dbeta(
      active,
      shape1 = location * common,
      shape2 = (1 - location) * common,
      log    = TRUE
    ),
    tolerance = 1e-12
  )
})


test_that("Chen weights dispatch by row-specific support", {

  supports <- matrix(
    c(0, 1,
      0, 1,
      1, 2,
      1, 2),
    ncol  = 2,
    byrow = TRUE
  )
  rownames(supports) <- as.character(seq_len(4L))
  calls <- list()

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(
        context, parameter, parameter_spec, active_rows, active_values,
        weight_rows, weight_values, support) {

      calls[[length(calls) + 1L]] <<- list(
        active_rows = active_rows,
        weight_rows = weight_rows,
        support     = support
      )

      list(
        log_weight = rep(sum(support), length(active_values)),
        method     = paste0("bounded_", paste(support, collapse = "_"))
      )
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(supports = supports),
    parameter      = "theta",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(4L),
    active_values  = c(.25, .75, 1.25, 1.75),
    weight_rows    = seq_len(4L),
    weight_values  = c(.30, .70, 1.30, 1.70),
    support        = c(0, 2)
  )

  expect_length(calls, 2L)
  expect_equal(calls[[1L]][["active_rows"]], 1:2)
  expect_equal(calls[[1L]][["weight_rows"]], 1:2)
  expect_equal(calls[[1L]][["support"]], c(0, 1))
  expect_equal(calls[[2L]][["active_rows"]], 3:4)
  expect_equal(calls[[2L]][["weight_rows"]], 3:4)
  expect_equal(calls[[2L]][["support"]], c(1, 2))
  expect_equal(weight[["log_weight"]], c(1, 1, 3, 3))
  expect_equal(weight[["method"]], "chen_mixed(bounded_0_1,bounded_1_2)")
  expect_length(weight[["partitions"]], 2L)
  expect_equal(
    vapply(weight[["partitions"]], `[[`, character(1), "method"),
    c("bounded_0_1", "bounded_1_2")
  )
  expect_equal(
    vapply(weight[["partitions"]], `[[`, integer(1), "n_eval_rows"),
    c(2L, 2L)
  )
})


test_that("Chen weights dispatch by bias branch before omega values", {

  supports <- matrix(
    rep(c(0, 1), 4L),
    ncol  = 2,
    byrow = TRUE
  )
  rownames(supports) <- as.character(seq_len(4L))
  samples <- data.frame(
    bias_indicator = c(1, 1, 2, 2),
    `omega[1]`     = c(.4, .4, .4, .4),
    check.names    = FALSE
  )
  calls   <- list()

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(
        context, parameter, parameter_spec, active_rows, active_values,
        weight_rows, weight_values, support) {

      calls[[length(calls) + 1L]] <<- list(
        active_rows = active_rows,
        weight_rows = weight_rows,
        support     = support
      )

      list(
        log_weight = rep(active_rows[[1L]], length(active_values)),
        method     = paste0("branch_", active_rows[[1L]])
      )
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(
      supports          = supports,
      posterior_samples = samples,
      indicator_names   = "bias_indicator",
      selection_spec    = list(jags_omega = "omega")
    ),
    parameter      = "theta",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(4L),
    active_values  = c(.25, .75, .25, .75),
    weight_rows    = seq_len(4L),
    weight_values  = c(.30, .70, .30, .70),
    support        = c(0, 1)
  )

  expect_length(calls, 2L)
  expect_equal(calls[[1L]][["active_rows"]], 1:2)
  expect_equal(calls[[1L]][["weight_rows"]], 1:2)
  expect_equal(calls[[2L]][["active_rows"]], 3:4)
  expect_equal(calls[[2L]][["weight_rows"]], 3:4)
  expect_equal(weight[["log_weight"]], c(1, 1, 3, 3))
  expect_equal(weight[["method"]], "chen_mixed(branch_1,branch_3)")
  expect_length(weight[["partitions"]], 2L)
  expect_equal(
    vapply(weight[["partitions"]], `[[`, character(1), "method"),
    c("branch_1", "branch_3")
  )
})


test_that("Chen weight dispatcher falls back when conditioning fails", {

  supports <- matrix(
    rep(c(-Inf, Inf), 3L),
    ncol  = 2L,
    byrow = TRUE
  )
  rownames(supports) <- as.character(seq_len(3L))

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_conditional_normal_log_weight = function(...) {

      .iwmde_chen_conditional_stop("conditioning failed")
    },
    .iwmde_chen_marginal_normal_log_weight = function(active_values,
                                                      weight_values) {

      list(log_weight = rep(log(.25), length(active_values)),
           method = "chen_marginal_normal")
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(supports = supports),
    parameter      = "mu",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(3L),
    active_values  = c(-1, 0, 1),
    weight_rows    = seq_len(3L),
    weight_values  = c(-1, 0, 1),
    support        = c(-Inf, Inf)
  )

  expect_equal(weight[["method"]], "chen_marginal_normal")
  expect_equal(weight[["log_weight"]], rep(log(.25), 3L))

  supports[,] <- rep(c(0, 1), each = 3L)
  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(...) {

      .iwmde_chen_conditional_stop("conditioning failed")
    },
    .iwmde_chen_beta_log_weight = function(active_values, weight_values,
                                           support) {

      list(log_weight = rep(log(.5), length(active_values)),
           method = "chen_beta")
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(supports = supports),
    parameter      = "rho",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(3L),
    active_values  = c(.25, .5, .75),
    weight_rows    = seq_len(3L),
    weight_values  = c(.25, .5, .75),
    support        = c(0, 1)
  )

  expect_equal(weight[["method"]], "chen_beta")
  expect_equal(weight[["log_weight"]], rep(log(.5), 3L))
})


test_that("Chen conditional-normal weights match bivariate normal oracle", {

  n <- 500L
  x <- seq(-2, 2, length.out = n)
  z <- .3 + 1.2 * x + stats::qnorm(stats::ppoints(n)) * .4
  active_values <- z + .15 * sin(seq_len(n) / 17)
  samples <- cbind(mu = z, PET = x)
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = NULL
  )

  weight <- .iwmde_chen_conditional_normal_log_weight(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(n),
    active_values  = active_values,
    weight_rows    = seq_len(n),
    weight_values  = z
  )

  x_center <- mean(x)
  x_scale  <- stats::sd(x)
  x_fit    <- (x - x_center) / x_scale
  cov_mat  <- stats::cov(cbind(z, x_fit))
  beta     <- cov_mat[1L, 2L] / cov_mat[2L, 2L]
  means    <- mean(z) + ((x - x_center) / x_scale) * beta
  variance <- sum((z - means)^2) / (length(z) - 1L)
  expected <- stats::dnorm(
    active_values,
    mean = means,
    sd   = sqrt(variance),
    log  = TRUE
  )

  expect_equal(weight[["method"]], "chen_conditional_normal")
  expect_equal(weight[["log_weight"]], expected, tolerance = 1e-12)
})


test_that("Chen Gaussian weights preserve represented small scales", {

  x      <- 2^-400 * seq(-2, 2, length.out = 200L)
  z      <- 2^-300 * (seq(-1, 1, length.out = 200L) +
    stats::qnorm(stats::ppoints(200L)) / 10)
  active <- z + 2^-305 * sin(seq_along(z))

  gaussian <- .iwmde_chen_conditional_gaussian(
    z_fit  = z,
    x_fit  = matrix(x),
    x_eval = matrix(x)
  )
  expected <- stats::dnorm(
    active / gaussian[["focal_scale"]],
    mean = gaussian[["means"]],
    sd   = gaussian[["sd"]],
    log  = TRUE
  ) - log(gaussian[["focal_scale"]])

  expect_true(all(is.finite(expected)))
  expect_true(gaussian[["sd"]] > 0)
})


test_that("Chen conditioning transforms do not clip support violations", {

  unit <- .iwmde_chen_transform_unit_interval(
    c(0, 1, -.Machine$double.eps, 1 + .Machine$double.eps),
    c(0, 1)
  )
  nonnegative <- .iwmde_chen_transform_nonnegative(
    c(-.Machine$double.xmin, 0, .Machine$double.xmin),
    0
  )

  expect_identical(unit[["fit"]], c(-Inf, Inf, NA_real_, NA_real_))
  expect_identical(
    nonnegative[["fit"]],
    c(NA_real_, 0, log1p(.Machine$double.xmin))
  )
})


test_that("logit conditional Chen weights are proper on original scale", {

  n       <- 500L
  x       <- seq(-2, 2, length.out = n)
  z       <- -.2 + .8 * x + stats::qnorm(stats::ppoints(n)) * .15
  rho     <- stats::plogis(z)
  samples <- cbind(rho = rho, mu = x)
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = NULL
  )

  grid <- seq(.001, .999, length.out = 2000)
  weight <- .iwmde_chen_logit_conditional_normal_log_weight(
    context        = context,
    parameter      = "rho",
    parameter_spec = list(type = "primitive"),
    active_rows    = rep(250L, length(grid)),
    active_values  = grid,
    weight_rows    = seq_len(n),
    weight_values  = rho,
    support        = c(0, 1)
  )

  expect_equal(weight[["method"]], "chen_logit_conditional_normal")
  expect_equal(
    .iwmde_trapz(grid, exp(weight[["log_weight"]])),
    1,
    tolerance = 5e-3
  )
})

test_that("IWMDE omega matrix collapse matches row-wise collapse", {

  omega       <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), nrow = 3, byrow = TRUE)
  global_cuts <- c(0, .025, .05, .50, 1)
  active_cuts <- c(0, .05, 1)

  fast <- .iwmde_collapse_omega_matrix(
    omega       = omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  )
  ref <- t(apply(omega, 1, .iwmde_collapse_omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  ))

  expect_equal(fast, ref)
  expect_equal(fast[, 1], c(1, 3, 5))
})


test_that("IWMDE omega collapse rejects unequal merged bins", {

  global_cuts <- c(0, .025, .05, .50, 1)
  active_cuts <- c(0, .05, .50, 1)
  omega       <- c(1, .5, .25, .25)

  collapsed <- .iwmde_collapse_omega(
    omega       = omega,
    global_cuts = global_cuts,
    active_cuts = active_cuts
  )

  expect_true(is.na(collapsed[1]))
  expect_equal(collapsed[2:3], c(.25, .25))
})


test_that("IWMDE structural identities are exact", {

  expect_error(
    .iwmde_indicator_index(1 + .Machine$double.eps, "model_indicator", 2L),
    "integer-valued"
  )
  expect_false(.iwmde_same_p_cuts(
    c(0, .5, 1),
    c(0, .5 + .Machine$double.eps, 1)
  ))

  collapsed <- .iwmde_collapse_omega(
    omega       = c(1, 1 + .Machine$double.eps),
    global_cuts = c(0, .25, .5),
    active_cuts = c(0, .5)
  )
  expect_true(is.na(collapsed[[1L]]))

  tiny <- .Machine$double.eps / 2
  expect_identical(
    .iwmde_linear_weights(c(mu = tiny, tau = 0)),
    c(mu = tiny)
  )

  points <- .iwmde_point_mass_table(
    c(1, 1 + .Machine$double.eps, 1),
    denominator = 3L
  )
  expect_identical(points[["x"]], c(1, 1 + .Machine$double.eps))
  expect_identical(points[["mass"]], c(2 / 3, 1 / 3))
})


test_that("IWMDE omega helpers use selection-specific JAGS names", {

  context <- list(
    selection_spec = list(
      n_bins     = 4L,
      p_cuts     = c(0, .025, .05, .50, 1),
      jags_omega = "custom.omega+beta"
    )
  )
  active_setup <- list(
    is_weightfunction = TRUE,
    selection_spec    = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega+beta"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  samples <- matrix(
    c(1, 1, .5, .5, 2, 10,
      2, 2, .4, .4, 2, 11),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      paste0("custom.omega+beta[", 1:4, "]"),
      "bias_indicator",
      "mu"
    ))
  )

  out <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  row_omega <- .iwmde_active_omega(
    context      = context,
    row          = samples[1, ],
    active_setup = list(
      selection_spec = context[["selection_spec"]],
      priors         = active_setup[["priors"]]
    )
  )

  expect_equal(samples[, "bias_indicator"], c(2, 2))
  expect_equal(out[, "bias_indicator"], c(1, 1))
  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(out))))
  expect_true(all(paste0("active.omega+beta[", 1:2, "]") %in% colnames(out)))
  expect_equal(out[, "active.omega+beta[1]"], samples[, "custom.omega+beta[1]"])
  expect_equal(out[, "active.omega+beta[2]"], samples[, "custom.omega+beta[3]"])
  expect_equal(out[, "mu"], samples[, "mu"])
  expect_equal(row_omega, as.numeric(samples[1, paste0("custom.omega+beta[", 1:4, "]")]))
})


test_that("IWMDE parameter filters use selection-specific JAGS names", {

  context <- list(
    selection_spec = list(jags_omega = "custom.omega+beta"),
    posterior_samples = matrix(
      0,
      nrow = 2,
      ncol = 4,
      dimnames = list(NULL, c(
        "mu",
        "custom.omega+beta[1]",
        "custom.omega+beta[2]",
        "eta[1]"
      ))
    ),
    indicator_names    = character(),
    flat_prior_list    = list(),
    focal_prior_cache  = new.env(parent = emptyenv()),
    support_cache      = new.env(parent = emptyenv())
  )

  expect_true(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "custom.omega+beta[1]",
    context   = context
  ))
  expect_true(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "eta[1]",
    context   = context
  ))
  expect_false(.iwmde_parameter_is_weightfunction_coordinate(
    parameter = "log_omega[1]",
    context   = context
  ))
  expect_true(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "custom.omega+beta[1]"
  ))
  expect_false(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "eta[1]"
  ))
  expect_false(.iwmde_chen_is_global_conditioning_column(
    context = context,
    column  = "log_omega[1]"
  ))
  expect_equal(
    .iwmde_parameter_spec(context, "custom.omega+beta[1]")[["status"]],
    "unsupported"
  )
})


test_that("IWMDE omega localization renames same-cut active branches", {

  context <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "custom.omega+beta"
    )
  )
  active_setup <- list(
    is_weightfunction = TRUE,
    selection_spec    = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega+beta"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  samples <- matrix(
    c(1, .5, 2, 10,
      2, .4, 2, 11),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      paste0("custom.omega+beta[", 1:2, "]"),
      "bias_indicator",
      "mu"
    ))
  )

  out <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  out_again <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = out,
    active_setup = active_setup
  )
  with_active <- cbind(
    samples,
    "active.omega+beta[1]" = c(.9, .8),
    "active.omega+beta[2]" = c(.7, .6)
  )
  active_first <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = with_active,
    active_setup = active_setup
  )

  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(out))))
  expect_true(all(paste0("active.omega+beta[", 1:2, "]") %in% colnames(out)))
  expect_equal(out[, "active.omega+beta[1]"], samples[, "custom.omega+beta[1]"])
  expect_equal(out[, "active.omega+beta[2]"], samples[, "custom.omega+beta[2]"])
  expect_equal(out[, "bias_indicator"], c(1, 1))
  expect_equal(out_again, out)
  expect_equal(active_first[, "active.omega+beta[1]"], c(.9, .8))
  expect_equal(active_first[, "active.omega+beta[2]"], c(.7, .6))
  expect_false(any(grepl("^custom\\.omega\\+beta\\[", colnames(active_first))))
})


test_that("IWMDE active omega does not truncate unexpected omega lengths", {

  context <- list(
    selection_spec = list(
      n_bins     = 4L,
      p_cuts     = c(0, .025, .05, .50, 1),
      jags_omega = "omega"
    )
  )
  active_setup <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "omega"
    ),
    priors = list(outcome = list(bias = NULL))
  )
  row <- c("omega[1]" = 1, "omega[2]" = .5, "omega[3]" = .25)

  omega <- .iwmde_active_omega(
    context      = context,
    row          = row,
    active_setup = active_setup
  )

  expect_true(all(is.na(omega)))
  expect_length(omega, 2L)
})

test_that("IWMDE active omega fails closed when weights are unavailable", {

  context <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "global.omega"
    )
  )
  active_setup <- list(
    selection_spec = list(
      n_bins     = 2L,
      p_cuts     = c(0, .05, 1),
      jags_omega = "active.omega"
    ),
    priors = list(outcome = list(bias = NULL))
  )

  omega <- .iwmde_active_omega(
    context      = context,
    row          = c(mu = 0),
    active_setup = active_setup
  )

  expect_equal(omega, rep(NA_real_, 2L))

  active_setup[["selection_spec"]][["jags_data"]] <- list(
    active.omega = c(1, .5)
  )
  omega <- .iwmde_active_omega(
    context      = context,
    row          = c(mu = 0),
    active_setup = active_setup
  )

  expect_equal(omega, c(1, .5))
})

test_that("IWMDE selection context uses localized active-branch indicators", {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .7, .35))
  )
  data <- list(
    outcome = data.frame(
      yi  = c(.1, .2),
      sei = c(.1, .2)
    )
  )
  attr(data, "effect_direction") <- "positive"

  context <- list(
    object = list(fit = list()),
    data   = data
  )
  active_setup <- list(
    priors            = list(outcome = list(bias = prior)),
    is_weightfunction = TRUE
  )
  samples <- matrix(
    c(2, 2),
    nrow     = 2,
    dimnames = list(NULL, "bias_indicator")
  )

  selection_context <- .iwmde_selection_context_active_branch(
    context           = context,
    active_setup      = active_setup,
    posterior_samples = samples
  )

  expect_equal(selection_context[["bias_indicator"]], c(1L, 1L))
  expect_false(any(selection_context[["use_normal"]]))
  expect_true(all(selection_context[["kernel_mode"]] != SELKERNEL_NORMAL))
  expect_equal(samples[, "bias_indicator"], c(2, 2))
})


test_that("negative-direction PET and PEESE location fast paths match flipped likelihoods", {

  yi  <- c(.20, -.10, .35)
  sei <- c(.25, .40, .55)
  data <- list(outcome = data.frame(yi = yi, sei = sei))
  attr(data, "outcome_type")      <- "norm"
  attr(data, "effect_direction")  <- "negative"

  mu <- matrix(
    c(.10, -.05, .20,
      -.15, .05, .25),
    nrow  = 2L,
    byrow = TRUE
  )
  tau <- matrix(
    c(.30, .35, .40,
      .32, .38, .44),
    nrow  = 2L,
    byrow = TRUE
  )
  current <- c(.15, -.20)
  values  <- c(-.30, .05, .40)
  active_setup <- list(is_weightfunction = FALSE)
  context <- list(data = data)
  setup <- list(
    yi               = yi,
    mu               = mu,
    tau_within       = tau,
    sei              = sei,
    weights          = NULL,
    effect_direction = "negative"
  )
  explicit_log_lik <- function(row, value, mu_basis) {

    candidate_mu <- mu[row, ] + (value - current[[row]]) * mu_basis
    sum(stats::dnorm(
      -yi,
      mean = -candidate_mu,
      sd   = sqrt(tau[row, ]^2 + sei^2),
      log  = TRUE
    ))
  }

  cases <- list(
    list(parameter = "PET",   mu_basis = -sei),
    list(parameter = "PEESE", mu_basis = -sei^2)
  )
  for (case in cases) {
    mu_basis <- case[["mu_basis"]]
    row_states <- lapply(seq_len(nrow(mu)), function(row) {
      list(
        active_setup     = active_setup,
        baseline_log_lik = explicit_log_lik(row, current[[row]], mu_basis)
      )
    })
    basis <- list(
      formula_mu      = FALSE,
      formula_logtau  = FALSE,
      scale_update    = "none",
      log_tau_basis   = NULL,
      mu_basis        = matrix(mu_basis, nrow = nrow(mu), ncol = length(sei),
                               byrow = TRUE),
      current         = current
    )

    testthat::local_mocked_bindings(
      .iwmde_predictor_log_prior = function(context, parameter, values,
                                            row_states, replacement) {

        rep(0, length(values) * length(row_states))
      },
      .package = "RoBMA"
    )

    fast <- .iwmde_log_q_grid_normal_location_group(
      context     = context,
      parameter   = case[["parameter"]],
      values      = values,
      row_states  = row_states,
      replacement = list(type = "scalar"),
      setup       = setup,
      basis       = basis
    )
    reference <- matrix(
      vapply(seq_along(row_states), function(row) {
        vapply(values, explicit_log_lik, numeric(1), row = row,
               mu_basis = mu_basis)
      }, numeric(length(values))),
      nrow = length(values)
    )

    expect_equal(fast, reference, tolerance = 1e-12)
  }
})


test_that("negative-direction selected-normal normalizer change matches matrix reference", {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .7, .35))
  )
  yi  <- c(.10, -.20, .30)
  sei <- c(.20, .30, .40)
  data <- list(outcome = data.frame(yi = yi, sei = sei))
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "negative"

  S      <- 2L
  G      <- 3L
  values <- c(-.30, .05, .40)
  current <- c(.10, -.15)

  setup <- list(
    mu                = matrix(c(.10, .15, .20, -.05, .00, .05),
                               nrow = S, byrow = TRUE),
    tau_within        = matrix(c(.25, .30, .35, .28, .32, .36),
                               nrow = S, byrow = TRUE),
    sei               = sei,
    weights           = c(1, .5, 2),
    posterior_samples = matrix(
      0,
      nrow     = S,
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    )
  )
  active_setup <- list(
    priors            = list(outcome = list(bias = prior)),
    is_weightfunction = TRUE
  )
  context <- list(
    object          = list(fit = list()),
    data            = data,
    predictor_cache = new.env(parent = emptyenv())
  )
  selection_context <- .iwmde_selection_context_active_branch(
    context           = context,
    active_setup      = active_setup,
    posterior_samples = setup[["posterior_samples"]]
  )
  basis <- list(
    formula_logtau = FALSE,
    scale_update   = "none",
    log_tau_basis  = NULL,
    mu_basis       = matrix(-sei, nrow = S, ncol = length(sei), byrow = TRUE),
    current        = current
  )
  row_states <- lapply(seq_len(S), function(row) {
    list(row_index = row, active_key = "negative-selection")
  })

  fast <- .iwmde_selected_normal_location_normalizer_change(
    context      = context,
    setup        = setup,
    basis        = basis,
    values       = values,
    row_states   = row_states,
    active_setup = active_setup
  )

  row_index  <- rep(seq_len(S), each = G)
  grid_index <- rep(seq_len(G), times = S)
  delta      <- values[grid_index] - current[row_index]
  sd <- sqrt(setup[["tau_within"]]^2 +
    matrix(sei^2, nrow = S, ncol = length(sei), byrow = TRUE))
  current_log_norm <- .selection_step_log_norm_matrix(
    mean              = setup[["mu"]],
    sd                = sd,
    sei               = sei,
    selection_context = selection_context
  )
  candidate_context <- BayesTools::selection_context_subset_rows(
    context = selection_context,
    rows    = row_index
  )
  candidate_log_norm <- .selection_step_log_norm_matrix(
    mean              = setup[["mu"]][row_index, , drop = FALSE] +
      basis[["mu_basis"]][row_index, , drop = FALSE] * delta,
    sd                = sd[row_index, , drop = FALSE],
    sei               = sei,
    selection_context = candidate_context
  )
  reference <- rowSums(.apply_log_lik_weights(
    candidate_log_norm - current_log_norm[row_index, , drop = FALSE],
    setup[["weights"]]
  ))

  expect_equal(fast, reference, tolerance = 1e-12)
})


test_that("multilevel weightfunction location fast path dispatches to cluster selected-normal grid", {

  yi  <- c(.10, .20)
  sei <- c(.20, .25)
  data <- list(outcome = data.frame(yi = yi, sei = sei, cluster = c(1L, 1L)))
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "positive"
  attr(data, "cluster")          <- TRUE

  values <- c(-.10, .30)
  setup <- list(
    yi                = yi,
    sei               = sei,
    mu                = matrix(c(.05, .10, .15, .20), nrow = 2L, byrow = TRUE),
    tau_within        = matrix(.20, nrow = 2L, ncol = 2L),
    tau_between       = matrix(.10, nrow = 2L, ncol = 2L),
    cluster           = list(`1` = 1:2),
    posterior_samples = matrix(0, nrow = 2L, ncol = 1L)
  )
  basis <- list(
    formula_mu      = FALSE,
    formula_logtau  = FALSE,
    scale_update    = "none",
    log_tau_basis   = NULL,
    mu_basis        = matrix(1, nrow = 2L, ncol = 2L),
    current         = c(0, .10)
  )
  row_states <- lapply(seq_len(2L), function(row) {
    list(
      row_index    = row,
      active_setup = list(is_weightfunction = TRUE)
    )
  })
  log_lik <- matrix(c(1, 2, 3, 4), nrow = length(values))
  log_prior <- c(.1, .2, .3, .4)
  captured <- NULL

  testthat::local_mocked_bindings(
    .iwmde_selection_context_active_branch = function(context, active_setup,
                                                      posterior_samples) {

      list(marker = TRUE)
    },
    .log_lik_cluster_selnorm_location_grid = function(setup, yi, sei, basis,
                                                      current, values,
                                                      selection_context) {

      captured <<- list(
        yi                = yi,
        sei               = sei,
        basis             = basis,
        current           = current,
        values            = values,
        selection_context = selection_context
      )
      log_lik
    },
    .iwmde_predictor_log_prior = function(context, parameter, values,
                                          row_states, replacement) {

      log_prior
    },
    .package = "RoBMA"
  )

  fast <- .iwmde_log_q_grid_normal_location_group(
    context     = list(data = data),
    parameter   = "mu",
    values      = values,
    row_states  = row_states,
    replacement = list(type = "scalar"),
    setup       = setup,
    basis       = basis
  )

  expect_equal(fast, log_lik + matrix(log_prior, nrow = length(values)))
  expect_equal(captured[["yi"]], yi)
  expect_equal(captured[["sei"]], sei)
  expect_equal(captured[["basis"]], basis[["mu_basis"]])
  expect_equal(captured[["current"]], basis[["current"]])
  expect_equal(captured[["values"]], values)
  expect_true(isTRUE(captured[["selection_context"]][["marker"]]))
})


test_that("multilevel weightfunction location fast path matches scalar fallback", {

  fit_name  <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  parameter <- "mu_Preregistered"
  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, NULL)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep <- vapply(row_states, function(state) {
    isTRUE(state[["active_setup"]][["is_weightfunction"]]) &&
      identical(state[["likelihood_mode"]], "marginal") &&
      is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  keys <- vapply(row_states, function(state) {
    .iwmde_state_active_key(context, state)
  }, character(1))
  active_key <- names(sort(table(keys), decreasing = TRUE))[[1L]]
  row_states <- head(row_states[keys == active_key], 3L)

  expect_gt(length(row_states), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = 3L),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement  <- .iwmde_replacement_spec(context, parameter, spec)
  active_setup <- row_states[[1L]][["active_setup"]]
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = "cluster"
  )
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )

  expect_false(is.null(basis))

  basis <- .iwmde_predictor_resolve_formula_basis(
    context      = context,
    active_setup = active_setup,
    setup        = setup,
    basis        = basis,
    row_states   = row_states
  )
  fast <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  scalar <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  finite_terms <- is.finite(fast) & is.finite(scalar)
  expect_true(any(finite_terms))
  expect_equal(fast[finite_terms], scalar[finite_terms], tolerance = 1e-8)
})


test_that("negative-direction multilevel selected-normal location grid matches quadrature reference", {

  skip_if_not(.has_native_selnorm_cluster_location_grid())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .65, .30))
  )
  yi  <- c(.12, -.18, .31, -.05)
  sei <- c(.20, .28, .24, .18)
  selection_context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "negative"
  )

  expect_equal(selection_context[["sign"]], -1L)

  S <- 2L
  K <- length(yi)
  setup <- list(
    S           = S,
    mu          = matrix(c(
       .05, .08, .12, .16,
      -.04, .02, .09, .13
    ), nrow = S, byrow = TRUE),
    tau_within  = matrix(c(
      .09, .11, .13, .15,
      .08, .10, .12, .14
    ), nrow = S, byrow = TRUE),
    tau_between = matrix(c(
      .05, .07, .09, .11,
      .04, .06, .08, .10
    ), nrow = S, byrow = TRUE),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L)),
    weights     = NULL
  )
  selection_context[["omega"]] <- matrix(c(
    1, .65, .30,
    1, .55, .25
  ), nrow = S, byrow = TRUE)
  selection_context[["alpha"]]       <- rep(0, S)
  selection_context[["phack_kind"]]  <- rep(0L, S)
  selection_context[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)

  basis <- matrix(c(
     .04, -.03, .05, -.02,
    -.02,  .06, .01,  .04
  ), nrow = S, byrow = TRUE)
  current <- c(-.10, .15)
  values  <- c(-.20, .05, .30)
  n_gamma <- 11L

  fast <- .log_lik_cluster_selnorm_location_grid(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    basis             = basis,
    current           = current,
    values            = values,
    selection_context = selection_context,
    n_gamma           = n_gamma
  )
  reference <- matrix(NA_real_, nrow = length(values), ncol = S)
  for (s in seq_len(S)) {
    context_s <- BayesTools::selection_context_subset_rows(
      context = selection_context,
      rows    = s
    )
    for (g in seq_along(values)) {
      candidate_setup <- setup
      delta <- values[[g]] - current[[s]]
      candidate_setup[["S"]] <- 1L
      candidate_setup[["mu"]] <- setup[["mu"]][s, , drop = FALSE] +
        basis[s, , drop = FALSE] * delta
      candidate_setup[["tau_within"]]  <- setup[["tau_within"]][s, , drop = FALSE]
      candidate_setup[["tau_between"]] <- setup[["tau_between"]][s, , drop = FALSE]
      reference[g, s] <- sum(.log_lik_cluster_norm_quadrature_r(
        setup             = candidate_setup,
        yi                = yi,
        sei               = sei,
        is_weightfunction = TRUE,
        selection_context = context_s,
        n_gamma           = n_gamma
      ))
    }
  }

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(reference))
  quadrature_change <- attr(fast, "quadrature_relative_change", exact = TRUE)
  expect_equal(dim(quadrature_change), dim(fast))
  expect_true(all(is.finite(quadrature_change)))
  attr(fast, "quadrature_relative_change") <- NULL
  expect_equal(fast, reference, tolerance = 1e-7)
})


test_that("IWMDE selected-normal native normalizer skips fallback candidates", {

  G <- 3L
  S <- 2L
  K <- 2L
  setup <- list(
    mu                = matrix(c(.1, .2, .3, .4), nrow = S),
    tau_within        = c(.05, .10),
    sei               = c(.20, .30),
    weights           = NULL,
    posterior_samples = matrix(0, nrow = S, ncol = 1)
  )
  basis <- list(
    log_tau_basis  = NULL,
    scale_update   = "none",
    formula_logtau = FALSE,
    mu_basis       = NULL,
    current        = c(.10, .30)
  )
  row_states <- list(
    list(row_index = 1L),
    list(row_index = 2L)
  )
  selection_context <- list(
    omega       = matrix(1, nrow = S, ncol = 1),
    alpha       = rep(NA_real_, S),
    phack_kind  = rep(NA_integer_, S),
    kernel_mode = rep(SELKERNEL_NORMAL, S)
  )

  testthat::local_mocked_bindings(
    .iwmde_selection_context_active_branch = function(context, active_setup,
                                                      posterior_samples) {

      selection_context
    },
    .iwmde_selected_normal_current_log_norm = function(context, setup, sd,
                                                       selection_context,
                                                       row_states) {

      matrix(0, nrow = S, ncol = K)
    },
    .has_native_selnorm_log_norm_delta = function() TRUE,
    .selnorm_kernel_log_norm_delta_grid = function(
        mean, sd, basis, current_log_norm, current, values, sei, weights,
        omega, selection_spec, alpha, phack_kind, kernel_mode) {

      expect_equal(mean, setup[["mu"]])
      expect_null(basis)
      matrix(.25, nrow = length(values), ncol = length(current))
    },
    .package = "RoBMA"
  )

  out <- .iwmde_selected_normal_location_normalizer_change(
    context      = list(),
    setup        = setup,
    basis        = basis,
    values       = seq(-1, 1, length.out = G),
    row_states   = row_states,
    active_setup = list()
  )

  expect_equal(out, rep(.25, G * S))
})


test_that("IWMDE active-branch samples reject uncollapsed prior mixtures", {

  prior <- BayesTools::prior_mixture(
    prior_list = list(
      BayesTools::prior("spike", parameters = list(location = 0)),
      BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
    ),
    is_null = c(TRUE, FALSE)
  )
  samples <- matrix(
    2,
    nrow     = 1,
    dimnames = list(NULL, "mu_indicator")
  )

  expect_error(
    .iwmde_likelihood_posterior_samples(
      context      = list(),
      samples      = samples,
      active_setup = list(
        priors            = list(outcome = list(mu = prior)),
        is_weightfunction = FALSE
      )
    ),
    "requires priors localized to a single mixture component"
  )
})

test_that("IWMDE batched q evaluation warns on invalid likelihood length", {

  context <- list(
    posterior_samples = matrix(
      0,
      nrow     = 1,
      dimnames = list(NULL, "mu")
    ),
    row_cache = new.env(parent = emptyenv())
  )
  row_states <- list(
    list(
      row_index       = 1L,
      active_key      = "same",
      active_setup    = list(),
      likelihood_mode = "marginal",
      prior_list      = list()
    )
  )

  expect_warning(
    out <- .iwmde_log_q_grid_from_samples(
      context         = context,
      parameter       = "mu",
      values          = c(0, 1),
      row_states      = row_states,
      replacement     = list(type = "scalar"),
      likelihood_mode = "marginal",
      log_lik_fun     = function(samples, active_setup) 0
    ),
    "invalid length"
  )
  expect_null(out)
})


test_that("IWMDE replacement keeps inverse-gamma auxiliary parameters synced", {

  prior <- BayesTools::prior(
    "invgamma",
    parameters = list(2, .5)
  )
  context <- list(flat_prior_list = list(tau = prior))
  samples <- matrix(
    c(.5, 99,
      2, 99),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("tau", "inv_tau"))
  )

  synced <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "tau"
  )
  bad <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "mu"
  )
  row <- .iwmde_sync_invgamma_auxiliary_row(
    context    = context,
    row        = c(tau = 4, inv_tau = 99),
    parameters = "tau"
  )
  invalid <- .iwmde_sync_invgamma_auxiliary_row(
    context    = context,
    row        = c(tau = 0, inv_tau = 99),
    parameters = "tau"
  )

  expect_false(.iwmde_can_use_focal_prior_delta(prior))
  expect_true(all(synced[["valid"]]))
  expect_equal(synced[["samples"]][, "inv_tau"], 1 / samples[, "tau"])
  expect_equal(bad[["samples"]], samples)
  expect_true(isTRUE(row[["valid"]]))
  expect_equal(row[["row"]][["inv_tau"]], .25)
  expect_false(invalid[["valid"]])
  expect_true(is.na(invalid[["row"]][["inv_tau"]]))
})


test_that("IWMDE inverse-gamma auxiliary matrix sync reuses active states", {

  inv_prior <- BayesTools::prior(
    "invgamma",
    parameters = list(2, .5)
  )
  normal_prior <- BayesTools::prior(
    "normal",
    parameters = list(0, 1)
  )
  context <- list(
    flat_prior_list = list(
      tau = BayesTools::prior_mixture(list(inv_prior, normal_prior))
    ),
    indicator_names = "tau_indicator"
  )
  samples <- matrix(
    c(.5, 99, 1,
      2, 99, 1,
      0, 99, 1,
      0, 99, 2,
      4, 99, 2,
      5, 99, 1),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(NULL, c("tau", "inv_tau", "tau_indicator"))
  )

  calls <- 0L
  testthat::local_mocked_bindings(
    .iwmde_row_uses_invgamma_prior = function(context, parameter, row) {
      calls <<- calls + 1L
      row[["tau_indicator"]] == 1
    },
    .package = "RoBMA"
  )
  synced <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = "tau"
  )

  expect_equal(calls, 2L)
  expect_equal(
    synced[["samples"]][c(1, 2, 6), "inv_tau"],
    1 / samples[c(1, 2, 6), "tau"]
  )
  expect_true(is.na(synced[["samples"]][3, "inv_tau"]))
  expect_equal(synced[["samples"]][c(4, 5), "inv_tau"], samples[c(4, 5), "inv_tau"])
  expect_equal(synced[["valid"]], c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
})


test_that("IWMDE batched q evaluation matches scalar fallback", {

  cases <- list(
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("bangertdrowns2004_location-scale", "log_tau_intercept", NULL),
    list("dat.lehmann2018-3PSM", "mu", NULL),
    list("dat.lehmann2018-3PSM_neg", "mu", NULL),
    list("dat.lehmann2018-PET_neg", "PET", NULL),
    list("dat.lehmann2018-PEESE_neg", "PEESE", NULL),
    list("dat.lehmann2018_RoBMA_3lvl_mods_scale", "mu_Preregistered", NULL),
    list("bcg_BMA.glmm_3lvl_location_scale", "mu_year", NULL),
    list("konstantopoulos2011_3lvl", "gamma[1]", NULL),
    list("bcg_glmm", "theta[1]", NULL),
    list("bcg_glmm", "pi[1]", NULL),
    list("nielweise2008_glmm", "phi[1]", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )

  for (case in cases) {
    .expect_iwmde_batch_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})

test_that("IWMDE predictor fast path matches scalar fallback", {

  cases <- list(
    list("bcg_meta-analysis", "mu", NULL),
    list("brma.mv_block_mvn_fixed_random_null", "mu", NULL),
    list("brma.mv_block_mvn", "mu", NULL),
    list("bcg_meta-analysis", "tau", NULL),
    list("konstantopoulos2011_3lvl", "rho", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
    list("dat.lehmann2018-PET_neg", "PET", NULL),
    list("dat.lehmann2018-PEESE_neg", "PEESE", NULL),
    list("dat.lehmann2018-3PSM_neg", "mu", NULL),
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("bangertdrowns2004_location-scale", "log_tau_intercept", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )

  for (case in cases) {
    .expect_iwmde_predictor_fast_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})

test_that("IWMDE normal quadratic fast path matches scalar fallback", {

  cases <- list(
    list("bcg_meta-analysis", "mu", NULL),
    list("brma.mv_block_mvn_fixed_random_null", "mu", NULL),
    list("brma.mv_block_mvn", "mu", NULL),
    list("dat.lehmann2018-3PSM", "mu", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
    list("dat.lehmann2018-PET_neg", "PET", NULL),
    list("dat.lehmann2018-PEESE_neg", "PEESE", NULL),
    list("dat.lehmann2018-3PSM_neg", "mu", NULL),
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("konstantopoulos2011_3lvl", "mu", NULL),
    list(
      "bcg_meta-regression",
      "mu_intercept + mu_ablat",
      list(type = "linear", weights = c(mu_intercept = 1, mu_ablat = 1))
    )
  )

  for (case in cases) {
    .expect_iwmde_normal_quadratic_equals_scalar(
      fit_name       = case[[1L]],
      parameter      = case[[2L]],
      parameter_spec = case[[3L]]
    )
  }
})


test_that("known-V normal q-grid factors invariant blocks once", {

  fit_name <- "brma.mv_block_mvn_fixed_random_null"
  .skip_if_missing_raw_fits(fit_name)

  context   <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  parameter <- "mu"
  spec      <- .iwmde_parameter_spec(context, parameter, NULL)
  values    <- .iwmde_parameter_values(context, parameter, spec)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & is.finite(values))
  rows      <- head(rows, 3L)
  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  active_setup <- row_states[[1L]][["active_setup"]]
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = "estimate"
  )
  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )
  grid <- seq(min(values[rows]), max(values[rows]), length.out = 20L)

  original_factor <- .known_v_chol_covariance
  factor_calls     <- 0L
  testthat::local_mocked_bindings(
    .known_v_chol_covariance = function(covariance, context) {

      factor_calls <<- factor_calls + 1L
      original_factor(covariance = covariance, context = context)
    },
    .package = "RoBMA"
  )

  out <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = grid,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  reused <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = rev(grid),
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )

  block_count <- length(.known_v_blocks(.data_known_v_data(context[["data"]])))
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(length(grid), length(row_states)))
  expect_equal(reused, out[nrow(out):1L, , drop = FALSE], tolerance = 1e-12)
  expect_equal(factor_calls, block_count)
})
