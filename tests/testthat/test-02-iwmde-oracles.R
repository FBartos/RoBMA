context("IWMDE numerical oracles")

source(testthat::test_path("common-functions.R"))


.iwmde_oracle_skip_missing_fits <- function(names) {

  skip_if_missing_fits(names)
}


.iwmde_oracle_point_estimate <- function(fit, parameter, value,
                                         density_method, density_control) {

  context  <- .iwmde_context(fit)
  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = parameter,
    density_method  = density_method,
    density_control = c(
      density_control,
      list(display_grid = "ordinate")
    ),
    outputs         = "ordinate",
    values          = value,
    parameter_spec  = list(
      type             = "primitive",
      conditional      = NULL,
      conditional_rule = "AND"
    ),
    metadata        = list(parameter = parameter),
    cache           = .iwmde_estimate_cache()
  )

  return(list(context = context, estimate = estimate))
}


# ============================================================================ #
# Test: Closed-form conditional-normal ordinary random-effects oracle
# ============================================================================ #

test_that("qCMDE matches the conditional-normal random-effects oracle", {

  skip_on_cran()
  .iwmde_oracle_skip_missing_fits("bcg_meta-analysis")

  fit <- load_fit("bcg_meta-analysis")
  out <- .iwmde_oracle_point_estimate(
    fit             = fit,
    parameter       = "mu",
    value           = 0,
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 40,
      max_samples          = 500,
      normalization_points = 101,
      normalization_prob   = .9999
    )
  )

  rows  <- out[["estimate"]][["plan"]][["rows"]][["estimator_rows"]]
  tau   <- out[["context"]][["posterior_samples"]][rows, "tau"]
  yi    <- fit[["data"]][["outcome"]][["yi"]]
  sei   <- fit[["data"]][["outcome"]][["sei"]]
  prior <- fit[["priors"]][["outcome"]][["mu"]][["parameters"]]

  observation_precision <- vapply(
    tau,
    function(tau_value) sum(1 / (sei^2 + tau_value^2)),
    numeric(1)
  )
  weighted_outcome       <- vapply(
    tau,
    function(tau_value) sum(yi / (sei^2 + tau_value^2)),
    numeric(1)
  )
  posterior_variance     <- 1 / (
    1 / prior[["sd"]]^2 + observation_precision
  )
  posterior_mean         <- posterior_variance * (
    prior[["mean"]] / prior[["sd"]]^2 + weighted_outcome
  )
  oracle_ordinate        <- mean(stats::dnorm(
    0,
    mean = posterior_mean,
    sd   = sqrt(posterior_variance)
  ))
  estimated_ordinate     <- out[["estimate"]][["posterior_ordinate"]][["ordinate"]]

  expect_equal(
    log(estimated_ordinate),
    log(oracle_ordinate),
    tolerance = 1e-3,
    info = "qCMDE ordinate must equal the average exact conditional density"
  )
})


# ============================================================================ #
# Test: Exact conjugate known-V oracle
# ============================================================================ #

test_that("known-V posterior ordinate and marginal likelihood are exact", {

  skip_on_cran()
  .iwmde_oracle_skip_missing_fits("brma.mv_block_mvn_fixed_random_null")

  fit   <- load_fit("brma.mv_block_mvn_fixed_random_null")
  V     <- .known_v_materialize(attr(fit[["data"]], "known_V_data"))
  yi    <- fit[["data"]][["outcome"]][["yi"]]
  prior <- fit[["priors"]][["outcome"]][["mu"]][["parameters"]]

  V_inverse         <- solve(V)
  intercept         <- rep(1, length(yi))
  posterior_variance <- 1 / (
    1 / prior[["sd"]]^2 +
      as.numeric(crossprod(intercept, V_inverse %*% intercept))
  )
  posterior_mean     <- posterior_variance * (
    prior[["mean"]] / prior[["sd"]]^2 +
      as.numeric(crossprod(intercept, V_inverse %*% yi))
  )
  oracle_ordinate    <- stats::dnorm(
    0,
    mean = posterior_mean,
    sd   = sqrt(posterior_variance)
  )
  oracle_log_marglik <- mvtnorm::dmvnorm(
    yi,
    mean  = rep(prior[["mean"]], length(yi)),
    sigma = V + prior[["sd"]]^2,
    log   = TRUE
  )
  log_kernel <- function(mu) {

    vapply(mu, function(value) {
      log_lik <- .log_lik_known_v_joint_sum_from_evaluated_predictors(
        fit                = fit[["fit"]],
        data               = fit[["data"]],
        priors             = fit[["priors"]],
        mu_samples         = matrix(value, nrow = 1L, ncol = length(yi)),
        tau_within_samples = matrix(0, nrow = 1L, ncol = length(yi)),
        unit                = "estimate"
      )
      as.numeric(log_lik) + stats::dnorm(
        value,
        mean = prior[["mean"]],
        sd   = prior[["sd"]],
        log  = TRUE
      )
    }, numeric(1))
  }
  mode <- stats::optimize(
    function(mu) -log_kernel(mu),
    interval = prior[["mean"]] + c(-8, 8) * prior[["sd"]]
  )
  integrated <- stats::integrate(
    function(mu) exp(log_kernel(mu) + mode[["objective"]]),
    lower         = -Inf,
    upper         = Inf,
    rel.tol       = 1e-10,
    subdivisions  = 1000L,
    stop.on.error = TRUE
  )[["value"]]
  package_log_marglik <- log(integrated) - mode[["objective"]]

  expect_equal(
    package_log_marglik,
    as.numeric(oracle_log_marglik),
    tolerance = 1e-7,
    info = "known-V likelihood integrated over the prior"
  )

  controls <- list(
    qCMDE = list(
      n_points             = 60,
      max_samples          = 40,
      normalization_points = 101,
      normalization_prob   = .9999
    ),
    IWMDE = list(
      n_points             = 60,
      max_samples          = 240,
      normalization_points = 200,
      normalization_prob   = .9999
    )
  )

  for (density_method in names(controls)) {
    out <- .iwmde_oracle_point_estimate(
      fit             = fit,
      parameter       = "mu",
      value           = 0,
      density_method  = density_method,
      density_control = controls[[density_method]]
    )
    estimated_ordinate <- out[["estimate"]][["posterior_ordinate"]][["ordinate"]]

    expect_equal(
      log(estimated_ordinate),
      log(oracle_ordinate),
      tolerance = 0.01,
      info = paste(density_method, "must match the exact known-V ordinate")
    )
  }
})


# ============================================================================ #
# Test: Bridge oracles for GLMM and active selection likelihoods
# ============================================================================ #

.iwmde_oracle_bridge_mcse <- function(object) {

  repetitions <- object[["marglik"]][["repetitions"]]
  included    <- repetitions[["success"]] & repetitions[["finite"]]
  if (!any(included)) {
    return(NA_real_)
  }

  return(max(repetitions[["mcse"]][included]))
}


.iwmde_oracle_expect_bridge <- function(result, null, full,
                                         hard_tolerance = .10,
                                         info = NULL) {

  raw_bf           <- as.numeric(attr(result, "raw_BF", exact = TRUE))[[1L]]
  bf_error_percent <- as.numeric(result[["BF_error"]])[[1L]]
  log_difference   <- abs(log(raw_bf) - (logml(full) - logml(null)))
  estimator_mcse   <- sqrt(log1p((bf_error_percent / 100)^2))
  bridge_mcse      <- sqrt(
    .iwmde_oracle_bridge_mcse(full)^2 +
      .iwmde_oracle_bridge_mcse(null)^2
  )

  expect_lte(log_difference, hard_tolerance, label = info)
  expect_lte(
    log_difference,
    3 * sqrt(estimator_mcse^2 + bridge_mcse^2),
    label = info
  )
}


test_that("qCMDE matches GLMM and both estimators match selection bridge factors", {

  skip_on_cran()
  pairs <- list(
    glmm = c(
      "nielweise2008_glmm_effect_null",
      "nielweise2008_glmm"
    ),
    selection = c(
      "dat.lehmann2018-3PSM_effect_null",
      "dat.lehmann2018-3PSM"
    )
  )
  .iwmde_oracle_skip_missing_fits(unique(unlist(pairs, use.names = FALSE)))

  for (pair_name in names(pairs)) {
    null <- load_fit(pairs[[pair_name]][[1L]])
    full <- load_fit(pairs[[pair_name]][[2L]])
    .expect_bridge_nesting(null, full, "mu")
    expect_true(is.finite(logml(null)))
    expect_true(is.finite(logml(full)))
    expect_lt(.iwmde_oracle_bridge_mcse(null), .05)
    expect_lt(.iwmde_oracle_bridge_mcse(full), .05)

    # GLMM local-state IWMDE has a high-dimensional fitted weight function;
    # qCMDE provides the stable fixed-runtime bridge certification here.
    density_methods <- if (identical(pair_name, "glmm")) {
      "qCMDE"
    } else {
      c("qCMDE", "IWMDE")
    }
    for (density_method in density_methods) {
      result <- hypothesis(
        full,
        "mu = 0",
        columns         = "all",
        density_method  = density_method,
        density_control = list(
          n_points             = 40,
          initial_samples      = 2000,
          max_samples          = 2000,
          normalization_points = 80
        ),
        n_samples = 1000
      )

      .iwmde_oracle_expect_bridge(
        result,
        null,
        full,
        info = paste(pair_name, density_method)
      )
      diagnostics <- density_diagnostics(result)
      expect_equal(diagnostics[["achieved_row_budget"]], 2000L)
      expect_equal(diagnostics[["status"]], "ok")
    }
  }
})


test_that("GLMM IWMDE point Bayes factors are explicitly unsupported", {

  skip_on_cran()
  names <- c("nielweise2008_glmm_effect_null", "nielweise2008_glmm")
  .iwmde_oracle_skip_missing_fits(names)

  null <- load_fit(names[[1L]])
  full <- load_fit(names[[2L]])
  .expect_bridge_nesting(null, full, "mu")

  expect_error(
    hypothesis(
      full,
      "mu = 0",
      density_method = "IWMDE"
    ),
    "do not meet the bridge-sampling certification tolerance"
  )
})


test_that("qCMDE and IWMDE match the known-V tau boundary bridge factor", {

  skip_on_cran()
  names <- c("iwmde_known_v_tau_null", "iwmde_known_v_tau_full")
  .iwmde_oracle_skip_missing_fits(names)

  null <- load_fit(names[[1L]])
  full <- load_fit(names[[2L]])
  .expect_bridge_nesting(null, full, "tau")
  expect_true(is.finite(logml(null)))
  expect_true(is.finite(logml(full)))
  expect_lt(.iwmde_oracle_bridge_mcse(null), .05)
  expect_lt(.iwmde_oracle_bridge_mcse(full), .05)

  for (density_method in c("qCMDE", "IWMDE")) {
    result <- hypothesis(
      full,
      "tau = 0",
      columns         = "all",
      density_method  = density_method,
      density_control = list(
        n_points             = 60,
        initial_samples      = 240,
        max_samples          = 240,
        normalization_points = 120,
        normalization_prob   = .9999
      ),
      n_samples = 1000
    )

    .iwmde_oracle_expect_bridge(
      result,
      null,
      full,
      info = paste("known-V tau boundary", density_method)
    )
    diagnostics <- density_diagnostics(result)
    expect_equal(diagnostics[["achieved_row_budget"]], 240L)
    expect_equal(diagnostics[["status"]], "ok")
  }
})
