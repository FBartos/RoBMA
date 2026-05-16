# ============================================================================ #
# test-02-iwmde-fast-paths.R
# ============================================================================ #

context("IWMDE fast-path equivalence")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


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
})

test_that("IWMDE omega matrix collapse matches row-wise collapse", {

  omega       <- matrix(seq_len(12), nrow = 3, byrow = TRUE)
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
})

test_that("IWMDE batched q evaluation matches scalar fallback", {

  cases <- list(
    list("bcg_meta-regression", "mu_ablat", NULL),
    list("bangertdrowns2004_location-scale", "log_tau_intercept", NULL),
    list("dat.lehmann2018-3PSM", "mu", NULL),
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
    list("bcg_meta-analysis", "tau", NULL),
    list("konstantopoulos2011_3lvl", "rho", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
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
    list("dat.lehmann2018-3PSM", "mu", NULL),
    list("dat.lehmann2018-PET", "PET", NULL),
    list("dat.lehmann2018-PEESE", "PEESE", NULL),
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
