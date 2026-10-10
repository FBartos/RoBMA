context("Cached fits for the bias-member design tests")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("bias-branch-design")

# Ensembles whose members differ in their bias terms (null, PET and PEESE
# members; no member holds both columns), and the data they were fitted to.
# test-02-bias-branch-design.R reads the draws of each design class from them.

.bbd_data <- function() {

  withr::local_seed(2)
  k <- 20
  d <- data.frame(sei = stats::runif(k, 0.08, 0.30), x = stats::rnorm(k))
  d[["yi"]] <- 0.25 + 0.1 * d[["x"]] + 0.5 * d[["sei"]] +
    stats::rnorm(k, 0, sqrt(0.02^2 + d[["sei"]]^2))
  d[["ni"]] <- 50

  return(d)
}

.bbd_fit_settings <- function() {
  list(
    chains = 1, burnin = 500, adapt = 300, sample = 1000, seed = 1,
    silent = TRUE, parallel = FALSE,
    convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  )
}

# 12 estimates in six studies with correlated sampling errors
.bbd_known_v_input <- function() {

  withr::local_seed(3)
  k_study <- 6L
  sei     <- stats::runif(2L * k_study, 0.08, 0.30)
  x       <- stats::rnorm(2L * k_study)
  d <- data.frame(
    sei = sei, x = x,
    yi  = 0.2 + 0.1 * x + 0.5 * sei + stats::rnorm(2L * k_study, 0, sei)
  )
  V <- matrix(0, 2L * k_study, 2L * k_study)
  for (s in seq_len(k_study)) {
    index <- 2L * s - 1:0
    V[index, index] <- 0.4 * tcrossprod(sei[index])
    diag(V[index, index]) <- sei[index]^2
  }

  return(list(data = d, V = V))
}


test_that("PP ensembles and a known-V ensemble are fitted for the bias-member design tests", {

  bbd_data <- .bbd_data()
  args <- c(
    list(
      yi = quote(yi), sei = quote(sei), ni = quote(ni), mods = ~ x,
      data = bbd_data, measure = "GEN", standardize_continuous_predictors = FALSE
    ),
    .bbd_fit_settings()
  )

  fit_mods <- suppressWarnings(do.call(RoBMA, c(args, list(model_type = "PP"))))
  fit_mods <- suppressWarnings(add_loo(fit_mods))
  fit_pet  <- suppressWarnings(do.call(bPET, args))
  intercept_args <- c(
    args[c("yi", "sei", "ni", "data", "measure")], .bbd_fit_settings()
  )
  fit_intercept <- suppressWarnings(
    do.call(RoBMA, c(intercept_args, list(model_type = "PP")))
  )

  known_v <- .bbd_known_v_input()
  fit_known_v <- suppressWarnings(RoBMA.mv(
    yi = yi, V = known_v[["V"]], mods = ~ x, data = known_v[["data"]],
    measure = "GEN",
    prior_bias = list(
      BayesTools::prior_PET("cauchy", list(0, 1), list(0, Inf)),
      BayesTools::prior_PEESE("cauchy", list(0, 1), list(0, Inf))
    ),
    prior_unit_information_sd = 1, standardize_continuous_predictors = FALSE,
    chains = 1, sample = 400, burnin = 200, adapt = 200, seed = 1, silent = TRUE,
    parallel = FALSE,
    convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  ))

  save_fit("bias_branch_PP_mods", fit_mods, info = list(data = bbd_data))
  save_fit("bias_branch_PP_intercept", fit_intercept, info = list(data = bbd_data))
  save_fit("bias_branch_bPET_mods", fit_pet, info = list(data = bbd_data))
  save_fit("bias_branch_known_v_ensemble", fit_known_v, info = known_v)

  expect_true(is_cached_fit_valid("bias_branch_PP_mods"))
  expect_true(is_cached_fit_valid("bias_branch_PP_intercept"))
  expect_true(is_cached_fit_valid("bias_branch_bPET_mods"))
  expect_true(is_cached_fit_valid("bias_branch_known_v_ensemble"))
})
