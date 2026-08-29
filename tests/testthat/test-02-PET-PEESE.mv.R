source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))

skip_if_missing_fits(c("bPET.mv_random", "bPEESE.mv_random"))
fit_pet_mv   <- load_fit("bPET.mv_random", validate = FALSE)
fit_peese_mv <- load_fit("bPEESE.mv_random", validate = FALSE)


test_that("multivariate PET and PEESE summaries retain both model parts", {

  pet_summary   <- summary(fit_pet_mv, include_mcmc_diagnostics = FALSE)
  peese_summary <- summary(fit_peese_mv, include_mcmc_diagnostics = FALSE)

  expect_match(pet_summary[["name"]], "Multivariate.*PET")
  expect_match(peese_summary[["name"]], "Multivariate.*PEESE")
  expect_true(nrow(pet_summary[["estimates_random"]]) > 0L)
  expect_true(nrow(peese_summary[["estimates_random"]]) > 0L)
  expect_true(nrow(pet_summary[["estimates_bias"]]) > 0L)
  expect_true(nrow(peese_summary[["estimates_bias"]]) > 0L)
  expect_true(all(c("random", "bias") %in%
                    as.data.frame(pet_summary)[["component"]]))
  expect_true(all(c("random", "bias") %in%
                    as.data.frame(peese_summary)[["component"]]))
})


test_that("multivariate PET and PEESE predictions use their row predictors", {

  specifications <- list(
    PET = list(
      object    = fit_pet_mv,
      parameter = "PET",
      predictor = function(sei) sei
    ),
    PEESE = list(
      object    = fit_peese_mv,
      parameter = "PEESE",
      predictor = function(sei) sei^2
    )
  )

  for (bias_type in names(specifications)) {
    specification <- specifications[[bias_type]]
    object         <- specification[["object"]]
    posterior      <- .get_posterior_samples(object[["fit"]])
    adjusted <- predict(
      object,
      type          = "terms",
      bias_adjusted = TRUE
    )
    observed <- predict(
      object,
      type          = "terms",
      bias_adjusted = FALSE
    )
    coefficient <- .posterior_or_fixed_scalar(
      posterior_samples = posterior,
      parameter         = specification[["parameter"]],
      fixed_value       = .fixed_bias_parameter_value(
        object[["priors"]],
        specification[["parameter"]]
      )
    )
    expected <- outer(
      coefficient,
      specification[["predictor"]](
        object[["data"]][["outcome"]][["sei"]]
      )
    )

    expect_equal(
      matrix(as.numeric(observed), nrow = nrow(observed)) -
        matrix(as.numeric(adjusted), nrow = nrow(adjusted)),
      expected,
      tolerance = 1e-12
    )
  }
})


test_that("multivariate PET and PEESE downstream targets remain available", {

  for (object in list(fit_pet_mv, fit_peese_mv)) {
    n       <- nobs(object)
    draws   <- .get_posterior_samples(object[["fit"]])
    log_lik <- log_lik(object, unit = "estimate")

    expect_identical(dim(log_lik), c(nrow(draws), n))
    expect_true(all(is.finite(log_lik)))
    expect_s3_class(loo(object, unit = "estimate"), "loo")
    object_with_waic <- suppressWarnings(
      add_waic(object, unit = "estimate")
    )
    expect_s3_class(waic(object_with_waic, unit = "estimate"), "waic")
    expect_true(is.finite(as.numeric(logml(object))))
    expect_brma_samples_matrix(
      predict(object, type = "estimate", conditioning_depth = "marginal"),
      n,
      "marginal estimate"
    )
    expect_brma_samples_matrix(
      predict(object, type = "estimate", conditioning_depth = "estimate"),
      n,
      "fitted estimate"
    )
    expect_brma_samples_matrix(
      predict(object, type = "response", bias_adjusted = FALSE),
      n,
      "observed response"
    )
    expect_brma_samples_matrix(
      ranef(object, component = "total", expand = TRUE),
      n,
      "fitted random effects"
    )
    expect_length(residuals(object), n)
    expect_equal(nrow(rstudent(object)), n)
    expect_length(hatvalues(object), n)
  }
})


test_that("multivariate PET and PEESE support correlated new responses", {

  newdata <- data.frame(
    study = factor(c("new study", "new study")),
    obs   = factor(c("new 1", "new 2")),
    x     = c(0, 1)
  )
  V_new <- matrix(c(0.02, 0.006, 0.006, 0.03), nrow = 2L)

  for (object in list(fit_pet_mv, fit_peese_mv)) {
    response <- predict(
      object,
      newdata       = newdata,
      V_new         = V_new,
      type          = "response",
      bias_adjusted = FALSE
    )

    expect_brma_samples_matrix(response, 2L, "new observed response")
  }
})
