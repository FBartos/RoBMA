context("Point-test eligibility of scalar formula coefficients")

# hypothesis_quantities() must list a point-test method for a scalar formula
# coefficient exactly when hypothesis() runs a point hypothesis on it with
# that method. With a standardized predictor, the original-scale intercept
# combines the intercept with the slope:
# - two Cauchy priors (intercept and slope): their sum is itself Cauchy,
#   which BayesTools evaluates exactly;
# - the same priors with a raw predictor keep the intercept as the fitted
#   coefficient itself (exact);
# - a normal intercept with a Cauchy slope is a Gaussian convolution, which
#   BayesTools evaluates exactly;
# - model-averaged fits combine the intercept and slope mixtures: their null
#   components put a point mass at 0, and point hypotheses elsewhere use the
#   continuous part, which is exact per component for the default normal
#   priors and for Cauchy alternatives.
.scalar_eligibility_cache <- new.env(parent = emptyenv())

.scalar_eligibility_fits <- function() {

  if (!is.null(.scalar_eligibility_cache[["fits"]])) {
    return(.scalar_eligibility_cache[["fits"]])
  }

  set.seed(1)
  k    <- 40L
  data <- data.frame(
    x   = stats::rnorm(k, 1, 2),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(k, 0.1 * data[["x"]], data[["sei"]])
  cauchy <- BayesTools::prior("cauchy", list(0, 0.5))
  fits <- lapply(c(standardized = TRUE, raw = FALSE), function(standardize) {
    suppressWarnings(brma(
      yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
      prior_effect = cauchy, prior_mods = list(x = cauchy),
      standardize_continuous_predictors = standardize,
      chains = 1, sample = 1000, burnin = 200, adapt = 100,
      seed = 1, silent = TRUE
    ))
  })
  fits[["convolution"]] <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    prior_mods = list(x = cauchy),
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  fits[["averaged"]] <- suppressWarnings(BMA.norm(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    prior_effect = cauchy, prior_mods = list(x = cauchy),
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  fits[["averaged_default"]] <- suppressWarnings(BMA.norm(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  .scalar_eligibility_cache[["fits"]] <- fits

  return(fits)
}


test_that("scalar point-test eligibility matches the hypotheses that run", {

  skip_on_cran()
  fits <- .scalar_eligibility_fits()

  for (fit_name in names(fits)) {
    fit        <- fits[[fit_name]]
    quantities <- hypothesis_quantities(fit)
    rows <- which(
      is.na(quantities[["bracket"]]) &
        quantities[["component"]] != "random" &
        !duplicated(quantities[["parameter"]])
    )
    for (i in rows) {
      parameter <- quantities[["parameter"]][[i]]
      listed    <- strsplit(
        quantities[["point_test_methods"]][[i]], ", ", fixed = TRUE
      )[[1L]]
      runs <- vapply(c("KDE", "qCMDE"), function(method) {
        result <- tryCatch(
          suppressWarnings(hypothesis(
            fit, paste0(parameter, " = 0.1"),
            density_method = method, n_samples = 2000, seed = 1
          )),
          error = function(error) error
        )
        !inherits(result, "error")
      }, logical(1))
      info <- paste(fit_name, parameter)
      expect_identical(unname(runs), c("KDE", "qCMDE") %in% listed, info = info)
      expect_identical(quantities[["point_test"]][[i]], any(runs), info = info)
    }
  }

  # A sum of Cauchy terms has an exact induced ordinate.
  standardized <- hypothesis_quantities(fits[["standardized"]])
  intercept    <- standardized[standardized[["alias"]] == "mu_intercept", , drop = FALSE]
  expect_true(intercept[["point_test"]])
  expect_identical(intercept[["point_test_methods"]], "KDE, qCMDE, IWMDE")
  expect_identical(intercept[["reason"]], "")
  expect_s3_class(suppressWarnings(hypothesis(
    fits[["standardized"]], "mu_intercept = 0", density_method = "KDE", seed = 1
  )), "data.frame")
  # Point hypotheses on the fitted-scale coefficient run as well.
  expect_s3_class(suppressWarnings(hypothesis(
    fits[["standardized"]], "mu_intercept = 0.1", density_method = "KDE",
    standardized_coefficients = TRUE, seed = 1
  )), "data.frame")

  # A normal intercept with a Cauchy slope has an exact induced ordinate.
  convolution <- hypothesis_quantities(fits[["convolution"]])
  convolution <- convolution[convolution[["alias"]] == "mu_intercept", , drop = FALSE]
  expect_true(convolution[["point_test"]])
  expect_identical(convolution[["point_test_methods"]], "KDE, qCMDE, IWMDE")
  expect_identical(convolution[["reason"]], "")
  expect_s3_class(suppressWarnings(hypothesis(
    fits[["convolution"]], "mu_intercept = 0", density_method = "KDE", seed = 1
  )), "data.frame")

  # The continuous ordinate of a model-averaged intercept (next to its atom at
  # 0) is exact with Cauchy alternatives and with the default normal priors.
  averaged <- hypothesis_quantities(fits[["averaged"]])
  averaged <- averaged[averaged[["alias"]] == "mu_intercept", , drop = FALSE]
  expect_true(averaged[["point_test"]])
  expect_identical(averaged[["point_test_methods"]], "KDE, qCMDE, IWMDE")
  expect_s3_class(suppressWarnings(hypothesis(
    fits[["averaged"]], "mu_intercept = 0.1", density_method = "KDE",
    standardized_coefficients = TRUE, seed = 1
  )), "data.frame")

  averaged_default <- hypothesis_quantities(fits[["averaged_default"]])
  averaged_default <- averaged_default[
    averaged_default[["alias"]] == "mu_intercept", , drop = FALSE
  ]
  expect_true(averaged_default[["point_test"]])
  expect_identical(averaged_default[["point_test_methods"]], "KDE, qCMDE, IWMDE")

  raw <- hypothesis_quantities(fits[["raw"]])
  expect_true(all(raw[["point_test"]]))
  expect_identical(unique(raw[["point_test_methods"]]), "KDE, qCMDE, IWMDE")
})
