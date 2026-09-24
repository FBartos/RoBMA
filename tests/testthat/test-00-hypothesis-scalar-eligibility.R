context("Point-test eligibility of scalar formula coefficients")

# hypothesis_quantities() must list a point-test method for a scalar formula
# coefficient exactly when hypothesis() runs a point hypothesis on it with
# that method. The same Cauchy slope prior is fitted with a standardized and
# with a raw predictor: the standardized fit's original-scale intercept
# combines the intercept with the Cauchy slope (no exact prior ordinate),
# while the raw fit's intercept is the fitted coefficient itself.
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
  fits <- lapply(c(standardized = TRUE, raw = FALSE), function(standardize) {
    suppressWarnings(brma(
      yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
      prior_mods = list(x = BayesTools::prior("cauchy", list(0, 0.5))),
      standardize_continuous_predictors = standardize,
      chains = 1, sample = 1000, burnin = 200, adapt = 100,
      seed = 1, silent = TRUE
    ))
  })
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

  standardized <- hypothesis_quantities(fits[["standardized"]])
  intercept    <- standardized[standardized[["alias"]] == "mu_intercept", , drop = FALSE]
  expect_false(intercept[["point_test"]])
  expect_identical(intercept[["point_test_methods"]], "")
  expect_identical(intercept[["reason"]], paste0(
    "Point hypotheses are not supported for 'mu_intercept': its induced prior ",
    "on the original scale has no exact ordinate. Point hypotheses on the ",
    "fitted-scale coefficient are available with ",
    "standardized_coefficients = TRUE."
  ))
  expect_error(
    suppressWarnings(hypothesis(
      fits[["standardized"]], "mu_intercept = 0", density_method = "KDE"
    )),
    "is not exact enough for a point-null Bayes factor",
    fixed = TRUE
  )
  # The alternative named in the reason runs.
  expect_s3_class(suppressWarnings(hypothesis(
    fits[["standardized"]], "mu_intercept = 0.1", density_method = "KDE",
    standardized_coefficients = TRUE, seed = 1
  )), "data.frame")

  raw <- hypothesis_quantities(fits[["raw"]])
  expect_true(all(raw[["point_test"]]))
  expect_identical(unique(raw[["point_test_methods"]]), "KDE, qCMDE, IWMDE")
})
