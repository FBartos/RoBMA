context("nobs and coef")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)


# ============================================================================ #
# Test: nobs
# ============================================================================ #

test_that("nobs returns correct number of observations", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  expect_equal(nobs(fit_brma), nrow(fit_brma$data$outcome),
               info = "nobs should match number of rows in outcome data")

  expect_equal(nobs(fit_brma), as.integer(nobs(fit_metafor)),
               info = "nobs should match metafor")
})

test_that("nobs works across model types", {

  for (name in names(fits)) {
    fit <- fits[[name]]
    expect_true(is.numeric(nobs(fit)) && nobs(fit) > 0,
                info = paste0("nobs should be a positive number for '", name, "'"))
  }
})

# ============================================================================ #
# Test: coef
# ============================================================================ #

test_that("coef returns named numeric vector", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  cf <- coef(fit_brma)

  expect_true(is.numeric(cf),
              info = "coef should return numeric vector")
  expect_true(!is.null(names(cf)),
              info = "coef should return named vector")
  expect_identical(cf, fit_brma$coefficients,
                   info = "coef should match $coefficients")
})

test_that("coef works across model types", {

  for (name in names(fits)) {
    fit <- fits[[name]]
    cf  <- coef(fit)
    expect_true(is.numeric(cf) && length(cf) > 0,
                info = paste0("coef should return non-empty numeric for '", name, "'"))
  }
})
