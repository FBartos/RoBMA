context("nobs and coef methods")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

test_that("nobs matches outcome data and metafor", {

  metafor_names <- list_fits(has_metafor = TRUE)
  skip_if_missing_fits(metafor_names)

  for (name in names(fits)) {
    fit <- fits[[name]]
    expect_true(is.numeric(nobs(fit)) && nobs(fit) > 0,
                info = paste0("nobs is positive for '", name, "'"))
    expect_equal(nobs(fit), nrow(fit$data$outcome),
                 info = paste0("nobs matches outcome rows for '", name, "'"))
  }

  for (name in metafor_names) {
    expect_equal(
      nobs(fits[[name]]),
      as.integer(nobs(info[[name]][["metafor"]])),
      info = paste0("nobs matches metafor for '", name, "'")
    )
  }
})

test_that("coef returns model coefficients across cached fits", {

  for (name in names(fits)) {
    fit <- fits[[name]]
    cf  <- coef(fit)
    expect_true(is.numeric(cf) && length(cf) > 0,
                info = paste0("coef is non-empty numeric for '", name, "'"))
    expect_true(!is.null(names(cf)),
                info = paste0("coef is named for '", name, "'"))
    expect_identical(cf, fit$coefficients,
                     info = paste0("coef matches $coefficients for '", name, "'"))
  }
})
