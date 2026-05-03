context("Autocompute options")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")

test_that("autocompute options are applied for brma()", {

  old_options <- list(
    autocompute.loo     = RoBMA.get_option("autocompute.loo"),
    autocompute.waic    = RoBMA.get_option("autocompute.waic"),
    autocompute.marglik = RoBMA.get_option("autocompute.marglik")
  )
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  # Using a very simple model to keep the test-01 fit cheap.
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )
  dat <- dat[1:3, ]

  RoBMA.options(
    autocompute.loo     = FALSE,
    autocompute.waic    = FALSE,
    autocompute.marglik = FALSE
  )
  fit_default <- suppressWarnings(brma(
    yi      = yi,
    vi      = vi,
    data    = dat,
    measure = "RR",
    seed    = 1,
    silent  = TRUE,
    chains  = 2,
    sample  = 250,
    burnin  = 50,
    adapt   = 100
  ))

  expect_null(fit_default$loo)
  expect_null(fit_default$waic)
  expect_null(fit_default$marglik)

  RoBMA.options(
    autocompute.loo     = TRUE,
    autocompute.waic    = TRUE,
    autocompute.marglik = TRUE
  )
  fit_auto <- suppressWarnings(brma(
    yi      = yi,
    vi      = vi,
    data    = dat,
    measure = "RR",
    seed    = 1,
    silent  = TRUE,
    chains  = 2,
    sample  = 250,
    burnin  = 50,
    adapt   = 100
  ))

  expect_s3_class(fit_auto$loo$estimate, "loo")
  expect_s3_class(fit_auto$waic$estimate, "waic")
  expect_s3_class(fit_auto$marglik, "bridge")
})
