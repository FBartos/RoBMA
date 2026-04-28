context("Autocompute options")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metafor")

test_that("autocompute options work for brma()", {

  # Using a very simple model to be fast
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  # subset to 3 studies for speed
  dat <- dat[1:3, ]

  # Ensure options are off by default
  RoBMA.options(autocompute.loo = FALSE, autocompute.waic = FALSE, autocompute.marglik = FALSE)
  fit_default <- suppressWarnings(brma(yi = yi, vi = vi, data = dat, measure = "RR",
                                       seed = 1, silent = TRUE, chains = 3, sample = 500, burnin = 100))

  expect_null(fit_default$loo)
  expect_null(fit_default$waic)
  expect_null(fit_default$marglik)


  # Turn on options
  RoBMA.options(autocompute.loo = TRUE, autocompute.waic = TRUE, autocompute.marglik = TRUE)
  on.exit(RoBMA.options(autocompute.loo = FALSE, autocompute.waic = FALSE, autocompute.marglik = FALSE))

  fit_auto <- suppressWarnings(brma(yi = yi, vi = vi, data = dat, measure = "RR",
                                    seed = 1, silent = TRUE, chains = 3, sample = 500, burnin = 100))

  expect_true(!is.null(fit_auto$loo))
  expect_true(!is.null(fit_auto$waic))
  expect_true(!is.null(fit_auto$marglik))

  expect_s3_class(fit_auto$loo$estimate, "loo")
  expect_s3_class(fit_auto$waic$estimate, "waic")
  expect_s3_class(fit_auto$marglik, "bridge")
})
