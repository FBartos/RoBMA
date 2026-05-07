context("CRAN smoke test")

on_cran <- get("on_cran", envir = asNamespace("testthat"), inherits = FALSE)

skip_if_not(on_cran(), "CRAN smoke test only")

test_that("package loads and fits a minimal brma model", {

  expect_true("RoBMA" %in% loadedNamespaces())

  yi  <- c(-0.12, 0.05, 0.18, 0.31)
  sei <- c(0.20, 0.18, 0.22, 0.19)

  fit <- brma(
    yi = yi, sei = sei, measure = "GEN",
    prior_unit_information_sd = 1,
    sample = 100, burnin = 50, adapt = 500,
    chains = 1, seed = 1, silent = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )

  fit_summary <- summary(fit)

  expect_s3_class(fit, "brma.norm")
  expect_s3_class(fit_summary, "summary.brma")
})
