context("Model fitting for VIF parity")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("vif-parity")

test_that("brma and diagonal-V brma.mv fit equivalent meta-regressions", {

  dat <- data.frame(
    yi = c(-0.20, 0.05, 0.12, 0.31, -0.08, 0.22, 0.40, 0.14),
    vi = c(0.03, 0.05, 0.02, 0.08, 0.04, 0.06, 0.03, 0.07),
    x  = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
    z  = c(0, 1, 0, 2, 1, 0, 2, 1)
  )
  checks <- set_convergence_checks(
    max_Rhat = NULL,
    min_ESS  = NULL
  )

  fit_brma <- brma(
    yi                        = yi,
    vi                        = vi,
    mods                      = ~ x + z,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    chains                    = 2,
    sample                    = 1000,
    burnin                    = 1000,
    adapt                     = 1000,
    seed                      = 701,
    silent                    = TRUE,
    convergence_checks        = checks
  )
  fit_brma_mv <- brma.mv(
    yi                        = yi,
    V                         = vi,
    mods                      = ~ x + z,
    random                    = NULL,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    known_v_parameterization  = "block_mvn",
    chains                    = 2,
    sample                    = 1000,
    burnin                    = 1000,
    adapt                     = 1000,
    seed                      = 702,
    silent                    = TRUE,
    convergence_checks        = checks
  )
  fit_brma    <- add_marglik(fit_brma)
  fit_brma    <- suppressWarnings(add_loo(fit_brma))
  fit_brma_mv <- suppressWarnings(add_loo(fit_brma_mv))

  save_fit(
    "vif_parity_brma",
    fit_brma,
    info = list(data = dat)
  )
  save_fit(
    "vif_parity_brma_mv",
    fit_brma_mv,
    info = list(data = dat)
  )

  expect_s3_class(fit_brma, "brma.norm")
  expect_s3_class(fit_brma_mv, "brma.mv")
})
