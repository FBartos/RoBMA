context("Fully fixed point-prior models")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("fixed-null-models")

.fixed_model_settings <- function() {

  list(
    chains             = 2,
    sample             = 100,
    burnin             = 50,
    adapt              = 100,
    seed               = 1,
    silent             = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
}

.add_fixed_model_diagnostics <- function(fit) {

  fit <- add_marglik(fit)
  fit <- suppressWarnings(add_loo(fit))

  return(fit)
}

test_that("all single-model classes fit fixed effect and heterogeneity priors", {

  settings      <- .fixed_model_settings()
  yi            <- c(-0.05, 0.08, 0.12)
  sei           <- c(0.12, 0.10, 0.15)
  effect_null   <- BayesTools::prior("spike", parameters = list(location = 0))
  tau_null      <- BayesTools::prior("spike", parameters = list(location = 0))
  effect_fixed  <- BayesTools::prior("spike", parameters = list(location = 0.20))
  tau_fixed     <- BayesTools::prior("spike", parameters = list(location = 0.15))

  fit_null <- brma.norm(
    yi                        = yi,
    sei                       = sei,
    measure                   = "GEN",
    prior_effect              = effect_null,
    prior_heterogeneity       = tau_null,
    prior_unit_information_sd = 1,
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  fit_null <- .add_fixed_model_diagnostics(fit_null)
  expect_s3_class(fit_null, "brma.norm")
  save_fit(
    "fixed_null_brma",
    fit_null,
    info = list(yi = yi, sei = sei, mu = 0, tau = 0)
  )

  fit_fixed <- brma.norm(
    yi                        = yi,
    sei                       = sei,
    measure                   = "GEN",
    prior_effect              = effect_fixed,
    prior_heterogeneity       = tau_fixed,
    prior_unit_information_sd = 1,
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  fit_fixed <- .add_fixed_model_diagnostics(fit_fixed)
  expect_s3_class(fit_fixed, "brma.norm")
  save_fit(
    "fixed_nonzero_brma",
    fit_fixed,
    info = list(yi = yi, sei = sei, mu = 0.20, tau = 0.15)
  )

  V <- matrix(
    c(
      0.0144, 0.0030, 0.0000,
      0.0030, 0.0100, 0.0020,
      0.0000, 0.0020, 0.0225
    ),
    nrow  = 3,
    byrow = TRUE
  )
  fit_mv <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = NULL,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_effect              = effect_fixed,
    prior_heterogeneity       = tau_fixed,
    prior_unit_information_sd = 1,
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  fit_mv <- .add_fixed_model_diagnostics(fit_mv)
  expect_s3_class(fit_mv, "brma.mv")
  save_fit(
    "fixed_nonzero_brma_mv",
    fit_mv,
    info = list(yi = yi, V = V, mu = 0.20, tau = 0.15)
  )

  fit_glmm <- brma.glmm(
    ai                        = c(3, 4, 5),
    n1i                       = c(10, 12, 14),
    ci                        = c(2, 3, 4),
    n2i                       = c(10, 12, 14),
    measure                   = "OR",
    prior_effect              = effect_null,
    prior_heterogeneity       = tau_null,
    prior_unit_information_sd = 1,
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  fit_glmm <- .add_fixed_model_diagnostics(fit_glmm)
  expect_s3_class(fit_glmm, "brma.glmm")
  save_fit(
    "fixed_null_brma_glmm",
    fit_glmm,
    info = list(mu = 0, tau = 0)
  )

  fit_PET <- bPET(
    yi                        = yi,
    sei                       = sei,
    measure                   = "GEN",
    prior_effect              = effect_null,
    prior_heterogeneity       = tau_null,
    prior_bias                = BayesTools::prior_PET(
      "spike",
      parameters = list(location = 0.30)
    ),
    prior_unit_information_sd = 1,
    effect_direction          = "positive",
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  fit_PET <- .add_fixed_model_diagnostics(fit_PET)
  expect_s3_class(fit_PET, "bPET")
  save_fit(
    "fixed_null_bPET",
    fit_PET,
    info = list(yi = yi, sei = sei, mu = 0, tau = 0, PET = 0.30)
  )

  fit_PEESE <- bPEESE(
    yi                        = yi,
    sei                       = sei,
    measure                   = "GEN",
    prior_effect              = effect_null,
    prior_heterogeneity       = tau_null,
    prior_bias                = BayesTools::prior_PEESE(
      "spike",
      parameters = list(location = 0.40)
    ),
    prior_unit_information_sd = 1,
    effect_direction          = "positive",
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  fit_PEESE <- .add_fixed_model_diagnostics(fit_PEESE)
  expect_s3_class(fit_PEESE, "bPEESE")
  save_fit(
    "fixed_null_bPEESE",
    fit_PEESE,
    info = list(yi = yi, sei = sei, mu = 0, tau = 0, PEESE = 0.40)
  )

  fit_selection <- bselmodel(
    yi                        = yi,
    sei                       = sei,
    measure                   = "GEN",
    steps                     = 0.05,
    prior_effect              = effect_null,
    prior_heterogeneity       = tau_null,
    prior_bias                = BayesTools::prior_weightfunction(
      "one-sided",
      0.05,
      BayesTools::wf_fixed(c(1, 0.80))
    ),
    prior_unit_information_sd = 1,
    effect_direction          = "positive",
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  fit_selection <- .add_fixed_model_diagnostics(fit_selection)
  expect_s3_class(fit_selection, "bselmodel")
  save_fit(
    "fixed_null_bselmodel",
    fit_selection,
    info = list(yi = yi, sei = sei, mu = 0, tau = 0)
  )
})
