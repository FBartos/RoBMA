context("Fixed point-prior monitoring")

.fixed_prior_monitors <- function(object) {

  fit_priors <- .create_fit_priors(
    data   = object[["data"]],
    priors = object[["priors"]]
  )

  return(BayesTools::JAGS_to_monitor(fit_priors))
}

test_that("single-model constructors monitor fixed mu and tau", {

  yi          <- c(0.10, 0.20)
  sei         <- c(0.10, 0.12)
  effect      <- BayesTools::prior("spike", parameters = list(location = 0))
  tau_nonzero <- BayesTools::prior("spike", parameters = list(location = 0.15))
  common      <- list(
    yi                        = yi,
    sei                       = sei,
    measure                   = "GEN",
    prior_effect              = effect,
    prior_heterogeneity       = tau_nonzero,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  objects <- list(
    brma.norm = do.call(brma.norm, common),
    brma.mv = brma.mv(
      yi                        = yi,
      V                         = matrix(c(0.01, 0.002, 0.002, 0.0144), nrow = 2),
      random                    = NULL,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_effect              = effect,
      prior_heterogeneity       = tau_nonzero,
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    brma.glmm = brma.glmm(
      ai                        = c(3, 4),
      n1i                       = c(10, 12),
      ci                        = c(2, 3),
      n2i                       = c(10, 12),
      measure                   = "OR",
      prior_effect              = effect,
      prior_heterogeneity       = tau_nonzero,
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    bPET = do.call(bPET, c(
      common,
      list(prior_bias = BayesTools::prior_PET(
        "spike",
        parameters = list(location = 0.30)
      ))
    )),
    bPEESE = do.call(bPEESE, c(
      common,
      list(prior_bias = BayesTools::prior_PEESE(
        "spike",
        parameters = list(location = 0.40)
      ))
    )),
    bselmodel = do.call(bselmodel, c(
      common,
      list(
        steps      = 0.05,
        prior_bias = BayesTools::prior_weightfunction(
          "one-sided",
          0.05,
          BayesTools::wf_fixed(c(1, 0.80))
        )
      )
    ))
  )

  for (name in names(objects)) {
    monitors <- .fixed_prior_monitors(objects[[name]])
    expect_true("mu" %in% monitors, info = name)
    expect_true("tau" %in% monitors, info = name)
  }

  expect_true("PET" %in% .fixed_prior_monitors(objects[["bPET"]]))
  expect_true("PEESE" %in% .fixed_prior_monitors(objects[["bPEESE"]]))
  expect_true("omega" %in% .fixed_prior_monitors(objects[["bselmodel"]]))
})

test_that("product-space constructors retain mu and tau monitors", {

  BMA_norm <- BMA.norm(
    yi                        = c(0.10, 0.20),
    sei                       = c(0.10, 0.12),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  BMA_glmm <- BMA.glmm(
    ai                        = c(3, 4),
    n1i                       = c(10, 12),
    ci                        = c(2, 3),
    n2i                       = c(10, 12),
    measure                   = "OR",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  RoBMA_fit <- RoBMA(
    yi                        = c(0.10, 0.20),
    sei                       = c(0.10, 0.12),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  for (object in list(BMA_norm, BMA_glmm, RoBMA_fit)) {
    monitors <- .fixed_prior_monitors(object)
    expect_true("mu" %in% monitors)
    expect_true("tau" %in% monitors)
  }
})
