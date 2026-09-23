context("Prior inputs for absent model components")

yi  <- c(-0.2, 0.1, 0.3, 0.0)
sei <- c(0.2, 0.3, 0.25, 0.2)

normal_prior <- BayesTools::prior(
  distribution = "normal",
  parameters   = list(mean = 0, sd = 1)
)
allocation_prior <- BayesTools::prior(
  distribution = "beta",
  parameters   = list(alpha = 2, beta = 2)
)
point_prior <- BayesTools::prior(
  distribution = "spike",
  parameters   = list(location = 0)
)

test_that("brma rejects explicit priors for absent components", {

  cases <- list(
    prior_mods                       = normal_prior,
    prior_scale                      = normal_prior,
    prior_heterogeneity_allocation   = allocation_prior
  )
  required_components <- c(
    prior_mods                     = "mods",
    prior_scale                    = "scale",
    prior_heterogeneity_allocation = "cluster"
  )

  for (argument in names(cases)) {
    call_args <- list(
      yi          = yi,
      sei         = sei,
      measure     = "SMD",
      only_priors = TRUE
    )
    call_args[[argument]] <- cases[[argument]]

    expect_error(
      do.call(brma, call_args),
      regexp = paste0("'", argument, "'.*'", required_components[[argument]], "'")
    )
  }

  expect_s3_class(
    brma(yi = yi, sei = sei, measure = "SMD", only_priors = TRUE),
    "brma"
  )
})

test_that("BMA rejects alternative and null priors for absent components", {

  cases <- list(
    prior_mods                           = normal_prior,
    prior_mods_null                      = NULL,
    prior_scale                          = normal_prior,
    prior_scale_null                     = NULL,
    prior_heterogeneity_allocation       = allocation_prior,
    prior_heterogeneity_allocation_null  = NULL
  )
  required_components <- c(
    prior_mods                          = "mods",
    prior_mods_null                     = "mods",
    prior_scale                         = "scale",
    prior_scale_null                    = "scale",
    prior_heterogeneity_allocation      = "cluster",
    prior_heterogeneity_allocation_null = "cluster"
  )

  for (argument in names(cases)) {
    call_args <- list(
      yi          = yi,
      sei         = sei,
      measure     = "SMD",
      only_priors = TRUE
    )
    call_args[argument] <- cases[argument]

    expect_error(
      do.call(BMA, call_args),
      regexp = paste0("'", argument, "'.*'", required_components[[argument]], "'")
    )
  }

  expect_s3_class(
    BMA(yi = yi, sei = sei, measure = "SMD", only_priors = TRUE),
    "BMA.norm"
  )
})

test_that("explicit priors remain valid for present components", {

  model_data <- data.frame(
    yi      = yi,
    sei     = sei,
    x       = c(-1, 0, 1, 2),
    z       = c(0, 1, 0, 1),
    cluster = c("a", "a", "b", "b")
  )

  fit <- BMA(
    yi                                    = yi,
    sei                                   = sei,
    mods                                  = ~ x,
    scale                                 = ~ z,
    cluster                               = cluster,
    data                                  = model_data,
    measure                               = "SMD",
    prior_mods                            = normal_prior,
    prior_mods_null                       = point_prior,
    prior_scale                           = normal_prior,
    prior_scale_null                      = point_prior,
    prior_heterogeneity_null              = NULL,
    prior_heterogeneity_allocation        = allocation_prior,
    prior_heterogeneity_allocation_null   = point_prior,
    only_priors                           = TRUE
  )

  expect_s3_class(fit, "BMA.norm")
  expect_false(is.null(fit[["priors"]][["mods"]]))
  expect_false(is.null(fit[["priors"]][["scale"]]))
  expect_true("rho" %in% names(fit[["priors"]][["outcome"]]))
})
