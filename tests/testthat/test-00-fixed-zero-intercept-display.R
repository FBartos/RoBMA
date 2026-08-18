context("Fixed-zero meta-regression intercept display")


.fixed_zero_intercept_test_data <- data.frame(
  yi    = c(-0.20, 0.05, 0.18, 0.31, 0.42),
  sei   = c(0.20, 0.18, 0.22, 0.19, 0.21),
  mod   = c(-1, -0.5, 0, 0.5, 1),
  study = factor(c(1, 1, 2, 2, 3))
)

.fixed_zero_intercept_point <- BayesTools::prior(
  "spike",
  parameters = list(0)
)

.fixed_zero_intercept_object <- function(mods, prior_effect = NULL) {

  args <- list(
    yi                        = .fixed_zero_intercept_test_data[["yi"]],
    sei                       = .fixed_zero_intercept_test_data[["sei"]],
    mods                      = mods,
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  if (!is.null(prior_effect)) {
    args[["prior_effect"]] <- prior_effect
  }

  do.call(brma.norm, args)
}


test_that("only fixed-zero intercepts accompanying moderators are omitted", {

  no_intercept <- .fixed_zero_intercept_object(~ 0 + mod)
  fixed_zero   <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  estimated    <- .fixed_zero_intercept_object(~ mod)
  intercept_only <- brma.norm(
    yi                        = yi,
    sei                       = sei,
    prior_effect              = .fixed_zero_intercept_point,
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_model <- brma.mv(
    yi                        = yi,
    V                         = sei^2,
    mods                      = ~ 0 + mod,
    random                    = ~ (1 | study),
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_true(.location_omit_fixed_zero_intercept(no_intercept))
  expect_true(.location_omit_fixed_zero_intercept(fixed_zero))
  expect_true(.location_omit_fixed_zero_intercept(random_model))
  expect_false(.location_omit_fixed_zero_intercept(estimated))
  expect_false(.location_omit_fixed_zero_intercept(intercept_only))
})


test_that("public parameter discovery follows fixed-zero intercept display", {

  fixed_zero <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  estimated <- .fixed_zero_intercept_object(~ mod)

  fixed_catalog     <- .brma_parameter_catalog_unfitted(fixed_zero)
  estimated_catalog <- .brma_parameter_catalog_unfitted(estimated)

  expect_false("mu_intercept" %in% fixed_catalog[["parameter"]])
  expect_true("mu_mod" %in% fixed_catalog[["parameter"]])
  expect_true("mu_intercept" %in% estimated_catalog[["parameter"]])
})


test_that("marginal means drop the intercept from structural metadata", {

  fixed_zero <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  estimated <- .fixed_zero_intercept_object(~ mod)
  formula   <- ~ mod
  terms     <- c("intercept", "mod")
  parameters <- c("mu_intercept", "mu_mod")

  fixed_setup <- .marginal_means_drop_fixed_zero_intercept(
    object     = fixed_zero,
    formula    = formula,
    terms      = terms,
    parameters = parameters
  )
  estimated_setup <- .marginal_means_drop_fixed_zero_intercept(
    object     = estimated,
    formula    = formula,
    terms      = terms,
    parameters = parameters
  )

  expect_identical(fixed_setup[["terms"]], "mod")
  expect_identical(fixed_setup[["parameters"]], "mu_mod")
  expect_identical(attr(stats::terms(fixed_setup[["formula"]]), "intercept"), 0L)
  expect_identical(estimated_setup[["terms"]], terms)
  expect_identical(estimated_setup[["parameters"]], parameters)
})


test_that("summary routes omit the fixed-zero location intercept", {

  object <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  object[["fit"]] <- structure(
    list(),
    class      = "BayesTools_fit",
    prior_list = list()
  )
  observed <- NULL

  testthat::local_mocked_bindings(
    .location_omit_fixed_zero_intercept = function(object) TRUE,
    .summary_estimates_pair = function(
        enabled, object, probs, include_mcmc_diagnostics, is_robma,
        conditional, main_args, conditional_args = main_args) {

      if (identical(main_args[["title"]], "Meta-Regression")) {
        observed <<- list(
          main        = main_args,
          conditional = conditional_args
        )
      }
      list(estimates = list(), conditional = list())
    },
    .package = "RoBMA"
  )

  summary.brma(object, include_mcmc_diagnostics = FALSE)

  expect_identical(observed[["main"]][["remove_parameters"]], "mu_intercept")
  expect_identical(
    observed[["conditional"]][["remove_parameters"]],
    "mu_intercept"
  )
})


test_that("stored summaries omit the fixed-zero location intercept", {

  observed <- NULL
  testthat::local_mocked_bindings(
    .location_omit_fixed_zero_intercept = function(object) TRUE,
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    JAGS_estimates_table = function(...) {

      observed <<- list(...)
      data.frame(Mean = 0, row.names = "mu_mod")
    },
    .package = "BayesTools"
  )

  .object_summary(list(fit = structure(list(), class = "BayesTools_fit")))

  expect_true("mu_intercept" %in% observed[["remove_parameters"]])
  expect_true(all(c("theta", "gamma", "sampling_z", "pi", "phi") %in%
                    observed[["remove_parameters"]]))
})
