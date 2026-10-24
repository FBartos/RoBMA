# Argument order, defaults, and entry points of RoBMA 4.0.0 that later
# development must keep. New arguments follow the released ones.

.released_arguments_norm <- c(
  "yi", "vi", "sei", "weights", "ni", "mods", "scale", "cluster", "data",
  "slab", "subset", "measure"
)

test_that("selection arguments follow the released RoBMA() and bselmodel() arguments", {

  released_RoBMA <- c(
    .released_arguments_norm, "effect_direction",
    "prior_effect", "prior_heterogeneity", "prior_mods", "prior_scale",
    "prior_heterogeneity_allocation", "prior_bias",
    "prior_effect_null", "prior_heterogeneity_null", "prior_mods_null",
    "prior_scale_null", "prior_heterogeneity_allocation_null", "prior_bias_null",
    "standardize_continuous_predictors", "set_contrast_factor_predictors",
    "prior_unit_information_sd", "rescale_priors",
    "prior_informed_field", "prior_informed_subfield", "model_type",
    "sample", "burnin", "adapt", "chains", "thin", "parallel", "autofit",
    "autofit_control", "convergence_checks", "seed", "silent"
  )
  released_bselmodel <- c(
    .released_arguments_norm,
    "prior_effect", "prior_heterogeneity", "prior_mods", "prior_scale",
    "prior_heterogeneity_allocation", "prior_bias",
    "standardize_continuous_predictors", "set_contrast_factor_predictors",
    "prior_unit_information_sd", "rescale_priors",
    "prior_informed_field", "prior_informed_subfield",
    "effect_direction", "steps",
    "sample", "burnin", "adapt", "chains", "thin", "parallel", "autofit",
    "autofit_control", "convergence_checks", "seed", "silent"
  )

  for (case in list(
    list(name = "RoBMA",     fun = RoBMA,     released = released_RoBMA),
    list(name = "bselmodel", fun = bselmodel, released = released_bselmodel)
  )) {
    arguments <- names(formals(case[["fun"]]))
    expect_identical(arguments, c(case[["released"]], "selection", "selection_control", "..."),
                     info = case[["name"]])
  }
})

test_that("parameter selectors of plot(), plot_prior(), and print_prior() have no default", {

  no_default <- function(fun, name) identical(formals(fun)[[name]], quote(expr = ))

  for (name in c("parameter", "parameter_mods", "parameter_scale")) {
    expect_true(no_default(plot.brma, name), info = paste("plot.brma", name))
    expect_true(no_default(print_prior.brma, name), info = paste("print_prior.brma", name))
  }
  expect_identical(formals(plot_prior.brma)[["parameter"]], "mu")
  expect_true(no_default(plot_prior.brma, "parameter_mods"))
  expect_true(no_default(plot_prior.brma, "parameter_scale"))

  prior_object <- RoBMA(
    yi = c(.1, .2, .3), sei = c(.1, .1, .1), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  expect_no_error(utils::capture.output(print_prior(prior_object)))
  expect_no_error(utils::capture.output(print_prior(prior_object, parameter = "tau")))
  expect_no_error(plot_prior(prior_object))
})
