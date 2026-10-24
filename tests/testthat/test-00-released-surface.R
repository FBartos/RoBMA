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
