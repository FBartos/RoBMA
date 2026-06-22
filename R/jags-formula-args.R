# ============================================================================ #
# jags-formula-args.R
# ============================================================================ #
#
# Shared construction of BayesTools formula arguments for JAGS fitting and
# bridge sampling.
#
# ============================================================================ #


.create_fit_scale_formula_prior_list <- function(priors, parameter) {

  if (inherits(priors[["scale"]], "RoBMA_scale_priors")) {
    return(priors[["scale"]][[parameter]])
  }

  priors[["scale"]]
}


.create_jags_formula_args <- function(data, priors) {

  formula_args <- list(
    formula_list                        = list(),
    formula_data_list                   = list(),
    formula_prior_list                  = list(),
    formula_scale_list                  = list(),
    formula_random_prior_list           = list(),
    formula_random_effects_compile_list = list(),
    add_parameters                      = character()
  )

  ### add effect regressions and random effects
  if (.is_data_mods(data) || .is_data_random(data)) {
    location_parameter <- if (.is_data_random(data)) "location" else "mods"

    formula_args[["formula_list"]][["mu"]]       <- .create_fit_formula_list(data = data, parameter = location_parameter)
    formula_args[["formula_data_list"]][["mu"]]  <- .create_fit_formula_data_list(data = data, parameter = location_parameter)
    formula_args[["formula_prior_list"]][["mu"]] <- .create_fit_formula_prior_list(priors = priors, parameter = location_parameter)
    formula_args[["formula_scale_list"]][["mu"]] <- .data_standardize_continuous_predictors(data)
    if (.is_data_random(data)) {
      formula_args[["formula_random_prior_list"]][["mu"]] <- priors[["random"]]
      random_effects_compile <- .data_random_effects_compile(data)
      if (!is.null(random_effects_compile)) {
        formula_args[["formula_random_effects_compile_list"]][["mu"]] <- random_effects_compile
      }
    }
  }

  ### add heterogeneity regressions
  if (.is_data_scale(data)) {
    for (scale_spec in .data_scale_component_specs(data)) {
      parameter <- scale_spec[["parameter"]]
      formula_args[["formula_list"]][[parameter]]       <- .create_fit_scale_formula(scale_spec[["formula"]])
      formula_args[["formula_data_list"]][[parameter]]  <- scale_spec[["data"]]
      formula_args[["formula_prior_list"]][[parameter]] <- .create_fit_scale_formula_prior_list(
        priors    = priors,
        parameter = parameter
      )
      formula_args[["formula_scale_list"]][[parameter]] <- .data_standardize_continuous_predictors(data)
    }
    if (.is_data_random(data)) {
      formula_args[["add_parameters"]] <- unique(c(
        formula_args[["add_parameters"]],
        .data_scale_formula_sources(data)
      ))
    }
  }

  return(formula_args)
}


.optional_jags_list <- function(x) {

  if (length(x) > 0L) {
    return(x)
  }

  NULL
}


.optional_jags_character <- function(x) {

  if (length(x) > 0L) {
    return(x)
  }

  NULL
}
