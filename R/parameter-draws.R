# ============================================================================ #
# parameter-draws.R
# ============================================================================ #
#
# Deferred extraction for quantities resolved through the fitted parameter
# catalog. BayesTools owns its quantities; RoBMA evaluates only provider-owned
# extraction keys.
#
# ============================================================================ #


#' @exportS3Method BayesTools::parameter_draws
#' @noRd
parameter_draws.brma <- function(object, selection, ...) {

  quantities <- selection[["quantities"]]
  if (!inherits(selection, "BayesTools_parameter_selection") ||
      !is.data.frame(quantities) || nrow(quantities) != 1L) {
    stop(
      "'selection' must contain one resolved BayesTools parameter quantity.",
      call. = FALSE
    )
  }

  provider <- quantities[["provider"]]
  if (identical(provider, "BayesTools")) {
    return(BayesTools::parameter_draws(object[["fit"]], selection, ...))
  }
  if (!identical(provider, "RoBMA")) {
    stop(
      "Parameter quantity provider '", provider, "' is not supported.",
      call. = FALSE
    )
  }

  key <- quantities[["extraction_key"]][[1L]]
  if (identical(key[["type"]], "robma_formula_group")) {
    return(BayesTools::JAGS_materialize_draws(
      fit        = object[["fit"]],
      parameters = key[["dependencies"]]
    ))
  }
  if (identical(key[["type"]], "robma_bias_component")) {
    stop(
      "Publication-bias components do not define one generic draw quantity; ",
      "use the corresponding RoBMA plot or summary method.",
      call. = FALSE
    )
  }

  stop(
    "Unsupported RoBMA parameter extraction key '", key[["type"]], "'.",
    call. = FALSE
  )
}
