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

  if (!is.list(selection) || !inherits(selection, "BayesTools_parameter_selection")) {
    stop("'selection' must contain one resolved BayesTools parameter quantity.",
         call. = FALSE)
  }
  quantities <- selection[["quantities"]]
  if (
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

  keys <- quantities[["extraction_key"]]
  if (!is.list(keys) || length(keys) != 1L || !is.list(keys[[1L]]) ||
      !is.character(keys[[1L]][["type"]]) || length(keys[[1L]][["type"]]) != 1L ||
      is.na(keys[[1L]][["type"]]) || !nzchar(keys[[1L]][["type"]])) {
    stop("RoBMA parameter quantity has no valid extraction key.", call. = FALSE)
  }
  key <- keys[[1L]]
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
