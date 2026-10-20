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
parameter_draws.brma <- function(object, selection, model_samples = NULL, ...) {

  if (!is.list(selection) || !inherits(selection, "BayesTools_parameter_selection")) {
    stop("'selection' must contain one resolved BayesTools parameter quantity.",
         call. = FALSE)
  }
  quantities <- selection[["quantities"]]
  if (!is.data.frame(quantities) || nrow(quantities) != 1L) {
    stop(
      "'selection' must contain one resolved BayesTools parameter quantity.",
      call. = FALSE
    )
  }

  provider <- quantities[["provider"]]
  if (identical(provider, "BayesTools")) {
    return(BayesTools::parameter_draws(
      object[["fit"]], selection, model_samples = model_samples, ...
    ))
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
  .check_unused_dots(list(...), allowed = character(),
                     caller = "parameter_draws.brma()")
  catalog <- .brma_parameter_catalog_metadata(object)[["catalog"]]
  current <- catalog[["quantities"]]
  current <- current[
    current[["quantity_id"]] %in% quantities[["quantity_id"]], , drop = FALSE
  ]
  if (nrow(current) != 1L) {
    stop("'selection' does not identify a current RoBMA parameter quantity.",
         call. = FALSE)
  }
  resolved <- BayesTools::parameter_catalog_resolve(
    catalog,
    alias     = current[["canonical_name"]],
    namespace = current[["namespace"]],
    component = current[["component"]]
  )
  rownames(current) <- NULL
  rownames(quantities) <- NULL
  rownames(resolved[["quantities"]]) <- NULL
  wrapper_fields <- c("schema_version", "parameter_map_version", "quantity_id")
  if (!identical(names(selection), names(resolved)) ||
      !all(vapply(wrapper_fields, function(field) {
        identical(selection[[field]], resolved[[field]])
      }, logical(1))) || !identical(quantities, resolved[["quantities"]])) {
    stop("'selection' does not identify a current RoBMA parameter quantity.",
         call. = FALSE)
  }
  key <- current[["extraction_key"]][[1L]]
  if (identical(key[["type"]], "robma_formula_group")) {
    fit <- object[["fit"]]
    if (!is.null(model_samples)) {
      if (!is.matrix(model_samples) || !is.numeric(model_samples) ||
          is.null(colnames(model_samples)) ||
          anyDuplicated(colnames(model_samples))) {
        stop("'model_samples' must be a numeric matrix with unique column names.",
             call. = FALSE)
      }
      missing <- setdiff(key[["dependencies"]], colnames(model_samples))
      if (length(missing) > 0L) {
        stop(
          "The supplied 'model_samples' are missing declared source coordinates: ",
          paste(missing, collapse = ", "), ".",
          call. = FALSE
        )
      }
      source <- if (inherits(model_samples, "mcmc")) {
        model_samples
      } else {
        coda::mcmc(model_samples)
      }
      fit <- BayesTools::JAGS_with_draws(fit, source)
    }
    return(BayesTools::JAGS_materialize_draws(
      fit        = fit,
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
