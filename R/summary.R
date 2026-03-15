
#' @title Summarize brma Object
#'
#' @description \code{summary.brma} creates summary tables for a
#' brma object.
#'
#' @param object a fitted brma object
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .50, .975)}
#' @param include_MCMC_diagnostics whether to include MCMC diagnostics in the output.
#' Defaults to \code{TRUE}.
#' @param standardized_coefficients whether to show standardized meta-regression coefficients.
#' Defaults to \code{FALSE}. When set to \code{TRUE}, standardized meta-regression
#' coefficients are returned for the intercept and continuous predictors. These coefficients
#' correspond to the standardized scale on which prior distributions are specified by default
#' (i.e., `standardize_continuous_predictors = TRUE`).
#' @param ... additional arguments
#'
#' @examples \dontrun{
#' }
#'
#'
#' @seealso [brma.norm()], [brma.glmm()]
#' @export
summary.brma       <- function(
    object, probs = c(.025, .50, .975),
    include_MCMC_diagnostics = TRUE, standardized_coefficients = FALSE, ...) {

  ### model information
  is_mods       <- .is_mods(object)
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)
  is_bias       <- .is_bias(object)
  outcome_type  <- .outcome_type(object)

  ### special cases handling
  # deal with `only_data` fit
  if (is.null(object[["priors"]]) && is.null(object[["fit"]])) {
    return(object[["data"]])
  }

  ### provide common estimates
  estimates_common <- BayesTools::JAGS_estimates_table(
    fit                = object[["fit"]],
    transform_factors  = TRUE,
    transform_scaled   = !standardized_coefficients,
    remove_diagnostics = !include_MCMC_diagnostics,
    keep_parameters    = c("mu", "tau", "rho"),
    probs              = probs,
    title              = if (is_mods || is_scale) "Model Results (Common Estimates):" else "Model Results:"
  )

  if (nrow(estimates_common) == 0) {
    estimates_common <- list()
  }

  ### provide regression estimates for the effect size meta-regression
  if (is_mods) {
    estimates_mods <- BayesTools::JAGS_estimates_table(
      fit                = object[["fit"]],
      transform_factors  = TRUE,
      transform_scaled   = !standardized_coefficients,
      remove_diagnostics = !include_MCMC_diagnostics,
      keep_formulas      = "mu",
      probs              = probs,
      formula_prefix     = FALSE,
      title              = if (is_scale) "Model Results (Location):" else "Model Results (Meta-Regression):"
    )
  } else {
    estimates_mods <- list()
  }

  ### provide regression estimates for the scale meta-regression
  if (is_scale) {
    estimates_scale <- BayesTools::JAGS_estimates_table(
      fit                = object[["fit"]],
      transform_factors  = TRUE,
      transform_scaled   = !standardized_coefficients,
      remove_diagnostics = !include_MCMC_diagnostics,
      keep_formulas      = "log_tau",
      probs              = probs,
      formula_prefix     = FALSE,
      title              = "Model Results (Scale):",
      footnotes          = "exp(Intercept) corresponds to the between-study heterogeneity tau; the meta-regression coefficients correspond to the multiplicative effects on log-scale."
    )
  } else {
    estimates_scale <- list()
  }

  ### provide publication bias estimates
  if (is_bias) {
    estimates_bias <- BayesTools::JAGS_estimates_table(
      fit                = object[["fit"]],
      remove_diagnostics = !include_MCMC_diagnostics,
      keep_parameters    = c("bias", "omega", "PET", "PEESE"),
      probs              = probs,
      title              = "Model Results (Publication Bias):",
      footnotes          = if (.is_weightfunction(object)) "P-value intervals for publication bias weights omega correspond to one-sided p-values."
    )
  } else {
    estimates_bias <- list()
  }

  out <- list(
    name            = .summary.brma_model_names(object),
    estimates       = estimates_common,
    estimates_mods  = estimates_mods,
    estimates_scale = estimates_scale,
    estimates_bias  = estimates_bias
  )

  class(out) <- "summary.brma"
  attr(out, "mods")         <- is_mods
  attr(out, "scale")        <- is_scale
  attr(out, "multilevel")   <- is_multilevel
  attr(out, "bias")         <- is_bias
  attr(out, "outcome_type") <- outcome_type

  return(out)
}

#' @export
print.summary.brma <- function(x, ...) {

  cat("\n")
  cat(x[["name"]])
  cat("\n")

  for (type in c("estimates", "estimates_mods", "estimates_scale", "estimates_bias")) {
    if (length(x[[type]]) > 0) {
      cat("\n")
      print(x[[type]])
    }
  }


  cat("\n")

  return(invisible(x))
}

#' @export
print.brma <- function(x, ...) {
  print(summary(x, ...))
}

.summary.brma_model_names <- function(object) {

  is_mods       <- .is_mods(object)
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)

  model_name <- "Bayesian"

  if (is_multilevel) {
    model_name <- paste(model_name, "Multilevel")
  }

  if (is_scale) {
    model_name <- paste(model_name, "Location-Scale")
  } else if (is_mods) {
    model_name <- paste(model_name, "Mixed-Effect")
  } else if (!is_mods && !is_scale) {
    model_name <- paste(model_name, "Random-Effects")
  }

  if ("bselmodel" %in% class(object)) {
    model_name <- paste(model_name, "Selection")
  } else if ("bPET" %in% class(object)) {
    model_name <- paste(model_name, "PET")
  } else if ("bPEESE" %in% class(object)) {
    model_name <- paste(model_name, "PEESE")
  }

  model_name <- paste(model_name, "Model")
  if (is_multilevel) {
    model_name <- paste(model_name, sprintf("(k = %1$i, lvls = %2$i)", nrow(object[["data"]][["outcome"]]), length(unique(object[["data"]][["outcome"]][["study_ids"]]))))
  } else {
    model_name <- paste(model_name, sprintf("(k = %1$i)", nrow(object[["data"]][["outcome"]])))
  }

  return(model_name)
}
