# ============================================================================ #
# brma.dfbetas.R
# ============================================================================ #
#
# This file contains the dfbetas method for brma model fits. DFBETAS values
# measure the influence of each observation on the estimated model coefficients.
#
# Instead of computationally expensive refitting (leave-one-out), we use
# Pareto Smoothed Importance Sampling (PSIS) weights from the LOO-CV
# approximation to estimate the leave-one-out coefficients.
#
# ============================================================================ #


#' @title DFBETAS for brma Objects
#'
#' @description Computes DFBETAS (Difference in BETAS, standardized) for a
#' fitted brma object. DFBETAS measures the influence of each observation on
#' the estimated model coefficients. Positive values indicate that deleting the
#' observation yields a smaller estimate, negative values indicate that deleting
#' the observation yields a larger estimate.
#'
#' @param model a fitted brma object.
#' @param type type of parameters to be summarized. Defaults to \code{"mods"}
#' (for the effect size and meta-regression coefficients). The other options are
#' \code{"scale"} (for the heterogeneity and scale-regression coefficients) and
#' \code{"bias"} (for omega, PET, and PEESE publication-bias parameters).
#' @param standardized_coefficients whether to show standardized meta-regression coefficients.
#' Defaults to \code{FALSE}. When set to \code{TRUE}, standardized meta-regression
#' coefficients are returned for the intercept and continuous predictors. These coefficients
#' correspond to the standardized scale on which prior distributions are specified by default
#' (i.e., `standardize_continuous_predictors = TRUE`).
#' @param transform_factors whether to transform factors to their original names.
#' Defaults to \code{TRUE}.
#' @param return_loo_estimates whether to return the leave-one-out coefficient
#' estimates used to compute DFBETAS instead of standardized DFBETAS values.
#' Defaults to \code{FALSE}.
#' @param component optional parameter namespace. Use \code{"random"} for
#'   semantic random-effect SD, correlation, and allocation quantities.
#' @param parameter optional semantic random-effect quantity. When omitted with
#'   \code{component = "random"}, all available random quantities are returned.
#' @param ... additional arguments (currently ignored).
#'
#' @details
#' This function computes DFBETAS values using the Leave-One-Out (LOO)
#' approximation based on Pareto Smoothed Importance Sampling (PSIS) weights.
#' Ideally, DFBETAS is defined as:
#' \deqn{DFBETAS_{ij} = \frac{\hat{\beta}_j - \hat{\beta}_{j(-i)}}{SE(\hat{\beta}_{j(-i)})}}
#' where \eqn{\hat{\beta}_j} is the estimate of the \eqn{j}-th coefficient using
#' the full data, \eqn{\hat{\beta}_{j(-i)}} is the estimate when observation \eqn{i}
#' is omitted, and \eqn{SE(\hat{\beta}_{j(-i)})} is the standard error of the
#' coefficient when observation \eqn{i} is omitted.
#'
#' In the Bayesian context using LOO approximation:
#' \itemize{
#'   \item \eqn{\hat{\beta}_{j(-i)}} is estimated as the importance sampling
#'     weighted mean of the posterior samples, using PSIS weights \eqn{w_{is}}.
#'   \item \eqn{SE(\hat{\beta}_{j(-i)})} is estimated as the importance sampling
#'     weighted standard deviation of the posterior samples.
#' }
#'
#' This approximation allows computing influence statistics without refitting
#' the model \eqn{K} times, making it computationally efficient.
#' For \code{brma.mv()} objects, DFBETAS use estimate-unit PSIS weights. With
#' correlated known-\code{V} data, this is influence under conditional estimate
#' deletion, not independent-study deletion.
#' For \code{type = "bias"}, fixed identification parameters (e.g., the reference
#' \eqn{\omega = 1} interval) can have zero LOO posterior standard deviation.
#' These parameters are retained in the output, but their DFBETAS values are
#' reported as \code{NaN} because the standardized diagnostic is undefined.
#'
#' Note: This function requires that LOO-CV has been computed for the model
#' using \code{\link{add_loo}}.
#'
#' @return If `return_loo_estimates = FALSE`, a data frame with \eqn{K} rows
#' (observations) and \eqn{P} columns (parameters), containing DFBETAS values.
#' If `return_loo_estimates = TRUE`, returns the corresponding leave-one-out
#' coefficient estimates. Row names correspond to study labels (if available)
#' or indices.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit <- add_loo(fit)
#'
#'   inf <- dfbetas(fit)
#'   plot(inf[, 1], type = "h")
#' }
#' }
#'
#' @seealso \code{\link{add_loo}}, \code{\link{loo_weights.brma}}
#' @importFrom stats dfbetas
#' @exportS3Method
dfbetas.brma <- function(model, type = "mods", standardized_coefficients = FALSE,
                         transform_factors = TRUE,
                         return_loo_estimates = FALSE, component = NULL,
                         parameter = NULL, ...) {

  dots <- list(...)
  .weights <- dots[[".weights"]]
  BayesTools::check_char(type, "type", allow_values = c("mods", "scale", "bias"))
  BayesTools::check_char(component, "component", check_length = 1,
                         allow_NULL = TRUE)
  BayesTools::check_char(parameter, "parameter", check_length = 1,
                         allow_NULL = TRUE)
  component <- .diagnostic_parameter_component(
    type          = type,
    component     = component,
    type_supplied = !missing(type),
    allow_bias    = TRUE
  )
  if (!identical(component, "random") && !is.null(parameter)) {
    stop("'parameter' is currently available only with component = 'random'.",
         call. = FALSE)
  }

  if (is.null(.weights)) {
    psis_context <- .diagnostic_psis_context(model)
    .diagnostic_check_loo(model, context = psis_context, unit = "estimate")
    weights <- psis_context[["psis_weights"]]
  } else {
    weights <- .diagnostic_psis_weights(model, .weights)
  }

  # determine whether to extract formula (for meta-regression) or parameter (for intercept-only)
  is_scale <- .is_scale(model)
  is_bias  <- .is_bias(model)

  if (component == "mods") {
    samples_table <- .diagnostic_location_parameter_samples(
      model                     = model,
      standardized_coefficients = standardized_coefficients,
      transform_factors         = transform_factors
    )
  } else if (component == "scale") {
    if (is_scale) {
      # scale-regression: tau is modeled via log_tau formula
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                    = model[["fit"]],
        keep_formulas          = "log_tau",
        random_effects_summary = "none",
        remove_diagnostics     = TRUE,
        transform_factors      = transform_factors,
        return_samples         = TRUE
      )
    } else {
      # random/fixed effects: tau is a parameter (intercept)
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                    = model[["fit"]],
        keep_parameters        = "tau",
        random_effects_summary = "none",
        remove_diagnostics     = TRUE,
        transform_factors      = transform_factors,
        return_samples         = TRUE
      )
    }
  } else if (component == "bias") {
    if (!is_bias) {
      stop("type = 'bias' is only available for models with publication bias adjustment.", call. = FALSE)
    }

    samples_table <- .dfbetas_bias_samples(model)
  } else if (component == "random") {
    samples_table <- .diagnostic_random_parameter_samples(
      model                     = model,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients
    )
  }


  # samples_table is a flattened data frame (S rows x P columns)
  # ensure it's a matrix for computation
  samples_mat <- as.matrix(samples_table)

  if (ncol(samples_mat) == 0) {
    stop("No parameters available for DFBETAS with the requested type.", call. = FALSE)
  }

  result <- .dfbetas_internal(samples_mat, weights)

  if (return_loo_estimates) {
    beta_loo_df <- as.data.frame(result[["loo_estimates"]])
    colnames(beta_loo_df) <- colnames(samples_mat)
    beta_loo_df <- .diagnostic_set_rownames(beta_loo_df, model)
    return(beta_loo_df)
  }

  dfbetas_val <- result[["values"]]
  undefined   <- result[["undefined"]]
  note        <- NULL
  if (any(undefined)) {
    dfbetas_val[undefined] <- NaN
    note <- .diagnostic_zero_variance_note(
      diagnostic = "DFBETAS",
      parameters = colnames(samples_mat)[apply(undefined, 2, any)]
    )
  }
  # convert to data frame
  dfbetas_df <- as.data.frame(dfbetas_val)
  colnames(dfbetas_df) <- colnames(samples_mat)
  dfbetas_df <- .diagnostic_set_rownames(dfbetas_df, model)

  if (!is.null(note)) {
    attr(dfbetas_df, "note") <- note
  }
  class(dfbetas_df) <- c("dfbetas.brma", class(dfbetas_df))

  return(dfbetas_df)
}


# ---------------------------------------------------------------------------- #
# .dfbetas_internal
# ---------------------------------------------------------------------------- #
#
# Compute PSIS DFBETAS in affine-standardized posterior coordinates.
#
# ---------------------------------------------------------------------------- #
.dfbetas_internal <- function(samples, weights) {

  summary <- .psis_influence_summary(
    samples     = samples,
    weights     = weights,
    fit_moments = "all",
    variance    = "all"
  )
  K       <- nrow(summary[["loo_fit"]])
  P       <- ncol(summary[["loo_fit"]])
  se_loo  <- sqrt(summary[["loo_var"]])

  beta_full <- matrix(
    summary[["full_fit"]],
    nrow  = K,
    ncol  = P,
    byrow = TRUE
  )
  undefined <- se_loo == 0
  values    <- matrix(
    NaN,
    nrow     = K,
    ncol     = P,
    dimnames = dimnames(summary[["loo_fit"]])
  )
  values[!undefined] <-
    (beta_full[!undefined] - summary[["loo_fit"]][!undefined]) /
    se_loo[!undefined]

  loo_estimates <- sweep(summary[["loo_fit"]], 2, summary[["scale"]], "*")
  loo_estimates <- sweep(loo_estimates, 2, summary[["origin"]], "+")

  return(list(
    values        = values,
    loo_estimates = loo_estimates,
    undefined     = undefined
  ))
}


#' @exportS3Method
print.dfbetas.brma <- function(x, ...) {

  note <- attr(x, "note")
  class(x) <- "data.frame"
  print(x, ...)
  .print_diagnostic_note(note)

  return(invisible(x))
}


#' @title Convert DFBETAS Results to a Data Frame
#'
#' @description Converts a \code{dfbetas.brma} object to a component-aware
#' long data frame while retaining its displayed coefficient columns.
#'
#' @param x a \code{dfbetas.brma} object.
#' @param row.names \code{NULL} or a character vector giving the row names.
#' @param optional logical; passed to the final data-frame coercion.
#' @param stringsAsFactors accepted for compatibility with \code{data.frame()}.
#' @param ... unused additional arguments.
#'
#' @return A plain \code{data.frame} with leading \code{component} and
#' \code{parameter} columns.
#'
#' @export
as.data.frame.dfbetas.brma <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  output <- .output_table_as_long_data_frame(
    table            = x,
    component        = "dfbetas",
    row.names        = row.names,
    optional         = optional,
    stringsAsFactors = stringsAsFactors
  )

  return(output)
}


# ---------------------------------------------------------------------------- #
# .dfbetas_bias_samples
# ---------------------------------------------------------------------------- #
#
# Extract publication-bias posterior samples without using JAGS_estimates_table.
# BayesTools summary tables always run formula/factor post-processing, which can
# fail for RoBMA mixture-prior objects when we only need raw omega/PET/PEESE draws.
#
# ---------------------------------------------------------------------------- #
.dfbetas_bias_samples <- function(model) {

  posterior_samples <- .get_posterior_samples(model[["fit"]])
  bias_samples      <- list()
  omega_cols        <- grep("^omega(\\[|$)", colnames(posterior_samples), value = TRUE)

  if (length(omega_cols) > 0L) {
    bias_samples[["omega"]] <- .extract_omega_samples(posterior_samples)
  }

  for (par in c("PET", "PEESE")) {
    if (par %in% colnames(posterior_samples)) {
      bias_samples[[par]] <- as.matrix(posterior_samples[, par, drop = FALSE])
    }
  }

  if (length(bias_samples) == 0L) {
    return(matrix(nrow = nrow(posterior_samples), ncol = 0))
  }

  return(do.call(cbind, bias_samples))
}
