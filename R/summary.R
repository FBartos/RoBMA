
#' @title Summarize brma Object
#'
#' @description \code{summary.brma} creates summary tables for a
#' brma object. For RoBMA objects, inclusion summaries are printed before
#' parameter estimates.
#'
#' @param object a fitted brma object
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .50, .975)}
#' @param include_mcmc_diagnostics whether to include MCMC diagnostics in the output.
#' Defaults to \code{TRUE}.
#' @param standardized_coefficients whether to show standardized meta-regression coefficients.
#' Defaults to \code{FALSE}. When set to \code{TRUE}, standardized meta-regression
#' coefficients are returned for the intercept and continuous predictors. These coefficients
#' correspond to the standardized scale on which prior distributions are specified by default
#' (i.e., `standardize_continuous_predictors = TRUE`).
#' @param conditional whether to include conditional estimates for RoBMA
#'   product-space objects. Defaults to \code{FALSE}.
#' @param logBF whether to show inclusion Bayes factors on the log scale.
#' Defaults to \code{FALSE}.
#' @param BF01 whether to show inverse inclusion Bayes factors. Defaults to
#' \code{FALSE}.
#' @param ... additional arguments
#'
#' @return A list of class `summary.brma` with model name, optional RoBMA
#' inclusion tables, common estimates, moderator estimates, scale estimates,
#' publication-bias estimates, and optional conditional estimates. The printed
#' form displays the non-empty tables.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- bPET(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary(fit)
#' }
#' }
#'
#'
#' @seealso [brma()], [brma.glmm()]
#' @export
summary.brma       <- function(
    object, probs = c(.025, .50, .975),
    include_mcmc_diagnostics  = TRUE,
    standardized_coefficients = FALSE,
    conditional               = FALSE,
    logBF                     = FALSE,
    BF01                      = FALSE, ...) {

  ### model information
  is_mods       <- .is_mods(object)
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)
  is_bias       <- .is_bias(object)
  is_robma      <- .is_RoBMA(object)
  outcome_type  <- .outcome_type(object)

  BayesTools::check_bool(include_mcmc_diagnostics, "include_mcmc_diagnostics")
  BayesTools::check_bool(standardized_coefficients, "standardized_coefficients")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(BF01, "BF01")
  if (conditional && !is_robma) {
    stop("'conditional' summaries are available only for RoBMA objects.",
         call. = FALSE)
  }

  ### special cases handling
  # deal with `only_data` fit
  if (is.null(object[["priors"]]) && is.null(object[["fit"]])) {
    return(object[["data"]])
  }

  ### provide common estimates
  common_parameters <- intersect(
    c("mu", "tau", "rho"),
    names(attr(object[["fit"]], "prior_list"))
  )
  estimates_common_pair <- .summary_estimates_pair(
    enabled                  = length(common_parameters) > 0L,
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      transform_factors = TRUE,
      transform_scaled  = !standardized_coefficients,
      keep_parameters   = common_parameters,
      title             = if (is_mods || is_scale) "Common Estimates" else "Estimates"
    ),
    conditional_args         = list(
      transform_factors = TRUE,
      transform_scaled  = !standardized_coefficients,
      keep_parameters   = common_parameters,
      title             = if (is_mods || is_scale) {
        "Conditional Common Estimates"
      } else {
        "Conditional Estimates"
      }
    )
  )
  estimates_common             <- estimates_common_pair[["estimates"]]
  estimates_common_conditional <- estimates_common_pair[["conditional"]]

  ### provide regression estimates for the effect size meta-regression
  estimates_mods_pair <- .summary_estimates_pair(
    enabled                  = is_mods,
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      transform_factors = TRUE,
      transform_scaled  = !standardized_coefficients,
      keep_formulas     = "mu",
      formula_prefix    = FALSE,
      title             = if (is_scale) "Location" else "Meta-Regression"
    ),
    conditional_args         = list(
      transform_factors = TRUE,
      transform_scaled  = !standardized_coefficients,
      keep_formulas     = "mu",
      formula_prefix    = FALSE,
      title             = if (is_scale) {
        "Conditional Location"
      } else {
        "Conditional Meta-Regression"
      }
    )
  )
  estimates_mods             <- estimates_mods_pair[["estimates"]]
  estimates_mods_conditional <- estimates_mods_pair[["conditional"]]

  ### provide regression estimates for the scale meta-regression
  scale_footnotes <- "exp(Intercept) corresponds to the between-study heterogeneity tau; the meta-regression coefficients correspond to the multiplicative effects on log-scale."
  estimates_scale_pair <- .summary_estimates_pair(
    enabled                  = is_scale,
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      transform_factors = TRUE,
      transform_scaled  = !standardized_coefficients,
      keep_formulas     = "log_tau",
      formula_prefix    = FALSE,
      title             = "Scale",
      footnotes         = scale_footnotes
    ),
    conditional_args         = list(
      transform_factors = TRUE,
      transform_scaled  = !standardized_coefficients,
      keep_formulas     = "log_tau",
      formula_prefix    = FALSE,
      title             = "Conditional Scale",
      footnotes         = scale_footnotes
    )
  )
  estimates_scale             <- estimates_scale_pair[["estimates"]]
  estimates_scale_conditional <- estimates_scale_pair[["conditional"]]

  ### provide publication bias estimates
  bias_footnotes <- if (.is_weightfunction(object)) {
    "P-value intervals for publication bias weights omega correspond to one-sided p-values."
  }
  estimates_bias_pair <- .summary_estimates_pair(
    enabled                  = is_bias,
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      keep_parameters = c("bias", "omega", "PET", "PEESE"),
      title           = "Publication Bias",
      footnotes       = bias_footnotes
    ),
    conditional_args         = list(
      keep_parameters = c("bias", "omega", "PET", "PEESE"),
      title           = "Conditional Publication Bias",
      footnotes       = bias_footnotes
    )
  )
  estimates_bias             <- estimates_bias_pair[["estimates"]]
  estimates_bias_conditional <- estimates_bias_pair[["conditional"]]

  ### provide RoBMA inclusion summaries
  if (is_robma) {
    inclusion <- .summary.RoBMA_inclusion_tables(
      object                   = object,
      include_mcmc_diagnostics = include_mcmc_diagnostics,
      logBF                    = logBF,
      BF01                     = BF01
    )
  } else {
    inclusion <- list(
      inclusion_components = list(),
      inclusion_mods       = list(),
      inclusion_scale      = list()
    )
  }

  out <- list(
    name                        = .summary.brma_model_names(object),
    inclusion_components        = inclusion[["inclusion_components"]],
    inclusion_mods              = inclusion[["inclusion_mods"]],
    inclusion_scale             = inclusion[["inclusion_scale"]],
    estimates                   = estimates_common,
    estimates_conditional       = estimates_common_conditional,
    estimates_mods              = estimates_mods,
    estimates_mods_conditional  = estimates_mods_conditional,
    estimates_scale             = estimates_scale,
    estimates_scale_conditional = estimates_scale_conditional,
    estimates_bias              = estimates_bias,
    estimates_bias_conditional  = estimates_bias_conditional
  )

  class(out) <- "summary.brma"
  attr(out, "mods")         <- is_mods
  attr(out, "scale")        <- is_scale
  attr(out, "multilevel")   <- is_multilevel
  attr(out, "bias")         <- is_bias
  attr(out, "RoBMA")        <- is_robma
  attr(out, "outcome_type") <- outcome_type

  return(out)
}

#' @rdname summary.brma
#' @param x a `summary.brma` or fitted `brma` object.
#' @export
print.summary.brma <- function(x, ...) {

  cat("\n")
  cat(x[["name"]])
  cat("\n")

  for (type in c(
    "inclusion_components", "inclusion_mods", "inclusion_scale",
    "estimates", "estimates_conditional",
    "estimates_mods", "estimates_mods_conditional",
    "estimates_scale", "estimates_scale_conditional",
    "estimates_bias", "estimates_bias_conditional"
  )) {
    if (length(x[[type]]) > 0) {
      cat("\n")
      print(x[[type]])
    }
  }


  cat("\n")

  return(invisible(x))
}

#' @rdname summary.brma
#' @export
print.brma <- function(x, ...) {
  print(summary(x, ...))
}

.summary_estimates_table <- function(object, probs, include_mcmc_diagnostics,
                                     is_robma, conditional = FALSE, ...) {

  args <- list(
    fit                = object[["fit"]],
    conditional        = conditional,
    remove_diagnostics = !include_mcmc_diagnostics,
    remove_inclusion   = if (conditional) TRUE else is_robma,
    probs              = probs,
    diagnostic_columns = .summary_estimates_diagnostic_columns(
      include_mcmc_diagnostics
    )
  )
  args <- c(args, list(...))

  return(do.call(BayesTools::JAGS_estimates_table, args))
}

.summary_estimates_pair <- function(enabled, object, probs,
                                    include_mcmc_diagnostics, is_robma,
                                    conditional, main_args,
                                    conditional_args = main_args) {

  if (!enabled) {
    return(list(estimates = list(), conditional = list()))
  }

  shared_args <- list(
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma
  )

  estimates <- do.call(
    .summary_estimates_table,
    c(shared_args, main_args)
  )

  if (conditional) {
    estimates_conditional <- do.call(
      .summary_estimates_table,
      c(shared_args, list(conditional = TRUE), conditional_args)
    )
  } else {
    estimates_conditional <- list()
  }

  return(list(
    estimates   = estimates,
    conditional = estimates_conditional
  ))
}

.summary.brma_model_names <- function(object) {

  is_mods       <- .is_mods(object)
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)

  if (.is_robust_RoBMA(object)) {
    model_name <- "Robust Bayesian Model-Averaged"
  } else if (.is_BMA(object)) {
    model_name <- "Bayesian Model-Averaged"
  } else {
    model_name <- "Bayesian"
  }

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

  n_effects <- nrow(object[["data"]][["outcome"]])

  model_name <- paste(model_name, "Model")
  if (is_multilevel) {
    n_clusters <- length(unique(object[["data"]][["outcome"]][["cluster"]]))
    model_name <- paste(model_name, sprintf("(k = %1$i, clusters = %2$i)", n_effects, n_clusters))
  } else {
    model_name <- paste(model_name, sprintf("(k = %1$i)", n_effects))
  }

  return(model_name)
}

# Split product-space inclusion inference into summary sections.
.summary.RoBMA_inclusion_tables <- function(object, include_mcmc_diagnostics,
                                            logBF, BF01) {

  args <- list(
    fit                   = object[["fit"]],
    formula_prefix        = TRUE,
    logBF                 = logBF,
    BF01                  = BF01,
    BF_diagnostic_columns = .summary_BF_diagnostic_columns(
      include_mcmc_diagnostics
    )
  )

  inclusion <- do.call(BayesTools::JAGS_inference_table, args)

  parameters <- attr(inclusion, "parameters")
  row_labels <- rownames(inclusion)

  core_map <- c(
    mu                = "Effect",
    mu_intercept      = "Effect",
    tau               = "Heterogeneity",
    log_tau_intercept = "Heterogeneity",
    bias              = "Publication Bias"
  )
  core_parameters <- names(core_map)[names(core_map) %in% parameters]
  core_indices    <- match(core_parameters, parameters)

  mods_indices <- grep("^mu_", parameters)
  mods_indices <- mods_indices[parameters[mods_indices] != "mu_intercept"]

  scale_indices <- grep("^log_tau_", parameters)
  scale_indices <- scale_indices[parameters[scale_indices] != "log_tau_intercept"]

  output <- list(
    inclusion_components = .summary.inclusion_subtable(
      table      = inclusion,
      indices    = core_indices,
      row_labels = unname(core_map[core_parameters]),
      title      = "Component Inclusion"
    ),
    inclusion_mods = .summary.inclusion_subtable(
      table      = inclusion,
      indices    = mods_indices,
      row_labels = .summary_parameter_label(
        sub("^\\(mu\\) ", "", row_labels[mods_indices])
      ),
      title      = if (.is_scale(object)) {
        "Location Inclusion"
      } else {
        "Meta-Regression Inclusion"
      }
    ),
    inclusion_scale = .summary.inclusion_subtable(
      table      = inclusion,
      indices    = scale_indices,
      row_labels = .summary_parameter_label(
        sub("^\\(log_tau\\) ", "", row_labels[scale_indices])
      ),
      title      = "Scale Inclusion"
    )
  )

  return(output)
}

.summary_estimates_diagnostic_columns <- function(include_mcmc_diagnostics) {

  if (isTRUE(include_mcmc_diagnostics)) {
    return(c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  }

  return("none")
}

.summary_BF_diagnostic_columns <- function(include_mcmc_diagnostics) {

  if (isTRUE(include_mcmc_diagnostics)) {
    return("BF_error_percent")
  }

  return("none")
}

# Convert internal interaction separators back to formula syntax.
.summary_parameter_label <- function(label) {

  return(gsub("__xXx__", ":", label, fixed = TRUE))
}

# Create a labelled BayesTools inclusion subtable.
.summary.inclusion_subtable <- function(table, indices, row_labels, title) {

  if (length(indices) == 0) {
    return(list())
  }

  output <- table[indices, , drop = FALSE]
  output <- .summary.inclusion_subtable_restore_BF_attributes(
    output  = output,
    table   = table,
    indices = indices
  )
  rownames(output) <- row_labels

  attr(output, "parameters") <- row_labels
  attr(output, "title")      <- title

  return(output)
}

# Restore row-level BayesTools_BF attributes after table subsetting.
.summary.inclusion_subtable_restore_BF_attributes <- function(output, table,
                                                              indices) {

  for (column in intersect(colnames(output), colnames(table))) {
    if (inherits(table[[column]], "BayesTools_BF")) {
      output[[column]] <- table[[column]][indices]
    }
  }

  return(output)
}
