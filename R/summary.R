
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
    conditional               = FALSE, ...) {

  ### model information
  is_mods       <- .is_mods(object)
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)
  is_bias       <- .is_bias(object)
  is_robma      <- .is_RoBMA(object)
  outcome_type  <- .outcome_type(object)

  BayesTools::check_bool(conditional, "conditional")
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
  if (length(common_parameters) > 0) {
    estimates_common <- BayesTools::JAGS_estimates_table(
      fit                = object[["fit"]],
      transform_factors  = TRUE,
      transform_scaled   = !standardized_coefficients,
      remove_diagnostics = !include_mcmc_diagnostics,
      remove_inclusion   = is_robma,
      keep_parameters    = common_parameters,
      probs              = probs,
      title              = if (is_mods || is_scale) "Common Estimates" else "Estimates"
    )
    if (conditional) {
      estimates_common_conditional <- BayesTools::JAGS_estimates_table(
        fit                = object[["fit"]],
        conditional        = TRUE,
        transform_factors  = TRUE,
        transform_scaled   = !standardized_coefficients,
        remove_diagnostics = !include_mcmc_diagnostics,
        remove_inclusion   = TRUE,
        keep_parameters    = common_parameters,
        probs              = probs,
        title              = if (is_mods || is_scale) {
          "Conditional Common Estimates"
        } else {
          "Conditional Estimates"
        }
      )
    } else {
      estimates_common_conditional <- list()
    }
  } else {
    estimates_common             <- list()
    estimates_common_conditional <- list()
  }

  ### provide regression estimates for the effect size meta-regression
  if (is_mods) {
    estimates_mods <- BayesTools::JAGS_estimates_table(
      fit                = object[["fit"]],
      transform_factors  = TRUE,
      transform_scaled   = !standardized_coefficients,
      remove_diagnostics = !include_mcmc_diagnostics,
      remove_inclusion   = is_robma,
      keep_formulas      = "mu",
      probs              = probs,
      formula_prefix     = FALSE,
      title              = if (is_scale) "Location" else "Meta-Regression"
    )
    if (conditional) {
      estimates_mods_conditional <- BayesTools::JAGS_estimates_table(
        fit                = object[["fit"]],
        conditional        = TRUE,
        transform_factors  = TRUE,
        transform_scaled   = !standardized_coefficients,
        remove_diagnostics = !include_mcmc_diagnostics,
        remove_inclusion   = TRUE,
        keep_formulas      = "mu",
        probs              = probs,
        formula_prefix     = FALSE,
        title              = if (is_scale) {
          "Conditional Location"
        } else {
          "Conditional Meta-Regression"
        }
      )
    } else {
      estimates_mods_conditional <- list()
    }
  } else {
    estimates_mods             <- list()
    estimates_mods_conditional <- list()
  }

  ### provide regression estimates for the scale meta-regression
  if (is_scale) {
    estimates_scale <- BayesTools::JAGS_estimates_table(
      fit                = object[["fit"]],
      transform_factors  = TRUE,
      transform_scaled   = !standardized_coefficients,
      remove_diagnostics = !include_mcmc_diagnostics,
      remove_inclusion   = is_robma,
      keep_formulas      = "log_tau",
      probs              = probs,
      formula_prefix     = FALSE,
      title              = "Scale",
      footnotes          = "exp(Intercept) corresponds to the between-study heterogeneity tau; the meta-regression coefficients correspond to the multiplicative effects on log-scale."
    )
    if (conditional) {
      estimates_scale_conditional <- BayesTools::JAGS_estimates_table(
        fit                = object[["fit"]],
        conditional        = TRUE,
        transform_factors  = TRUE,
        transform_scaled   = !standardized_coefficients,
        remove_diagnostics = !include_mcmc_diagnostics,
        remove_inclusion   = TRUE,
        keep_formulas      = "log_tau",
        probs              = probs,
        formula_prefix     = FALSE,
        title              = "Conditional Scale",
        footnotes          = "exp(Intercept) corresponds to the between-study heterogeneity tau; the meta-regression coefficients correspond to the multiplicative effects on log-scale."
      )
    } else {
      estimates_scale_conditional <- list()
    }
  } else {
    estimates_scale             <- list()
    estimates_scale_conditional <- list()
  }

  ### provide publication bias estimates
  if (is_bias) {
    estimates_bias <- BayesTools::JAGS_estimates_table(
      fit                = object[["fit"]],
      remove_diagnostics = !include_mcmc_diagnostics,
      remove_inclusion   = is_robma,
      keep_parameters    = c("bias", "omega", "PET", "PEESE"),
      probs              = probs,
      title              = "Publication Bias",
      footnotes          = if (.is_weightfunction(object)) "P-value intervals for publication bias weights omega correspond to one-sided p-values."
    )
    if (conditional) {
      estimates_bias_conditional <- BayesTools::JAGS_estimates_table(
        fit                = object[["fit"]],
        conditional        = TRUE,
        remove_diagnostics = !include_mcmc_diagnostics,
        remove_inclusion   = TRUE,
        keep_parameters    = c("bias", "omega", "PET", "PEESE"),
        probs              = probs,
        title              = "Conditional Publication Bias",
        footnotes          = if (.is_weightfunction(object)) "P-value intervals for publication bias weights omega correspond to one-sided p-values."
      )
    } else {
      estimates_bias_conditional <- list()
    }
  } else {
    estimates_bias             <- list()
    estimates_bias_conditional <- list()
  }

  ### provide RoBMA inclusion summaries
  if (is_robma) {
    inclusion <- .summary.RoBMA_inclusion_tables(object)
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

  model_name <- paste(model_name, "Model")
  if (is_multilevel) {
    model_name <- paste(model_name, sprintf("(k = %1$i, lvls = %2$i)", nrow(object[["data"]][["outcome"]]), length(unique(object[["data"]][["outcome"]][["cluster"]]))))
  } else {
    model_name <- paste(model_name, sprintf("(k = %1$i)", nrow(object[["data"]][["outcome"]])))
  }

  return(model_name)
}

# Split product-space inclusion inference into summary sections.
.summary.RoBMA_inclusion_tables <- function(object) {

  inclusion <- BayesTools::JAGS_inference_table(
    fit            = object[["fit"]],
    formula_prefix = TRUE
  )

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
  rownames(output) <- row_labels

  attr(output, "parameters") <- row_labels
  attr(output, "title")      <- title

  return(output)
}
