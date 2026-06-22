
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
#' @param standardized_coefficients whether to show standardized formula
#' coefficients. Defaults to \code{FALSE}. When set to \code{TRUE},
#' meta-regression coefficients and random-component summaries are returned on
#' the standardized scale on which prior distributions are specified by default
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
  is_random     <- .is_random(object)
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
      transform_factors      = TRUE,
      transform_scaled       = !standardized_coefficients,
      keep_parameters        = common_parameters,
      random_effects_summary = "none",
      title                  = if (is_mods || is_scale) "Common Estimates" else "Estimates"
    ),
    conditional_args         = list(
      transform_factors      = TRUE,
      transform_scaled       = !standardized_coefficients,
      keep_parameters        = common_parameters,
      random_effects_summary = "none",
      title                  = if (is_mods || is_scale) {
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
    enabled                  = is_mods || is_random,
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      transform_factors      = TRUE,
      transform_scaled       = !standardized_coefficients,
      keep_formulas          = "mu",
      random_effects_summary = "none",
      formula_prefix         = FALSE,
      title                  = if (is_scale || is_random) "Location" else "Meta-Regression"
    ),
    conditional_args         = list(
      transform_factors      = TRUE,
      transform_scaled       = !standardized_coefficients,
      keep_formulas          = "mu",
      random_effects_summary = "none",
      formula_prefix         = FALSE,
      title                  = if (is_scale || is_random) {
        "Conditional Location"
      } else {
        "Conditional Meta-Regression"
      }
    )
  )
  estimates_mods             <- estimates_mods_pair[["estimates"]]
  estimates_mods_conditional <- estimates_mods_pair[["conditional"]]

  ### provide regression estimates for the scale meta-regression
  scale_footnotes <- .summary_scale_footnotes(object)
  scale_formulas       <- .summary_scale_formula_parameters(object)
  scale_formula_prefix <- length(scale_formulas) > 1L
  estimates_scale_pair <- .summary_estimates_pair(
    enabled                  = is_scale,
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      transform_factors      = TRUE,
      transform_scaled       = !standardized_coefficients,
      keep_formulas          = scale_formulas,
      random_effects_summary = "none",
      formula_prefix         = scale_formula_prefix,
      title                  = "Scale",
      footnotes              = scale_footnotes
    ),
    conditional_args         = list(
      transform_factors      = TRUE,
      transform_scaled       = !standardized_coefficients,
      keep_formulas          = scale_formulas,
      random_effects_summary = "none",
      formula_prefix         = scale_formula_prefix,
      title                  = "Conditional Scale",
      footnotes              = scale_footnotes
    )
  )
  estimates_scale             <- .summary_scale_repair_row_labels(
    estimates = estimates_scale_pair[["estimates"]],
    object    = object
  )
  estimates_scale_conditional <- .summary_scale_repair_row_labels(
    estimates = estimates_scale_pair[["conditional"]],
    object    = object
  )

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
      keep_parameters        = c("bias", "omega", "PET", "PEESE"),
      random_effects_summary = "none",
      title                  = "Publication Bias",
      footnotes              = bias_footnotes
    ),
    conditional_args         = list(
      keep_parameters        = c("bias", "omega", "PET", "PEESE"),
      random_effects_summary = "none",
      title                  = "Conditional Publication Bias",
      footnotes              = bias_footnotes
    )
  )
  estimates_bias             <- estimates_bias_pair[["estimates"]]
  estimates_bias_conditional <- estimates_bias_pair[["conditional"]]

  ### provide random-effect component estimates
  estimates_random_pair <- .summary_estimates_pair(
    enabled                  = .summary_random_components_enabled(object),
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      keep_parameters         = "random_effects",
      transform_scaled        = !standardized_coefficients,
      random_effects_summary  = "standard",
      random_effects_metadata = TRUE,
      formula_prefix          = FALSE,
      title                   = "Random"
    ),
    conditional_args         = list(
      keep_parameters         = "random_effects",
      transform_scaled        = !standardized_coefficients,
      random_effects_summary  = "standard",
      random_effects_metadata = TRUE,
      formula_prefix          = FALSE,
      title                   = "Conditional Random"
    )
  )
  estimates_random             <- estimates_random_pair[["estimates"]]
  estimates_random_conditional <- estimates_random_pair[["conditional"]]

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
    known_v_backend             = .brma_mv_known_v_backend_metadata(object),
    inclusion_components        = inclusion[["inclusion_components"]],
    inclusion_mods              = inclusion[["inclusion_mods"]],
    inclusion_scale             = inclusion[["inclusion_scale"]],
    estimates                   = estimates_common,
    estimates_conditional       = estimates_common_conditional,
    estimates_mods              = estimates_mods,
    estimates_mods_conditional  = estimates_mods_conditional,
    estimates_scale             = estimates_scale,
    estimates_scale_conditional = estimates_scale_conditional,
    estimates_random             = estimates_random,
    estimates_random_conditional = estimates_random_conditional,
    estimates_bias              = estimates_bias,
    estimates_bias_conditional  = estimates_bias_conditional
  )

  class(out) <- "summary.brma"
  attr(out, "mods")         <- is_mods
  attr(out, "scale")        <- is_scale
  attr(out, "random")       <- is_random
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

  x_print <- .summary_brma_prepare_print_sections(x)

  cat("\n")
  cat(x_print[["name"]])
  cat("\n")

  for (type in c(
    "inclusion_components", "inclusion_mods", "inclusion_scale",
    "estimates", "estimates_conditional",
    "estimates_mods", "estimates_mods_conditional",
    "estimates_scale", "estimates_scale_conditional",
    "estimates_random", "estimates_random_conditional",
    "estimates_bias", "estimates_bias_conditional"
  )) {
    if (length(x_print[[type]]) > 0) {
      cat("\n")
      print(x_print[[type]])
    }
  }


  cat("\n")

  return(invisible(x))
}

.summary_brma_prepare_print_sections <- function(x) {

  x[["estimates_random"]] <- .summary_brma_random_section_for_print(
    random = x[["estimates_random"]],
    title  = "Random"
  )
  x[["estimates_random_conditional"]] <- .summary_brma_random_section_for_print(
    random = x[["estimates_random_conditional"]],
    title  = "Conditional Random"
  )

  return(x)
}

.summary_brma_random_section_for_print <- function(random, title) {

  if (length(random) == 0L) {
    return(random)
  }

  out <- .summary_brma_random_estimates_for_print(random)
  if (length(out) == 0L) {
    return(out)
  }
  attr(out, "title") <- title
  attr(out, "rownames") <- attr(random, "rownames", exact = TRUE)

  return(out)
}

.summary_brma_random_estimates_for_print <- function(random) {

  metadata_cols <- c("Random name", "Random grouping", "Random structure")
  estimate_cols <- setdiff(colnames(random), metadata_cols)
  if (length(estimate_cols) == 0L) {
    return(list())
  }

  random_names <- if ("Random name" %in% colnames(random)) {
    as.character(random[["Random name"]])
  } else {
    rep("", nrow(random))
  }

  out <- random[, estimate_cols, drop = FALSE]
  rownames(out) <- .summary_brma_random_print_labels(
    labels       = rownames(random),
    random_names = random_names
  )

  return(out)
}

.summary_brma_random_print_labels <- function(labels, random_names) {

  has_name <- !is.na(random_names) & nzchar(random_names)
  if (length(unique(random_names[has_name])) <= 1L) {
    return(labels)
  }

  for (i in seq_along(labels)) {
    if (!has_name[[i]]) {
      next
    }
    label <- labels[[i]]
    label <- sub("^sd\\((.*) \\| .*\\)$", "sd(\\1)", label)
    label <- sub("^rho\\(.*\\)$", "rho", label)
    labels[[i]] <- paste0(random_names[[i]], ": ", label)
  }

  return(labels)
}

#' @rdname summary.brma
#' @export
print.brma <- function(x, ...) {
  print(summary(x, ...))
}

.summary_scale_formula_parameters <- function(object) {

  parameters <- .data_scale_formula_parameters(object[["data"]])
  if (length(parameters) == 0L) {
    return("log_tau")
  }

  parameters
}

.summary_scale_footnotes <- function(object) {

  if (inherits(object, "brma.mv") && .is_random(object)) {
    return(
      "exp(Intercept) corresponds to the targeted random-effect SD; the meta-regression coefficients correspond to multiplicative effects on log-scale."
    )
  }

  "exp(Intercept) corresponds to the between-study heterogeneity tau; the meta-regression coefficients correspond to the multiplicative effects on log-scale."
}

.summary_scale_repair_row_labels <- function(estimates, object) {

  if (length(estimates) == 0L || is.null(rownames(estimates))) {
    return(estimates)
  }

  scale_names <- .summary_scale_display_names(object)
  for (parameter in names(scale_names)) {
    rownames(estimates) <- sub(
      pattern     = paste0("^\\(", parameter, "\\)"),
      replacement = paste0("(", scale_names[[parameter]], ")"),
      x           = rownames(estimates)
    )
  }
  rownames(estimates) <- sub(
    pattern     = "\\) intercept$",
    replacement = ") exp(intercept)",
    x           = rownames(estimates)
  )

  estimates
}

.summary_scale_display_names <- function(object) {

  scale_specs <- .data_scale_component_specs(object[["data"]])
  if (length(scale_specs) == 0L) {
    return(character(0))
  }

  out <- vapply(scale_specs, `[[`, character(1), "display_name")
  stats::setNames(out, vapply(scale_specs, `[[`, character(1), "parameter"))
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

.summary_random_components_enabled <- function(object) {

  if (!.is_random(object)) {
    return(FALSE)
  }

  design <- .fitted_formula_design(object, "mu", required = FALSE)
  if (is.null(design) || length(design[["random_effects"]]) == 0L) {
    return(FALSE)
  }
  if (!.is_scale(object)) {
    return(TRUE)
  }

  any(vapply(design[["random_effects"]], function(term) {
    binding <- term[["sd_binding"]]
    if (is.null(binding)) {
      return(TRUE)
    }
    length(binding[["factors_by_column"]]) > 0L ||
      length(binding[["allocations"]]) > 0L ||
      is.null(binding[["source"]])
  }, logical(1)))
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

  if (inherits(object, "brma.mv")) {
    model_name <- paste(model_name, "Multivariate")
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
