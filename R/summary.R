
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
#' form displays the non-empty tables. In meta-regressions with moderators, a
#' location intercept fixed at zero is omitted; intercept-only models retain it.
#' The random table reports the quantities aligned with prior specification;
#' use [summary_heterogeneity()] for aggregate variances and the complete family
#' of deterministic allocation transforms.
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
  location_remove_parameters <- if (
    .location_omit_fixed_zero_intercept(object)
  ) {
    "mu_intercept"
  } else {
    NULL
  }
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
      remove_parameters      = location_remove_parameters,
      random_effects_summary = "none",
      formula_prefix         = FALSE,
      title                  = if (is_scale || is_random) "Location" else "Meta-Regression"
    ),
    conditional_args         = list(
      transform_factors      = TRUE,
      transform_scaled       = !standardized_coefficients,
      keep_formulas          = "mu",
      remove_parameters      = location_remove_parameters,
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
    likelihood_note <- if (.is_data_exact_selection(object[["data"]])) {
      paste0(
        "The exact finite-vector likelihood applies selection jointly after ",
        "analytically marginalizing Gaussian random effects."
      )
    } else {
      paste0(
        "The approximate likelihood applies selection row by row conditional ",
        "on sampled random effects."
      )
    }
    paste(
      "P-value intervals for publication bias weights omega correspond to one-sided p-values.",
      likelihood_note
    )
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
  random_footnotes <- .summary_random_footnotes(object, conditional = FALSE)
  random_conditional_footnotes <- .summary_random_footnotes(
    object,
    conditional = TRUE
  )
  estimates_random_pair <- .summary_estimates_pair(
    enabled                  = .summary_random_components_enabled(object),
    object                   = object,
    probs                    = probs,
    include_mcmc_diagnostics = include_mcmc_diagnostics,
    is_robma                 = is_robma,
    conditional              = conditional,
    main_args                = list(
      keep_parameters         = "random",
      transform_scaled        = !standardized_coefficients,
      random_effects_summary  = "standard",
      random_effects_metadata = TRUE,
      formula_prefix          = FALSE,
      title                   = "Random",
      footnotes               = random_footnotes
    ),
    conditional_args         = list(
      keep_parameters         = "random",
      transform_scaled        = !standardized_coefficients,
      random_effects_summary  = "standard",
      random_effects_metadata = TRUE,
      formula_prefix          = FALSE,
      title                   = "Conditional Random",
      footnotes               = random_conditional_footnotes
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
      inclusion_scale      = list(),
      inclusion_random     = list()
    )
  }

  out <- list(
    name                        = .summary.brma_model_names(object),
    known_v_backend             = .brma_mv_known_v_backend_metadata(object),
    inclusion_components        = inclusion[["inclusion_components"]],
    inclusion_mods              = inclusion[["inclusion_mods"]],
    inclusion_scale             = inclusion[["inclusion_scale"]],
    inclusion_random            = inclusion[["inclusion_random"]],
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
    "inclusion_random",
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


#' @title Convert brma Summaries to a Data Frame
#'
#' @description Converts the non-empty tables displayed by a fitted
#' \code{brma} object or its \code{summary.brma} result to one plain,
#' long-form data frame. The leading
#' \code{component} column identifies the summary section and
#' \code{parameter} retains the displayed row label. Printed quantile labels
#' are returned as syntactic \code{CI_} column names.
#'
#' @param x a fitted \code{brma} or \code{summary.brma} object.
#' @param row.names \code{NULL} or a character vector giving the row names.
#' @param optional logical; passed to the final data-frame coercion.
#' @param stringsAsFactors accepted for compatibility with \code{data.frame()}.
#' @param ... for a fitted object, additional arguments passed to
#' \code{summary.brma()}; otherwise unused.
#'
#' @return A plain \code{data.frame} containing all displayed summary rows.
#'
#' @export
as.data.frame.brma <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  output <- as.data.frame.summary.brma(
    x                = summary(x, ...),
    row.names        = row.names,
    optional         = optional,
    stringsAsFactors = stringsAsFactors
  )

  return(output)
}


#' @rdname as.data.frame.brma
#' @export
as.data.frame.summary.brma <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  x_data <- .summary_brma_prepare_print_sections(x)
  components <- c(
    inclusion_components         = "inclusion",
    inclusion_mods               = "inclusion location",
    inclusion_scale              = "inclusion scale",
    inclusion_random             = "inclusion random",
    estimates                    = "common",
    estimates_conditional        = "conditional common",
    estimates_mods               = "location",
    estimates_mods_conditional   = "conditional location",
    estimates_scale              = "scale",
    estimates_scale_conditional  = "conditional scale",
    estimates_random             = "random",
    estimates_random_conditional = "conditional random",
    estimates_bias               = "bias",
    estimates_bias_conditional   = "conditional bias"
  )

  tables <- lapply(names(components), function(section) {
    table <- x_data[[section]]
    if (length(table) == 0L) {
      return(NULL)
    }

    .output_table_as_long_data_frame(
      table            = table,
      component        = components[[section]],
      stringsAsFactors = stringsAsFactors
    )
  })

  output <- .output_bind_long_data_frames(
    tables    = tables,
    row.names = row.names,
    optional  = optional
  )

  return(output)
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

  out <- random[, estimate_cols, drop = FALSE]

  return(out)
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

.summary_random_footnotes <- function(object, conditional) {

  design      <- .fitted_formula_design(object, "mu", required = FALSE)
  allocations <- design[["random_allocations"]]
  gated_roots <- .random_gated_root_allocations(allocations)
  if (length(gated_roots) == 0L) {
    return(NULL)
  }

  component_clause <- if (conditional) {
    "Component SDs condition on their own inclusion gate."
  } else {
    "Component SDs include their excluded zero branches."
  }
  has_gated_aggregate <- any(vapply(
    gated_roots,
    function(allocation) length(allocation[["terms"]]) > 1L,
    logical(1)
  ))
  if (has_gated_aggregate) {
    return(paste0(
      "sd_total and var_prop(...) describe the slab allocation before ",
      "independent component gates. ", component_clause
    ))
  }

  component_clause
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
    simplify_names     = TRUE,
    remove_diagnostics = !include_mcmc_diagnostics,
    remove_inclusion   = if (conditional) TRUE else is_robma,
    remove_spike_0     = FALSE,
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
    selection_type <- if (.is_data_exact_selection(object[["data"]])) {
      "Exact Selection"
    } else {
      "Approximate Selection"
    }
    model_name <- paste(model_name, selection_type)
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
  parameter_roles <- attr(inclusion, "parameter_roles", exact = TRUE)
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

  random_indices      <- which(parameter_roles == "random_inclusion")
  random_slab_indices <- which(parameter_roles == "random_slab")

  mods_indices <- grep("^mu_", parameters)
  mods_indices <- mods_indices[parameters[mods_indices] != "mu_intercept"]
  mods_indices <- setdiff(
    mods_indices,
    c(random_indices, random_slab_indices)
  )

  scale_indices <- grep("^log_tau_", parameters)
  scale_indices <- scale_indices[parameters[scale_indices] != "log_tau_intercept"]
  scale_indices <- setdiff(
    scale_indices,
    c(random_indices, random_slab_indices)
  )

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
    ),
    inclusion_random = .summary.inclusion_subtable(
      table      = inclusion,
      indices    = random_indices,
      row_labels = sub(
        pattern     = "^.*inclusion\\((.*)\\)$",
        replacement = "\\1",
        x           = row_labels[random_indices]
      ),
      title      = "Random-Effect Inclusion"
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
