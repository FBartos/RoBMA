.check_and_select_plot_parameter <- function(parameter, parameter_mods, parameter_scale, object) {

  # Check for model characteristics
  is_mods           <- .is_mods(object)
  is_scale          <- .is_scale(object)
  is_multilevel     <- .is_multilevel(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)

  # Determine which arguments were provided
  has_parameter       <- !missing(parameter)       && !is.null(parameter)
  has_parameter_mods  <- !missing(parameter_mods)  && !is.null(parameter_mods)
  has_parameter_scale <- !missing(parameter_scale) && !is.null(parameter_scale)

  # Count how many were specified

  n_specified <- sum(c(has_parameter, has_parameter_mods, has_parameter_scale))

  # Error if more than one specified
  if (n_specified > 1) {
    stop("Only one of 'parameter', 'parameter_mods', or 'parameter_scale' can be specified.", call. = FALSE)
  }

  # Default behavior if none specified
  if (n_specified == 0) {
    if (is_mods) {
      # Default to intercept for meta-regression models
      return("mu_intercept")
    } else {
      # Default to mu for simple models
      return("mu")
    }
  }

  # Process based on which argument was specified
  if (has_parameter) {

    # Validate base parameters
    BayesTools::check_char(parameter, "parameter")

    # Handle mu parameter
    if (parameter == "mu") {
      if (is_mods) {
        # mu becomes mu_intercept when mods are present
        return("mu_intercept")
      } else {
        return("mu")
      }
    }

    # Handle tau parameter
    if (parameter == "tau") {
      if (is_scale) {
        # tau becomes log_tau_intercept when scale is present
        return("log_tau_intercept")
      } else {
        return("tau")
      }
    }

    # Handle rho parameter (only for multilevel models)
    if (parameter == "rho") {
      if (!is_multilevel) {
        stop("The 'rho' parameter is only available for multilevel models.", call. = FALSE)
      }
      return("rho")
    }

    ### Handle publication bias parameters
    if (parameter == "PET") {
      if (!is_PET) {
        stop("The 'PET' parameter is only available for PET models.", call. = FALSE)
      }
      return("PET")
    }

    if (parameter == "PEESE") {
      if (!is_PEESE) {
        stop("The 'PEESE' parameter is only available for PEESE models.", call. = FALSE)
      }
      return("PEESE")
    }

    if (parameter == "omega" || parameter == "weightfunction") {
      if (!is_weightfunction) {
        stop("The 'omega'/'weightfunction' parameter is only available for selection models.", call. = FALSE)
      }
      return("omega")
    }

    # Unknown base parameter
    # Build list of valid parameters dynamically
    valid_params <- c("'mu'", "'tau'")
    if (is_multilevel)     valid_params <- c(valid_params, "'rho'")
    if (is_PET)            valid_params <- c(valid_params, "'PET'")
    if (is_PEESE)          valid_params <- c(valid_params, "'PEESE'")
    if (is_weightfunction) valid_params <- c(valid_params, "'omega'/'weightfunction'")

    stop(sprintf(
      "Unknown parameter '%s'. Valid base parameters are: %s.",
      parameter,
      paste(valid_params, collapse = ", ")
    ), call. = FALSE)

  } else if (has_parameter_mods) {

    # Validate that mods are present
    if (!is_mods) {
      stop("The 'parameter_mods' argument can only be used for models with moderators.", call. = FALSE)
    }

    BayesTools::check_char(parameter_mods, "parameter_mods")

    available_terms <- .fitted_formula_terms(
      object            = object,
      parameter         = "mu",
      include_intercept = TRUE,
      display           = TRUE
    )

    # Validate the specified term exists
    if (!parameter_mods %in% available_terms) {
      stop(sprintf(
        "The specified 'parameter_mods' term '%s' is not available. Available terms are: %s.",
        parameter_mods,
        paste0("'", available_terms, "'", collapse = ", ")
      ), call. = FALSE)
    }

    return(BayesTools::JAGS_parameter_names(
      parameters        = parameter_mods,
      formula_parameter = "mu"
    ))

  } else if (has_parameter_scale) {

    # Validate that scale is present
    if (!is_scale) {
      stop("The 'parameter_scale' argument can only be used for location-scale models.", call. = FALSE)
    }

    BayesTools::check_char(parameter_scale, "parameter_scale")

    available_terms <- .fitted_formula_terms(
      object            = object,
      parameter         = "log_tau",
      include_intercept = TRUE,
      display           = TRUE
    )

    # Validate the specified term exists
    if (!parameter_scale %in% available_terms) {
      stop(sprintf(
        "The specified 'parameter_scale' term '%s' is not available. Available terms are: %s.",
        parameter_scale,
        paste0("'", available_terms, "'", collapse = ", ")
      ), call. = FALSE)
    }

    return(BayesTools::JAGS_parameter_names(
      parameters        = parameter_scale,
      formula_parameter = "log_tau"
    ))
  }
}

.set_dots_plot        <- function(..., n_levels = 1) {

  dots <- list(...)
  if (is.null(dots[["col"]]) & n_levels == 1) {
    dots[["col"]]      <- "black"
  }else if (is.null(dots[["col"]]) & n_levels > 1) {
    dots[["col"]]      <- .plot_level_palette(n_levels)
  }
  if (is.null(dots[["col.fill"]])) {
    dots[["col.fill"]] <- "#4D4D4D4C" # scales::alpha("grey30", .30)
  }

  return(dots)
}

.plot_dots_allowed <- function() {

  return(c(
    "lwd", "lty", "col", "col.fill", "xlab", "ylab", "main",
    "xlim", "ylim", "par_name", "legend", "legend_title",
    "legend_labels", "legend_position", "cex", "cex.axis",
    "cex.lab", "cex.main", "col.axis", "col.lab", "col.main",
    "las", "bty", "xaxs", "yaxs", "axes", "xaxt", "yaxt",
    "pch", "bg", "border", "width", "scale_y2", "color",
    "colour", "fill", "alpha", "size", "linewidth", "linetype"
  ))
}

.iwmde_diagnostics_dots_allowed <- function() {

  return(c(
    "hist_col", "hist_border", "kde_col", "kde_lwd",
    "iwmde_col", "iwmde_lwd"
  ))
}

.plot_level_palette   <- function(n_levels) {

  if (n_levels <= 0L) {
    return(character())
  }
  if (n_levels == 1L) {
    return("black")
  }

  okabe_levels <- min(n_levels + 1L, 9L)
  colors       <- grDevices::palette.colors(
    n       = okabe_levels,
    palette = "Okabe-Ito"
  )[-1L]

  if (length(colors) < n_levels) {
    colors <- c(
      colors,
      grDevices::hcl.colors(
        n       = n_levels - length(colors),
        palette = "Dark 3"
      )
    )
  }

  return(colors[seq_len(n_levels)])
}
.set_dots_prior       <- function(dots_prior) {

  if (is.null(dots_prior)) {
    dots_prior <- list()
  }

  if (is.null(dots_prior[["col"]])) {
    dots_prior[["col"]]      <- "grey60"
  }
  if (is.null(dots_prior[["lty"]])) {
    dots_prior[["lty"]]      <- 1
  }
  if (is.null(dots_prior[["col.fill"]])) {
    dots_prior[["col.fill"]] <- "#B3B3B34C" # scales::alpha("grey70", .30)
  }

  return(dots_prior)
}
.set_dots_diagnostics <- function(..., type, chains) {

  dots <- list(...)
  if (is.null(dots[["col"]])) {
    dots[["col"]] <- if(type == "autocorrelation") "black" else rev(scales::viridis_pal()(chains))
  }

  return(dots)
}
.get_samples_n_levels <- function(samples, parameter) {

  if (inherits(samples[[parameter]], "mixed_posteriors.factor")) {
    if (attr(samples[[parameter]], "orthonormal") || attr(samples[[parameter]], "meandif")) {
      n_levels <- length(attr(samples[[parameter]], "level_names"))
    } else if (attr(samples[[parameter]], "treatment")) {
      n_levels <- length(attr(samples[[parameter]], "level_names")) - 1
    } else {
      n_levels <- 1
    }
  } else {
    n_levels <- 1
  }

  return(n_levels)
}

.as_mixed_posteriors_parameters <- function(object, parameters) {

  fit_priors <- attr(object[["fit"]], "prior_list")

  if (!is.null(fit_priors[["bias"]])) {
    parameters[parameters %in% c("omega", "PET", "PEESE")] <- "bias"
    parameters <- unique(parameters)
  }

  return(parameters)
}

.plot_parameter_label <- function(parameter, effect_transform = NULL) {

  if (.is_effect_location_parameter(parameter)) {
    if (.effect_output_active(effect_transform)) {
      return(paste0("Effect Size (", effect_transform[["label"]], ")"))
    }

    return("Effect Size")
  }

  if (grepl("^mu_", parameter)) {
    return(paste0("Effect Size: ", .summary_parameter_label(sub("^mu_", "", parameter))))
  }

  if (parameter %in% c("tau", "log_tau_intercept")) {
    return("Heterogeneity")
  }

  if (grepl("^log_tau_", parameter)) {
    return(paste0("Heterogeneity: ", .summary_parameter_label(sub("^log_tau_", "", parameter))))
  }

  label <- switch(
    parameter,
    "rho"   = "Heterogeneity Allocation",
    "PET"   = "PET",
    "PEESE" = "PEESE",
    "omega" = "Publication Bias",
    "bias"  = "Publication Bias",
    parameter
  )

  return(label)
}

.select_plot_prior_parameter <- function(
    object, parameter, parameter_mods, parameter_scale, allow_mixed_bias = FALSE) {

  n_specified <- sum(!vapply(
    list(parameter, parameter_mods, parameter_scale),
    is.null,
    logical(1)
  ))

  if (n_specified == 0) {
    parameter <- "mu"
  }

  priors <- object[["priors"]]

  if (!is.null(parameter_mods)) {
    .check_plot_prior_name(parameter_mods, "parameter_mods")
    return(.select_plot_prior_term(
      prior_list = priors[["mods"]],
      term       = parameter_mods,
      argument   = "parameter_mods",
      prefix     = "mu",
      source     = "mods"
    ))
  }

  if (!is.null(parameter_scale)) {
    .check_plot_prior_name(parameter_scale, "parameter_scale")
    return(.select_plot_prior_term(
      prior_list = priors[["scale"]],
      term       = parameter_scale,
      argument   = "parameter_scale",
      prefix     = "log_tau",
      source     = "scale"
    ))
  }

  .check_plot_prior_name(parameter, "parameter")

  parameter <- switch(
    parameter,
    "effect"        = "mu",
    "heterogeneity" = "tau",
    "weightfunction" = "omega",
    parameter
  )

  if (parameter == "mu" && is.null(priors[["outcome"]][["mu"]]) && !is.null(priors[["mods"]][["intercept"]])) {
    return(list(
      prior  = priors[["mods"]][["intercept"]],
      label  = "mu_intercept",
      source = "mods",
      term   = "intercept"
    ))
  }

  if (parameter == "tau" && is.null(priors[["outcome"]][["tau"]]) && !is.null(priors[["scale"]][["intercept"]])) {
    return(list(
      prior  = priors[["scale"]][["intercept"]],
      label  = "log_tau_intercept",
      source = "scale",
      term   = "intercept"
    ))
  }

  if (parameter %in% c("omega", "PET", "PEESE", "bias") && !is.null(priors[["outcome"]][["bias"]])) {
    prior <- .select_plot_prior_bias(
      prior            = priors[["outcome"]][["bias"]],
      parameter        = parameter,
      allow_mixed_bias = allow_mixed_bias
    )
    return(list(prior = prior, label = parameter, source = "bias", term = parameter))
  }

  if (parameter %in% names(priors[["outcome"]]) && !is.null(priors[["outcome"]][[parameter]])) {
    return(list(prior = priors[["outcome"]][[parameter]], label = parameter, source = "outcome", term = parameter))
  }

  in_mods  <- !is.null(priors[["mods"]])  && parameter %in% names(priors[["mods"]])
  in_scale <- !is.null(priors[["scale"]]) && parameter %in% names(priors[["scale"]])

  if (in_mods && in_scale) {
    stop(sprintf(
      "The term '%s' is available in both moderator and scale priors. Use 'parameter_mods' or 'parameter_scale'.",
      parameter
    ), call. = FALSE)
  }

  if (in_mods) {
    return(list(
      prior  = priors[["mods"]][[parameter]],
      label  = paste0("mu_", parameter),
      source = "mods",
      term   = parameter
    ))
  }

  if (in_scale) {
    return(list(
      prior  = priors[["scale"]][[parameter]],
      label  = paste0("log_tau_", parameter),
      source = "scale",
      term   = parameter
    ))
  }

  stop(sprintf(
    "Unknown prior parameter '%s'. Available parameters are: %s.",
    parameter,
    .collapse_plot_prior_names(priors)
  ), call. = FALSE)
}

.select_plot_prior_term <- function(prior_list, term, argument, prefix, source) {

  if (is.null(prior_list)) {
    stop(sprintf("The '%s' argument can only be used when the object contains the corresponding priors.", argument), call. = FALSE)
  }

  if (!term %in% names(prior_list)) {
    stop(sprintf(
      "The specified '%s' term '%s' is not available. Available terms are: %s.",
      argument,
      term,
      paste0("'", names(prior_list), "'", collapse = ", ")
    ), call. = FALSE)
  }

  return(list(
    prior  = prior_list[[term]],
    label  = paste0(prefix, "_", term),
    source = source,
    term   = term
  ))
}

.select_plot_prior_bias <- function(prior, parameter, allow_mixed_bias = FALSE) {

  if (.prior_has_phacking(prior)) {
    .selection_stop_phacking_deferred()
  }

  if (!BayesTools::is.prior.mixture(prior)) {
    if (parameter == "bias" ||
        (parameter == "omega" && .prior_has_selection(prior)) ||
        (parameter == "PET" && BayesTools::is.prior.PET(prior)) ||
        (parameter == "PEESE" && BayesTools::is.prior.PEESE(prior))) {
      return(prior)
    }

    stop(sprintf("The publication-bias prior does not contain a '%s' component.", parameter), call. = FALSE)
  }

  has_weightfunction <- any(vapply(prior, .prior_has_selection, logical(1)))
  has_petpeese       <- any(vapply(prior, function(x) BayesTools::is.prior.PET(x) || BayesTools::is.prior.PEESE(x), logical(1)))

  if (parameter == "bias") {
    if (has_weightfunction && has_petpeese) {
      if (!allow_mixed_bias) {
        stop(
          "The publication-bias prior mixes weight-function and PET-PEESE components. ",
          "Use parameter = 'omega', 'PET', or 'PEESE' to plot one component type.",
          call. = FALSE
        )
      }
    }
    return(prior)
  }

  keep <- switch(
    parameter,
    "omega" = vapply(prior, .prior_has_selection, logical(1)),
    "PET"   = vapply(prior, BayesTools::is.prior.PET, logical(1)),
    "PEESE" = vapply(prior, BayesTools::is.prior.PEESE, logical(1))
  )

  if (!any(keep)) {
    stop(sprintf("The publication-bias prior does not contain a '%s' component.", parameter), call. = FALSE)
  }

  selected <- unclass(prior)[keep]
  if (length(selected) == 1) {
    return(selected[[1]])
  }

  class(selected) <- c("prior", "prior.mixture")
  attr(selected, "components")    <- rep("alternative", length(selected))
  attr(selected, "prior_weights") <- rep(1, length(selected))

  return(selected)
}

.append_print_prior_selection <- function(selected, prior, label, source, term) {

  if (is.null(prior)) {
    return(selected)
  }

  selected[[label]] <- list(
    prior  = prior,
    label  = label,
    source = source,
    term   = term
  )

  return(selected)
}

.select_print_prior_all <- function(object) {

  priors  <- object[["priors"]]
  outcome <- priors[["outcome"]]
  mods    <- priors[["mods"]]
  scale   <- priors[["scale"]]

  selected <- list()

  selected <- .append_print_prior_selection(
    selected = selected,
    prior    = outcome[["mu"]],
    label    = "mu",
    source   = "outcome",
    term     = "mu"
  )

  if (!is.null(mods)) {
    selected <- .append_print_prior_selection(
      selected = selected,
      prior    = mods[["intercept"]],
      label    = "mu_intercept",
      source   = "mods",
      term     = "intercept"
    )

    for (term in setdiff(names(mods), "intercept")) {
      selected <- .append_print_prior_selection(
        selected = selected,
        prior    = mods[[term]],
        label    = paste0("mu_", term),
        source   = "mods",
        term     = term
      )
    }
  }

  selected <- .append_print_prior_selection(
    selected = selected,
    prior    = outcome[["tau"]],
    label    = "tau",
    source   = "outcome",
    term     = "tau"
  )

  if (!is.null(scale)) {
    selected <- .append_print_prior_selection(
      selected = selected,
      prior    = scale[["intercept"]],
      label    = "log_tau_intercept",
      source   = "scale",
      term     = "intercept"
    )

    for (term in setdiff(names(scale), "intercept")) {
      selected <- .append_print_prior_selection(
        selected = selected,
        prior    = scale[[term]],
        label    = paste0("log_tau_", term),
        source   = "scale",
        term     = term
      )
    }
  }

  for (parameter in setdiff(names(outcome), c("mu", "tau"))) {
    selected <- .append_print_prior_selection(
      selected = selected,
      prior    = outcome[[parameter]],
      label    = parameter,
      source   = if (parameter == "bias") "bias" else "outcome",
      term     = parameter
    )
  }

  return(selected)
}

.print_prior_object <- function(x, label = NULL, ...) {

  dots <- list(...)
  silent <- isTRUE(dots[["silent"]])
  dots[["silent"]] <- TRUE

  output <- do.call(print, c(list(x = x), dots))

  if (!silent) {
    if (!is.null(label)) {
      output <- paste(output, collapse = "\n")
      output <- unlist(strsplit(output, "\n", fixed = TRUE), use.names = FALSE)
      output <- paste0("  ", output)
      output <- paste(output, collapse = "\n")
      cat(label, ":\n", sep = "")
    }
    cat(output, "\n", sep = "")
  }

  return(invisible(x))
}

.plot_prior_unstandardized <- function(object, selected, plot_type, dots) {

  if (!(selected[["source"]] %in% c("mods", "scale"))) {
    return(NULL)
  }

  if (!isTRUE(.standardize_continuous_predictors(object))) {
    return(NULL)
  }

  formula_info <- .plot_prior_formula_info(
    object = object,
    source = selected[["source"]]
  )

  if (is.null(formula_info) || is.null(formula_info[["formula_scale"]])) {
    return(NULL)
  }

  parameter <- gsub(":", "__xXx__", selected[["label"]], fixed = TRUE)
  if (!parameter %in% formula_info[["column_names"]]) {
    return(NULL)
  }

  n_points <- if (!is.null(dots[["n_points"]])) dots[["n_points"]] else 1000
  BayesTools::check_int(n_points, "n_points", lower = 2)

  par_name <- dots[["par_name"]]
  plot_dots <- dots
  plot_dots[c(
    "par_name", "n_points", "n_samples", "force_samples", "x_range_quant",
    "individual", "show_figures", "rescale_x", "transformation",
    "transformation_arguments", "transformation_settings"
  )] <- NULL

  plot <- do.call(
    BayesTools::plot_transformed_prior,
    c(
      list(
        prior_list               = formula_info[["prior_list"]],
        column_names             = formula_info[["column_names"]],
        formula_scale            = formula_info[["formula_scale"]],
        parameter                = parameter,
        n_points                 = n_points,
        x_range                  = dots[["xlim"]],
        transformation           = dots[["transformation"]],
        transformation_arguments = dots[["transformation_arguments"]],
        transformation_settings  = if (!is.null(dots[["transformation_settings"]])) {
          dots[["transformation_settings"]]
        } else {
          FALSE
        },
        plot_type                = plot_type,
        par_name                 = par_name
      ),
      plot_dots
    )
  )

  return(plot)
}

.plot_prior_formula_info <- function(object, source) {

  parameter <- .fitted_formula_parameter(source)
  design    <- .fitted_formula_design(
    object    = object,
    parameter = parameter,
    required  = TRUE
  )

  if (is.null(design[["formula_scale"]])) {
    return(NULL)
  }

  formula_scale <- list(design[["formula_scale"]])
  names(formula_scale) <- parameter

  return(list(
    prior_list    = design[["prior_list"]],
    column_names  = names(design[["prior_list"]]),
    formula_scale = formula_scale
  ))
}

.use_plot_prior_list_dispatch <- function(prior) {

  if (.is_prior_phacking(prior) || .prior_has_phacking(prior)) {
    .selection_stop_phacking_deferred()
  }

  return(
    (
      BayesTools::is.prior.factor(prior) &&
        !BayesTools::is.prior.mixture(prior) &&
        (BayesTools::is.prior.meandif(prior) || BayesTools::is.prior.orthonormal(prior))
    ) ||
      .is_prior_bias_kernel(prior)
  )
}

.prior_mixture_plot_list <- function(prior) {

  prior_list  <- unclass(prior)
  prior_names <- names(prior_list)

  attributes(prior_list) <- NULL
  names(prior_list)     <- prior_names

  return(prior_list)
}

.check_plot_prior_name <- function(x, argument) {

  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop(sprintf("The '%s' argument must be a single non-empty character string.", argument), call. = FALSE)
  }

  return(invisible(TRUE))
}

.collapse_plot_prior_names <- function(priors) {

  names_available <- names(priors[["outcome"]])

  if (!is.null(priors[["outcome"]][["bias"]])) {
    names_available <- unique(c(names_available, "omega", "PET", "PEESE"))
  }

  if (!is.null(priors[["mods"]])) {
    names_available <- c(names_available, paste0("parameter_mods = '", names(priors[["mods"]]), "'"))
  }

  if (!is.null(priors[["scale"]])) {
    names_available <- c(names_available, paste0("parameter_scale = '", names(priors[["scale"]]), "'"))
  }

  return(paste0("'", names_available, "'", collapse = ", "))
}
