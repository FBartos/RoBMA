.check_and_select_plot_parameter <- function(parameter, parameter_mods,
                                             parameter_scale, object,
                                             component = "auto") {

  component           <- .parameter_component_normalize(component)
  has_parameter       <- !missing(parameter)       && !is.null(parameter)
  has_parameter_mods  <- !missing(parameter_mods)  && !is.null(parameter_mods)
  has_parameter_scale <- !missing(parameter_scale) && !is.null(parameter_scale)

  n_specified <- sum(c(has_parameter, has_parameter_mods, has_parameter_scale))
  if (n_specified > 1) {
    stop("Only one of 'parameter', 'parameter_mods', or 'parameter_scale' can be specified.", call. = FALSE)
  }

  if (has_parameter_mods) {
    .parameter_component_check_compatible(component, "mods", "parameter_mods")
    return(.brma_parameter_select(
      object    = object,
      parameter = parameter_mods,
      component = "mods",
      argument  = "parameter_mods"
    ))
  }

  if (has_parameter_scale) {
    .parameter_component_check_compatible(component, "scale", "parameter_scale")
    return(.brma_parameter_select(
      object    = object,
      parameter = parameter_scale,
      component = "scale",
      argument  = "parameter_scale"
    ))
  }

  if (!has_parameter) {
    parameter <- .brma_parameter_default(component, object)
  }

  return(.brma_parameter_select(
    object    = object,
    parameter = parameter,
    component = component,
    argument  = "parameter"
  ))
}

.component_values <- function(allow_auto = FALSE, allow_all = FALSE,
                              allow_outcome = FALSE, allow_bias = FALSE) {

  values <- c("mods", "location", "scale")
  if (allow_outcome) {
    values <- c("outcome", values)
  }
  if (allow_bias) {
    values <- c(values, "bias")
  }
  if (allow_all) {
    values <- c(values, "all")
  }
  if (allow_auto) {
    values <- c("auto", values)
  }

  return(values)
}

.component_normalize <- function(component, argument = "component",
                                 allow_auto = FALSE, allow_all = FALSE,
                                 allow_outcome = FALSE, allow_bias = FALSE,
                                 null = "auto",
                                 location_value = c("mods", "location")) {

  location_value <- match.arg(location_value)
  if (is.null(component)) {
    component <- null
  }

  BayesTools::check_char(component, argument, check_length = 1, allow_NA = FALSE)
  component <- match.arg(
    component,
    .component_values(
      allow_auto    = allow_auto,
      allow_all     = allow_all,
      allow_outcome = allow_outcome,
      allow_bias    = allow_bias
    )
  )

  if (component %in% c("mods", "location")) {
    return(location_value)
  }

  return(component)
}

.parameter_component_normalize <- function(component) {

  .component_normalize(
    component      = component,
    allow_auto     = TRUE,
    allow_bias     = TRUE,
    location_value = "mods"
  )
}

.fitted_component_normalize <- function(component) {

  .component_normalize(
    component      = component,
    allow_all      = TRUE,
    null           = "location",
    location_value = "location"
  )
}

.parameter_component_check_compatible <- function(component, selected_component,
                                                  argument) {

  if (!identical(component, "auto") &&
      !identical(component, selected_component)) {
    stop(
      "The '", argument, "' argument selects component = '", selected_component,
      "' but 'component' was set to '", component, "'.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}

.brma_parameter_default <- function(component, object) {

  if (identical(component, "mods")) {
    if (.is_mods(object)) {
      return("intercept")
    }
    return("mu")
  }
  if (identical(component, "scale")) {
    if (.is_scale(object)) {
      return("intercept")
    }
    return("tau")
  }
  if (identical(component, "bias")) {
    stop("Specify 'parameter' when component = 'bias'.", call. = FALSE)
  }

  if (.is_mods(object)) {
    return("mu")
  }

  return("mu")
}

.brma_parameter_select <- function(object, parameter,
                                   component = "auto",
                                   argument = "parameter") {

  component <- .parameter_component_normalize(component)
  BayesTools::check_char(parameter, argument, check_length = 1, allow_NA = FALSE)

  catalog <- .brma_parameter_catalog(object)
  matches <- catalog[catalog[["alias"]] == parameter, , drop = FALSE]
  if (!identical(component, "auto")) {
    matches <- matches[matches[["component"]] == component, , drop = FALSE]
  }

  if (nrow(matches) == 0L) {
    stop(
      "Unknown ", argument, " '", parameter, "'. Available quantities are: ",
      .brma_parameter_available(catalog, component), ".",
      call. = FALSE
    )
  }

  selected <- unique(matches[["parameter"]])
  if (length(selected) == 1L) {
    return(selected)
  }

  components <- unique(matches[["component"]])
  stop(
    "Parameter '", parameter, "' is ambiguous across components: ",
    paste0("'", components, "'", collapse = ", "),
    ". Set 'component' to 'mods'/'location', 'scale', or 'bias'.",
    call. = FALSE
  )
}

.brma_parameter_catalog <- function(object) {

  rows <- list()
  add <- function(parameter, component, term, aliases) {
    aliases <- unique(aliases[!is.na(aliases) & nzchar(aliases)])
    if (length(aliases) == 0L) {
      return(invisible(NULL))
    }
    for (alias in aliases) {
      rows[[length(rows) + 1L]] <<- data.frame(
        alias     = alias,
        parameter = parameter,
        component = component,
        term      = term,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    return(invisible(NULL))
  }

  if (.is_mods(object)) {
    add("mu_intercept", "mods", "intercept",
        c("mu_intercept", "mu", "intercept"))
    .brma_parameter_catalog_terms(
      object            = object,
      model_parameter   = "mu",
      component         = "mods",
      formula_parameter = "mu",
      add               = add
    )
  } else {
    add("mu", "mods", "mu", c("mu", "effect"))
  }

  if (.is_scale(object)) {
    add("log_tau_intercept", "scale", "intercept",
        c("log_tau_intercept", "tau", "intercept"))
    .brma_parameter_catalog_terms(
      object            = object,
      model_parameter   = "log_tau",
      component         = "scale",
      formula_parameter = "log_tau",
      add               = add
    )
  } else {
    add("tau", "scale", "tau", c("tau", "heterogeneity"))
  }

  if (.is_multilevel(object)) {
    add("rho", "scale", "rho", "rho")
  }
  if (.is_PET(object)) {
    add("PET", "bias", "PET", "PET")
  }
  if (.is_PEESE(object)) {
    add("PEESE", "bias", "PEESE", "PEESE")
  }
  if (.is_weightfunction(object)) {
    add("omega", "bias", "omega", c("omega", "weightfunction"))
  }

  fit_priors <- attr(object[["fit"]], "prior_list")
  if (!is.null(fit_priors[["bias"]])) {
    add("bias", "bias", "bias", "bias")
  }

  out <- do.call(rbind, rows)
  out <- unique(out)
  rownames(out) <- NULL
  return(out)
}

.brma_parameter_catalog_terms <- function(object, model_parameter, component,
                                          formula_parameter, add) {

  terms <- .fitted_formula_terms(
    object            = object,
    parameter         = model_parameter,
    include_intercept = FALSE,
    display           = TRUE,
    required          = FALSE
  )
  if (is.null(terms) || length(terms) == 0L) {
    return(invisible(NULL))
  }

  for (term in terms) {
    parameter <- BayesTools::JAGS_parameter_names(
      parameters        = term,
      formula_parameter = formula_parameter
    )
    add(parameter, component, term, c(parameter, term))
  }

  return(invisible(NULL))
}

.brma_parameter_available <- function(catalog, component = "auto") {

  component <- .parameter_component_normalize(component)
  if (!identical(component, "auto")) {
    catalog <- catalog[catalog[["component"]] == component, , drop = FALSE]
  }

  available <- unique(catalog[["alias"]])
  available <- sort(available[nzchar(available)])

  paste0("'", available, "'", collapse = ", ")
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
    object, parameter, parameter_mods, parameter_scale,
    component = "auto", allow_mixed_bias = FALSE) {

  component <- .parameter_component_normalize(component)
  n_specified <- sum(!vapply(
    list(parameter, parameter_mods, parameter_scale),
    is.null,
    logical(1)
  ))

  if (n_specified == 0 && identical(component, "auto")) {
    parameter <- "mu"
  }

  priors <- object[["priors"]]

  if (!is.null(parameter_mods)) {
    .parameter_component_check_compatible(component, "mods", "parameter_mods")
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
    .parameter_component_check_compatible(component, "scale", "parameter_scale")
    .check_plot_prior_name(parameter_scale, "parameter_scale")
    return(.select_plot_prior_term(
      prior_list = priors[["scale"]],
      term       = parameter_scale,
      argument   = "parameter_scale",
      prefix     = "log_tau",
      source     = "scale"
    ))
  }

  if (!identical(component, "auto")) {
    return(.select_plot_prior_component(
      priors           = priors,
      parameter        = parameter,
      component        = component,
      allow_mixed_bias = allow_mixed_bias
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
      "The term '%s' is available in both moderator and scale priors. Set 'component' to 'mods'/'location' or 'scale'.",
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
    "Unknown prior parameter '%s'. Available selections are: %s.",
    parameter,
    .collapse_plot_prior_names(priors)
  ), call. = FALSE)
}

.select_plot_prior_component <- function(priors, parameter, component,
                                         allow_mixed_bias) {

  if (is.null(parameter)) {
    parameter <- switch(
      component,
      "mods"    = "intercept",
      "scale"   = "intercept",
      "bias"    = "bias"
    )
  }

  .check_plot_prior_name(parameter, "parameter")

  if (identical(component, "mods")) {
    parameter <- .select_plot_prior_component_term(parameter, "mods")
    if (!is.null(priors[["mods"]]) && parameter %in% names(priors[["mods"]])) {
      return(.select_plot_prior_term(
        prior_list = priors[["mods"]],
        term       = parameter,
        argument   = "parameter",
        prefix     = "mu",
        source     = "mods"
      ))
    }
    if (identical(parameter, "intercept") &&
        !is.null(priors[["outcome"]][["mu"]])) {
      return(list(
        prior  = priors[["outcome"]][["mu"]],
        label  = "mu",
        source = "outcome",
        term   = "mu"
      ))
    }
    stop(sprintf("Unknown location prior parameter '%s'.", parameter), call. = FALSE)
  }

  if (identical(component, "scale")) {
    parameter <- .select_plot_prior_component_term(parameter, "scale")
    if (!is.null(priors[["scale"]]) && parameter %in% names(priors[["scale"]])) {
      return(.select_plot_prior_term(
        prior_list = priors[["scale"]],
        term       = parameter,
        argument   = "parameter",
        prefix     = "log_tau",
        source     = "scale"
      ))
    }
    if (identical(parameter, "intercept") &&
        !is.null(priors[["outcome"]][["tau"]])) {
      return(list(
        prior  = priors[["outcome"]][["tau"]],
        label  = "tau",
        source = "outcome",
        term   = "tau"
      ))
    }
    if (!is.null(priors[["outcome"]][[parameter]]) &&
        parameter %in% c("tau", "rho")) {
      return(list(
        prior  = priors[["outcome"]][[parameter]],
        label  = parameter,
        source = "outcome",
        term   = parameter
      ))
    }
    stop(sprintf("Unknown scale prior parameter '%s'.", parameter), call. = FALSE)
  }

  parameter <- switch(
    parameter,
    "effect"         = "mu",
    "heterogeneity"  = "tau",
    "weightfunction" = "omega",
    parameter
  )

  if (identical(component, "bias")) {
    if (parameter %in% c("omega", "PET", "PEESE", "bias") &&
        !is.null(priors[["outcome"]][["bias"]])) {
      prior <- .select_plot_prior_bias(
        prior            = priors[["outcome"]][["bias"]],
        parameter        = parameter,
        allow_mixed_bias = allow_mixed_bias
      )
      return(list(prior = prior, label = parameter, source = "bias", term = parameter))
    }
    stop(sprintf("Unknown bias prior parameter '%s'.", parameter), call. = FALSE)
  }

  stop(sprintf("Unknown prior component '%s'.", component), call. = FALSE)
}

.select_plot_prior_component_term <- function(parameter, component) {

  if (identical(component, "mods")) {
    if (parameter %in% c("mu", "mu_intercept")) {
      return("intercept")
    }
    return(.formula_design_display_names(sub("^mu_", "", parameter)))
  }

  if (identical(component, "scale")) {
    if (parameter %in% c("tau", "log_tau_intercept")) {
      return("intercept")
    }
    return(.formula_design_display_names(sub("^log_tau_", "", parameter)))
  }

  return(parameter)
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

  selections <- paste0("parameter = '", names(priors[["outcome"]]), "'")

  if (!is.null(priors[["outcome"]][["bias"]])) {
    selections <- unique(c(
      selections,
      paste0("parameter = '", c("omega", "PET", "PEESE"), "', component = 'bias'")
    ))
  }

  if (!is.null(priors[["mods"]])) {
    selections <- c(
      selections,
      paste0("parameter = '", names(priors[["mods"]]), "', component = 'mods'")
    )
  }

  if (!is.null(priors[["scale"]])) {
    selections <- c(
      selections,
      paste0("parameter = '", names(priors[["scale"]]), "', component = 'scale'")
    )
  }

  return(paste(unique(selections), collapse = ", "))
}
