# Prior plots and prints of the random-effect quantities of fitted objects:
# plot_prior() and print_prior() of selectors such as 'tau', 'tau_total',
# 'tau2' or 'tau2_prop(study)'.

# The selected prior of a random-effect catalog entry of a fitted object, from
# the persisted provenance of the quantity (never from posterior draws). A
# quantity that is exactly one declared prior under the identity transform
# selects that prior. Any other quantity selects the prior law BayesTools
# derives from the declared priors, parameter_prior_density(), which
# plot(prior = TRUE) overlays too, with the derivation of the quantity from
# its declared source priors. 'standardized_coefficients' selects the fitted
# (TRUE) or the original predictor scale of that law.
.select_plot_prior_random <- function(object, entry,
                                      standardized_coefficients = TRUE) {

  quantities <- entry[["selection"]][["quantities"]]
  provenance <- .brma_random_parameter_design_spec(
    object,
    as.list(.brma_random_parameter_specs(quantities)[1L, , drop = FALSE])
  )
  spec <- provenance[["spec"]]

  selected <- list(
    prior             = .brma_random_parameter_declared_prior(object, spec),
    label             = spec[["label"]],
    source            = entry[["source"]],
    term              = entry[["term"]],
    formula_parameter = entry[["formula_parameter"]]
  )
  if (!is.null(selected[["prior"]])) {
    return(selected)
  }

  fit <- object[["fit"]]
  if (standardized_coefficients) {
    attr(fit, "formula_scale") <- list()
  }
  law <- BayesTools::parameter_prior_density(fit, entry[["selection"]])
  if (is.null(law)) {
    reason <- paste0(
      "The prior of the random-effect quantity '", spec[["label"]],
      "' is unavailable: the fitted model declares no prior density for it."
    )
    if (!standardized_coefficients && length(attr(fit, "formula_scale")) > 0L) {
      reason <- paste0(
        reason, " Its prior on the fitted scale may be available with ",
        "'standardized_coefficients = TRUE'."
      )
    }
    .hypothesis_stop(.hypothesis_refusal(reason, "target"))
  }

  selected[["prior"]]      <- law
  selected[["derivation"]] <- .brma_random_prior_derivation(
    object, spec, provenance
  )

  return(selected)
}


# The derivation of a random-effect quantity from its declared source priors:
# the header line (the quantity as a function of the declared priors, in the
# notation of print_prior(object, "random")) and the named source priors it
# refers to. A quantity whose derivation is not described here is derived from
# the random-effect specification as a whole.
.brma_random_prior_derivation <- function(object, spec, provenance) {

  fit        <- object[["fit"]]
  prior_list <- attr(fit, "prior_list", exact = TRUE)
  label      <- spec[["label"]]
  quantity   <- spec[["quantity"]]

  allocation <- provenance[["allocation_definition"]]
  factors    <- NULL
  if (is.null(allocation) && isTRUE(spec[["allocation_derived"]])) {
    term <- provenance[["term"]]
    if (!is.null(term) && .marginalized_random_effect_has_allocation(term)) {
      allocation <- term[["sd_binding"]][["allocations"]][[1L]]
      column     <- .brma_random_parameter_component_column(term, spec)
      factors    <- tryCatch(
        .marginalized_random_effect_allocation_factors(term, column = column),
        error = function(error) NULL
      )
    }
  }
  chain    <- .brma_random_allocation_chain(
    attr(fit, "formula_design", exact = TRUE), spec, allocation
  )
  symbols  <- .brma_random_allocation_symbols(chain)
  product  <- function(factors) .brma_random_allocation_terms(factors, symbols)
  square   <- function(expression) paste0("(", expression, ")^2")
  # the SD symbol of a variance: its label without the block prefix
  sd_label <- sub("^.*: ", "", sub("tau2", "tau", label, fixed = TRUE))

  # the declared priors that the expression refers to
  sd_source  <- FALSE
  weights    <- character()
  indicators <- character()
  sources    <- list()
  refer_factors <- function(factors) {
    weights    <<- c(weights, unlist(lapply(factors, `[[`, "weight_name")))
    indicators <<- c(indicators, unlist(lapply(factors, `[[`, "inclusion_name")))
  }
  refer_allocation <- function(allocation, gates = TRUE) {
    weights <<- c(weights, allocation[["weight_name"]])
    if (gates) {
      indicators <<- c(indicators, unlist(lapply(
        allocation[["inclusion"]], `[[`, "indicator_name"
      )))
    }
  }

  expression <- NULL
  one_to_one <- identical(spec[["source_type"]], "one_to_one_transform") &&
    identical(spec[["source_transform"]], "square") &&
    !isTRUE(spec[["allocation_derived"]])
  if (quantity %in% c("var", "var_total", "var_common") && one_to_one) {
    # the square of one declared SD prior
    expression <- if (grepl("[^A-Za-z0-9_.]", sd_label)) {
      paste0("(", sd_label, ")^2")
    } else {
      paste0(sd_label, "^2")
    }
    source_prior <- .brma_random_parameter_source_prior(object, spec)
    if (!is.null(source_prior)) {
      sources[[sd_label]] <- source_prior
    }
  } else if (quantity %in% c("sd", "var") && !is.null(allocation) &&
             !is.null(factors)) {
    # the SD (variance) of a component of a variance allocation
    expression <- paste(c("sigma_total", product(factors)), collapse = " * ")
    if (identical(quantity, "var")) {
      expression <- square(expression)
    }
    sd_source <- TRUE
    refer_factors(factors)
  } else if (quantity %in% c("sd_total", "var_total") && !is.null(allocation)) {
    # the total of an allocation: the SD source times the factors of its
    # parents and the share of its active components
    gated <- length(allocation[["inclusion"]]) > 0L
    expression <- paste(c(
      "sigma_total",
      product(allocation[["parent_factors"]]),
      if (gated) {
        paste0(
          "sqrt(sum of I_c * ", symbols[["weight"]](allocation[["weight_name"]]),
          "[c] over the components c)"
        )
      }
    ), collapse = " * ")
    if (identical(quantity, "var_total")) {
      expression <- square(expression)
    }
    sd_source <- TRUE
    refer_factors(allocation[["parent_factors"]])
    if (gated) {
      refer_allocation(allocation)
    }
  } else if (identical(quantity, "var_prop") && !is.null(allocation) &&
             is.character(allocation[["weight_name"]])) {
    # the share of a component of the variance of its allocation
    index  <- spec[["allocation_index"]]
    weight <- symbols[["weight"]](allocation[["weight_name"]])
    gates  <- allocation[["inclusion"]]
    gate   <- Filter(function(gate) {
      identical(as.integer(gate[["index"]]), index)
    }, gates)
    of <- paste0(" of allocation '", allocation[["label"]], "'")
    expression <- if (length(gates) == 0L) {
      paste0(weight, "[", index, "], the share of the total variance", of)
    } else {
      paste0(
        if (length(gate) == 1L) paste0("I_", gate[[1L]][["component"]], " * "),
        weight, "[", index, "] / sum of I_c * ", weight,
        "[c] over the components c, the share of the total variance",
        of, " of the active components"
      )
    }
    refer_allocation(allocation)
  } else if (quantity %in% c("sd_mult", "var_mult") && !is.null(allocation) &&
             is.character(allocation[["weight_name"]])) {
    # the multiplier of the SD (variance) of a component of an allocation
    weight <- paste0(
      symbols[["weight"]](allocation[["weight_name"]]), "[",
      spec[["allocation_index"]], "]"
    )
    if (identical(allocation[["scale"]], "mean_variance")) {
      weight <- paste0(as.integer(allocation[["n_targets"]]), " * ", weight)
    }
    expression <- if (identical(quantity, "sd_mult")) {
      paste0("sqrt(", weight, ")")
    } else {
      weight
    }
    refer_allocation(allocation, gates = FALSE)
  }

  add_source <- function(name, prior_name) {
    if (is.character(prior_name) && length(prior_name) == 1L &&
        !is.na(prior_name) && prior_name %in% names(prior_list) &&
        BayesTools::is.prior(prior_list[[prior_name]])) {
      sources[[name]] <<- prior_list[[prior_name]]
    }
  }
  if (sd_source && length(chain) > 0L) {
    add_source("sigma_total", chain[[1L]][["source"]][["name"]])
  }
  for (allocation_i in chain) {
    weight_name <- allocation_i[["weight_name"]]
    if (is.character(weight_name) && weight_name %in% weights) {
      add_source(symbols[["weight"]](weight_name), weight_name)
    }
    for (gate in allocation_i[["inclusion"]]) {
      if (gate[["indicator_name"]] %in% indicators) {
        add_source(paste0("p_", gate[["component"]]), gate[["prob_name"]])
      }
    }
  }

  header <- if (is.null(expression)) {
    paste0(label, ": prior law derived from the declared random-effect priors")
  } else {
    paste0(label, " = ", expression)
  }
  if (length(sources) == 0L &&
      inherits(object[["priors"]][["random"]], "prior_random")) {
    sources[["random"]] <- object[["priors"]][["random"]]
  }

  return(list(header = header, sources = sources))
}


# The variance allocation of a random-effect quantity and its ancestors, root
# first (an empty list without an allocation).
.brma_random_allocation_chain <- function(formula_design, spec, allocation) {

  chain <- list()
  while (!is.null(allocation) && length(chain) < 16L) {
    chain  <- c(list(allocation), chain)
    parent <- allocation[["parent"]][["allocation"]]
    allocation <- if (is.character(parent) && length(parent) == 1L &&
                      !is.na(parent) && nzchar(parent)) {
      .brma_random_parameter_design_allocation(
        formula_design,
        list(allocation = parent, formula_parameter = spec[["formula_parameter"]])
      )
    }
  }

  return(chain)
}


# The symbols of the weights ('w', or 'w_<allocation>' when several
# allocations have weights) and of the inclusion indicators ('I_<component>')
# of an allocation chain, as print_prior(object, "random") names them.
.brma_random_allocation_symbols <- function(chain) {

  weights <- unique(unlist(lapply(chain, `[[`, "weight_name"), use.names = FALSE))

  list(
    weight    = function(name) {
      if (length(weights) <= 1L) {
        return("w")
      }
      owner <- Filter(function(allocation) {
        identical(allocation[["weight_name"]], name)
      }, chain)
      paste0("w_", owner[[1L]][["label"]])
    },
    indicator = function(name) {
      for (allocation in chain) {
        for (gate in allocation[["inclusion"]]) {
          if (identical(gate[["indicator_name"]], name)) {
            return(paste0("I_", gate[["component"]]))
          }
        }
      }
      "I"
    }
  )
}


# The factors of an allocated SD as the terms of its product: the inclusion
# indicators and the square roots of the (scaled) Dirichlet weights.
.brma_random_allocation_terms <- function(factors, symbols) {

  unlist(lapply(factors, function(factor) {
    gate   <- factor[["inclusion_name"]]
    weight <- factor[["weight_name"]]
    c(
      if (is.character(gate)) symbols[["indicator"]](gate),
      if (is.character(weight)) {
        term <- paste0(
          symbols[["weight"]](weight), "[", as.integer(factor[["index"]]), "]"
        )
        if (identical(factor[["scale"]], "mean_variance")) {
          term <- paste0(as.integer(factor[["n_targets"]]), " * ", term)
        }
        paste0("sqrt(", term, ")")
      }
    )
  }), use.names = FALSE)
}


# Plot the prior law of a random-effect quantity (a prior_linear_density of
# BayesTools::parameter_prior_density()) without a posterior, as the
# declared-prior path does: 'xlim' and 'n_points' select the plotted range and
# grid, 'par_name' the axis label, and the other arguments are graphical
# arguments. Continuous parts are drawn as curves and point masses as arrows
# on the probability axis. The sampling controls of plot() of a prior are
# no-ops for a law and dropped; 'x_seq' and 'x_range_quant' cannot be honored
# and stop.
.plot_prior_density_law <- function(law, plot_type, dots) {

  for (argument in c("x_seq", "x_range_quant")) {
    if (!is.null(dots[[argument]])) {
      stop(
        "The '", argument, "' argument is not available for the prior of a ",
        "derived random-effect quantity; use 'xlim' and 'n_points'.",
        call. = FALSE
      )
    }
  }
  n_points <- if (!is.null(dots[["n_points"]])) dots[["n_points"]] else 1000
  BayesTools::check_int(n_points, "n_points", lower = 2)
  BayesTools::check_real(dots[["xlim"]], "xlim", check_length = 2, allow_NULL = TRUE)

  # as for declared priors, the default range is the pretty range of the
  # support the law is built on
  xlim <- dots[["xlim"]]
  if (is.null(xlim) && is.null(dots[["transformation"]])) {
    support <- range(c(law[["density"]][["x"]], law[["points"]][["x"]]))
    if (all(is.finite(support)) && support[[2L]] > support[[1L]]) {
      xlim <- range(pretty(support))
    }
  }

  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    law,
    n_points                 = n_points,
    x_range                  = xlim,
    transformation           = dots[["transformation"]],
    transformation_arguments = dots[["transformation_arguments"]],
    transformation_settings  = isTRUE(dots[["transformation_settings"]])
  )
  graphical <- dots
  graphical[c(
    "xlim", "n_points", "n_samples", "force_samples", "individual",
    "show_figures", "rescale_x", "transformation", "transformation_arguments",
    "transformation_settings"
  )] <- NULL

  plot <- do.call(
    BayesTools:::.plot_prior_list.both,
    c(list(plot_data = plot_data, plot_type = plot_type), graphical)
  )

  if (plot_type == "ggplot") {
    return(plot)
  }

  return(invisible())
}


# Print a selected prior: a prior object under its label, and a derived
# random-effect prior as its derivation line followed by the declared source
# priors under their symbols.
.print_prior_selected <- function(selected, ...) {

  derivation <- selected[["derivation"]]
  if (is.null(derivation)) {
    return(.print_prior_object(
      selected[["prior"]], label = selected[["label"]], ...
    ))
  }

  if (!isTRUE(list(...)[["silent"]])) {
    cat(derivation[["header"]], "\n", sep = "")
    for (name in names(derivation[["sources"]])) {
      .print_prior_object(
        derivation[["sources"]][[name]], label = name, indent = "  ", ...
      )
    }
  }

  return(invisible(selected[["prior"]]))
}
