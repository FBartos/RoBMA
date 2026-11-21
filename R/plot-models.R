#' Plot posterior model probabilities
#'
#' Display the posterior probabilities of the individual models in a fitted
#' model-averaged meta-analysis. Each column is one full model and its width is
#' its posterior probability. Rows describe the model's fitted components.
#'
#' @param x a fitted model-averaged object, such as a fit from [RoBMA()] or
#'   [BMA()].
#' @param ... additional arguments for methods.
#' @param parameter an ordered character vector of model component, term, or
#'   factor-level selectors. Use `NULL` to show the model components, including
#'   random-effect rows, followed by the location and scale terms.
#'   Random-effect inclusion rows use the labels printed by [summary.brma()]
#'   and appear before publication bias.
#'   Use qualified selectors such as `"mu_year"`, `"log_tau_year"`, or
#'   `"mu_delivery[online]"` to disambiguate terms or levels.
#' @param component component to select from: `"auto"`, `"mods"` (also
#'   `"location"`), `"scale"`, `"bias"`, or `"random"`.
#' @param factor_view `"levels"` shows a row per fitted factor level and
#'   `"terms"` shows one neutral inclusion row per factor term. Supply a single
#'   setting or named overrides for individual terms. Unspecified terms use
#'   `"levels"`.
#' @param mixture_view `"inclusion"` displays each member's declared null or
#'   alternative status. `"model"` shows one categorical row for the mixture.
#'   `"individual"` shows one row for each declared member, including null
#'   members. Supply a single setting or named overrides for individual terms.
#'   Unspecified terms use `"inclusion"`.
#' @param group_labels whether to display group headings beside the rows.
#' @param labels named character vector of row label overrides. Use stable
#'   selectors such as `"mu"`, `"mu_year"`, `"mu_delivery[online]"`, or
#'   `"bias[PET]"` (when the mixture member is named `PET`).
#'   Overrides also update the mixture labels in the legend.
#' @param decreasing whether to draw models in decreasing posterior
#'   probability order.
#' @param col named colors for `positive`, `negative`, `included`, and
#'   `excluded`. Defaults are blue, red, black, and white, respectively.
#' @param col_bias named colors for `bias` and `no_bias` or for publication-bias
#'   mixture members, using their names from [summary_models()] or qualified member
#'   selectors of the form `"bias[PET]"`. Member-specific colors override
#'   bias and no-bias colors in categorical and individual views.
#'   Use qualified selectors for members named `bias` or `no_bias`.
#'   This argument requires a displayed publication-bias component.
#' @param legend whether to show the color legend.
#' @param legend_position location of the legend, `"bottom"` or `"top"`.
#' @param plot_type `"base"` for base graphics or `"ggplot"` for ggplot2.
#'
#' @details
#' The column widths come from the same joint indicator draws as
#' `summary_models(x, type = "individual")`. Choosing rows changes the
#' displayed information while keeping the fitted models and their widths.
#' Ungated coefficients can appear alongside the model components.
#' Select the main components by `"mu"`, `"tau"`, `"bias"`, and `"rho"`
#' for an available heterogeneity allocation mixture. Use qualified selectors
#' for moderators with those names.
#'
#' Inclusion uses the null or alternative assignment of the active prior
#' member, or the fitted inclusion gate. In the inclusion view, directional
#' colors indicate the sign of the posterior median within that full model on
#' the fitted predictor scale. If the predictors were centered or standardized
#' during fitting, intercept and moderator directions refer to that
#' parameterization. Heterogeneity, heterogeneity allocation, publication
#' bias, random-effect gates, whole factors, and coefficients with a median of
#' zero use neutral inclusion colors. Factor level rows share the factor's
#' inclusion gate.
#' Treatment contrasts show coefficients relative to the reference level and
#' omit the reference row. Mean-difference and orthonormal contrasts show
#' deviations from the equally weighted mean of the factor levels.
#'
#' @return Invisibly `NULL` for base graphics or a ggplot object for ggplot2.
#'   The ggplot object has an integer `n_rows` attribute giving the number of
#'   displayed rows after selection and expansion, excluding the group headings
#'   and legend.
#'
#' @examples
#' \dontrun{
#' fit <- RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1)
#' p <- plot_models(fit, parameter = c("mu", "bias"),
#'                  mixture_view = c(bias = "individual"), plot_type = "ggplot")
#' attr(p, "n_rows")
#' p
#' }
#' @seealso [summary_models()], [plot.brma()]
#' @md
#' @export
plot_models <- function(x, ...) {

  UseMethod("plot_models")
}


#' @export
plot_models.default <- function(x, ...) {

  stop("'plot_models' requires a fitted model-averaged RoBMA object.", call. = FALSE)
}


#' @rdname plot_models
#' @export
plot_models.RoBMA <- function(x, parameter = NULL, component = "auto",
                              factor_view = "levels", mixture_view = "inclusion",
                              group_labels = FALSE, labels = NULL, decreasing = TRUE,
                              col = NULL, col_bias = NULL, legend = FALSE,
                              legend_position = "bottom", plot_type = "base", ...) {

  .check_unused_dots(list(...), allowed = character(), caller = "plot_models()")
  BayesTools::check_bool(group_labels, "group_labels", allow_NA = FALSE)
  if (!is.null(col)) {
    .plot_models_check_named(col, "col")
    .plot_models_check_keys(col, c("positive", "negative", "included", "excluded"), "col")
  }
  data <- .plot_models_data(x, parameter, component, factor_view, mixture_view,
                           group_labels, labels, col_bias)
  BayesTools::plot_model_probabilities(
    model_probabilities = data[["probabilities"]],
    cells               = data[["cells"]],
    group_labels        = group_labels,
    decreasing          = decreasing,
    col                 = c(col, data[["colors"]]),
    legend              = legend,
    legend_position     = legend_position,
    state_labels        = data[["state_labels"]],
    plot_type           = plot_type
  )
}


.plot_models_data <- function(x, parameter = NULL, component = "auto",
                              factor_view = "levels", mixture_view = "inclusion",
                              group_labels = FALSE, labels = NULL, col_bias = NULL) {

  if (!.is_RoBMA(x) || is.null(x[["fit"]]) || length(x[["fit"]]) == 0L) {
    stop("'plot_models' requires a fitted model-averaged RoBMA object.", call. = FALSE)
  }
  component  <- .parameter_component_normalize(component)
  available_components <- .summary_models_components(x, include_fixed = TRUE)
  components <- Filter(function(info) is.null(info[["fixed"]]), available_components)
  if (length(components) == 0L) {
    stop("No mixture-prior model components found.", call. = FALSE)
  }
  terms <- .plot_models_terms(x, available_components)
  if (!identical(component, "auto")) {
    terms <- Filter(function(term) term[["component"]] == component, terms)
  }
  if (length(terms) == 0L) {
    stop("No model components are available for component = '", component, "'.",
         call. = FALSE)
  }
  factor_view <- .plot_models_view(
    factor_view, c("levels", "terms"), "factor_view",
    Filter(function(term) length(term[["levels"]]) > 0L, terms)
  )
  mixture_view <- .plot_models_view(
    mixture_view, c("inclusion", "model", "individual"), "mixture_view",
    Filter(function(term) is.null(term[["info"]][["fixed"]]), terms)
  )
  for (owner in setdiff(names(terms), names(mixture_view))) {
    mixture_view[[owner]] <- "inclusion"
  }

  if (is.null(parameter)) {
    selected <- names(terms)
  } else {
    BayesTools::check_char(parameter, "parameter", check_length = FALSE, allow_NA = FALSE)
    levels <- do.call(c, unname(lapply(terms, `[[`, "levels")))
    available <- c(lapply(terms, `[[`, "aliases"), lapply(levels, `[[`, "aliases"))
    selected <- .plot_models_match(parameter, available, "parameter")
  }
  if (length(selected) == 0L) {
    stop("'parameter' must select at least one model component or factor level.",
         call. = FALSE)
  }
  level_owners <- do.call(c, unname(lapply(terms, function(term) {
    ids <- names(term[["levels"]])
    stats::setNames(rep(term[["id"]], length(ids)), ids)
  })))
  owners <- ifelse(selected %in% names(terms), selected, level_owners[selected])
  if (!is.null(col_bias)) {
    .plot_models_check_named(col_bias, "col_bias")
    if (!"bias" %in% owners) {
      stop("'col_bias' requires a displayed publication-bias row.", call. = FALSE)
    }
    members <- terms[["bias"]][["info"]][["names"]]
    allowed <- c("bias", "no_bias", members,
      paste0("bias[", members, "]", recycle0 = TRUE),
      paste0("bias{", members, "}", recycle0 = TRUE))
    .plot_models_check_keys(col_bias, allowed, "col_bias")
  }
  groups <- vapply(terms[unique(owners)], `[[`, character(1), "group")
  draw_groups <- length(unique(groups)) > 1L && !group_labels

  posterior <- .get_posterior_samples(x[["fit"]])
  models    <- .plot_models_joint(components, posterior)
  rows <- lapply(seq_along(selected), function(i) {
    owner <- owners[[i]]
    .plot_models_rows(
      terms[[owner]], selected[[i]], factor_view[[owner]],
      mixture_view[[owner]], models, draw_groups, col_bias
    )
  })
  rows <- unlist(rows, recursive = FALSE, use.names = FALSE)
  row_ids <- vapply(rows, `[[`, character(1), "id")
  rows <- stats::setNames(rows[!duplicated(row_ids)], unique(row_ids))
  signed <- Filter(function(row) {
    !is.null(row[["estimate"]]) && any(row[["state"]] == "included")
  }, rows)
  if (length(signed) > 0L) {
    parameters <- unique(vapply(signed, `[[`, character(1), "parameter"))
    estimates <- .plot_models_estimate_draws(x, posterior, parameters)
    for (id in names(signed)) {
      rows[[id]][["state"]] <- .plot_models_sign(
        signed[[id]][["state"]], estimates, models[["draws"]], signed[[id]][["estimate"]]
      )
    }
  }
  if (!is.null(labels)) {
    .plot_models_check_named(labels, "labels")
    label_rows <- .plot_models_match(names(labels), lapply(rows, `[[`, "aliases"), "labels")
    for (i in seq_along(labels)) {
      id <- label_rows[[i]]
      row <- rows[[id]]
      row[["label"]] <- labels[[i]]
      if (!is.null(row[["members"]])) {
        row[["state_labels"]][] <- paste0(row[["label"]], ": ", row[["members"]])
      } else if (isTRUE(row[["individual"]])) {
        row[["state_labels"]][] <- row[["label"]]
      }
      rows[[id]] <- row
    }
  }
  row_cells <- lapply(rows, function(row) {
    data.frame(
      model = names(models[["probabilities"]]),
      row   = row[["id"]],
      label = row[["label"]],
      group = row[["group"]],
      state = row[["state"]],
      stringsAsFactors = FALSE
    )
  })
  colors <- unlist(lapply(unname(rows), `[[`, "colors"), use.names = TRUE)
  state_labels <- unlist(lapply(unname(rows), `[[`, "state_labels"), use.names = TRUE)
  list(
    probabilities = models[["probabilities"]],
    cells         = do.call(rbind, unname(row_cells)),
    colors        = colors[!duplicated(names(colors))],
    state_labels  = state_labels[!duplicated(names(state_labels))]
  )
}


.plot_models_terms <- function(x, components) {

  main <- c("Effect", "Heterogeneity", "Heterogeneity Allocation")
  location <- names(components)[startsWith(names(components), "Location: ")]
  scale <- names(components)[startsWith(names(components), "Scale: ")]
  random <- setdiff(names(components), c(main, "Publication Bias", location, scale))
  order <- c(intersect(c(main, random, "Publication Bias"), names(components)), location, scale)
  gates <- Filter(function(info) !is.null(info[["selector"]]), components[random])
  gate_labels <- if (length(gates) > 0L) {
    selectors <- vapply(gates, `[[`, character(1), "selector")
    labels <- .summary_random_inclusion_labels(x, selectors, names(gates))
    as.list(stats::setNames(.summary_random_component_labels(labels), names(gates)))
  }
  metadata <- if (length(c(location, scale))) .brma_parameter_catalog_metadata(x) else NULL
  catalog <- metadata[["catalog"]]
  quantities <- catalog[["quantities"]]
  catalog_aliases <- catalog[["aliases"]]
  scale_names <- .summary_scale_display_names(x)

  terms <- lapply(components[order], function(info) {
    label <- info[["component"]]
    parameter <- info[["parameter"]]
    component <- switch(label, Effect = "mods", Heterogeneity = "scale",
                        "Publication Bias" = "bias", "Heterogeneity Allocation" = "scale",
                        "random")
    term <- list(
      id        = switch(label, Effect = "mu", Heterogeneity = "tau",
                         "Publication Bias" = "bias", parameter),
      label     = if (is.null(gate_labels[[label]])) label else gate_labels[[label]],
      group     = "Components",
      component = component,
      aliases   = unique(c(parameter, label, gate_labels[[label]], info[["aliases"]], switch(label,
        Effect = "effect", Heterogeneity = "heterogeneity",
        "Heterogeneity Allocation" = "heterogeneity_allocation",
        "Publication Bias" = "publication_bias"))),
      signed    = identical(label, "Effect"),
      estimate  = if (parameter == "mu_intercept") "(mu) intercept" else parameter,
      levels    = list(),
      info      = info
    )
    if (!label %in% c(location, scale)) {
      return(term)
    }

    entry <- metadata[["entries"]][metadata[["entries"]][["parameter"]] == parameter, ,
                                    drop = FALSE]
    if (nrow(entry) != 1L) {
      .stop_refit_required("Fitted formula metadata for '", parameter,
                           "' are unavailable. Refit the model.")
    }
    is_location <- label %in% location
    formula <- entry[["formula_parameter"]]
    prefix  <- if (is_location) "mu" else unname(scale_names[[formula]])
    term[["label"]] <- if (entry[["term"]] == "intercept" && !is_location) {
      "exp(intercept)"
    } else {
      entry[["term"]]
    }
    term[["prefix"]] <- prefix
    term[["keep_prefix"]] <- !is_location && length(scale_names) > 1L
    term[["group"]] <- if (is_location) "Location" else "Scale"
    term[["component"]] <- entry[["component"]]
    term[["estimate"]] <- paste0("(", formula, ") ", term[["label"]])
    term[["aliases"]] <- unique(c(term[["aliases"]], entry[["aliases"]][[1L]],
      term[["estimate"]], paste0("(", prefix, ") ", term[["label"]])))

    level_rows <- quantities[match(entry[["member_quantity_ids"]][[1L]],
                                  quantities[["quantity_id"]]), , drop = FALSE]
    level_rows <- level_rows[level_rows[["source_type"]] != "structural_zero", , drop = FALSE]
    term[["signed"]] <- nrow(level_rows) == 0L && entry[["term"]] != "intercept"
    centered <- any(vapply(info[["prior"]], function(prior) {
      BayesTools::is.prior.meandif(prior) || BayesTools::is.prior.orthonormal(prior)
    }, logical(1)))
    term[["levels"]] <- lapply(seq_len(nrow(level_rows)), function(i) {
      quantity <- level_rows[i, , drop = FALSE]
      parts <- quantity[["label_parts"]][[1L]]
      if (centered) {
        parts[["transformation"]][[1L]] <- "dif"
      }
      label <- BayesTools::parameter_labels(list(parts), style = "table", formula_prefix = FALSE)
      list(
        id       = quantity[["canonical_name"]],
        label    = label,
        estimate = BayesTools::parameter_labels(list(parts), style = "table"),
        aliases  = setdiff(c(label, paste0("(", prefix, ") ", label),
          catalog_aliases[["alias"]][catalog_aliases[["quantity_id"]] == quantity[["quantity_id"]]]),
          term[["aliases"]])
      )
    })
    names(term[["levels"]]) <- level_rows[["canonical_name"]]
    term
  })
  stats::setNames(terms, vapply(terms, `[[`, character(1), "id"))
}


.plot_models_check_named <- function(values, argument) {

  BayesTools::check_char(values, argument, check_length = FALSE, allow_NA = FALSE)
  keys <- names(values)
  if (length(values) && (is.null(keys) || anyNA(keys) ||
                         any(!nzchar(keys)) || anyDuplicated(keys) ||
                         any(!nzchar(values)))) {
    stop("'", argument, "' must be a named character vector with unique ",
         "non-empty names and values.", call. = FALSE)
  }
}


.plot_models_check_keys <- function(values, allowed, argument) {

  unknown <- setdiff(names(values), allowed)
  if (length(unknown) > 0L) {
    stop("Unknown '", argument, "' ", if (length(unknown) == 1L) "key" else "keys",
         ": ", paste(sprintf("'%s'", unknown), collapse = ", "), ".", call. = FALSE)
  }
}


.plot_models_match <- function(selectors, aliases, argument) {

  if (any(!nzchar(selectors))) {
    stop("'", argument, "' selectors must be non-empty.", call. = FALSE)
  }
  vapply(selectors, function(selector) {
    if (selector %in% names(aliases)) {
      return(selector)
    }
    hits <- names(aliases)[vapply(aliases, function(keys) selector %in% keys, logical(1))]
    if (length(hits) != 1L) {
      problem <- if (length(hits) == 0L) {
        " is unavailable."
      } else {
        paste0(" is ambiguous. Use ", paste(sprintf("'%s'", hits), collapse = " or "), ".")
      }
      stop("'", argument, "' selector '", selector, "'", problem, call. = FALSE)
    }
    hits
  }, character(1), USE.NAMES = FALSE)
}


.plot_models_view <- function(view, choices, argument, terms) {

  BayesTools::check_char(view, argument, check_length = FALSE, allow_NA = FALSE)
  if (any(!view %in% choices)) {
    stop("'", argument, "' must use ", paste(sprintf("'%s'", choices), collapse = ", "),
         ".", call. = FALSE)
  }
  if (length(view) == 1L && is.null(names(view))) {
    return(as.list(stats::setNames(rep(view, length(terms)), names(terms))))
  }
  .plot_models_check_named(view, argument)
  overrides <- .plot_models_match(names(view), lapply(terms, `[[`, "aliases"), argument)
  result <- stats::setNames(rep(choices[[1L]], length(terms)), names(terms))
  result[overrides] <- view
  as.list(result)
}


.plot_models_joint <- function(components, posterior) {

  indicators <- lapply(components, function(info) {
    .summary_models_indicators(
      posterior_samples = posterior,
      parameter         = info[["parameter"]],
      prior             = info[["prior"]],
      column            = info[["indicator"]],
      offset            = info[["indicator_offset"]]
    )
  })
  keys <- do.call(paste, c(unname(indicators), sep = ":"))
  first <- which(!duplicated(keys))
  # Use the same tie order as the component grid in summary_models().
  first <- first[do.call(order, lapply(rev(unname(indicators)), `[`, first))]
  ids <- keys[first]
  index <- match(keys, ids)
  list(
    probabilities = stats::setNames(tabulate(index, length(ids)) / length(index), ids),
    indicators    = lapply(indicators, `[`, first),
    draws         = split(seq_along(index), index)
  )
}


.plot_models_estimate_draws <- function(x, posterior, parameters) {

  # Convert the exact source rows that supplied the joint indicators.
  fit <- BayesTools::JAGS_with_draws(x[["fit"]], posterior)
  draws <- BayesTools::JAGS_estimates_table(
    fit                    = fit,
    keep_parameters        = parameters,
    return_samples         = TRUE,
    transform_scaled       = FALSE,
    transform_factors      = TRUE,
    remove_spike_0         = FALSE,
    remove_inclusion       = TRUE,
    remove_diagnostics     = TRUE,
    random_effects_summary = "none",
    simplify_names         = TRUE
  )
  if (nrow(draws) != nrow(posterior)) {
    stop("Model indicators and coefficient draws have different row counts.", call. = FALSE)
  }
  draws
}


.plot_models_rows <- function(term, selected, factor_view, mixture_view,
                              models, draw_groups, col_bias) {

  info <- term[["info"]]
  indicator <- models[["indicators"]][[info[["component"]]]]
  included <- if (is.null(info[["fixed"]])) {
    info[["hypothesis"]][indicator] == "Alternative"
  } else {
    rep(info[["fixed"]], length(models[["probabilities"]]))
  }
  is_bias <- identical(term[["id"]], "bias")
  prefix <- if (draw_groups || isTRUE(term[["keep_prefix"]])) term[["prefix"]]
  row <- function(id, label, aliases, state, colors = NULL, state_labels = NULL,
                  estimate = NULL, members = NULL, individual = FALSE) {
    list(id = id, label = label, aliases = aliases, group = term[["group"]],
         state = state, colors = colors, state_labels = state_labels,
         estimate = estimate, parameter = info[["parameter"]],
         members = members, individual = individual)
  }
  term_label <- if (is.null(prefix)) term[["label"]] else paste0("(", prefix, ") ", term[["label"]])
  direct_level <- selected %in% names(term[["levels"]])

  if (mixture_view != "inclusion") {
    if (direct_level) {
      stop("Select the whole term '", term[["id"]],
           "' to display its mixture members.", call. = FALSE)
    }
    members <- info[["names"]]
    keys <- paste0(term[["id"]], " mixture member ", seq_along(members))
    member_ids <- paste0(term[["id"]], "{", members, "}")
    member_aliases <- paste0(term[["id"]], "[", members, "]")
    palette <- grDevices::hcl.colors(length(members), "Dark 3")
    if (is_bias) {
      types <- ifelse(info[["hypothesis"]] == "Null", "no_bias", "bias")
      replace <- types %in% names(col_bias)
      palette[replace] <- col_bias[types[replace]]
      qualified <- c(member_aliases, member_ids)
      bare_members <- ifelse(members %in% c("bias", "no_bias", qualified),
                              NA_character_, members)
      for (aliases in list(bare_members, member_aliases, member_ids)) {
        replace <- aliases %in% names(col_bias)
        palette[replace] <- col_bias[aliases[replace]]
      }
    }
    palette <- stats::setNames(palette, keys)
    legend_labels <- stats::setNames(paste0(term_label, ": ", members), keys)
    if (mixture_view == "model") {
      return(list(row(term[["id"]], term_label, term[["aliases"]],
                      keys[indicator], palette, legend_labels, members = members)))
    }
    return(lapply(seq_along(members), function(i) {
      row(member_ids[[i]], paste0(term_label, ": ", members[[i]]),
          c(member_aliases[[i]], paste0(term_label, ": ", members[[i]])),
          ifelse(indicator == i, keys[[i]], "excluded"),
          palette[i], legend_labels[i], individual = TRUE)
    }))
  }

  state <- ifelse(included, "included", "excluded")
  colors <- state_labels <- NULL
  if (is_bias) {
    custom <- intersect(names(col_bias), c("bias", "no_bias"))
    state[included & "bias" %in% custom] <- "bias"
    state[!included & "no_bias" %in% custom] <- "no_bias"
    colors <- col_bias[custom]
    state_labels <- c(bias = "Bias", no_bias = "No bias")[custom]
  }
  levels <- if (direct_level) {
    term[["levels"]][selected]
  } else if (identical(factor_view, "levels")) {
    term[["levels"]]
  }
  if (length(levels) > 0L) {
    return(lapply(levels, function(level) {
      label <- if (is.null(prefix)) level[["label"]] else paste0("(", prefix, ") ", level[["label"]])
      row(level[["id"]], label, level[["aliases"]],
          state, estimate = level[["estimate"]])
    }))
  }
  estimate <- if (isTRUE(term[["signed"]])) term[["estimate"]]
  list(row(term[["id"]], term_label, term[["aliases"]], state, colors, state_labels,
           estimate = estimate))
}


.plot_models_sign <- function(state, estimates, draw_rows, parameter) {

  included <- which(state == "included")
  if (length(included) == 0L) return(state)
  column <- which(colnames(estimates) == parameter)
  if (length(column) != 1L) {
    stop("Expected one posterior coefficient for '", parameter, "'.", call. = FALSE)
  }
  values <- estimates[, column]
  for (i in included) {
    draws <- values[draw_rows[[i]]]
    if (any(!is.finite(draws))) {
      stop("Non-finite posterior coefficient for '", parameter, "'.", call. = FALSE)
    }
    median <- stats::median(draws)
    if (median > 0) {
      state[[i]] <- "positive"
    } else if (median < 0) {
      state[[i]] <- "negative"
    }
  }
  state
}
