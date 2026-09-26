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
      argument  = "parameter_mods",
      allow_factor_cells = TRUE
    ))
  }

  if (has_parameter_scale) {
    .parameter_component_check_compatible(component, "scale", "parameter_scale")
    return(.brma_parameter_select(
      object    = object,
      parameter = parameter_scale,
      component = "scale",
      argument  = "parameter_scale",
      allow_factor_cells = TRUE
    ))
  }

  if (!has_parameter) {
    parameter <- .brma_parameter_default(component, object)
  }

  return(.brma_parameter_select(
    object    = object,
    parameter = parameter,
    component = component,
    argument  = "parameter",
    allow_factor_cells = TRUE
  ))
}

.component_values <- function(allow_auto = FALSE, allow_all = FALSE,
                              allow_outcome = FALSE, allow_bias = FALSE,
                              allow_random = FALSE) {

  values <- c("mods", "location", "scale")
  if (allow_outcome) {
    values <- c("outcome", values)
  }
  if (allow_bias) {
    values <- c(values, "bias")
  }
  if (allow_random) {
    values <- c(values, "random")
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
                                 allow_random = FALSE,
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
      allow_bias    = allow_bias,
      allow_random  = allow_random
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
    allow_random   = TRUE,
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
                                                  argument, class = NULL) {

  if (!identical(component, "auto") &&
      !identical(component, selected_component)) {
    .stop_component_mismatch(
      paste0(
        "The '", argument, "' argument selects component = '",
        selected_component, "' but 'component' was set to '", component, "'."
      ),
      class = class
    )
  }

  return(invisible(TRUE))
}


# A selection whose component differs from the requested 'component': the
# error of class "RoBMA_component_mismatch" at every entry point (plot(),
# the prior functions, hypothesis()); 'class' prefixes the classes of the
# entry point (hypothesis() adds "RoBMA_hypothesis_statement").
.stop_component_mismatch <- function(message, class = NULL) {

  stop(structure(
    class = c(class, "RoBMA_component_mismatch", "error", "condition"),
    list(message = message, call = NULL)
  ))
}
.brma_parameter_default <- function(component, object) {

  if (identical(component, "mods")) {
    if (.is_mods(object) || .is_random(object)) {
      return(.brma_parameter_default_formula(object, "mods"))
    }
    return("mu")
  }
  if (identical(component, "scale")) {
    if (.is_scale(object)) {
      scale_specs <- .data_scale_component_specs(object[["data"]])
      if (length(scale_specs) > 1L) {
        stop(
          "Specify 'parameter' when component = 'scale' for component-specific scale models.",
          call. = FALSE
        )
      }
      return(.brma_parameter_default_formula(object, "scale"))
    }
    if (!is.null(object[["priors"]][["outcome"]][["tau"]])) {
      return("tau")
    }
    stop("The object does not contain scale priors.", call. = FALSE)
  }
  if (identical(component, "bias")) {
    stop("Specify 'parameter' when component = 'bias'.", call. = FALSE)
  }
  if (identical(component, "random")) {
    stop("Specify 'parameter' when component = 'random'.", call. = FALSE)
  }

  if (.is_mods(object) || .is_random(object)) {
    return(.brma_parameter_default_formula(object, "mods"))
  }

  return("mu")
}

.brma_parameter_default_formula <- function(object, component) {

  catalog <- .brma_parameter_catalog(object)
  rows    <- catalog[catalog[["component"]] == component, , drop = FALSE]
  if (any(rows[["alias"]] == "intercept")) {
    return("intercept")
  }

  parameters <- unique(rows[["parameter"]])
  if (length(parameters) == 1L) {
    return(parameters)
  }

  stop(
    "Specify 'parameter' when component = '", component,
    "' because the fitted formula has no intercept.",
    call. = FALSE
  )
}

.brma_parameter_select <- function(object, parameter,
                                   component = "auto",
                                   argument = "parameter",
                                   allow_factor_cells = FALSE) {

  entry <- .brma_parameter_select_entry(
    object    = object,
    parameter = parameter,
    component = component,
    argument  = argument,
    allow_factor_cells = allow_factor_cells
  )

  return(entry[["parameter"]])
}

.brma_parameter_select_entry <- function(object, parameter,
                                         component = "auto",
                                         argument = "parameter",
                                         allow_factor_cells = FALSE) {

  component <- .parameter_component_normalize(component)
  BayesTools::check_char(parameter, argument, check_length = 1, allow_NA = FALSE)

  if (is.null(object[["fit"]])) {
    catalog <- .brma_parameter_catalog_unfitted(object)
    matches <- catalog[catalog[["alias"]] == parameter, , drop = FALSE]
    if (!identical(component, "auto")) {
      matches <- matches[matches[["component"]] == component, , drop = FALSE]
    }
    if (nrow(matches) == 0L) {
      stop(
        "The specified ", argument, " '", parameter, "' is not available. ",
        "Available quantities are: ",
        .brma_parameter_available(catalog, component), ".",
        call. = FALSE
      )
    }
    selected <- unique(matches[["parameter"]])
    if (length(selected) == 1L) {
      return(as.list(matches[1L, , drop = FALSE]))
    }
    quantities <- unique(paste0(
      "'", matches[["parameter"]], "' (", matches[["component"]], ")"
    ))
    stop(
      "Parameter '", parameter, "' is ambiguous across quantities: ",
      paste(quantities, collapse = ", "),
      ". Set 'component' or use the full parameter name.",
      call. = FALSE
    )
  }

  metadata <- .brma_parameter_catalog_metadata(object)
  selection <- tryCatch(
    BayesTools::parameter_catalog_resolve(
      catalog   = metadata[["catalog"]],
      alias     = parameter,
      component = if (identical(component, "auto")) NULL else component,
      simplify_names = TRUE
    ),
    BayesTools_parameter_not_found = function(error) {
      if (identical(component, "auto")) stop(error)
      parents <- metadata[["entries"]]
      parents <- parents[parents[["component"]] == component, , drop = FALSE]
      allowed <- unique(c(parents[["quantity_id"]],
        unlist(parents[["member_quantity_ids"]], use.names = FALSE)))
      catalog <- metadata[["catalog"]]
      quantities <- catalog[["quantities"]]
      aliases <- catalog[["aliases"]]
      candidates <- unique(c(
        quantities[["quantity_id"]][quantities[["canonical_name"]] == parameter],
        aliases[["quantity_id"]][aliases[["alias"]] == parameter]
      ))
      candidates <- intersect(candidates, allowed)
      if (length(candidates) != 1L) stop(error)
      quantity <- quantities[quantities[["quantity_id"]] == candidates, , drop = FALSE]
      BayesTools::parameter_catalog_resolve(catalog,
        alias = quantity[["canonical_name"]], namespace = quantity[["namespace"]],
        component = quantity[["component"]], simplify_names = TRUE)
    },
    BayesTools_parameter_ambiguous = function(error) {
      ambiguity <- .brma_parameter_catalog_group_ambiguity(
        entries      = metadata[["entries"]],
        quantity_ids = error[["candidates"]][["quantity_id"]],
        component    = component
      )
      if (!is.null(ambiguity[["component"]])) {
        return(BayesTools::parameter_catalog_resolve(
          catalog   = metadata[["catalog"]],
          alias     = parameter,
          component = ambiguity[["component"]],
          simplify_names = TRUE
        ))
      }
      stop(error)
    }
  )
  entry <- metadata[["entries"]][
    metadata[["entries"]][["quantity_id"]] == selection[["quantity_id"]],
    ,
    drop = FALSE
  ]
  if (nrow(entry) == 0L) {
    parent <- .brma_parameter_catalog_entries_for_quantities(
      metadata[["entries"]], selection[["quantity_id"]]
    )
    if (nrow(parent) == 1L &&
        identical(parent[["role"]], "formula_coefficient_group")) {
      if (!allow_factor_cells) {
        stop("Individual factor-cell selection is unavailable for this method. ",
          "Use 'parameter = \"", parent[["term"]], "\"' to select the whole factor term.",
          call. = FALSE)
      }
      out <- as.list(parent[1L, setdiff(names(parent), "aliases"), drop = FALSE])
      quantity <- selection[["quantities"]]
      out[["parent_parameter"]] <- out[["parameter"]]
      out[["parameter"]] <- quantity[["canonical_name"]]
      for (field in c("quantity_id", "role", "status", "fixed_value")) {
        out[[field]] <- quantity[[field]]
      }
      out[["selection"]] <- selection
      return(out)
    }
  }
  if (nrow(entry) != 1L) {
    coefficients <- .brma_contrast_coefficient_quantities(
      metadata, selection[["quantity_id"]]
    )
    if (nrow(coefficients) > 0L) {
      .brma_stop_contrast_coefficient(
        metadata    = metadata,
        selector    = parameter,
        coefficient = coefficients[1L, , drop = FALSE],
        hypothesis  = FALSE
      )
    }
    .stop_refit_required(
      "Resolved parameter metadata are unavailable. Refit the model with the ",
      "current RoBMA/BayesTools build."
    )
  }
  out <- as.list(entry[1L, setdiff(names(entry), "aliases"), drop = FALSE])
  out[["selection"]] <- selection
  return(out)
}

# The fitted coordinate that a catalog quantity is identical to: the coordinate
# of a coordinate key, or the only dependency of a factor-level key with unit
# weight (a direct level cell or a contrast coefficient '<parameter>{j}').
# Combinations of coordinates and structural levels have none.
.brma_catalog_key_coordinate <- function(key) {

  if (!is.list(key) || length(key[["dependencies"]]) != 1L) {
    return(NULL)
  }
  if (identical(key[["type"]], "coordinate")) {
    return(key[["dependencies"]][[1L]])
  }
  if (identical(key[["type"]], "factor_level") &&
      identical(unname(as.numeric(key[["weights"]])), 1)) {
    return(key[["dependencies"]][[1L]])
  }

  return(NULL)
}

# The fitted coordinate that a factor level (or cell) structurally is: a
# direct level cell of treatment, independent, or ordered coding. BayesTools
# labels every coordinate that is not structurally a level cell as a contrast
# coefficient '<parameter>{j}' holding that coordinate, so a level whose
# unit-weight coordinate a coefficient also holds is not a level cell. This
# keeps a mean-difference or orthonormal level whose design row is a unit
# vector in floating point (the first of four mean-difference levels) from
# being taken for that coordinate.
.brma_catalog_level_coordinate <- function(quantity, quantities) {

  key        <- quantity[["extraction_key"]][[1L]]
  coordinate <- .brma_catalog_key_coordinate(key)
  if (is.null(coordinate) || !identical(key[["type"]], "factor_level")) {
    return(coordinate)
  }
  holders <- vapply(quantities[["extraction_key"]], function(other) {
    is.list(other) && identical(other[["type"]], "factor_level") &&
      identical(.brma_catalog_key_coordinate(other), coordinate)
  }, logical(1))
  if (sum(holders) != 1L) {
    return(NULL)
  }

  return(coordinate)
}

# Catalog quantities that are factor contrast coefficients ('g{1}'): their
# label parts name a contrast coefficient rather than a level cell.
.brma_catalog_contrast_coefficient <- function(quantities) {

  vapply(quantities[["label_parts"]], function(parts) {
    !is.null(parts) && !is.na(parts[["coefficient"]])
  }, logical(1))
}

# Contrast coefficients '<term>{j}' of mean-difference, orthonormal, and
# ordered factors are catalog quantities that no RoBMA entry covers: RoBMA
# addresses factor terms through their level labels.
.brma_contrast_coefficient_quantities <- function(metadata, quantity_ids) {

  quantities <- metadata[["catalog"]][["quantities"]]
  if (!is.data.frame(quantities) || length(quantity_ids) == 0L) {
    return(data.frame())
  }
  rows <- quantities[
    quantities[["quantity_id"]] %in% quantity_ids,
    ,
    drop = FALSE
  ]
  entries <- metadata[["entries"]]
  covered <- rows[["quantity_id"]] %in% c(
    entries[["quantity_id"]],
    unlist(entries[["member_quantity_ids"]], use.names = FALSE)
  )
  factor_level <- vapply(rows[["extraction_key"]], function(key) {
    is.list(key) && identical(key[["type"]], "factor_level")
  }, logical(1))
  coefficient <- .brma_catalog_contrast_coefficient(rows)

  return(rows[(factor_level | coefficient) & !covered, , drop = FALSE])
}

# The public label of a factor term and the labels of its non-structural
# levels (the treatment reference level is structural).
.brma_factor_term_selectors <- function(metadata, entry) {

  aliases <- as.list(rep(entry[["parameter"]], length(entry[["aliases"]][[1L]])))
  names(aliases) <- entry[["aliases"]][[1L]]
  quantities <- metadata[["catalog"]][["quantities"]]
  levels     <- quantities[
    quantities[["quantity_id"]] %in%
      unlist(entry[["member_quantity_ids"]], use.names = FALSE) &
      quantities[["status"]] != "structural",
    ,
    drop = FALSE
  ]

  return(list(
    label  = .hypothesis_brma_alias_label(aliases, entry[["parameter"]]),
    levels = levels[["component"]]
  ))
}

# Stop for a selector of a factor contrast coefficient ('g{1}' or, for terms
# without coefficient quantities, any '{j}' selector), naming the level-label
# form that RoBMA supports. 'term_alias' identifies the term by one of its
# aliases when the selector did not resolve to a catalog quantity.
.brma_stop_contrast_coefficient <- function(metadata, selector,
                                            coefficient = NULL,
                                            term_alias = NULL,
                                            hypothesis = TRUE) {

  entries <- metadata[["entries"]]
  groups  <- entries[
    entries[["role"]] == "formula_coefficient_group",
    ,
    drop = FALSE
  ]
  group <- if (!is.null(coefficient)) {
    groups[
      groups[["term"]] == coefficient[["term"]] &
        groups[["formula_parameter"]] == coefficient[["formula_parameter"]],
      ,
      drop = FALSE
    ]
  } else {
    groups[groups[["parameter"]] %in% term_alias |
      vapply(groups[["aliases"]], function(aliases) {
        term_alias %in% aliases
      }, logical(1)), , drop = FALSE]
  }
  example <- NULL
  label   <- NULL
  if (nrow(group) == 1L) {
    selectors <- .brma_factor_term_selectors(metadata, group)
    label     <- selectors[["label"]]
    if (!is.null(coefficient)) {
      selector <- paste0(label, coefficient[["component"]])
    }
    if (length(selectors[["levels"]]) > 0L) {
      example <- paste0(label, "[", selectors[["levels"]][[1L]], "]")
    }
  }
  example <- if (is.null(example)) {
    "their labels"
  } else {
    paste0("their labels, such as '", example, "'")
  }

  # hypothesis(): a statement to restate on factor levels
  # ("RoBMA_hypothesis_statement").
  if (hypothesis) {
    .hypothesis_stop(list(
      reason = paste0(
        "Hypotheses on factor contrast coefficients such as '", selector,
        "' are not supported. State them on factor levels by ", example, "."
      ),
      class  = "RoBMA_hypothesis_statement"
    ))
  }
  stop(
    "Factor contrast coefficients such as '", selector, "' cannot be ",
    "selected. Select factor levels by ", example,
    if (!is.null(label)) paste0(", or the whole term '", label, "'"), ".",
    call. = FALSE
  )
}

.brma_parameter_catalog <- function(object) {

  if (is.null(object[["fit"]])) {
    return(.brma_parameter_catalog_unfitted(object))
  }

  entries <- .brma_parameter_catalog_metadata(object)[["entries"]]
  rows <- lapply(seq_len(nrow(entries)), function(i) {
    aliases <- entries[["aliases"]][[i]]
    data.frame(
      alias             = aliases,
      parameter         = entries[["parameter"]][[i]],
      component         = entries[["component"]][[i]],
      term              = entries[["term"]][[i]],
      source            = entries[["source"]][[i]],
      formula_parameter = entries[["formula_parameter"]][[i]],
      quantity_id       = entries[["quantity_id"]][[i]],
      stringsAsFactors  = FALSE,
      check.names       = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- unique(out)
  rownames(out) <- NULL
  return(out)
}


.brma_parameter_catalog_metadata <- function(object) {

  .brma_validate_fit_contract(
    object,
    requires = c(
      "name_encoding",
      "formula_name_map",
      "formula_design",
      "parameter_map"
    )
  )

  # The metadata is derived from the fitted map, but also from the object's
  # data and priors, so those have to be part of the cache key: the map alone
  # would let a differently specified object read back another one's entries.
  BayesTools::parameter_map_cache(
    map      = BayesTools::parameter_map(object[["fit"]]),
    provider = "RoBMA",
    key      = list(
      data   = object[["data"]],
      priors = object[["priors"]]
    ),
    compute  = function() .brma_parameter_catalog_metadata_compute(object)
  )
}

.brma_parameter_catalog_metadata_compute <- function(object) {

  catalog    <- BayesTools::parameter_catalog(object[["fit"]])
  quantities <- catalog[["quantities"]]
  entries    <- list()
  aliases    <- list()
  extension_quantities <- list()

  add_entry <- function(quantity, parameter, component, term, source,
                        formula_parameter, entry_aliases,
                        entry_alias_simplified = rep(FALSE,
                          length(entry_aliases)),
                        member_quantity_ids = character()) {
    keep_aliases <- !is.na(entry_aliases) & nzchar(entry_aliases)
    alias_metadata <- unique(data.frame(
      alias      = entry_aliases[keep_aliases],
      simplified = entry_alias_simplified[keep_aliases],
      stringsAsFactors = FALSE
    ))
    entry_aliases <- unique(alias_metadata[["alias"]])
    entry <- data.frame(
      quantity_id       = quantity[["quantity_id"]],
      parameter         = parameter,
      component         = component,
      term              = term,
      source            = source,
      formula_parameter = formula_parameter,
      role              = quantity[["role"]],
      quantity          = quantity[["quantity"]],
      status            = quantity[["status"]],
      fixed_value       = quantity[["fixed_value"]],
      stringsAsFactors  = FALSE,
      check.names       = FALSE
    )
    entry[["aliases"]] <- I(list(entry_aliases))
    entry[["member_quantity_ids"]] <- I(list(unique(member_quantity_ids)))
    entries[[length(entries) + 1L]] <<- entry
    alias_rows <- data.frame(
      alias       = alias_metadata[["alias"]],
      quantity_id = rep(quantity[["quantity_id"]], nrow(alias_metadata)),
      namespace   = rep(component, nrow(alias_metadata)),
      component   = rep(component, nrow(alias_metadata)),
      simplified  = alias_metadata[["simplified"]],
      stringsAsFactors = FALSE
    )
    aliases[[length(aliases) + 1L]] <<- alias_rows[
      , c("alias", "quantity_id", "namespace", "component", "simplified"),
      drop = FALSE
    ]
    invisible(NULL)
  }

  public <- !quantities[["internal"]]
  add_ordinary <- function(parameter, component, term, source, entry_aliases) {
    rows <- which(public & quantities[["canonical_name"]] == parameter)
    if (length(rows) == 1L) {
      add_entry(
        quantity          = as.list(quantities[rows, , drop = FALSE]),
        parameter         = parameter,
        component         = component,
        term              = term,
        source            = source,
        formula_parameter = NA_character_,
        entry_aliases     = entry_aliases
      )
    }
  }
  if (!.is_mods(object) && !.is_random(object)) {
    add_ordinary("mu", "mods", "mu", "outcome", c("mu", "effect"))
  }
  if (!.is_scale(object)) {
    add_ordinary(
      "tau", "scale", "tau", "outcome", c("tau", "heterogeneity")
    )
  }
  if (.is_multilevel(object)) {
    add_ordinary("rho", "scale", "rho", "outcome", "rho")
  }

  formula_specs <- list()
  if (.is_mods(object) || .is_random(object)) {
    formula_specs[["mu"]] <- list(
      component = "mods",
      source    = if (.is_random(object)) "location" else "mods"
    )
  }
  if (.is_scale(object)) {
    scale_specs <- .data_scale_component_specs(object[["data"]])
    for (scale_spec in scale_specs) {
      formula_specs[[scale_spec[["parameter"]]]] <- list(
        component  = "scale",
        source     = "scale",
        scale_spec = scale_spec,
        single     = length(scale_specs) == 1L
      )
    }
  }
  for (formula_parameter in names(formula_specs)) {
    spec     <- formula_specs[[formula_parameter]]
    name_map <- .fitted_formula_name_map(
      object    = object,
      parameter = formula_parameter,
      required  = TRUE
    )
    fixed_rows <- name_map[
      name_map[["kind"]] == "fixed" &
        name_map[["role"]] == "coefficient",
      ,
      drop = FALSE
    ]
    formula_design <- BayesTools::JAGS_formula_design(
      object[["fit"]], formula_parameter
    )
    omit_location_intercept <-
      identical(formula_parameter, "mu") &&
      .location_omit_fixed_zero_intercept(object)
    if (!is.null(formula_design) &&
        (!.fitted_formula_has_intercept(
          object, formula_parameter, required = TRUE
        ) || omit_location_intercept)) {
      fixed_rows <- fixed_rows[
        fixed_rows[["term"]] != "intercept",
        ,
        drop = FALSE
      ]
    }
    for (row in seq_len(nrow(fixed_rows))) {
      map_row <- fixed_rows[row, , drop = FALSE]
      term    <- map_row[["term"]]
      entry_aliases <- c(map_row[["jags_name"]], term)
      if (identical(formula_parameter, "mu") && identical(term, "intercept")) {
        entry_aliases <- c(entry_aliases, "mu", "effect", "intercept")
      }
      if (identical(spec[["component"]], "scale") &&
          identical(term, "intercept")) {
        scale_spec <- spec[["scale_spec"]]
        if (isTRUE(spec[["single"]]) &&
            identical(formula_parameter, "log_tau")) {
          entry_aliases <- c(
            entry_aliases, "log_tau_intercept", "tau", "heterogeneity",
            "intercept"
          )
        } else {
          entry_aliases <- c(
            entry_aliases,
            scale_spec[["display_name"]],
            scale_spec[["aliases"]],
            paste0(scale_spec[["aliases"]], "_intercept")
          )
        }
      }
      entry_aliases <- unique(entry_aliases)
      # The term's coefficients, or the level cells of a factor term (its
      # contrast coefficients are no level cells).
      coordinate_rows <- which(
        public &
          quantities[["formula_parameter"]] == formula_parameter &
          quantities[["role"]] == "fixed_coefficient" &
          quantities[["term"]] == term &
          !.brma_catalog_contrast_coefficient(quantities)
      )
      if (length(coordinate_rows) == 0L) {
        .stop_refit_required(
          "Fitted coefficient metadata for formula parameter '",
          formula_parameter, "' and term '", term,
          "' are incomplete. Refit the model with the current ",
          "RoBMA/BayesTools build."
        )
      }
      if (length(coordinate_rows) == 1L) {
        quantity <- quantities[coordinate_rows, , drop = FALSE]
        member_quantity_ids <- character()
      } else {
        quantity <- .brma_parameter_catalog_formula_quantity(
          catalog            = catalog,
          map_row            = map_row,
          coordinates        = quantities[coordinate_rows, , drop = FALSE],
          semantic_component = spec[["component"]]
        )
        extension_quantities[[length(extension_quantities) + 1L]] <- quantity
        member_quantity_ids <- quantities[["quantity_id"]][coordinate_rows]
      }
      quantity <- as.list(quantity)
      add_entry(
        quantity          = quantity,
        parameter         = map_row[["jags_name"]],
        component         = spec[["component"]],
        term              = term,
        source            = spec[["source"]],
        formula_parameter = formula_parameter,
        entry_aliases     = entry_aliases,
        member_quantity_ids = member_quantity_ids
      )
    }
  }

  if (.is_random(object)) {
    rows <- which(
      public & quantities[["status"]] != "unavailable" &
        startsWith(quantities[["role"]], "random_") &
        quantities[["quantity"]] %in%
          .brma_random_parameter_supported_quantities()
    )
    for (row in rows) {
      quantity <- as.list(quantities[row, , drop = FALSE])
      # RoBMA selects random-effect quantities by its own quantity names,
      # rendered from the catalog label parts.
      parameter  <- .brma_random_parameter_io_labels(
        quantities[row, , drop = FALSE],
        "selector"
      )
      io_aliases <- .brma_random_parameter_io_aliases(
        catalog     = catalog,
        quantity_id = quantity[["quantity_id"]]
      )
      add_entry(
        quantity          = quantity,
        parameter         = parameter,
        component         = "random",
        term              = parameter,
        source            = "random",
        formula_parameter = quantity[["formula_parameter"]],
        entry_aliases     = c(parameter, io_aliases[["alias"]]),
        entry_alias_simplified = c(FALSE, io_aliases[["simplified"]])
      )
    }
  }

  bias_parameters <- .brma_parameter_catalog_bias_parameters(object)
  for (parameter in bias_parameters) {
    quantity_rows <- which(
      public & quantities[["canonical_name"]] == parameter
    )
    if (length(quantity_rows) > 1L) {
      .stop_refit_required(
        "Fitted publication-bias metadata for '", parameter,
        "' are ambiguous. Refit the model with the current ",
        "RoBMA/BayesTools build."
      )
    }
    if (length(quantity_rows) == 1L) {
      quantity <- quantities[quantity_rows, , drop = FALSE]
    } else {
      quantity <- .brma_parameter_catalog_bias_quantity(
        catalog   = catalog,
        parameter = parameter
      )
      extension_quantities[[length(extension_quantities) + 1L]] <- quantity
    }
    add_entry(
      quantity          = as.list(quantity),
      parameter         = parameter,
      component         = "bias",
      term              = parameter,
      source            = "bias",
      formula_parameter = NA_character_,
      entry_aliases     = switch(
        parameter,
        omega = c("omega", "weightfunction"),
        parameter
      )
    )
  }

  entries <- do.call(rbind, entries)
  rownames(entries) <- NULL
  extension_aliases <- do.call(rbind, aliases)
  extension_quantities <- if (length(extension_quantities) == 0L) {
    catalog[["quantities"]][FALSE, , drop = FALSE]
  } else {
    do.call(rbind, extension_quantities)
  }
  random_quantity_ids <- entries[["quantity_id"]][
    entries[["component"]] == "random"
  ]
  catalog[["aliases"]] <- catalog[["aliases"]][
    !catalog[["aliases"]][["quantity_id"]] %in% random_quantity_ids,
    ,
    drop = FALSE
  ]
  catalog <- BayesTools::parameter_catalog_extend(
    catalog    = catalog,
    quantities = extension_quantities,
    aliases    = extension_aliases,
    provider   = "RoBMA"
  )

  return(list(catalog = catalog, entries = entries))
}


.brma_parameter_catalog_entries_for_quantities <- function(entries,
                                                            quantity_ids) {

  direct <- entries[["quantity_id"]] %in% quantity_ids
  grouped <- vapply(
    entries[["member_quantity_ids"]],
    function(members) any(members %in% quantity_ids),
    logical(1)
  )
  return(entries[direct | grouped, , drop = FALSE])
}


.brma_parameter_catalog_group_ambiguity <- function(entries, quantity_ids,
                                                     component) {

  candidates <- .brma_parameter_catalog_entries_for_quantities(
    entries      = entries,
    quantity_ids = quantity_ids
  )
  grouped <- candidates[
    candidates[["role"]] == "formula_coefficient_group",
    ,
    drop = FALSE
  ]
  retry_component <- if (identical(component, "auto") &&
                         nrow(candidates) == 1L && nrow(grouped) == 1L) {
    grouped[["component"]]
  } else {
    NULL
  }

  return(list(candidates = candidates, component = retry_component))
}


.brma_parameter_catalog_formula_quantity <- function(
    catalog, map_row, coordinates, semantic_component) {

  fitted_scale  <- unique(coordinates[["fitted_scale"]])
  display_scale <- unique(coordinates[["display_scale"]])
  if (length(fitted_scale) != 1L || length(display_scale) != 1L) {
    .stop_refit_required(
      "Grouped formula coefficient coordinates have inconsistent scale ",
      "metadata. Refit the model with the current RoBMA/BayesTools build."
    )
  }
  dependencies <- unique(unlist(lapply(
    coordinates[["extraction_key"]],
    function(key) key[["dependencies"]]
  ), use.names = FALSE))
  if (!is.character(dependencies) || length(dependencies) == 0L ||
      anyNA(dependencies) || any(!nzchar(dependencies))) {
    .stop_refit_required(
      "Grouped formula coefficient cells have invalid BayesTools extraction ",
      "metadata. Refit the model with the current RoBMA/BayesTools build."
    )
  }

  out <- data.frame(
    quantity_id       = paste0("RoBMA::formula_", map_row[["encoded_name"]]),
    canonical_name    = map_row[["jags_name"]],
    provider          = "RoBMA",
    namespace         = map_row[["formula_parameter"]],
    role              = "formula_coefficient_group",
    formula_parameter = map_row[["formula_parameter"]],
    owner_type        = "",
    owner_name        = "",
    quantity          = "",
    scale_role        = "",
    parent_quantity_id = "",
    term              = map_row[["term"]],
    component         = semantic_component,
    display_label     = map_row[["term"]],
    fitted_scale      = fitted_scale,
    display_scale     = display_scale,
    status            = "derived",
    fixed_value       = NA_real_,
    internal          = FALSE,
    source_type       = "composite",
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
  out[["arguments"]] <- I(list(character()))
  out[["extraction_key"]] <- I(list(list(
    type              = "robma_formula_group",
    dependencies      = dependencies,
    formula_parameter = map_row[["formula_parameter"]],
    term              = map_row[["term"]]
  )))
  # A coefficient group has no scalar support of its own; its display label
  # is RoBMA's term label.
  out[["support"]]     <- I(list(NULL))
  out[["definedness"]] <- "always"
  out[["label_parts"]] <- I(list(NULL))
  out <- out[, names(catalog[["quantities"]]), drop = FALSE]
  return(out)
}


.brma_parameter_catalog_bias_parameters <- function(object) {

  prior <- object[["priors"]][["outcome"]][["bias"]]
  if (is.null(prior)) {
    return(character())
  }
  priors <- if (BayesTools::is.prior.mixture(prior)) unclass(prior) else list(prior)
  out <- "bias"
  if (any(vapply(priors, .prior_has_selection, logical(1)))) {
    out <- c(out, "omega")
  }
  if (any(vapply(priors, BayesTools::is.prior.PET, logical(1)))) {
    out <- c(out, "PET")
  }
  if (any(vapply(priors, BayesTools::is.prior.PEESE, logical(1)))) {
    out <- c(out, "PEESE")
  }
  return(out)
}


.brma_parameter_catalog_bias_quantity <- function(catalog, parameter) {

  out <- data.frame(
    quantity_id       = paste0("RoBMA::bias_", tolower(parameter)),
    canonical_name    = parameter,
    provider          = "RoBMA",
    namespace         = "bias",
    role              = "publication_bias_component",
    formula_parameter = "",
    owner_type        = "",
    owner_name        = "",
    quantity          = "",
    scale_role        = "",
    parent_quantity_id = "",
    term              = parameter,
    component         = "bias",
    display_label     = parameter,
    fitted_scale      = "fitted_original",
    display_scale     = "original",
    status            = "derived",
    fixed_value       = NA_real_,
    internal          = FALSE,
    source_type       = "composite",
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
  out[["arguments"]] <- I(list(character()))
  out[["extraction_key"]] <- I(list(list(
    type         = "robma_bias_component",
    dependencies = "bias",
    parameter    = parameter
  )))
  # The publication-bias component mixes several priors; its support is not
  # declared here, and its display label is the component name.
  out[["support"]]     <- I(list(NULL))
  out[["definedness"]] <- "always"
  out[["label_parts"]] <- I(list(NULL))
  out <- out[, names(catalog[["quantities"]]), drop = FALSE]
  return(out)
}


.brma_parameter_catalog_unfitted <- function(object) {

  rows <- list()
  outcome_priors <- object[["priors"]][["outcome"]]
  if (is.null(outcome_priors)) {
    outcome_priors <- list()
  }
  add <- function(parameter, component, term, aliases,
                 source, formula_parameter = NA_character_) {
    aliases <- unique(aliases[!is.na(aliases) & nzchar(aliases)])
    if (length(aliases) == 0L) {
      return(invisible(NULL))
    }
    for (alias in aliases) {
      rows[[length(rows) + 1L]] <<- data.frame(
        alias             = alias,
        parameter         = parameter,
        component         = component,
        term              = term,
        source            = source,
        formula_parameter = formula_parameter,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    return(invisible(NULL))
  }

  if (.is_mods(object) || .is_random(object)) {
    location_source <- if (.is_random(object)) "location" else "mods"
    if (.fitted_formula_has_intercept(object, "mu", required = FALSE) &&
        !.location_omit_fixed_zero_intercept(object)) {
      add("mu_intercept", "mods", "intercept",
          c("mu_intercept", "mu", "intercept"),
          source = location_source, formula_parameter = "mu")
    }
    .brma_parameter_catalog_terms(
      object            = object,
      model_parameter   = "mu",
      component         = "mods",
      formula_parameter = "mu",
      source            = location_source,
      add               = add
    )
  } else {
    add("mu", "mods", "mu", c("mu", "effect"), source = "outcome")
  }

  if (.is_scale(object)) {
    scale_specs <- .data_scale_component_specs(object[["data"]])
    single_scale <- length(scale_specs) == 1L
    for (scale_spec in scale_specs) {
      formula_parameter <- scale_spec[["parameter"]]
      if (single_scale && identical(formula_parameter, "log_tau")) {
        intercept_parameter <- "log_tau_intercept"
        intercept_aliases   <- c("log_tau_intercept", "tau", "intercept")
      } else {
        component_aliases <- scale_spec[["aliases"]]
        intercept_parameter <- paste0(formula_parameter, "_intercept")
        intercept_aliases   <- c(
          intercept_parameter,
          scale_spec[["display_name"]],
          component_aliases,
          paste0(component_aliases, "_intercept")
        )
      }
      if (.fitted_formula_has_intercept(
          object, formula_parameter, required = FALSE)) {
        add(intercept_parameter, "scale", "intercept",
            intercept_aliases, source = "scale",
            formula_parameter = formula_parameter)
      }
      .brma_parameter_catalog_terms(
        object            = object,
        model_parameter   = formula_parameter,
        component         = "scale",
        formula_parameter = formula_parameter,
        source            = "scale",
        add               = add
      )
    }
  } else {
    if (!is.null(outcome_priors[["tau"]])) {
      add("tau", "scale", "tau", c("tau", "heterogeneity"), source = "outcome")
    }
  }

  if (.is_multilevel(object) && !is.null(outcome_priors[["rho"]])) {
    add("rho", "scale", "rho", "rho", source = "outcome")
  }
  if (.is_PET(object)) {
    add("PET", "bias", "PET", "PET", source = "bias")
  }
  if (.is_PEESE(object)) {
    add("PEESE", "bias", "PEESE", "PEESE", source = "bias")
  }
  if (.is_weightfunction(object)) {
    add("omega", "bias", "omega", c("omega", "weightfunction"), source = "bias")
  }

  if (!is.null(outcome_priors[["bias"]])) {
    add("bias", "bias", "bias", "bias", source = "bias")
  }

  if (.is_random(object)) {
    random_specs <- .brma_random_parameter_bundle(object)[["specs"]]
    if (nrow(random_specs) > 0L) {
      for (i in seq_len(nrow(random_specs))) {
        spec <- random_specs[i, , drop = FALSE]
        add(
          parameter         = spec[["parameter"]],
          component         = "random",
          term              = spec[["label"]],
          aliases           = c(
            spec[["parameter"]],
            spec[["label"]]
          ),
          source            = "random",
          formula_parameter = spec[["formula_parameter"]]
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  out <- unique(out)
  rownames(out) <- NULL
  return(out)
}

.brma_parameter_catalog_terms <- function(object, model_parameter, component,
                                          formula_parameter, source, add) {

  if (!is.null(object[["fit"]])) {
    name_map <- .fitted_formula_name_map(
      object    = object,
      parameter = model_parameter,
      required  = TRUE
    )
    rows <- name_map[
      name_map[["kind"]] == "fixed" &
        name_map[["role"]] == "coefficient" &
        name_map[["term"]] != "intercept",
      ,
      drop = FALSE
    ]
    for (i in seq_len(nrow(rows))) {
      add(
        parameter         = rows[["jags_name"]][[i]],
        component         = component,
        term              = rows[["term"]][[i]],
        aliases           = c(rows[["jags_name"]][[i]], rows[["term"]][[i]]),
        source            = source,
        formula_parameter = formula_parameter
      )
    }
    return(invisible(NULL))
  }

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
    add(parameter, component, term, c(parameter, term),
        source = source, formula_parameter = formula_parameter)
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
    "lwd", "lty", "col", "col.fill", "xlab", "ylab", "ylab2", "main",
    "xlim", "ylim", "ylim2", "par_name", "legend", "legend_title",
    "legend_labels", "legend_position", "cex", "cex.axis",
    "cex.lab", "cex.main", "col.axis", "col.lab", "col.main",
    "las", "bty", "xaxs", "yaxs", "axes", "xaxt", "yaxt",
    "pch", "bg", "border", "width", "scale_y2", "color",
    "colour", "fill", "alpha", "size", "linewidth", "linetype"
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
    if (isTRUE(attr(samples[[parameter]], "orthonormal")) ||
        isTRUE(attr(samples[[parameter]], "meandif")) ||
        isTRUE(attr(samples[[parameter]], "independent"))) {
      n_levels <- length(attr(samples[[parameter]], "level_names"))
    } else if (isTRUE(attr(samples[[parameter]], "treatment"))) {
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

.brma_as_mixed_posteriors <- function(object, parameters, conditional = NULL,
                                      ...) {

  fit        <- object[["fit"]]
  prior_list <- attr(fit, "prior_list", exact = TRUE)

  if (!is.null(prior_list)) {
    protected <- unique(c(parameters, conditional))
    drop <- vapply(prior_list, function(prior) {
      prior_attributes <- names(attributes(prior))
      BayesTools::is.prior.vector(prior) &&
        any(startsWith(prior_attributes, "random_"))
    }, logical(1))
    drop[names(prior_list) %in% protected] <- FALSE

    if (any(drop)) {
      attr(fit, "prior_list") <- prior_list[!drop]
    }
  }

  BayesTools::as_mixed_posteriors(
    model       = fit,
    parameters  = parameters,
    conditional = conditional,
    ...
  )
}

# Axis label of a plotted parameter. Formula coefficients are labelled by the
# formula term of their catalog entry.
.plot_parameter_label <- function(parameter, effect_transform = NULL,
                                  entry = NULL, object = NULL) {

  if (.is_effect_location_parameter(parameter)) {
    if (.effect_output_active(effect_transform)) {
      return(paste0("Effect Size (", effect_transform[["label"]], ")"))
    }

    return("Effect Size")
  }

  formula_parameter <- entry[["formula_parameter"]]
  if (!is.character(formula_parameter) || length(formula_parameter) != 1L ||
      is.na(formula_parameter)) {
    formula_parameter <- ""
  }

  if (identical(formula_parameter, "mu")) {
    label <- paste0("Effect Size: ", entry[["term"]])
    if (.effect_output_active(effect_transform)) {
      label <- paste0(label, " (", effect_transform[["label"]], ")")
    }
    return(label)
  }

  if (parameter %in% c("tau", "log_tau_intercept")) {
    label <- "Heterogeneity"
    if (.effect_output_active(effect_transform)) {
      label <- paste0(label, " (", effect_transform[["label"]], ")")
    }
    return(label)
  }

  if (nzchar(formula_parameter) && !is.null(object) &&
      formula_parameter %in% .summary_scale_formula_parameters(object)) {
    label <- paste0(
      "Heterogeneity: ",
      .summary_formula_term_label(object, formula_parameter, entry[["term"]])
    )
    if (.effect_output_active(effect_transform)) {
      label <- paste0(label, " (", effect_transform[["label"]], ")")
    }
    return(label)
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

  if (!is.null(parameter_mods)) {
    .parameter_component_check_compatible(component, "mods", "parameter_mods")
    .check_plot_prior_name(parameter_mods, "parameter_mods")
    entry <- .brma_parameter_select_entry(
      object    = object,
      parameter = parameter_mods,
      component = "mods",
      argument  = "parameter_mods"
    )
    return(.select_plot_prior_entry(object, entry, allow_mixed_bias))
  }

  if (!is.null(parameter_scale)) {
    .parameter_component_check_compatible(component, "scale", "parameter_scale")
    .check_plot_prior_name(parameter_scale, "parameter_scale")
    entry <- .brma_parameter_select_entry(
      object    = object,
      parameter = parameter_scale,
      component = "scale",
      argument  = "parameter_scale"
    )
    return(.select_plot_prior_entry(object, entry, allow_mixed_bias))
  }

  if (is.null(parameter)) {
    parameter <- .brma_parameter_default(component, object)
  }

  .check_plot_prior_name(parameter, "parameter")

  entry <- .brma_parameter_select_entry(
    object    = object,
    parameter = parameter,
    component = component,
    argument  = "parameter"
  )

  return(.select_plot_prior_entry(object, entry, allow_mixed_bias))
}

.select_print_prior_parameter <- function(
    object, parameter, parameter_mods, parameter_scale,
    component = "auto", allow_mixed_bias = FALSE) {

  component <- .parameter_component_normalize(component)
  if (is.null(parameter_mods) && is.null(parameter_scale) &&
      !is.null(parameter) && identical(parameter, "random")) {
    catalog <- .brma_parameter_catalog(object)
    matches <- catalog[catalog[["alias"]] == parameter, , drop = FALSE]
    if (!identical(component, "auto")) {
      matches <- matches[matches[["component"]] == component, , drop = FALSE]
    }
    if (nrow(matches) == 0L && identical(component, "auto")) {
      return(.select_print_prior_random(object))
    }
  }

  .select_plot_prior_parameter(
    object           = object,
    parameter        = parameter,
    parameter_mods   = parameter_mods,
    parameter_scale  = parameter_scale,
    component        = component,
    allow_mixed_bias = allow_mixed_bias
  )
}

.select_plot_prior_entry <- function(object, entry, allow_mixed_bias = FALSE) {

  priors            <- object[["priors"]]
  parameter         <- entry[["parameter"]]
  source            <- entry[["source"]]
  term              <- entry[["term"]]
  formula_parameter <- entry[["formula_parameter"]]

  if (identical(source, "outcome")) {
    prior <- priors[["outcome"]][[parameter]]
    if (is.null(prior)) {
      stop(sprintf("Unknown outcome prior parameter '%s'.", parameter), call. = FALSE)
    }
    return(list(
      prior             = prior,
      label             = parameter,
      source            = source,
      term              = term,
      formula_parameter = formula_parameter
    ))
  }

  if (source %in% c("mods", "location")) {
    return(.select_plot_prior_term(
      prior_list        = priors[[source]],
      term              = term,
      argument          = "parameter",
      prefix            = "mu",
      source            = source,
      label             = parameter,
      formula_parameter = formula_parameter
    ))
  }

  if (identical(source, "scale")) {
    scale_priors <- priors[["scale"]]
    if (inherits(scale_priors, "RoBMA_scale_priors")) {
      prior_list <- scale_priors[[formula_parameter]]
    } else {
      prior_list <- scale_priors
    }
    return(.select_plot_prior_term(
      prior_list        = prior_list,
      term              = term,
      argument          = "parameter",
      prefix            = formula_parameter,
      source            = source,
      label             = parameter,
      formula_parameter = formula_parameter
    ))
  }

  if (identical(source, "bias")) {
    prior <- .select_plot_prior_bias(
      prior            = priors[["outcome"]][["bias"]],
      parameter        = parameter,
      allow_mixed_bias = allow_mixed_bias
    )
    return(list(
      prior             = prior,
      label             = parameter,
      source            = source,
      term              = term,
      formula_parameter = formula_parameter
    ))
  }

  stop(sprintf("Unknown prior source '%s'.", source), call. = FALSE)
}

.select_plot_prior_term <- function(prior_list, term, argument, prefix, source,
                                    label = NULL,
                                    formula_parameter = NA_character_) {

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
    prior             = prior_list[[term]],
    label             = if (!is.null(label)) label else paste0(prefix, "_", term),
    source            = source,
    term              = term,
    formula_parameter = formula_parameter
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

.select_print_prior_all <- function(object) {

  catalog  <- .brma_parameter_catalog(object)
  catalog  <- catalog[catalog[["source"]] != "random", , drop = FALSE]
  catalog  <- catalog[!duplicated(catalog[["parameter"]]), , drop = FALSE]
  if (any(catalog[["source"]] == "bias" & catalog[["parameter"]] == "bias")) {
    catalog <- catalog[
      !(catalog[["source"]] == "bias" & catalog[["parameter"]] != "bias"),
      ,
      drop = FALSE
    ]
  }
  selected <- list()

  for (i in seq_len(nrow(catalog))) {
    entry          <- as.list(catalog[i, , drop = FALSE])
    selected_entry <- .select_plot_prior_entry(
      object           = object,
      entry            = entry,
      allow_mixed_bias = TRUE
    )
    if (is.null(selected_entry[["prior"]])) {
      next
    }
    selected[[selected_entry[["label"]]]] <- selected_entry
  }

  if (.has_print_prior_random(object)) {
    selected[["random"]] <- .select_print_prior_random(object)
  }

  return(selected)
}

.has_print_prior_random <- function(object) {

  inherits(object[["priors"]][["random"]], "prior_random")
}

.select_print_prior_random <- function(object) {

  if (!.has_print_prior_random(object)) {
    stop(
      "The specified parameter 'random' is not available. ",
      "The object does not contain random-effect priors.",
      call. = FALSE
    )
  }

  list(
    prior             = object[["priors"]][["random"]],
    label             = "random",
    source            = "random",
    term              = "random",
    formula_parameter = NA_character_
  )
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

  if (!(selected[["source"]] %in% c("mods", "location", "scale"))) {
    return(NULL)
  }

  if (!isTRUE(.standardize_continuous_predictors(object))) {
    return(NULL)
  }

  formula_info <- .plot_prior_formula_info(
    object   = object,
    selected = selected
  )

  if (is.null(formula_info) || is.null(formula_info[["formula_scale"]])) {
    return(NULL)
  }

  # The formula columns are named as BayesTools names formula terms.
  parameter <- BayesTools::JAGS_parameter_names(selected[["label"]])
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

.plot_prior_formula_info <- function(object, selected) {

  source    <- selected[["source"]]
  parameter <- selected[["formula_parameter"]]
  if (is.null(parameter) || is.na(parameter) || !nzchar(parameter)) {
    parameter <- .fitted_formula_parameter(source)
  }
  design    <- .fitted_formula_design(
    object    = object,
    parameter = parameter,
    required  = TRUE
  )

  if (is.null(design[["formula_scale"]])) {
    return(NULL)
  }

  # The formula scaling that carries the fitted design transforms the
  # coefficients to the original predictor scale.
  formula_scale <- .object_formula_scale(object, parameter)

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
