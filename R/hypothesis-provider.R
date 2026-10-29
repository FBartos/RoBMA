# ============================================================================ #
# hypothesis-provider.R
# ============================================================================ #
#
# The providers that describe RoBMA models to the hypothesis engine of
# BayesTools (BayesTools::hypothesis_BF_model()): the fitted brma objects and
# the stored marginal means of marginal_means(). A provider starts from the
# default provider of the BayesTools fit of the model and replaces what is
# specific to RoBMA: the parameter catalog and its component vocabulary, the
# draws of the tested quantities (the product-space ensemble), the
# random-effect quantities (R/random-parameters.R), the wording of refusals, and
# the density backends that evaluate point hypotheses (hypothesis-backends.R).
# The engine owns planning, evaluation, and binding. See
# BayesTools/.agents/instructions/hypothesis.md.
#
# ============================================================================ #


#' @noRd
#' @exportS3Method BayesTools::hypothesis_provider
hypothesis_provider.brma <- function(object, ...) {

  if (...length() > 0L) {
    stop("Unused argument(s) in hypothesis_provider(): ",
         paste0("'", names(list(...)), "'", collapse = ", "), ".",
         call. = FALSE)
  }
  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0L) {
    stop("'hypothesis' requires a fitted brma object.", call. = FALSE)
  }
  .hypothesis_density_backends_ensure()
  # A fit without the current BayesTools metadata is a request to refit.
  .brma_parameter_catalog_metadata(object)
  fit      <- object[["fit"]]
  provider <- utils::modifyList(BayesTools::hypothesis_provider(fit), list(
    build                    = "RoBMA/BayesTools",
    refit_class              = "RoBMA_refit_required",
    caller                   = "hypothesis()",
    ordinate_label           = "qCMDE/IWMDE",
    transform_scaled_default = TRUE,
    random_point_mass_advice = paste0(
      "Use a region or directional hypothesis; the Component Inclusion table ",
      "from 'summary(object)' or 'summary_models(object)' compares the ",
      "exclusion and inclusion of gated components."
    ),
    density_methods          = c("KDE", "normal", .hypothesis_density_backend_names()),
    components               = list(
      values   = c("auto", "mods", "location", "scale", "bias", "random"),
      synonyms = c(location = "mods"),
      default  = "auto",
      refused  = c(bias = .hypothesis_brma_bias_reason()),
      random   = "random",
      advice   = "Set 'component' to 'mods'/'location', 'scale', or 'random'."
    ),
    metadata                 = function() .brma_parameter_catalog_metadata(object),
    uncovered                = function(resolution, quantity_ids, symbol,
                                        metadata) {
      .hypothesis_brma_uncovered(
        object, metadata, quantity_ids, symbol
      )
    },
    sample_parameters        = function(parameter) {
      .as_mixed_posteriors_parameters(object, parameter)
    },
    mixed_posteriors         = function(parameters, conditional,
                                        conditional_rule, transform_scaled,
                                        n_prior_samples) {
      .brma_as_mixed_posteriors(
        object           = object,
        parameters       = parameters,
        conditional      = conditional,
        conditional_rule = conditional_rule,
        transform_scaled = transform_scaled,
        n_prior_samples  = n_prior_samples
      )
    },
    density_parameter        = function(samples, parameter, sample_parameters) {
      .plot_brma_density_sample_parameter(
        samples          = samples,
        parameter        = parameter,
        sample_parameter = sample_parameters
      )
    },
    formula_name_map         = function(formula_parameter) {
      .fitted_formula_name_map(object, formula_parameter)
    },
    random                   = .hypothesis_brma_random_callbacks(object),
    full_ensemble            = function(parameter) {
      .hypothesis_brma_full_ensemble(object, parameter)
    },
    combination_hint         = function(plan) {
      if (!.is_RoBMA(object) || isTRUE(plan[["conditional"]])) {
        return(NULL)
      }
      paste0(
        " The levels of '", plan[["label"]], "' share the null ",
        "component of the model-averaged prior; test the combination within ",
        "the models that include the term with conditional = TRUE."
      )
    },
    check_density_control    = function(density_method, density_control) {
      # The refusal of the settings of the methods that use none has the
      # wording of the settings of qCMDE and IWMDE.
      if (!density_method %in% .hypothesis_density_backend_names()) {
        .density_control_normalize(
          density_method  = density_method,
          density_control = density_control,
          allow_normal    = TRUE,
          purpose         = "ordinate"
        )
      }
      invisible(NULL)
    }
  ))
  # The model objects are not merged with the defaults of the fit.
  provider[["object"]]            <- object
  provider[["fit"]]               <- fit
  provider[["capability_object"]] <- object
  provider[["point_mass_reason"]] <- if (.is_RoBMA(object)) {
    function() {
      paste0(
        "There is a point mass in the prior at the exact null hypothesis value. ",
        "The Savage-Dickey density ratio is invalid. This parameter has a null ",
        "component, so its evidence against the null is the inclusion Bayes ",
        "factor reported by 'summary()' and 'summary_models()'."
      )
    }
  }

  provider
}


# Whether omitting 'conditional' tests the model-averaged 'parameter' within
# the full ensemble, including the models in which it is fixed by its null
# component: a RoBMA ensemble whose prior of the parameter has a null and an
# alternative component.
.hypothesis_brma_full_ensemble <- function(object, parameter) {

  prior_list <- attr(object[["fit"]], "prior_list", exact = TRUE)
  if (!.is_RoBMA(object) || !parameter %in% names(prior_list) ||
      !BayesTools::is.prior.mixture(prior_list[[parameter]])) {
    return(FALSE)
  }
  prior_components <- attr(prior_list[[parameter]], "components", exact = TRUE)

  any(prior_components == "null") && any(prior_components == "alternative")
}


# The callbacks of the random-effect quantities (R/random-parameters.R);
# 'transform_scaled' is the opposite of the 'standardized_coefficients' of
# hypothesis().
.hypothesis_brma_random_callbacks <- function(object) {

  list(
    select           = function(parameter, transform_scaled) {
      .brma_random_parameter_select(
        object                    = object,
        parameter                 = parameter,
        standardized_coefficients = !transform_scaled
      )
    },
    mixed_posterior  = function(parameter, transform_scaled, conditional,
                                selected) {
      .brma_random_parameter_mixed_posterior(
        object                    = object,
        parameter                 = parameter,
        standardized_coefficients = !transform_scaled,
        conditional               = conditional,
        selected                  = selected
      )
    },
    support          = function(selected) {
      .brma_random_parameter_support(selected)
    },
    density_target   = function(parameter, transform_scaled, operation) {
      .brma_random_parameter_density_target(
        object,
        parameter,
        operation                 = operation,
        standardized_coefficients = !transform_scaled
      )
    },
    sd_pair          = function(selected) {
      .brma_random_parameter_sd_pair(object, selected)
    },
    zero_boundary    = function(selected) {
      .brma_random_parameter_zero_boundary(object, selected)
    },
    defined_share    = function(values, selected, conditional) {
      .brma_random_parameter_defined_share(
        object, values, selected = selected, conditional = conditional
      )
    },
    defined_footnote = .brma_random_parameter_defined_footnote,
    defined_draws    = .brma_random_parameter_defined_draws,
    prior_draws      = function(parameter, selected, transform_scaled,
                                n_samples, seed) {
      .brma_random_parameter_select(
        object                    = object,
        parameter                 = parameter,
        standardized_coefficients = !transform_scaled,
        prior                     = TRUE,
        n_prior_samples           = n_samples,
        seed                      = seed
      )
    }
  )
}


# Resolved catalog quantities without a RoBMA parameter entry ('quantity_ids')
# are no hypothesis targets of a current fit, so a refit cannot help: contrast
# coefficients of formula factor terms are statements to restate on factor
# levels, quantities of the publication-bias prior have the publication-bias
# refusal, and other quantities (e.g. inclusion indicators, 'inclusion(<component>)'
# of variance allocations, and latent cluster effects) are refused as targets,
# named as the statement references them (.brma_stop_uncovered_quantity(),
# which the parameter selection of plot() shares). The refusal is returned for
# the engine to raise.
.hypothesis_brma_uncovered <- function(object, metadata, quantity_ids, symbol) {

  tryCatch(
    .brma_stop_uncovered_quantity(
      object       = object,
      metadata     = metadata,
      quantity_ids = quantity_ids,
      symbol       = symbol,
      hypothesis   = TRUE
    ),
    BayesTools_hypothesis_statement   = .hypothesis_brma_refusal_of,
    BayesTools_hypothesis_unavailable = .hypothesis_brma_refusal_of
  )
}


# The refusal (reason, complete classes, and further fields) a condition raised
# by BayesTools::hypothesis_stop() states.
.hypothesis_brma_refusal_of <- function(condition) {

  fields <- unclass(condition)
  list(
    reason = conditionMessage(condition),
    class  = setdiff(class(condition), c("error", "condition")),
    fields = fields[setdiff(names(fields), c("message", "call"))]
  )
}


# ---------------------------------------------------------------------------- #
# Marginal means
# ---------------------------------------------------------------------------- #

#' @noRd
#' @exportS3Method BayesTools::hypothesis_provider
hypothesis_provider.marginal_means.brma <- function(object, ...) {

  if (...length() > 0L) {
    stop("Unused argument(s) in hypothesis_provider(): ",
         paste0("'", names(list(...)), "'", collapse = ", "), ".",
         call. = FALSE)
  }
  .hypothesis_density_backends_ensure()
  source_object <- object[["source_object"]]
  source_valid  <- !is.null(source_object) && inherits(source_object, "brma") &&
    !is.null(source_object[["fit"]])

  provider <- list(
    kind                  = "marginal_means",
    object                = object,
    build                 = "RoBMA/BayesTools",
    refit_class           = "RoBMA_refit_required",
    caller                = "hypothesis()",
    ordinate_label        = "qCMDE/IWMDE",
    density_methods       = c("KDE", .hypothesis_density_backend_names()),
    inference             = object[["inference"]],
    model_averaged        = isTRUE(object[["model_averaged"]]) ||
      inherits(source_object, "RoBMA"),
    capability_object     = if (source_valid) source_object,
    source_missing_reason = paste0(
      "The marginal-means object does not contain the source fitted brma ",
      "object needed to compute qCMDE/IWMDE ordinates."
    ),
    select                = function(hypothesis, parameter) {
      .hypothesis_marginal_means_select_parameter(object, hypothesis, parameter)
    },
    density_method_default = function() {
      .marginal_means_density_method(object, NULL)
    },
    density_control_default = function(density_method) {
      .hypothesis_density_control_fill(
        .marginal_means_iwmde_settings_control(object, density_method, "ordinate")
      )
    },
    check_density_control = function(density_method, density_control) {
      if (!is.null(density_control) &&
          !density_method %in% .hypothesis_density_backend_names()) {
        stop("'density_control' is only used when 'density_method' is ",
             "'qCMDE' or 'IWMDE'.", call. = FALSE)
      }
      invisible(NULL)
    },
    metadata              = function() list(term_map = object[["term_map"]])
  )
  class(provider) <- c("BayesTools_hypothesis_provider", "list")

  provider
}


# The marginal-means parameter that a statement of the marginal means of
# 'object' refers to, and the aliases that name it ('hypothesis' is a parsed
# statement).
.hypothesis_marginal_means_select_parameter <- function(object, hypothesis,
                                                         parameter) {

  # References that name no marginal-means parameter are statement errors
  # ("BayesTools_hypothesis_statement"), with the classes with which the
  # fitted-model path refuses them: a statement without parameter symbols
  # (also with 'parameter') has the classes of BayesTools' refusal
  # "BayesTools_hypothesis_no_parameters" and its parent
  # "BayesTools_parameter_resolution_error"; an unknown name has those of
  # BayesTools' "BayesTools_parameter_not_found".
  roots <- BayesTools::hypothesis_symbols(hypothesis)
  roots <- unique(roots[nzchar(roots)])
  if (length(roots) == 0L) {
    BayesTools::hypothesis_stop(
      "Hypothesis must reference a marginal-means parameter.",
      type  = "statement",
      class = c("BayesTools_hypothesis_no_parameters",
                "BayesTools_parameter_resolution_error")
    )
  }

  alias_catalog <- .hypothesis_marginal_means_alias_catalog(object)
  if (!is.null(parameter)) {
    selected <- .marginal_means_select_parameter(object, parameter)[["parameter"]]
    .hypothesis_marginal_means_check_parameter(hypothesis, alias_catalog, selected)
    aliases <- .hypothesis_marginal_means_aliases(alias_catalog, selected)
    return(list(parameter = selected, aliases = aliases))
  }

  matches <- lapply(roots, function(root) {

    unique(alias_catalog[["parameter"]][alias_catalog[["alias"]] == root])
  })
  # An alias naming several parameters is ambiguous, as for fitted objects;
  # 'parameter' selects one of them.
  ambiguous <- lengths(matches) > 1L
  if (any(ambiguous)) {
    root <- roots[ambiguous][[1L]]
    BayesTools::hypothesis_stop(
      paste0(
        "Marginal-means alias '", root, "' is ambiguous (",
        paste(matches[[which(ambiguous)[[1L]]]], collapse = ", "),
        "). Specify 'parameter' using an internal parameter name."
      ),
      type = "ambiguous"
    )
  }
  known <- unlist(matches[lengths(matches) == 1L], use.names = FALSE)
  known <- unique(known)

  if (length(known) == 1L) {
    aliases <- .hypothesis_marginal_means_aliases(alias_catalog, known)
    return(list(parameter = known, aliases = aliases))
  }

  # A statement is tested on the marginal means of one parameter: hypotheses
  # that name several parameters are tested by statement (the method resolves
  # each statement on its own), and 'parameter' cannot resolve a statement
  # that names several.
  if (length(known) > 1L) {
    BayesTools::hypothesis_stop(
      .hypothesis_marginal_means_multiple_parameters_message(known),
      type = "ambiguous"
    )
  }

  # The fields of BayesTools' condition: the unresolved alias and the
  # available ones.
  available <- sort(unique(alias_catalog[["alias"]]))
  stop(structure(
    class = c("BayesTools_hypothesis_statement",
              "BayesTools_parameter_not_found",
              "BayesTools_parameter_resolution_error", "error", "condition"),
    list(
      message   = paste0(
        "Could not infer a marginal-means parameter from the hypothesis. ",
        if (length(roots) > 1L) {
          paste0("The names ", paste0("'", roots, "'", collapse = ", "),
                 " are not marginal-means quantities. ")
        } else {
          paste0("The name '", roots, "' is not a marginal-means quantity. ")
        },
        "Available quantities are: ",
        paste0("'", available, "'", collapse = ", "), "."
      ),
      call      = NULL,
      alias     = roots[[1L]],
      available = available
    )
  ))
}


# The message of a statement that names the marginal means of several
# parameters, with and without 'parameter'.
.hypothesis_marginal_means_multiple_parameters_message <- function(parameters) {

  paste0(
    "Hypothesis references multiple marginal-means parameters (",
    paste(parameters, collapse = ", "),
    "), but one statement must refer to the marginal means of one ",
    "parameter. Use separate statements for different parameters."
  )
}


# The refusal of a statement that names the marginal means of another
# parameter than the 'selected' one: the roots of every statement of
# 'hypothesis' that are aliases of marginal-means parameters must include the
# selected one. A root that is no alias of any parameter is an unknown name,
# which the planning refuses as at the other entry points.
.hypothesis_marginal_means_check_parameter <- function(hypothesis, alias_catalog,
                                                       selected) {

  statement <- BayesTools::hypothesis_render(hypothesis)
  roots <- BayesTools::hypothesis_symbols(hypothesis)
  roots <- unique(roots[nzchar(roots)])
  mismatch <- function(message) {
    stop(structure(
      class = c("BayesTools_hypothesis_parameter_mismatch",
                "BayesTools_hypothesis_statement", "error", "condition"),
      list(message = message, call = NULL)
    ))
  }

  # A statement that names the marginal means of several parameters has the
  # refusal it has without 'parameter' (omitting 'parameter' does not help),
  # with the class of the mismatch with the selected parameter.
  named <- lapply(roots, function(root) {
    unique(alias_catalog[["parameter"]][alias_catalog[["alias"]] == root])
  })
  named <- named[lengths(named) > 0L]
  if (length(named) > 1L && length(Reduce(intersect, named)) == 0L) {
    mismatch(.hypothesis_marginal_means_multiple_parameters_message(
      unique(unlist(named))
    ))
  }

  for (root in roots) {
    parameters <- unique(alias_catalog[["parameter"]][alias_catalog[["alias"]] == root])
    if (length(parameters) > 0L && !selected %in% parameters) {
      mismatch(paste0(
        "Statement '", statement, "' refers to the marginal means of '",
        paste(parameters, collapse = "', '"), "', but 'parameter' selects '",
        selected, "'. A statement must refer to the marginal means of one ",
        "parameter; omit 'parameter' to test statements on different ",
        "parameters in one call."
      ))
    }
  }

  invisible(NULL)
}


.hypothesis_marginal_means_alias_catalog <- function(object) {

  term_map <- object[["term_map"]]
  rows     <- list()

  for (i in seq_len(nrow(term_map))) {
    parameter <- term_map[["parameter"]][[i]]
    term      <- term_map[["term"]][[i]]

    aliases <- parameter
    if (!is.null(term) && nzchar(term)) {
      aliases <- c(aliases, term)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      alias      = unique(aliases),
      parameter  = parameter,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  intercept <- term_map[["parameter"]][term_map[["term"]] == "intercept"]
  if (length(intercept) == 1L) {
    rows[[length(rows) + 1L]] <- data.frame(
      alias      = "mu",
      parameter  = intercept,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  aliases <- unique(do.call(rbind, rows))
  rownames(aliases) <- NULL
  return(aliases)
}


.hypothesis_marginal_means_aliases <- function(alias_catalog, parameter) {

  rows <- alias_catalog[alias_catalog[["parameter"]] == parameter, , drop = FALSE]
  aliases <- as.list(rows[["parameter"]])
  names(aliases) <- rows[["alias"]]
  return(aliases)
}


# ---------------------------------------------------------------------------- #
# hypothesis_quantities() of an unfitted object
# ---------------------------------------------------------------------------- #

# The eligibility table of a fitted object is rendered by BayesTools from the
# provider (hypothesis_provider.brma()); an object without a fit lists the
# names of its catalog with the refusal that hypothesis() needs a fit.
#' @export
hypothesis_quantities.brma <- function(object, ...) {

  if (!is.null(object[["fit"]]) && length(object[["fit"]]) > 0L) {
    return(NextMethod())
  }
  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis_quantities()"
  )
  catalog <- .brma_parameter_catalog(object)
  catalog <- catalog[catalog[["component"]] != "bias", , drop = FALSE]
  out  <- catalog[, c("alias", "parameter", "component", "term"), drop = FALSE]
  out[["bracket"]]             <- NA_character_
  out[["point_test"]]            <- FALSE
  out[["direction_test"]]        <- FALSE
  out[["contrast_test"]]         <- NA
  out[["point_test_methods"]]    <- ""
  out[["contrast_test_methods"]] <- NA_character_
  out[["reason"]]                <- "'hypothesis' requires a fitted brma object."
  out <- out[, c(
    "alias", "parameter", "component", "term", "bracket", "point_test",
    "direction_test", "contrast_test", "point_test_methods",
    "contrast_test_methods", "reason"
  ), drop = FALSE]
  rownames(out) <- NULL

  return(out)
}
