# ============================================================================ #
# formula-design.R
# ============================================================================ #
#
# Shared accessors for BayesTools formula-design metadata.
#
# ============================================================================ #


# Validate the versioned BayesTools metadata consumed by a fitted-object path.
.brma_validate_fit_contract <- function(object, requires) {

  fit <- if (inherits(object, "BayesTools_fit")) object else object[["fit"]]
  if (is.null(fit) || !inherits(fit, "BayesTools_fit")) {
    stop(
      "Current BayesTools fitted metadata are unavailable. Refit the model ",
      "with the current RoBMA/BayesTools build.",
      call. = FALSE
    )
  }

  BayesTools::JAGS_validate_fit_contract(fit, requires = requires)
  return(invisible(TRUE))
}


# Return the persisted semantic-to-backend map for one fitted formula.
.fitted_formula_name_map <- function(object, parameter, required = TRUE) {

  if (is.null(object[["fit"]])) {
    if (required) {
      stop(
        "Fitted formula name-map metadata for parameter '", parameter,
        "' are unavailable. Fit the model before requesting fitted ",
        "parameter identities.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  .brma_validate_fit_contract(
    object,
    requires = c(
      "name_encoding",
      "formula_name_map",
      "formula_design",
      "parameter_map"
    )
  )
  name_map <- BayesTools::JAGS_formula_name_map(object[["fit"]], parameter)
  if (is.null(name_map) && required) {
    stop(
      "Fitted formula name-map metadata for parameter '", parameter,
      "' are missing. Refit the model with the current RoBMA/BayesTools build.",
      call. = FALSE
    )
  }

  return(name_map)
}


# Return the BayesTools formula design for a brma formula parameter.
.fitted_formula_design <- function(object, parameter, required = TRUE) {

  design <- NULL

  if (!is.null(object[["fit"]])) {
    designs <- BayesTools::JAGS_formula_design(object[["fit"]])
    design  <- designs[[parameter]]
  }

  if (is.null(design)) {
    design <- .object_formula_design(object, parameter)
  }

  if (is.null(design) && required) {
    stop(
      "Fitted formula design metadata for parameter '", parameter,
      "' is missing. Refit the model with the current RoBMA/BayesTools build.",
      call. = FALSE
    )
  }

  return(design)
}


# Return formula design metadata stored on or reconstructed from a brma object.
.object_formula_design <- function(object, parameter) {

  if (!is.null(object[["formula_design"]][[parameter]])) {
    return(object[["formula_design"]][[parameter]])
  }

  source <- .fitted_formula_source(
    parameter = parameter,
    data      = object[["data"]]
  )
  if (is.null(source) ||
      is.null(object[["data"]]) ||
      is.null(object[["data"]][[source]])) {
    return(NULL)
  }

  if (!is.null(object[["priors"]]) &&
      !is.null(object[["priors"]][[source]])) {
    return(.object_bayestools_formula_design(
      object    = object,
      parameter = parameter,
      source    = source
    ))
  }

  if (source == "location" && .is_data_random(object[["data"]])) {
    return(NULL)
  }

  return(.object_data_formula_design(
    object    = object,
    parameter = parameter,
    source    = source
  ))
}


# Reconstruct BayesTools formula design for objects that store data and priors
# but do not have a fitted JAGS object, e.g. only_priors objects.
.object_bayestools_formula_design <- function(
    object, parameter, source,
    random_effects_compile = .object_formula_random_effects_compile(object, source)) {

  if (identical(source, "scale")) {
    scale_spec <- .fitted_scale_spec(
      data      = object[["data"]],
      parameter = parameter
    )
    formula_design <- BayesTools::JAGS_formula(
      formula       = .create_fit_scale_formula(scale_spec[["formula"]]),
      parameter     = parameter,
      data          = scale_spec[["data"]],
      prior_list    = .create_fit_scale_formula_prior_list(
        priors    = object[["priors"]],
        parameter = parameter
      ),
      formula_scale = .data_standardize_continuous_predictors(object[["data"]])
    )[["formula_design"]]

    return(formula_design)
  }

  formula_design <- BayesTools::JAGS_formula(
    formula       = .create_fit_formula_list(
      data      = object[["data"]],
      parameter = source
    ),
    parameter     = parameter,
    data          = .create_fit_formula_data_list(
      data      = object[["data"]],
      parameter = source
    ),
    prior_list    = .create_fit_formula_prior_list(
      priors    = object[["priors"]],
      parameter = source
    ),
    formula_scale = .data_standardize_continuous_predictors(object[["data"]]),
    prior_random  = .object_formula_prior_random(
      object = object,
      source = source
    ),
    random_effects_compile = random_effects_compile
  )[["formula_design"]]

  return(formula_design)
}

.object_formula_prior_random <- function(object, source) {

  if (source != "location" || !.is_data_random(object[["data"]])) {
    return(NULL)
  }

  return(object[["priors"]][["random"]])
}


# Build minimal formula design metadata for data-only objects. These objects do
# not have priors, so BayesTools::JAGS_formula() cannot be used.
.object_data_formula_design <- function(object, parameter, source) {

  if (identical(source, "scale")) {
    scale_spec <- .fitted_scale_spec(
      data      = object[["data"]],
      parameter = parameter
    )
    formula <- .create_fit_scale_formula(scale_spec[["formula"]])
    data    <- scale_spec[["data"]]
  } else {
    formula <- .create_fit_formula_list(data = object[["data"]], parameter = source)
    data    <- object[["data"]][[source]]
  }

  terms <- attr(data, "terms")

  if (is.null(terms)) {
    terms <- stats::terms(formula, data = data)
  }

  model_matrix <- stats::model.matrix(terms, data = data)
  column_names <- colnames(model_matrix)
  assign       <- attr(model_matrix, "assign")
  term_labels  <- attr(terms, "term.labels")
  term_labels  <- gsub(":", "__xXx__", term_labels, fixed = TRUE)

  if (is.null(assign)) {
    assign <- seq_len(ncol(model_matrix)) - 1L
  }

  qr_x    <- qr(model_matrix)
  aliased <- rep(FALSE, ncol(model_matrix))
  if (qr_x[["rank"]] < ncol(model_matrix)) {
    aliased[qr_x[["pivot"]][seq(qr_x[["rank"]] + 1L, ncol(model_matrix))]] <- TRUE
  }
  names(aliased) <- column_names

  predictors <- all.vars(formula)
  predictors <- predictors[predictors %in% names(data)]
  predictor_types <- stats::setNames(
    vapply(data[predictors], .formula_design_predictor_type, character(1)),
    predictors
  )

  model_terms <- term_labels
  if (isTRUE(attr(terms, "intercept") == 1L)) {
    model_terms <- c("intercept", model_terms)
  }
  model_terms_type <- stats::setNames(
    vapply(model_terms, .formula_design_term_type, character(1), predictor_types = predictor_types),
    model_terms
  )

  xlevels <- lapply(data[vapply(data, is.factor, logical(1))], levels)

  design <- list(
    parameter         = parameter,
    formula           = formula,
    model_frame       = data,
    model_matrix      = model_matrix,
    column_names      = column_names,
    raw_column_names  = column_names,
    assign            = assign,
    terms             = terms,
    contrasts         = attr(model_matrix, "contrasts"),
    xlevels           = xlevels,
    predictors        = predictors,
    predictor_types   = predictor_types,
    model_terms       = model_terms,
    model_terms_type  = model_terms_type,
    prior_list        = NULL,
    formula_scale     = NULL,
    rank              = qr_x[["rank"]],
    qr_pivot          = qr_x[["pivot"]],
    aliased           = aliased,
    transformed_terms = list(),
    random_effects    = list(),
    jags_data_names   = list()
  )
  class(design) <- c("BayesTools_formula_design", "list")

  return(design)
}


# Classify a formula predictor for minimal data-only designs.
.formula_design_predictor_type <- function(x) {

  if (is.factor(x) || is.character(x)) {
    return("factor")
  }

  return("continuous")
}


# Classify a formula term from the predictor types.
.formula_design_term_type <- function(term, predictor_types) {

  if (identical(term, "intercept")) {
    return("continuous")
  }

  components <- unlist(strsplit(term, "__xXx__", fixed = TRUE), use.names = FALSE)
  if (any(predictor_types[components] == "factor", na.rm = TRUE)) {
    return("factor")
  }

  return("continuous")
}


# Map a BayesTools formula parameter to the RoBMA data/prior source.
.fitted_formula_source <- function(parameter, data = NULL) {

  if (identical(parameter, "mu")) {
    return(if (!is.null(data) && .is_data_random(data)) {
      "location"
    } else {
      "mods"
    })
  }
  if (identical(parameter, "log_tau")) {
    return("scale")
  }
  if (!is.null(data) && .is_data_scale(data) &&
      parameter %in% .data_scale_formula_parameters(data)) {
    return("scale")
  }

  NULL
}


# Map a RoBMA formula source name to the BayesTools formula parameter.
.fitted_formula_parameter <- function(source) {

  switch(
    source,
    "mods"     = "mu",
    "location" = "mu",
    "scale"    = "log_tau",
    source
  )
}

.fitted_scale_spec <- function(data, parameter) {

  specs <- .data_scale_component_specs(data)
  for (scale_spec in specs) {
    if (identical(scale_spec[["parameter"]], parameter)) {
      return(scale_spec)
    }
  }

  stop(
    "Scale formula design metadata for parameter '", parameter,
    "' is missing.",
    call. = FALSE
  )
}


# Return fitted model terms, optionally converted to user-facing interaction names.
.fitted_formula_terms <- function(object, parameter, include_intercept = TRUE,
                                  display = TRUE, required = TRUE) {

  design <- .fitted_formula_design(
    object    = object,
    parameter = parameter,
    required  = required
  )

  if (is.null(design)) {
    return(NULL)
  }

  terms <- design[["model_terms"]]
  if (!include_intercept) {
    terms <- terms[terms != "intercept"]
  }
  if (display) {
    terms <- .formula_design_display_names(terms)
  }

  return(terms)
}


# Return whether the user-facing fitted formula contains an intercept.
.fitted_formula_has_intercept <- function(object, parameter, required = TRUE) {

  design <- .fitted_formula_design(
    object    = object,
    parameter = parameter,
    required  = required
  )
  if (is.null(design)) {
    return(FALSE)
  }

  source_terms <- attr(design[["source_data"]], "terms", exact = TRUE)
  if (!is.null(source_terms)) {
    return(identical(attr(source_terms, "intercept", exact = TRUE), 1L))
  }

  "intercept" %in% design[["model_terms"]]
}


# Return whether a fixed-zero location intercept should be omitted publicly.
# Keep an intercept-only model visible even when its prior is a point at zero.
.location_omit_fixed_zero_intercept <- function(object) {

  if (!.is_mods(object)) {
    return(FALSE)
  }

  moderator_terms <- .fitted_formula_terms(
    object            = object,
    parameter         = "mu",
    include_intercept = FALSE,
    display           = FALSE,
    required          = FALSE
  )
  if (length(moderator_terms) == 0L) {
    return(FALSE)
  }

  source          <- if (.is_random(object)) "location" else "mods"
  intercept_prior <- object[["priors"]][[source]][["intercept"]]
  if (is.null(intercept_prior) && !is.null(object[["fit"]])) {
    prior_list      <- attr(object[["fit"]], "prior_list", exact = TRUE)
    intercept_prior <- prior_list[["mu_intercept"]]
  }

  !is.null(intercept_prior) &&
    BayesTools::is.prior.point(intercept_prior) &&
    identical(as.numeric(mean(intercept_prior)), 0)
}


# Return a fitted formula with an evaluation environment suitable for model.frame().
.fitted_formula_evaluable <- function(object, design, source) {

  formula <- design[["formula"]]
  if (is.null(formula)) {
    return(NULL)
  }

  stored_formula <- attr(object[["data"]][[source]], "formula")
  if (!is.null(stored_formula) && !is.null(environment(stored_formula))) {
    environment(formula) <- environment(stored_formula)
  } else {
    environment(formula) <- baseenv()
  }

  return(formula)
}


# Return fitted formula column display names keyed by internal matrix names.
.fitted_formula_column_display_map <- function(design) {

  if (is.null(design[["column_names"]]) || is.null(design[["raw_column_names"]])) {
    return(NULL)
  }

  if (length(design[["column_names"]]) != length(design[["raw_column_names"]])) {
    return(NULL)
  }

  return(stats::setNames(
    .formula_design_display_names(design[["raw_column_names"]]),
    design[["column_names"]]
  ))
}


# Convert BayesTools' JAGS-safe interaction separator to formula syntax.
.formula_design_display_names <- function(x) {

  return(gsub("__xXx__", ":", x, fixed = TRUE))
}


# ---------------------------------------------------------------------------- #
# .get_model_matrix
# ---------------------------------------------------------------------------- #
#
# Extract the fitted location model matrix X from a brma object.
#
# ---------------------------------------------------------------------------- #
.get_model_matrix <- function(object) {

  is_mods  <- .is_mods(object)
  is_PET   <- .is_PET(object)
  is_PEESE <- .is_PEESE(object)

  if (!is_mods) {
    K <- nrow(object[["data"]][["outcome"]])
    X <- matrix(1, nrow = K, ncol = 1)
    colnames(X) <- "(Intercept)"
    attr(X, "assign") <- 0L
  } else {
    design <- .fitted_formula_design(object, "mu", required = TRUE)
    X      <- design[["model_matrix"]]

    attr(X, "assign")       <- design[["assign"]]
    attr(X, "rank")         <- design[["rank"]]
    attr(X, "aliased")      <- design[["aliased"]]
    attr(X, "term_labels")  <- .fitted_formula_terms(
      object            = object,
      parameter         = "mu",
      include_intercept = FALSE,
      display           = TRUE,
      required          = TRUE
    )

    raw_colnames <- .fitted_formula_column_display_map(design)
    if (!is.null(raw_colnames)) {
      attr(X, "raw_colnames") <- raw_colnames
    }
  }

  if (is_PET || is_PEESE) {
    outcome_data <- object[["data"]][["outcome"]]
    effect_direction <- .effect_direction(object)
    assign       <- attr(X, "assign")
    raw_colnames <- attr(X, "raw_colnames")
    aliased      <- attr(X, "aliased")
    term_labels  <- attr(X, "term_labels")
    if (is.null(assign)) {
      assign <- seq_len(ncol(X)) - 1L
    }
    next_assign  <- max(assign) + 1L

    bias_parameters <- c(if (is_PET) "PET", if (is_PEESE) "PEESE")
    for (parameter in bias_parameters) {
      predictor <- .bias_regression_predictor(
        sei              = outcome_data[["sei"]],
        parameter        = parameter,
        effect_direction = effect_direction
      )
      X <- cbind(X, predictor)
      colnames(X)[ncol(X)] <- parameter
      assign       <- c(assign, next_assign)
      raw_colnames <- c(
        raw_colnames,
        stats::setNames(parameter, parameter)
      )
      aliased      <- c(aliased, stats::setNames(FALSE, parameter))
      term_labels  <- c(term_labels, parameter)
      next_assign  <- next_assign + 1L
    }

    attr(X, "assign") <- assign
    if (!is.null(raw_colnames) && length(raw_colnames) == ncol(X)) {
      attr(X, "raw_colnames") <- raw_colnames
    }
    if (!is.null(aliased) && length(aliased) == ncol(X)) {
      attr(X, "aliased") <- aliased
    }
    if (!is.null(term_labels)) {
      attr(X, "term_labels") <- term_labels
    }
  }

  attr(X, "rank") <- qr(X)[["rank"]]

  return(X)
}
