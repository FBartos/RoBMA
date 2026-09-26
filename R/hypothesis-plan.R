# ============================================================================ #
# hypothesis-plan.R
# ============================================================================ #
#
# One plan per hypothesis statement. A plan records what the statement tests
# and how each density method can evaluate it:
#   - kind: "scalar" (a scalar parameter or coefficient), "linear" (factor
#     levels of every contrast, ordered levels, and linear combinations of
#     the levels of one term: level contrasts and derived linear
#     expressions), "random" (random-effect quantities), "exp_affine"
#     (exponentiated affine original-scale coefficients), or "marginal_mean";
#   - targets: one per point value, with the target's weights on the fitted
#     coordinates, its prior density built by BayesTools, its declared
#     atoms, and the classification of the value by
#     BayesTools::prior_ordinate_status();
#   - methods: the status of KDE, qCMDE, IWMDE, and the normal
#     approximation, NULL when the method is available and a refusal (the
#     reason and the condition classes a refused call stops with) otherwise.
# hypothesis(), hypothesis.marginal_means.brma(), and hypothesis_quantities()
# execute or render plans only: the eligibility of a statement has one
# source. Region statements are method-free: their plans refuse only
# statements on fixed quantities or unsupported targets.
#
# ============================================================================ #


.hypothesis_plan_methods <- function() {

  c("KDE", "qCMDE", "IWMDE", "normal")
}


# The methods hypothesis_quantities() lists; the normal approximation is a
# rough check of near-normal interior ordinates and is not advertised.
.hypothesis_plan_advertised_methods <- function() {

  c("KDE", "qCMDE", "IWMDE")
}


# A refusal: the reason and the condition classes of the error a refused
# statement stops with. BayesTools refusals keep their classes ('condition',
# e.g. "BayesTools_point_mass_at_null"). RoBMA's own refusals of a test that
# cannot be computed have the class "RoBMA_hypothesis_<type>" with the
# parent "RoBMA_hypothesis_unavailable": "fixed" (the quantity is fixed by
# the fitted model), "target" (no supported target or prior density), and
# "method" (the density method is unavailable for the target). A statement
# that has to be restated, or an argument that has to change, is no
# unavailable test: "statement" (the statement's form is unsupported) has
# the class "RoBMA_hypothesis_statement" only, and "ambiguous" (the
# statement names several parameters; 'component' or, for marginal means,
# 'parameter' selects one) has it as its parent
# (.hypothesis_ambiguous_class()).
# Method refusals by the qCMDE/IWMDE capability of the fitted model
# (.iwmde_capability()), and the qCMDE/IWMDE refusals of the causes for
# which plot() refuses these methods (conditional random-effect statements,
# random-effect quantities without a scalar coordinate, nonlinear
# original-scale coefficients), also have the classes of that refusal:
# "RoBMA_density_method_<cause>" and the parent
# "RoBMA_density_method_unavailable" (.hypothesis_refusal_density_method()).
.hypothesis_refusal <- function(reason, type = NULL, condition = NULL) {

  class <- if (!is.null(condition)) {
    unique(condition)
  } else if (identical(type, "ambiguous")) {
    .hypothesis_ambiguous_class()
  } else if (identical(type, "statement")) {
    "RoBMA_hypothesis_statement"
  } else {
    c(paste0("RoBMA_hypothesis_", type), "RoBMA_hypothesis_unavailable")
  }

  list(reason = reason, class = class)
}


# The classes of an ambiguous statement, at every entry point of
# hypothesis(): its own class and the parent "RoBMA_hypothesis_statement".
.hypothesis_ambiguous_class <- function() {

  c("RoBMA_hypothesis_ambiguous", "RoBMA_hypothesis_statement")
}


# A hypothesis() refusal whose cause plot() (and, for the capability of the
# fitted model, marginal_means()) refuses as a qCMDE/IWMDE request
# ('unavailable', a refusal of .iwmde_unavailable() or .iwmde_capability()):
# the refusal of 'type' with the reason of that refusal, followed by its
# classes ("RoBMA_density_method_<cause>" and
# "RoBMA_density_method_unavailable"), so that one cause has the same
# classes at every entry point.
.hypothesis_refusal_density_method <- function(unavailable, type = "method") {

  refusal <- .hypothesis_refusal(unavailable[["reason"]], type)
  refusal[["class"]] <- c(refusal[["class"]], unavailable[["class"]])

  refusal
}


# An unknown reference is a statement error with the classes with which
# BayesTools refuses it ("BayesTools_parameter_not_found" and its parent
# "BayesTools_parameter_resolution_error"), after "RoBMA_hypothesis_statement":
# the classes of an unknown name or level at every entry point of
# hypothesis().
.hypothesis_not_found_class <- function() {

  c("RoBMA_hypothesis_statement", "BayesTools_parameter_not_found",
    "BayesTools_parameter_resolution_error")
}


# The refusal of a statement that references an unknown level of a factor
# term (the resolver refuses unknown levels of fitted terms first; marginal
# means are refused by their plans).
.hypothesis_unknown_level_refusal <- function(level, parameter) {

  .hypothesis_refusal(
    paste0(
      "Hypothesis references unknown level '", level, "' for parameter '",
      parameter, "'."
    ),
    condition = .hypothesis_not_found_class()
  )
}


.hypothesis_refusal_from_condition <- function(condition, reason = NULL) {

  list(
    reason = if (is.null(reason)) conditionMessage(condition) else reason,
    class  = setdiff(class(condition), c("error", "condition"))
  )
}


.hypothesis_stop <- function(refusal) {

  stop(structure(
    class = c(refusal[["class"]], "error", "condition"),
    list(message = refusal[["reason"]], call = NULL)
  ))
}


# Fitted metadata that a hypothesis needs are missing or unsupported (a fit
# of an older RoBMA/BayesTools build): the error of class
# "RoBMA_refit_required" (.refit_required_class()), whose message asks for
# a refit. BayesTools refuses such fits with unclassed errors, so the class
# has no BayesTools parent.
.stop_refit_required <- function(...) {

  stop(structure(
    class = c(.refit_required_class(), "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}


.refit_required_class <- function() {

  "RoBMA_refit_required"
}


# The status of a density method in a plan: NULL when the method evaluates
# the statement, otherwise its refusal.
.hypothesis_plan_status <- function(plan, method) {

  plan[["methods"]][[method]]
}


.hypothesis_plan_check <- function(plan, method) {

  refusal <- .hypothesis_plan_status(plan, method)
  if (!is.null(refusal)) {
    .hypothesis_stop(refusal)
  }

  invisible(plan)
}


# A cache shared by the plans of one call: mixed posteriors, marginal
# posteriors, and random-effect selections are built once per target.
.hypothesis_plan_cache <- function() {

  new.env(parent = emptyenv())
}


.hypothesis_plan_cached <- function(cache, key, compute) {

  if (is.null(cache)) {
    return(compute())
  }
  if (!exists(key, envir = cache, inherits = FALSE)) {
    assign(key, compute(), envir = cache)
  }

  get(key, envir = cache, inherits = FALSE)
}


# ---------------------------------------------------------------------------- #
# Plans of fitted brma objects
# ---------------------------------------------------------------------------- #

# The plans of the statements of a hypothesis, one per statement.
.hypothesis_plans <- function(object, hypothesis, component = "auto",
                              standardized = FALSE, conditional = FALSE,
                              metadata = NULL, n_samples = 10000L,
                              cache = .hypothesis_plan_cache()) {

  if (is.null(metadata)) {
    metadata <- .brma_parameter_catalog_metadata(object)
  }
  ast        <- .hypothesis_brma_ast(hypothesis, metadata[["catalog"]])
  statements <- BayesTools::hypothesis_render(ast)

  lapply(statements, function(statement) {
    .hypothesis_plan(
      object       = object,
      statement    = statement,
      component    = component,
      standardized = standardized,
      conditional  = conditional,
      metadata     = metadata,
      n_samples    = n_samples,
      cache        = cache
    )
  })
}


# The plan of one statement of a fitted brma object on the fitted
# ('standardized') or original coefficient scale, conditional on the
# inclusion of the tested parameter ('conditional', model-averaged objects).
.hypothesis_plan <- function(object, statement, component = "auto",
                             standardized = FALSE, conditional = FALSE,
                             metadata = NULL, n_samples = 10000L,
                             cache = .hypothesis_plan_cache()) {

  if (is.null(metadata)) {
    metadata <- .brma_parameter_catalog_metadata(object)
  }
  statement <- .hypothesis_brma_ast(statement, metadata[["catalog"]])
  if (length(statement[["statements"]]) != 1L) {
    stop("Internal error: a hypothesis plan covers one statement.",
         call. = FALSE)
  }
  selected <- .hypothesis_brma_select_parameter(
    object     = object,
    hypothesis = statement,
    component  = component,
    metadata   = metadata
  )
  .hypothesis_brma_check_supported_component(selected[["component"]])
  parameter <- selected[["parameter"]]

  plan <- .hypothesis_plan_new(
    statement    = statement,
    hypothesis   = .hypothesis_brma_rewrite(
      hypothesis = statement,
      aliases    = selected[["aliases"]],
      parameter  = parameter
    ),
    parameter    = parameter,
    label        = .hypothesis_brma_alias_label(selected[["aliases"]], parameter),
    component    = selected[["component"]],
    selected     = selected,
    standardized = standardized,
    conditional  = conditional,
    n_samples    = n_samples
  )
  entry <- selected[["entry"]]
  plan  <- if (identical(selected[["component"]], "random")) {
    .hypothesis_plan_random(plan, object, cache)
  } else if (identical(entry[["role"]], "formula_coefficient_group")) {
    .hypothesis_plan_linear(plan, object, metadata, cache)
  } else {
    .hypothesis_plan_scalar(plan, object, cache)
  }

  .hypothesis_plan_finish(plan, object)
}


.hypothesis_plan_new <- function(statement, hypothesis, parameter, label,
                                 component, selected = NULL,
                                 standardized = FALSE, conditional = FALSE,
                                 n_samples = 10000L) {

  sides <- unlist(lapply(hypothesis[["statements"]], function(x) {
    c(x[["left"]][["type"]], x[["right"]][["type"]])
  }), use.names = FALSE)
  refs  <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )

  list(
    kind            = NA_character_,
    route           = NA_character_,
    statement       = statement,
    hypothesis      = hypothesis,
    parameter       = parameter,
    label           = label,
    component       = component,
    selected        = selected,
    standardized    = standardized,
    conditional     = conditional,
    n_samples       = n_samples,
    point           = any(sides %in% c("point", "not_point")),
    region          = any(sides == "region"),
    refs            = refs,
    targets         = list(),
    refusal         = NULL,
    method_refusals = list(),
    support         = NULL,
    group           = NA_character_
  )
}


# The label of a point side of the statement, as the user wrote it.
.hypothesis_plan_side_label <- function(plan, side) {

  statement <- plan[["statement"]][["statements"]][[1L]]
  label     <- statement[[side]][["label"]]
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    label <- statement[["source"]]
  }

  label
}


# A point target: the value, its label, the prior density, the declared
# atoms, the classification of the value, and the qCMDE/IWMDE parameter spec
# (NULL when the method refusals of the plan cover the target).
.hypothesis_plan_target <- function(value, label, prior_density,
                                    spec = NULL, atoms = NULL, level = NA,
                                    refusal = NULL) {

  status <- if (is.null(refusal) && !is.null(prior_density)) {
    BayesTools::prior_ordinate_status(prior_density, value, labels = label)
  }

  list(
    value         = value,
    label         = label,
    level         = level,
    prior_density = prior_density,
    atoms         = atoms,
    status        = status,
    spec          = spec,
    refusal       = refusal
  )
}


# The refusal of a point target: its own refusal, or its prior ordinate
# classification when the value is not eligible (BayesTools classes).
.hypothesis_plan_target_refusal <- function(target) {

  if (!is.null(target[["refusal"]])) {
    return(target[["refusal"]])
  }
  if (is.null(target[["prior_density"]])) {
    return(.hypothesis_refusal(
      paste0(
        "The prior density of '", target[["label"]], "' is unavailable, so ",
        "the point hypothesis has no Savage-Dickey Bayes factor. Use a region ",
        "or directional hypothesis."
      ),
      "target"
    ))
  }
  status <- target[["status"]]
  if (is.null(status) || all(status[["eligible"]])) {
    return(NULL)
  }
  i         <- which(!status[["eligible"]])[[1L]]
  condition <- status[["condition"]][[i]]
  reason    <- status[["reason"]][[i]]
  if (identical(condition, "BayesTools_point_mass_at_null") &&
      !is.null(target[["point_mass_reason"]])) {
    reason <- target[["point_mass_reason"]]
  }

  refusal <- .hypothesis_refusal(
    reason    = reason,
    condition = c(condition, "BayesTools_hypothesis_ordinate")
  )
  refusal[["value_specific"]] <- TRUE

  refusal
}


.hypothesis_plan_finish <- function(plan, object) {

  methods <- .hypothesis_plan_methods()
  plan[["methods"]] <- stats::setNames(lapply(methods, function(method) {
    .hypothesis_plan_method_refusal(plan, object, method)
  }), methods)

  plan
}


# The refusal of one density method: a refusal of the whole plan (fixed
# quantities, unsupported targets and statements), then the refusals of the
# point targets (their values), then the refusal of an object without the
# fitted model the method needs (a marginal-means object without its source
# fit), then the qCMDE/IWMDE capability of the fitted model, then the
# quantity's own refusal of the method. The capability precedes the
# quantity's refusal as in plot(), which checks it before selecting the
# quantity, so that a qCMDE/IWMDE request with both causes reports the
# capability at both entry points. Region statements need no density method.
.hypothesis_plan_method_refusal <- function(plan, object, method) {

  if (!is.null(plan[["refusal"]])) {
    return(plan[["refusal"]])
  }
  if (!plan[["point"]]) {
    return(NULL)
  }
  for (target in plan[["targets"]]) {
    refusal <- .hypothesis_plan_target_refusal(target)
    if (!is.null(refusal)) {
      return(refusal)
    }
  }
  if (!is.null(plan[["source_refusals"]][[method]])) {
    return(plan[["source_refusals"]][[method]])
  }
  if (method %in% c("qCMDE", "IWMDE")) {
    capability_object <- if (is.null(plan[["capability_object"]])) {
      object
    } else {
      plan[["capability_object"]]
    }
    capability <- .iwmde_capability(
      object         = capability_object,
      density_method = method
    )
    if (!capability[["available"]]) {
      # A capability refusal also has the classes of the capability refusal
      # outside hypothesis(): its cause and the parent class.
      return(.hypothesis_refusal_density_method(capability, "method"))
    }
  }

  plan[["method_refusals"]][[method]]
}


# A quantity fixed by the fitted model (catalog status or fixed value).
.hypothesis_plan_fixed_refusal <- function(label) {

  .hypothesis_refusal(
    paste0(
      "The quantity '", label, "' is fixed by the fitted model; posterior ",
      "hypothesis tests are undefined."
    ),
    "fixed"
  )
}


.hypothesis_plan_quantity_fixed <- function(status, fixed_value) {

  identical(as.character(status), "structural") ||
    (length(fixed_value) == 1L && is.finite(fixed_value))
}


.hypothesis_plan_direct_refusal <- function(plan) {

  refs <- plan[["refs"]]
  if (nrow(refs) == 0L || all(refs[["direct"]])) {
    return(NULL)
  }

  .hypothesis_refusal(
    paste0(
      "Point-null hypotheses require a direct parameter or level reference; ",
      "unsupported point expression in: '",
      refs[["hypothesis"]][!refs[["direct"]]][[1L]], "'."
    ),
    "statement"
  )
}


# The draws of a fixed-effect target: the mixed posteriors of the tested
# parameter and of the parameters owning the coordinates that its prior
# density combines, on the fitted ('standardized') or original scale, and the
# marginal posterior of the tested parameter. 'prior_density' is set as the
# parameter's prior density before the marginal posterior is built.
.hypothesis_plan_draws <- function(plan, object, cache,
                                   extra_parameters = character(),
                                   prior_density = NULL, key = "") {

  parameter        <- plan[["parameter"]]
  sample_parameter <- unique(c(
    .as_mixed_posteriors_parameters(object, parameter),
    extra_parameters
  ))
  cache_key <- paste(
    "draws", parameter, paste(sort(sample_parameter), collapse = ","),
    plan[["standardized"]], plan[["conditional"]], plan[["n_samples"]], key,
    sep = "\r"
  )

  .hypothesis_plan_cached(cache, cache_key, function() {
    samples <- .brma_as_mixed_posteriors(
      object           = object,
      parameters       = sample_parameter,
      conditional      = if (plan[["conditional"]]) parameter else NULL,
      conditional_rule = "AND",
      transform_scaled = !plan[["standardized"]],
      n_prior_samples  = plan[["n_samples"]]
    )
    density <- if (is.function(prior_density)) {
      prior_density(samples)
    } else {
      prior_density
    }
    if (!is.null(density)) {
      prior_densities <- BayesTools::posterior_metadata(samples, "prior_densities")
      prior_densities[[parameter]] <- density
      BayesTools::posterior_metadata(samples, "prior_densities") <- prior_densities
    }
    density_parameter <- .plot_brma_density_sample_parameter(
      samples          = samples,
      parameter        = parameter,
      sample_parameter = sample_parameter
    )
    posterior <- BayesTools::marginal_posterior(
      samples       = samples,
      parameter     = density_parameter,
      prior_samples = TRUE,
      use_formula   = FALSE,
      n_samples     = plan[["n_samples"]]
    )
    if (!is.null(density)) {
      BayesTools::posterior_metadata(posterior, "prior_density") <- density
    }

    list(
      samples           = samples,
      posterior         = posterior,
      sample_parameter  = sample_parameter,
      density_parameter = density_parameter,
      prior_density     = density
    )
  })
}


# The route of an original-scale formula coefficient as BayesTools declares
# it in the fitted coefficient transform: the map type (identity, affine,
# exp_affine, or unsupported), the support of the map, and the target's
# weights on the fitted coordinates. An unsupported route has a reason and
# its cause: "metadata" (the transform lacks the certified metadata),
# "fixed" (the coefficient is structurally fixed), or "nonlinear" (a
# nonlinear joint transform). NULL for parameters that are no formula
# coefficients.
.brma_formula_coefficient_route <- function(object, selected) {

  target_info <- .hypothesis_brma_formula_coefficient_target(
    object   = object,
    selected = selected
  )
  if (is.null(target_info)) {
    return(NULL)
  }
  transform <- target_info[["transform"]]
  target    <- target_info[["target"]]
  targets   <- transform[["targets"]]
  row       <- if (is.data.frame(targets) &&
                   all(c("target", "map_type", "support") %in% names(targets))) {
    match(target, targets[["target"]])
  } else {
    NA_integer_
  }
  if (!identical(transform[["target_scale"]], "original") || is.na(row)) {
    return(c(target_info, list(
      type   = "unsupported",
      reason = paste0(
        "The fitted coefficient transform for '", target,
        "' lacks the certified structural metadata required for hypothesis testing."
      ),
      cause  = "metadata"
    )))
  }
  weights <- stats::setNames(
    as.numeric(transform[["matrix"]][target_info[["target_i"]], , drop = FALSE]),
    colnames(transform[["matrix"]])
  )
  weights  <- weights[weights != 0]
  map_type <- targets[["map_type"]][[row]]
  route <- c(target_info, list(
    type    = map_type,
    weights = weights,
    support = targets[["support"]][[row]]
  ))
  if (length(weights) == 0L) {
    route[["type"]]   <- "unsupported"
    route[["reason"]] <- paste0(
      "The fitted coefficient '", target,
      "' is structurally fixed and has no posterior hypothesis route."
    )
    route[["cause"]]  <- "fixed"
  } else if (!map_type %in% c("identity", "affine", "exp_affine")) {
    route[["type"]]   <- "unsupported"
    route[["reason"]] <- paste0(
      "The fitted nonlinear joint coefficient transform for '", target,
      "' is not supported by hypothesis()."
    )
    route[["cause"]]  <- "nonlinear"
  }

  route
}


# Point values outside the open support of an original-scale formula target.
.hypothesis_plan_support_refusal <- function(value, support, description) {

  if (is.null(support)) {
    return(NULL)
  }
  outside <- !is.finite(value) ||
    (is.finite(support[[1L]]) && value <= support[[1L]]) ||
    (is.finite(support[[2L]]) && value >= support[[2L]])
  if (!outside) {
    return(NULL)
  }

  refusal <- .hypothesis_refusal(
    paste0(
      "Point-null value ", value, " is outside or on the boundary of the ",
      "open support for ", description, "."
    ),
    "target"
  )
  refusal[["value_specific"]] <- TRUE

  refusal
}


# The reason of a prior point mass at a tested value of a model-averaged
# parameter names the inclusion Bayes factor.
.hypothesis_plan_point_mass_reason <- function(object) {

  if (!.is_RoBMA(object)) {
    return(NULL)
  }

  paste0(
    "There is a point mass in the prior at the exact null hypothesis value. ",
    "The Savage-Dickey density ratio is invalid. This parameter has a null ",
    "component, so its evidence against the null is the inclusion Bayes ",
    "factor reported by 'summary()' and 'summary_models()'."
  )
}


# Whether BayesTools evaluates the region probabilities of a prior density
# (BayesTools::prior_density_has_provenance()); a density without that
# provenance, e.g. the plotting grid of a nested variance allocation, has none.
.hypothesis_plan_density_has_provenance <- function(prior_density) {

  !is.null(prior_density) && BayesTools::prior_density_has_provenance(prior_density)
}


# The exact support bounds BayesTools declares for posterior draws, or NULL.
.hypothesis_plan_support_bounds <- function(x) {

  support <- BayesTools::posterior_metadata(x, "support")
  if (is.list(support) && !inherits(support, "BayesTools_posterior_support") &&
      length(support) == 1L) {
    support <- support[[1L]]
  }
  if (!inherits(support, "BayesTools_posterior_support") ||
      !isTRUE(support[["exact"]])) {
    return(NULL)
  }

  as.numeric(support[["bounds"]])
}


# ---------------------------------------------------------------------------- #
# Scalar and exp(affine) targets
# ---------------------------------------------------------------------------- #

.hypothesis_plan_scalar <- function(plan, object, cache) {

  entry <- plan[["selected"]][["entry"]]
  plan[["kind"]]  <- "scalar"
  plan[["route"]] <- "scalar"
  plan[["group"]] <- paste("scalar", plan[["parameter"]], sep = "\r")
  if (.hypothesis_plan_quantity_fixed(entry[["status"]], entry[["fixed_value"]])) {
    plan[["refusal"]] <- .hypothesis_plan_fixed_refusal(plan[["label"]])
    return(plan)
  }
  plan[["refusal"]] <- .hypothesis_plan_direct_refusal(plan)
  if (!is.null(plan[["refusal"]])) {
    return(plan)
  }

  route <- if (!plan[["standardized"]]) {
    .brma_formula_coefficient_route(object, plan[["selected"]])
  }
  if (identical(route[["type"]], "unsupported")) {
    plan[["refusal"]] <- .hypothesis_refusal(route[["reason"]], "target")
    return(plan)
  }
  if (identical(route[["type"]], "exp_affine")) {
    return(.hypothesis_plan_exp_affine(plan, object, route, cache))
  }

  formula_density <- if (!is.null(route)) {
    function(samples) {
      BayesTools::JAGS_formula_prior_density(
        fit          = object[["fit"]],
        parameter    = route[["formula_parameter"]],
        target       = route[["target"]],
        target_scale = "original",
        context      = BayesTools::posterior_metadata(samples, "prior_context")
      )
    }
  }
  draws <- .hypothesis_plan_draws(
    plan             = plan,
    object           = object,
    cache            = cache,
    extra_parameters = if (!is.null(route)) {
      .hypothesis_brma_target_prior_parameters(object, route[["weights"]])
    } else {
      character()
    },
    prior_density    = formula_density,
    key              = "scalar"
  )
  prior_density <- BayesTools::posterior_metadata(draws[["posterior"]], "prior_density")
  spec <- if (is.null(route) || identical(route[["type"]], "identity")) {
    list(type = "primitive", prior_density = prior_density)
  } else {
    list(type = "linear", weights = route[["weights"]], prior_density = prior_density)
  }
  plan[["draws"]]         <- draws
  plan[["route_info"]]    <- route
  plan[["prior_density"]] <- prior_density
  plan[["support"]]       <- if (is.null(route[["support"]])) {
    .hypothesis_plan_support_bounds(draws[["posterior"]])
  } else {
    route[["support"]]
  }
  plan[["weights"]]       <- if (is.null(route)) {
    stats::setNames(1, plan[["parameter"]])
  } else {
    route[["weights"]]
  }
  plan[["targets"]] <- .hypothesis_plan_point_targets(
    plan          = plan,
    object        = object,
    prior_density = prior_density,
    spec          = spec,
    atoms         = BayesTools::posterior_metadata(draws[["posterior"]], "atoms"),
    support       = route[["support"]],
    description   = paste0("transformed coefficient '", plan[["parameter"]], "'")
  )
  if (length(plan[["targets"]]) > 0L &&
      is.list(draws[["posterior"]])) {
    plan[["method_refusals"]] <- .hypothesis_plan_precomputed_refusals(paste0(
      "qCMDE/IWMDE point hypotheses for factor parameters must specify ",
      "a level, e.g. '", plan[["label"]], "[level] = ",
      plan[["targets"]][[1L]][["value"]], "'."
    ))
  }

  plan
}


# One point target per point side of the statement of a scalar quantity.
.hypothesis_plan_point_targets <- function(plan, object, prior_density,
                                           spec = NULL, atoms = NULL,
                                           support = NULL, description = NULL) {

  refs <- plan[["refs"]]
  lapply(seq_len(nrow(refs)), function(i) {
    label  <- .hypothesis_plan_side_label(plan, refs[["side"]][[i]])
    target <- .hypothesis_plan_target(
      value         = refs[["value"]][[i]],
      label         = label,
      prior_density = prior_density,
      spec          = spec,
      atoms         = atoms,
      refusal       = if (!is.null(support)) {
        .hypothesis_plan_support_refusal(
          refs[["value"]][[i]], support, description
        )
      }
    )
    target[["point_mass_reason"]] <- .hypothesis_plan_point_mass_reason(object)
    target
  })
}


.hypothesis_plan_precomputed_refusals <- function(reason, type = "method") {

  refusal <- .hypothesis_refusal(reason, type)
  list(qCMDE = refusal, IWMDE = refusal)
}


# An exponentiated affine original-scale coefficient: KDE point and region
# hypotheses on the coefficient's own scale for a certified atom-free,
# unconditional target, with the prior density of the fitted coefficient
# transform.
.hypothesis_plan_exp_affine <- function(plan, object, route, cache) {

  plan[["kind"]]  <- "exp_affine"
  plan[["route"]] <- "scalar"
  target <- route[["target"]]
  if (plan[["conditional"]]) {
    plan[["refusal"]] <- .hypothesis_refusal(
      paste0(
        "Nonlinear fitted-scale KDE hypotheses are unavailable for conditional ",
        "product-space posteriors."
      ),
      "target"
    )
    return(plan)
  }
  draws <- .hypothesis_plan_exp_affine_draws(plan, object, route, cache)
  plan[["refusal"]] <- .hypothesis_plan_exp_affine_certify(
    sample        = draws[["samples"]][[target]],
    prior_density = draws[["prior_density"]]
  )
  if (!is.null(plan[["refusal"]])) {
    return(plan)
  }

  plan[["draws"]]         <- draws
  plan[["route_info"]]    <- route
  plan[["prior_density"]] <- draws[["prior_density"]]
  plan[["support"]]       <- route[["support"]]
  plan[["weights"]]       <- route[["weights"]]
  plan[["targets"]]       <- .hypothesis_plan_point_targets(
    plan          = plan,
    object        = object,
    prior_density = draws[["prior_density"]],
    atoms         = BayesTools::posterior_metadata(draws[["posterior"]], "atoms"),
    support       = route[["support"]],
    description   = paste0("transformed coefficient '", target, "'")
  )
  reason <- paste0(
    "The requested nonlinear fitted-scale hypothesis for '",
    plan[["parameter"]], "' is supported only with density_method = 'KDE'. ",
    "qCMDE/IWMDE ordinates support only direct parameter or level point ",
    "hypotheses with an exact linear fitted-scale map."
  )
  # qCMDE/IWMDE are refused for the cause for which plot() refuses them
  # ("original_scale"); the normal approximation is no qCMDE/IWMDE request.
  precomputed <- .hypothesis_refusal_density_method(
    .iwmde_unavailable(reason, "original_scale"),
    "method"
  )
  plan[["method_refusals"]] <- list(
    qCMDE  = precomputed,
    IWMDE  = precomputed,
    normal = .hypothesis_refusal(reason, "method")
  )

  plan
}


# The mixed posteriors of an exp(affine) target and of the parameters owning
# the fitted coordinates its prior combines, and the target's marginal
# posterior on its original scale (BayesTools::marginal_posterior()): its
# draws with the prior density of the fitted coefficient transform, the exact
# support of the map, and the declared atoms and conditioning of the mixed
# posterior.
.hypothesis_plan_exp_affine_draws <- function(plan, object, route, cache) {

  parameter        <- plan[["parameter"]]
  sample_parameter <- unique(c(
    .as_mixed_posteriors_parameters(object, parameter),
    .hypothesis_brma_target_prior_parameters(object, route[["weights"]])
  ))
  .hypothesis_plan_cached(
    cache,
    paste("exp_affine", parameter, plan[["n_samples"]], sep = "\r"),
    function() {
      samples <- .brma_as_mixed_posteriors(
        object           = object,
        parameters       = sample_parameter,
        transform_scaled = TRUE,
        n_prior_samples  = plan[["n_samples"]]
      )
      posterior <- BayesTools::marginal_posterior(
        samples       = samples,
        parameter     = route[["target"]],
        prior_samples = TRUE,
        use_formula   = FALSE,
        n_samples     = plan[["n_samples"]]
      )

      list(
        samples       = samples,
        posterior     = posterior,
        prior_density = BayesTools::posterior_metadata(posterior, "prior_density")
      )
    }
  )
}


# Structural certification of an exp(affine) target: the averaged,
# atom-free posterior of a scalar target whose prior density is atom-free.
.hypothesis_plan_exp_affine_certify <- function(sample, prior_density) {

  refusal <- function(reason) .hypothesis_refusal(reason, "target")
  if (is.null(sample) || !inherits(sample, "mixed_posteriors.simple")) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require a certified scalar ",
      "mixed posterior."
    )))
  }
  condition <- BayesTools::posterior_metadata(sample, "condition")
  if (!isTRUE(condition[["averaged"]])) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence for ",
      "an unconditional posterior."
    )))
  }
  if (!BayesTools::posterior_atoms_free(sample)) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence that ",
      "the posterior is atom-free."
    )))
  }
  points <- prior_density[["points"]]
  if (!inherits(prior_density, "prior_density") || !is.data.frame(points) ||
      !all(c("x", "p") %in% names(points)) || nrow(points) > 0L) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence that ",
      "the prior is atom-free."
    )))
  }

  NULL
}


# ---------------------------------------------------------------------------- #
# Linear targets: factor levels and linear combinations of levels
# ---------------------------------------------------------------------------- #

# Levels of every contrast are linear in the fitted coordinates. A statement
# whose point sides reference levels directly ('g[a] = 0', 'g[a] = 0 vs
# g[b] = 0.1') and a region statement ('g[a] > g[b]') are evaluated on the
# marginal posterior of the factor term ("levels" route); a statement whose
# point sides are linear combinations of levels (level contrasts such as
# 'g[a] = g[b]' and derived linear expressions such as '2 * g[a] = 0.1')
# is evaluated on the scalar linear target of BayesTools
# hypothesis_linear_target() ("combination" route). On the original scale of
# formula coefficients, a target's prior density is
# BayesTools::JAGS_formula_prior_density() at the target's weights on the
# original coefficients; on the fitted scale it is the density BayesTools
# attaches to the level or linear target.
.hypothesis_plan_linear <- function(plan, object, metadata, cache) {

  plan[["kind"]] <- "linear"
  entry     <- plan[["selected"]][["entry"]]
  parameter <- plan[["parameter"]]
  levels    <- .hypothesis_plan_term_levels(metadata, entry)
  occurrences <- BayesTools::hypothesis_symbols(
    plan[["hypothesis"]],
    occurrences = TRUE
  )
  referenced <- unique(occurrences[["level"]][occurrences[["parameter"]] == parameter])
  whole_term <- any(is.na(referenced))
  referenced <- referenced[!is.na(referenced)]
  unknown <- setdiff(referenced, levels[["level"]])
  if (length(unknown) > 0L) {
    plan[["refusal"]] <- .hypothesis_unknown_level_refusal(
      level     = unknown[[1L]],
      parameter = plan[["label"]]
    )
    return(plan)
  }
  # A statement tests a level the contrast fixes (the treatment reference
  # level) when it states a point on it directly, references the whole term
  # (a point or region event on every level jointly), or references no other
  # level; linear combinations and comparisons with other levels remain
  # defined.
  fixed_levels <- levels[["level"]][levels[["fixed"]]]
  refs         <- plan[["refs"]]
  tested       <- unique(c(
    refs[["level"]][refs[["direct"]] & !is.na(refs[["level"]])],
    if (whole_term) levels[["level"]],
    if (length(referenced) > 0L && all(referenced %in% fixed_levels)) referenced
  ))
  fixed <- intersect(tested, fixed_levels)
  if (length(fixed) > 0L) {
    plan[["refusal"]] <- .hypothesis_plan_fixed_refusal(
      paste0(plan[["label"]], "[", fixed[[1L]], "]")
    )
    return(plan)
  }

  scale <- .hypothesis_plan_linear_scale(
    plan        = plan,
    object      = object,
    entry       = entry,
    coordinates = unique(unlist(levels[["dependencies"]], use.names = FALSE))
  )
  draws <- .hypothesis_plan_draws(
    plan             = plan,
    object           = object,
    cache            = cache,
    extra_parameters = scale[["prior_parameters"]],
    key              = "linear"
  )
  posterior <- draws[["posterior"]]
  if (!is.list(posterior)) {
    stop("Internal error: factor levels have no factor marginal posterior.",
         call. = FALSE)
  }
  plan[["draws"]] <- draws
  plan[["scale"]] <- scale

  combination <- nrow(plan[["refs"]]) > 0L && any(!plan[["refs"]][["direct"]])
  if (combination) {
    return(.hypothesis_plan_combination(plan, object))
  }

  plan[["route"]] <- "levels"
  plan[["group"]] <- paste("levels", parameter, sep = "\r")
  refs <- plan[["refs"]]
  if (length(referenced) > 0L) {
    plan[["support"]] <- .hypothesis_plan_support_bounds(
      posterior[[referenced[[1L]]]]
    )
  }
  plan[["targets"]] <- unlist(lapply(seq_len(nrow(refs)), function(i) {
    ref_levels <- if (is.na(refs[["level"]][[i]])) {
      levels[["level"]]
    } else {
      refs[["level"]][[i]]
    }
    label <- .hypothesis_plan_side_label(plan, refs[["side"]][[i]])
    lapply(ref_levels, function(level) {
      .hypothesis_plan_level_target(
        plan   = plan,
        object = object,
        level  = level,
        value  = refs[["value"]][[i]],
        label  = if (is.na(refs[["level"]][[i]])) {
          paste0(plan[["label"]], "[", level, "] = ", refs[["value"]][[i]])
        } else {
          label
        }
      )
    })
  }), recursive = FALSE)
  if (whole_term && plan[["point"]]) {
    plan[["method_refusals"]] <- .hypothesis_plan_precomputed_refusals(paste0(
      "qCMDE/IWMDE point hypotheses for factor parameters must specify ",
      "a level, e.g. '", plan[["label"]], "[level] = ",
      refs[["value"]][[1L]], "'."
    ))
  }

  plan
}


# The levels of a factor term: their labels and whether the contrast fixes
# them (the treatment reference level), from the catalog.
.hypothesis_plan_term_levels <- function(metadata, entry) {

  quantities <- metadata[["catalog"]][["quantities"]]
  members    <- quantities[
    quantities[["quantity_id"]] %in%
      unlist(entry[["member_quantity_ids"]], use.names = FALSE),
    ,
    drop = FALSE
  ]

  out <- data.frame(
    level       = members[["component"]],
    quantity_id = members[["quantity_id"]],
    fixed       = vapply(seq_len(nrow(members)), function(i) {
      .hypothesis_plan_quantity_fixed(
        members[["status"]][[i]],
        members[["fixed_value"]][[i]]
      )
    }, logical(1)),
    stringsAsFactors = FALSE
  )
  out[["dependencies"]] <- lapply(members[["extraction_key"]], function(key) {
    as.character(key[["dependencies"]])
  })

  out
}


# The coefficient scale of the targets of a factor term: on the original
# scale of formula coefficients, the fitted coefficient transform (its
# matrix maps original-coefficient weights to fitted-coordinate weights) and
# the prior-list entries owning the fitted coordinates it combines, whose
# priors the prior densities combine.
.hypothesis_plan_linear_scale <- function(plan, object, entry, coordinates) {

  formula_parameter <- entry[["formula_parameter"]]
  if (plan[["standardized"]] || !is.character(formula_parameter) ||
      length(formula_parameter) != 1L || is.na(formula_parameter) ||
      !nzchar(formula_parameter)) {
    return(list(original = FALSE, prior_parameters = character()))
  }
  transform <- BayesTools::JAGS_formula_coefficient_transform(
    fit          = object[["fit"]],
    parameter    = formula_parameter,
    target_scale = "original"
  )
  term_coordinates <- intersect(coordinates, rownames(transform[["matrix"]]))
  dependencies <- colnames(transform[["matrix"]])[
    colSums(abs(transform[["matrix"]][term_coordinates, , drop = FALSE])) > 0
  ]

  list(
    original          = TRUE,
    formula_parameter = formula_parameter,
    transform         = transform,
    prior_parameters  = if (length(dependencies) > 0L) {
      .hypothesis_brma_target_prior_parameters(
        object,
        stats::setNames(rep(1, length(dependencies)), dependencies)
      )
    } else {
      character()
    }
  )
}


# A linear target's prior density and its weights on the fitted coordinates
# from its weights 'weights' on the coefficients of its posterior (original
# coefficients on the original scale, fitted coordinates otherwise).
.hypothesis_plan_linear_density <- function(plan, object, weights, density) {

  weights <- weights[weights != 0]
  scale   <- plan[["scale"]]
  if (!isTRUE(scale[["original"]])) {
    return(list(prior_density = density, weights = weights))
  }
  matrix <- scale[["transform"]][["matrix"]]
  fitted <- drop(
    matrix(weights, nrow = 1L) %*% matrix[names(weights), , drop = FALSE]
  )
  names(fitted) <- colnames(matrix)

  list(
    prior_density = BayesTools::JAGS_formula_prior_density(
      fit          = object[["fit"]],
      parameter    = scale[["formula_parameter"]],
      weights      = weights,
      target_scale = "original",
      context      = BayesTools::posterior_metadata(
        plan[["draws"]][["samples"]],
        "prior_context"
      )
    ),
    weights = fitted[fitted != 0]
  )
}


.hypothesis_plan_level_target <- function(plan, object, level, value, label) {

  sample  <- plan[["draws"]][["posterior"]][[level]]
  weights <- .iwmde_linear_weights(
    BayesTools::posterior_metadata(sample, "linear_weights")
  )
  if (length(weights) == 0L) {
    # Missing fitted metadata: the target refusal also requires a refit.
    refusal <- .hypothesis_refusal(
      paste0(
        "The linear weights of factor level '", plan[["label"]], "[", level,
        "]' on the fitted coefficients are unavailable. Refit the model ",
        "with the current RoBMA/BayesTools build."
      ),
      "target"
    )
    refusal[["class"]] <- c(refusal[["class"]], .refit_required_class())
    return(.hypothesis_plan_target(
      value         = value,
      label         = label,
      level         = level,
      prior_density = NULL,
      refusal       = refusal
    ))
  }
  density <- .hypothesis_plan_linear_density(
    plan    = plan,
    object  = object,
    weights = weights,
    density = BayesTools::posterior_metadata(sample, "prior_density")
  )
  target <- .hypothesis_plan_target(
    value         = value,
    label         = label,
    level         = level,
    prior_density = density[["prior_density"]],
    atoms         = BayesTools::posterior_metadata(sample, "atoms"),
    spec          = list(
      type          = "linear",
      weights       = density[["weights"]],
      prior_density = density[["prior_density"]]
    )
  )
  target[["point_mass_reason"]] <- .hypothesis_plan_point_mass_reason(object)

  target
}


# A statement on one linear combination of the levels: the scalar linear
# target of BayesTools, which certifies the target's structural atom-freeness
# and joint prior context (classed refusals otherwise).
.hypothesis_plan_combination <- function(plan, object) {

  plan[["route"]] <- "combination"
  target <- tryCatch(
    BayesTools::hypothesis_linear_target(
      posterior  = plan[["draws"]][["posterior"]],
      hypothesis = plan[["hypothesis"]],
      parameter  = plan[["parameter"]]
    ),
    error = .hypothesis_linear_target_condition
  )
  if (inherits(target, "condition")) {
    reason <- conditionMessage(target)
    if (inherits(target, "BayesTools_linear_target_unavailable") &&
        identical(target[["reason"]], "posterior_atoms") &&
        .is_RoBMA(object) && !plan[["conditional"]]) {
      reason <- paste0(
        reason, " The levels of '", plan[["label"]], "' share the null ",
        "component of the model-averaged prior; test the combination within ",
        "the models that include the term with conditional = TRUE."
      )
    }
    plan[["refusal"]] <- .hypothesis_refusal_from_condition(target, reason)
    plan[["group"]] <- paste("combination", plan[["parameter"]], sep = "\r")
    return(plan)
  }
  weights <- .iwmde_linear_weights(target[["weights"]])
  density <- .hypothesis_plan_linear_density(
    plan    = plan,
    object  = object,
    weights = weights,
    density = BayesTools::posterior_metadata(target[["posterior"]], "prior_density")
  )
  BayesTools::posterior_metadata(target[["posterior"]], "prior_density") <-
    density[["prior_density"]]
  plan[["linear_target"]] <- target
  plan[["group"]] <- paste(
    "combination", plan[["parameter"]],
    paste(names(weights), format(weights, digits = 15L), collapse = ","),
    sep = "\r"
  )

  refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = target[["hypothesis"]],
    allow_compound = TRUE
  )
  spec <- list(
    type          = "linear",
    weights       = density[["weights"]],
    prior_density = density[["prior_density"]]
  )
  plan[["targets"]] <- lapply(seq_len(nrow(refs)), function(i) {
    .hypothesis_plan_target(
      value         = refs[["value"]][[i]],
      label         = .hypothesis_plan_side_label(plan, refs[["side"]][[i]]),
      prior_density = density[["prior_density"]],
      atoms         = BayesTools::posterior_metadata(target[["posterior"]], "atoms"),
      spec          = spec
    )
  })

  plan
}


# ---------------------------------------------------------------------------- #
# Random-effect quantities
# ---------------------------------------------------------------------------- #

# Random-effect quantities are tested on their BayesTools mixed posterior
# (catalog support, declared inclusion- and allocation-gate atoms, exact
# prior density). Point statements on a variance, and the point parts of
# statements comparing a variance point with a region, are evaluated through
# its standard deviation with the square display transform, so that a
# variance and its standard deviation give the same Bayes factor with every
# method.
.hypothesis_plan_random <- function(plan, object, cache) {

  parameter <- plan[["parameter"]]
  plan[["kind"]]  <- "random"
  plan[["route"]] <- "random"
  plan[["group"]] <- paste("random", parameter, sep = "\r")
  selected <- .hypothesis_plan_cached(
    cache,
    paste("random_select", parameter, plan[["standardized"]], sep = "\r"),
    function() {
      .brma_random_parameter_select(
        object                    = object,
        parameter                 = parameter,
        standardized_coefficients = plan[["standardized"]]
      )
    }
  )
  label <- selected[["spec"]][["label"]]
  plan[["random"]]  <- selected
  plan[["support"]] <- .brma_random_parameter_support(selected)
  if (identical(selected[["spec"]][["status"]], "structural")) {
    plan[["refusal"]] <- .hypothesis_plan_fixed_refusal(label)
    return(plan)
  }
  samples <- .hypothesis_plan_cached(
    cache,
    paste("random_samples", parameter, plan[["standardized"]],
          plan[["conditional"]], sep = "\r"),
    function() {
      .brma_random_parameter_mixed_posterior(
        object                    = object,
        parameter                 = parameter,
        standardized_coefficients = plan[["standardized"]],
        conditional               = plan[["conditional"]],
        selected                  = selected
      )
    }
  )
  prior_density <- BayesTools::posterior_metadata(samples[[parameter]], "prior_density")
  plan[["draws"]]         <- list(samples = samples, posterior = samples[[parameter]])
  plan[["prior_density"]] <- prior_density
  # Region probabilities need a prior density with deterministic provenance;
  # without one (no density, or a plotting grid such as the SD components of
  # nested allocations) they come from prior draws.
  plan[["prior_draws"]] <- !.hypothesis_plan_density_has_provenance(prior_density)
  if (plan[["prior_draws"]] && plan[["conditional"]]) {
    plan[["refusal"]] <- .hypothesis_refusal(
      paste0(
        "Conditional hypotheses are unavailable for random-effect quantity '",
        label, "' because its prior density is unavailable."
      ),
      "target"
    )
    return(plan)
  }
  if (!plan[["point"]]) {
    return(plan)
  }
  refs <- plan[["refs"]]
  if (any(!refs[["direct"]])) {
    plan[["refusal"]] <- .hypothesis_refusal(
      paste0(
        "Point-null tests for random-effect quantities require a direct ",
        "scalar parameter reference."
      ),
      "statement"
    )
    return(plan)
  }
  if (is.null(prior_density)) {
    plan[["refusal"]] <- .hypothesis_refusal(
      paste0(
        "Point-null Bayes factors are unavailable for random-effect quantity '",
        label, "' because its prior density is unavailable. Use a region or ",
        "directional hypothesis."
      ),
      "target"
    )
    return(plan)
  }

  evaluation <- .hypothesis_plan_random_evaluation(plan, object, cache, selected)
  if (!is.null(evaluation[["refusal"]])) {
    plan[["refusal"]] <- evaluation[["refusal"]]
    return(plan)
  }
  plan[["evaluation"]] <- evaluation
  if (!is.null(evaluation[["parameter"]])) {
    plan[["group"]] <- paste(
      "random", parameter, if (plan[["region"]]) "sd_region" else "sd",
      sep = "\r"
    )
  }
  zero_alternative <- if (any(refs[["value"]] == 0)) {
    .brma_random_parameter_zero_boundary_alternative(object, selected)
  }
  point_mass_reason <- paste0(
    "Point-null Bayes factors at %s are unavailable for random-effect ",
    "quantity '", label, "' because its prior has a point mass at %s. Use a ",
    "region or directional hypothesis; the Component Inclusion table from ",
    "'summary(object)' or 'summary_models(object)' compares the exclusion ",
    "and inclusion of gated components."
  )
  plan[["targets"]] <- lapply(seq_len(nrow(refs)), function(i) {
    value  <- refs[["value"]][[i]]
    label_value <- paste0(label, " = ", format(value, digits = 15L, trim = TRUE))
    # A prior point mass at 0 is refused as such; otherwise 0 is a
    # nonregular product boundary of an allocation-derived SD.
    point_mass <- !is.null(zero_alternative) && value == 0 && identical(
      BayesTools::prior_ordinate_status(prior_density, value, labels = label_value)[["condition"]],
      "BayesTools_point_mass_at_null"
    )
    target <- .hypothesis_plan_target(
      value         = value,
      label         = label_value,
      prior_density = prior_density,
      atoms         = BayesTools::posterior_metadata(samples[[parameter]], "atoms"),
      refusal       = if (!is.null(zero_alternative) && value == 0 && !point_mass) {
        c(
          .hypothesis_refusal(
            paste0(
              "Point-null Bayes factors are unavailable for allocation-derived ",
              "random-effect quantity '", label, "' at 0 because zero is a ",
              "nonregular product boundary of the common scale and allocation ",
              "weight. Test '", zero_alternative, "' to compare omission of ",
              "this component."
            ),
            "target"
          ),
          list(value_specific = TRUE)
        )
      }
    )
    value_text <- format(value, digits = 15L, trim = TRUE)
    target[["point_mass_reason"]] <- sprintf(point_mass_reason, value_text, value_text)
    target
  })
  precomputed <- .hypothesis_plan_random_precomputed(
    plan       = plan,
    object     = object,
    evaluation = evaluation,
    label      = label
  )
  plan[["method_refusals"]] <- precomputed[["refusals"]]
  plan[["density_target"]]  <- precomputed[["density_target"]]

  plan
}


# How point statements on a random-effect quantity are evaluated: variances
# with a standard-deviation counterpart in the catalog are evaluated on the
# standard deviation ('parameter', values mapped by the square root); other
# quantities on themselves (NULL 'parameter'). Of a statement comparing a
# variance point with a region, only the point part is evaluated on the
# standard deviation; the region part is evaluated on the variance draws.
.hypothesis_plan_random_evaluation <- function(plan, object, cache, selected) {

  pair <- .brma_random_parameter_sd_pair(object, selected)
  # At a variance of 0 the square map is singular: the variance is evaluated
  # on its own draws there (its prior ordinate is classified on its own
  # scale).
  if (is.null(pair) || any(plan[["refs"]][["value"]] <= 0)) {
    return(list(parameter = NULL))
  }
  if (length(unique(plan[["refs"]][["value"]])) > 1L) {
    return(list(refusal = .hypothesis_refusal(
      paste0(
        "Point hypotheses on random-effect variance '",
        selected[["spec"]][["label"]], "' are evaluated through its standard ",
        "deviation and must compare one point value per statement."
      ),
      "statement"
    )))
  }
  pair_samples <- .hypothesis_plan_cached(
    cache,
    paste("random_samples", pair, plan[["standardized"]],
          plan[["conditional"]], sep = "\r"),
    function() {
      .brma_random_parameter_mixed_posterior(
        object                    = object,
        parameter                 = pair,
        standardized_coefficients = plan[["standardized"]],
        conditional               = plan[["conditional"]]
      )
    }
  )

  list(
    parameter = pair,
    samples   = pair_samples
  )
}


# The standard deviation whose square is the selected variance quantity (the
# catalog quantity with the same owner and component, and the standard
# deviation evaluator of the variance's extraction key), or NULL.
.brma_random_parameter_sd_pair <- function(object, selected) {

  quantity <- selected[["entry"]][["selection"]][["quantities"]]
  key      <- quantity[["extraction_key"]][[1L]]
  sd_evaluator <- switch(
    as.character(key[["evaluator"]]),
    sd_variance    = "sd",
    allocation_var = "allocation_sd",
    NULL
  )
  if (is.null(sd_evaluator)) {
    return(NULL)
  }
  quantities <- BayesTools::parameter_catalog(object[["fit"]])[["quantities"]]
  same <- function(field) {
    vapply(quantities[[field]], identical, logical(1), y = quantity[[field]][[1L]])
  }
  candidates <- which(
    !quantities[["internal"]] &
      same("formula_parameter") & same("owner_type") & same("owner_name") &
      same("component") &
      vapply(quantities[["extraction_key"]], function(other) {
        is.list(other) && identical(other[["type"]], key[["type"]]) &&
          identical(other[["evaluator"]], sd_evaluator) &&
          identical(other[["random_block"]], key[["random_block"]]) &&
          identical(other[["dependencies"]], key[["dependencies"]])
      }, logical(1))
  )
  if (length(candidates) != 1L) {
    return(NULL)
  }

  .brma_random_parameter_io_labels(quantities[candidates, , drop = FALSE], "selector")
}


# qCMDE/IWMDE refusals of a random-effect point statement: conditional
# statements, quantities without a supported scalar random-component
# coordinate, and values where the public transformation of that coordinate
# is singular. The first two are the causes for which plot() refuses
# qCMDE/IWMDE ("conditional_random", "random_target"): their refusals also
# have its classes.
.hypothesis_plan_random_precomputed <- function(plan, object, evaluation,
                                                label) {

  refused <- function(reason, type = "method") {
    list(refusals = .hypothesis_plan_precomputed_refusals(reason, type))
  }
  refused_density_method <- function(unavailable, type) {
    refusal <- .hypothesis_refusal_density_method(unavailable, type)
    list(refusals = list(qCMDE = refusal, IWMDE = refusal))
  }
  if (plan[["conditional"]]) {
    return(refused_density_method(
      .iwmde_unavailable(
        "Conditional random-effect hypotheses support 'density_method = \"KDE\"' only.",
        "conditional_random"
      ),
      "method"
    ))
  }
  parameter <- if (is.null(evaluation[["parameter"]])) {
    plan[["parameter"]]
  } else {
    evaluation[["parameter"]]
  }
  target <- .brma_random_parameter_density_target(
    object,
    parameter,
    operation = "point hypotheses"
  )
  if (is.null(target[["parameter"]])) {
    # The refusal names the tested quantity (a variance rather than the
    # standard deviation it is evaluated through).
    if (!identical(parameter, plan[["parameter"]])) {
      own <- .brma_random_parameter_density_target(
        object,
        plan[["parameter"]],
        operation = "point hypotheses"
      )
      if (!is.null(own[["reason"]])) {
        target[["reason"]] <- own[["reason"]]
      }
    }
    return(refused_density_method(target, "target"))
  }
  # Values outside the support are refused by their prior ordinates.
  values <- .hypothesis_plan_random_source_values(plan, evaluation)
  values <- values[is.finite(values)]
  display_transform <- target[["display_transform"]]
  if (!is.null(display_transform) && length(values) > 0L) {
    source_values <- BayesTools::parameter_transform_inverse(values, display_transform)
    jacobian      <- BayesTools::parameter_transform_jacobian(
      source_values,
      display_transform
    )
    singular <- !is.finite(source_values) | !is.finite(jacobian) | jacobian <= 0
    if (any(singular)) {
      return(refused(paste0(
        "Point-null Bayes factors are unavailable for random-effect ",
        "quantity '", label, "' at ", values[which(singular)[[1L]]],
        " because its public transformation is singular at that support ",
        "boundary. Use the corresponding directly modeled scale or a ",
        "region hypothesis."
      ), "target"))
    }
  }

  list(refusals = list(), density_target = target)
}


# The tested values on the evaluated quantity's scale: the square roots of
# the values of a variance evaluated through its standard deviation.
.hypothesis_plan_random_source_values <- function(plan, evaluation) {

  values <- plan[["refs"]][["value"]]
  if (is.null(evaluation[["parameter"]])) {
    return(values)
  }

  ifelse(values >= 0, sqrt(pmax(values, 0)), NA_real_)
}


# ---------------------------------------------------------------------------- #
# Marginal means
# ---------------------------------------------------------------------------- #

# The plan of one statement on a marginal-means object. Model-averaged point
# statements use the alternative-conditioned marginal means (their levels
# need not share one conditioning event, so a point statement tests one
# level); other statements use the averaged marginal means. A level's prior
# density is the one BayesTools attached to the marginal mean; linear
# combinations of levels are the scalar linear targets of
# BayesTools::hypothesis_linear_target(). qCMDE/IWMDE ordinates are computed
# from the stored source model.
.hypothesis_plan_marginal_means <- function(object, statement,
                                            parameter = NULL, cache = NULL) {

  if (!inherits(statement, "BayesTools_hypothesis_ast")) {
    statement <- BayesTools::hypothesis_parse(statement)
  }
  selected <- .hypothesis_marginal_means_select_parameter(
    object     = object,
    hypothesis = statement,
    parameter  = parameter
  )
  parameter <- selected[["parameter"]]
  plan <- .hypothesis_plan_new(
    statement  = statement,
    hypothesis = .hypothesis_brma_rewrite(
      hypothesis = statement,
      aliases    = selected[["aliases"]],
      parameter  = parameter
    ),
    parameter  = parameter,
    label      = .hypothesis_brma_alias_label(selected[["aliases"]], parameter),
    component  = "marginal_means",
    selected   = selected
  )
  plan[["kind"]]              <- "marginal_mean"
  plan[["route"]]             <- "levels"
  plan[["group"]]             <- paste("levels", parameter, sep = "\r")
  plan[["capability_object"]] <- object[["source_object"]]
  model_averaged <- isTRUE(object[["model_averaged"]]) ||
    inherits(object[["source_object"]], "RoBMA")
  plan[["inference_type"]] <- if (model_averaged && plan[["point"]]) {
    "conditional"
  } else {
    "averaged"
  }
  if (model_averaged && plan[["point"]] && plan[["region"]]) {
    plan[["refusal"]] <- .hypothesis_refusal(
      paste0(
        "Model-averaged marginal-means hypotheses cannot mix point and ",
        "region events because they require different posterior ",
        "conditioning. Use a pure point-null or a pure region hypothesis."
      ),
      "statement"
    )
    return(.hypothesis_plan_finish(plan, object))
  }
  if (model_averaged && plan[["point"]] && length(unique(
    .hypothesis_marginal_means_ast_levels(plan[["hypothesis"]], parameter)
  )) > 1L) {
    plan[["refusal"]] <- .hypothesis_refusal(
      paste0(
        "Point hypotheses spanning multiple marginal-means levels are not ",
        "supported because the levels need not share one conditioning event."
      ),
      "statement"
    )
    return(.hypothesis_plan_finish(plan, object))
  }
  posterior <- object[["inference"]][[plan[["inference_type"]]]][[parameter]]
  if (is.null(posterior)) {
    plan[["refusal"]] <- .hypothesis_refusal(
      paste0(
        "'marginal_means' object does not contain ",
        plan[["inference_type"]], " marginal means."
      ),
      "target"
    )
    return(.hypothesis_plan_finish(plan, object))
  }
  plan[["draws"]] <- list(posterior = posterior)

  plan[["method_refusals"]] <- list(normal = .hypothesis_refusal(
    "'normal' density_method is not supported for marginal-means hypotheses.",
    "method"
  ))
  # Without its fitted source model, qCMDE/IWMDE are refused before that
  # model's capability is checked, as the marginal-means plot() does.
  if (is.null(object[["source_object"]]) ||
      !inherits(object[["source_object"]], "brma") ||
      is.null(object[["source_object"]][["fit"]])) {
    plan[["source_refusals"]] <- .hypothesis_plan_precomputed_refusals(paste0(
      "The marginal-means object does not contain the source fitted brma ",
      "object needed to compute qCMDE/IWMDE ordinates."
    ), "target")
  }

  refs <- plan[["refs"]]
  plan <- if (nrow(refs) > 0L && any(!refs[["direct"]])) {
    .hypothesis_plan_marginal_combination(plan)
  } else {
    .hypothesis_plan_marginal_levels(plan)
  }

  .hypothesis_plan_finish(plan, object)
}


# The refusal of a model-averaged marginal-means request as a whole, or NULL.
# Its statements are evaluated on one posterior conditioning (point
# statements on the alternative-conditioned marginal mean of one level,
# region statements on the averaged marginal means), so a request cannot
# combine point and region statements or point statements on several levels.
.hypothesis_plan_marginal_means_request_refusal <- function(plans, object,
                                                            parameter) {

  model_averaged <- isTRUE(object[["model_averaged"]]) ||
    inherits(object[["source_object"]], "RoBMA")
  if (!model_averaged) {
    return(NULL)
  }
  point <- vapply(plans, `[[`, logical(1), "point")
  if (any(point) && !all(point)) {
    return(.hypothesis_refusal(
      paste0(
        "A model-averaged marginal-means hypothesis request cannot mix ",
        "point-null and region statements because they require different ",
        "posterior conditioning."
      ),
      "statement"
    ))
  }
  levels <- unique(.hypothesis_marginal_means_ast_levels(
    hypothesis = .hypothesis_plan_group_ast(plans, "hypothesis"),
    parameter  = parameter
  ))
  if (all(point) && length(levels) > 1L) {
    return(.hypothesis_refusal(
      paste0(
        "Point hypotheses spanning multiple marginal-means levels are not ",
        "supported because the levels need not share one conditioning event."
      ),
      "statement"
    ))
  }

  NULL
}


# A marginal mean the fitted model fixes: its declared atoms carry the whole
# posterior mass.
.hypothesis_plan_draws_fixed <- function(sample) {

  atoms <- BayesTools::posterior_metadata(sample, "atoms")
  mass  <- atoms[["mass"]]
  isTRUE(atoms[["declared"]]) && length(mass) > 0L &&
    abs(sum(mass) - 1) <= sqrt(.Machine[["double.eps"]])
}


.hypothesis_plan_marginal_label <- function(plan, level) {

  if (is.na(level)) {
    return(plan[["label"]])
  }

  paste0(plan[["label"]], "[", level, "]")
}


.hypothesis_plan_marginal_levels <- function(plan) {

  posterior <- plan[["draws"]][["posterior"]]
  parameter <- plan[["parameter"]]
  levels    <- if (is.list(posterior)) names(posterior) else NULL
  occurrences <- BayesTools::hypothesis_symbols(
    plan[["hypothesis"]],
    occurrences = TRUE
  )
  referenced <- unique(occurrences[["level"]][occurrences[["parameter"]] == parameter])
  referenced <- referenced[!is.na(referenced)]
  unknown <- setdiff(referenced, levels)
  if (length(unknown) > 0L) {
    plan[["refusal"]] <- .hypothesis_unknown_level_refusal(
      level     = unknown[[1L]],
      parameter = parameter
    )
    return(plan)
  }
  sample_of <- function(level) {
    if (is.na(level)) posterior else posterior[[level]]
  }
  tested <- if (length(referenced) > 0L) referenced else if (!is.list(posterior)) NA_character_
  fixed  <- tested[vapply(tested, function(level) {
    .hypothesis_plan_draws_fixed(sample_of(level))
  }, logical(1))]
  if (length(fixed) > 0L) {
    plan[["refusal"]] <- .hypothesis_plan_fixed_refusal(
      .hypothesis_plan_marginal_label(plan, fixed[[1L]])
    )
    return(plan)
  }

  refs <- plan[["refs"]]
  # A factor marginal mean with one level resolves statements without a
  # level to that level.
  if (is.list(posterior) && length(posterior) == 1L) {
    refs[["level"]][is.na(refs[["level"]])] <- levels
  }
  if (is.list(posterior) && any(is.na(refs[["level"]]))) {
    plan[["method_refusals"]] <- c(
      plan[["method_refusals"]],
      .hypothesis_plan_precomputed_refusals(paste0(
        "qCMDE/IWMDE point hypotheses for marginal-means factor parameters ",
        "must specify a level, e.g. '", plan[["label"]], "[level] = ",
        refs[["value"]][is.na(refs[["level"]])][[1L]], "'."
      ))
    )
  }
  plan[["targets"]] <- unlist(lapply(seq_len(nrow(refs)), function(i) {
    ref_levels <- if (is.na(refs[["level"]][[i]]) && is.list(posterior)) {
      levels
    } else {
      refs[["level"]][[i]]
    }
    label <- .hypothesis_plan_side_label(plan, refs[["side"]][[i]])
    lapply(ref_levels, function(level) {
      sample  <- sample_of(level)
      density <- BayesTools::posterior_metadata(sample, "prior_density")
      .hypothesis_plan_target(
        value         = refs[["value"]][[i]],
        label         = label,
        level         = level,
        prior_density = density,
        atoms         = BayesTools::posterior_metadata(sample, "atoms"),
        refusal       = if (.hypothesis_plan_draws_fixed(sample)) {
          .hypothesis_plan_fixed_refusal(.hypothesis_plan_marginal_label(plan, level))
        } else if (is.null(density)) {
          .hypothesis_refusal(
            paste0(
              "The prior density of marginal mean '",
              .hypothesis_plan_marginal_label(plan, level),
              "' is unavailable, so point hypotheses on it have no ",
              "Savage-Dickey Bayes factor. Use a region or directional ",
              "hypothesis."
            ),
            "target"
          )
        }
      )
    })
  }), recursive = FALSE)
  if (length(plan[["targets"]]) > 0L) {
    plan[["support"]] <- .hypothesis_plan_support_bounds(
      sample_of(plan[["targets"]][[1L]][["level"]])
    )
  }

  plan
}


# The condition with which BayesTools::hypothesis_linear_target() refuses a
# statement: its unclassed errors refuse the form of the statement (e.g. a
# nonlinear expression of the levels), which is restated
# ("RoBMA_hypothesis_statement"); an unknown level (a BayesTools resolution
# error such as "BayesTools_parameter_not_found") is a statement error with
# BayesTools' classes, as at the other entry points
# (.hypothesis_not_found_class()); other classed BayesTools conditions
# (linear targets BayesTools cannot certify,
# "BayesTools_linear_target_unavailable", or refusals of the fitted
# metadata) keep their classes.
.hypothesis_linear_target_condition <- function(condition) {

  if (inherits(condition, "BayesTools_parameter_resolution_error")) {
    class(condition) <- unique(c("RoBMA_hypothesis_statement", class(condition)))
    return(condition)
  }
  if (any(startsWith(class(condition), "BayesTools_"))) {
    return(condition)
  }

  structure(
    class = c("RoBMA_hypothesis_statement", "error", "condition"),
    list(message = conditionMessage(condition), call = NULL)
  )
}


.hypothesis_plan_marginal_combination <- function(plan) {

  plan[["route"]] <- "combination"
  plan[["group"]] <- paste("combination", plan[["parameter"]], sep = "\r")
  target <- tryCatch(
    BayesTools::hypothesis_linear_target(
      posterior  = plan[["draws"]][["posterior"]],
      hypothesis = plan[["hypothesis"]],
      parameter  = plan[["parameter"]]
    ),
    error = .hypothesis_linear_target_condition
  )
  if (inherits(target, "condition")) {
    plan[["refusal"]] <- .hypothesis_refusal_from_condition(target)
    return(plan)
  }
  weights <- .iwmde_linear_weights(target[["weights"]])
  density <- BayesTools::posterior_metadata(target[["posterior"]], "prior_density")
  plan[["linear_target"]] <- target
  plan[["group"]] <- paste(
    "combination", plan[["parameter"]],
    paste(names(weights), format(weights, digits = 15L), collapse = ","),
    sep = "\r"
  )
  refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = target[["hypothesis"]],
    allow_compound = TRUE
  )
  plan[["targets"]] <- lapply(seq_len(nrow(refs)), function(i) {
    .hypothesis_plan_target(
      value         = refs[["value"]][[i]],
      label         = .hypothesis_plan_side_label(plan, refs[["side"]][[i]]),
      prior_density = density,
      atoms         = BayesTools::posterior_metadata(target[["posterior"]], "atoms"),
      spec          = list(
        type          = "linear",
        weights       = weights,
        prior_density = density
      )
    )
  })

  plan
}
