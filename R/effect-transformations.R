# ============================================================================ #
# effect-transformations.R
# ============================================================================ #
#
# Effect-size transformations for posterior draws and plotting.
#
# Measures follow metafor naming: OR, RR, HR, and IRR are log-scale measures.
# Use transform = "EXP" to exponentiate a log-scale measure for display.
#
# Transformation objects use the BayesTools convention:
# - fun: forward transformation
# - inv: inverse transformation
# - jac: derivative of fun(x) with respect to x
#
# ============================================================================ #


.effect_known_measures <- function() {

  return(c("SMD", "COR", "ZCOR", "OR", "RR", "HR", "IRR", "RD", "GEN"))
}

.effect_core_measures <- function() {

  return(c("SMD", "COR", "ZCOR", "OR"))
}

.effect_log_measures <- function() {

  return(c("OR", "RR", "HR", "IRR"))
}

.effect_measure_aliases <- function() {

  return(c(
    "D"          = "SMD",
    "COHENS_D"   = "SMD",
    "COHENSD"    = "SMD",
    "R"          = "COR",
    "CORRELATION" = "COR",
    "FISHERS_Z"  = "ZCOR",
    "FISHERSZ"   = "ZCOR",
    "FISHER_Z"   = "ZCOR",
    "FISHERZ"    = "ZCOR",
    "Z"          = "ZCOR",
    "LOGOR"      = "OR",
    "LOGRR"      = "RR",
    "LOGHR"      = "HR",
    "LOGIRR"     = "IRR",
    "NONE"       = "GEN"
  ))
}

.normalize_effect_measure <- function(measure, argument, allow_NULL = FALSE) {

  if (is.null(measure)) {
    if (allow_NULL) {
      return(NULL)
    }
    stop("The '", argument, "' argument must not be NULL.", call. = FALSE)
  }

  BayesTools::check_char(measure, argument, check_length = 1, allow_NA = FALSE)

  measure_key <- toupper(gsub("[^A-Za-z0-9]", "", measure))
  aliases     <- .effect_measure_aliases()

  if (measure_key %in% names(aliases)) {
    measure_key <- aliases[[measure_key]]
  }

  if (!measure_key %in% .effect_known_measures()) {
    stop(
      "Unknown effect-size measure in '", argument, "'. Available measures are: ",
      paste(.effect_known_measures(), collapse = ", "), ".",
      call. = FALSE
    )
  }

  return(measure_key)
}

.normalize_effect_transform <- function(transform, allow_log = FALSE) {

  if (is.null(transform)) {
    return("identity")
  }

  BayesTools::check_char(transform, "transform", check_length = 1, allow_NA = FALSE)

  transform_key <- toupper(gsub("[^A-Za-z0-9]", "", transform))
  transform_map <- c(
    NONE     = "identity",
    IDENTITY = "identity",
    ID       = "identity",
    EXP      = "EXP"
  )
  if (allow_log) {
    transform_map <- c(transform_map, LOG = "LOG")
  }
  if (!transform_key %in% names(transform_map)) {
    available <- if (allow_log) {
      "'EXP', 'LOG', and 'identity'"
    } else {
      "'EXP' and 'identity'"
    }
    stop(
      "Unknown 'transform'. Available options are ", available, ".",
      call. = FALSE
    )
  }

  return(unname(transform_map[[transform_key]]))
}

.identity_effect_transformation <- function() {

  return(list(
    fun = function(x) x,
    inv = function(x) x,
    jac = function(x) rep(1, length(x))
  ))
}

.exp_effect_transformation <- function() {

  return(list(
    fun = exp,
    inv = log,
    jac = exp
  ))
}

.log_plot_transformation <- function() {

  return(list(
    fun = log,
    inv = exp,
    jac = function(x) 1 / x
  ))
}

.d_to_cor <- function(x) {

  out <- x / sqrt(x^2 + 4)

  out[is.infinite(x) & x > 0] <-  1
  out[is.infinite(x) & x < 0] <- -1

  return(out)
}

.cor_to_z <- function(x) {

  if (!is.numeric(x) || any(!is.finite(x)) || any(x < -1 | x > 1)) {
    stop(
      "Correlation values must be numeric, finite, and within [-1, 1].",
      call. = FALSE
    )
  }

  return(atanh(x))
}

.effect_basic_transformations <- function() {

  return(list(
    "SMD->COR" = list(
      fun = .d_to_cor,
      inv = function(x) 2 * x / sqrt(1 - x^2),
      jac = function(x) 4 / (x^2 + 4)^(3 / 2)
    ),
    "COR->SMD" = list(
      fun = function(x) 2 * x / sqrt(1 - x^2),
      inv = .d_to_cor,
      jac = function(x) 2 / (1 - x^2)^(3 / 2)
    ),
    "COR->ZCOR" = list(
      fun = .cor_to_z,
      inv = tanh,
      jac = function(x) 1 / (1 - x^2)
    ),
    "ZCOR->COR" = list(
      fun = tanh,
      inv = .cor_to_z,
      jac = function(x) 1 - tanh(x)^2
    ),
    "SMD->OR" = list(
      fun = function(x) x * pi / sqrt(3),
      inv = function(x) x * sqrt(3) / pi,
      jac = function(x) rep(pi / sqrt(3), length(x))
    ),
    "OR->SMD" = list(
      fun = function(x) x * sqrt(3) / pi,
      inv = function(x) x * pi / sqrt(3),
      jac = function(x) rep(sqrt(3) / pi, length(x))
    )
  ))
}

.compose_effect_transformations <- function(first, second) {

  force(first)
  force(second)

  return(list(
    fun = function(x) second[["fun"]](first[["fun"]](x)),
    inv = function(x) first[["inv"]](second[["inv"]](x)),
    jac = function(x) {

      first[["jac"]](x) * second[["jac"]](first[["fun"]](x))
    }
  ))
}

.effect_neighbor_edges <- function(edges, node) {

  edge_names <- names(edges)
  keep       <- startsWith(edge_names, paste0(node, "->"))

  return(edge_names[keep])
}

.effect_measure_transformation <- function(from, to) {

  if (from == to) {
    return(.identity_effect_transformation())
  }

  if (!all(c(from, to) %in% .effect_core_measures())) {
    stop(
      "Effect-size transformation from ", from, " to ", to,
      " is not available. Only SMD, COR, ZCOR, and OR can be converted. ",
      "RR, HR, IRR, RD, and GEN can only be shown on their fitted measure; ",
      "use transform = 'EXP' for RR, HR, IRR, and OR ratios.",
      call. = FALSE
    )
  }

  edges   <- .effect_basic_transformations()
  queue   <- list(list(
    node           = from,
    transformation = .identity_effect_transformation()
  ))
  visited <- from

  while (length(queue) > 0L) {

    current <- queue[[1]]
    queue   <- queue[-1]

    for (edge_name in .effect_neighbor_edges(edges, current[["node"]])) {

      next_node <- sub("^.*->", "", edge_name)
      if (next_node %in% visited) {
        next
      }

      next_transformation <- .compose_effect_transformations(
        first  = current[["transformation"]],
        second = edges[[edge_name]]
      )

      if (next_node == to) {
        return(next_transformation)
      }

      visited <- c(visited, next_node)
      queue[[length(queue) + 1L]] <- list(
        node           = next_node,
        transformation = next_transformation
      )
    }
  }

  stop("No effect-size transformation path found from ", from, " to ", to, ".",
       call. = FALSE)
}

.effect_output_setup <- function(object, output_measure = NULL, transform = NULL) {

  return(.effect_output_setup_measure(
    input_measure  = .measure(object),
    output_measure = output_measure,
    transform      = transform
  ))
}

.effect_output_setup_measure <- function(input_measure, output_measure = NULL,
                                         transform = NULL) {

  input_measure  <- .normalize_effect_measure(input_measure, "input_measure")
  output_measure <- .normalize_effect_measure(output_measure, "output_measure", allow_NULL = TRUE)
  transform      <- .normalize_effect_transform(transform)

  requested <- !is.null(output_measure) || transform != "identity"

  if (is.null(output_measure)) {
    output_measure <- input_measure
  }

  if (transform == "EXP" && !output_measure %in% .effect_log_measures()) {
    stop(
      "transform = 'EXP' is only available for log-scale measures ",
      "OR, RR, HR, and IRR.",
      call. = FALSE
    )
  }

  measure_transformation <- .effect_measure_transformation(
    from = input_measure,
    to   = output_measure
  )

  if (transform == "EXP") {
    transformation <- .compose_effect_transformations(
      first  = measure_transformation,
      second = .exp_effect_transformation()
    )
  } else {
    transformation <- measure_transformation
  }

  active <- input_measure != output_measure || transform != "identity"

  output <- list(
    input_measure  = input_measure,
    output_measure = output_measure,
    transform      = transform,
    requested      = requested,
    active         = active,
    transformation = transformation,
    label          = .effect_output_label(output_measure, transform),
    note           = .effect_output_note(input_measure, output_measure, transform)
  )

  return(output)
}

.plot_output_setup <- function(object, parameter, parameter_entry,
                               output_measure = NULL, transform = NULL) {

  transform <- .normalize_effect_transform(transform, allow_log = TRUE)

  formula_parameter     <- parameter_entry[["formula_parameter"]]
  has_formula_parameter <- (
    is.character(formula_parameter) &&
    length(formula_parameter) == 1L &&
    !is.na(formula_parameter) &&
    nzchar(formula_parameter)
  )
  is_formula_coef       <- (
    identical(parameter_entry[["role"]], "fixed_coefficient") ||
    identical(parameter_entry[["role"]], "formula_coefficient_group")
  )

  is_effect_intercept <- .is_effect_location_parameter(parameter)
  is_location_coef    <- identical(parameter_entry[["component"]], "mods")
  is_scale_coef       <- (
    identical(parameter_entry[["component"]], "scale") &&
    is_formula_coef &&
    has_formula_parameter &&
    !identical(parameter_entry[["term"]], "intercept")
  )
  is_scale_intercept  <- (
    identical(parameter_entry[["component"]], "scale") &&
    identical(parameter_entry[["term"]], "intercept")
  )
  random_quantity <- parameter_entry[["quantity"]]
  is_random_multiplier <- (
    identical(parameter_entry[["component"]], "random") &&
    is.character(random_quantity) && length(random_quantity) == 1L &&
    !is.na(random_quantity) &&
    random_quantity %in% c("var_mult", "sd_mult")
  )

  if (!is.null(output_measure) && !is_effect_intercept) {
    stop(
      "'output_measure' is only available for effect-size location parameters ",
      "('mu' or the meta-regression intercept).",
      call. = FALSE
    )
  }

  if (identical(transform, "EXP")) {
    if (is_scale_coef) {
      return(list(
        input_measure  = .measure(object),
        output_measure = .measure(object),
        transform      = transform,
        requested      = TRUE,
        active         = TRUE,
        transformation = .exp_effect_transformation(),
        label          = "multiplicative scale",
        note           = NULL
      ))
    }

    if (!is_effect_intercept && !is_location_coef) {
      stop(
        "transform = 'EXP' is only available for effect-size location, ",
        "location meta-regression coefficients, and scale-regression ",
        "coefficients.",
        call. = FALSE
      )
    }

    return(.effect_output_setup(
      object         = object,
      output_measure = output_measure,
      transform      = transform
    ))
  }

  if (identical(transform, "LOG")) {
    if (!is_scale_intercept && !is_random_multiplier) {
      stop(
        "transform = 'LOG' is only available for positive heterogeneity ",
        "intercepts and random-effect variance or SD multipliers.",
        call. = FALSE
      )
    }

    return(list(
      input_measure  = .measure(object),
      output_measure = .measure(object),
      transform      = transform,
      requested      = TRUE,
      active         = TRUE,
      transformation = .log_plot_transformation(),
      label          = "log scale",
      note           = NULL
    ))
  }

  return(.effect_output_setup(
    object         = object,
    output_measure = output_measure,
    transform      = transform
  ))
}

.effect_output_requested <- function(effect_transform) {

  if (is.null(effect_transform)) {
    return(FALSE)
  }

  return(isTRUE(effect_transform[["requested"]]))
}

.effect_output_active <- function(effect_transform) {

  if (is.null(effect_transform)) {
    return(FALSE)
  }

  return(isTRUE(effect_transform[["active"]]))
}

.effect_output_label <- function(output_measure, transform) {

  if (transform == "EXP") {
    return(switch(
      output_measure,
      "OR"  = "odds ratio",
      "RR"  = "risk ratio",
      "HR"  = "hazard ratio",
      "IRR" = "incidence rate ratio"
    ))
  }

  return(switch(
    output_measure,
    "SMD"  = "standardized mean difference",
    "COR"  = "correlation",
    "ZCOR" = "Fisher's z",
    "OR"   = "log odds ratio",
    "RR"   = "log risk ratio",
    "HR"   = "log hazard ratio",
    "IRR"  = "log incidence rate ratio",
    "RD"   = "risk difference",
    "GEN"  = "effect size"
  ))
}

.effect_measure_label <- function(measure) {

  return(.effect_output_label(measure, "identity"))
}

.effect_output_note <- function(input_measure, output_measure, transform) {

  if (input_measure == output_measure && transform == "identity") {
    return(NULL)
  }

  input_label  <- .effect_measure_label(input_measure)
  output_label <- .effect_output_label(output_measure, transform)

  if (transform == "EXP" && input_measure == output_measure) {
    return(paste0(
      "Effect estimates are summarized on the ", output_label,
      " scale using EXP on the ", input_label, " measure."
    ))
  }

  if (transform == "EXP") {
    return(paste0(
      "Effect estimates are transformed from ", input_label, " to ",
      .effect_measure_label(output_measure), " and summarized on the ",
      output_label, " scale using EXP."
    ))
  }

  return(paste0(
    "Effect estimates are transformed from ", input_label,
    " to ", output_label, "."
  ))
}

.effect_output_title <- function(title, effect_transform) {

  if (!.effect_output_active(effect_transform)) {
    return(title)
  }

  suffix <- paste0(" (", effect_transform[["label"]], ")")

  if (grepl(":$", title)) {
    title <- sub(":$", paste0(suffix, ":"), title)
  } else {
    title <- paste0(title, suffix)
  }

  return(title)
}

.effect_plot_transformation <- function(effect_transform) {

  if (!.effect_output_active(effect_transform)) {
    return(NULL)
  }

  return(effect_transform[["transformation"]])
}

.is_effect_location_parameter <- function(parameter) {

  return(parameter %in% c("mu", "mu_intercept"))
}

.transform_effect_matrix <- function(samples, effect_transform) {

  if (!.effect_output_active(effect_transform)) {
    return(samples)
  }

  dim_names <- dimnames(samples)
  samples   <- effect_transform[["transformation"]][["fun"]](samples)

  if (!is.matrix(samples)) {
    samples <- as.matrix(samples)
  }

  dimnames(samples) <- dim_names

  return(samples)
}

.transform_effect_vector <- function(samples, effect_transform) {

  if (!.effect_output_active(effect_transform)) {
    return(samples)
  }

  return(effect_transform[["transformation"]][["fun"]](samples))
}

.new_effect_brma_samples <- function(samples, n_chains, n_iter, title,
                                     component = "location",
                                     probs = c(.025, .975), data = NULL,
                                     effect_transform = NULL,
                                     prediction_samples = NULL) {

  if (!is.null(effect_transform)) {
    samples <- .transform_effect_matrix(samples, effect_transform)
    if (!is.null(prediction_samples)) {
      prediction_samples <- .transform_effect_matrix(
        prediction_samples,
        effect_transform
      )
    }
    title <- .effect_output_title(title, effect_transform)
  }

  return(.new_brma_samples(
    samples            = samples,
    n_chains           = n_chains,
    n_iter             = n_iter,
    title              = title,
    component          = component,
    probs              = probs,
    data               = data,
    effect_transform   = effect_transform,
    prediction_samples = prediction_samples
  ))
}

.transform_marginal_samples_effect <- function(samples, effect_transform) {

  if (!.effect_output_active(effect_transform)) {
    return(samples)
  }

  for (parameter in names(samples)) {
    samples[[parameter]] <- .transform_marginal_posterior_effect(
      samples          = samples[[parameter]],
      effect_transform = effect_transform
    )
  }

  return(samples)
}

.transform_marginal_posterior_effect <- function(samples, effect_transform) {

  sample_attributes <- attributes(samples)

  for (i in seq_along(samples)) {
    samples[[i]] <- .transform_effect_vector(
      samples           = samples[[i]],
      effect_transform  = effect_transform
    )
  }

  attributes(samples) <- sample_attributes

  return(samples)
}
