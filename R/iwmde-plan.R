# ============================================================================ #
# IWMDE Estimate Plan Helpers
# ============================================================================ #

.iwmde_plan <- function(context, parameter, density_method, density_control,
                        outputs = c("density", "ordinate"), values = NULL,
                        parameter_spec = NULL, metadata = NULL) {

  context        <- .iwmde_context_ensure_caches(context)
  density_method <- .density_method_normalize_precomputed(density_method)
  outputs <- unique(match.arg(
    outputs,
    c("density", "ordinate"),
    several.ok = TRUE
  ))
  density_control <- .iwmde_density_control_resolve(
    density_method  = density_method,
    density_control = density_control,
    purpose         = if ("ordinate" %in% outputs) {
      "ordinate"
    } else {
      "density"
    }
  )
  if (is.null(density_control[["normalization_points"]])) {
    density_control[["normalization_points"]] <- max(
      50L,
      density_control[["n_points"]]
    )
  }
  row_budget <- if ("ordinate" %in% outputs) {
    density_control[["samples"]]
  } else {
    density_control[["max_samples"]]
  }

  values <- as.numeric(values)
  values <- values[is.finite(values)]
  ordinate_values <- if ("ordinate" %in% outputs && length(values) > 0L) {
    values[[1L]]
  } else {
    NULL
  }

  parameter_spec <- .iwmde_prepare_prior_ordinates(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    values         = ordinate_values
  )
  prior_ordinates <- parameter_spec[["prior_ordinates"]]
  parameter_spec[["prior_ordinates"]] <- NULL
  ordinate_warnings <- .iwmde_ordinate_prior_warnings(
    parameter       = parameter,
    prior_ordinates = prior_ordinates
  )
  parameter_spec        <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  method                <- .density_method_iwmde_estimator(density_method)
  target_key            <- .iwmde_target_key(parameter, parameter_spec)
  stored_parameter_spec <- .iwmde_plan_parameter_spec(parameter_spec)

  plan <- list(
    target = .iwmde_plan_target(
      parameter      = parameter,
      parameter_spec = parameter_spec,
      metadata       = metadata,
      target_key     = target_key
    ),
    parameter_spec = stored_parameter_spec,
    execution_spec = .iwmde_plan_execution_spec(parameter_spec),
    outputs = list(
      need_density     = "density" %in% outputs,
      need_ordinate    = "ordinate" %in% outputs,
      requested_values = values,
      display_grid     = density_control[["display_grid"]]
    ),
    control = density_control,
    row_budget = row_budget,
    method  = method,
    density_method = density_method,
    source_fingerprint = .iwmde_source_fingerprint(context),
    prior_ordinates     = prior_ordinates,
    ordinate_warnings   = ordinate_warnings
  )
  plan <- .iwmde_plan_prepare_contract(context, plan)
  plan <- .iwmde_new_plan(plan)
  plan[["plan_key"]] <- .iwmde_hash("iwmde_plan", .iwmde_plan_key_payload(plan))

  return(plan)
}


.iwmde_plan_prepare_contract <- function(context, plan) {

  parameter      <- plan[["target"]][["parameter"]]
  parameter_spec <- plan[["parameter_spec"]]

  unavailable_reason <- .iwmde_context_unavailable_reason(context)
  if (!is.null(unavailable_reason)) {
    plan[["status"]] <- "unsupported"
    plan[["reason"]] <- unavailable_reason
    plan[["rows"]]   <- .iwmde_plan_empty_rows()
    return(plan)
  }

  if (identical(parameter_spec[["status"]], "unsupported")) {
    plan[["status"]] <- "unsupported"
    plan[["reason"]] <- parameter_spec[["reason"]]
    plan[["rows"]]   <- .iwmde_plan_empty_rows()
    return(plan)
  }

  boundary_reason <- .iwmde_known_v_tau_zero_boundary_reason(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    values         = plan[["outputs"]][["requested_values"]]
  )
  if (!is.null(boundary_reason)) {
    plan[["status"]] <- "unsupported"
    plan[["reason"]] <- boundary_reason
    plan[["rows"]]   <- .iwmde_plan_empty_rows()
    return(plan)
  }

  posterior_values <- tryCatch(
    .iwmde_parameter_values(
      context        = context,
      parameter      = parameter,
      parameter_spec = parameter_spec
    ),
    error = function(e) {
      if (inherits(e, "iwmde_construction_error")) {
        stop(e)
      }
      .iwmde_stop_construction_failure(
        estimator = plan[["method"]],
        parameter = parameter,
        rows      = integer(),
        stage     = "target-draw evaluation",
        detail    = conditionMessage(e)
      )
    }
  )
  if (!is.numeric(posterior_values) ||
      length(posterior_values) != nrow(context[["posterior_samples"]])) {
    .iwmde_stop_construction_failure(
      estimator = plan[["method"]],
      parameter = parameter,
      rows      = integer(),
      stage     = "target-draw evaluation",
      detail    = "the target draw vector has an invalid length or type"
    )
  }
  finite <- is.finite(posterior_values)
  condition_rows <- .iwmde_parameter_condition_rows(
    context        = context,
    parameter_spec = parameter_spec
  )
  if (!is.logical(condition_rows) ||
      length(condition_rows) != length(posterior_values) ||
      anyNA(condition_rows)) {
    .iwmde_stop_construction_failure(
      estimator = plan[["method"]],
      parameter = parameter,
      rows      = integer(),
      stage     = "target-condition validation",
      detail    = "the target conditioning mask has an invalid length or type"
    )
  }
  invalid_target_rows <- which(condition_rows & !finite)
  if (length(invalid_target_rows) > 0L) {
    .iwmde_stop_construction_failure(
      estimator = plan[["method"]],
      parameter = parameter,
      rows      = invalid_target_rows,
      stage     = "target-draw validation",
      detail    = "a conditioned target draw was not finite"
    )
  }
  finite_rows <- finite & condition_rows
  if (!any(finite_rows)) {
    plan[["status"]] <- "unsupported"
    plan[["reason"]] <- "posterior samples are not finite"
    plan[["rows"]]   <- .iwmde_plan_empty_rows()
    return(plan)
  }

  component <- .iwmde_parameter_components(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec
  )
  component <- .iwmde_restrict_parameter_component(component, finite_rows)
  if (!any(component[["active"]])) {
    plan[["status"]] <- "point_only"
    plan[["rows"]] <- .iwmde_plan_rows(
      posterior_values = posterior_values,
      finite_rows      = finite_rows,
      component        = component
    )
    return(plan)
  }

  continuous_rows <- which(component[["active"]] & finite_rows)
  chain_coverage <- .iwmde_plan_chain_coverage(
    context       = context,
    eligible_rows = continuous_rows,
    selected_rows = continuous_rows
  )
  missing_chains <- chain_coverage[["expected_chain_ids"]][
    chain_coverage[["eligible_chain_counts"]] == 0L
  ]
  if (length(missing_chains) > 0L) {
    .iwmde_stop_construction_failure(
      estimator = plan[["method"]],
      parameter = parameter,
      rows      = integer(),
      stage     = "chain-coverage validation",
      detail    = paste0(
        "no eligible continuous row remains in fitted chain(s): ",
        paste(missing_chains, collapse = ", ")
      ),
      chain_coverage = chain_coverage
    )
  }
  if (length(continuous_rows) < 20L) {
    plan[["status"]] <- "unsupported"
    plan[["reason"]] <- "fewer than 20 continuous active samples"
    plan[["rows"]]   <- .iwmde_plan_empty_rows()
    return(plan)
  }

  candidate_rows <- continuous_rows
  if (length(candidate_rows) > plan[["row_budget"]]) {
    candidate_rows <- .nested_srs_rows(
      rows        = candidate_rows,
      max_samples = plan[["row_budget"]]
    )
  }
  chain_coverage <- .iwmde_plan_chain_coverage(
    context       = context,
    eligible_rows = continuous_rows,
    selected_rows = candidate_rows
  )

  candidate_values <- posterior_values[candidate_rows]
  if (all(candidate_values == candidate_values[[1L]])) {
    plan[["status"]] <- "unsupported"
    plan[["reason"]] <- "active samples have zero variance"
    plan[["rows"]]   <- .iwmde_plan_empty_rows()
    return(plan)
  }

  baseline_contract <- .iwmde_plan_baseline_contract(
    context          = context,
    plan             = plan,
    candidate_rows   = candidate_rows,
    candidate_values = candidate_values
  )

  support <- .iwmde_parameter_support(
    context        = context,
    parameter      = parameter,
    rows           = candidate_rows,
    parameter_spec = parameter_spec
  )
  transform <- .iwmde_parameter_transform(support)
  xlim <- .iwmde_plot_range(posterior_values[finite_rows], support)
  if (!all(is.finite(xlim)) || xlim[1] >= xlim[2]) {
    plan[["status"]] <- "unsupported"
    plan[["reason"]] <- "could not construct a finite plotting range"
    plan[["rows"]]   <- .iwmde_plan_empty_rows()
    return(plan)
  }

  requested_values <- plan[["outputs"]][["requested_values"]]
  requested_values <- requested_values[
    requested_values >= support[1] &
      requested_values <= support[2]
  ]
  evaluation_values <- requested_values

  continuous_values <- posterior_values[continuous_rows]
  tail_probabilities <- NULL
  tail_values        <- NULL
  if (isTRUE(plan[["outputs"]][["need_density"]])) {
    tail_probabilities <- .iwmde_density_tail_probabilities()
    tail_values        <- .iwmde_density_tail_values(continuous_values)
  }
  density_xlim <- .iwmde_include_plot_values(
    xlim    = xlim,
    values  = requested_values,
    support = support
  )
  display_grid <- NULL
  if (isTRUE(plan[["outputs"]][["need_density"]])) {
    display_grid <- .iwmde_display_grid(
      xlim        = density_xlim,
      n_points    = plan[["control"]][["n_points"]],
      transform   = transform,
      values      = c(continuous_values, requested_values, tail_values),
      grid_method = plan[["control"]][["display_grid"]]
    )
    display_grid <- .iwmde_include_display_values(
      grid    = display_grid,
      values  = c(requested_values, tail_values),
      xlim    = density_xlim,
      support = support
    )
  }

  normalization_grid <- .iwmde_normalization_grid(
    values               = continuous_values,
    display_grid         = numeric(),
    support              = support,
    transform            = transform,
    normalization_points = plan[["control"]][["normalization_points"]],
    normalization_prob   = plan[["control"]][["normalization_prob"]]
  )

  plan[["status"]] <- "ok"
  plan[["rows"]] <- .iwmde_plan_rows(
    posterior_values  = posterior_values,
    finite_rows       = finite_rows,
    component         = component,
    continuous_rows   = continuous_rows,
    candidate_rows    = candidate_rows,
    candidate_values  = candidate_values,
    continuous_values = continuous_values,
    baseline_contract = baseline_contract,
    chain_coverage     = chain_coverage
  )
  plan[["support"]] <- list(
    x_lower   = support[1],
    x_upper   = support[2],
    support   = support,
    transform = transform,
    xlim      = xlim,
    density_xlim = density_xlim
  )
  plan[["replacement"]] <- .iwmde_replacement_spec(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec
  )
  plan[["grids"]] <- list(
    display_grid       = display_grid,
    requested_values   = requested_values,
    evaluation_values  = evaluation_values,
    tail_probabilities = tail_probabilities,
    tail_values        = tail_values,
    normalization_grid = normalization_grid
  )
  return(plan)
}


# Retain only fields used to execute or identify a plan.
.iwmde_plan_parameter_spec <- function(parameter_spec) {

  fields <- c(
    "type", "parameter", "weights", "conditional", "conditional_rule",
    "condition_key", "status", "reason"
  )

  return(parameter_spec[intersect(fields, names(parameter_spec))])
}


.iwmde_plan_baseline_contract <- function(context, plan, candidate_rows,
                                          candidate_values) {

  row_states <- .iwmde_row_states(
    context        = context,
    rows           = candidate_rows,
    parameter      = plan[["target"]][["parameter"]],
    parameter_spec = plan[["execution_spec"]],
    estimator      = plan[["method"]]
  )

  baseline_log_q <- tryCatch(
    vapply(row_states, function(state) {
      state[["baseline_log_q"]]
    }, numeric(1)),
    error = function(e) {
      if (inherits(e, "iwmde_construction_error")) {
        stop(e)
      }
      .iwmde_stop_construction_failure(
        estimator = plan[["method"]],
        parameter = plan[["target"]][["parameter"]],
        rows      = candidate_rows,
        stage     = "baseline joint-density evaluation",
        detail    = conditionMessage(e)
      )
    }
  )
  finite_baseline <- is.finite(baseline_log_q)
  if (!all(finite_baseline)) {
    .iwmde_stop_construction_failure(
      estimator = plan[["method"]],
      parameter = plan[["target"]][["parameter"]],
      rows      = candidate_rows[!finite_baseline],
      stage     = "baseline joint-density evaluation",
      detail    = "the baseline joint log density was not finite"
    )
  }
  .iwmde_validate_row_states(row_states)

  return(list(
    row_states       = row_states,
    baseline_log_q   = baseline_log_q,
    estimator_rows   = candidate_rows,
    estimator_values = candidate_values,
    baseline_rows_hash = .iwmde_hash(
      "iwmde_baseline_rows",
      candidate_rows
    )
  ))
}


.iwmde_plan_empty_rows <- function() {

  return(list(
    posterior_values  = numeric(),
    finite_rows       = logical(),
    samples           = numeric(),
    component         = NULL,
    continuous_rows   = integer(),
    candidate_rows    = integer(),
    candidate_values  = numeric(),
    continuous_values = numeric(),
    population_rows   = integer(),
    estimator_rows    = integer(),
    estimator_values  = numeric(),
    row_states        = list(),
    baseline_log_q    = numeric(),
    n_denominator_rows = 0L,
    n_estimator_rows  = 0L,
    baseline_rows_hash = .iwmde_hash(
      "iwmde_baseline_rows",
      integer()
    ),
    active_mass       = NA_real_,
    point_masses      = data.frame(x = numeric(), mass = numeric()),
    point_mass_total  = 0,
    n_candidate_rows  = 0L,
    n_total           = 0L,
    row_thinning_policy = list(
      method             = "nested_srswor",
      selected_rows_hash = .iwmde_hash("iwmde_rows", integer())
    )
  ))
}


.iwmde_plan_rows <- function(posterior_values, finite_rows, component,
                             continuous_rows = integer(),
                             candidate_rows = integer(),
                             candidate_values = numeric(),
                             continuous_values = numeric(),
                             baseline_contract = NULL,
                             chain_coverage = NULL) {

  point_masses <- component[["point_masses"]]
  if (is.null(point_masses)) {
    point_masses <- data.frame(x = numeric(), mass = numeric())
  }
  if (is.null(baseline_contract)) {
    baseline_contract <- list(
      row_states       = list(),
      baseline_log_q   = rep(NA_real_, length(candidate_rows)),
      estimator_rows   = candidate_rows,
      estimator_values = candidate_values,
      baseline_rows_hash = .iwmde_hash(
        "iwmde_baseline_rows",
        candidate_rows
      )
    )
  }

  return(list(
    posterior_values  = posterior_values,
    finite_rows       = finite_rows,
    samples           = posterior_values[finite_rows],
    component         = component,
    continuous_rows   = continuous_rows,
    candidate_rows    = candidate_rows,
    candidate_values  = candidate_values,
    continuous_values = continuous_values,
    population_rows   = which(finite_rows),
    estimator_rows    = baseline_contract[["estimator_rows"]],
    estimator_values  = baseline_contract[["estimator_values"]],
    row_states        = baseline_contract[["row_states"]],
    baseline_log_q    = baseline_contract[["baseline_log_q"]],
    n_denominator_rows = length(candidate_rows),
    n_estimator_rows  = length(baseline_contract[["estimator_rows"]]),
    baseline_rows_hash = baseline_contract[["baseline_rows_hash"]],
    chain_coverage     = chain_coverage,
    active_mass       = mean(component[["active"]][finite_rows]),
    point_masses      = point_masses,
    point_mass_total  = sum(point_masses[["mass"]]),
    n_candidate_rows  = length(candidate_rows),
    n_total           = sum(finite_rows),
    row_thinning_policy = list(
      method              = "nested_srswor",
      selected_rows_hash = .iwmde_hash("iwmde_rows", candidate_rows),
      n_continuous_rows  = length(continuous_rows),
      n_selected_rows    = length(candidate_rows),
      n_estimator_rows   = length(baseline_contract[["estimator_rows"]]),
      baseline_rows_hash = baseline_contract[["baseline_rows_hash"]]
    )
  ))
}


.iwmde_plan_chain_coverage <- function(context, eligible_rows,
                                       selected_rows = eligible_rows) {

  chain_id <- context[["chain_id"]]
  if (is.null(chain_id)) {
    chain_id <- rep(1L, nrow(context[["posterior_samples"]]))
  }
  if (length(chain_id) != nrow(context[["posterior_samples"]]) ||
      anyNA(chain_id)) {
    stop("Invalid fitted-chain metadata for IWMDE.", call. = FALSE)
  }
  expected_chain_ids <- unique(chain_id)
  chain_counts <- function(rows) {
    counts <- tabulate(
      match(chain_id[rows], expected_chain_ids),
      nbins = length(expected_chain_ids)
    )
    names(counts) <- as.character(expected_chain_ids)
    counts
  }

  return(list(
    expected_chain_ids    = expected_chain_ids,
    eligible_chain_counts = chain_counts(eligible_rows),
    selected_chain_counts = chain_counts(selected_rows)
  ))
}


.iwmde_plan_execution_spec <- function(parameter_spec) {

  keep <- c(
    "type",
    "parameter",
    "weights",
    "conditional",
    "conditional_rule",
    "condition_key"
  )
  parameter_spec <- parameter_spec[intersect(keep, names(parameter_spec))]

  return(.iwmde_compact_nulls(parameter_spec))
}


.iwmde_plan_target <- function(parameter, parameter_spec, metadata, target_key) {

  target <- list(
    type             = parameter_spec[["type"]],
    parameter        = parameter,
    target_key       = target_key,
    status           = parameter_spec[["status"]],
    reason           = parameter_spec[["reason"]],
    conditional      = parameter_spec[["conditional"]],
    conditional_rule = parameter_spec[["conditional_rule"]],
    condition_key    = parameter_spec[["condition_key"]],
    metadata         = .iwmde_target_provenance(metadata)
  )
  if (identical(parameter_spec[["type"]], "linear")) {
    target[["weights"]] <- .iwmde_plan_weights(parameter_spec[["weights"]])
  }

  return(.iwmde_compact_nulls(target))
}


.iwmde_plan_weights <- function(weights) {

  weights <- .iwmde_linear_weights(weights)
  if (is.null(weights) || length(weights) == 0L) {
    return(NULL)
  }
  weights <- weights[order(names(weights))]

  return(weights)
}


.iwmde_plan_key_payload <- function(plan) {

  return(.iwmde_compact_nulls(list(
    schema_version     = .iwmde_schema_version(),
    algorithm_version  = .iwmde_algorithm_version(),
    target             = plan[["target"]],
    outputs            = plan[["outputs"]],
    control            = .iwmde_density_control_provenance(plan[["control"]]),
    row_budget         = plan[["row_budget"]],
    method             = plan[["method"]],
    density_method     = plan[["density_method"]],
    status             = plan[["status"]],
    reason             = plan[["reason"]],
    rows               = .iwmde_plan_rows_provenance(plan),
    support            = .iwmde_plan_support_provenance(plan),
    grids              = .iwmde_plan_grid_provenance(plan),
    prior_ordinates    = plan[["prior_ordinates"]],
    ordinate_warnings  = plan[["ordinate_warnings"]],
    replacement        = .iwmde_hash(
      "iwmde_replacement",
      plan[["replacement"]]
    ),
    source_fingerprint = plan[["source_fingerprint"]]
  )))
}


.iwmde_source_fingerprint <- function(context) {

  if (!is.null(context[["source_fingerprint"]])) {
    return(context[["source_fingerprint"]])
  }

  return(.iwmde_compute_source_fingerprint(context))
}


.iwmde_compute_source_fingerprint <- function(context) {

  samples <- context[["posterior_samples"]]

  return(.iwmde_compact_nulls(list(
    object_class    = class(context[["object"]]),
    fit_class       = class(context[["object"]][["fit"]]),
    posterior_dim   = dim(samples),
    posterior_names = .iwmde_hash("iwmde_columns", colnames(samples)),
    posterior_values = .iwmde_hash("iwmde_posterior_draws", samples),
    chain_id          = .iwmde_hash("iwmde_chain_id", context[["chain_id"]]),
    data             = .iwmde_hash("iwmde_data", context[["data"]]),
    priors           = .iwmde_hash("iwmde_priors", context[["priors"]]),
    selection_spec   = .iwmde_hash(
      "iwmde_selection",
      context[["selection_spec"]]
    ),
    likelihood_family = context[["object"]][["likelihood"]][["family"]],
    measure           = context[["object"]][["data"]][["measure"]]
  )))
}


.iwmde_plan_provenance <- function(plan) {

  return(.iwmde_compact_nulls(list(
    schema_version      = .iwmde_schema_version(),
    algorithm_version   = .iwmde_algorithm_version(),
    provenance_level    = "iwmde_plan",
    plan_key            = plan[["plan_key"]],
    density_method      = plan[["density_method"]],
    method              = plan[["method"]],
    internal_method     = plan[["method"]],
    density_control     = .iwmde_density_control_provenance(plan[["control"]]),
    row_budget          = plan[["row_budget"]],
    target              = plan[["target"]],
    status              = plan[["status"]],
    reason              = plan[["reason"]],
    rows                = .iwmde_plan_rows_provenance(plan),
    support             = .iwmde_plan_support_provenance(plan),
    grids               = .iwmde_plan_grid_provenance(plan),
    prior_ordinates     = plan[["prior_ordinates"]],
    ordinate_warnings   = plan[["ordinate_warnings"]],
    source_fingerprint  = plan[["source_fingerprint"]]
  )))
}


.iwmde_plan_rows_provenance <- function(plan) {

  rows <- plan[["rows"]]
  if (is.null(rows)) {
    return(NULL)
  }

  return(.iwmde_compact_nulls(list(
    n_total              = rows[["n_total"]],
    n_candidate_rows     = rows[["n_candidate_rows"]],
    n_denominator_rows   = rows[["n_denominator_rows"]],
    n_estimator_rows     = rows[["n_estimator_rows"]],
    n_continuous_rows    = length(rows[["continuous_rows"]]),
    active_mass          = rows[["active_mass"]],
    point_mass_total     = rows[["point_mass_total"]],
    chain_coverage       = rows[["chain_coverage"]],
    selected_rows_hash   =
      rows[["row_thinning_policy"]][["selected_rows_hash"]],
    baseline_rows_hash   = rows[["baseline_rows_hash"]],
    row_thinning_policy  = rows[["row_thinning_policy"]]
  )))
}


.iwmde_plan_support_provenance <- function(plan) {

  support <- plan[["support"]]
  if (is.null(support)) {
    return(NULL)
  }

  return(.iwmde_compact_nulls(list(
    x_lower   = support[["x_lower"]],
    x_upper   = support[["x_upper"]],
    transform = support[["transform"]],
    xlim      = support[["xlim"]],
    density_xlim = support[["density_xlim"]]
  )))
}


.iwmde_plan_grid_provenance <- function(plan) {

  grids <- plan[["grids"]]
  if (is.null(grids)) {
    return(NULL)
  }

  normalization_grid <- grids[["normalization_grid"]]
  return(.iwmde_compact_nulls(list(
    display_grid_hash = .iwmde_hash(
      "iwmde_display_grid",
      grids[["display_grid"]]
    ),
    requested_values = grids[["requested_values"]],
    evaluation_values = grids[["evaluation_values"]],
    tail_probabilities = grids[["tail_probabilities"]],
    tail_values = grids[["tail_values"]],
    normalization_grid = if (is.null(normalization_grid)) {
      NULL
    } else {
      list(
        n_points = length(normalization_grid[["x"]]),
        range    = range(normalization_grid[["x"]]),
        z_range  = range(normalization_grid[["z"]])
      )
    }
  )))
}


.iwmde_request_provenance <- function(context, parameter, density_method,
                                      density_control,
                                      attribute = c("density", "ordinate"),
                                      value = NULL, parameter_spec = NULL,
                                      metadata = NULL) {

  context        <- .iwmde_context_ensure_caches(context)
  attribute      <- match.arg(attribute)
  density_method <- .density_method_normalize_precomputed(density_method)
  density_control <- .iwmde_density_control_resolve(
    density_method  = density_method,
    density_control = density_control,
    purpose         = attribute
  )
  if (is.null(density_control[["normalization_points"]])) {
    density_control[["normalization_points"]] <- max(
      50L,
      density_control[["n_points"]]
    )
  }

  if (identical(attribute, "ordinate")) {
    value <- suppressWarnings(as.numeric(value))
    value <- value[is.finite(value)]
    if (length(value) > 1L) {
      value <- value[[1L]]
    }
  } else {
    value <- NULL
  }
  ordinate_values <- if (identical(attribute, "ordinate")) value else NULL
  parameter_spec <- .iwmde_prepare_prior_ordinates(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    values         = ordinate_values
  )
  prior_ordinates <- if (identical(attribute, "ordinate")) {
    parameter_spec[["prior_ordinates"]]
  } else {
    NULL
  }
  parameter_spec[["prior_ordinates"]] <- NULL
  parameter_spec <- .iwmde_parameter_spec(context, parameter, parameter_spec)

  return(.iwmde_provenance_request(
    density_method     = density_method,
    method             = .density_method_iwmde_estimator(density_method),
    metadata           = metadata,
    density_control    = density_control,
    value              = value,
    attribute          = attribute,
    target_key         = .iwmde_target_key(parameter, parameter_spec),
    source_fingerprint = .iwmde_source_fingerprint(context),
    prior_ordinates    = prior_ordinates
  ))
}
