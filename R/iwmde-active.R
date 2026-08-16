# ============================================================================ #
# IWMDE Active Priors, Rows, and Formula Inputs
# ============================================================================ #

.iwmde_active_setup <- function(context, row, active_key = NULL) {

  key <- if (is.null(active_key)) .iwmde_active_key(context, row) else active_key
  if (exists(key, envir = context[["active_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["active_cache"]]))
  }

  priors <- .iwmde_active_nested_priors(context, row)
  fit_priors <- .create_fit_priors(
    data   = context[["data"]],
    priors = priors
  )
  fit_data <- .create_fit_data(data = context[["data"]], priors = priors)
  fit_data <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = priors,
    effect_direction = .data_effect_direction(context[["data"]])
  )

  setup <- list(
    priors            = priors,
    fit_priors        = fit_priors,
    fit_data          = fit_data,
    selection_spec    = .iwmde_selection_spec(context[["data"]], priors),
    is_PET            = .is_priors_PET(priors),
    is_PEESE          = .is_priors_PEESE(priors),
    is_weightfunction = .is_priors_weightfunction(priors)
  )

  assign(key, setup, envir = context[["active_cache"]])
  return(setup)
}


.iwmde_active_key <- function(context, row) {

  indicator_names <- context[["indicator_names"]]
  if (length(indicator_names) == 0L) {
    return("all")
  }

  values <- vapply(indicator_names, function(name) {
    if (name %in% names(row)) {
      as.character(.iwmde_indicator_index(row[[name]], name))
    } else {
      "NA"
    }
  }, character(1))

  return(paste(paste(indicator_names, values, sep = "="), collapse = "|"))
}


.iwmde_active_nested_priors <- function(context, row) {

  priors <- context[["priors"]]

  if (!is.null(priors[["outcome"]])) {
    for (name in names(priors[["outcome"]])) {
      priors[["outcome"]][[name]] <- .iwmde_select_prior_component(
        prior     = priors[["outcome"]][[name]],
        row       = row,
        indicator = .iwmde_nested_indicator_name("outcome", name)
      )
    }
  }
  if (!is.null(priors[["mods"]])) {
    for (name in names(priors[["mods"]])) {
      priors[["mods"]][[name]] <- .iwmde_select_prior_component(
        prior     = priors[["mods"]][[name]],
        row       = row,
        indicator = paste0(BayesTools::JAGS_parameter_names(name, "mu"), "_indicator")
      )
    }
  }
  if (!is.null(priors[["scale"]])) {
    for (name in names(priors[["scale"]])) {
      priors[["scale"]][[name]] <- .iwmde_select_prior_component(
        prior     = priors[["scale"]][[name]],
        row       = row,
        indicator = paste0(BayesTools::JAGS_parameter_names(name, "log_tau"), "_indicator")
      )
    }
  }

  return(priors)
}


.iwmde_nested_indicator_name <- function(section, name) {

  if (section == "outcome" && name == "bias") {
    return("bias_indicator")
  }

  return(paste0(name, "_indicator"))
}


.iwmde_select_prior_component <- function(prior, row, indicator) {

  if (!BayesTools::is.prior.mixture(prior) || !indicator %in% names(row)) {
    return(prior)
  }

  index <- .iwmde_indicator_index(row[[indicator]], indicator, length(prior))

  return(prior[[index]])
}


.iwmde_indicator_index <- function(value, name, max_index = NULL) {

  index <- .as_exact_model_indicator(value, name, scalar = TRUE)
  if (index < 1L || (!is.null(max_index) && index > max_index)) {
    stop("'", name, "' is outside the available prior components.", call. = FALSE)
  }

  return(index)
}


.iwmde_active_flat_prior_list <- function(context, row, parameter = NULL,
                                          state_scope = c("local", "global")) {

  state_scope <- match.arg(state_scope)
  prior_list <- context[["flat_prior_list"]]
  if (length(prior_list) == 0L) {
    return(list())
  }

  out <- list()
  for (name in names(prior_list)) {
    if (identical(state_scope, "global") &&
        .iwmde_prior_name_is_local_latent(name)) {
      next
    }

    prior <- .iwmde_select_prior_component(
      prior     = prior_list[[name]],
      row       = row,
      indicator = .iwmde_indicator_name(name)
    )
    if (!BayesTools::is.prior.none(prior)) {
      if (BayesTools::is.prior.weightfunction(prior)) {
        if (.iwmde_parameter_is_eta(parameter) &&
            length(BayesTools::JAGS_indexed_parameter_vector(row, "eta")) > 0L) {
          out[[name]] <- prior
        }
      } else {
        out[[name]] <- prior
      }
    }
  }

  return(out)
}


.iwmde_selection_spec <- function(data, priors) {

  if (.data_outcome_type(data) != "norm" || !.is_priors_weightfunction(priors)) {
    return(NULL)
  }

  yi <- data[["outcome"]][["yi"]]
  if (.data_effect_direction(data) == "negative") {
    yi <- -yi
  }

  return(.selection_spec(
    priors           = priors,
    yi               = yi,
    sei              = data[["outcome"]][["sei"]],
    effect_direction = .data_effect_direction(data),
    signed_data      = TRUE
  ))
}


.iwmde_row_parameters <- function(context, row, active_setup,
                                  state_scope = c("local", "global")) {

  state_scope <- match.arg(state_scope)
  row <- .resolve_fixed_prior_row(row, active_setup[["fit_priors"]])
  parameters <- .iwmde_marglik_row_parameters(
    row        = row,
    fit_priors = active_setup[["fit_priors"]]
  )

  data <- context[["data"]]
  K    <- nrow(data[["outcome"]])

  if (.is_data_mods(data)) {
    parameters[["mu"]] <- .iwmde_evaluate_formula_row(
      context   = context,
      row       = row,
      source    = "mods",
      parameter = "mu"
    )
  } else if ("mu" %in% names(row)) {
    parameters[["mu"]] <- row[["mu"]]
  } else if (is.null(parameters[["mu"]])) {
    parameters[["mu"]] <- 0
  }

  if (.is_data_scale(data)) {
    parameters[["log_tau"]] <- .iwmde_evaluate_formula_row(
      context   = context,
      row       = row,
      source    = "scale",
      parameter = "log_tau"
    )
  } else if ("tau" %in% names(row)) {
    parameters[["tau"]] <- row[["tau"]]
  } else if (is.null(parameters[["tau"]])) {
    fixed_tau <- .fixed_tau_prior_value(active_setup[["priors"]])
    parameters[["tau"]] <- if (is.null(fixed_tau)) 0 else fixed_tau
  }

  if (.is_data_multilevel(data)) {
    if ("rho" %in% names(row)) {
      parameters[["rho"]] <- row[["rho"]]
    }
    if (identical(state_scope, "local")) {
      gamma_cols <- paste0(
        "gamma[",
        seq_len(length(unique(data[["outcome"]][["cluster"]]))),
        "]"
      )
      if (all(gamma_cols %in% names(row))) {
        parameters[["gamma"]] <- as.numeric(row[gamma_cols])
      }
    }
  }

  if (active_setup[["is_PET"]]) {
    if ("PET" %in% names(row)) {
      parameters[["PET"]] <- row[["PET"]]
    } else if (is.null(parameters[["PET"]])) {
      fixed_PET <- .fixed_bias_parameter_value(active_setup[["priors"]], "PET")
      parameters[["PET"]] <- if (is.null(fixed_PET)) 0 else fixed_PET
    }
  }
  if (active_setup[["is_PEESE"]]) {
    if ("PEESE" %in% names(row)) {
      parameters[["PEESE"]] <- row[["PEESE"]]
    } else if (is.null(parameters[["PEESE"]])) {
      fixed_PEESE <- .fixed_bias_parameter_value(active_setup[["priors"]], "PEESE")
      parameters[["PEESE"]] <- if (is.null(fixed_PEESE)) 0 else fixed_PEESE
    }
  }
  if (active_setup[["is_weightfunction"]] &&
      !.iwmde_omega_is_usable(parameters[["omega"]], active_setup)) {
    parameters[["omega"]] <- .iwmde_active_omega(context, row, active_setup)
  }

  if (identical(state_scope, "local") &&
      .data_outcome_type(data) == "bin") {
    pi_cols    <- paste0("pi[", seq_len(K), "]")
    theta_cols <- paste0("theta[", seq_len(K), "]")
    if (all(pi_cols %in% names(row))) {
      parameters[["pi"]] <- as.numeric(row[pi_cols])
    }
    if (all(theta_cols %in% names(row))) {
      parameters[["theta"]] <- as.numeric(row[theta_cols])
    }
  } else if (identical(state_scope, "local") &&
             .data_outcome_type(data) == "pois") {
    phi_cols   <- paste0("phi[", seq_len(K), "]")
    theta_cols <- paste0("theta[", seq_len(K), "]")
    if (all(phi_cols %in% names(row))) {
      parameters[["phi"]] <- as.numeric(row[phi_cols])
    }
    if (all(theta_cols %in% names(row))) {
      parameters[["theta"]] <- as.numeric(row[theta_cols])
    }
  }

  if (identical(state_scope, "global")) {
    parameters[c("gamma", "theta", "pi", "phi")] <- NULL
  }

  parameters <- .iwmde_attach_marginalized_random_parameters(
    parameters = parameters,
    data       = data,
    row        = row
  )

  return(parameters)
}


.iwmde_marglik_row_parameters <- function(row, fit_priors) {

  if (length(fit_priors) == 0L) {
    return(list())
  }

  parameters <- list()
  for (i in seq_along(fit_priors)) {
    if (.iwmde_direct_marglik_parameter_missing(
      row   = row,
      prior = fit_priors[[i]],
      name  = names(fit_priors)[[i]]
    )) {
      next
    }
    current <- tryCatch(
      BayesTools::JAGS_marglik_parameters(row, fit_priors[i]),
      error = function(e) {

        if (.iwmde_marglik_parameters_missing(e)) {
          return(list())
        }
        stop(e)
      }
    )
    parameters <- c(parameters, current)
  }

  return(parameters)
}


.iwmde_direct_marglik_parameter_missing <- function(row, prior, name) {

  if (!BayesTools::is.prior.simple(prior) ||
      BayesTools::is.prior.point(prior)) {
    return(FALSE)
  }

  return(!name %in% names(row))
}


.iwmde_marglik_parameters_missing <- function(error) {

  message <- conditionMessage(error)
  return(startsWith(message, "'samples' does not contain all monitored "))
}


.iwmde_attach_marginalized_random_parameters <- function(parameters, data,
                                                         row) {

  if (!.data_has_marginalized_random_effects(data)) {
    return(parameters)
  }

  terms <- .data_marginalized_random_effects(data)
  for (term in terms) {
    parameter <- .iwmde_marginalized_random_parameter_columns(term)
    for (name in parameter) {
      if (name %in% names(row) && is.null(parameters[[name]])) {
        parameters[[name]] <- row[[name]]
      }
    }
  }

  return(parameters)
}


.iwmde_marginalized_random_parameter_columns <- function(term) {

  columns <- term[["sd_parameter_names"]]

  sources <- .predict_known_v_marginalized_sd_sources(term)
  for (source in sources) {
    if (identical(source[["shape"]], "scalar")) {
      columns <- c(columns, source[["name"]])
    }
  }

  if (.marginalized_random_effect_has_allocation(term)) {
    factors <- .marginalized_random_effect_allocation_factors(term)
    columns <- c(columns, vapply(factors, function(factor) {
      paste0(factor[["weight_name"]], "[", factor[["index"]], "]")
    }, character(1)))
  }

  unique(columns[!is.na(columns) & nzchar(columns)])
}


.iwmde_omega_is_usable <- function(omega, active_setup) {

  if (is.null(omega) || is.null(active_setup[["selection_spec"]])) {
    return(FALSE)
  }
  if (length(omega) != active_setup[["selection_spec"]][["n_bins"]]) {
    return(FALSE)
  }

  return(all(is.finite(omega)) && all(omega >= 0))
}


.iwmde_active_omega <- function(context, row, active_setup) {

  active_spec <- active_setup[["selection_spec"]]
  if (is.null(active_spec)) {
    return(NULL)
  }

  omega <- as.numeric(BayesTools::JAGS_indexed_parameter_vector(
    row       = row,
    parameter = active_spec[["jags_omega"]]
  ))
  if (length(omega) == 0L) {
    fixed_omega <- .iwmde_fixed_weightfunction_omega(
      active_setup[["priors"]][["outcome"]][["bias"]]
    )
    if (length(fixed_omega) > 0L) {
      return(.iwmde_validate_omega(fixed_omega, active_spec[["n_bins"]]))
    }
  }
  if (length(omega) == 0L) {
    jags_data <- active_spec[["jags_data"]]
    if (is.list(jags_data)) {
      fixed_omega <- jags_data[[active_spec[["jags_omega"]]]]
      if (!is.null(fixed_omega)) {
        return(.iwmde_validate_omega(fixed_omega, active_spec[["n_bins"]]))
      }
    }

    return(rep(NA_real_, active_spec[["n_bins"]]))
  }

  global_spec <- context[["selection_spec"]]
  if (length(omega) == active_spec[["n_bins"]]) {
    return(.iwmde_validate_omega(omega, active_spec[["n_bins"]]))
  }

  if (is.null(global_spec) || length(omega) != global_spec[["n_bins"]]) {
    return(rep(NA_real_, active_spec[["n_bins"]]))
  }

  omega <- .iwmde_collapse_omega(
    omega        = omega,
    global_cuts  = global_spec[["p_cuts"]],
    active_cuts  = active_spec[["p_cuts"]]
  )

  return(.iwmde_validate_omega(omega, active_spec[["n_bins"]]))
}


.iwmde_same_p_cuts <- function(x, y) {

  if (length(x) != length(y)) {
    return(FALSE)
  }

  return(identical(as.numeric(x), as.numeric(y)))
}


.iwmde_fixed_weightfunction_omega <- function(prior) {

  if (is.null(prior) || !BayesTools::is.prior.weightfunction(prior)) {
    return(numeric())
  }
  weights <- prior[["weights"]]
  if (is.null(weights) || !identical(weights[["type"]], "fixed")) {
    return(numeric())
  }

  return(as.numeric(weights[["omega"]]))
}


.iwmde_collapse_omega <- function(omega, global_cuts, active_cuts) {

  bins <- .iwmde_collapse_omega_bins(global_cuts, active_cuts)
  out  <- rep(NA_real_, length(bins))
  for (i in seq_along(bins)) {
    index <- bins[[i]]
    if (length(index) == 0L || any(index > length(omega))) {
      next
    }
    if (length(index) == 1L) {
      out[i] <- omega[index]
      next
    }
    values <- as.numeric(omega[index])
    if (all(is.finite(values)) &&
        all(values == values[1L])) {
      out[i] <- values[1L]
    }
  }

  return(out)
}


.iwmde_collapse_omega_matrix <- function(omega, global_cuts, active_cuts) {

  bins <- .iwmde_collapse_omega_bins(global_cuts, active_cuts)
  out  <- matrix(NA_real_, nrow = nrow(omega), ncol = length(bins))
  for (i in seq_along(bins)) {
    index <- bins[[i]]
    if (length(index) == 0L || any(index > ncol(omega))) {
      next
    }
    if (length(index) == 1L) {
      out[, i] <- omega[, index]
      next
    }
    values <- omega[, index, drop = FALSE]
    finite <- rowSums(is.finite(values)) == ncol(values)
    keep <- finite
    for (j in seq.int(2L, ncol(values))) {
      keep <- keep & values[, j] == values[, 1L]
    }
    if (any(keep)) {
      out[keep, i] <- values[keep, 1L]
    }
  }

  return(out)
}


.iwmde_collapse_omega_bins <- function(global_cuts, active_cuts) {

  out <- vector("list", length(active_cuts) - 1L)

  for (i in seq_along(out)) {
    lower <- active_cuts[i]
    upper <- active_cuts[i + 1L]
    bins  <- which(
      global_cuts[-length(global_cuts)] >= lower &
        global_cuts[-1L] <= upper
    )
    if (length(bins) == 0L) {
      out[[i]] <- integer()
    } else {
      out[[i]] <- bins
    }
  }

  return(out)
}


.iwmde_validate_omega <- function(omega, n_bins) {

  if (length(omega) != n_bins || any(!is.finite(omega)) || any(omega < 0)) {
    return(rep(NA_real_, n_bins))
  }

  return(omega)
}


.iwmde_formula_inputs <- function(data, priors) {

  out <- list()
  if (.is_data_mods(data)) {
    out[["mods"]] <- .iwmde_formula_input(
      data      = data,
      priors    = priors,
      source    = "mods",
      parameter = "mu"
    )
  }
  if (.is_data_scale(data)) {
    out[["scale"]] <- .iwmde_formula_input(
      data      = data,
      priors    = priors,
      source    = "scale",
      parameter = "log_tau"
    )
  }

  return(out)
}


.iwmde_formula_input <- function(data, priors, source, parameter) {

  return(list(
    formula    = .create_fit_formula_list(data = data, parameter = source),
    data       = .create_fit_formula_data_list(data = data, parameter = source),
    prior_list = .repair_formula_prior_list(
      prior_list = .create_fit_formula_prior_list(
        priors    = priors,
        parameter = source
      ),
      parameter = parameter
    )
  ))
}


.iwmde_evaluate_formula_row <- function(context, row, source, parameter) {

  input <- context[["formula_inputs"]][[source]]

  return(.iwmde_evaluate_formula_row_input(
    context   = context,
    row       = row,
    input     = input,
    parameter = parameter
  ))
}


.iwmde_evaluate_formula_row_input <- function(context, row, input,
                                              parameter) {

  posterior_row <- matrix(row, nrow = 1)
  colnames(posterior_row) <- names(row)

  values <- BayesTools::JAGS_evaluate_formula(
    fit            = .posterior_formula_fit(context[["formula_fit"]], posterior_row),
    formula        = input[["formula"]],
    parameter      = parameter,
    data           = input[["data"]],
    prior_list     = input[["prior_list"]],
    formula_target = "fixed"
  )

  return(as.numeric(values[, 1]))
}
