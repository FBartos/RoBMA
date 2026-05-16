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

  if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
    stop("'", name, "' must be a finite scalar model indicator.", call. = FALSE)
  }

  rounded <- round(value)
  if (abs(value - rounded) > 1e-8) {
    stop("'", name, "' must be integer-valued.", call. = FALSE)
  }

  index <- as.integer(rounded)
  if (index < 1L || (!is.null(max_index) && index > max_index)) {
    stop("'", name, "' is outside the available prior components.", call. = FALSE)
  }

  return(index)
}


.iwmde_active_flat_prior_list <- function(context, row, parameter = NULL) {

  prior_list <- context[["flat_prior_list"]]
  if (length(prior_list) == 0L) {
    return(list())
  }

  out <- list()
  for (name in names(prior_list)) {
    prior <- .iwmde_select_prior_component(
      prior     = prior_list[[name]],
      row       = row,
      indicator = .iwmde_indicator_name(name)
    )
    if (!BayesTools::is.prior.none(prior)) {
      if (BayesTools::is.prior.weightfunction(prior)) {
        if (.iwmde_parameter_is_eta(parameter) &&
            length(.iwmde_row_indexed_parameter(row, "eta")) > 0L) {
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


.iwmde_row_parameters <- function(context, row, active_setup) {

  parameters <- list()
  base_parameters <- try(
    BayesTools::JAGS_marglik_parameters(row, active_setup[["fit_priors"]]),
    silent = TRUE
  )
  if (!inherits(base_parameters, "try-error")) {
    parameters <- c(parameters, base_parameters)
  }

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
    parameters[["tau"]] <- 0
  }

  if (.is_data_multilevel(data)) {
    if ("rho" %in% names(row)) {
      parameters[["rho"]] <- row[["rho"]]
    }
    gamma_cols <- paste0("gamma[", seq_len(length(unique(data[["outcome"]][["cluster"]]))), "]")
    if (all(gamma_cols %in% names(row))) {
      parameters[["gamma"]] <- as.numeric(row[gamma_cols])
    }
  }

  if (active_setup[["is_PET"]]) {
    parameters[["PET"]] <- if ("PET" %in% names(row)) row[["PET"]] else 0
  }
  if (active_setup[["is_PEESE"]]) {
    parameters[["PEESE"]] <- if ("PEESE" %in% names(row)) row[["PEESE"]] else 0
  }
  if (active_setup[["is_weightfunction"]] &&
      !.iwmde_omega_is_usable(parameters[["omega"]], active_setup)) {
    parameters[["omega"]] <- .iwmde_active_omega(context, row, active_setup)
  }

  if (.data_outcome_type(data) == "bin") {
    pi_cols    <- paste0("pi[", seq_len(K), "]")
    theta_cols <- paste0("theta[", seq_len(K), "]")
    if (all(pi_cols %in% names(row))) {
      parameters[["pi"]] <- as.numeric(row[pi_cols])
    }
    if (all(theta_cols %in% names(row))) {
      parameters[["theta"]] <- as.numeric(row[theta_cols])
    }
  } else if (.data_outcome_type(data) == "pois") {
    phi_cols   <- paste0("phi[", seq_len(K), "]")
    theta_cols <- paste0("theta[", seq_len(K), "]")
    if (all(phi_cols %in% names(row))) {
      parameters[["phi"]] <- as.numeric(row[phi_cols])
    }
    if (all(theta_cols %in% names(row))) {
      parameters[["theta"]] <- as.numeric(row[theta_cols])
    }
  }

  return(parameters)
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

  omega <- .iwmde_row_indexed_parameter(row, "omega")
  if (length(omega) == 0L) {
    log_omega <- .iwmde_row_indexed_parameter(row, "log_omega")
    if (length(log_omega) > 0L) {
      omega <- exp(log_omega)
    }
  }
  if (length(omega) == 0L) {
    fixed_omega <- .iwmde_fixed_weightfunction_omega(
      active_setup[["priors"]][["outcome"]][["bias"]]
    )
    if (length(fixed_omega) > 0L) {
      return(.iwmde_validate_omega(fixed_omega, active_spec[["n_bins"]]))
    }
  }
  if (length(omega) == 0L) {
    return(.iwmde_validate_omega(rep(1, active_spec[["n_bins"]]), active_spec[["n_bins"]]))
  }

  global_spec <- context[["selection_spec"]]
  if (length(omega) == active_spec[["n_bins"]] &&
      (is.null(global_spec) ||
        .iwmde_same_p_cuts(global_spec[["p_cuts"]], active_spec[["p_cuts"]]))) {
    return(.iwmde_validate_omega(omega, active_spec[["n_bins"]]))
  }

  if (is.null(global_spec) || length(omega) != global_spec[["n_bins"]]) {
    return(.iwmde_validate_omega(
      omega[seq_len(active_spec[["n_bins"]])],
      active_spec[["n_bins"]]
    ))
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

  return(all(abs(as.numeric(x) - as.numeric(y)) <= 1e-12))
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


.iwmde_row_indexed_parameter <- function(row, parameter) {

  cols <- grep(paste0("^", parameter, "\\["), names(row), value = TRUE)
  if (length(cols) == 0L) {
    return(numeric())
  }

  index <- as.integer(sub(
    paste0("^", parameter, "\\[([0-9]+)\\]$"),
    "\\1",
    cols
  ))
  cols <- cols[order(index)]

  return(as.numeric(row[cols]))
}


.iwmde_collapse_omega <- function(omega, global_cuts, active_cuts) {

  index <- .iwmde_collapse_omega_index(global_cuts, active_cuts)
  out   <- rep(NA_real_, length(index))
  keep  <- !is.na(index)
  out[keep] <- omega[index[keep]]

  return(out)
}


.iwmde_collapse_omega_matrix <- function(omega, global_cuts, active_cuts) {

  index <- .iwmde_collapse_omega_index(global_cuts, active_cuts)
  out   <- matrix(NA_real_, nrow = nrow(omega), ncol = length(index))
  keep  <- !is.na(index)
  if (any(keep)) {
    out[, keep] <- omega[, index[keep], drop = FALSE]
  }

  return(out)
}


.iwmde_collapse_omega_index <- function(global_cuts, active_cuts) {

  out <- rep(NA_integer_, length(active_cuts) - 1L)

  for (i in seq_along(out)) {
    lower <- active_cuts[i]
    upper <- active_cuts[i + 1L]
    bins  <- which(
      global_cuts[-length(global_cuts)] >= lower - 1e-12 &
        global_cuts[-1L] <= upper + 1e-12
    )
    if (length(bins) == 0L) {
      out[i] <- NA_integer_
    } else {
      out[i] <- bins[1L]
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
    fit        = .posterior_formula_fit(context[["formula_fit"]], posterior_row),
    formula    = input[["formula"]],
    parameter  = parameter,
    data       = input[["data"]],
    prior_list = input[["prior_list"]]
  )

  return(as.numeric(values[, 1]))
}
