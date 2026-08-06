# Fitting functions -----
.create_fit_priors <- function(data, priors) {

  .check_glmm_no_bias_priors(data, priors)

  # extract the common prior list
  prior_list <- priors[["outcome"]]

  # add levels to the baseline prior for outcomes
  if (.data_outcome_type(data) == "bin") {
    # encode number of levels for the baserate and random-effects prior
    attr(prior_list[["pi"]], "levels")    <- nrow(data[["outcome"]])
    attr(prior_list[["theta"]], "levels") <- nrow(data[["outcome"]])
  } else if (.data_outcome_type(data) == "pois") {
    # encode number of levels for the baserate and random-effects prior
    attr(prior_list[["phi"]], "levels")   <- nrow(data[["outcome"]])
    attr(prior_list[["theta"]], "levels") <- nrow(data[["outcome"]])
  }

  # add cluster-level indicators
  if (.is_data_multilevel(data)) {
    # encode number of levels for the random-effects prior
    attr(prior_list[["gamma"]], "levels") <- length(unique(data[["outcome"]][["cluster"]]))
  }
  # add known-V sampling dependency latent factors
  if (.is_data_known_v_backend(data, "latent") && .data_known_v_rank(data) > 0L) {
    prior_list[["sampling_z"]] <- BayesTools::prior_factor(
      "normal",
      parameters = list("mean" = 0, "sd" = 1),
      contrast   = "independent"
    )
    attr(prior_list[["sampling_z"]], "levels") <- .data_known_v_rank(data)
  }

  ### deal with non-prior mixture distributions (bPET, bPEESE, and bselmodel)
  # the prior name must match the output parameter for sample extraction in plots etc
  if (.is_priors_bias(priors)) {
    if (!BayesTools::is.prior.mixture(prior_list[["bias"]])) {
      if (BayesTools::is.prior.weightfunction(prior_list[["bias"]])) {
        names(prior_list)[names(prior_list) == "bias"] <- "omega"
      } else if (.prior_has_phacking(prior_list[["bias"]])) {
        .selection_stop_phacking_deferred()
      } else if (BayesTools::is.prior.PET(prior_list[["bias"]])) {
        names(prior_list)[names(prior_list) == "bias"] <- "PET"
      } else if (BayesTools::is.prior.PEESE(prior_list[["bias"]])) {
        names(prior_list)[names(prior_list) == "bias"] <- "PEESE"
      }
    }
  }


  return(prior_list)
}
.create_fit_data   <- function(data, priors) {

  .check_glmm_no_bias_priors(data, priors)

  ### add outcome specific data
  if (.data_outcome_type(data) == "norm") {
    # flip effect size direction (needed for selection models, PET, and PEESE models)
    # (done for everything for consistency)
    yi                <- data[["outcome"]][["yi"]]
    effect_direction  <- .data_effect_direction(data)
    if (effect_direction == "negative") {
      yi <- -1 * yi
    }

    known_v_whitened  <- .is_data_known_v_backend(data, "whitened")
    known_v_block_mvn <- .is_data_known_v_backend(data, "block_mvn")
    if (known_v_whitened) {
      known_V          <- .data_known_v_data(data)
      whitening_blocks <- .known_v_backend_blocks(known_V, "whitened")
      independent <- .known_v_independent_indices(known_V)
      fit_data <- list()
      if (length(independent) > 0L) {
        fit_data[["known_v_independent_n"]]     <- length(independent)
        fit_data[["known_v_independent_index"]] <- independent
        fit_data[["known_v_independent_y"]]     <- yi[independent]
        fit_data[["known_v_independent_var"]]   <-
          .known_v_diagonal(known_V)[independent]
      }
      for (b in seq_along(whitening_blocks)) {
        block <- whitening_blocks[[b]]
        index <- block[["index"]]
        fit_data[[paste0("whitening_y_", b)]] <- as.vector(
          block[["rotation"]] %*% yi[index]
        )
        fit_data[[paste0("whitening_var_", b)]] <- block[["variance"]]
        fit_data[[paste0("whitening_matrix_", b)]] <- block[["rotation"]]
      }
    } else if (known_v_block_mvn) {
      known_V <- .data_known_v_data(data)
      fit_data <- list()
      block_mvn_blocks <- .known_v_backend_blocks(known_V, "block_mvn")

      independent <- .known_v_independent_indices(known_V)
      if (length(independent) > 0L) {
        fit_data[["known_v_independent_n"]]     <- length(independent)
        fit_data[["known_v_independent_index"]] <- independent
        fit_data[["known_v_independent_y"]]     <- yi[independent]
        fit_data[["known_v_independent_var"]]   <- .known_v_diagonal(known_V)[independent]
      }

      for (b in seq_along(block_mvn_blocks)) {
        block <- block_mvn_blocks[[b]]
        if (block[["size"]] > 1L) {
          fit_data[[paste0("known_v_y_", b)]]     <- yi[block[["index"]]]
          fit_data[[paste0("known_v_lower_", b)]] <- block[["v_lower"]]
        }
      }
    } else {
      # ordinary normal and latent known-V bias graphs reference sei
      fit_data <- list(
        yi = yi
      )
      if (!.is_data_known_v(data) || .is_priors_bias(priors)) {
        fit_data[["sei"]] <- data[["outcome"]][["sei"]]
      }
    }

    # add selection-kernel data for selection models
    if (.is_priors_weightfunction(priors)) {
      selection_spec <- .selection_spec(
        priors           = priors,
        yi               = yi,
        sei              = fit_data[["sei"]],
        effect_direction = effect_direction,
        signed_data      = TRUE
      )

      fit_data <- c(fit_data, selection_spec[["jags_data"]])

    }

    if (.is_data_known_v_backend(data, c("latent", "diagonal"))) {
      known_V <- .data_known_v_data(data)
      fit_data[["sampling_var"]] <- .known_v_residual_variance(known_V)
      if (.known_v_rank(known_V) > 0L) {
        latent_blocks <- .known_v_backend_blocks(known_V, "latent")
        independent   <- .known_v_independent_indices(known_V)
        if (length(independent) > 0L) {
          fit_data[["known_v_independent_n"]]     <- length(independent)
          fit_data[["known_v_independent_index"]] <- independent
        }
        for (b in seq_along(latent_blocks)) {
          fit_data[[paste0("sampling_B_", b)]] <-
            latent_blocks[[b]][["B"]]
        }
      }
    }

    fit_data <- .add_marginalized_random_effect_row_multiplier_data(
      fit_data = fit_data,
      data     = data
    )

  } else if (.data_outcome_type(data) == "bin") {
    # always include ai, ci, n1i, and n2i
    fit_data <- list(
      ai  = data[["outcome"]][["ai"]],
      ci  = data[["outcome"]][["ci"]],
      n1i = data[["outcome"]][["n1i"]],
      n2i = data[["outcome"]][["n2i"]]
    )
  } else if (.data_outcome_type(data) == "pois") {
    # always include x1i, x2i, t1i, and t2i
    fit_data <- list(
      x1i = data[["outcome"]][["x1i"]],
      x2i = data[["outcome"]][["x2i"]],
      t1i = data[["outcome"]][["t1i"]],
      t2i = data[["outcome"]][["t2i"]]
    )
  }

  ### add common data
  # add number of estimates
  fit_data[["K"]] <- nrow(data[["outcome"]])

  # add cluster for 3lvl models
  if (.is_data_multilevel(data)) {
    fit_data[["cluster"]] <- data[["outcome"]][["cluster"]]
  }

  # add weights for weighted models
  if (.is_data_weights(data)) {
    fit_data[["weight"]] <- data[["outcome"]][["weights"]]
  }



  return(fit_data)
}
.create_fit_formula_list       <- function(data, parameter) {

  formula <- attr(data[[parameter]], "formula")

  # for scale regression - modify formula for exponential parameterization
  # (required for the exponential parameterization trick in BayesTools::JAGS_fit)
  if (parameter == "scale") {
    formula <- .create_fit_scale_formula(formula)
  }

  return(formula)
}
.create_fit_scale_formula      <- function(formula) {

  # check whether the intercept was removed from the formula
  # if so, warn the user and return it back (intercepts cannot be omitted from scale models)
  if (attr(terms(formula), "intercept") == 0) {
    warning("Intercept cannot be omitted from scale models (the regression estimates a multiplicative constant for the intercept). The intercept removal term has been ignored.", call. = FALSE)
    formula <- BayesTools::formula_add_intercept(formula)
  }

  # add the corresponding attribute
  attr(formula, "log(intercept)") <- TRUE

  return(formula)
}
.create_fit_formula_data_list  <- function(data, parameter) {
  return(data[[parameter]])
}
.create_fit_formula_prior_list <- function(priors, parameter) {
  return(priors[[parameter]])
}
.create_model_syntax           <- function(data, priors) {

  .check_glmm_no_bias_priors(data, priors)

  ### extract structural information about the model
  is_mods           <- .is_data_mods(data)
  is_scale          <- .is_data_scale(data)
  is_random         <- .is_data_random(data)
  is_multilevel     <- .is_data_multilevel(data)
  is_weights        <- .is_data_weights(data)
  is_PET            <- .is_priors_PET(priors)
  is_PEESE          <- .is_priors_PEESE(priors)
  is_weightfunction <- .is_priors_weightfunction(priors)
  outcome_type      <- .data_outcome_type(data)
  effect_direction  <- .data_effect_direction(data)
  is_known_v                   <- .is_data_known_v(data)
  known_v_rank                 <- .data_known_v_rank(data)
  known_v_parameterization     <- .data_known_v_parameterization(data)
  known_v_backend              <- .data_known_v_effective_backend(data)
  is_known_v_latent            <- is_known_v && known_v_backend %in% c("latent", "diagonal")
  is_known_v_whitened          <- is_known_v && known_v_backend == "whitened"
  is_known_v_block_mvn         <- is_known_v && known_v_backend == "block_mvn"
  has_marginalized_random      <- .data_has_marginalized_random_effects(data)

  if (is_known_v_whitened && is_scale) {
    stop(
      "known_v_parameterization = 'whitened' is currently available only without scale regression.",
      call. = FALSE
    )
  }

  ### create the model syntax
  model_syntax <- "model{\n"

  selection_spec <- NULL
  if (is_weightfunction) {
    selection_spec <- .selection_spec(
      priors           = priors,
      yi               = data[["outcome"]][["yi"]],
      sei              = data[["outcome"]][["sei"]],
      effect_direction = effect_direction,
      signed_data      = TRUE
    )
  }

  ### the main model parameters are created automatically via BayesTools::JAGS_fit
  # - mu  (!is_mods and !is_random)  / mu[i]  (is_mods or is_random)
  # - tau (!is_scale) / tau[i] (is_scale)
  # - rho (is_multilevel)
  # for publication bias
  # - PET            (is_PET)
  # - PEESE          (is_PEESE)
  # - weightfunction (is_weightfunction)

  ### prepare multilevel heterogeneity parameter (non-regression = defined outside of the for loop)
  # for multilevel: specify heterogeneity allocation
  # tau_within  - within-cluster variance (estimate-level variance)
  # tau_between - between-cluster variance (cluster-level variance)
  if (is_multilevel & !is_scale) {
    model_syntax <- paste0(model_syntax, "tau_within  = tau * sqrt(1-rho)\n")
    model_syntax <- paste0(model_syntax, "tau_between = tau * sqrt(rho)\n")
  }

  tau_total_node   <- if (is_scale) "tau[i]" else "tau"
  tau_within_node  <- if (is_multilevel) {
    if (is_scale) "tau_within[i]" else "tau_within"
  } else {
    tau_total_node
  }
  tau_between_node <- if (is_scale) "tau_between[i]" else "tau_between"

  if (is_weightfunction && isTRUE(selection_spec[["jags_use_step_switch"]])) {
    model_syntax <- paste0(
      model_syntax,
      "sel_kernel_mode_active = ",
      selection_spec[["jags_kernel_mode_expr"]],
      "\n"
    )
  }

  if (is_known_v_latent && known_v_rank > 0L) {
    known_V       <- .data_known_v_data(data)
    latent_blocks <- .known_v_backend_blocks(known_V, "latent")
    independent   <- .known_v_independent_indices(known_V)
    if (length(independent) > 0L) {
      model_syntax <- paste0(
        model_syntax,
        "for(j in 1:known_v_independent_n){\n",
        "  sampling_dependency[known_v_independent_index[j]] = 0\n",
        "}\n"
      )
    }
    for (b in seq_along(latent_blocks)) {
      block <- latent_blocks[[b]]
      for (j in seq_along(block[["index"]])) {
        model_syntax <- paste0(
          model_syntax,
          "sampling_dependency[", block[["index"]][[j]], "] = inprod(",
          "sampling_B_", b, "[", j, ",1:", block[["rank"]], "],",
          "sampling_z[", block[["z_start"]], ":", block[["z_end"]], "])\n"
        )
      }
    }
  }

  ### enter the main block
  model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")

  ### prepare effect size parameter
  # flip effect size direction (needed for selection models, PET, and PEESE models)
  # (done for everything for consistency, data are flipped within `.create_fit_data()`)
  if (is_mods || is_random) {
    mu_estimate <- ifelse(effect_direction == "negative", paste0("- mu[i]"), "mu[i]")
  } else {
    mu_estimate <- ifelse(effect_direction == "negative", "- mu", "mu")
  }
  # add cluster-level effects
  if (is_multilevel) {
    mu_estimate <- paste0(mu_estimate, ifelse(effect_direction == "negative", " - ", " + "),  "gamma[cluster[i]] * ", tau_between_node)
  }
  # add known-V sampling dependency latent factors
  if (is_known_v_latent && known_v_rank > 0L) {
    mu_estimate <- paste0(mu_estimate, " + sampling_dependency[i]")
  }
  # add PET/PEESE
  if (is_PET) {
    mu_estimate <- paste0(mu_estimate, " + PET * sei[i]")
  }
  if (is_PEESE) {
    mu_estimate <- paste0(mu_estimate, " + PEESE * pow(sei[i],2)")
  }
  if (is_known_v_whitened || is_known_v_block_mvn) {
    model_syntax <- paste0(model_syntax, "  mu_observed[i] = ", mu_estimate, "\n")
  }


  ### prepare heterogeneity parameter (regression)
  # in the case of a regression: the scale model is specified on multiplicative scale:
  # tau = intercept * exp(sum(beta_i * x_i))
  # the model uses a trick to separately pass the tau_intercept and multiply it with
  # log_tau scale regression with zero intercept
  if (is_scale) {
    for (scale_spec in .data_scale_component_specs(data)) {
      model_syntax <- paste0(
        model_syntax,
        "  ", scale_spec[["source"]], "[i] = exp(",
        scale_spec[["parameter"]], "[i])\n"
      )
    }
  }
  # for multilevel: specify heterogeneity allocation & dispatch estimate-specific/common parameter
  # tau_within  - within-cluster variance (estimate-level variance)
  # tau_between - between-cluster variance (cluster-level variance)
  # for rest: dispatch estimate-specific/common parameter
  if (is_multilevel) {
    if (is_scale) {
      model_syntax <- paste0(model_syntax, "  tau_within[i]  = tau[i] * sqrt(1-rho)\n")
      model_syntax <- paste0(model_syntax, "  tau_between[i] = tau[i] * sqrt(rho)\n")
      tau_estimate <- tau_within_node
    } else {
      # the variance allocation performed outside of the loop
      tau_estimate <- tau_within_node
    }
  } else if (is_random) {
    tau_estimate <- "0"
  } else {
    tau_estimate <- tau_total_node
  }
  # compute the total conditional variance of the observed estimate
  sampling_var_node          <- if (is_known_v_latent) "sampling_var[i]" else "pow(sei[i],2)"
  marginalized_random_var    <- .data_marginalized_random_variance_expression(data, row_index = "i")
  random_likelihood_var_expr <- if (identical(marginalized_random_var, "0")) {
    sampling_var_node
  } else {
    paste0(sampling_var_node, " + ", marginalized_random_var)
  }
  total_var_expr             <- if (is_random) {
    paste0("( ", random_likelihood_var_expr, " )")
  } else {
    paste0("( ", sampling_var_node, " + pow(", tau_estimate, ",2) )")
  }
  if (is_known_v_block_mvn) {
    tau2_expr <- if (is_random) marginalized_random_var else paste0("pow(", tau_estimate, ",2)")
    model_syntax <- paste0(model_syntax, "  tau2_observed[i] = ", tau2_expr, "\n")
  }


  ### specify model likelihood
  # separate likelihoods for normal / selection / binomial models
  # subversion for weighted/unweighted likelihoods
  if (outcome_type == "norm") {

    # selection and normal models for norm data outcome
    if (is_weightfunction) {
      likelihood_weight_expr <- if (is_weights) "weight[i]" else "1"
      total_sd_node          <- "sel_total_sd[i]"
      model_syntax           <- paste0(
        model_syntax,
        "  ", total_sd_node, " = ", paste0("sqrt", total_var_expr), "\n"
      )
      if (identical(selection_spec[["mode"]], "step") &&
          isTRUE(selection_spec[["jags_use_step_switch"]])) {
        model_syntax <- paste0(
          model_syntax,
          "  yi[i] ~ dselnorm_step_switch(",
          mu_estimate, ",", total_sd_node, ",",
          "sei[i],", likelihood_weight_expr, ",",
          selection_spec[["jags_omega"]], ",",
          "sel_z_lower,sel_z_upper,sel_obs_bin[i],sel_sign,",
          selection_spec[["jags_kernel_mode"]], ",",
          "sel_telescope_probabilities)\n"
        )
      } else if (identical(selection_spec[["mode"]], "step")) {
        model_syntax <- paste0(
          model_syntax,
          "  yi[i] ~ dselnorm_step(",
          mu_estimate, ",", total_sd_node, ",",
          "sei[i],", likelihood_weight_expr, ",",
          selection_spec[["jags_omega"]], ",",
          "sel_z_lower,sel_z_upper,sel_obs_bin[i],sel_sign,",
          "sel_telescope_probabilities)\n"
        )
      } else {
        model_syntax <- paste0(
          model_syntax,
          "  yi[i] ~ dselnorm_kernel(",
          mu_estimate, ",", total_sd_node, ",",
          mu_estimate, ",", total_sd_node, ",",
          "sei[i],", likelihood_weight_expr, ",",
          selection_spec[["jags_omega"]], ",",
          "sel_z_lower,sel_z_upper,sel_obs_bin[i],sel_sign,",
          selection_spec[["jags_alpha"]], ",",
          selection_spec[["jags_phack_kind"]], ",",
          "phack_z_source,phack_z_dest,",
          "sel_segment_bounds,sel_segment_step_bin,sel_segment_phack_region,",
          "sel_kernel_mode)\n"
        )
      }
    } else if (!is_known_v_whitened && !is_known_v_block_mvn) {
      if (is_weights) {
        model_syntax <- paste0(model_syntax, "  yi[i] ~ dwnorm(", mu_estimate, ",", "1/", total_var_expr, ", weight[i])\n")
      } else {
        model_syntax <- paste0(model_syntax, "  yi[i] ~ dnorm(",  mu_estimate, ",", "1/", total_var_expr, ")\n")
      }
    }

  } else if (outcome_type %in% c("bin", "pois")) {

    # estimate-level variance is not marginalized in binomial-normal/poisson-normal models
    mu_estimate <- paste0(mu_estimate, " + theta[i] * ", tau_estimate)

    # specify appropriate link for each type of outcome
    if (outcome_type == "bin") {
      # transform the parameters to the probability scale
      model_syntax <- paste0(model_syntax, "  logit(p1[i]) = logit(pi[i]) + 0.5 * (", mu_estimate, ")\n")
      model_syntax <- paste0(model_syntax, "  logit(p2[i]) = logit(pi[i]) - 0.5 * (", mu_estimate, ")\n")
    } else if (outcome_type == "pois") {
      # transform the parameters to the probability scale
      model_syntax <- paste0(model_syntax, "  log(r1[i]) = phi[i] + 0.5 * (", mu_estimate, ") + log(t1i[i])\n")
      model_syntax <- paste0(model_syntax, "  log(r2[i]) = phi[i] - 0.5 * (", mu_estimate, ") + log(t2i[i])\n")
    }

    # specify appropriate likelihoods
    if (outcome_type == "bin") {
      # the observed data
      if (is_weights) {
        model_syntax <- paste0(model_syntax, "  ai[i] ~ dwbinom(p1[i], n1i[i], weight[i])\n")
        model_syntax <- paste0(model_syntax, "  ci[i] ~ dwbinom(p2[i], n2i[i], weight[i])\n")
      } else {
        model_syntax <- paste0(model_syntax, "  ai[i] ~ dbinom(p1[i], n1i[i])\n")
        model_syntax <- paste0(model_syntax, "  ci[i] ~ dbinom(p2[i], n2i[i])\n")
      }
    } else if (outcome_type == "pois") {
      # the observed data
      if (is_weights) {
        model_syntax <- paste0(model_syntax, "  x1i[i] ~ dwpois(r1[i], weight[i])\n")
        model_syntax <- paste0(model_syntax, "  x2i[i] ~ dwpois(r2[i], weight[i])\n")
      } else {
        model_syntax <- paste0(model_syntax, "  x1i[i] ~ dpois(r1[i])\n")
        model_syntax <- paste0(model_syntax, "  x2i[i] ~ dpois(r2[i])\n")
      }
    }



  }


  model_syntax <- paste0(model_syntax, "}\n")
  if (outcome_type == "norm" && is_known_v_whitened) {
    known_V          <- .data_known_v_data(data)
    whitening_blocks <- .known_v_backend_blocks(known_V, "whitened")
    independent <- .known_v_independent_indices(known_V)
    if (length(independent) > 0L) {
      independent_extra <- if (is_random) {
        .data_marginalized_random_variance_expression(
          data,
          row_index = "known_v_independent_index[j]"
        )
      } else {
        "pow(tau,2)"
      }
      model_syntax <- paste0(
        model_syntax,
        "for(j in 1:known_v_independent_n){\n",
        "  known_v_independent_y[j] ~ dnorm(",
        "mu_observed[known_v_independent_index[j]], 1/( ",
        "known_v_independent_var[j] + ", independent_extra, " ))\n",
        "}\n"
      )
    }
    for (b in seq_along(whitening_blocks)) {
      block <- whitening_blocks[[b]]
      index <- block[["index"]]
      for (j in seq_along(index)) {
        model_syntax <- paste0(
          model_syntax,
          "whitening_mu_source_", b, "[", j, "] = mu_observed[", index[[j]], "]\n"
        )
      }
      model_syntax <- paste0(
        model_syntax,
        "for(j in 1:", block[["size"]], "){\n",
        "  whitening_mu_", b, "[j] = inprod(whitening_matrix_", b,
        "[j,1:", block[["size"]], "],whitening_mu_source_", b,
        "[1:", block[["size"]], "])\n"
      )
      extra <- if (is_random) {
        if (has_marginalized_random) {
          .data_marginalized_random_variance_expression(
            data,
            row_index = as.character(index[[1L]])
          )
        } else {
          "0"
        }
      } else {
        "pow(tau,2)"
      }
      model_syntax <- paste0(
        model_syntax,
        "  whitening_y_", b, "[j] ~ dnorm(whitening_mu_", b,
        "[j], 1/( whitening_var_", b, "[j] + ", extra, " ))\n",
        "}\n"
      )
    }
  }
  if (outcome_type == "norm" && is_known_v_block_mvn) {
    known_V <- .data_known_v_data(data)
    block_mvn_blocks <- .known_v_backend_blocks(known_V, "block_mvn")

    independent <- .known_v_independent_indices(known_V)
    if (length(independent) > 0L) {
      model_syntax <- paste0(
        model_syntax,
        "for(j in 1:known_v_independent_n){\n",
        "  known_v_independent_y[j] ~ dnorm(",
        "mu_observed[known_v_independent_index[j]], 1/( ",
        "known_v_independent_var[j] + ",
        "tau2_observed[known_v_independent_index[j]] ))\n",
        "}\n"
      )
    }

    for (b in seq_along(block_mvn_blocks)) {
      block <- block_mvn_blocks[[b]]
      idx   <- block[["index"]]

      for (j in seq_along(idx)) {
        model_syntax <- paste0(
          model_syntax,
          "known_v_mu_", b, "[", j, "] = mu_observed[", idx[[j]], "]\n",
          "known_v_tau2_", b, "[", j, "] = tau2_observed[", idx[[j]], "]\n"
        )
      }
      model_syntax <- paste0(
        model_syntax,
        "known_v_y_", b, "[1:", block[["size"]], "] ~ dknown_v_mnorm(",
        "known_v_mu_", b, "[1:", block[["size"]], "],",
        "known_v_tau2_", b, "[1:", block[["size"]], "],",
        "known_v_lower_", b, "[1:", length(block[["v_lower"]]), "])\n"
      )
    }
  }
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}

.fit <- function(object, extend = FALSE) {

  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  data               <- object[["data"]]
  priors             <- object[["priors"]]

  errors   <- NULL
  warnings <- NULL

  ### create model base
  fit_priors       <- .create_fit_priors(data = data, priors = priors)
  fit_data         <- .create_fit_data(data = data, priors = priors)
  fit_formula_args <- .create_jags_formula_args(data = data, priors = priors)

  ### generate the model syntax
  model_syntax <- .create_model_syntax(data = data, priors = priors)

  ### fit the model
  if (!extend || length(object[["fit"]]) == 0) {

    fit <- BayesTools::JAGS_fit(
      model_syntax          = model_syntax,
      data                  = fit_data,
      prior_list            = fit_priors,
      formula_list          = .optional_jags_list(fit_formula_args[["formula_list"]]),
      formula_data_list     = .optional_jags_list(fit_formula_args[["formula_data_list"]]),
      formula_prior_list    = .optional_jags_list(fit_formula_args[["formula_prior_list"]]),
      formula_scale_list    = .optional_jags_list(fit_formula_args[["formula_scale_list"]]),
      formula_random_prior_list           = .optional_jags_list(fit_formula_args[["formula_random_prior_list"]]),
      formula_random_effects_compile_list = .optional_jags_list(fit_formula_args[["formula_random_effects_compile_list"]]),
      add_parameters                      = .optional_jags_character(fit_formula_args[["add_parameters"]]),
      chains                = fit_control[["chains"]],
      adapt                 = fit_control[["adapt"]],
      burnin                = fit_control[["burnin"]],
      sample                = fit_control[["sample"]],
      thin                  = fit_control[["thin"]],
      autofit               = fit_control[["autofit"]],
      autofit_control       = autofit_control,
      parallel              = fit_control[["parallel"]],
      cores                 = fit_control[["cores"]],
      silent                = fit_control[["silent"]],
      seed                  = fit_control[["seed"]],
      required_packages     = "RoBMA",
      is_JASP               = object[["is_JASP"]],
      is_JASP_prefix        = object[["is_JASP_prefix"]]
    )

  } else {

    fit <- BayesTools::JAGS_extend(
      fit                = object[["fit"]],
      autofit_control    = autofit_control,
      parallel           = fit_control[["parallel"]],
      cores              = fit_control[["cores"]],
      silent             = fit_control[["silent"]],
      seed               = fit_control[["seed"]]
    )

  }


  # assess the model fit and deal with errors
  if (inherits(fit, "error")) {

    if(grepl("Unknown function", fit$message))
      stop("The RoBMA module could not be loaded. Check whether the RoBMA package was installed correctly and whether 'RoBMA::RoBMA.private$module_location' contains path to the RoBMA JAGS module.", call. = FALSE)

    fit_error_message       <- fit$message
    fit                     <- list()
    attr(fit, "prior_list") <- fit_priors

    converged      <- FALSE
    has_posterior  <- FALSE
    errors         <- c(errors, fit_error_message)

  } else {

    has_posterior <- TRUE
    check_fit     <- BayesTools::JAGS_check_convergence(
      fit            = fit,
      prior_list     = attr(fit, "prior_list"),
      add_parameters = .convergence_structural_parameters(priors),
      max_Rhat             = convergence_checks[["max_Rhat"]],
      min_ESS              = convergence_checks[["min_ESS"]],
      max_error            = convergence_checks[["max_error"]],
      max_SD_error         = convergence_checks[["max_SD_error"]],
      check_indicators     = isTRUE(convergence_checks[["check_indicators"]]),
      monitor              = convergence_checks[["monitor"]],
      allow_not_assessable = isTRUE(convergence_checks[["allow_not_assessable"]])
    )
    warnings  <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
    converged <- check_fit

  }

  # add results
  fit$errors        <- errors
  fit$warnings      <- warnings
  fit$converged     <- converged
  fit$has_posterior <- has_posterior

  return(fit)
}


# Stop before post-processing when fitting produced no posterior.
.stop_fit_errors <- function(fit) {

  fit_errors <- fit[["errors"]]
  fit_failed <- identical(fit[["has_posterior"]], FALSE) ||
    length(fit_errors) > 0L

  if (!fit_failed) {
    return(invisible(TRUE))
  }
  if (length(fit_errors) > 0L) {
    stop(paste(fit_errors, collapse = "\n"), call. = FALSE)
  }

  stop(
    "Model fitting failed before posterior samples were produced.",
    call. = FALSE
  )
}


.is_priors_PET            <- function(priors) {

  if (is.null(priors[["outcome"]][["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["outcome"]][["bias"]]))
    return(any(sapply(priors[["outcome"]][["bias"]], is.prior.PET)))

  return(is.prior.PET(priors[["outcome"]][["bias"]]))
}
.is_priors_PEESE          <- function(priors) {

  if (is.null(priors[["outcome"]][["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["outcome"]][["bias"]]))
    return(any(sapply(priors[["outcome"]][["bias"]], is.prior.PEESE)))

  return(is.prior.PEESE(priors[["outcome"]][["bias"]]))
}
.is_priors_weightfunction <- function(priors) {

  if (is.null(priors[["outcome"]][["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["outcome"]][["bias"]]))
    return(any(sapply(priors[["outcome"]][["bias"]], .prior_is_selection_kernel)))

  return(.prior_is_selection_kernel(priors[["outcome"]][["bias"]]))
}
.convergence_structural_parameters <- function(priors) {

  if (!.is_priors_weightfunction(priors)) {
    return(character())
  }

  bias_priors <- .selection_bias_priors(priors)
  has_selection <- vapply(
    bias_priors,
    .prior_has_selection,
    logical(1)
  )
  p_cuts <- BayesTools::weightfunctions_mapping(
    prior_list = bias_priors[has_selection],
    cuts_only  = TRUE,
    one_sided  = TRUE
  )

  public_reference <- paste0(
    "omega[", p_cuts[[1L]], ",", p_cuts[[2L]], "]"
  )
  return(c("omega[1]", public_reference))
}
.is_priors_bias           <- function(priors) {
  return(.is_priors_PET(priors) || .is_priors_PEESE(priors) || .is_priors_weightfunction(priors))
}
.check_glmm_no_bias_priors <- function(data, priors) {

  if (.data_outcome_type(data) %in% c("bin", "pois") && .is_priors_bias(priors)) {
    stop("Publication-bias priors are not supported for GLMM outcomes.",
         call. = FALSE)
  }

  invisible(TRUE)
}
.is_data_multilevel       <- function(data) {
  return(isTRUE(attr(data, "cluster")))
}
.is_data_mods             <- function(data) {
  return(isTRUE(attr(data, "mods")))
}
.is_data_scale            <- function(data) {
  return(isTRUE(attr(data, "scale")))
}
.data_scale_components    <- function(data) {

  if (!.is_data_scale(data)) {
    return(list())
  }

  scale <- data[["scale"]]
  if (inherits(scale, "RoBMA_scale_components")) {
    return(scale)
  }

  out <- list(tau = scale)
  class(out) <- c("RoBMA_scale_components", "list")

  return(out)
}
.data_scale_component_specs <- function(data) {

  components <- .data_scale_components(data)
  if (length(components) == 0L) {
    return(list())
  }

  out <- lapply(names(components), function(name) {
    component <- components[[name]]
    source    <- attr(component, "source")
    parameter <- attr(component, "parameter")
    if (is.null(source)) {
      source <- if (identical(name, "tau")) "tau" else paste0("tau_", name)
    }
    if (is.null(parameter)) {
      parameter <- if (identical(name, "tau")) "log_tau" else paste0("log_tau_", name)
    }
    component_name <- attr(component, "component_name", exact = TRUE)
    scale_name     <- attr(component, "scale_name", exact = TRUE)
    aliases        <- attr(component, "aliases", exact = TRUE)
    display_name   <- component_name
    if (is.null(display_name) || length(display_name) != 1L ||
        is.na(display_name) || !nzchar(display_name)) {
      display_name <- scale_name
    }
    if (is.null(display_name) || length(display_name) != 1L ||
        is.na(display_name) || !nzchar(display_name)) {
      display_name <- name
    }
    if (is.null(aliases)) {
      aliases <- unique(c(name, component_name, scale_name))
      aliases <- aliases[!is.na(aliases) & nzchar(aliases)]
    }
    list(
      name           = name,
      display_name   = display_name,
      component_name = component_name,
      scale_name     = scale_name,
      aliases        = aliases,
      source         = source,
      parameter      = parameter,
      data           = component,
      formula        = attr(component, "formula")
    )
  })
  names(out) <- names(components)

  return(out)
}
.data_scale_formula_parameters <- function(data) {

  vapply(.data_scale_component_specs(data), `[[`, character(1), "parameter")
}
.data_scale_formula_sources <- function(data) {

  vapply(.data_scale_component_specs(data), `[[`, character(1), "source")
}
.is_data_random           <- function(data) {
  return(isTRUE(attr(data, "random")))
}
.is_data_weights          <- function(data) {
  return(isTRUE(attr(data, "weights")))
}
.is_data_known_v          <- function(data) {
  return(isTRUE(attr(data, "known_V")))
}
.data_known_v_data        <- function(data) {
  known_V <- attr(data, "known_V_data")
  if (.is_data_known_v(data)) {
    if (is.null(known_V)) {
      stop("Internal error: known-V metadata are missing.", call. = FALSE)
    }
    .validate_prepared_known_v(known_V)
  }
  return(known_V)
}
.data_known_v_parameterization  <- function(data) {
  known_V <- .data_known_v_data(data)
  if (is.null(known_V)) {
    return("latent")
  }
  return(.known_v_parameterization(known_V))
}
.data_known_v_effective_backend <- function(data) {
  known_V <- .data_known_v_data(data)
  if (is.null(known_V)) {
    return(.data_known_v_parameterization(data))
  }
  return(.known_v_effective_backend(known_V))
}
.is_data_known_v_backend <- function(data, backend) {
  return(
    .is_data_known_v(data) &&
      .data_known_v_effective_backend(data) %in% backend
  )
}
.is_data_known_v_parameterization <- function(data, parameterization) {
  return(.is_data_known_v(data) && .data_known_v_parameterization(data) == parameterization)
}
.data_known_v_rank        <- function(data) {
  known_V <- .data_known_v_data(data)
  if (is.null(known_V)) {
    return(0L)
  }
  return(.known_v_rank(known_V))
}
.data_known_v_correlated  <- function(data) {
  known_V <- .data_known_v_data(data)
  if (is.null(known_V)) {
    return(FALSE)
  }
  return(.known_v_is_correlated(known_V))
}
.data_outcome_type        <- function(data) {
  return(attr(data, "outcome_type"))
}
.data_measure             <- function(data) {
  return(attr(data, "measure"))
}
.data_effect_direction    <- function(data) {
  return(attr(data, "effect_direction"))
}
.data_standardize_continuous_predictors <- function(data) {
  return(attr(data, "standardize_continuous_predictors"))
}

.is_PET            <- function(object) .is_priors_PET(object[["priors"]])
.is_PEESE          <- function(object) .is_priors_PEESE(object[["priors"]])
.is_weightfunction <- function(object) .is_priors_weightfunction(object[["priors"]])
.is_bias           <- function(object) .is_priors_bias(object[["priors"]])
.is_RoBMA          <- function(object) inherits(object, "RoBMA")
.is_BMA            <- function(object) inherits(object, "BMA.norm") || inherits(object, "BMA.glmm")
.is_robust_RoBMA   <- function(object) .is_RoBMA(object) && !.is_BMA(object)
.is_multilevel     <- function(object) .is_data_multilevel(object[["data"]])
.is_mods           <- function(object) .is_data_mods(object[["data"]])
.is_scale          <- function(object) .is_data_scale(object[["data"]])
.is_random         <- function(object) .is_data_random(object[["data"]])
.is_weights        <- function(object) .is_data_weights(object[["data"]])
.outcome_type      <- function(object) .data_outcome_type(object[["data"]])
.measure           <- function(object) .data_measure(object[["data"]])
.effect_direction  <- function(object) .data_effect_direction(object[["data"]])
.standardize_continuous_predictors  <- function(object) .data_standardize_continuous_predictors(object[["data"]])
