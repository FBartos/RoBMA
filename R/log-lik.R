# ============================================================================ #
# Estimate-Unit Log-Likelihood Dispatch
# ============================================================================ #



# ---------------------------------------------------------------------------- #
# .log_lik.brma
# ---------------------------------------------------------------------------- #
#
# Compute the log-likelihood matrix for a brma object.
#
# @param object brma object.
# @param unit   character; output/deletion unit.
#
# @return S x K or S x G matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.log_lik.brma <- function(object, unit = "estimate",
                          caller = ".log_lik.brma()") {

  unit <- .normalize_unit(unit)
  .check_log_lik_target_available(object, unit, caller)

  if (unit == "estimate") {
    return(.log_lik_estimate.brma(object))
  } else {
    return(.log_lik_cluster.brma(object))
  }
}



# ---------------------------------------------------------------------------- #
# .pdf.brma
# ---------------------------------------------------------------------------- #
#
# Back-end wrapper kept for internal callers that use the old helper name.
#
# @param object brma object.
# @param unit   character; output/deletion unit.
#
# @return S x K or S x G matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.pdf.brma <- function(object, unit = "estimate") {

  return(.log_lik.brma(object = object, unit = unit))
}



# ---------------------------------------------------------------------------- #
# .log_lik_from_posterior_samples
# ---------------------------------------------------------------------------- #
#
# Compute pointwise log-likelihoods from an explicit posterior sample matrix.
#
# This is the shared evaluator used by regular LOO/logLik paths and by IWMDE
# candidate rows. It owns only posterior-row evaluation and dispatch; the actual
# likelihood math stays in the existing outcome and cluster helpers.
#
# ---------------------------------------------------------------------------- #
.log_lik_from_posterior_samples <- function(fit, posterior_samples, data, priors,
                                             unit = "estimate",
                                             add_metadata = FALSE,
                                             data_hash = NULL) {

  unit              <- .normalize_unit(unit)
  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  setup             <- .log_lik_posterior_setup(
    fit                  = fit,
    posterior_samples    = posterior_samples,
    data                 = data,
    priors               = priors,
    unit                 = unit,
    data_hash            = data_hash
  )

  log_lik <- if (unit == "estimate") {
    .log_lik_estimate_from_setup(setup)
  } else {
    .log_lik_cluster_from_setup(setup)
  }

  if (!add_metadata) {
    return(log_lik)
  }

  if (unit == "estimate") {
    colnames(log_lik) <- paste0("log_lik[", seq_len(setup[["K"]]), "]")
    attr(log_lik, "RoBMA_target") <- .estimate_log_lik_target_metadata(
      setup     = setup,
      data_hash = data_hash
    )
  } else {
    log_lik <- .add_cluster_log_lik_metadata(
      log_lik         = log_lik,
      cluster_indices = setup[["cluster"]],
      data_hash       = data_hash
    )
  }

  return(log_lik)
}



.log_lik_from_evaluated_predictors <- function(fit, data, priors,
                                               mu_samples,
                                               tau_within_samples,
                                               tau_between_samples = NULL,
                                               posterior_samples = NULL,
                                               unit = "estimate",
                                               add_metadata = FALSE,
                                               data_hash = NULL) {

  unit  <- .normalize_unit(unit)
  setup <- .log_lik_evaluated_setup(
    fit                  = fit,
    data                 = data,
    priors               = priors,
    unit                 = unit,
    data_hash            = data_hash,
    mu_samples           = mu_samples,
    tau_within_samples   = tau_within_samples,
    tau_between_samples  = tau_between_samples,
    posterior_samples    = posterior_samples
  )

  log_lik <- if (unit == "estimate") {
    .log_lik_estimate_from_setup(setup)
  } else {
    .log_lik_cluster_from_setup(setup)
  }

  if (!add_metadata) {
    return(log_lik)
  }

  if (unit == "estimate") {
    colnames(log_lik) <- paste0("log_lik[", seq_len(setup[["K"]]), "]")
    attr(log_lik, "RoBMA_target") <- .estimate_log_lik_target_metadata(
      setup     = setup,
      data_hash = data_hash
    )
  } else {
    log_lik <- .add_cluster_log_lik_metadata(
      log_lik         = log_lik,
      cluster_indices = setup[["cluster"]],
      data_hash       = data_hash
    )
  }

  return(log_lik)
}


.log_lik_from_evaluated_predictors_sum <- function(fit, data, priors,
                                                   mu_samples,
                                                   tau_within_samples,
                                                   tau_between_samples = NULL,
                                                   posterior_samples = NULL,
                                                   unit = "estimate",
                                                   data_hash = NULL) {

  unit  <- .normalize_unit(unit)
  setup <- .log_lik_evaluated_setup(
    fit                  = fit,
    data                 = data,
    priors               = priors,
    unit                 = unit,
    data_hash            = data_hash,
    mu_samples           = mu_samples,
    tau_within_samples   = tau_within_samples,
    tau_between_samples  = tau_between_samples,
    posterior_samples    = posterior_samples
  )

  if (unit == "estimate") {
    return(.log_lik_estimate_sum_from_setup(setup))
  }

  return(.log_lik_cluster_sum_from_setup(setup))
}


.log_lik_glmm_conditional_sum_from_evaluated_predictors <- function(fit, data,
                                                                    posterior_samples,
                                                                    mu_samples,
                                                                    tau_within_samples,
                                                                    tau_between_samples = NULL) {

  outcome_type <- .data_outcome_type(data)
  if (!outcome_type %in% c("bin", "pois")) {
    stop("Conditional GLMM log-likelihood requires binomial or Poisson outcomes.",
         call. = FALSE)
  }

  K           <- nrow(data[["outcome"]])
  mu_samples  <- as.matrix(mu_samples)
  tau_within  <- as.matrix(tau_within_samples)

  if (ncol(mu_samples) != K || ncol(tau_within) != K) {
    stop("Evaluated predictor matrices must have one column per estimate.",
         call. = FALSE)
  }
  if (nrow(tau_within) != nrow(mu_samples)) {
    stop("Evaluated predictor matrices must have matching row counts.",
         call. = FALSE)
  }

  if (.is_data_multilevel(data)) {
    if (is.null(tau_between_samples)) {
      stop("Conditional multilevel GLMM log-likelihood requires 'tau_between_samples'.",
           call. = FALSE)
    }
    tau_between <- as.matrix(tau_between_samples)
    if (nrow(tau_between) != nrow(mu_samples) || ncol(tau_between) != K) {
      stop("Evaluated predictor matrices must have matching dimensions.",
           call. = FALSE)
    }
    mu_samples <- mu_samples + .evaluate.brma.cluster_effects(
      fit               = fit,
      tau_between       = tau_between,
      cluster           = data[["outcome"]][["cluster"]],
      same_data         = TRUE,
      effect_direction  = .data_effect_direction(data),
      posterior_samples = posterior_samples
    )
  }

  mu_samples <- mu_samples + .evaluate.brma.theta.glmm(
    fit               = fit,
    tau_within        = tau_within,
    same_data         = TRUE,
    K                 = K,
    posterior_samples = posterior_samples
  )
  weights <- if (.is_data_weights(data)) {
    data[["outcome"]][["weights"]]
  } else {
    NULL
  }

  if (outcome_type == "bin") {
    return(.outcome_pdf_sum.binom_conditional(
      ai             = data[["outcome"]][["ai"]],
      ci             = data[["outcome"]][["ci"]],
      n1i            = data[["outcome"]][["n1i"]],
      n2i            = data[["outcome"]][["n2i"]],
      mu_samples     = mu_samples,
      logit_baserate = .evaluate.brma.baserate(
        fit               = fit,
        K                 = K,
        posterior_samples = posterior_samples
      ),
      weights        = weights
    ))
  }

  return(.outcome_pdf_sum.pois_conditional(
    x1i        = data[["outcome"]][["x1i"]],
    x2i        = data[["outcome"]][["x2i"]],
    t1i        = data[["outcome"]][["t1i"]],
    t2i        = data[["outcome"]][["t2i"]],
    mu_samples = mu_samples,
    log_phi    = .evaluate.brma.lograte(
      fit               = fit,
      K                 = K,
      posterior_samples = posterior_samples
    ),
    weights    = weights
  ))
}


.log_lik_from_posterior_samples_sum <- function(fit, posterior_samples,
                                                 data, priors,
                                                 unit = "estimate",
                                                 data_hash = NULL) {

  unit  <- .normalize_unit(unit)
  setup <- .log_lik_posterior_setup(
    fit                  = fit,
    posterior_samples    = posterior_samples,
    data                 = data,
    priors               = priors,
    unit                 = unit,
    data_hash            = data_hash
  )

  if (unit == "estimate") {
    return(.log_lik_estimate_sum_from_setup(setup))
  }

  return(.log_lik_cluster_sum_from_setup(setup))
}



.log_lik_evaluated_setup <- function(fit, data, priors, unit, data_hash,
                                      mu_samples, tau_within_samples,
                                      tau_between_samples, posterior_samples) {

  is_multilevel <- .is_data_multilevel(data)
  is_random     <- .is_data_random(data)
  K             <- nrow(data[["outcome"]])

  if (is_random && !.is_data_known_v(data)) {
    stop(
      "Log-likelihood evaluation is not implemented for brma.mv() ",
      "random-formula models yet.",
      call. = FALSE
    )
  }

  if (unit == "cluster" && !is_multilevel) {
    stop("Cluster-unit log-likelihood is only available for multilevel models.",
         call. = FALSE)
  }

  mu_samples         <- as.matrix(mu_samples)
  tau_within_samples <- as.matrix(tau_within_samples)
  if (is.null(tau_between_samples)) {
    tau_between_samples <- matrix(0, nrow = nrow(mu_samples), ncol = K)
  } else {
    tau_between_samples <- as.matrix(tau_between_samples)
  }

  if (ncol(mu_samples) != K ||
      ncol(tau_within_samples) != K ||
      ncol(tau_between_samples) != K) {
    stop("Evaluated predictor matrices must have one column per estimate.",
         call. = FALSE)
  }
  if (nrow(tau_within_samples) != nrow(mu_samples) ||
      nrow(tau_between_samples) != nrow(mu_samples)) {
    stop("Evaluated predictor matrices must have matching row counts.",
         call. = FALSE)
  }

  if (is.null(posterior_samples)) {
    posterior_samples <- matrix(numeric(0), nrow = nrow(mu_samples), ncol = 0L)
  } else {
    posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  }
  if (nrow(posterior_samples) != nrow(mu_samples)) {
    stop("Posterior samples must match evaluated predictor rows.",
         call. = FALSE)
  }

  active_object <- list(
    fit    = fit,
    data   = data,
    priors = priors
  )
  if (.is_data_known_v(data) && !.known_v_estimate_target_uses_backend(data)) {
    mu_samples <- mu_samples + .evaluate.brma.sampling_dependency(
      fit               = fit,
      data              = data,
      posterior_samples = posterior_samples
    )
  }
  zero_mu <- matrix(0, nrow = nrow(mu_samples), ncol = K)
  setup_sei <- .log_lik_setup_sei(active_object, data)

  return(list(
    fit                  = fit,
    priors               = priors,
    data                 = data,
    yi                   = .outcome_data_yi(active_object),
    sei                  = setup_sei[["likelihood_sei"]],
    selection_sei        = setup_sei[["selection_sei"]],
    K                    = K,
    S                    = nrow(mu_samples),
    mu                   = mu_samples,
    mu_random            = zero_mu,
    tau_within           = tau_within_samples,
    tau_between          = tau_between_samples,
    cluster              = if (is_multilevel) split(seq_len(K), data[["outcome"]][["cluster"]]) else NULL,
    weights              = if (.is_data_weights(data)) data[["outcome"]][["weights"]] else NULL,
    data_hash            = data_hash,
    is_weightfunction    = .is_priors_weightfunction(priors),
    outcome_type         = .data_outcome_type(data),
    effect_direction     = .data_effect_direction(data),
    posterior_samples    = posterior_samples
  ))
}



.log_lik_posterior_setup <- function(fit, posterior_samples, data, priors,
                                     unit, data_hash) {

  .estimate_likelihood_setup_from_parts(
    fit               = fit,
    data              = data,
    priors            = priors,
    posterior_samples = posterior_samples,
    unit              = unit,
    data_hash         = data_hash,
    bias_adjusted     = FALSE,
    caller            = ".log_lik_posterior_setup()"
  )
}



.log_lik_estimate_from_setup <- function(setup) {

  data                <- setup[["data"]]
  priors              <- setup[["priors"]]
  yi                  <- setup[["yi"]]
  sei                 <- setup[["sei"]]
  selection_sei       <- setup[["selection_sei"]]
  K                   <- setup[["K"]]
  mu_samples          <- setup[["mu"]]
  tau_within_samples  <- setup[["tau_within"]]
  is_weightfunction   <- setup[["is_weightfunction"]]
  outcome_type        <- setup[["outcome_type"]]
  effect_direction    <- setup[["effect_direction"]]
  data_weights        <- setup[["weights"]]

  if (.known_v_estimate_target_uses_backend(data)) {
    return(.log_lik_known_v_estimate_target_from_setup(setup))
  }

  if (is.null(selection_sei)) {
    selection_sei <- sei
  }

  if (outcome_type == "norm") {

    if (effect_direction == "negative") {
      mu_samples_pdf <- -mu_samples
      yi_pdf         <- -yi
    } else {
      mu_samples_pdf <- mu_samples
      yi_pdf         <- yi
    }

    if (is_weightfunction) {
      selection_context <- .selection_context_from_parts(
        fit                  = setup[["fit"]],
        data                 = data,
        priors               = priors,
        posterior_samples    = setup[["posterior_samples"]],
        effect_direction     = effect_direction
      )

      log_lik <- .outcome_pdf.selnorm(
        yi                = yi,
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = sei,
        selection_sei     = selection_sei,
        selection_context = selection_context
      )
    } else {
      log_lik <- .outcome_pdf.norm(
        yi         = yi_pdf,
        mu_samples = mu_samples_pdf,
        tau_within = tau_within_samples,
        sei        = sei
      )
    }

    return(.apply_log_lik_weights(log_lik, data_weights))
  }

  if (outcome_type == "bin") {
    return(.outcome_pdf.binom(
      ai         = data[["outcome"]][["ai"]],
      ci         = data[["outcome"]][["ci"]],
      n1i        = data[["outcome"]][["n1i"]],
      n2i        = data[["outcome"]][["n2i"]],
      mu_samples = mu_samples,
      tau_within = tau_within_samples,
      prior_pi   = priors[["outcome"]][["pi"]],
      weights    = data_weights
    ))
  }

  if (outcome_type == "pois") {
    return(.outcome_pdf.pois(
      x1i        = data[["outcome"]][["x1i"]],
      x2i        = data[["outcome"]][["x2i"]],
      t1i        = data[["outcome"]][["t1i"]],
      t2i        = data[["outcome"]][["t2i"]],
      mu_samples = mu_samples,
      tau_within = tau_within_samples,
      prior_phi  = priors[["outcome"]][["phi"]],
      weights    = data_weights
    ))
  }

  stop("Unsupported outcome type for estimate-unit log-likelihood.",
       call. = FALSE)
}


.log_lik_estimate_sum_from_setup <- function(setup) {

  data                <- setup[["data"]]
  priors              <- setup[["priors"]]
  yi                  <- setup[["yi"]]
  sei                 <- setup[["sei"]]
  selection_sei       <- setup[["selection_sei"]]
  mu_samples          <- setup[["mu"]]
  tau_within_samples  <- setup[["tau_within"]]
  is_weightfunction   <- setup[["is_weightfunction"]]
  outcome_type        <- setup[["outcome_type"]]
  effect_direction    <- setup[["effect_direction"]]
  data_weights        <- setup[["weights"]]

  if (.known_v_estimate_target_uses_backend(data)) {
    return(rowSums(.log_lik_known_v_estimate_target_from_setup(setup)))
  }

  if (is.null(selection_sei)) {
    selection_sei <- sei
  }

  if (outcome_type == "norm") {
    if (effect_direction == "negative") {
      mu_samples_pdf <- -mu_samples
      yi_pdf         <- -yi
    } else {
      mu_samples_pdf <- mu_samples
      yi_pdf         <- yi
    }

    if (is_weightfunction) {
      selection_context <- .selection_context_from_parts(
        fit                  = setup[["fit"]],
        data                 = data,
        priors               = priors,
        posterior_samples    = setup[["posterior_samples"]],
        effect_direction     = effect_direction
      )

      return(.outcome_pdf_sum.selnorm(
        yi                = yi,
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = sei,
        selection_sei     = selection_sei,
        selection_context = selection_context,
        weights           = data_weights
      ))
    }

    return(.outcome_pdf_sum.norm(
      yi         = yi_pdf,
      mu_samples = mu_samples_pdf,
      tau_within = tau_within_samples,
      sei        = sei,
      weights    = data_weights
    ))
  }

  if (outcome_type == "bin") {
    return(.outcome_pdf_sum.binom(
      ai         = data[["outcome"]][["ai"]],
      ci         = data[["outcome"]][["ci"]],
      n1i        = data[["outcome"]][["n1i"]],
      n2i        = data[["outcome"]][["n2i"]],
      mu_samples = mu_samples,
      tau_within = tau_within_samples,
      prior_pi   = priors[["outcome"]][["pi"]],
      weights    = data_weights
    ))
  }

  if (outcome_type == "pois") {
    return(.outcome_pdf_sum.pois(
      x1i        = data[["outcome"]][["x1i"]],
      x2i        = data[["outcome"]][["x2i"]],
      t1i        = data[["outcome"]][["t1i"]],
      t2i        = data[["outcome"]][["t2i"]],
      mu_samples = mu_samples,
      tau_within = tau_within_samples,
      prior_phi  = priors[["outcome"]][["phi"]],
      weights    = data_weights
    ))
  }

  stop("Unsupported outcome type for estimate-unit log-likelihood.",
       call. = FALSE)
}



# ---------------------------------------------------------------------------- #
# .estimate_likelihood_setup.brma
# ---------------------------------------------------------------------------- #
#
# Extract common posterior quantities for estimate-unit predictive likelihoods.
#
# The estimate-unit target conditions on fitted cluster effects for multilevel
# models and integrates over estimate-level heterogeneity.
#
# @param object brma object.
#
# @return list with model metadata and posterior matrices.
#
# ---------------------------------------------------------------------------- #
.estimate_likelihood_setup_from_parts <- function(fit, data, priors,
                                                  posterior_samples,
                                                  unit = "estimate",
                                                  data_hash = NULL,
                                                  bias_adjusted = FALSE,
                                                  caller = ".estimate_likelihood_setup_from_parts()") {

  unit              <- .normalize_unit(unit)
  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  is_mods           <- .is_data_mods(data)
  is_scale          <- .is_data_scale(data)
  is_random         <- .is_data_random(data)
  is_multilevel     <- .is_data_multilevel(data)
  is_PET            <- .is_priors_PET(priors)
  is_PEESE          <- .is_priors_PEESE(priors)
  is_weightfunction <- .is_priors_weightfunction(priors)
  outcome_type      <- .data_outcome_type(data)
  effect_direction  <- .data_effect_direction(data)
  K                 <- nrow(data[["outcome"]])

  .check_glmm_no_bias_priors(data, priors)

  if (is_random && !.is_data_known_v(data)) {
    stop(
      caller,
      " is not implemented for brma.mv() random-formula models yet.",
      call. = FALSE
    )
  }

  if (unit == "cluster" && !is_multilevel) {
    stop("Cluster-unit log-likelihood is only available for multilevel models.",
         call. = FALSE)
  }

  tau_result <- if (is_random) {
    zero_tau <- matrix(0, nrow = nrow(posterior_samples), ncol = K)
    list(
      tau_total   = zero_tau,
      tau_within  = zero_tau,
      tau_between = zero_tau
    )
  } else {
    .evaluate.brma.tau(
      fit               = fit,
      scale_data        = data[["scale"]],
      scale_formula     = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
      scale_priors      = priors[["scale"]],
      is_scale          = is_scale,
      is_multilevel     = is_multilevel,
      K                 = K,
      posterior_samples = posterior_samples,
      allow_missing_tau = .fixed_tau_prior_value(priors)
    )
  }

  mu_samples <- .evaluate.brma.mu(
    fit               = fit,
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = if (is_random) priors[["location"]] else priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = bias_adjusted,
    K                 = K,
    posterior_samples = posterior_samples,
    priors            = priors
  )
  mu_random_samples <- matrix(0, nrow = nrow(mu_samples), ncol = K)

  fit_data <- NULL
  if (is_multilevel || is_weightfunction) {
    fit_data <- .create_fit_data(data = data, priors = priors)
  }

  cluster <- if (is_multilevel) {
    split(seq_len(K), data[["outcome"]][["cluster"]])
  } else {
    NULL
  }

  if (unit == "estimate" && is_multilevel) {
    cluster_effects <- .evaluate.brma.cluster_effects(
      fit               = fit,
      tau_between       = tau_result[["tau_between"]],
      cluster           = if (is.null(fit_data)) data[["outcome"]][["cluster"]] else fit_data[["cluster"]],
      same_data         = TRUE,
      effect_direction  = effect_direction,
      posterior_samples = posterior_samples
    )
    mu_samples         <- mu_samples + cluster_effects
    mu_random_samples <- mu_random_samples + cluster_effects
  }
  if (unit == "estimate" && is_random) {
    random_effects <- .evaluate.brma.random_effects(
      fit               = fit,
      data              = data,
      priors            = priors,
      posterior_samples = posterior_samples,
      same_data         = TRUE,
      required          = TRUE
    )
    mu_samples         <- mu_samples + random_effects
    mu_random_samples <- mu_random_samples + random_effects
  }

  active_object <- list(
    fit    = fit,
    data   = data,
    priors = priors
  )
  if (.is_data_known_v(data) && !.known_v_estimate_target_uses_backend(data)) {
    mu_samples <- mu_samples + .evaluate.brma.sampling_dependency(
      fit               = fit,
      data              = data,
      posterior_samples = posterior_samples
    )
  }
  setup_sei <- .log_lik_setup_sei(active_object, data)

  return(list(
    fit               = fit,
    priors            = priors,
    data              = data,
    fit_data          = fit_data,
    yi                = .outcome_data_yi(active_object),
    sei               = setup_sei[["likelihood_sei"]],
    selection_sei     = setup_sei[["selection_sei"]],
    K                 = K,
    S                 = nrow(mu_samples),
    mu                = mu_samples,
    mu_random         = mu_random_samples,
    tau_within        = tau_result[["tau_within"]],
    tau_between       = tau_result[["tau_between"]],
    cluster           = cluster,
    weights           = if (.is_data_weights(data)) data[["outcome"]][["weights"]] else NULL,
    data_hash         = data_hash,
    is_weightfunction = is_weightfunction,
    outcome_type      = outcome_type,
    effect_direction  = effect_direction,
    posterior_samples = posterior_samples
  ))
}


.estimate_likelihood_setup.brma <- function(object, bias_adjusted = FALSE,
                                            posterior_samples = NULL) {

  .estimate_likelihood_setup_from_parts(
    fit               = object[["fit"]],
    data              = object[["data"]],
    priors            = object[["priors"]],
    posterior_samples = .get_posterior_samples(object[["fit"]], posterior_samples),
    unit              = "estimate",
    data_hash         = .get_outcome_hash(object),
    bias_adjusted     = bias_adjusted,
    caller            = ".estimate_likelihood_setup.brma()"
  )
}



# ---------------------------------------------------------------------------- #
# .log_lik_estimate.brma
# ---------------------------------------------------------------------------- #
#
# Estimate-unit likelihood, one contribution per effect-size estimate.
#
# For normal multilevel models this is the model likelihood factor conditional
# on the cluster effect and marginal over estimate-level heterogeneity.
#
# @param object brma object.
#
# @return S x K matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.log_lik_estimate.brma <- function(object) {

  return(.log_lik_from_posterior_samples(
    fit                  = object[["fit"]],
    posterior_samples    = .get_posterior_samples(object[["fit"]]),
    data                 = object[["data"]],
    priors               = object[["priors"]],
    unit                 = "estimate",
    add_metadata         = TRUE,
    data_hash            = .get_outcome_hash(object)
  ))
}


.log_lik_setup_sei <- function(object, data) {

  selection_sei <- .outcome_data_sei(object)
  likelihood_sei <- if (.known_v_estimate_target_uses_backend(data)) {
    selection_sei
  } else {
    .outcome_data_likelihood_sei(object)
  }

  list(
    likelihood_sei = likelihood_sei,
    selection_sei  = selection_sei
  )
}


.estimate_log_lik_target_metadata <- function(setup, data_hash) {

  known_v_estimate_backend <- .known_v_estimate_target_uses_backend(setup[["data"]])
  known_v_schur            <- .known_v_estimate_target_uses_schur_conditioning(setup[["data"]])
  known_V <- if (.is_data_known_v(setup[["data"]])) {
    .data_known_v_data(setup[["data"]])
  } else {
    NULL
  }
  component_sizes <- if (known_v_estimate_backend && !is.null(known_V)) {
    block_indices <- .known_v_dependency_blocks(
      data = setup[["data"]],
      K    = setup[["K"]]
    )
    vapply(block_indices, length, integer(1))
  } else {
    rep(1L, setup[["K"]])
  }
  random_effect_terms <- .data_random_effect_terms(setup[["data"]])
  has_sampled_random <- length(random_effect_terms) > 0L && any(vapply(
    random_effect_terms,
    function(term) identical(term[["compile_mode"]], "sampled"),
    logical(1)
  ))

  list(
    unit                 = "estimate",
    conditioning_depth   = "estimate",
    target               = if (known_v_estimate_backend) "known_v_estimate" else "factorized_estimate",
    n                    = setup[["K"]],
    targets              = seq_len(setup[["K"]]),
    data_hash            = data_hash,
    known_v              = .is_data_known_v(setup[["data"]]),
    known_v_estimate_backend = known_v_estimate_backend,
    known_v_parameterization = if (is.null(known_V)) NULL else known_V[["parameterization"]],
    known_v_parameterization_requested = if (is.null(known_V)) {
      NULL
    } else {
      known_V[["parameterization_requested"]]
    },
    dependency_component_sizes = component_sizes,
    known_v_schur       = known_v_schur,
    random_effects       = if (has_sampled_random) "conditioned" else "none",
    estimate_level_random = if (.data_has_marginalized_random_effects(setup[["data"]])) {
      "marginalized"
    } else if (.data_has_sampled_estimate_level_random_effects(setup[["data"]])) {
      "conditioned"
    } else {
      "none"
    }
  )
}
