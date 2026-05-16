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
.log_lik.brma <- function(object, unit = "estimate") {

  unit <- .normalize_unit(unit)

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
                                            data_hash = NULL,
                                            honor_bias_indicator = FALSE) {

  unit              <- .normalize_unit(unit)
  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  setup             <- .log_lik_posterior_setup(
    fit                  = fit,
    posterior_samples    = posterior_samples,
    data                 = data,
    priors               = priors,
    unit                 = unit,
    data_hash            = data_hash,
    honor_bias_indicator = honor_bias_indicator
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
    attr(log_lik, "RoBMA_target") <- list(
      unit               = "estimate",
      conditioning_depth = "estimate",
      n                  = setup[["K"]],
      targets            = seq_len(setup[["K"]]),
      data_hash          = data_hash
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
                                               data_hash = NULL,
                                               honor_bias_indicator = FALSE) {

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
    posterior_samples    = posterior_samples,
    honor_bias_indicator = honor_bias_indicator
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
    attr(log_lik, "RoBMA_target") <- list(
      unit               = "estimate",
      conditioning_depth = "estimate",
      n                  = setup[["K"]],
      targets            = seq_len(setup[["K"]]),
      data_hash          = data_hash
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
                                                   data_hash = NULL,
                                                   honor_bias_indicator = FALSE) {

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
    posterior_samples    = posterior_samples,
    honor_bias_indicator = honor_bias_indicator
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
                                                data_hash = NULL,
                                                honor_bias_indicator = TRUE) {

  unit  <- .normalize_unit(unit)
  setup <- .log_lik_posterior_setup(
    fit                  = fit,
    posterior_samples    = posterior_samples,
    data                 = data,
    priors               = priors,
    unit                 = unit,
    data_hash            = data_hash,
    honor_bias_indicator = honor_bias_indicator
  )

  if (unit == "estimate") {
    return(.log_lik_estimate_sum_from_setup(setup))
  }

  return(.log_lik_cluster_sum_from_setup(setup))
}



.log_lik_evaluated_setup <- function(fit, data, priors, unit, data_hash,
                                     mu_samples, tau_within_samples,
                                     tau_between_samples, posterior_samples,
                                     honor_bias_indicator) {

  is_multilevel <- .is_data_multilevel(data)
  K             <- nrow(data[["outcome"]])

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

  return(list(
    fit                  = fit,
    priors               = priors,
    data                 = data,
    yi                   = .outcome_data_yi(active_object),
    sei                  = .outcome_data_sei(active_object),
    K                    = K,
    S                    = nrow(mu_samples),
    mu                   = mu_samples,
    tau_within           = tau_within_samples,
    tau_between          = tau_between_samples,
    cluster              = if (is_multilevel) split(seq_len(K), data[["outcome"]][["cluster"]]) else NULL,
    weights              = if (.is_data_weights(data)) data[["outcome"]][["weights"]] else NULL,
    data_hash            = data_hash,
    is_weightfunction    = .is_priors_weightfunction(priors),
    outcome_type         = .data_outcome_type(data),
    effect_direction     = .data_effect_direction(data),
    posterior_samples    = posterior_samples,
    honor_bias_indicator = honor_bias_indicator
  ))
}



.log_lik_posterior_setup <- function(fit, posterior_samples, data, priors,
                                     unit, data_hash,
                                     honor_bias_indicator) {

  is_mods           <- .is_data_mods(data)
  is_scale          <- .is_data_scale(data)
  is_multilevel     <- .is_data_multilevel(data)
  is_PET            <- .is_priors_PET(priors)
  is_PEESE          <- .is_priors_PEESE(priors)
  is_weightfunction <- .is_priors_weightfunction(priors)
  outcome_type      <- .data_outcome_type(data)
  effect_direction  <- .data_effect_direction(data)
  K                 <- nrow(data[["outcome"]])

  if (unit == "cluster" && !is_multilevel) {
    stop("Cluster-unit log-likelihood is only available for multilevel models.",
         call. = FALSE)
  }

  mu_samples <- .evaluate.brma.mu(
    fit               = fit,
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = FALSE,
    K                 = K,
    posterior_samples = posterior_samples
  )

  tau_result <- .evaluate.brma.tau(
    fit               = fit,
    scale_data        = data[["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
    scale_priors      = priors[["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples
  )

  cluster <- if (is_multilevel) {
    split(seq_len(K), data[["outcome"]][["cluster"]])
  } else {
    NULL
  }

  if (unit == "estimate" && is_multilevel) {
    mu_samples <- mu_samples + .evaluate.brma.cluster_effects(
      fit               = fit,
      tau_between       = tau_result[["tau_between"]],
      cluster           = data[["outcome"]][["cluster"]],
      same_data         = TRUE,
      effect_direction  = effect_direction,
      posterior_samples = posterior_samples
    )
  }

  active_object <- list(
    fit    = fit,
    data   = data,
    priors = priors
  )

  return(list(
    fit                  = fit,
    priors               = priors,
    data                 = data,
    yi                   = .outcome_data_yi(active_object),
    sei                  = .outcome_data_sei(active_object),
    K                    = K,
    S                    = nrow(mu_samples),
    mu                   = mu_samples,
    tau_within           = tau_result[["tau_within"]],
    tau_between          = tau_result[["tau_between"]],
    cluster              = cluster,
    weights              = if (.is_data_weights(data)) data[["outcome"]][["weights"]] else NULL,
    data_hash            = data_hash,
    is_weightfunction    = is_weightfunction,
    outcome_type         = outcome_type,
    effect_direction     = effect_direction,
    posterior_samples    = posterior_samples,
    honor_bias_indicator = honor_bias_indicator
  ))
}



.log_lik_estimate_from_setup <- function(setup) {

  data                <- setup[["data"]]
  priors              <- setup[["priors"]]
  yi                  <- setup[["yi"]]
  sei                 <- setup[["sei"]]
  K                   <- setup[["K"]]
  mu_samples          <- setup[["mu"]]
  tau_within_samples  <- setup[["tau_within"]]
  is_weightfunction   <- setup[["is_weightfunction"]]
  outcome_type        <- setup[["outcome_type"]]
  effect_direction    <- setup[["effect_direction"]]
  data_weights        <- setup[["weights"]]

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
        effect_direction     = effect_direction,
        honor_bias_indicator = setup[["honor_bias_indicator"]]
      )

      log_lik <- .outcome_pdf.selnorm(
        yi                = yi,
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = sei,
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
  mu_samples          <- setup[["mu"]]
  tau_within_samples  <- setup[["tau_within"]]
  is_weightfunction   <- setup[["is_weightfunction"]]
  outcome_type        <- setup[["outcome_type"]]
  effect_direction    <- setup[["effect_direction"]]
  data_weights        <- setup[["weights"]]

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
        effect_direction     = effect_direction,
        honor_bias_indicator = setup[["honor_bias_indicator"]]
      )

      return(.outcome_pdf_sum.selnorm(
        yi                = yi,
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = sei,
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
.estimate_likelihood_setup.brma <- function(object) {

  priors            <- object[["priors"]]
  data              <- object[["data"]]
  is_multilevel     <- .is_multilevel(object)
  is_mods           <- .is_mods(object)
  is_scale          <- .is_scale(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)
  posterior_samples <- .get_posterior_samples(object[["fit"]])

  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  K   <- length(yi)

  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = data[["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
    scale_priors      = priors[["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples
  )

  mu_samples <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = FALSE,
    K                 = K,
    posterior_samples = posterior_samples
  )

  fit_data <- NULL
  if (is_multilevel || is_weightfunction) {
    fit_data <- .create_fit_data(data = data, priors = priors)
  }

  if (is_multilevel) {
    mu_samples <- mu_samples + .evaluate.brma.cluster_effects(
      fit               = object[["fit"]],
      tau_between       = tau_result[["tau_between"]],
      cluster           = fit_data[["cluster"]],
      same_data         = TRUE,
      effect_direction  = effect_direction,
      posterior_samples = posterior_samples
    )
  }

  return(list(
    priors            = priors,
    data              = data,
    fit_data          = fit_data,
    yi                = yi,
    sei               = sei,
    K                 = K,
    S                 = nrow(mu_samples),
    mu                = mu_samples,
    tau_within        = tau_result[["tau_within"]],
    tau_between       = tau_result[["tau_between"]],
    is_weightfunction = is_weightfunction,
    outcome_type      = outcome_type,
    effect_direction  = effect_direction,
    posterior_samples = posterior_samples
  ))
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
    data_hash            = .get_outcome_hash(object),
    honor_bias_indicator = TRUE
  ))
}
