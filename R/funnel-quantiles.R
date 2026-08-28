# ============================================================================ #
# funnel-quantiles.R
# ============================================================================ #

# .get_radial_tau
# ---------------------------------------------------------------------------- #
#
# Extract mean scalar tau (heterogeneity) estimate for the radial plot.
#
# For multilevel models, returns total tau (combining within and between components).
#
# @param object       brma object
# @return numeric scalar representing the mean tau estimate
#
# ---------------------------------------------------------------------------- #
.get_radial_tau <- function(object) {

  if (.is_data_known_v(object[["data"]])) {
    return(.get_radial_tau_known_v(object))
  }

  posterior_samples <- .get_posterior_samples(object[["fit"]])
  tau_samples <- .radial_heterogeneity_samples(
    object            = object,
    posterior_samples = posterior_samples
  )

  return(mean(tau_samples))
}


.get_radial_tau_known_v <- function(object) {

  extra_variance <- .known_v_extra_variance_samples(object)
  tau_samples    <- sqrt(rowMeans(extra_variance))

  return(mean(tau_samples))
}


.radial_heterogeneity_samples <- function(object, posterior_samples) {

  if (inherits(object, "brma.mv")) {
    components <- .brma_mv_heterogeneity_components(
      object                         = object,
      posterior_samples              = posterior_samples,
      include_known_group_covariance = TRUE
    )
    total <- .total_brma_mv_heterogeneity_samples(components)
    return(.brma_mv_rms_sd_samples(total))
  }

  return(.radial_tau_samples(object, posterior_samples))
}


# Determine whether the fitted model has one row-invariant marginal
# heterogeneity distribution. For random-formula brma.mv models this evaluates
# the same row-marginal ZGZ' variances used by prediction and bridge likelihoods.
.funnel_common_heterogeneity <- function(object, posterior_samples = NULL) {

  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    posterior_samples
  )
  tau <- .funnel_row_heterogeneity_samples(object, posterior_samples)

  reference <- matrix(
    tau[, 1L],
    nrow = nrow(tau),
    ncol = ncol(tau)
  )
  common <- all(tau == reference)

  return(list(
    common            = common,
    posterior_samples = posterior_samples,
    tau                = if (common) tau[, 1L] else NULL
  ))
}


.funnel_row_heterogeneity_samples <- function(object, posterior_samples) {

  if (.is_data_known_v(object[["data"]])) {
    return(sqrt(.known_v_extra_variance_samples(
      object            = object,
      posterior_samples = posterior_samples
    )))
  }

  if (inherits(object, "brma.mv")) {
    components <- .brma_mv_heterogeneity_components(
      object                         = object,
      posterior_samples              = posterior_samples,
      include_known_group_covariance = TRUE
    )
    return(.total_brma_mv_heterogeneity_samples(components))
  }

  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (.is_scale(object)) {
      .create_fit_formula_list(data = object[["data"]], "scale")
    } else {
      NULL
    },
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = .is_scale(object),
    is_multilevel     = .is_multilevel(object),
    K                 = nrow(object[["data"]][["outcome"]]),
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )

  return(.expand_brma_mv_heterogeneity_samples(
    samples = tau_result[["tau_total"]],
    S       = nrow(posterior_samples),
    K       = nrow(object[["data"]][["outcome"]])
  ))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles
# ---------------------------------------------------------------------------- #
#
# Compute quantiles for funnel plot contours based on the sampling distribution.
#
# Plug-in contours use one conditional posterior-mean row per complete joint
# model. Posterior-predictive contours use posterior rows directly. In both
# cases the active bias model supplies the row CDF and the mixture is inverted.
#
# @param x                      brma object
# @param se_sequence            numeric vector of SE values for funnel contours
# @param sampling_heterogeneity logical; whether to incorporate tau
# @param sampling_bias          logical; whether to incorporate bias into sampling dist
# @param effect_direction       character; "positive" or "negative"
# @param max_samples             posterior-row budget for Bayesian contours
# @param estimand                plug-in or posterior-predictive construction
# @param common_heterogeneity    validated common-heterogeneity setup
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles <- function(x, se_sequence,
                                  sampling_heterogeneity, sampling_bias,
                                  effect_direction, max_samples, estimand,
                                  common_heterogeneity) {

  setup <- .funnel_sampling_setup(
    x                      = x,
    sampling_heterogeneity = sampling_heterogeneity,
    sampling_bias          = sampling_bias,
    max_samples            = max_samples,
    estimand               = estimand,
    common_heterogeneity   = common_heterogeneity
  )

  return(.get_funnel_quantiles_from_setup(
    se_sequence      = se_sequence,
    setup            = setup,
    effect_direction = effect_direction
  ))
}


.get_funnel_quantiles_from_setup <- function(se_sequence, setup,
                                             effect_direction) {

  if (length(setup[["mu"]]) == 1L &&
      !any(setup[["is_weightfunction"]])) {
    location <- vapply(
      se_sequence,
      .funnel_row_location,
      numeric(1),
      setup            = setup,
      effect_direction = effect_direction
    )
    total_sd <- .root_sum_squares(se_sequence, setup[["tau"]])
    lower    <- stats::qnorm(0.025, mean = location, sd = total_sd)
    upper    <- stats::qnorm(0.975, mean = location, sd = total_sd)
    lower[total_sd == 0] <- location[total_sd == 0]
    upper[total_sd == 0] <- location[total_sd == 0]

    return(list(lower = lower, upper = upper, mid = location))
  }

  return(.get_funnel_quantiles_model_averaged(
    se_sequence      = se_sequence,
    setup            = setup,
    effect_direction = effect_direction
  ))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles_model_averaged
# ---------------------------------------------------------------------------- #
#
# Compute mixture funnel contours.
#
# Each plug-in model row or posterior row dispatches to its active bias model:
# - no-bias rows use a normal CDF
# - PET/PEESE rows use a normal CDF with the row-specific bias offset
# - selection-model rows use the selected-normal CDF
#
# Quantiles are obtained from the weighted row CDF.
#
# @return list with 'lower', 'upper', and 'mid' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_model_averaged <- function(se_sequence, setup,
                                                 effect_direction) {

  return(.funnel_model_averaged_quantiles_native(
    se_sequence      = se_sequence,
    setup            = setup,
    effect_direction = effect_direction
  ))
}

.funnel_model_averaged_quantiles_native <- function(se_sequence, setup,
                                                    effect_direction) {

  S <- length(setup[["mu"]])
  L <- length(se_sequence)
  location <- vapply(
    se_sequence,
    .funnel_row_location,
    numeric(S),
    setup            = setup,
    effect_direction = effect_direction
  )
  location <- matrix(location, nrow = S, ncol = L)
  total_sd <- vapply(
    se_sequence,
    function(se) .root_sum_squares(setup[["tau"]], se),
    numeric(S)
  )
  total_sd <- matrix(total_sd, nrow = S, ncol = L)

  quantiles <- .plot_mixture_quantiles_native(
    mean_samples      = location,
    sd_samples        = total_sd,
    se                = se_sequence,
    probs             = c(0.025, 0.975, 0.5),
    weights           = .funnel_setup_weights(setup),
    selected_rows     = setup[["is_weightfunction"]],
    selection_context = setup[["selection"]],
    caller            = ".funnel_model_averaged_quantiles()"
  )

  return(list(
    lower = quantiles[, 1L],
    upper = quantiles[, 2L],
    mid   = quantiles[, 3L]
  ))
}


.funnel_setup_weights <- function(setup) {

  weights <- setup[["weights"]]
  if (is.null(weights)) {
    return(rep(1 / length(setup[["mu"]]), length(setup[["mu"]])))
  }

  return(weights / sum(weights))
}


.funnel_row_location <- function(se, setup, effect_direction) {

  direction <- ifelse(effect_direction == "negative", -1, 1)

  return(
    setup[["mu"]] +
      direction * setup[["PET"]] * se +
      direction * setup[["PEESE"]] * se^2
  )
}


# ---------------------------------------------------------------------------- #
# .funnel_sampling_setup
# ---------------------------------------------------------------------------- #
#
# Prepare plug-in model rows or posterior rows for funnel contours.
#
# ---------------------------------------------------------------------------- #
.funnel_sampling_setup <- function(x, sampling_heterogeneity, sampling_bias,
                                   max_samples, estimand,
                                   common_heterogeneity) {

  estimand          <- match.arg(estimand, c("plugin", "posterior_predictive"))
  posterior_samples <- common_heterogeneity[["posterior_samples"]]
  tau_samples       <- common_heterogeneity[["tau"]]

  if (estimand == "posterior_predictive") {
    selected_rows <- .nested_srs_rows(
      rows        = seq_len(nrow(posterior_samples)),
      max_samples = max_samples
    )
    posterior_samples <- posterior_samples[selected_rows, , drop = FALSE]
    tau_samples       <- tau_samples[selected_rows]

    return(.funnel_setup_from_samples(
      x                      = x,
      posterior_samples      = posterior_samples,
      tau_samples            = tau_samples,
      sampling_heterogeneity = sampling_heterogeneity,
      sampling_bias          = sampling_bias,
      weights                = NULL
    ))
  }

  groups     <- .funnel_joint_model_groups(x, posterior_samples)
  group_rows <- split(seq_len(nrow(posterior_samples)), groups)
  model_samples <- t(vapply(group_rows, function(rows) {

    colMeans(posterior_samples[rows, , drop = FALSE])
  }, numeric(ncol(posterior_samples))))
  colnames(model_samples) <- colnames(posterior_samples)

  model_tau <- vapply(group_rows, function(rows) {

    mean(tau_samples[rows])
  }, numeric(1))
  model_weights <- lengths(group_rows) / nrow(posterior_samples)

  return(.funnel_setup_from_samples(
    x                      = x,
    posterior_samples      = model_samples,
    tau_samples            = model_tau,
    sampling_heterogeneity = sampling_heterogeneity,
    sampling_bias          = sampling_bias,
    weights                = model_weights
  ))
}


.funnel_joint_model_groups <- function(x, posterior_samples) {

  components <- .summary_models_components(x)
  if (length(components) == 0L) {
    return(rep(1L, nrow(posterior_samples)))
  }

  indicators <- lapply(components, function(info) {

    .summary_models_indicators(
      posterior_samples = posterior_samples,
      parameter         = info[["parameter"]],
      prior             = info[["prior"]],
      column            = info[["indicator"]],
      offset            = info[["indicator_offset"]]
    )
  })
  keys <- do.call(paste, c(indicators, sep = ":"))

  return(match(keys, unique(keys)))
}


.funnel_setup_from_samples <- function(x, posterior_samples, tau_samples,
                                       sampling_heterogeneity, sampling_bias,
                                       weights) {

  priors_bias <- x[["priors"]][["outcome"]][["bias"]]

  if (!BayesTools::is.prior.mixture(priors_bias)) {
    priors_bias <- list(priors_bias)
  }

  branch_is_PET   <- vapply(priors_bias, BayesTools::is.prior.PET,   logical(1))
  branch_is_PEESE <- vapply(priors_bias, BayesTools::is.prior.PEESE, logical(1))
  bias_indicator  <- .extract_bias_indicator(x, posterior_samples = posterior_samples)

  if (any(is.na(bias_indicator)) ||
      any(bias_indicator < 1L | bias_indicator > length(priors_bias))) {
    stop("Invalid 'bias_indicator' values in posterior samples.", call. = FALSE)
  }

  S <- nrow(posterior_samples)

  mu_samples <- .funnel_mu_samples(x, posterior_samples)

  if (!sampling_heterogeneity) {
    tau_samples <- rep(0, S)
  }

  PET_samples   <- .funnel_posterior_column(posterior_samples, "PET",   S)
  PEESE_samples <- .funnel_posterior_column(posterior_samples, "PEESE", S)

  PET_samples[!branch_is_PET[bias_indicator]]     <- 0
  PEESE_samples[!branch_is_PEESE[bias_indicator]] <- 0

  if (!sampling_bias) {
    PET_samples[]   <- 0
    PEESE_samples[] <- 0
  }

  selection <- if (sampling_bias) {
    .selection_context(
      object            = x,
      posterior_samples = posterior_samples
    )
  } else {
    NULL
  }
  use_normal <- if (is.null(selection)) {
    rep(TRUE, S)
  } else {
    selection[["use_normal"]]
  }

  return(list(
    mu                    = mu_samples,
    tau                   = tau_samples,
    PET                   = PET_samples,
    PEESE                 = PEESE_samples,
    bias_indicator        = bias_indicator,
    is_weightfunction     = !use_normal,
    selection             = selection,
    weights               = weights
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_mu_samples
# ---------------------------------------------------------------------------- #
#
# Extract pooled location samples without PET/PEESE bias offsets.
#
# ---------------------------------------------------------------------------- #
.funnel_mu_samples <- function(x, posterior_samples) {

  data   <- x[["data"]]
  priors <- x[["priors"]]
  mu_samples <- .evaluate.brma.mu(
    fit               = x[["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (.is_mods(x)) {
      .create_fit_formula_list(data = data, "mods")
    } else {
      NULL
    },
    mods_priors       = if (.is_random(x)) {
      priors[["location"]]
    } else {
      priors[["mods"]]
    },
    priors            = priors,
    is_mods           = .is_mods(x),
    is_PET            = .is_PET(x),
    is_PEESE          = .is_PEESE(x),
    effect_direction  = .effect_direction(x),
    bias_adjusted     = TRUE,
    K                 = nrow(data[["outcome"]]),
    posterior_samples = posterior_samples
  )

  return(rowMeans(mu_samples))
}


# ---------------------------------------------------------------------------- #
# .radial_tau_samples
# ---------------------------------------------------------------------------- #
#
# Extract RMS total heterogeneity samples for the radial plot.
#
# ---------------------------------------------------------------------------- #
.radial_tau_samples <- function(x, posterior_samples) {

  tau_result <- .evaluate.brma.tau(
    fit               = x[["fit"]],
    scale_data        = x[["data"]][["scale"]],
    scale_formula     = if (.is_scale(x)) .create_fit_formula_list(data = x[["data"]], "scale") else NULL,
    scale_priors      = x[["priors"]][["scale"]],
    is_scale          = .is_scale(x),
    is_multilevel     = .is_multilevel(x),
    K                 = nrow(x[["data"]][["outcome"]]),
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(x[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(x[["priors"]])
  )

  return(.brma_mv_rms_sd_samples(tau_result[["tau_total"]]))
}


# ---------------------------------------------------------------------------- #
# .funnel_posterior_column
# ---------------------------------------------------------------------------- #
#
# Extract an optional posterior column, returning zeros when absent.
#
# ---------------------------------------------------------------------------- #
.funnel_posterior_column <- function(posterior_samples, column, S) {

  if (column %in% colnames(posterior_samples)) {
    return(as.numeric(posterior_samples[, column]))
  }

  return(rep(0, S))
}
