# ============================================================================ #
# funnel-quantiles.R
# ============================================================================ #

# .get_funnel_tau
# ---------------------------------------------------------------------------- #
#
# Extract mean tau (heterogeneity) estimate from brma object for funnel plot.
#
# For multilevel models, returns total tau (combining within and between components).
#
# @param object       brma object
# @return numeric scalar representing the mean tau estimate
#
# ---------------------------------------------------------------------------- #
.get_funnel_tau <- function(object) {

  if (.is_data_known_v(object[["data"]])) {
    return(.get_funnel_tau_known_v(object))
  }

  # use pooled_heterogeneity to get mean tau
  tau_samples <- pooled_heterogeneity(object)
  tau_summary <- summary(tau_samples)
  return(tau_summary["tau", "Mean"])
}


.get_funnel_tau_known_v <- function(object) {

  extra_variance <- .known_v_extra_variance_samples(object)
  tau_samples    <- sqrt(rowMeans(extra_variance))

  return(mean(tau_samples))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles
# ---------------------------------------------------------------------------- #
#
# Compute quantiles for funnel plot contours based on the sampling distribution.
#
# When sampling_bias = TRUE, posterior rows dispatch to their active bias model
# and the contour is obtained by inverting the model-averaged CDF.
#
# When sampling_bias = FALSE: uses standard normal quantiles.
#
# @param x                      brma object
# @param se_sequence            numeric vector of SE values for funnel contours
# @param mu                     numeric scalar of pooled effect estimate
# @param tau                    numeric scalar of heterogeneity estimate
# @param sampling_heterogeneity logical; whether to incorporate tau
# @param sampling_bias          logical; whether to incorporate bias into sampling dist
# @param is_weightfunction      logical; whether model has selection model
# @param is_PET                 logical; whether model has PET adjustment
# @param is_PEESE               logical; whether model has PEESE adjustment
# @param effect_direction       character; "positive" or "negative"
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles <- function(x, se_sequence, mu, tau,
                                  sampling_heterogeneity, sampling_bias,
                                  is_weightfunction, is_PET, is_PEESE,
                                  effect_direction, max_samples) {
  n_se   <- length(se_sequence)
  sd_seq <- .root_sum_squares(se_sequence, tau)

  # default: standard normal quantiles centered at mu
  if (!sampling_bias || (!is_weightfunction && !is_PET && !is_PEESE)) {
    return(list(
      lower = stats::qnorm(0.025, mean = mu, sd = sd_seq),
      upper = stats::qnorm(0.975, mean = mu, sd = sd_seq),
      mid   = rep(mu, n_se)
    ))
  }

  if (!BayesTools::is.prior.mixture(x[["priors"]][["outcome"]][["bias"]])) {
    if (is_PET || is_PEESE) {
      return(.get_funnel_quantiles_PETPEESE_plugin(
        x                = x,
        se_sequence      = se_sequence,
        mu               = mu,
        sd_seq           = sd_seq,
        is_PET           = is_PET,
        is_PEESE         = is_PEESE,
        effect_direction = effect_direction
      ))
    }
  }

  return(.get_funnel_quantiles_model_averaged(
    x                      = x,
    se_sequence            = se_sequence,
    sampling_heterogeneity = sampling_heterogeneity,
    effect_direction       = effect_direction,
    max_samples            = max_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles_PETPEESE_plugin
# ---------------------------------------------------------------------------- #
#
# Fast plug-in normal contours for single PET/PEESE priors.
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_PETPEESE_plugin <- function(x, se_sequence, mu, sd_seq,
                                                  is_PET, is_PEESE,
                                                  effect_direction) {

  n_se              <- length(se_sequence)
  posterior_samples <- .get_posterior_samples(x[["fit"]])
  bias_shift        <- rep(0, n_se)
  direction         <- ifelse(effect_direction == "negative", -1, 1)

  if (is_PET && "PET" %in% colnames(posterior_samples)) {
    PET_mean   <- mean(posterior_samples[, "PET"])
    bias_shift <- bias_shift + direction * PET_mean * se_sequence
  }

  if (is_PEESE && "PEESE" %in% colnames(posterior_samples)) {
    PEESE_mean <- mean(posterior_samples[, "PEESE"])
    bias_shift <- bias_shift + direction * PEESE_mean * se_sequence^2
  }

  mid   <- mu + bias_shift
  lower <- stats::qnorm(0.025, mean = mid, sd = sd_seq)
  upper <- stats::qnorm(0.975, mean = mid, sd = sd_seq)

  return(list(lower = lower, upper = upper, mid = mid))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles_model_averaged
# ---------------------------------------------------------------------------- #
#
# Compute model-averaged funnel contours for publication-bias models.
#
# Each posterior row dispatches to its active bias model:
# - no-bias rows use a normal CDF
# - PET/PEESE rows use a normal CDF with the row-specific bias offset
# - selection-model rows use the selected-normal CDF
#
# Quantiles are obtained from the averaged posterior-row CDF.
#
# @return list with 'lower', 'upper', and 'mid' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_model_averaged <- function(x, se_sequence,
                                                 sampling_heterogeneity,
                                                 effect_direction,
                                                 max_samples = Inf) {

  setup <- .funnel_model_averaged_setup(
    x                      = x,
    sampling_heterogeneity = sampling_heterogeneity,
    max_samples            = max_samples
  )

  if (.has_native_funnel_model_averaged_quantiles(setup)) {
    return(.funnel_model_averaged_quantiles_native(
      se_sequence      = se_sequence,
      setup            = setup,
      effect_direction = effect_direction
    ))
  }

  lower <- vapply(
    se_sequence,
    .funnel_model_averaged_quantile,
    numeric(1),
    p                = 0.025,
    setup            = setup,
    effect_direction = effect_direction
  )
  upper <- vapply(
    se_sequence,
    .funnel_model_averaged_quantile,
    numeric(1),
    p                = 0.975,
    setup            = setup,
    effect_direction = effect_direction
  )
  mid <- vapply(
    se_sequence,
    .funnel_model_averaged_quantile,
    numeric(1),
    p                = 0.5,
    setup            = setup,
    effect_direction = effect_direction
  )

  return(list(lower = lower, upper = upper, mid = mid))
}

.has_native_funnel_model_averaged_quantiles <- function(setup) {

  if (is.null(setup[["selection"]])) {
    return(FALSE)
  }

  return(is.loaded("RoBMA_funnel_model_averaged_quantiles", PACKAGE = "RoBMA"))
}

.funnel_model_averaged_quantiles_native <- function(se_sequence, setup,
                                                    effect_direction) {

  selection <- setup[["selection"]]
  if (any(setup[["is_weightfunction"]])) {
    .selection_require_step_evaluable(
      selection,
      ".funnel_model_averaged_quantiles()"
    )
  }

  direction <- ifelse(effect_direction == "negative", -1L, 1L)
  native_static <- BayesTools::selection_native_static_args(selection)

  return(.Call(
    "RoBMA_funnel_model_averaged_quantiles",
    .native_numeric_vector(se_sequence),
    .native_numeric_vector(setup[["mu"]]),
    .native_numeric_vector(setup[["tau"]]),
    .native_numeric_vector(setup[["PET"]]),
    .native_numeric_vector(setup[["PEESE"]]),
    .native_integer_vector(setup[["is_weightfunction"]]),
    .native_numeric_matrix(selection[["omega"]]),
    .native_numeric_vector(selection[["alpha"]]),
    .native_integer_vector(selection[["phack_kind"]]),
    .native_integer_vector(selection[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    .native_integer_vector(direction),
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_model_averaged_setup
# ---------------------------------------------------------------------------- #
#
# Prepare posterior-row quantities for model-averaged funnel contours.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_setup <- function(x, sampling_heterogeneity,
                                         max_samples = Inf) {

  posterior_samples <- .get_posterior_samples(x[["fit"]])
  priors_bias       <- x[["priors"]][["outcome"]][["bias"]]

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

  selected_rows <- .funnel_subsample_rows(
    bias_indicator = bias_indicator,
    max_samples    = max_samples
  )
  if (!is.null(selected_rows)) {
    posterior_samples <- posterior_samples[selected_rows, , drop = FALSE]
    bias_indicator    <- bias_indicator[selected_rows]
  }
  S <- nrow(posterior_samples)

  mu_samples <- .funnel_mu_samples(x, posterior_samples)

  if (sampling_heterogeneity) {
    tau_samples <- .funnel_tau_samples(x, posterior_samples)
  } else {
    tau_samples <- rep(0, S)
  }

  PET_samples   <- .funnel_posterior_column(posterior_samples, "PET",   S)
  PEESE_samples <- .funnel_posterior_column(posterior_samples, "PEESE", S)

  PET_samples[!branch_is_PET[bias_indicator]]     <- 0
  PEESE_samples[!branch_is_PEESE[bias_indicator]] <- 0

  selection <- .selection_context(
    object            = x,
    posterior_samples = posterior_samples
  )
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
    selection             = selection
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_subsample_rows
# ---------------------------------------------------------------------------- #
#
# Deterministically thin posterior rows for plotting while preserving bias-model
# proportions. Model-averaged funnel contours are visual summaries, and using
# every row from large MCMC fits makes root finding unnecessarily expensive.
#
# ---------------------------------------------------------------------------- #
.funnel_subsample_rows <- function(bias_indicator, max_samples) {

  return(.thin_sample_rows_by_group(
    group       = bias_indicator,
    max_samples = max_samples
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

  if (!.is_mods(x)) {
    return(as.numeric(posterior_samples[, "mu"]))
  }

  mu_samples <- pooled_effect.brma(
    object             = x,
    bias_adjusted      = TRUE,
    .posterior_samples = posterior_samples
  )

  return(as.numeric(as.matrix(mu_samples)[, 1]))
}


# ---------------------------------------------------------------------------- #
# .funnel_tau_samples
# ---------------------------------------------------------------------------- #
#
# Extract total heterogeneity samples for outcome-mode funnel contours.
#
# ---------------------------------------------------------------------------- #
.funnel_tau_samples <- function(x, posterior_samples) {

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
