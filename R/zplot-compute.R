# ============================================================================ #
# Internal Helper Functions
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# .zplot_fun.brma
# ---------------------------------------------------------------------------- #
#
# Core computation function for zplot estimates and densities.
#
# This function operates in two modes:
# 1. EDR mode (z_threshold set): Computes Expected Discovery Rate
# 2. Density mode (z_sequence set): Computes density over z-values
#
# The extrapolate parameter controls bias adjustment:
# - extrapolate=TRUE:  Removes publication bias (PET/PEESE/weights) to estimate
#                      "true" power distribution
# - extrapolate=FALSE: Includes all bias adjustments for fitted distribution
#
# @param object      brma object with fit and priors
# @param z_threshold z-value threshold for significance (EDR mode)
# @param z_sequence  vector of z-values for density evaluation (density mode)
# @param max_samples maximum posterior samples for computation
# @param extrapolate whether to remove bias adjustments
#
# @return EDR mode:     list with EDR (S-vector) and weights (S-vector)
#         Density mode: S x length(z_sequence) matrix of densities
#
.zplot_fun.brma <- function(object, z_threshold = NULL, z_sequence = NULL,
                             max_samples = 10000, extrapolate = FALSE) {

  max_samples <- .normalize_max_samples(max_samples, "max_samples")

  ### extract model info
  is_multilevel     <- .is_multilevel(object)
  is_weightfunction <- .is_weightfunction(object)
  effect_direction  <- .effect_direction(object)

  ### get outcome data
  outcome_data <- object[["data"]][["outcome"]]
  yi           <- outcome_data[["yi"]]
  sei          <- outcome_data[["sei"]]
  K            <- length(yi)

  ### determine mu type
  # "cluster" for multilevel (mu + gamma), "terms" for marginal/single-level
  mu_type <- if (is_multilevel) "cluster" else "terms"

  ### 1. Thin posterior samples before derived quantities
  posterior_samples <- .get_posterior_samples(object[["fit"]])
  selected_ind      <- .thin_sample_rows(nrow(posterior_samples), max_samples)
  if (!is.null(selected_ind)) {
    posterior_samples <- posterior_samples[selected_ind, , drop = FALSE]
  }

  ### 2. Predict mu (location) samples
  # predict.brma handles PET/PEESE via bias_adjusted argument
  # extrapolate=TRUE -> bias_adjusted=TRUE (unbiased predictions)
  mu_samples <- predict.brma(
    object             = object,
    type               = mu_type,
    bias_adjusted      = extrapolate,
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  ### 3. Get tau (heterogeneity) samples
  # specific handling to get tau_within
  # .evaluate.brma.tau handles ML splitting logic
  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (.is_scale(object)) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = .is_scale(object),
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples
  )
  tau_within <- tau_result[["tau_within"]]

  ### 4. Prepare for Computation
  selection <- .zplot_selection_context(
    object            = object,
    posterior_samples = posterior_samples,
    is_weightfunction = is_weightfunction
  )


  ### 5. Dispatch: EDR (Threshold) vs Density (Sequence)

  if (!is.null(z_threshold)) {
    return(.zplot_threshold_vectorized(
      z_threshold      = z_threshold,
      mu_samples       = mu_samples,
      tau_within       = tau_within,
      sei              = sei,
      selection        = selection,
      extrapolate      = extrapolate,
      effect_direction = effect_direction
    ))
  }

  if (!is.null(z_sequence)) {
    return(.zplot_density_vectorized(
      z_sequence       = z_sequence,
      mu_samples       = mu_samples,
      tau_within       = tau_within,
      sei              = sei,
      selection        = selection,
      extrapolate      = extrapolate,
      effect_direction = effect_direction
    ))
  }

  return(NULL)
}


# ---------------------------------------------------------------------------- #
# .zplot_threshold_vectorized
# ---------------------------------------------------------------------------- #
#
# Vectorized EDR computation for normal and selected-normal rows.
#
# ---------------------------------------------------------------------------- #
.zplot_threshold_vectorized <- function(z_threshold, mu_samples, tau_within,
                                         sei, selection, extrapolate,
                                         effect_direction) {

  S         <- nrow(mu_samples)
  K         <- ncol(mu_samples)
  sei_mat   <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd  <- sqrt(tau_within^2 + sei_mat^2)

  if (!is.null(selection) && .has_native_zplot_threshold()) {
    return(.zplot_selnorm_threshold_summary(
      z_threshold       = z_threshold,
      mean              = mu_samples,
      sd                = total_sd,
      sei               = sei,
      selection_context = selection,
      extrapolate       = extrapolate
    ))
  }

  q_upper   <- z_threshold * sei
  q_lower   <- -z_threshold * sei

  thresholds <- stats::pnorm(
    matrix(q_upper, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu_samples,
    sd         = total_sd,
    lower.tail = FALSE
  ) + stats::pnorm(
    matrix(q_lower, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu_samples,
    sd         = total_sd,
    lower.tail = TRUE
  )
  weights       <- .zplot_inverse_selection_weights(
    mean      = mu_samples,
    sd        = total_sd,
    sei       = sei,
    selection = selection
  )
  weighted_rows <- if (is.null(selection)) integer(0) else which(!selection[["use_normal"]])
  if (length(weighted_rows) > 0) {
    selection_weight <- BayesTools::selection_context_subset_rows(
      context = selection,
      rows    = weighted_rows
    )
    mean_weight      <- mu_samples[weighted_rows, , drop = FALSE]
    sd_weight        <- total_sd[weighted_rows, , drop = FALSE]

    if (!extrapolate) {
      prob_upper <- .selection_step_cdf_matrix(
        q                 = q_upper,
        mean              = mean_weight,
        sd                = sd_weight,
        sei               = sei,
        selection_context = selection_weight,
        lower.tail        = FALSE
      )
      prob_lower <- .selection_step_cdf_matrix(
        q                 = q_lower,
        mean              = mean_weight,
        sd                = sd_weight,
        sei               = sei,
        selection_context = selection_weight,
        lower.tail        = TRUE
      )

      thresholds[weighted_rows, ] <- prob_upper + prob_lower
    }
  }

  if (!is.null(selection) && extrapolate) {
    EDR <- rowSums(thresholds * weights) / rowSums(weights)
  } else {
    EDR <- rowMeans(thresholds)
  }

  return(list(
    EDR     = EDR,
    weights = rowMeans(weights)
  ))
}


# ---------------------------------------------------------------------------- #
# .zplot_selnorm_threshold_summary
# ---------------------------------------------------------------------------- #
#
# Native EDR and inverse-selection-weight reductions for zplot thresholds.
#
# ---------------------------------------------------------------------------- #
.has_native_zplot_threshold <- function() {

  return(is.loaded(
    "RoBMA_selnorm_zcurve_threshold_summary",
    PACKAGE = "RoBMA"
  ))
}

.zplot_selnorm_threshold_summary <- function(z_threshold, mean, sd, sei,
                                             selection_context, extrapolate) {

  .selection_require_step_evaluable(selection_context, ".zplot_threshold_vectorized()")
  native_static <- BayesTools::selection_native_static_args(selection_context)

  return(.Call(
    "RoBMA_selnorm_zcurve_threshold_summary",
    .native_numeric_vector(z_threshold),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    as.logical(extrapolate),
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}


# ---------------------------------------------------------------------------- #
# .zplot_density_vectorized
# ---------------------------------------------------------------------------- #
#
# Vectorized z-density computation for normal and selected-normal rows.
#
# ---------------------------------------------------------------------------- #
.has_native_zplot_density <- function(selection = FALSE) {

  symbols <- "RoBMA_zcurve_normal_density_matrix"
  if (selection) {
    symbols <- c(symbols, "RoBMA_selnorm_zcurve_density_matrix")
  }

  return(all(vapply(symbols, is.loaded, logical(1), PACKAGE = "RoBMA")))
}

.zplot_normal_density_matrix <- function(z_sequence, mean, sd, sei) {

  if (!.has_native_zplot_density(selection = FALSE)) {
    stop("The native zplot density kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_zcurve_normal_density_matrix",
    .native_numeric_vector(z_sequence),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    PACKAGE = "RoBMA"
  ))
}

.zplot_selnorm_density_matrix <- function(z_sequence, mean, sd, sei,
                                           selection_context, extrapolate) {

  .selection_require_step_evaluable(selection_context, ".zplot_density_vectorized()")

  if (!.has_native_zplot_density(selection = TRUE)) {
    stop("The native selected-normal zplot density kernel is not loaded.", call. = FALSE)
  }
  native_static <- BayesTools::selection_native_static_args(selection_context)

  return(.Call(
    "RoBMA_selnorm_zcurve_density_matrix",
    .native_numeric_vector(z_sequence),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    as.logical(extrapolate),
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.zplot_density_vectorized <- function(z_sequence, mu_samples, tau_within,
                                       sei, selection, extrapolate,
                                       effect_direction) {

  S           <- nrow(mu_samples)
  K           <- ncol(mu_samples)
  sei_mat     <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd    <- sqrt(tau_within^2 + sei_mat^2)

  if (is.null(selection)) {
    return(.zplot_normal_density_matrix(
      z_sequence = z_sequence,
      mean       = mu_samples,
      sd         = total_sd,
      sei        = sei
    ))
  }

  density <- .zplot_selnorm_density_matrix(
    z_sequence        = z_sequence,
    mean              = mu_samples,
    sd                = total_sd,
    sei               = sei,
    selection_context = selection,
    extrapolate       = extrapolate
  )

  return(density)
}


# ---------------------------------------------------------------------------- #
# .zplot_inverse_selection_weights
# ---------------------------------------------------------------------------- #
#
# Expected attempted studies represented by each observed study under the
# selection model. Normal/no-bias branches have weight one.
#
# ---------------------------------------------------------------------------- #
.zplot_inverse_selection_weights <- function(mean, sd, sei, selection) {

  weights <- matrix(1, nrow = nrow(mean), ncol = ncol(mean))

  if (is.null(selection)) {
    return(weights)
  }

  weighted_rows <- which(!selection[["use_normal"]])
  if (length(weighted_rows) == 0L) {
    return(weights)
  }

  selection_weight <- BayesTools::selection_context_subset_rows(
    context = selection,
    rows    = weighted_rows
  )
  log_norm <- .selection_step_log_norm_matrix(
    mean              = mean[weighted_rows, , drop = FALSE],
    sd                = sd[weighted_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_weight
  )

  weights[weighted_rows, ] <- exp(-log_norm)
  return(weights)
}


# ---------------------------------------------------------------------------- #
# .zplot_weighted_rows / .zplot_constant_rows
# ---------------------------------------------------------------------------- #
#
# Row selectors for fitted and extrapolated selection-model zplot densities.
#
# ---------------------------------------------------------------------------- #
.zplot_weighted_rows <- function(selection, extrapolate) {

  if (is.null(selection) || extrapolate) {
    return(integer(0))
  }

  return(which(!selection[["use_normal"]]))
}

.zplot_constant_rows <- function(selection, extrapolate) {

  return(integer(0))
}


# ---------------------------------------------------------------------------- #
# .zplot_selection_context
# ---------------------------------------------------------------------------- #
#
# Prepare posterior-row selection metadata for branch-aware zplot evaluation.
#
# ---------------------------------------------------------------------------- #
.zplot_selection_context <- function(object, posterior_samples,
                                     is_weightfunction) {

  if (!is_weightfunction) {
    return(NULL)
  }

  return(.selection_context(
    object            = object,
    posterior_samples = posterior_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .zplot_selection_args
# ---------------------------------------------------------------------------- #
#
# Extract active omega and cutpoints for one posterior row and estimate.
#
# ---------------------------------------------------------------------------- #
.zplot_selection_args <- function(selection, row, estimate, n = 1L) {

  .selection_require_step_evaluable(selection, ".zplot_selection_args()")

  omega <- selection[["omega"]][row, , drop = FALSE]

  if (n > 1L) {
    omega <- matrix(as.numeric(omega), nrow = n, ncol = ncol(omega), byrow = TRUE)
  }

  return(list(
    omega   = omega,
    crit_yi = stats::qnorm(rev(selection[["p_cuts"]])[-c(1L, length(selection[["p_cuts"]]))],
                           lower.tail = FALSE) *
      selection[["sei"]][estimate]
  ))
}


# ============================================================================ #
# Graphical Helper Functions
# ============================================================================ #
#
# These functions extract and set default graphical parameters for zplot
# plotting components. They handle the translation between base graphics
# and ggplot2 parameter naming conventions.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .get_dots_hist_zplot
# ---------------------------------------------------------------------------- #
#
# Extracts histogram graphical parameters with defaults.
#
# @param dots      list of user-supplied graphical parameters
# @param plot_type "base" or "ggplot"
# @param max_density maximum density value for setting ylim
#
# @return list of graphical parameters appropriate for plot_type
#
