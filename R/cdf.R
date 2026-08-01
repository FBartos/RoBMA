# ============================================================================ #
# Outcome CDF Functions (Cumulative Distribution)
# ============================================================================ #
#
# These functions compute pointwise CDF values F(yi | theta) for each
# observation and posterior sample. The target-specific estimate-unit CDF is
# used for LOO-PIT residuals via probability integral transformation.
#
# Parallels the structure of pdf.R but returns CDF values instead of
# density values.
#
# Note: For binomial and Poisson models, each "observation" consists of a pair
# of data points (ai+ci or x1i+x2i) that together define a single effect size
# estimate. The CDF is computed using the implied normal approximation for the
# effect size (log-OR or log-IRR).
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .outcome_cdf.norm
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for normal outcome models.
#
# For normal outcome models, the observed effect y_i has distribution:
#   y_i ~ N(mu_i, tau_within^2 + se_i^2)
#
# This function computes F(y_i | mu_i, sigma_i) for each observation at each
# posterior sample.
#
# @param yi               numeric vector of length K; observed effect sizes
# @param mu_samples       S x K matrix of location samples (with cluster effects if multilevel)
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
# @param sei              numeric vector of length K; standard errors
# @param lower.tail       logical; return P(Y <= yi) if TRUE, P(Y > yi) if FALSE
#
# @return S x K matrix of CDF values in [0, 1]
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.norm <- function(yi, mu_samples, tau_within, sei,
                              lower.tail = TRUE) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate yi and sei across samples: matrix(vec, S, K, byrow = TRUE)
  yi_mat  <- matrix(yi,  nrow = S, ncol = K, byrow = TRUE)
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

  # compute total SD: sqrt(tau^2 + se^2)
  total_sd <- .root_sum_squares(tau_within, sei_mat)

  # compute CDF value for each cell
  cdf_vals <- stats::pnorm(
    yi_mat,
    mean       = mu_samples,
    sd         = total_sd,
    lower.tail = lower.tail
  )

  return(cdf_vals)
}


# ---------------------------------------------------------------------------- #
# .outcome_cdf.selnorm
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for the selected-normal kernel.
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.selnorm <- function(yi, mu_samples, tau_within, sei,
                                 selection_context, lower.tail = TRUE) {

  S         <- nrow(mu_samples)
  K         <- ncol(mu_samples)
  sei_mat   <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd  <- .root_sum_squares(tau_within, sei_mat)

  cdf_vals <- .selection_step_cdf_matrix(
    q                 = yi,
    mean              = mu_samples,
    sd                = total_sd,
    sei               = sei,
    selection_context = selection_context,
    lower.tail        = lower.tail
  )

  return(cdf_vals)
}


# ---------------------------------------------------------------------------- #
# .outcome_cdf.binom
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for binomial outcome models.
#
# For binomial outcome models, we use a normal approximation based on the
# implied log-odds ratio effect size and its approximate sampling variance.
# This approach is consistent with how metafor computes residuals for GLMM.
#
# The scalar residual is computed on the approximate log-odds-ratio scale:
#   y_i ~ N(mu_i, tau_within_i^2 + sigma_i^2)
# where sigma_i is the approximate sampling SE from the cell counts.
#
# @param yi               numeric vector of length K; approximate log-OR effect sizes
# @param sei              numeric vector of length K; approximate sampling SEs
# @param mu_samples       S x K matrix of log-odds ratio samples
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
#
# @return S x K matrix of CDF values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.binom <- function(yi, sei, mu_samples, tau_within) {
  return(.outcome_cdf.norm(
    yi         = yi,
    mu_samples = mu_samples,
    tau_within = tau_within,
    sei        = sei
  ))
}


# ---------------------------------------------------------------------------- #
# .outcome_cdf.pois
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for Poisson outcome models.
#
# For Poisson outcome models, we use a normal approximation based on the
# implied log incidence rate ratio effect size and its approximate sampling
# variance. This is consistent with how metafor computes residuals for GLMM.
#
# The scalar residual is computed on the approximate log-incidence-rate-ratio
# scale:
#   y_i ~ N(mu_i, tau_within_i^2 + sigma_i^2)
# where sigma_i is the approximate sampling SE from the counts.
#
# @param yi               numeric vector of length K; approximate log-IRR effect sizes
# @param sei              numeric vector of length K; approximate sampling SEs
# @param mu_samples       S x K matrix of log-IRR samples
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
#
# @return S x K matrix of CDF values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.pois <- function(yi, sei, mu_samples, tau_within) {
  return(.outcome_cdf.norm(
    yi         = yi,
    mu_samples = mu_samples,
    tau_within = tau_within,
    sei        = sei
  ))
}


# ---------------------------------------------------------------------------- #
# .cdf.brma
# ---------------------------------------------------------------------------- #
#
# Compute the full CDF matrix for a brma object.
#
# This function computes pointwise CDF values F(yi | mu, tau) for each
# observation and posterior sample. It uses predict.brma to obtain the
# appropriate mu samples at the requested conditioning depth.
#
# @param object brma object
# @param conditioning_depth character; conditioning depth. Options are:
#                           - "marginal" (default): Fixed effects only (mu).
#                             CDF: yi ~ N(mu_i, tau^2 + sei^2)
#                           - "cluster": Fixed effects + cluster-level random
#                             effects (mu + gamma). Only available for
#                             multilevel models. CDF:
#                             yi ~ N(mu_i + gamma_j, tau_within^2 + sei^2)
#                           - "estimate": True estimate effects
#                             (mu + gamma + theta). CDF: yi ~ N(theta_i, sei^2)
#
# @return S x K matrix of CDF values in [0, 1]
#
# ---------------------------------------------------------------------------- #
.cdf.brma <- function(object, conditioning_depth = "marginal") {

  ### input validation
  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)
  if (.is_data_known_v(object[["data"]])) {
    if (conditioning_depth == "estimate") {
      cdf_vals <- .cdf_lik_estimate.brma(object)
      colnames(cdf_vals) <- paste0("cdf[", seq_len(ncol(cdf_vals)), "]")
      return(cdf_vals)
    }
    stop(
      ".cdf.brma() with known-V data is available only with ",
      "conditioning_depth = 'estimate'.",
      call. = FALSE
    )
  }
  if (.is_random(object)) {
    .check_random_formula_postfit_deferred(object, ".cdf.brma()")
  }

  ### extract structural information about the model
  priors            <- object[["priors"]]
  data              <- object[["data"]]
  is_multilevel     <- .is_multilevel(object)
  is_scale          <- .is_scale(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)
  posterior_samples <- .get_posterior_samples(object[["fit"]])

  .check_unit_conditioning_depth(
    object             = object,
    unit               = "estimate",
    conditioning_depth = conditioning_depth,
    caller             = ".cdf.brma()"
  )

  ### obtain observed effect sizes and sampling SEs
  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  K   <- length(yi)

  ### get mu samples at the appropriate conditioning depth using predict.brma
  # map CDF conditioning depths to predict.brma types
  predict_type <- switch(conditioning_depth,
    "marginal"  = "terms",
    "cluster"   = "cluster",
    "estimate"  = "estimate"
  )

  mu_samples <- predict.brma(
    object  = object,
    newdata = NULL,
    type    = predict_type,
    quiet   = TRUE
  )
  S <- nrow(mu_samples)

  ### determine tau component for CDF computation
  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors      = priors[["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )

  if (conditioning_depth == "estimate") {
    tau_within_samples <- matrix(0, nrow = S, ncol = K)
  } else if (conditioning_depth == "cluster") {
    tau_within_samples <- tau_result[["tau_within"]]
  } else {
    tau_within_samples <- tau_result[["tau_total"]]
  }

  ### compute CDF based on outcome type
  if (outcome_type == "norm") {

    # flip for negative effect direction
    if (effect_direction == "negative") {
      mu_samples_cdf <- -mu_samples
      yi_cdf         <- -yi
      lower_tail     <- FALSE
    } else {
      mu_samples_cdf <- mu_samples
      yi_cdf         <- yi
      lower_tail     <- TRUE
    }

    # dispatch between weighted and standard normal
    if (is_weightfunction) {
      selection_context <- .selection_context(
        object            = object,
        posterior_samples = posterior_samples
      )

      cdf_vals <- .outcome_cdf.selnorm(
        yi                = yi,
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = sei,
        selection_context = selection_context,
        lower.tail        = TRUE
      )

    } else {

      # standard normal CDF
      cdf_vals <- .outcome_cdf.norm(
        yi         = yi_cdf,
        mu_samples = mu_samples_cdf,
        tau_within = tau_within_samples,
        sei        = sei,
        lower.tail = lower_tail
      )

    }

  } else if (outcome_type == "bin") {

    # binomial CDF using normal approximation
    cdf_vals <- .outcome_cdf.binom(
      yi         = yi,
      sei        = sei,
      mu_samples = mu_samples,
      tau_within = tau_within_samples
    )

  } else if (outcome_type == "pois") {

    # Poisson CDF using normal approximation
    cdf_vals <- .outcome_cdf.pois(
      yi         = yi,
      sei        = sei,
      mu_samples = mu_samples,
      tau_within = tau_within_samples
    )

  } else {

    stop("Unsupported outcome type for CDF computation.", call. = FALSE)

  }

  # add column names for observations
  colnames(cdf_vals) <- paste0("cdf[", seq_len(K), "]")

  return(cdf_vals)
}


# ---------------------------------------------------------------------------- #
# .cdf_lik_estimate.brma
# ---------------------------------------------------------------------------- #
#
# Compute CDF values for the estimate-unit LOO target.
#
# This mirrors `.log_lik_estimate.brma()`: fixed effects plus fitted cluster
# effects for multilevel models, marginal over estimate-level heterogeneity.
# It is used by LOO-PIT residuals so the PSIS weights and CDF target match.
#
# @param object brma object.
#
# @return S x K matrix of CDF values in [0, 1]
#
# ---------------------------------------------------------------------------- #
.cdf_lik_estimate.brma <- function(object, setup = NULL) {

  if (is.null(setup)) {
    setup <- .estimate_likelihood_setup.brma(object)
  }

  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  K                 <- setup[["K"]]
  mu_samples        <- setup[["mu"]]
  tau_within        <- setup[["tau_within"]]
  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  effect_direction  <- setup[["effect_direction"]]
  posterior_samples <- setup[["posterior_samples"]]

  if (.known_v_estimate_target_uses_backend(setup[["data"]])) {
    cdf_vals <- .cdf_known_v_estimate_target_from_setup(setup)
    colnames(cdf_vals) <- paste0("cdf_lik[", seq_len(K), "]")
    return(cdf_vals)
  }

  if (outcome_type == "norm") {

    if (effect_direction == "negative") {
      mu_samples_cdf <- -mu_samples
      yi_cdf         <- -yi
      lower_tail     <- FALSE
    } else {
      mu_samples_cdf <- mu_samples
      yi_cdf         <- yi
      lower_tail     <- TRUE
    }

    if (is_weightfunction) {
      selection_context <- .selection_context(
        object            = object,
        posterior_samples = posterior_samples
      )

      cdf_vals <- .outcome_cdf.selnorm(
        yi                = yi,
        mu_samples        = mu_samples,
        tau_within        = tau_within,
        sei               = sei,
        selection_context = selection_context,
        lower.tail        = TRUE
      )

    } else {

      cdf_vals <- .outcome_cdf.norm(
        yi         = yi_cdf,
        mu_samples = mu_samples_cdf,
        tau_within = tau_within,
        sei        = sei,
        lower.tail = lower_tail
      )
    }

  } else if (outcome_type == "bin") {

    cdf_vals <- .outcome_cdf.binom(
      yi         = yi,
      sei        = sei,
      mu_samples = mu_samples,
      tau_within = tau_within
    )

  } else if (outcome_type == "pois") {

    cdf_vals <- .outcome_cdf.pois(
      yi         = yi,
      sei        = sei,
      mu_samples = mu_samples,
      tau_within = tau_within
    )

  } else {

    stop("Unsupported outcome type for estimate-unit CDF computation.",
         call. = FALSE)
  }

  colnames(cdf_vals) <- paste0("cdf_lik[", seq_len(K), "]")

  return(cdf_vals)
}
