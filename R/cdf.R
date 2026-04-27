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
# @return S x K matrix of CDF values in (0, 1)
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
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

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
# .outcome_cdf.wnorm
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for weighted normal distribution (selection models).
#
# For selection models, the observed effect y_i follows a weighted normal:
#   y_i ~ f(y) * omega(y) / Z
# where f(y) is normal density and omega(y) is the selection weight function.
#
# This uses the .pwnorm_fast.ss function that handles the weighted CDF.
#
# Performance optimization: When use_normal is provided, samples with
# use_normal[i] = TRUE use the fast normal path (pnorm), skipping the
# expensive weighted CDF computation. This is correct because these samples
# have omega = all 1s (set by JAGS for non-weightfunction models).
#
# @param yi               numeric vector of length K; observed effect sizes
# @param mu_samples       S x K matrix of location samples
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
# @param sei              numeric vector of length K; standard errors
# @param omega            S x W matrix of omega (weight) samples
# @param crit_yi          W - 1 x K matrix of critical values for each observation
# @param use_normal       optional logical vector of length S; TRUE if sample should
#                         use fast normal path (for samples with omega = all 1s)
# @param bias_indicator   optional integer vector of length S; active bias branch
#                         for branch-specific cutpoint mapping
# @param crit_yi_mapping  optional cutpoint mapping matrix, cuts x branches
# @param crit_yi_mapping_max optional active cutoff count per branch
# @param lower.tail       logical; return P(Y <= yi) if TRUE, P(Y > yi) if FALSE
#
# @return S x K matrix of CDF values from weighted distribution
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.wnorm <- function(yi, mu_samples, tau_within, sei, omega, crit_yi,
                               use_normal = NULL, bias_indicator = NULL,
                               crit_yi_mapping = NULL,
                               crit_yi_mapping_max = NULL,
                               lower.tail = TRUE) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # pre-compute total SD for each observation
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # use the mapped dwnorm_mix backend when branch-specific cutpoint mapping is available
  has_mapping <- .selection_has_mapping(
    bias_indicator      = bias_indicator,
    crit_yi_mapping     = crit_yi_mapping,
    crit_yi_mapping_max = crit_yi_mapping_max
  )

  if (has_mapping) {
    return(.wnorm_mix_cdf_matrix(
      q                   = yi,
      mean                = mu_samples,
      sd                  = total_sd,
      omega               = omega,
      crit_yi             = crit_yi,
      bias_indicator      = bias_indicator,
      crit_yi_mapping     = crit_yi_mapping,
      crit_yi_mapping_max = crit_yi_mapping_max,
      lower.tail          = lower.tail,
      log.p               = FALSE
    ))
  }

  cdf_vals <- matrix(NA_real_, nrow = S, ncol = K)
  for (k in seq_len(K)) {
    cdf_vals[, k] <- .pwnorm_fast.ss(
      q          = rep(yi[k], S),
      mean       = mu_samples[, k],
      sd         = total_sd[, k],
      omega      = omega,
      crit_x     = crit_yi[, k],
      lower.tail = lower.tail,
      use_normal = use_normal  # <-- Pass through for internal subdispatch
    )
  }

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
# @return S x K matrix of CDF values in (0, 1)
#
# ---------------------------------------------------------------------------- #
.cdf.brma <- function(object, conditioning_depth = "marginal") {

  ### input validation
  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)

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
    posterior_samples = posterior_samples
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

    # flip for negative effect direction (applies to normal and weighted normal)
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

      # extract omega samples for weight function
      omega_samples <- .extract_omega_samples(posterior_samples)

      # get fit_data for crit_yi
      fit_data <- .create_fit_data(data = data, priors = priors)

      # compute use_normal indicator for performance optimization
      # this identifies which samples come from non-weightfunction bias models
      use_normal <- .extract_use_normal(object, posterior_samples = posterior_samples)

      cdf_vals <- .outcome_cdf.wnorm(
        yi                  = yi_cdf,
        mu_samples          = mu_samples_cdf,
        tau_within          = tau_within_samples,
        sei                 = sei,
        omega               = omega_samples,
        crit_yi             = fit_data$crit_yi,
        lower.tail          = lower_tail,
        use_normal          = use_normal,
        bias_indicator      = .extract_bias_indicator(object, posterior_samples = posterior_samples),
        crit_yi_mapping     = fit_data[["crit_yi_mapping"]],
        crit_yi_mapping_max = fit_data[["crit_yi_mapping_max"]]
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
# @return S x K matrix of CDF values in (0, 1)
#
# ---------------------------------------------------------------------------- #
.cdf_lik_estimate.brma <- function(object) {

  setup             <- .estimate_likelihood_setup.brma(object)
  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  K                 <- setup[["K"]]
  mu_samples        <- setup[["mu"]]
  tau_within        <- setup[["tau_within"]]
  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  effect_direction  <- setup[["effect_direction"]]
  posterior_samples <- setup[["posterior_samples"]]

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
      omega_samples <- .extract_omega_samples(posterior_samples)
      use_normal    <- .extract_use_normal(object, posterior_samples = posterior_samples)

      cdf_vals <- .outcome_cdf.wnorm(
        yi                  = yi_cdf,
        mu_samples          = mu_samples_cdf,
        tau_within          = tau_within,
        sei                 = sei,
        omega               = omega_samples,
        crit_yi             = setup[["fit_data"]][["crit_yi"]],
        lower.tail          = lower_tail,
        use_normal          = use_normal,
        bias_indicator      = .extract_bias_indicator(object, posterior_samples = posterior_samples),
        crit_yi_mapping     = setup[["fit_data"]][["crit_yi_mapping"]],
        crit_yi_mapping_max = setup[["fit_data"]][["crit_yi_mapping_max"]]
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
