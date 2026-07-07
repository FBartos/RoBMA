# ============================================================================ #
# evaluate.R
# ============================================================================ #
#
# This file contains modular helper functions for evaluating posterior samples
# from brma model fits. These functions extract and transform MCMC samples for:
# - heterogeneity (tau) parameters
# - location (mu) parameters
# - true effects (theta)
# - GLMM-specific parameters (baserate, lograte)
#
# The functions are designed for:
# - Prediction (predict.brma)
# - Simulation from the posterior predictive distribution
# - Likelihood evaluation
# - Other downstream tasks requiring posterior sample manipulation
#
# Design principles:
# - Decomposed parameters: minimize memory by passing only required components
# - Vectorized operations: use outer(), sweep(), matrix ops instead of loops
# - Consistent return structures: always return lists with predictable components
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# Posterior Extraction Helpers
# ---------------------------------------------------------------------------- #
#
# Keep posterior extraction centralized so helper functions can share an already
# materialized matrix in hot paths.
#
# ---------------------------------------------------------------------------- #
.get_posterior_samples <- function(fit, posterior_samples = NULL) {

  if (!is.null(posterior_samples)) {
    return(posterior_samples)
  }
  if (is.null(fit)) {
    stop(
      "Posterior samples are required for this operation; refit the model ",
      "or supply '.posterior_samples'.",
      call. = FALSE
    )
  }

  return(suppressWarnings(coda::as.mcmc(fit)))
}

.is_zero_point_prior <- function(prior) {

  value <- .point_prior_value(prior)

  return(!is.null(value) && isTRUE(all.equal(value, 0)))
}

.point_prior_value <- function(prior) {

  if (is.null(prior) || !BayesTools::is.prior.point(prior)) {
    return(NULL)
  }

  value <- prior[["parameters"]][["location"]]
  if (length(value) == 1L && is.finite(value)) {
    return(as.numeric(value))
  }

  return(NULL)
}

.fixed_tau_prior_value <- function(priors) {

  if (is.null(priors) || is.null(priors[["outcome"]])) {
    return(NULL)
  }

  return(.point_prior_value(priors[["outcome"]][["tau"]]))
}

.fixed_mu_prior_value <- function(priors) {

  if (is.null(priors) || is.null(priors[["outcome"]])) {
    return(NULL)
  }

  return(.point_prior_value(priors[["outcome"]][["mu"]]))
}

.fixed_bias_parameter_value <- function(priors, parameter) {

  if (is.null(priors) ||
      is.null(priors[["outcome"]]) ||
      is.null(priors[["outcome"]][["bias"]])) {
    return(NULL)
  }

  prior <- priors[["outcome"]][["bias"]]
  if (parameter == "PET" && !BayesTools::is.prior.PET(prior)) {
    return(NULL)
  }
  if (parameter == "PEESE" && !BayesTools::is.prior.PEESE(prior)) {
    return(NULL)
  }

  return(.point_prior_value(prior))
}

.fixed_prior_list_values <- function(prior_list) {

  if (is.null(prior_list) || length(prior_list) == 0L) {
    return(numeric())
  }

  values <- vapply(prior_list, function(prior) {
    value <- .point_prior_value(prior)
    if (is.null(value)) NA_real_ else value
  }, numeric(1))
  values <- values[is.finite(values)]

  return(values)
}

.resolve_fixed_prior_row <- function(row, prior_list) {

  fixed_values <- .fixed_prior_list_values(prior_list)
  if (length(fixed_values) == 0L) {
    return(row)
  }

  for (parameter in names(fixed_values)) {
    row[[parameter]] <- fixed_values[[parameter]]
  }

  return(row)
}

.resolve_fixed_prior_sample_columns <- function(posterior_samples, prior_list) {

  fixed_values <- .fixed_prior_list_values(prior_list)
  if (length(fixed_values) == 0L) {
    return(posterior_samples)
  }

  columns <- colnames(posterior_samples)
  if (is.null(columns)) {
    return(posterior_samples)
  }

  fixed_columns <- intersect(names(fixed_values), columns)
  for (column in fixed_columns) {
    posterior_samples[, column] <- fixed_values[[column]]
  }

  missing_values <- fixed_values[!names(fixed_values) %in% columns]
  if (length(missing_values) == 0L) {
    return(posterior_samples)
  }

  fixed_samples <- matrix(
    rep(missing_values, each = nrow(posterior_samples)),
    nrow     = nrow(posterior_samples),
    ncol     = length(missing_values),
    dimnames = list(NULL, names(missing_values))
  )

  return(cbind(posterior_samples, fixed_samples))
}

.has_fixed_zero_tau_prior <- function(priors) {

  if (is.null(priors) || is.null(priors[["outcome"]])) {
    return(FALSE)
  }

  return(.is_zero_point_prior(priors[["outcome"]][["tau"]]))
}

.resolve_missing_tau_value <- function(allow_missing_tau) {

  if (isTRUE(allow_missing_tau)) {
    return(0)
  }
  if (is.numeric(allow_missing_tau) &&
      length(allow_missing_tau) == 1L &&
      is.finite(allow_missing_tau)) {
    return(as.numeric(allow_missing_tau))
  }

  return(NULL)
}

.posterior_or_fixed_scalar <- function(posterior_samples, parameter,
                                       fixed_value = NULL,
                                       required = TRUE) {

  if (!is.null(fixed_value)) {
    return(rep(fixed_value, nrow(posterior_samples)))
  }
  if (parameter %in% colnames(posterior_samples)) {
    return(as.numeric(posterior_samples[, parameter]))
  }
  if (isTRUE(required)) {
    stop("Missing posterior ", parameter, " columns.", call. = FALSE)
  }

  return(NULL)
}


# ---------------------------------------------------------------------------- #
# .extract_posterior_matrix
# ---------------------------------------------------------------------------- #
#
# Extract vectorized JAGS parameters such as theta[1], ..., theta[K] as an
# S x K matrix, preserving the requested column order.
#
# ---------------------------------------------------------------------------- #
.extract_posterior_matrix <- function(posterior_samples, parameter, K) {

  column_names <- paste0(parameter, "[", seq_len(K), "]")
  missing_cols <- setdiff(column_names, colnames(posterior_samples))

  if (length(missing_cols) > 0) {
    stop("Missing posterior column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  return(as.matrix(posterior_samples[, column_names, drop = FALSE]))
}


# ---------------------------------------------------------------------------- #
# .extract_indexed_parameter_samples
# ---------------------------------------------------------------------------- #
#
# Extract scalar/vectorized posterior parameter columns, sorting fully indexed
# names such as `omega[10]` by numeric index instead of column order.
#
# ---------------------------------------------------------------------------- #
.extract_indexed_parameter_samples <- function(posterior_samples, parameter,
                                               n_expected = NULL,
                                               required = TRUE) {

  parameter_matrix <- BayesTools::JAGS_indexed_parameter_matrix(
    samples   = posterior_samples,
    parameter = parameter
  )
  if (is.null(parameter_matrix) && parameter %in% colnames(posterior_samples)) {
    parameter_matrix <- as.matrix(posterior_samples[, parameter, drop = FALSE])
  }

  if (is.null(parameter_matrix)) {
    if (required) {
      stop("Missing posterior ", parameter, " columns.", call. = FALSE)
    }
    return(NULL)
  }

  if (!is.null(n_expected) && ncol(parameter_matrix) != n_expected) {
    stop(
      "Expected ", n_expected, " posterior ", parameter, " column(s), found ",
      ncol(parameter_matrix), ".",
      call. = FALSE
    )
  }

  return(parameter_matrix)
}


# ---------------------------------------------------------------------------- #
# .extract_omega_samples
# ---------------------------------------------------------------------------- #
#
# Extract selection-model omega columns from posterior samples in numeric order.
#
# @param posterior_samples posterior sample matrix.
#
# @return S x W matrix of omega samples.
#
# ---------------------------------------------------------------------------- #
.extract_omega_samples <- function(posterior_samples) {

  return(.extract_indexed_parameter_samples(posterior_samples, "omega"))
}


# ---------------------------------------------------------------------------- #
# .posterior_formula_fit
# ---------------------------------------------------------------------------- #
#
# Build a BayesTools::JAGS_evaluate_formula() input from already selected
# posterior rows while preserving formula scaling metadata from the JAGS fit.
#
# ---------------------------------------------------------------------------- #
.posterior_formula_fit <- function(fit, posterior_samples,
                                   formula_design = TRUE) {

  formula_fit <- if (inherits(posterior_samples, "mcmc")) {
    posterior_samples
  } else {
    coda::as.mcmc(as.matrix(posterior_samples))
  }

  attr(formula_fit, "formula_scale") <- attr(fit, "formula_scale")
  if (isTRUE(formula_design)) {
    attr(formula_fit, "formula_design") <- attr(fit, "formula_design")
  } else {
    attr(formula_fit, "formula_design") <- NULL
  }

  return(formula_fit)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.tau
# ---------------------------------------------------------------------------- #
#
# Extract and compute heterogeneity (tau) posterior samples from a brma fit.
#
# This function handles:
# - Scale regression models: evaluates log_tau formula, then exponentiates
# - Simple models: extracts tau column and replicates to K columns
# - Multilevel models: splits tau into within/between components via rho
#
# @param fit         runjags fit object containing posterior samples
# @param scale_data  data.frame of scale predictors for formula evaluation
#                    (NULL if not a scale regression model)
# @param scale_formula formula object for scale regression (from attr(data$scale, "formula"))
# @param scale_priors list of priors for scale parameters
# @param is_scale    logical; whether model uses scale regression
# @param is_multilevel logical; whether model is 3-level (clustered)
# @param K           integer; number of observations (determines output columns)
#
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return A list with three components (all S x K matrices):
#   - tau_total: total heterogeneity
#   - tau_within: estimate-level heterogeneity component
#   - tau_between: cluster-level heterogeneity component
#   For non-multilevel models: tau_within = tau (total), tau_between = 0 matrix
#   Total tau can be reconstructed as: sqrt(tau_within^2 + tau_between^2)
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.tau <- function(fit, scale_data, scale_formula, scale_priors,
                               is_scale, is_multilevel, K,
                               posterior_samples = NULL,
                               allow_missing_tau = FALSE) {

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  S <- nrow(posterior_samples)  # number of posterior samples

  ### compute tau samples based on model type
  if (is_scale) {

    scale_priors <- .repair_formula_prior_list(
      prior_list = scale_priors,
      parameter  = "log_tau"
    )

    # scale regression: evaluate log_tau formula then exponentiate
    # BayesTools::JAGS_evaluate_formula returns K x S matrix, we need S x K
    log_tau_samples <- t(BayesTools::JAGS_evaluate_formula(
      fit            = .posterior_formula_fit(fit, posterior_samples),
      formula        = scale_formula,
      parameter      = "log_tau",
      data           = scale_data,
      prior_list     = scale_priors,
      formula_target = "fixed"
    ))
    tau_samples <- exp(log_tau_samples)

  } else {

    fixed_tau <- if (is.numeric(allow_missing_tau)) {
      .resolve_missing_tau_value(allow_missing_tau)
    } else {
      NULL
    }

    if (!is.null(fixed_tau)) {
      tau_samples <- matrix(fixed_tau, nrow = S, ncol = K)
    } else {
      tau_parameter <- .extract_indexed_parameter_samples(
        posterior_samples = posterior_samples,
        parameter         = "tau",
        required          = FALSE
      )
      if (is.null(tau_parameter)) {
        missing_tau <- .resolve_missing_tau_value(allow_missing_tau)
        if (is.null(missing_tau) &&
            !any(grepl("__xRE", colnames(posterior_samples), fixed = TRUE))) {
          stop("Missing posterior tau columns.", call. = FALSE)
        }
        if (is.null(missing_tau)) {
          missing_tau <- 0
        }
        tau_samples <- matrix(missing_tau, nrow = S, ncol = K)
      } else if (ncol(tau_parameter) == 1L) {
        # simple model: extract tau column and replicate to K columns
        # matrix(vec, nrow = S, ncol = K) replicates vec across columns
        # equivalent to: for (k in 1:K) result[, k] <- tau_parameter[, 1]
        tau_samples <- matrix(tau_parameter[, 1L], nrow = S, ncol = K)
      } else if (ncol(tau_parameter) == K) {
        tau_samples <- tau_parameter
      } else {
        stop(
          "Posterior tau samples must be scalar or row-wise.",
          call. = FALSE
        )
      }
    }

  }

  ### split tau into within/between components for multilevel models
  if (is_multilevel) {

    # extract rho (proportion of variance at cluster-level)
    rho_samples <- posterior_samples[, "rho"]

    # clamp rho to [0, 1] to handle JAGS numerical precision issues
    # pmin/pmax are vectorized min/max: faster than rho[rho > 1] <- 1
    rho_samples <- pmin(pmax(rho_samples, 0), 1)

    # tau_within = tau * sqrt(1 - rho)  (estimate-level heterogeneity)
    # tau_between = tau * sqrt(rho)     (cluster-level heterogeneity)
    # multiplication by vector rho_samples broadcasts across columns
    tau_within_samples  <- tau_samples * sqrt(1 - rho_samples)
    tau_between_samples <- tau_samples * sqrt(rho_samples)

  } else {

    # non-multilevel: all heterogeneity is at estimate-level
    # tau_between is zero matrix for consistent interface
    tau_within_samples  <- tau_samples
    tau_between_samples <- matrix(0, nrow = S, ncol = K)

  }

  return(list(
    tau_total   = tau_samples,
    tau_within  = tau_within_samples,
    tau_between = tau_between_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.mu
# ---------------------------------------------------------------------------- #
#
# Extract and compute location (mu) posterior samples from a brma fit.
#
# This function handles:
# - Meta-regression models: evaluates mu formula with moderators
# - Simple models: extracts mu column and replicates to K columns
# - PET/PEESE bias shifts (added when bias_adjusted = FALSE)
#
# The returned mu is on the original effect-size scale. For
# effect_direction = "negative", the JAGS likelihood already uses -mu against
# internally flipped outcomes, so this helper must not flip mu itself.
#
# @param fit              runjags fit object containing posterior samples
# @param outcome_data     data.frame with outcome info (must contain 'sei')
# @param mods_data        data.frame of moderators for formula evaluation
#                         (NULL if not a meta-regression model)
# @param mods_formula     formula object for meta-regression
# @param mods_priors      list of priors for moderator parameters
# @param is_mods          logical; whether model is meta-regression
# @param is_PET           logical; whether model includes PET adjustment
# @param is_PEESE         logical; whether model includes PEESE adjustment
# @param effect_direction character; "positive" or "negative" - direction of true effect
# @param bias_adjusted    logical; if TRUE, PET/PEESE adjustments are skipped
#                         (returns bias-adjusted mu); if FALSE, adds bias terms
# @param K                integer; number of observations
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return S x K matrix of mu (location) posterior samples
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.mu <- function(fit, outcome_data, mods_data, mods_formula, mods_priors,
                              is_mods, is_PET, is_PEESE, effect_direction,
                              bias_adjusted, K, posterior_samples = NULL,
                              priors = NULL) {

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  S <- nrow(posterior_samples)

  ### compute base mu samples
  if (is_mods) {

    mods_priors <- .repair_formula_prior_list(
      prior_list = mods_priors,
      parameter  = "mu"
    )

    # meta-regression: evaluate mu formula with moderators
    # returns K x S, transpose to S x K
    mu_samples <- t(BayesTools::JAGS_evaluate_formula(
      fit            = .posterior_formula_fit(fit, posterior_samples),
      formula        = mods_formula,
      parameter      = "mu",
      data           = mods_data,
      prior_list     = mods_priors,
      formula_target = "fixed"
    ))

  } else {

    fixed_mu <- .fixed_mu_prior_value(priors)
    if (!is.null(fixed_mu)) {

      # point prior: mu is fixed and therefore absent from posterior samples
      mu_samples <- matrix(fixed_mu, nrow = S, ncol = K)

    } else if ("mu" %in% colnames(posterior_samples)) {

      # simple model: replicate mu column to K columns
      mu_samples <- matrix(posterior_samples[, "mu"], nrow = S, ncol = K)

    } else {

      location_data <- data.frame(row.names = seq_len(K))
      mu_samples <- t(BayesTools::JAGS_evaluate_formula(
        fit            = .posterior_formula_fit(
          fit               = fit,
          posterior_samples = posterior_samples,
          formula_design    = FALSE
        ),
        formula        = stats::as.formula("~ 1"),
        parameter      = "mu",
        data           = location_data,
        prior_list     = .repair_formula_prior_list(
          prior_list = mods_priors,
          parameter  = "mu"
        ),
        formula_target = "fixed"
      ))
    }

  }

  # NOTE: No effect direction flipping needed here.
  # The JAGS model uses `-mu` in the likelihood for negative effects, so `mu`
  # remains on the original effect-size scale. The data flip in
  # .create_fit_data() is matched by the likelihood sign.

  ### add PET adjustment when NOT incorporating publication bias adjustment
  # (i.e., when we want to show the biased predictions)
  # PET model in JAGS: yi_flipped ~ N(-mu + PET * sei, tau) for negative effects
  #                    yi ~ N(mu + PET * sei, tau) for positive effects
  # To get biased effect in original scale:
  #   - positive: E[yi] = mu + PET * sei
  #   - negative: E[yi_original] = -E[yi_flipped] = -(-mu + PET * sei) = mu - PET * sei
  # So the sign of PET adjustment depends on effect_direction
  if (is_PET && !bias_adjusted) {

    PET_samples <- .posterior_or_fixed_scalar(
      posterior_samples = posterior_samples,
      parameter         = "PET",
      fixed_value       = .fixed_bias_parameter_value(priors, "PET")
    )
    sei_vec     <- outcome_data[["sei"]]

    # vectorized: outer(PET_samples, sei_vec) creates S x K matrix
    # outer(a, b) computes a[i] * b[j] for all i,j pairs
    # direction multiplier: +1 for positive, -1 for negative
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * outer(PET_samples, sei_vec)

  }

  ### add PEESE adjustment when NOT incorporating publication bias adjustment
  # PEESE model: Same logic as PET but with sei^2
  if (is_PEESE && !bias_adjusted) {

    PEESE_samples <- .posterior_or_fixed_scalar(
      posterior_samples = posterior_samples,
      parameter         = "PEESE",
      fixed_value       = .fixed_bias_parameter_value(priors, "PEESE")
    )
    sei_sq_vec    <- outcome_data[["sei"]]^2

    # direction multiplier: +1 for positive, -1 for negative
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * outer(PEESE_samples, sei_sq_vec)

  }

  return(mu_samples)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.random_effects
# ---------------------------------------------------------------------------- #
#
# Evaluate the conditional random-effect contribution for same-data
# random-formula models. This intentionally returns only the group-level
# contribution; callers add it to `.evaluate.brma.mu()` explicitly.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.random_effects <- function(fit, data, priors,
                                          posterior_samples = NULL,
                                          same_data = TRUE,
                                          required = FALSE,
                                          formula_target = "conditional",
                                          blocks = NULL,
                                          object = NULL) {

  formula_target <- match.arg(formula_target, c("conditional", "marginal"))

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  S                 <- nrow(posterior_samples)
  K                 <- nrow(data[["outcome"]])

  if (is.null(object)) {
    object <- list(
      fit    = fit,
      data   = data,
      priors = priors
    )
  }

  if (!.is_data_random(data)) {
    if (required) {
      stop("Random-effect predictions require a random-formula model.",
           call. = FALSE)
    }
    return(matrix(0, nrow = S, ncol = K))
  }

  location_priors <- attr(fit, "prior_list")
  if (is.null(location_priors)) {
    location_priors <- priors[["location"]]
  }
  if (is.null(location_priors)) {
    stop("Random-effect prediction requires location prior metadata.",
         call. = FALSE)
  }

  formula_design <- if (.is_scale(object)) {
    .predict_known_v_formula_design_with_row_source_values(
      object = object,
      data   = data
    )
  } else {
    .fitted_formula_design(object, "mu", required = TRUE)
  }
  if (is.null(formula_design)) {
    stop(
      "Conditional random-effect evaluation requires fitted formula-design metadata.",
      call. = FALSE
    )
  }

  if (identical(formula_target, "conditional")) {
    sampled_blocks <- .formula_design_sampled_random_effect_blocks(formula_design)
    if (is.null(blocks)) {
      blocks <- sampled_blocks
    }
    if (length(blocks) == 0L) {
      return(matrix(0, nrow = S, ncol = K))
    }
  }
  has_known_group_covariance <- .formula_design_blocks_have_known_group_covariance(
    formula_design = formula_design,
    blocks         = blocks
  )

  formula_fit <- .posterior_formula_fit(
    fit               = fit,
    posterior_samples = posterior_samples,
    formula_design    = TRUE
  )
  attr(formula_fit, "formula_design") <- list(mu = formula_design)

  call_args <- list(
    fit            = formula_fit,
    formula        = .create_fit_formula_list(data = data, parameter = "location"),
    parameter      = "mu",
    data           = data[["location"]],
    prior_list     = location_priors,
    formula_target = formula_target
  )

  if (identical(formula_target, "conditional")) {
    call_args[["blocks"]]     <- blocks
    call_args[["new_levels"]] <- if (isTRUE(same_data) || has_known_group_covariance) {
      "error"
    } else {
      "sample"
    }
  } else {
    if (!is.null(blocks)) {
      call_args[["blocks"]] <- blocks
    }
    call_args[["marginal_method"]] <- if (has_known_group_covariance) {
      "covariance"
    } else {
      "sample"
    }
    call_args[["new_levels"]]      <- if (has_known_group_covariance) "error" else "sample"
  }

  prediction     <- do.call(BayesTools::JAGS_predict_formula, call_args)
  random_samples <- prediction[["random"]]
  if (is.null(random_samples)) {
    random_samples <- .random_effects_from_marginal_vcov(
      prediction = prediction,
      S          = S,
      K          = K
    )
  }
  if (is.null(random_samples)) {
    return(matrix(0, nrow = S, ncol = K))
  }

  return(t(random_samples))
}


.random_effects_from_marginal_vcov <- function(prediction, S, K) {

  vcov <- prediction[["vcov"]]
  if (is.null(vcov) || is.null(vcov[["samples"]])) {
    return(NULL)
  }

  covariance_samples <- vcov[["samples"]]
  if (length(dim(covariance_samples)) != 3L ||
      dim(covariance_samples)[1L] != S ||
      dim(covariance_samples)[2L] != K ||
      dim(covariance_samples)[3L] != K) {
    stop("Random-effect covariance samples have inconsistent dimensions.",
         call. = FALSE)
  }

  zero_mu <- matrix(0, nrow = S, ncol = K)
  return(t(.outcome_rng.norm_known_v_covariance(
    mu_samples         = zero_mu,
    covariance_samples = covariance_samples
  )))
}


.evaluate.brma.random_effects_components <- function(object, data = object[["data"]],
                                                     posterior_samples = NULL,
                                                     same_data = TRUE) {

  posterior_samples <- .get_posterior_samples(object[["fit"]], posterior_samples)
  formula_design    <- if (.is_scale(object)) {
    .predict_known_v_formula_design_with_row_source_values(
      object = object,
      data   = data
    )
  } else {
    .fitted_formula_design(object, "mu", required = TRUE)
  }
  sampled_blocks <- .formula_design_sampled_random_effect_blocks(formula_design)
  if (length(sampled_blocks) == 0L) {
    return(list())
  }

  out <- lapply(sampled_blocks, function(block) {
    .evaluate.brma.random_effects(
      fit               = object[["fit"]],
      data              = data,
      priors            = object[["priors"]],
      posterior_samples = posterior_samples,
      same_data         = same_data,
      required          = TRUE,
      formula_target    = "conditional",
      blocks            = block,
      object            = object
    )
  })
  names(out) <- sampled_blocks

  return(out)
}


.evaluate.brma.mv_marginalized_random_blup.norm <- function(object, mu_samples,
                                                            sampled_random_samples = NULL,
                                                            posterior_samples = NULL,
                                                            bias_offset = NULL) {

  data  <- object[["data"]]
  terms <- .data_marginalized_random_effects(data)
  if (length(terms) == 0L) {
    return(list())
  }
  if (!.is_data_known_v(data)) {
    stop("Marginalized brma.mv random-effect BLUPs require known-V metadata.",
         call. = FALSE)
  }

  posterior_samples <- .get_posterior_samples(object[["fit"]], posterior_samples)
  S                 <- nrow(posterior_samples)
  K                 <- nrow(data[["outcome"]])

  if (!identical(dim(mu_samples), c(S, K))) {
    stop("'mu_samples' must have dimensions posterior draw x observation.",
         call. = FALSE)
  }
  if (is.null(sampled_random_samples)) {
    sampled_random_samples <- matrix(0, nrow = S, ncol = K)
  } else if (!identical(dim(sampled_random_samples), c(S, K))) {
    stop("'sampled_random_samples' must have dimensions posterior draw x observation.",
         call. = FALSE)
  }
  if (is.null(bias_offset)) {
    bias_offset <- matrix(0, nrow = S, ncol = K)
  } else if (!identical(dim(bias_offset), c(S, K))) {
    stop("'bias_offset' must have dimensions posterior draw x observation.",
         call. = FALSE)
  }

  known_V <- .data_known_v_data(data)
  V       <- known_V[["V"]]
  if (!is.matrix(V) || !identical(dim(V), c(K, K))) {
    stop("Known-V covariance dimensions do not match fitted rows.",
         call. = FALSE)
  }

  block_indices <- known_V[["block_indices"]]
  if (is.null(block_indices)) {
    block_indices <- list(seq_len(K))
  }

  source_samples <- .predict_known_v_newdata_marginalized_source_samples(
    object            = object,
    data              = data,
    posterior_samples = posterior_samples
  )
  total_variance <- .evaluate_marginalized_random_variance(
    data              = data,
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = source_samples
  )
  residual <- matrix(data[["outcome"]][["yi"]], nrow = S, ncol = K, byrow = TRUE) -
    bias_offset - mu_samples - sampled_random_samples

  out <- lapply(terms, function(term) {
    term_sd <- .marginalized_random_effect_sd_samples(
      term              = term,
      posterior_samples = posterior_samples,
      K                 = K,
      source_samples    = source_samples
    )
    term_variance <- .expand_random_effect_variance(term_sd^2, S = S, K = K)
    contribution  <- matrix(0, nrow = S, ncol = K)

    for (idx in block_indices) {
      V_block <- V[idx, idx, drop = FALSE]
      for (s in seq_len(S)) {
        covariance <- V_block
        diag(covariance) <- diag(covariance) + total_variance[s, idx]
        chol_covariance <- try(chol(covariance), silent = TRUE)
        if (inherits(chol_covariance, "try-error")) {
          stop(
            "Cannot solve marginalized random-effect BLUP covariance block; ",
            "covariance is not positive definite.",
            call. = FALSE
          )
        }
        weights <- backsolve(
          chol_covariance,
          forwardsolve(t(chol_covariance), residual[s, idx])
        )
        contribution[s, idx] <- term_variance[s, idx] * weights
      }
    }

    contribution
  })
  names(out) <- vapply(terms, `[[`, character(1), "block_name")

  return(out)
}


.expand_random_effect_variance <- function(variance, S, K) {

  if (nrow(variance) != S) {
    stop("Random-effect variance draw count does not match posterior samples.",
         call. = FALSE)
  }
  if (ncol(variance) == 1L) {
    return(matrix(variance[, 1L], nrow = S, ncol = K))
  }
  if (ncol(variance) == K) {
    return(variance)
  }

  stop("Random-effect variance must be scalar or row-wise.",
       call. = FALSE)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.scale_terms
# ---------------------------------------------------------------------------- #
#
# Evaluate row-wise scale formulas. For component-specific random-formula scale
# models this returns one block of columns per random component source.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.scale_terms <- function(fit, data, priors,
                                       posterior_samples = NULL,
                                       as_list = FALSE) {

  if (!.is_data_scale(data)) {
    stop("Scale predictions are not available because the model has no scale formula.",
         call. = FALSE)
  }

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  scale_specs       <- .data_scale_component_specs(data)
  scale_samples     <- lapply(scale_specs, function(scale_spec) {

    scale_priors <- .repair_formula_prior_list(
      prior_list = .create_fit_scale_formula_prior_list(
        priors    = priors,
        parameter = scale_spec[["parameter"]]
      ),
      parameter  = scale_spec[["parameter"]]
    )
    log_scale_samples <- t(BayesTools::JAGS_evaluate_formula(
      fit        = .posterior_formula_fit(
        fit               = fit,
        posterior_samples = posterior_samples,
        formula_design    = FALSE
      ),
      formula    = .create_fit_scale_formula(scale_spec[["formula"]]),
      parameter  = scale_spec[["parameter"]],
      data       = scale_spec[["data"]],
      prior_list = scale_priors
    ))
    scale_samples <- exp(log_scale_samples)
    if (as_list) {
      colnames(scale_samples) <- paste0("tau[", seq_len(ncol(scale_samples)), "]")
    } else {
      colnames(scale_samples) <- paste0(
        scale_spec[["source"]], "[",
        seq_len(ncol(scale_samples)), "]"
      )
    }
    scale_samples
  })

  if (as_list) {
    names(scale_samples) <- .evaluate.brma.scale_component_names(scale_specs)
    return(scale_samples)
  }

  return(do.call(cbind, scale_samples))
}


.evaluate.brma.aggregate_scale_terms_list <- function(scale_samples) {

  lapply(scale_samples, function(component_samples) {
    out <- matrix(rowMeans(component_samples), ncol = 1L)
    colnames(out) <- "tau"
    out
  })
}


.evaluate.brma.scale_component_names <- function(scale_specs) {

  vapply(scale_specs, `[[`, character(1), "display_name")
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.known_v_blup.norm
# ---------------------------------------------------------------------------- #
#
# Exact same-data BLUP means for independent latent effects observed through a
# known covariance matrix V: theta | y, mu, tau.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.known_v_blup.norm <- function(mu_samples, tau_within, yi,
                                             known_V, bias_offset = NULL) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  if (is.null(bias_offset)) {
    bias_offset <- matrix(0, nrow = S, ncol = K)
  } else if (!identical(dim(bias_offset), c(S, K))) {
    stop("'bias_offset' must have the same dimensions as 'mu_samples'.",
         call. = FALSE)
  }
  if (is.null(known_V) || is.null(known_V[["V"]])) {
    stop("Known-V BLUP requires known-V metadata.", call. = FALSE)
  }

  V <- known_V[["V"]]
  if (!is.matrix(V) || !identical(dim(V), c(K, K))) {
    stop("Known-V covariance dimensions do not match prediction rows.",
         call. = FALSE)
  }

  block_indices <- known_V[["block_indices"]]
  if (is.null(block_indices)) {
    block_indices <- list(seq_len(K))
  }

  true_effects_samples <- mu_samples
  for (idx in block_indices) {
    n_block      <- length(idx)
    V_block      <- V[idx, idx, drop = FALSE]
    tau_block    <- tau_within[, idx, drop = FALSE]
    residual     <- matrix(yi[idx], nrow = S, ncol = n_block, byrow = TRUE) -
      bias_offset[, idx, drop = FALSE] - mu_samples[, idx, drop = FALSE]
    constant_tau <- n_block == 1L ||
      isTRUE(all(tau_block == tau_block[, 1L]))

    if (n_block == 1L) {
      tau2        <- tau_block[, 1L]^2
      denominator <- tau2 + V_block[1L, 1L]
      if (any(!is.finite(denominator) | denominator <= 0)) {
        stop("Cannot solve known-V BLUP covariance block; covariance is not positive definite.",
             call. = FALSE)
      }
      true_effects_samples[, idx] <- mu_samples[, idx] +
        residual[, 1L] * tau2 / denominator
      next
    }

    if (constant_tau) {
      eigen_v     <- eigen(V_block, symmetric = TRUE)
      tau2        <- tau_block[, 1L]^2
      denominator <- outer(tau2, eigen_v[["values"]], "+")
      if (any(!is.finite(denominator) | denominator <= 0)) {
        stop("Cannot solve known-V BLUP covariance block; covariance is not positive definite.",
             call. = FALSE)
      }
      solved <- (residual %*% eigen_v[["vectors"]] / denominator) %*%
        t(eigen_v[["vectors"]])
      true_effects_samples[, idx] <- mu_samples[, idx, drop = FALSE] +
        solved * tau2
      next
    }

    for (s in seq_len(S)) {
      tau2    <- tau_block[s, ]^2
      M_block <- V_block
      diag(M_block) <- diag(M_block) + tau2
      chol_m  <- try(chol(M_block), silent = TRUE)
      if (inherits(chol_m, "try-error")) {
        stop("Cannot solve known-V BLUP covariance block; covariance is not positive definite.",
             call. = FALSE)
      }
      solved <- backsolve(chol_m, forwardsolve(t(chol_m), residual[s, ]))
      true_effects_samples[s, idx] <- mu_samples[s, idx] + tau2 * solved
    }
  }

  return(true_effects_samples)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.bias_offset
# ---------------------------------------------------------------------------- #
#
# Compute the PET/PEESE bias offset on the original effect-size scale.
#
# For observed-data BLUPs with bias_adjusted = TRUE, the residual must be formed
# from the bias-corrected estimate:
#   yi - bias_offset - mu
# rather than:
#   yi - mu
#
# @param fit              runjags fit object containing posterior samples
# @param outcome_data     data.frame with outcome info (must contain 'sei')
# @param is_PET           logical; whether model includes PET adjustment
# @param is_PEESE         logical; whether model includes PEESE adjustment
# @param effect_direction character; "positive" or "negative" - direction of true effect
# @param K                integer; number of observations
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return S x K matrix of additive bias offsets.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.bias_offset <- function(fit, outcome_data, is_PET, is_PEESE,
                                       effect_direction, K,
                                       posterior_samples = NULL,
                                       priors = NULL) {

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  S                 <- nrow(posterior_samples)
  bias_offset       <- matrix(0, nrow = S, ncol = K)
  direction         <- if (effect_direction == "negative") -1 else 1

  if (is_PET) {
    PET_samples <- .posterior_or_fixed_scalar(
      posterior_samples = posterior_samples,
      parameter         = "PET",
      fixed_value       = .fixed_bias_parameter_value(priors, "PET")
    )
    bias_offset <- bias_offset +
      direction * outer(PET_samples, outcome_data[["sei"]])
  }

  if (is_PEESE) {
    PEESE_samples <- .posterior_or_fixed_scalar(
      posterior_samples = posterior_samples,
      parameter         = "PEESE",
      fixed_value       = .fixed_bias_parameter_value(priors, "PEESE")
    )
    bias_offset <- bias_offset +
      direction * outer(PEESE_samples, outcome_data[["sei"]]^2)
  }

  return(bias_offset)
}


# ---------------------------------------------------------------------------- #
# .get_multilevel_block_indices
# ---------------------------------------------------------------------------- #
#
# Create a reusable block index list for multilevel covariance calculations.
#
# @param cluster integer vector of cluster identifiers
#
# @return A named list of integer index vectors, one per cluster.
#
# ---------------------------------------------------------------------------- #
.get_multilevel_block_indices <- function(cluster) {

  split(seq_along(cluster), cluster)
}


# ---------------------------------------------------------------------------- #
# .build_multilevel_marginal_covariance
# ---------------------------------------------------------------------------- #
#
# Construct the marginal covariance matrix for a 3-level normal model.
#
# The covariance decomposes into:
# - sampling variance: diag(vi)
# - estimate-level heterogeneity: diag(tau_within^2)
# - cluster-level heterogeneity: block-wise tcrossprod(tau_between)
#
# @param tau_within    numeric vector of length K with estimate-level SDs
# @param tau_between   numeric vector of length K with cluster-level SDs
# @param vi            numeric vector of length K with sampling variances
# @param block_indices list of observation indices for each cluster
#
# @return A K x K marginal covariance matrix.
#
# ---------------------------------------------------------------------------- #
.build_multilevel_marginal_covariance <- function(tau_within, tau_between, vi,
                                                  block_indices) {

  K <- length(vi)

  marginal_covariance <- diag(vi + tau_within^2, nrow = K, ncol = K)

  for (idx in block_indices) {
    marginal_covariance[idx, idx] <- marginal_covariance[idx, idx] +
      tcrossprod(tau_between[idx])
  }

  return(marginal_covariance)
}


# ---------------------------------------------------------------------------- #
# .solve_diagonal_rank_one_block
# ---------------------------------------------------------------------------- #
#
# Solve (diag(diagonal) + rank_one %*% t(rank_one))^-1 residual. The analytic
# Sherman-Morrison path is O(K) per cluster and falls back to Cholesky when a
# diagonal element is non-positive.
#
# ---------------------------------------------------------------------------- #
.solve_diagonal_rank_one_block <- function(diagonal, rank_one, residual) {

  if (!all(is.finite(diagonal)) ||
      !all(is.finite(rank_one)) ||
      !all(is.finite(residual))) {
    stop("Cannot solve known-V BLUP covariance block with non-finite inputs.",
         call. = FALSE)
  }

  if (all(is.finite(diagonal)) && all(diagonal > 0)) {
    inv_diag_residual <- residual / diagonal
    inv_diag_rank_one <- rank_one / diagonal
    denom             <- 1 + sum(rank_one * inv_diag_rank_one)

    if (is.finite(denom) && denom > .Machine$double.eps) {
      correction <- inv_diag_rank_one * sum(rank_one * inv_diag_residual) / denom
      return(inv_diag_residual - correction)
    }
  }

  covariance <- diag(diagonal, nrow = length(diagonal), ncol = length(diagonal)) +
    tcrossprod(rank_one)
  chol_m <- try(chol(covariance), silent = TRUE)

  if (inherits(chol_m, "try-error")) {
    stop("Cannot solve known-V BLUP covariance block; covariance is not positive definite.",
         call. = FALSE)
  }
  weights <- chol2inv(chol_m)

  return(as.vector(weights %*% residual))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.multilevel_blup.norm
# ---------------------------------------------------------------------------- #
#
# Compute exact same-data BLUP contributions for 3-level normal models.
#
# The conditional predictor is obtained from the marginal mixed-model identity:
#   b_hat = G M^{-1} (y - X beta)
# where G is decomposed into cluster-level and estimate-level components.
#
# @param mu_samples   S x K matrix of fixed-effect predictions
# @param tau_within   S x K matrix of estimate-level SDs
# @param tau_between  S x K matrix of cluster-level SDs
# @param yi           numeric vector of observed outcomes
# @param vi           numeric vector of sampling variances
# @param cluster      integer vector of cluster identifiers
#
# @return A list with S x K matrices `cluster` and `estimate`.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.multilevel_blup.norm <- function(mu_samples, tau_within, tau_between,
                                                yi, vi, cluster,
                                                bias_offset = NULL) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  if (is.null(bias_offset)) {
    bias_offset <- matrix(0, nrow = S, ncol = K)
  } else if (!identical(dim(bias_offset), c(S, K))) {
    stop("'bias_offset' must have the same dimensions as 'mu_samples'.",
         call. = FALSE)
  }

  block_indices <- .get_multilevel_block_indices(cluster)

  cluster_contribution  <- matrix(0, nrow = S, ncol = K)
  estimate_contribution <- matrix(0, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    for (idx in block_indices) {
      weighted_residual <- .solve_diagonal_rank_one_block(
        diagonal = vi[idx] + tau_within[s, idx]^2,
        rank_one = tau_between[s, idx],
        residual = yi[idx] - bias_offset[s, idx] - mu_samples[s, idx]
      )

      estimate_contribution[s, idx] <- tau_within[s, idx]^2 * weighted_residual

      cluster_scale <- sum(tau_between[s, idx] * weighted_residual)
      cluster_contribution[s, idx] <- tau_between[s, idx] * cluster_scale
    }
  }

  return(list(
    cluster  = cluster_contribution,
    estimate = estimate_contribution
  ))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.cluster_effects
# ---------------------------------------------------------------------------- #
#
# Extract or sample cluster-level (gamma) random effects for multilevel models.
#
# For multilevel models, gamma[cluster] represents the standardized cluster-level
# random effect (i.e., gamma ~ N(0, 1)). The actual contribution to mu is
# gamma * tau_between.
#
# This function handles:
# - Same data: extracts fitted gamma samples from posterior
# - New data: samples new gamma from N(0, 1) (marginalizes over cluster effects)
#
# @param fit              runjags fit object (needed to extract gamma if same_data)
# @param tau_between      S x K matrix of cluster-level heterogeneity samples
# @param cluster          integer vector of length K; cluster ID for each observation
# @param same_data        logical; TRUE if predicting on original data (use fitted gamma)
# @param effect_direction character; "positive" or "negative" (currently unused but kept for interface)
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return S x K matrix: contribution from cluster-level random effects
#         (gamma[cluster[k]] * tau_between[,k])
#         Can be added directly to mu samples.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.cluster_effects <- function(fit, tau_between, cluster,
                                           same_data, effect_direction,
                                           posterior_samples = NULL) {

  S <- nrow(tau_between)
  K <- ncol(tau_between)

  # NOTE: No direction flipping needed for cluster effects!
  # The JAGS model uses: "-gamma*tau_between" for negative effects, "+gamma*tau_between" for positive
  # But when converting to original scale:
  # E[yi_original] = -E[yi_flipped] = -(-mu - gamma*tau_between) = mu + gamma*tau_between
  # So the contribution to the original-scale effect is always +gamma*tau_between

  if (same_data) {

    # extract fitted gamma samples from posterior
    posterior_samples <- .get_posterior_samples(fit, posterior_samples)
    n_clusters        <- max(cluster)

    # extract all gamma columns at once: S x n_clusters matrix
    gamma_samples <- .extract_posterior_matrix(
      posterior_samples = posterior_samples,
      parameter         = "gamma",
      K                 = n_clusters
    )

    # gamma_samples[, cluster] reorders columns to match observations
    # this is S x K after reordering, element-wise multiply with tau_between
    cluster_contribution <- gamma_samples[, cluster, drop = FALSE] * tau_between

  } else {

    # new data: sample fresh gamma ~ N(0, 1) for each supplied cluster
    cluster_levels <- sort(unique(cluster))
    n_clusters     <- length(cluster_levels)
    gamma_new      <- matrix(stats::rnorm(S * n_clusters), nrow = S, ncol = n_clusters)
    cluster_map    <- match(cluster, cluster_levels)

    cluster_contribution <- gamma_new[, cluster_map, drop = FALSE] * tau_between

  }

  return(cluster_contribution)
}


.evaluate.brma.sampling_dependency <- function(fit, data, posterior_samples = NULL) {

  known_V <- .data_known_v_data(data)
  if (is.null(known_V) || !.is_data_known_v_parameterization(data, "latent") ||
      known_V[["rank"]] == 0L) {
    K <- nrow(data[["outcome"]])
    posterior_samples <- .get_posterior_samples(fit, posterior_samples)
    return(matrix(0, nrow = nrow(posterior_samples), ncol = K))
  }

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  z_samples <- .extract_posterior_matrix(
    posterior_samples = posterior_samples,
    parameter         = "sampling_z",
    K                 = known_V[["rank"]]
  )

  sampling_dependency <- z_samples %*% t(known_V[["B"]])
  if (.data_effect_direction(data) == "negative") {
    sampling_dependency <- -sampling_dependency
  }

  return(sampling_dependency)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.true_effects.norm
# ---------------------------------------------------------------------------- #
#
# Compute posterior samples of true-effect summaries for normal models.
# For same-data selection models, the selection weight is constant after
# conditioning on yi, so the estimate-level shrinkage remains Gaussian
# conditional on the fitted location and heterogeneity draw.
#
# For same_data = TRUE: Uses empirical Bayes shrinkage to return the
# conditional BLUP mean, E(theta_i | y_i, mu_i, tau_i), for each posterior row:
#   theta_i = lambda_i * y_i + (1 - lambda_i) * mu_i
# where:
#   lambda_i = tau_within^2 / (tau_within^2 + se_i^2)
#
# If bias_offset is supplied, the observed estimate is first corrected by the
# posterior-row PET/PEESE offset:
#   theta_i = mu_i + lambda_i * (y_i - bias_offset_i - mu_i)
#
# For same_data = FALSE: Samples from the marginal distribution of true effects:
#   theta_i = mu_i + epsilon_i * tau_within_i, where epsilon_i ~ N(0, 1)
#
# IMPORTANT: For multilevel models, mu_samples must already include the
# cluster-level contribution (gamma * tau_between) before calling this function.
#
# @param mu_samples       S x K matrix of location samples (must include
#                         gamma * tau_between contribution for multilevel models)
# @param tau_within       S x K matrix of estimate-level heterogeneity
# @param yi               numeric vector of length K; observed effect sizes (used only if same_data = TRUE)
# @param sei              numeric vector of length K; standard errors (used only if same_data = TRUE)
# @param same_data        logical; TRUE for BLUP estimates, FALSE for marginal sampling
# @param bias_offset      optional S x K matrix of PET/PEESE offsets to subtract
#                         from yi before BLUP shrinkage.
#
# @return S x K matrix of true-effect BLUP means or new-effect posterior samples
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.true_effects.norm <- function(mu_samples, tau_within, yi, sei,
                                             same_data, bias_offset = NULL) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  if (same_data) {

    # BLUP: empirical Bayes conditional means for existing observations
    # lambda = tau^2 / (tau^2 + se^2) ranges from 0 (strong shrinkage)
    # to 1 (weak shrinkage).
    lambda <- tau_within^2
    lambda <- lambda / sweep(lambda, 2, sei^2, "+")

    if (is.null(bias_offset)) {
      bias_offset <- matrix(0, nrow = S, ncol = K)
    } else if (!identical(dim(bias_offset), c(S, K))) {
      stop("'bias_offset' must have the same dimensions as 'mu_samples'.",
           call. = FALSE)
    }

    # BLUP: weighted average of observed effect and model prediction
    # high tau -> lambda -> 1 -> trust data more
    # low tau -> lambda -> 0 -> trust model more
    true_effects_samples <- mu_samples +
      lambda * (sweep(mu_samples, 2, yi, function(mu, yi_i) yi_i - mu) - bias_offset)

  } else {

    # Marginal sampling: sample new theta from N(mu, tau_within)
    # epsilon ~ N(0, 1), then theta = mu + epsilon * tau_within
    epsilon <- matrix(stats::rnorm(S * K), nrow = S, ncol = K)
    true_effects_samples <- mu_samples + epsilon * tau_within

  }

  return(true_effects_samples)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.true_effects.glmm
# ---------------------------------------------------------------------------- #
#
# Compute posterior samples of true effects (theta) for GLMM models.
#
# For GLMM models (binomial or Poisson), the estimate-level random effects
# (theta) are directly sampled in JAGS (not marginalized as in normal models).
# The true effect is:
#   true_effect_i = mu_i + theta_i * tau_within_i
#
# For same_data: extracts theta[k] from posterior samples
# For new_data: samples new theta ~ N(0, 1) to marginalize over random effects
#
# IMPORTANT: For multilevel models, mu_samples must already include the
# cluster-level contribution (gamma * tau_between) before calling this function.
# This is handled by .evaluate.brma.cluster_effects() in predict.R.
#
# @param fit              runjags fit object (needed to extract theta if same_data = TRUE)
# @param mu_samples       S x K matrix of location samples (must include
#                         gamma * tau_between contribution for multilevel models)
# @param tau_within       S x K matrix of estimate-level heterogeneity
# @param same_data        logical; TRUE if predicting on original data
# @param K                integer; number of observations
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return S x K matrix of true effect (theta) posterior samples
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.true_effects.glmm <- function(fit, mu_samples, tau_within, same_data, K,
                                             posterior_samples = NULL) {

  # add the estimate-level random effects (theta * tau_within) to mu
  theta_contribution <- .evaluate.brma.theta.glmm(
    fit               = fit,
    tau_within        = tau_within,
    same_data         = same_data,
    K                 = K,
    posterior_samples = posterior_samples
  )
  true_effects_samples <- mu_samples + theta_contribution

  return(true_effects_samples)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.baserate
# ---------------------------------------------------------------------------- #
#
# Extract base-rate (pi) samples for binomial GLMM models.
#
# For binomial outcomes, logit(pi[i]) is the midpoint of the two arm logits.
# The log-odds ratio effect size is applied symmetrically around this midpoint.
#
# @param fit              runjags fit object containing posterior samples
# @param K                integer; number of observations
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return A matrix of logit(pi) samples for each observation
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.baserate <- function(fit, K, posterior_samples = NULL) {

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  pi_samples        <- .extract_posterior_matrix(
    posterior_samples = posterior_samples,
    parameter         = "pi",
    K                 = K
  )

  return(.logit(pi_samples))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.baserate_newdata
# ---------------------------------------------------------------------------- #
#
# Sample new-study base-rate values for binomial GLMM response prediction.
#
# @param prior_pi BayesTools prior object on pi.
# @param S        integer; number of posterior rows.
# @param K        integer; number of new observations.
#
# @return S x K matrix of logit(pi) samples.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.baserate_newdata <- function(prior_pi, S, K) {

  pi_samples <- .draw_prior_samples_matrix(prior = prior_pi, S = S, K = K)
  pi_samples <- pmin(pmax(pi_samples, .Machine$double.eps), 1 - .Machine$double.eps)

  return(.logit(pi_samples))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.lograte
# ---------------------------------------------------------------------------- #
#
# Extract log-rate (phi) samples for Poisson GLMM models.
#
# For Poisson outcomes, phi[i] is the midpoint of the two arm log incidence
# rates. The incidence rate ratio effect size is applied symmetrically around
# this midpoint.
#
# @param fit              runjags fit object containing posterior samples
# @param K                integer; number of observations
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return A matrix of log-rate (phi) samples for each observation
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.lograte <- function(fit, K, posterior_samples = NULL) {

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)

  return(.extract_posterior_matrix(
    posterior_samples = posterior_samples,
    parameter         = "phi",
    K                 = K
  ))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.lograte_newdata
# ---------------------------------------------------------------------------- #
#
# Sample new-study log-rate values for Poisson GLMM response prediction.
#
# @param prior_phi BayesTools prior object on log(phi).
# @param S         integer; number of posterior rows.
# @param K         integer; number of new observations.
#
# @return S x K matrix of log-rate samples.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.lograte_newdata <- function(prior_phi, S, K) {

  return(.draw_prior_samples_matrix(prior = prior_phi, S = S, K = K))
}


# ---------------------------------------------------------------------------- #
# .draw_prior_samples_matrix
# ---------------------------------------------------------------------------- #
#
# Draw independent prior samples for GLMM newdata nuisance parameters.
#
# @param prior BayesTools prior object.
# @param S     integer; number of posterior rows.
# @param K     integer; number of new observations.
#
# @return S x K matrix of samples.
#
# ---------------------------------------------------------------------------- #
.draw_prior_samples_matrix <- function(prior, S, K) {

  values <- .draw_prior_samples(prior = prior, n = S * K)

  return(matrix(values, nrow = S, ncol = K))
}


# ---------------------------------------------------------------------------- #
# .draw_prior_samples
# ---------------------------------------------------------------------------- #
#
# Thin wrapper around BayesTools prior RNG.
#
# @param prior BayesTools prior object.
# @param n     integer; number of draws.
#
# @return numeric vector.
#
# ---------------------------------------------------------------------------- #
.draw_prior_samples <- function(prior, n) {

  return(as.numeric(BayesTools::rng(prior, n = n)))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.theta.glmm
# ---------------------------------------------------------------------------- #
#
# Extract or sample estimate-level random effects (theta) for GLMM models.
#
# For GLMMs, theta[i] represents the standardized estimate-level random effect
# (i.e., theta ~ N(0, 1)). The actual random effect is theta * tau_within.
#
# @param fit              runjags fit object (needed to extract theta if same_data)
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
# @param same_data        logical; TRUE if predicting on original data
# @param K                integer; number of observations
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return S x K matrix: mu contribution from estimate-level random effects
#         (theta[k] * tau_within[,k])
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.theta.glmm <- function(fit, tau_within, same_data, K,
                                      posterior_samples = NULL) {

  S <- nrow(tau_within)

  if (same_data) {

    posterior_samples  <- .get_posterior_samples(fit, posterior_samples)
    theta_contribution <- .extract_posterior_matrix(
      posterior_samples = posterior_samples,
      parameter         = "theta",
      K                 = K
    ) * tau_within

  } else {

    # new data: sample fresh theta ~ N(0, 1) for each observation
    theta_new          <- matrix(stats::rnorm(S * K), nrow = S, ncol = K)
    theta_contribution <- theta_new * tau_within

  }

  return(theta_contribution)
}


.logit <- function(p) {
  return(stats::qlogis(p))
}
.inv_logit <- function(x) {
  return(stats::plogis(x))
}


# ---------------------------------------------------------------------------- #
# .extract_use_normal
# ---------------------------------------------------------------------------- #
#
# Extract bias indicator and compute use_normal logical vector.
#
# For RoBMA objects with mixture bias priors, `bias_indicator` tracks which
# bias model generated each posterior sample. For brma objects with a single
# bias prior, no indicator column exists (all samples use the same model).
#
# This function is used for performance optimization in PDF/CDF/RNG functions.
# When `use_normal[i] = TRUE`, the sample uses a non-weightfunction bias model
# (PET, PEESE, or no bias), so the fast normal path can be used instead of
# the selected-normal computation.
#
# @param object brma or RoBMA object.
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return A logical vector of length S (number of posterior samples):
#   - TRUE if sample uses non-weightfunction bias model (PET, PEESE, or no bias)
#   - FALSE if sample uses weightfunction bias model
#
# ---------------------------------------------------------------------------- #
.extract_use_normal <- function(object, posterior_samples = NULL) {

  posterior_samples <- .get_posterior_samples(object[["fit"]], posterior_samples)
  routing           <- .selection_row_routing(
    priors               = object[["priors"]],
    posterior_samples    = posterior_samples
  )

  return(routing[["use_normal"]])
}


# ---------------------------------------------------------------------------- #
# .extract_bias_indicator
# ---------------------------------------------------------------------------- #
#
# Extract posterior bias branch indicators. Single-prior models do not monitor
# `bias_indicator`, but their branch is always 1.
#
# @param object brma object.
# @param posterior_samples optional posterior sample matrix; avoids repeated
#                          coda extraction when supplied by callers.
#
# @return integer vector of length S.
#
# ---------------------------------------------------------------------------- #
.extract_bias_indicator <- function(object, posterior_samples = NULL) {

  posterior_samples <- .get_posterior_samples(object[["fit"]], posterior_samples)
  routing           <- .selection_row_routing(
    priors               = object[["priors"]],
    posterior_samples    = posterior_samples
  )

  return(routing[["bias_indicator"]])
}
