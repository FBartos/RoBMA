#' @title Predict From brma Object
#'
#' @description \code{predict.brma} predicts values
#'
#' @inheritParams summary.brma
#' @param newdata specification for prediction data. Defaults to \code{NULL} which
#' corresponds to prediction for the observed data. Alternatives are:
#' \itemize{
#'   \item{\code{TRUE} returns a single aggregated prediction (either the scalar
#'   parameter for models without moderators/scale, or the average across the
#'   model matrix for regression models). Not available for \code{type = "response"}
#'   since observation-level sampling variance is required.}
#'   \item{A data.frame (for meta-regression) or a named list with effect size
#'   measure and variability metrics (for meta-analysis) for new studies. The input
#'   must contain the variables used by the requested prediction type. Outcome
#'   columns are optional for \code{type = "terms"} unless PET/PEESE bias terms
#'   are included via \code{bias_adjusted = FALSE}.}
#' }
#' @param type type of prediction to be performed. Options are:
#' \itemize{
#'   \item{\code{"terms"} (alias: \code{"marginal"}): Fixed-effect parameters only (mu).
#'   Produces the mean effect size estimate at the given predictor levels,
#'   not accounting for random effects.}
#'   \item{\code{"cluster"}: Fixed effects plus cluster-level random effects (mu + gamma).
#'   Only available for multilevel (3-level) models.}
#'   \item{\code{"estimate"} (aliases: \code{"effect"}, \code{"blup"}): True
#'   effects (mu + gamma + theta). For existing normal data, returns posterior
#'   draws of conditional BLUP means, not simulated latent-effect draws. For
#'   new data, samples estimate-level random effects.}
#'   \item{\code{"response"} (alias: \code{"outcome"}): Predicted observed values (yi).
#'   Incorporates both heterogeneity and sampling variability.}
#'   \item{\code{"terms.scale"}: Scale parameter (tau), incorporating scale
#'   regression if present.}
#' }
#' @param as_measure logical; whether to convert GLMM response predictions from
#' simulated counts to continuity-corrected effect-size estimators (logOR for
#' binomial, logIRR for Poisson). Defaults to \code{TRUE}. Only relevant for
#' GLMM models with \code{type = "response"}. When \code{FALSE}, returns raw
#' frequency data (counts). Use \code{type = "estimate"} for latent logOR/logIRR
#' predictions.
#' @param output_measure effect-size measure for location/effect predictions.
#' Defaults to the fitted measure. Supported conversions are among \code{"SMD"},
#' \code{"COR"}, \code{"ZCOR"}, and \code{"OR"}; \code{"RR"}, \code{"HR"},
#' \code{"IRR"}, \code{"RD"}, and \code{"GEN"} can only be returned on their
#' fitted measure. Use \code{transform = "EXP"} for ratio-scale output from
#' log-scale measures.
#' @param transform optional display transformation. Currently \code{"EXP"}
#' exponentiates log-scale measures \code{"OR"}, \code{"RR"}, \code{"HR"},
#' and \code{"IRR"}.
#' @param bias_adjusted whether predictions should adjust for publication bias.
#' Defaults to \code{FALSE}. When \code{TRUE}:
#' \itemize{
#'   \item{PET/PEESE terms are NOT added to the mean parameter (mu), returning
#'   the bias-corrected effect estimate.}
#'   \item{For \code{type = "response"} with selection models, samples from
#'   the ordinary normal predictive distribution instead of the selected-normal
#'   distribution, simulating what would be observed without publication bias.}
#' }
#' When \code{FALSE}:
#' \itemize{
#'   \item{PET/PEESE terms ARE added to mu, returning predictions that include
#'   the expected bias (i.e., what we expect to observe given publication bias).}
#'   \item{For \code{type = "response"} with selection models, samples from
#'   the selected-normal distribution reflecting the selective publishing
#'   process.}
#' }
#' @param conditional whether to return conditional posterior predictions for
#' RoBMA product-space objects. For location predictions, samples are conditioned
#' on the effect component; for \code{type = "terms.scale"}, samples are
#' conditioned on the heterogeneity component.
#'
#' @details
#' \strong{Type hierarchy:}
#' \itemize{
#'   \item{\code{"terms"}: mu (fixed effects only)}
#'   \item{\code{"cluster"}: mu + gamma (adds cluster-level random effect)}
#'   \item{\code{"estimate"}: mu + gamma + theta (adds estimate-level random effect)}
#'   \item{\code{"response"}: mu + gamma + theta + epsilon (adds sampling error)}
#' }
#' For existing normal observations, \code{type = "estimate"} reports the
#' conditional BLUP mean \eqn{E(\theta_i \mid y_i, \mu_i, \tau_i)} for each
#' posterior draw. It is therefore an empirical-Bayes summary, not a draw from
#' the full latent-effect posterior \eqn{\theta_i \mid y_i}. For new data,
#' estimate-level random effects are sampled from their model distribution.
#'
#' Note that in contrast to \link[metafor]{predict}, the \code{type = "response"} produces
#' predictions for the new effect size estimates. To obtain results corresponding to
#' metafor's predict function, use \code{type = "terms"} for the mean effect size
#' and \code{type = "estimate"} for true effects (prediction interval).
#'
#' @examples \dontrun{
#' }
#'
#' @return A \code{brma_samples} object containing posterior samples. When printed,
#' displays a summary table via \code{BayesTools::ensemble_estimates_table}. The
#' underlying samples matrix can be accessed directly (the object inherits from matrix)
#' or via \code{summary()} to obtain the summary table. The samples can also be converted
#' to \pkg{posterior} draws formats using \code{as_draws()} and related functions.
#' @seealso [pooled_effect()], [pooled_heterogeneity()], [blup()]
#' @export
predict.brma <- function(object, newdata = NULL,
                         type = "terms",
                         as_measure = TRUE,
                         output_measure = NULL,
                         transform = NULL,
                         probs = c(.025, .975),
                         bias_adjusted = FALSE,
                         quiet = FALSE,
                         conditional = FALSE,
                         ...){

  # normalize type aliases
  type <- match.arg(type, c("terms", "marginal", "cluster", "estimate", "effect", "blup",
                            "response", "outcome", "terms.scale"))
  type <- switch(type,
    "marginal" = "terms",
    "effect"   = "estimate",
    "blup"     = "estimate",
    "outcome"  = "response",
    type  # default: keep as is
  )

  # input validation
  BayesTools::check_bool(as_measure, "as_measure")
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(quiet, "quiet")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' predictions are available only for RoBMA objects.",
         call. = FALSE)
  }

  # check newdata: NULL, TRUE, or data.frame/list
  if (!is.null(newdata) && !isTRUE(newdata)) {
    if (!is.data.frame(newdata) && !is.list(newdata)) {
      stop("'newdata' must be NULL, TRUE, a data.frame, or a named list.", call. = FALSE)
    }
    if (is.list(newdata) && !is.data.frame(newdata) && is.null(names(newdata))) {
      stop("'newdata' list must be named.", call. = FALSE)
    }
  }

  # check incompatible options: aggregate predictions not available for response
  if (isTRUE(newdata) && type == "response") {
    stop("Aggregated predictions (newdata = TRUE) are not available for type = 'response' ",
         "because observation-level sampling variance is required.", call. = FALSE)
  }

  ### types of predictions
  # terms:       fixed effects terms for the overall effect (mu) / incorporating mods if present
  # cluster:     terms + cluster-level random effects (mu + gamma) - multilevel only
  # terms.scale: fixed effects terms for the overall heterogeneity (tau) / incorporating scale if present
  # estimate:    incorporating heterogeneity into terms to obtain the true effects
  #              (via empirical Bayes for existing data, new random effect sampled for new data)
  # response:    incorporating heterogeneity and sampling variability

  ### dispatch between prediction on the current data vs. new data
  if (is.null(newdata)) {

    # existing data are used
    same_data <- TRUE
    aggregate <- FALSE
    new_data  <- object[["data"]]

  } else if (isTRUE(newdata)) {

    # aggregated prediction: use original data but marginalize random effects
    # and average across observations (for mods/scale models)
    same_data <- FALSE
    aggregate <- TRUE
    new_data  <- object[["data"]]

  } else {

    # prepare newdata using the same settings as the original fit
    same_data <- FALSE
    aggregate <- FALSE
    new_data  <- .prepare_newdata(
      object        = object,
      newdata       = newdata,
      type          = type,
      bias_adjusted = bias_adjusted
    )

  }

  ### extract priors and structural information about the model
  priors            <- object[["priors"]]
  is_mods           <- .is_mods(object)
  is_multilevel     <- .is_multilevel(object)
  is_scale          <- .is_scale(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)
  is_weights        <- .is_weights(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)
  posterior_samples <- .get_posterior_samples(object[["fit"]])
  effect_transform  <- .effect_output_setup(
    object         = object,
    output_measure = output_measure,
    transform      = transform
  )

  if (type == "terms.scale" && .effect_output_requested(effect_transform)) {
    stop(
      "'output_measure' and 'transform' are only available for effect-size ",
      "predictions, not for type = 'terms.scale'.",
      call. = FALSE
    )
  }
  if (type == "response" && !isTRUE(as_measure) &&
      .effect_output_requested(effect_transform)) {
    stop(
      "'output_measure' and 'transform' require type = 'response' predictions ",
      "to be returned as effect-size measures (as_measure = TRUE).",
      call. = FALSE
    )
  }

  # check: cluster type requires multilevel model
  if (type == "cluster" && !is_multilevel) {
    stop("type = 'cluster' is only available for multilevel (3-level) models. ",
         "Use type = 'terms' for non-multilevel models.", call. = FALSE)
  }

  ### extract outcome data and fit data for convenience
  outcome_data      <- new_data[["outcome"]]
  fit_data          <- .create_fit_data(data = new_data, priors = priors)

  # outcome dimensions
  K_original <- nrow(outcome_data)
  K          <- K_original

  ### obtain tau samples using helper function
  # returns list(tau_within, tau_between) - all S x K matrices
  # see .evaluate.brma.tau() in evaluate.R for details
  tau_result          <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = new_data[["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = new_data, "scale") else NULL,
    scale_priors      = priors[["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K_original,
    posterior_samples = posterior_samples
  )
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  ### aggregate tau samples if requested (for terms.scale only at this point)
  ### mu aggregation happens after mu computation
  if (aggregate && type == "terms.scale") {
    # average tau across observations (rows are samples, columns are observations)
    # for non-scale models, all columns are identical so this is a no-op
    # for scale models, this averages across the model matrix
    if (is_scale && K_original > 1 && !quiet) {
      message("Aggregated prediction averages tau across the scale model matrix (K = ", K_original, " observations).")
    }
    tau_within_samples  <- matrix(rowMeans(tau_within_samples), ncol = 1)
    tau_between_samples <- matrix(rowMeans(tau_between_samples), ncol = 1)
    K <- 1L
  }

  ### extract MCMC chain info for brma_samples construction
  n_chains <- length(object[["fit"]][["mcmc"]])
  n_iter   <- object[["fit"]][["sample"]]

  ### return only tau samples if type = "terms.scale" is selected
  if (type == "terms.scale") {
    # reconstruct total tau from components: tau = sqrt(tau_within^2 + tau_between^2)
    tau_samples <- sqrt(tau_within_samples^2 + tau_between_samples^2)

    # rename samples
    colnames(tau_samples) <- if (aggregate) "tau" else paste0("tau[", seq_len(K), "]")

    out <- .new_brma_samples(
      samples  = tau_samples,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = if (aggregate) "Aggregated Scale Term Posterior Prediction:" else "Scale Term Posterior Prediction:",
      probs    = probs,
      data     = if (aggregate) NULL else new_data
    )
    return(.condition_prediction_samples(
      object      = object,
      samples     = out,
      conditional = conditional,
      parameters  = .conditional_heterogeneity_parameters(object)
    ))
  }

  ### get the base mu samples using helper function
  # returns S x K matrix of location samples
  # see .evaluate.brma.mu() in evaluate.R for details
  mu_samples <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = outcome_data,
    mods_data         = new_data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = new_data, "mods") else NULL,
    mods_priors       = priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = bias_adjusted,
    K                 = K_original,
    posterior_samples = posterior_samples
  )

  ### aggregate mu and tau samples if requested (for terms/cluster/estimate types)
  if (aggregate && type %in% c("terms", "cluster", "estimate")) {
    # average mu across observations (rows are samples, columns are observations)
    # for non-mods models, all columns are identical so this is a no-op
    # for mods models, this averages across the model matrix
    if (is_mods && K_original > 1 && !quiet) {
      message("Aggregated prediction averages mu across the moderator model matrix (K = ", K_original, " observations).")
    }
    mu_samples          <- matrix(rowMeans(mu_samples), ncol = 1)
    tau_within_samples  <- matrix(rowMeans(tau_within_samples), ncol = 1)
    tau_between_samples <- matrix(rowMeans(tau_between_samples), ncol = 1)
    K <- 1L
  }

  ### return only mu samples if type = "terms" is selected
  # terms incorporate fixed effects only (i.e., random effects are not incorporated)
  if (type == "terms") {
    # rename samples
    colnames(mu_samples) <- if (aggregate) "mu" else paste0("mu[", seq_len(K), "]")

    out <- .new_effect_brma_samples(
      samples          = mu_samples,
      n_chains         = n_chains,
      n_iter           = n_iter,
      title            = if (aggregate) "Aggregated Location Term Posterior Prediction:" else "Location Term Posterior Prediction:",
      probs            = probs,
      data             = if (aggregate) NULL else new_data,
      effect_transform = effect_transform
    )
    return(.condition_prediction_samples(
      object      = object,
      samples     = out,
      conditional = conditional,
      parameters  = .conditional_effect_parameters(object)
    ))
  }

  ### include 3-level (cluster-level) random-effects using helper function
  # returns contribution matrix gamma[cluster] * tau_between for multilevel models
  # see .evaluate.brma.cluster_effects() in evaluate.R for details
  # for aggregated predictions (same_data = FALSE, K = 1), function samples new gamma ~ N(0,1)
  blup_vi <- NULL
  if (outcome_type == "norm") {
    blup_vi <- outcome_data[["sei"]]^2
    if (is_weights) {
      blup_vi <- blup_vi / outcome_data[["weights"]]
    }
  }

  blup_bias_offset <- NULL
  if (outcome_type == "norm" && same_data && bias_adjusted && (is_PET || is_PEESE)) {
    blup_bias_offset <- .evaluate.brma.bias_offset(
      fit               = object[["fit"]],
      outcome_data      = outcome_data,
      is_PET            = is_PET,
      is_PEESE          = is_PEESE,
      effect_direction  = effect_direction,
      K                 = K,
      posterior_samples = posterior_samples
    )
  }

  # Exact block BLUP applies when the Gaussian residual model matches the target.
  # For weightfunction models, gamma is sampled under the selected-data
  # likelihood; keep those draws and apply estimate-level shrinkage conditional
  # on them instead of replacing them by the ordinary Gaussian block shortcut.
  multilevel_blup     <- NULL
  use_multilevel_blup <- is_multilevel && outcome_type == "norm" && same_data &&
    !is_weightfunction

  if (use_multilevel_blup) {
    multilevel_blup <- .evaluate.brma.multilevel_blup.norm(
      mu_samples  = mu_samples,
      tau_within  = tau_within_samples,
      tau_between = tau_between_samples,
      yi          = outcome_data[["yi"]],
      vi          = blup_vi,
      cluster     = fit_data[["cluster"]],
      bias_offset = blup_bias_offset
    )
  }

  if (is_multilevel) {
    if (!is.null(multilevel_blup)) {
      cluster_contribution <- multilevel_blup[["cluster"]]
    } else {
      cluster <- if (aggregate) seq_len(K) else fit_data[["cluster"]]
      cluster_contribution <- .evaluate.brma.cluster_effects(
        fit               = object[["fit"]],
        tau_between       = tau_between_samples,
        cluster           = cluster,
        same_data         = same_data,
        effect_direction  = effect_direction,
        posterior_samples = posterior_samples
      )
    }
    mu_samples <- mu_samples + cluster_contribution
  }

  ### return cluster-level predictions if type = "cluster" is selected
  # cluster incorporates fixed effects + cluster-level random effects (mu + gamma)
  if (type == "cluster") {
    # rename samples
    colnames(mu_samples) <- if (aggregate) "mu_cluster" else paste0("mu_cluster[", seq_len(K), "]")

    out <- .new_effect_brma_samples(
      samples          = mu_samples,
      n_chains         = n_chains,
      n_iter           = n_iter,
      title            = if (aggregate) "Aggregated Cluster-Level Posterior Prediction:" else "Cluster-Level Posterior Prediction:",
      probs            = probs,
      data             = if (aggregate) NULL else new_data,
      effect_transform = effect_transform
    )
    return(.condition_prediction_samples(
      object      = object,
      samples     = out,
      conditional = conditional,
      parameters  = .conditional_effect_parameters(object)
    ))
  }

  ### create true-effect prediction
  # dispatches between normal and GLMM approaches
  # both helpers handle same_data dispatch internally:
  # - same_data = TRUE: use fitted values (BLUPs for normal, extracted theta for GLMM)
  # - same_data = FALSE: sample from marginal distribution N(mu, tau_within)
  if (type == "estimate") {

    if (!is.null(multilevel_blup)) {
      true_effects_samples <- mu_samples + multilevel_blup[["estimate"]]
    } else if (outcome_type == "norm") {
      true_effects_samples <- .evaluate.brma.true_effects.norm(
        mu_samples  = mu_samples,
        tau_within  = tau_within_samples,
        yi          = outcome_data[["yi"]],
        sei         = sqrt(blup_vi),
        same_data   = same_data,
        bias_offset = blup_bias_offset
      )
    } else {
      true_effects_samples <- .evaluate.brma.true_effects.glmm(
        fit               = object[["fit"]],
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        same_data         = same_data,
        K                 = K,
        posterior_samples = posterior_samples
      )
    }

    # rename samples
    colnames(true_effects_samples) <- if (aggregate) "theta" else paste0("theta[", seq_len(K), "]")

    out <- .new_effect_brma_samples(
      samples          = true_effects_samples,
      n_chains         = n_chains,
      n_iter           = n_iter,
      title            = if (same_data && outcome_type == "norm") {
        "True Effects (BLUP Means):"
      } else if (aggregate) {
        "Aggregated True Effect Posterior Prediction:"
      } else {
        "True Effect Posterior Prediction:"
      },
      probs            = probs,
      data             = if (aggregate) NULL else new_data,
      effect_transform = effect_transform
    )
    return(.condition_prediction_samples(
      object      = object,
      samples     = out,
      conditional = conditional,
      parameters  = .conditional_effect_parameters(object)
    ))

  }

  ### create observed effects prediction
  if (type == "response") {

    # different model types have different output structures
    if (is.element(outcome_type, c("bin", "pois"))) {

      if (as_measure) {
        .check_glmm_response_as_measure(outcome_type, outcome_data)
      }

      # the estimate-level effects are not marginalized for GLMMs
      # include the sampled random effects or marginalize over them for new data
      # see .evaluate.brma.theta.glmm() in evaluate.R for details
      theta_contribution <- .evaluate.brma.theta.glmm(
        fit               = object[["fit"]],
        tau_within        = tau_within_samples,
        same_data         = same_data,
        K                 = K,
        posterior_samples = posterior_samples
      )
      mu_samples <- mu_samples + theta_contribution

      if (outcome_type == "bin") {

        ### incorporate with base-rate using helper function
        # Existing data use fitted pi[k]. Newdata samples new pi from the
        # estimate-specific prior; reusing fitted pi[k] would index old studies.
        if (same_data) {
          if (K != nrow(object[["data"]][["outcome"]])) {
            stop("Existing-data GLMM response prediction has inconsistent K.",
                 call. = FALSE)
          }
          logit_baserate <- .evaluate.brma.baserate(
            fit               = object[["fit"]],
            K                 = K,
            posterior_samples = posterior_samples
          )
        } else {
          logit_baserate <- .evaluate.brma.baserate_newdata(
            prior_pi = priors[["outcome"]][["pi"]],
            S        = nrow(mu_samples),
            K        = K
          )
        }

        ### sample outcome using RNG helper
        outcome_samples <- .outcome_rng.binom(
          mu_samples     = mu_samples,
          logit_baserate = logit_baserate,
          n1i            = outcome_data[["n1i"]],
          n2i            = outcome_data[["n2i"]]
        )

      } else if (outcome_type == "pois") {

        ### incorporate with log-rate using helper function
        # Existing data use fitted phi[k]. Newdata samples new phi from the
        # estimate-specific prior; reusing fitted phi[k] would index old studies.
        if (same_data) {
          if (K != nrow(object[["data"]][["outcome"]])) {
            stop("Existing-data GLMM response prediction has inconsistent K.",
                 call. = FALSE)
          }
          log_phi <- .evaluate.brma.lograte(
            fit               = object[["fit"]],
            K                 = K,
            posterior_samples = posterior_samples
          )
        } else {
          log_phi <- .evaluate.brma.lograte_newdata(
            prior_phi = priors[["outcome"]][["phi"]],
            S         = nrow(mu_samples),
            K         = K
          )
        }

        ### sample outcome using RNG helper
        outcome_samples <- .outcome_rng.pois(
          mu_samples = mu_samples,
          log_phi    = log_phi,
          t1i        = outcome_data[["t1i"]],
          t2i        = outcome_data[["t2i"]]
        )

      }

      # convert to effect size measure if requested (default)
      if (as_measure) {
        outcome_samples <- .glmm_response_counts_to_measure(
          outcome_samples = outcome_samples,
          outcome_data    = outcome_data,
          outcome_type    = outcome_type
        )

      }


    } else if (outcome_type == "norm") {

      # normal outcome models dispatch between:
      # - bias_adjusted = TRUE: sample from unweighted normal (as if no bias)
      # - bias_adjusted = FALSE with weightfunction: sample from selected-normal
      # - bias_adjusted = FALSE without weightfunction: sample from unweighted normal

      if (bias_adjusted || !is_weightfunction) {

        # sample from unweighted normal distribution
        # y[s,k] ~ N(mu[s,k], sqrt(tau_within[s,k]^2 + sei[k]^2))
        outcome_samples <- .outcome_rng.norm(
          mu_samples = mu_samples,
          tau_within = tau_within_samples,
          sei        = outcome_data[["sei"]]
        )

      } else {

        selection_context <- .selection_context(
          object            = object,
          posterior_samples = posterior_samples,
          newdata           = new_data
        )

        outcome_samples <- .outcome_rng.selnorm(
          mu_samples        = mu_samples,
          tau_within        = tau_within_samples,
          sei               = outcome_data[["sei"]],
          selection_context = selection_context
        )

      }

      # rename samples
      colnames(outcome_samples) <- paste0("yi[", seq_len(K), "]")
    }

    out <- .new_effect_brma_samples(
      samples          = outcome_samples,
      n_chains         = n_chains,
      n_iter           = n_iter,
      title            = .response_prediction_title(outcome_type, as_measure),
      probs            = probs,
      data             = new_data,
      effect_transform = effect_transform
    )
    return(.condition_prediction_samples(
      object      = object,
      samples     = out,
      conditional = conditional,
      parameters  = .conditional_effect_parameters(object)
    ))
  }
}


# ---------------------------------------------------------------------------- #
# .check_glmm_response_as_measure
# ---------------------------------------------------------------------------- #
#
# Check whether simulated GLMM counts can be converted to effect-size estimator
# draws. Raw count predictions can allow zero exposure/person-time, but the
# derived logOR/logIRR estimator cannot.
#
# ---------------------------------------------------------------------------- #
.check_glmm_response_as_measure <- function(outcome_type, outcome_data) {

  if (outcome_type == "bin" &&
      (any(outcome_data[["n1i"]] <= 0, na.rm = TRUE) ||
       any(outcome_data[["n2i"]] <= 0, na.rm = TRUE))) {
    stop(
      "GLMM response predictions with as_measure = TRUE require positive ",
      "'n1i' and 'n2i'. Use as_measure = FALSE for raw count predictions ",
      "or type = 'estimate' for latent logOR predictions.",
      call. = FALSE
    )
  }

  if (outcome_type == "pois" &&
      (any(outcome_data[["t1i"]] <= 0, na.rm = TRUE) ||
       any(outcome_data[["t2i"]] <= 0, na.rm = TRUE))) {
    stop(
      "Poisson response predictions with as_measure = TRUE require positive ",
      "'t1i' and 't2i'. Use as_measure = FALSE for raw count predictions ",
      "or type = 'estimate' for latent logIRR predictions.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


# ---------------------------------------------------------------------------- #
# .glmm_response_counts_to_measure
# ---------------------------------------------------------------------------- #
#
# Convert simulated GLMM response counts to escalc-style effect-size estimator
# draws. This is not the latent effect posterior; it is the sampling distribution
# of the continuity-corrected estimator applied to posterior predictive counts.
#
# ---------------------------------------------------------------------------- #
.glmm_response_counts_to_measure <- function(outcome_samples, outcome_data,
                                             outcome_type) {

  K <- nrow(outcome_data)

  if (outcome_type == "bin") {
    ai <- outcome_samples[, seq(1, 2 * K, by = 2), drop = FALSE]
    ci <- outcome_samples[, seq(2, 2 * K, by = 2), drop = FALSE]

    n1i_mat <- matrix(outcome_data[["n1i"]], nrow = nrow(ai), ncol = K, byrow = TRUE)
    n2i_mat <- matrix(outcome_data[["n2i"]], nrow = nrow(ai), ncol = K, byrow = TRUE)
    bi      <- n1i_mat - ai
    di      <- n2i_mat - ci
    zero    <- ai == 0 | bi == 0 | ci == 0 | di == 0

    ai_adj <- ai + 0.5 * zero
    bi_adj <- bi + 0.5 * zero
    ci_adj <- ci + 0.5 * zero
    di_adj <- di + 0.5 * zero

    out <- log((ai_adj * di_adj) / (bi_adj * ci_adj))

  } else if (outcome_type == "pois") {
    x1i <- outcome_samples[, seq(1, 2 * K, by = 2), drop = FALSE]
    x2i <- outcome_samples[, seq(2, 2 * K, by = 2), drop = FALSE]

    t1i_mat <- matrix(outcome_data[["t1i"]], nrow = nrow(x1i), ncol = K, byrow = TRUE)
    t2i_mat <- matrix(outcome_data[["t2i"]], nrow = nrow(x1i), ncol = K, byrow = TRUE)
    zero    <- x1i == 0 | x2i == 0

    x1i_adj <- x1i + 0.5 * zero
    x2i_adj <- x2i + 0.5 * zero

    out <- log((x1i_adj / t1i_mat) / (x2i_adj / t2i_mat))

  } else {
    stop("Unsupported GLMM response outcome type.", call. = FALSE)
  }

  colnames(out) <- paste0("yi[", seq_len(K), "]")

  return(out)
}


# ---------------------------------------------------------------------------- #
# .response_prediction_title
# ---------------------------------------------------------------------------- #
#
# Human-readable prediction title matching the response estimand.
#
# ---------------------------------------------------------------------------- #
.response_prediction_title <- function(outcome_type, as_measure) {

  if (is.element(outcome_type, c("bin", "pois"))) {
    if (as_measure) {
      return("Observed Effect-Size Estimator Posterior Prediction:")
    }
    return("Observed Counts Posterior Prediction:")
  }

  return("Observations Posterior Prediction:")
}

.condition_prediction_samples <- function(object, samples, conditional,
                                          parameters) {

  if (!conditional) {
    return(samples)
  }

  keep <- .conditional_parameter_rows(
    object     = object,
    parameters = parameters,
    rule       = "OR"
  )
  sample_matrix <- as.matrix(samples)[keep, , drop = FALSE]

  return(.new_brma_samples(
    samples          = sample_matrix,
    n_chains         = attr(samples, "nchains"),
    n_iter           = attr(samples, "niter"),
    title            = paste("Conditional", attr(samples, "title")),
    probs            = attr(samples, "probs"),
    data             = attr(samples, "data"),
    effect_transform = attr(samples, "effect_transform")
  ))
}

.conditional_effect_parameters <- function(object) {

  prior_names <- names(attr(object[["fit"]], "prior_list"))

  if ("mu" %in% prior_names) {
    return("mu")
  }
  if ("mu_intercept" %in% prior_names) {
    return("mu_intercept")
  }

  parameters <- grep("^mu_", prior_names, value = TRUE)
  if (length(parameters) > 0L) {
    return(parameters)
  }

  stop("No location parameter is available for conditional prediction.",
       call. = FALSE)
}

.conditional_heterogeneity_parameters <- function(object) {

  prior_names <- names(attr(object[["fit"]], "prior_list"))

  if ("tau" %in% prior_names) {
    return("tau")
  }
  if ("log_tau_intercept" %in% prior_names) {
    return("log_tau_intercept")
  }

  parameters <- grep("^log_tau_", prior_names, value = TRUE)
  if (length(parameters) > 0L) {
    return(parameters)
  }

  stop("No heterogeneity parameter is available for conditional prediction.",
       call. = FALSE)
}

.conditional_parameter_rows <- function(object, parameters,
                                        posterior_samples = NULL,
                                        rule = "OR") {

  BayesTools::check_char(parameters, "parameters", check_length = 0)
  BayesTools::check_char(rule, "rule", allow_values = c("AND", "OR"))

  posterior_samples <- .get_posterior_samples(object[["fit"]], posterior_samples)
  conditions <- vapply(
    parameters,
    .conditional_parameter_rows_single,
    logical(nrow(posterior_samples)),
    object            = object,
    posterior_samples = posterior_samples
  )

  if (rule == "AND") {
    keep <- apply(conditions, 1, all)
  } else {
    keep <- apply(conditions, 1, any)
  }

  if (!any(keep)) {
    warning("No samples left after conditioning.", call. = FALSE,
            immediate. = TRUE)
  }

  return(keep)
}

.conditional_parameter_rows_single <- function(parameter, object,
                                               posterior_samples) {

  prior_list <- attr(object[["fit"]], "prior_list")
  if (!parameter %in% names(prior_list)) {
    stop("Missing prior for conditional parameter '", parameter, "'.",
         call. = FALSE)
  }

  prior <- prior_list[[parameter]]
  S     <- nrow(posterior_samples)

  if (BayesTools::is.prior.spike_and_slab(prior)) {
    indicator <- .conditional_indicator_column(posterior_samples, parameter)
    return(indicator == 1)
  }

  if (BayesTools::is.prior.mixture(prior)) {
    components <- attr(prior, "components")
    if (!all(components %in% c("null", "alternative"))) {
      stop(
        "Conditional mixture posterior distributions are available only ",
        "for 'null' and 'alternative' components.",
        call. = FALSE
      )
    }
    indicator <- .conditional_indicator_column(posterior_samples, parameter)
    return(indicator %in% which(components == "alternative"))
  }

  warning(
    "The parameter '", parameter, "' is not a conditional parameter. ",
    "All samples are assumed to come from the conditional posterior distribution.",
    call. = FALSE,
    immediate. = TRUE
  )
  return(rep(TRUE, S))
}

.conditional_indicator_column <- function(posterior_samples, parameter) {

  column <- paste0(parameter, "_indicator")
  if (!column %in% colnames(posterior_samples)) {
    stop("Missing posterior model indicator: '", column, "'.",
         call. = FALSE)
  }

  return(as.integer(posterior_samples[, column]))
}
