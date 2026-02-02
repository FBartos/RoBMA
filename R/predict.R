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
#'   must correspond to the format and naming used in the original fit.}
#' }
#' @param type type of prediction to be performed. Options are:
#' \itemize{
#'   \item{\code{"terms"} (alias: \code{"marginal"}): Fixed-effect parameters only (mu).
#'   Produces the mean effect size estimate at the given predictor levels,
#'   not accounting for random effects.}
#'   \item{\code{"study"}: Fixed effects plus study-level random effects (mu + gamma).
#'   Only available for multilevel (3-level) models.}
#'   \item{\code{"estimate"} (aliases: \code{"effect"}, \code{"blup"}): True study
#'   effects (mu + gamma + theta). Incorporates all random effects including estimate-level
#'   heterogeneity. For existing data, returns BLUPs (empirical Bayes estimates).}
#'   \item{\code{"response"} (alias: \code{"outcome"}): Predicted observed values (yi).
#'   Incorporates both heterogeneity and sampling variability.}
#'   \item{\code{"terms.scale"}: Scale parameter (tau), incorporating scale
#'   regression if present.}
#' }
#' @param as_measure logical; whether to return GLMM response predictions as
#' effect size measures (logOR for binomial, logIRR for Poisson). Defaults to
#' \code{TRUE}. Only relevant for GLMM models with \code{type = "response"}.
#' When \code{FALSE}, returns raw frequency data (counts).
#' @param bias_adjusted whether predictions should adjust for publication bias.
#' Defaults to \code{FALSE}. When \code{TRUE}:
#' \itemize{
#'   \item{PET/PEESE terms are NOT added to the mean parameter (mu), returning
#'   the bias-corrected effect estimate.}
#'   \item{For \code{type = "response"} with selection models, samples from
#'   unweighted normal distribution instead of weighted distribution, simulating
#'   what would be observed without publication bias.}
#' }
#' When \code{FALSE}:
#' \itemize{
#'   \item{PET/PEESE terms ARE added to mu, returning predictions that include
#'   the expected bias (i.e., what we expect to observe given publication bias).}
#'   \item{For \code{type = "response"} with selection models, samples from
#'   weighted distribution reflecting the selective publishing process.}
#' }
#'
#' @details
#' \strong{Type hierarchy:}
#' \itemize{
#'   \item{\code{"terms"}: mu (fixed effects only)}
#'   \item{\code{"study"}: mu + gamma (adds study-level random effect)}
#'   \item{\code{"estimate"}: mu + gamma + theta (adds estimate-level random effect)}
#'   \item{\code{"response"}: mu + gamma + theta + epsilon (adds sampling error)}
#' }
#'
#' Note that in contrast to \link[metafor]{predict}, the \code{type = "response"} produces
#' predictions for the new effect size estimates. To obtain results corresponding to
#' metafor's predict function, use \code{type = "terms"} for the mean effect size
#' and \code{type = "estimate"} for true study effects (prediction interval).
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
                         probs = c(.025, .975),
                         bias_adjusted = FALSE,
                         quiet = FALSE,
                         ...){

  # normalize type aliases
  type <- match.arg(type, c("terms", "marginal", "study", "estimate", "effect", "blup",
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
  BayesTools::check_bool(quiet, "quiet")

  # check newdata: NULL, TRUE, or data.frame/list
  if (!is.null(newdata) && !isTRUE(newdata)) {
    if (!is.data.frame(newdata) && !is.list(newdata)) {
      stop("'newdata' must be NULL, TRUE, a data.frame, or a named list.", call. = FALSE)
    }
  }

  # check incompatible options: aggregate predictions not available for response
  if (isTRUE(newdata) && type == "response") {
    stop("Aggregated predictions (newdata = TRUE) are not available for type = 'response' ",
         "because observation-level sampling variance is required.", call. = FALSE)
  }

  ### types of predictions
  # terms:       fixed effects terms for the overall effect (mu) / incorporating mods if present
  # study:       terms + study-level random effects (mu + gamma) - multilevel only
  # terms.scale: fixed effects terms for the overall heterogeneity (tau) / incorporating scale if present
  # estimate:    incorporating between-study heterogeneity into terms to obtain the true study effects
  #              (via empirical Bayes for existing data, new random effect sampled for new data)
  # response:    incorporating between-study heterogeneity and sampling variability

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
      type          = type
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
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)

  # check: study type requires multilevel model
  if (type == "study" && !is_multilevel) {
    stop("type = 'study' is only available for multilevel (3-level) models. ",
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
  # see .evaluate.brma.tau() in brma.evaluate.R for details
  tau_result          <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = new_data[["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = new_data, "scale") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K_original
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

    return(.new_brma_samples(
      samples  = tau_samples,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = if (aggregate) "Aggregated Scale Term Posterior Prediction:" else "Scale Term Posterior Prediction:",
      probs    = probs,
      data     = if (aggregate) NULL else new_data
    ))
  }

  ### get the base mu samples using helper function
  # returns S x K matrix of location samples
  # see .evaluate.brma.mu() in brma.evaluate.R for details
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
    K                 = K_original
  )

  ### aggregate mu and tau samples if requested (for terms/study/estimate types)
  if (aggregate && type %in% c("terms", "study", "estimate")) {
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

    return(.new_brma_samples(
      samples  = mu_samples,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = if (aggregate) "Aggregated Location Term Posterior Prediction:" else "Location Term Posterior Prediction:",
      probs    = probs,
      data     = if (aggregate) NULL else new_data
    ))
  }

  ### include 3-level (study-level) random-effects using helper function
  # returns contribution matrix gamma[study_id] * tau_between for multilevel models
  # see .evaluate.brma.study_effects() in brma.evaluate.R for details
  # for aggregated predictions (same_data = FALSE, K = 1), function samples new gamma ~ N(0,1)
  if (is_multilevel) {
    study_contribution <- .evaluate.brma.study_effects(
      fit              = object[["fit"]],
      tau_between      = tau_between_samples,
      study_ids        = fit_data[["study_ids"]],
      same_data        = same_data,
      effect_direction = effect_direction
    )
    mu_samples <- mu_samples + study_contribution
  }

  ### return study-level predictions if type = "study" is selected
  # study incorporates fixed effects + study-level random effects (mu + gamma)
  if (type == "study") {
    # rename samples
    colnames(mu_samples) <- if (aggregate) "mu_study" else paste0("mu_study[", seq_len(K), "]")

    return(.new_brma_samples(
      samples  = mu_samples,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = if (aggregate) "Aggregated Study-Level Posterior Prediction:" else "Study-Level Posterior Prediction:",
      probs    = probs,
      data     = if (aggregate) NULL else new_data
    ))
  }

  ### create true-effect prediction
  # dispatches between normal and GLMM approaches
  # both helpers handle same_data dispatch internally:
  # - same_data = TRUE: use fitted values (BLUPs for normal, extracted theta for GLMM)
  # - same_data = FALSE: sample from marginal distribution N(mu, tau_within)
  if (type == "estimate") {

    if (outcome_type == "norm") {
      true_effects_samples <- .evaluate.brma.true_effects.norm(
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        yi         = outcome_data[["yi"]],
        sei        = outcome_data[["sei"]],
        same_data  = same_data
      )
    } else {
      true_effects_samples <- .evaluate.brma.true_effects.glmm(
        fit        = object[["fit"]],
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        same_data  = same_data,
        K          = K
      )
    }

    # rename samples
    colnames(true_effects_samples) <- if (aggregate) "theta" else paste0("theta[", seq_len(K), "]")

    return(.new_brma_samples(
      samples  = true_effects_samples,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = if (aggregate) "Aggregated True Effect Posterior Prediction:" else "True Effect Posterior Prediction:",
      probs    = probs,
      data     = if (aggregate) NULL else new_data
    ))

  }

  ### create observed effects prediction
  if (type == "response") {

    # different model types have different output structures
    if (is.element(outcome_type, c("bin", "pois"))) {

      # the estimate-level effects are not marginalized for GLMMs
      # include the sampled random effects or marginalize over them for new data
      # see .evaluate.brma.theta.glmm() in brma.evaluate.R for details
      theta_contribution <- .evaluate.brma.theta.glmm(
        fit        = object[["fit"]],
        tau_within = tau_within_samples,
        same_data  = same_data,
        K          = K
      )
      mu_samples <- mu_samples + theta_contribution

      if (outcome_type == "bin") {

        ### incorporate with base-rate using helper function
        # .evaluate.brma.baserate() returns S x K matrix of logit(pi) samples
        logit_baserate <- .evaluate.brma.baserate(fit = object[["fit"]], K = K)

        ### sample outcome using RNG helper
        outcome_samples <- .outcome_rng.binom(
          mu_samples     = mu_samples,
          logit_baserate = logit_baserate,
          n1i            = outcome_data[["n1i"]],
          n2i            = outcome_data[["n2i"]]
        )

      } else if (outcome_type == "pois") {

        ### incorporate with log-rate using helper function
        # .evaluate.brma.lograte() returns S x K matrix of log(phi) samples
        log_phi <- .evaluate.brma.lograte(fit = object[["fit"]], K = K)

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

        if (outcome_type == "bin") {
          # convert counts to log-odds ratio using escalc-like formula
          # logOR = log((ai/bi) / (ci/di)) where bi = n1i - ai, di = n2i - ci
          # extract ai and ci from interleaved outcome_samples
          ai <- outcome_samples[, seq(1, 2 * K, by = 2), drop = FALSE]
          ci <- outcome_samples[, seq(2, 2 * K, by = 2), drop = FALSE]
          # compute complementary counts
          n1i_mat <- matrix(outcome_data[["n1i"]], nrow = nrow(ai), ncol = K, byrow = TRUE)
          n2i_mat <- matrix(outcome_data[["n2i"]], nrow = nrow(ai), ncol = K, byrow = TRUE)
          bi <- n1i_mat - ai
          di <- n2i_mat - ci
          # apply 0.5 continuity correction for zero cells
          ai_adj <- ai + 0.5 * (ai == 0 | bi == 0 | ci == 0 | di == 0)
          bi_adj <- bi + 0.5 * (ai == 0 | bi == 0 | ci == 0 | di == 0)
          ci_adj <- ci + 0.5 * (ai == 0 | bi == 0 | ci == 0 | di == 0)
          di_adj <- di + 0.5 * (ai == 0 | bi == 0 | ci == 0 | di == 0)
          # compute log-odds ratio
          outcome_samples <- log((ai_adj * di_adj) / (bi_adj * ci_adj))
          colnames(outcome_samples) <- paste0("yi[", seq_len(K), "]")

        } else if (outcome_type == "pois") {
          # convert counts to log incidence rate ratio
          # logIRR = log((x1i/t1i) / (x2i/t2i))
          x1i <- outcome_samples[, seq(1, 2 * K, by = 2), drop = FALSE]
          x2i <- outcome_samples[, seq(2, 2 * K, by = 2), drop = FALSE]
          t1i_mat <- matrix(outcome_data[["t1i"]], nrow = nrow(x1i), ncol = K, byrow = TRUE)
          t2i_mat <- matrix(outcome_data[["t2i"]], nrow = nrow(x1i), ncol = K, byrow = TRUE)
          # apply 0.5 continuity correction for zero counts
          x1i_adj <- x1i + 0.5 * (x1i == 0 | x2i == 0)
          x2i_adj <- x2i + 0.5 * (x1i == 0 | x2i == 0)
          # compute log incidence rate ratio
          outcome_samples <- log((x1i_adj / t1i_mat) / (x2i_adj / t2i_mat))
          colnames(outcome_samples) <- paste0("yi[", seq_len(K), "]")
        }

      }


    } else if (outcome_type == "norm") {

      # normal outcome models dispatch between:
      # - bias_adjusted = TRUE: sample from unweighted normal (as if no bias)
      # - bias_adjusted = FALSE with weightfunction: sample from weighted normal
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

        # sample from weighted normal distribution for selection models
        # extract omega samples for weight function
        posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
        omega_samples     <- posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE]

        # for weighted distributions, crit_yi is computed in "positive" space
        # (yi flipped for negative effect direction in .create_fit_data)
        # so we need to flip mu_samples to match, sample, then flip back
        if (effect_direction == "negative") {
          # flip mu to positive space for sampling
          mu_samples_for_wnorm <- -mu_samples
        } else {
          mu_samples_for_wnorm <- mu_samples
        }

        outcome_samples <- .outcome_rng.wnorm(
          mu_samples = mu_samples_for_wnorm,
          tau_within = tau_within_samples,
          sei        = outcome_data[["sei"]],
          omega      = omega_samples,
          crit_yi    = fit_data$crit_yi
        )

        # flip samples back to original space
        if (effect_direction == "negative") {
          outcome_samples <- -outcome_samples
        }

      }

      # rename samples
      colnames(outcome_samples) <- paste0("yi[", seq_len(K), "]")
    }

    return(.new_brma_samples(
      samples  = outcome_samples,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = "Observations Posterior Prediction:",
      probs    = probs,
      data     = new_data
    ))
  }
}
