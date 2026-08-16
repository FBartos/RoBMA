#' @title Predict From brma Object
#'
#' @description \code{predict.brma} predicts values
#'
#' @inheritParams summary.brma
#' @param probs quantiles of the posterior samples to be displayed when the
#' returned `brma_samples` object is printed. Defaults to `c(.025, .975)`.
#' @param newdata prediction data. Defaults to \code{NULL}, which predicts the
#'   observed design rows. Marginal new-effect predictions may use an explicit data frame or
#'   named list with one row per target effect and every variable used by the
#'   requested fixed, scale, and random-effect designs. Missing ordinary
#'   grouping-only identifiers are synthesized as row-specific new levels;
#'   known-\code{R} grouping identifiers must be supplied. Outcome columns are
#'   optional for \code{type = "terms"} unless PET/PEESE bias terms are included
#'   via \code{bias_adjusted = FALSE}. Missing values are rejected; prediction
#'   never drops rows silently. Binomial GLMM response predictions require arm
#'   totals \code{n1i} and \code{n2i}; supply the totals alone, optionally with
#'   event counts \code{ai} and \code{ci}, or supply the four cell counts
#'   \code{ai}, \code{bi}, \code{ci}, and \code{di}. Supplied cells and totals
#'   must agree. Poisson GLMM response predictions require exposures \code{t1i}
#'   and \code{t2i}; event counts \code{x1i} and \code{x2i} are optional. For
#'   known-\code{V} \code{brma.mv()} response prediction, explicit
#'   \code{newdata} also requires a new sampling covariance matrix.
#' @param type type of prediction to be performed. Options are:
#' \itemize{
#'   \item{\code{"terms"} (alias: \code{"marginal"}): fixed location
#'   parameters only, \eqn{X\beta}.}
#'   \item{\code{"estimate"} (alias: \code{"effect"}): latent true effects.
#'   The random-effect distribution is selected by \code{conditioning_depth}.}
#'   \item{\code{"response"} (alias: \code{"outcome"}): Predicted observed values (yi).
#'   Adds the outcome sampling distribution to the corresponding latent-effect
#'   target.}
#'   \item{\code{"terms.scale"}: Scale parameter (tau), incorporating scale
#'   regression if present. For \code{brma.mv()} random-formula models with
#'   scale formulas, returns a named list of component-specific
#'   \code{brma_samples} matrices.}
#' }
#' The released \code{"cluster"} and \code{"blup"} spellings are retained as
#' compatibility shortcuts for conditional location samples: they fix
#' \code{conditioning_depth = "cluster"} and \code{"estimate"}, respectively.
#' Use \code{fitted()} or \code{blup()} when conditional means are the intended
#' quantity.
#' @param conditioning_depth random-effect conditioning depth for
#'   \code{type = "estimate"} and \code{"response"}. \code{"marginal"}
#'   (default) predicts new latent effects, \code{"cluster"} conditions on the
#'   fitted cluster effect and predicts a new estimate within that cluster, and
#'   \code{"estimate"} predicts the fitted latent effect. Cluster depth is
#'   available only for legacy multilevel models. Non-marginal prediction is
#'   currently restricted to \code{newdata = NULL}, where fitted identities are
#'   unambiguous.
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
#' conditioned on the heterogeneity component. Conditional samples are flattened
#' to one chain after subsetting posterior rows.
#' @param quiet logical; whether to suppress informational messages about
#' prediction scale and bias adjustment.
#' @param V_new optional sampling covariance matrix for explicit
#' \code{newdata} response predictions from known-\code{V} \code{brma.mv()}
#' models. May be a square matrix or a list of block covariance matrices whose
#' total dimension matches \code{nrow(newdata)}. Cross-covariance with observed
#' rows is not supported.
#'
#' @details
#' Prediction has two independent axes. \code{type} selects the quantity:
#' fixed location (\code{"terms"}), latent true effect (\code{"estimate"}), or
#' observed outcome (\code{"response"}). \code{conditioning_depth} selects how
#' much of the fitted random-effect hierarchy is retained. In the legacy model
#' \eqn{y_{ij} = X_{ij}\beta + u_j + v_{ij} + \epsilon_{ij}}, the targets are:
#' \itemize{
#'   \item{\code{"marginal"}: condition on neither \eqn{u_j} nor \eqn{v_{ij}};}
#'   \item{\code{"cluster"}: condition on fitted \eqn{u_j}, but draw a new
#'   \eqn{v_{ij}};}
#'   \item{\code{"estimate"}: draw the fitted latent effect from its posterior.}
#' }
#' \code{type = "response"} adds \eqn{\epsilon} to the selected latent-effect
#' target. For normal fitted effects this includes the conditional latent
#' uncertainty, not only the BLUP mean. \code{blup()} and \code{fitted(...,
#' conditioning_depth = "estimate")} remain the conditional-mean interfaces.
#'
#' \code{newdata} supplies design values and grouping identities; it does not
#' select an estimand. Consequently, marginal prediction has the same law for
#' \code{newdata = NULL} and an equivalent explicit design. Matching fitted
#' grouping labels retains the requested dependence among new marginal draws
#' but does not reuse fitted random effects. Cluster- and estimate-depth
#' predictions require fitted identities and therefore currently reject
#' explicit \code{newdata}.
#'
#' For binomial and Poisson GLMM response predictions, estimate depth uses the
#' fitted posterior nuisance base rate \eqn{\pi_i} or rate \eqn{\phi_i}.
#' Marginal and cluster depths describe a new estimate and therefore draw a new
#' nuisance rate independently from its model prior. Those predictions contain
#' a prior-predictive nuisance-rate component, which can dominate their width
#' when the nuisance prior is diffuse; fitted \eqn{\pi_i}/\eqn{\phi_i} values are
#' never silently reused for a marginal or cluster target.
#'
#' For RoBMA product-space objects, conditional posterior predictions subset
#' posterior rows according to model indicators. This removes the original
#' chain structure, so returned \code{brma_samples} objects are intentionally
#' stored as one flattened chain.
#'
#' Note that in contrast to \link[metafor]{predict}, the \code{type = "response"} produces
#' predictions for the new effect size estimates. To obtain results corresponding to
#' metafor's predict function, use \code{type = "terms"} for the mean effect size
#' and \code{type = "estimate"} for true effects (prediction interval).
#'
#' \strong{Likelihood weights:} If the model was fitted with \code{weights},
#' the weights affect the posterior fit, log-likelihood, LOO/WAIC, and
#' existing-data conditional diagnostics such as BLUP shrinkage and leverage.
#' They do not change the observation-level sampling error used by
#' \code{type = "response"} for normal models: response predictions simulate
#' raw future effect-size estimates with the supplied \code{sei}, not
#' \code{sei / sqrt(weight)}.
#'
#' For known-\code{V} \code{brma.mv()} models, explicit \code{newdata}
#' response prediction requires \code{V_new}. Random-formula models preserve
#' the joint random-effect structure selected by \code{conditioning_depth} and
#' then add \code{V} or \code{V_new} sampling noise. Cross-covariance with
#' observed rows is not supported. Explicit
#' \code{newdata} estimate/response predictions are not
#' available for \code{brma.mv()} objects with marginalized known-\code{R}
#' blocks, because their row multipliers are fitted-row metadata.
#' \code{brma.mv()} prediction results carry a
#' \code{brma_mv_prediction_target} attribute recording the formula,
#' random-effect, and response-generation target used.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE) &&
#'     requireNamespace("metafor", quietly = TRUE)) {
#'   data(dat.bcg, package = "metadat")
#'   dat <- metafor::escalc(
#'     measure = "RR",
#'     ai      = tpos,
#'     bi      = tneg,
#'     ci      = cpos,
#'     di      = cneg,
#'     data    = dat.bcg
#'   )
#'
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     mods    = ~ ablat + year,
#'     data    = dat,
#'     measure = "RR",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   predict(fit, type = "terms")
#'   predict(fit, newdata = data.frame(ablat = 50, year = 1970),
#'           type = "estimate")
#'   predict(fit, type = "estimate")
#' }
#' }
#'
#' @return A \code{brma_samples} object containing posterior samples. For
#' \code{brma.mv()} random-formula \code{type = "terms.scale"} predictions,
#' returns a named list of component-specific \code{brma_samples} objects.
#' When printed, each \code{brma_samples} object displays a summary table via
#' \code{BayesTools::ensemble_estimates_table}. The underlying samples matrix
#' can be accessed directly because the object inherits from matrix. Use
#' \code{summary()} or \code{as.data.frame()} to obtain the summary table. For
#' multi-component results, the data frame contains component and parameter
#' identifiers; use \code{as.data.frame(format = "list")} to retain separate
#' tables. The samples can also be converted to \pkg{posterior} draws formats
#' using \code{as_draws()} and related functions.
#' @seealso [pooled_effect()], [pooled_heterogeneity()], [blup()]
#' @export
predict.brma <- function(object, newdata = NULL, type = "terms",
                         as_measure = TRUE,
                         output_measure = NULL,
                         transform = NULL,
                         probs = c(.025, .975),
                         bias_adjusted = FALSE,
                         quiet = FALSE,
                         conditional = FALSE,
                         V_new = NULL,
                         conditioning_depth = "marginal",
                         ...){

  conditioning_depth_specified <- !missing(conditioning_depth)

  context <- .predict_brma_context(
    object                       = object,
    newdata                      = newdata,
    V_new                        = V_new,
    type                         = type,
    conditioning_depth           = conditioning_depth,
    conditioning_depth_specified = conditioning_depth_specified,
    as_measure                   = as_measure,
    output_measure               = output_measure,
    transform                    = transform,
    probs                        = probs,
    bias_adjusted                = bias_adjusted,
    quiet                        = quiet,
    conditional                  = conditional,
    dots                         = list(...)
  )

  type <- context[["type"]]

  if (type == "terms") {
    return(.predict_brma_terms(context))
  }

  scale_state <- .predict_brma_scale_state(context)
  if (type == "terms.scale") {
    return(.predict_brma_scale(context, scale_state))
  }
  location_state <- .predict_brma_location_state(context, scale_state)

  if (type == "location") {
    mu_samples <- location_state[["mu"]]
    prefix <- if (context[["conditioning_depth"]] == "cluster") {
      "mu_cluster"
    } else {
      "theta_blup"
    }
    colnames(mu_samples) <- paste0(prefix, "[", seq_len(context[["K"]]), "]")

    return(.predict_brma_finalize(
      context    = context,
      samples    = mu_samples,
      title      = if (context[["conditioning_depth"]] == "cluster") {
        "Cluster-Level Fitted Location:"
      } else {
        "True Effects (BLUP Means):"
      },
      parameters = .conditional_effect_parameters(context[["object"]])
    ))
  }

  if (type == "estimate") {
    return(.predict_brma_estimate(context, location_state, scale_state))
  }

  ### create observed effects prediction
  if (type == "response") {
    return(.predict_brma_response(
      context        = context,
      location_state = location_state,
      scale_state    = scale_state
    ))
  }
}


.predict_brma_context <- function(
    object, newdata, V_new, type, conditioning_depth,
    conditioning_depth_specified, as_measure, output_measure, transform,
    probs, bias_adjusted, quiet, conditional, dots) {

  .check_unused_dots(
    dots    = dots,
    allowed = ".posterior_samples",
    caller  = "predict.brma()"
  )

  requested_type <- type
  type <- match.arg(type, c("terms", "marginal", "cluster", "estimate", "effect", "blup",
                            "response", "outcome", "terms.scale"))
  type <- switch(type,
    "marginal" = "terms",
    "effect"   = "estimate",
    "cluster"  = "location",
    "blup"     = "location",
    "outcome"  = "response",
    type  # default: keep as is
  )
  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)
  if (requested_type == "cluster") {
    if (conditioning_depth_specified && conditioning_depth != "cluster") {
      stop(
        "type = 'cluster' fixes conditioning_depth = 'cluster'.",
        call. = FALSE
      )
    }
    conditioning_depth <- "cluster"
  }
  if (requested_type == "blup") {
    if (conditioning_depth_specified && conditioning_depth != "estimate") {
      stop(
        "type = 'blup' fixes conditioning_depth = 'estimate'.",
        call. = FALSE
      )
    }
    conditioning_depth <- "estimate"
  }

  # input validation
  BayesTools::check_bool(as_measure, "as_measure")
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(quiet, "quiet")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' predictions are available only for RoBMA objects.",
         call. = FALSE)
  }

  is_brma_mv_object <- inherits(object, "brma.mv")
  is_random_object  <- .is_random(object)
  # check newdata: NULL or an explicit data.frame/list
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata) && !is.list(newdata)) {
      stop("'newdata' must be NULL, a data.frame, or a named list.", call. = FALSE)
    }
    if (is.list(newdata) && !is.data.frame(newdata) && is.null(names(newdata))) {
      stop("'newdata' list must be named.", call. = FALSE)
    }
  }

  if (type %in% c("terms", "terms.scale") &&
      conditioning_depth_specified && conditioning_depth != "marginal") {
    stop(
      "'conditioning_depth' is only available for type = 'estimate' or ",
      "type = 'response'.",
      call. = FALSE
    )
  }
  if (!is.null(newdata) && conditioning_depth != "marginal" &&
      type %in% c("location", "estimate", "response")) {
    stop(
      "Non-marginal predictions require fitted observation identities and ",
      "are currently available only with newdata = NULL.",
      call. = FALSE
    )
  }

  known_v_newdata_response <- .is_data_known_v(object[["data"]]) &&
    !is.null(newdata) && type == "response"
  known_V_new <- NULL
  if (!is.null(V_new)) {
    if (!known_v_newdata_response) {
      stop(
        "'V_new' is only available for known-V brma.mv() response ",
        "predictions with explicit 'newdata'.",
        call. = FALSE
      )
    }
    newdata_frame <- .prepare_newdata_as_data_frame(newdata)
    known_V_new   <- .known_v_newdata_prepare(
      V_new,
      k = nrow(newdata_frame)
    )
    newdata <- .predict_known_v_newdata_add_vi(
      newdata = newdata_frame,
      V_new   = .known_v_diagonal(known_V_new)
    )
  }

  # check incompatible options
  if (.is_data_known_v(object[["data"]]) && !is.null(newdata) &&
      type == "response" && is.null(known_V_new)) {
    stop(
      "Newdata response predictions for brma.mv() known-V models require ",
      "a supplied 'V_new' sampling covariance matrix.",
      call. = FALSE
    )
  }
  if (is_brma_mv_object &&
      !is.null(newdata) &&
      type %in% c("estimate", "response") &&
      .data_has_marginalized_known_group_covariance(object[["data"]])) {
    stop(
      "Explicit prediction rows are not available for marginalized known-R ",
      "blocks because their row multipliers are tied to the fitted layout. ",
      "Support for a supplied 'R_new' is deferred.",
      call. = FALSE
    )
  }

  ### types of predictions
  # terms:       fixed location
  # location:    conditional location mean selected by conditioning depth
  # terms.scale: fitted heterogeneity parameters
  # estimate:    latent true-effect draws selected by conditioning depth
  # response:    estimate target plus outcome sampling variability

  ### dispatch between prediction on the current data vs. new data
  if (is.null(newdata)) {

    # existing data are used
    same_data <- TRUE
    new_data  <- object[["data"]]

  } else {

    # prepare newdata using the same settings as the original fit
    same_data <- FALSE
    new_data  <- .prepare_newdata(
      object        = object,
      newdata       = newdata,
      type          = type,
      bias_adjusted  = bias_adjusted,
      include_scale  = type != "terms",
      include_random = is_random_object && is_brma_mv_object &&
        type %in% c("estimate", "response", "terms.scale")
    )

  }

  ### extract priors and structural information about the model
  priors            <- object[["priors"]]
  is_mods           <- .is_mods(object)
  is_multilevel     <- .is_multilevel(object)
  is_scale          <- .is_scale(object)
  is_random         <- is_random_object
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)
  is_weights        <- .is_weights(object)
  is_known_v        <- .is_data_known_v(new_data) || !is.null(known_V_new)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)

  if (type == "response" && outcome_type == "norm" && is_known_v &&
      is_weightfunction && !bias_adjusted) {
    stop(
      "Bias-unadjusted response predictions for known-V weightfunction ",
      "models are not supported because selected-normal sampling does not ",
      "support known-V covariance.",
      call. = FALSE
    )
  }

  if (conditioning_depth == "cluster" && !is_multilevel) {
    stop(
      "conditioning_depth = 'cluster' is only available for legacy ",
      "multilevel (3-level) models.",
      call. = FALSE
    )
  }
  if (conditioning_depth == "cluster" && is_brma_mv_object && is_random_object) {
    stop(
      "conditioning_depth = 'cluster' is unavailable for brma.mv() ",
      "random-formula models because their hierarchy need not define one ",
      "distinguished cluster level.",
      call. = FALSE
    )
  }

  posterior_samples <- .get_posterior_samples(object[["fit"]], dots[[".posterior_samples"]])
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

  ### extract outcome data and fit data for convenience
  outcome_data      <- new_data[["outcome"]]
  fit_data          <- .create_fit_data(data = new_data, priors = priors)

  # outcome dimensions
  K_original <- nrow(outcome_data)
  K          <- K_original

  ### extract MCMC chain info for brma_samples construction
  chain_info <- .brma_samples_chain_info(
    fit       = object[["fit"]],
    n_samples = nrow(posterior_samples)
  )
  n_chains <- chain_info[["n_chains"]]
  n_iter   <- chain_info[["n_iter"]]

  random_mv <- is_brma_mv_object && is_random

  list(
    object              = object,
    type                = type,
    requested_type      = requested_type,
    conditioning_depth  = conditioning_depth,
    as_measure          = as_measure,
    probs               = probs,
    bias_adjusted       = bias_adjusted,
    quiet               = quiet,
    conditional         = conditional,
    same_data           = same_data,
    new_data            = new_data,
    known_V_new         = known_V_new,
    priors              = priors,
    is_mods             = is_mods,
    is_multilevel       = is_multilevel,
    is_scale            = is_scale,
    is_random           = is_random,
    is_PET              = is_PET,
    is_PEESE            = is_PEESE,
    is_weightfunction   = is_weightfunction,
    is_weights          = is_weights,
    is_known_v          = is_known_v,
    outcome_type        = outcome_type,
    effect_direction    = effect_direction,
    posterior_samples   = posterior_samples,
    effect_transform    = effect_transform,
    outcome_data        = outcome_data,
    fit_data            = fit_data,
    K_original          = K_original,
    K                   = K,
    n_chains            = n_chains,
    n_iter              = n_iter,
    random_mv           = random_mv
  )
}


.predict_brma_finalize <- function(
    context, samples, title, parameters, effect = TRUE) {

  constructor <- if (effect) .new_effect_brma_samples else .new_brma_samples
  arguments   <- list(
    samples  = samples,
    n_chains = context[["n_chains"]],
    n_iter   = context[["n_iter"]],
    title    = title,
    probs    = context[["probs"]],
    data     = context[["new_data"]]
  )
  if (effect) {
    arguments[["effect_transform"]] <- context[["effect_transform"]]
  }
  out <- do.call(constructor, arguments)
  out <- .predict_brma_attach_mv_metadata(
    samples            = out,
    object             = context[["object"]],
    type               = context[["type"]],
    conditioning_depth = context[["conditioning_depth"]],
    same_data          = context[["same_data"]],
    random_mv          = context[["random_mv"]],
    known_V_new        = context[["known_V_new"]]
  )

  .condition_prediction_samples(
    object            = context[["object"]],
    samples           = out,
    conditional       = context[["conditional"]],
    parameters        = parameters,
    posterior_samples = context[["posterior_samples"]],
    quiet             = context[["quiet"]]
  )
}


.predict_brma_terms <- function(context) {

  object            <- context[["object"]]
  new_data          <- context[["new_data"]]
  priors            <- context[["priors"]]
  posterior_samples <- context[["posterior_samples"]]
  K                 <- context[["K"]]

  mu_samples <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = context[["outcome_data"]],
    mods_data         = new_data[["mods"]],
    mods_formula      = if (context[["is_mods"]]) {
      .create_fit_formula_list(data = new_data, "mods")
    } else {
      NULL
    },
    mods_priors       = if (context[["is_random"]]) {
      priors[["location"]]
    } else {
      priors[["mods"]]
    },
    is_mods           = context[["is_mods"]],
    is_PET            = context[["is_PET"]],
    is_PEESE          = context[["is_PEESE"]],
    effect_direction  = context[["effect_direction"]],
    bias_adjusted     = context[["bias_adjusted"]],
    K                 = context[["K_original"]],
    posterior_samples = posterior_samples,
    priors            = priors
  )
  colnames(mu_samples) <- paste0("mu[", seq_len(K), "]")

  .predict_brma_finalize(
    context    = context,
    samples    = mu_samples,
    title      = "Location Term Posterior Prediction:",
    parameters = .conditional_effect_parameters(object)
  )
}


.predict_brma_scale_state <- function(context) {

  posterior_samples <- context[["posterior_samples"]]
  K                  <- context[["K_original"]]
  if (context[["random_mv"]]) {
    return(list(
      within  = matrix(0, nrow = nrow(posterior_samples), ncol = K),
      between = matrix(0, nrow = nrow(posterior_samples), ncol = K)
    ))
  }

  new_data <- context[["new_data"]]
  priors   <- context[["priors"]]
  result   <- .evaluate.brma.tau(
    fit               = context[["object"]][["fit"]],
    scale_data        = new_data[["scale"]],
    scale_formula     = if (context[["is_scale"]]) {
      .create_fit_formula_list(data = new_data, "scale")
    } else {
      NULL
    },
    scale_priors      = priors[["scale"]],
    is_scale          = context[["is_scale"]],
    is_multilevel     = context[["is_multilevel"]],
    K                 = K,
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(priors),
    fixed_rho         = .fixed_rho_prior_value(priors)
  )

  list(
    within  = result[["tau_within"]],
    between = result[["tau_between"]]
  )
}


.predict_brma_scale <- function(context, scale_state) {

  object            <- context[["object"]]
  new_data          <- context[["new_data"]]
  posterior_samples <- context[["posterior_samples"]]
  if (context[["random_mv"]]) {
    scale_samples <- .brma_mv_random_heterogeneity_components(
      object            = object,
      posterior_samples = posterior_samples,
      K                 = context[["K_original"]],
      data              = new_data,
      new_levels        = if (context[["same_data"]]) NULL else "sample"
    )
    scale_samples <- lapply(scale_samples, function(samples) {
      colnames(samples) <- paste0("tau[", seq_len(ncol(samples)), "]")
      samples
    })
    out <- lapply(names(scale_samples), function(component) {
      .new_brma_samples(
        samples  = scale_samples[[component]],
        n_chains = context[["n_chains"]],
        n_iter   = context[["n_iter"]],
        title    = paste0("Scale Term Posterior Prediction (", component, "):"),
        probs    = context[["probs"]],
        data     = new_data
      )
    })
    names(out) <- names(scale_samples)
    out <- .predict_brma_attach_mv_metadata(
      samples            = out,
      object             = object,
      type               = context[["type"]],
      conditioning_depth = context[["conditioning_depth"]],
      same_data          = context[["same_data"]],
      random_mv          = TRUE,
      known_V_new        = context[["known_V_new"]]
    )
    out <- .new_brma_samples_list(out)
    return(.condition_prediction_samples(
      object            = object,
      samples           = out,
      conditional       = context[["conditional"]],
      parameters        = .conditional_heterogeneity_parameters(object),
      posterior_samples = posterior_samples,
      quiet             = context[["quiet"]]
    ))
  }

  tau_samples <- .root_sum_squares(
    scale_state[["within"]],
    scale_state[["between"]]
  )
  colnames(tau_samples) <- paste0("tau[", seq_len(context[["K"]]), "]")
  .predict_brma_finalize(
    context    = context,
    samples    = tau_samples,
    title      = "Scale Term Posterior Prediction:",
    parameters = .conditional_heterogeneity_parameters(object),
    effect     = FALSE
  )
}


.predict_brma_location_state <- function(context, scale_state) {

  object              <- context[["object"]]
  conditioning_depth  <- context[["conditioning_depth"]]
  same_data           <- context[["same_data"]]
  new_data            <- context[["new_data"]]
  known_V_new         <- context[["known_V_new"]]
  priors              <- context[["priors"]]
  is_mods             <- context[["is_mods"]]
  is_multilevel       <- context[["is_multilevel"]]
  is_random           <- context[["is_random"]]
  is_PET              <- context[["is_PET"]]
  is_PEESE            <- context[["is_PEESE"]]
  is_weightfunction   <- context[["is_weightfunction"]]
  is_weights          <- context[["is_weights"]]
  is_known_v          <- context[["is_known_v"]]
  outcome_type        <- context[["outcome_type"]]
  effect_direction    <- context[["effect_direction"]]
  posterior_samples   <- context[["posterior_samples"]]
  outcome_data        <- context[["outcome_data"]]
  fit_data            <- context[["fit_data"]]
  K_original          <- context[["K_original"]]
  K                   <- context[["K"]]
  random_mv           <- context[["random_mv"]]
  bias_adjusted       <- context[["bias_adjusted"]]
  tau_within_samples  <- scale_state[["within"]]
  tau_between_samples <- scale_state[["between"]]

  fixed_mu <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = outcome_data,
    mods_data         = new_data[["mods"]],
    mods_formula      = if (is_mods) {
      .create_fit_formula_list(data = new_data, "mods")
    } else {
      NULL
    },
    mods_priors       = if (is_random) priors[["location"]] else priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = bias_adjusted,
    K                 = K_original,
    posterior_samples = posterior_samples,
    priors            = priors
  )

  blup_vi <- NULL
  if (outcome_type == "norm") {
    blup_vi <- outcome_data[["sei"]]^2
    if (!is.null(known_V_new)) {
      blup_vi <- .known_v_residual_variance(known_V_new)
    } else if (is_known_v) {
      blup_vi <- .data_known_v_data(new_data)[["residual_variance"]]
    }
    if (is_weights) {
      blup_vi <- blup_vi / outcome_data[["weights"]]
    }
  }

  blup_bias_offset <- NULL
  if (outcome_type == "norm" && same_data && bias_adjusted &&
      (is_PET || is_PEESE)) {
    blup_bias_offset <- .evaluate.brma.bias_offset(
      fit               = object[["fit"]],
      outcome_data      = outcome_data,
      is_PET            = is_PET,
      is_PEESE          = is_PEESE,
      effect_direction  = effect_direction,
      K                 = K,
      posterior_samples = posterior_samples,
      priors            = priors
    )
  }

  mu_samples       <- fixed_mu
  cluster_mu       <- fixed_mu
  multilevel_blup  <- NULL
  use_known_v_blup <- outcome_type == "norm" && same_data && is_known_v &&
    !random_mv

  needs_sampled_cluster <- is_multilevel &&
    (conditioning_depth == "cluster" ||
     (conditioning_depth == "estimate" &&
      (outcome_type != "norm" || is_weightfunction)))
  if (needs_sampled_cluster) {
    cluster_contribution <- .evaluate.brma.cluster_effects(
      fit               = object[["fit"]],
      tau_between       = tau_between_samples,
      cluster           = fit_data[["cluster"]],
      same_data         = TRUE,
      effect_direction  = effect_direction,
      posterior_samples = posterior_samples
    )
    cluster_mu <- fixed_mu + cluster_contribution
  }
  if (conditioning_depth == "cluster") {
    mu_samples <- cluster_mu
  }

  if (conditioning_depth == "estimate") {
    if (outcome_type == "norm" && context[["type"]] != "location") {
      mu_samples <- cluster_mu
    } else if (random_mv && outcome_type == "norm" && same_data && is_known_v) {
      random_blup <- .evaluate.brma.mv_random_blup.norm(
        object            = object,
        mu_samples        = fixed_mu,
        posterior_samples = posterior_samples,
        bias_offset       = blup_bias_offset
      )
      mu_samples <- fixed_mu + random_blup
    } else if (is_multilevel && outcome_type == "norm" && same_data &&
               !is_weightfunction) {
      multilevel_blup <- .evaluate.brma.multilevel_blup.norm(
        mu_samples  = fixed_mu,
        tau_within  = tau_within_samples,
        tau_between = tau_between_samples,
        yi          = outcome_data[["yi"]],
        vi          = blup_vi,
        cluster     = fit_data[["cluster"]],
        bias_offset = blup_bias_offset
      )
      mu_samples <- fixed_mu + multilevel_blup[["cluster"]] +
        multilevel_blup[["estimate"]]
    } else {
      mu_samples <- cluster_mu

      if (outcome_type == "norm" && use_known_v_blup) {
        mu_samples <- .evaluate.brma.known_v_blup.norm(
          mu_samples  = mu_samples,
          tau_within  = tau_within_samples,
          yi          = outcome_data[["yi"]],
          known_V     = .data_known_v_data(new_data),
          bias_offset = blup_bias_offset
        )
      } else if (outcome_type == "norm") {
        mu_samples <- .evaluate.brma.true_effects.norm(
          mu_samples  = mu_samples,
          tau_within  = tau_within_samples,
          yi          = outcome_data[["yi"]],
          sei         = sqrt(blup_vi),
          same_data   = TRUE,
          bias_offset = blup_bias_offset
        )
      } else {
        mu_samples <- .evaluate.brma.true_effects.glmm(
          fit               = object[["fit"]],
          mu_samples        = mu_samples,
          tau_within        = tau_within_samples,
          same_data         = TRUE,
          K                 = K,
          posterior_samples = posterior_samples
        )
      }
    }
  }


  list(
    mu               = mu_samples,
    fixed_mu          = fixed_mu,
    cluster_mu        = cluster_mu,
    blup_vi          = blup_vi,
    blup_bias_offset = blup_bias_offset,
    multilevel_blup  = multilevel_blup,
    use_known_v_blup = use_known_v_blup
  )
}


.predict_brma_estimate <- function(context, location_state, scale_state) {

  true_effects_samples <- .predict_brma_estimate_draws(
    context        = context,
    location_state = location_state,
    scale_state    = scale_state
  )
  colnames(true_effects_samples) <- paste0(
    "theta[", seq_len(context[["K"]]), "]"
  )

  .predict_brma_finalize(
    context    = context,
    samples    = true_effects_samples,
    title      = "True Effect Posterior Prediction:",
    parameters = .conditional_effect_parameters(context[["object"]])
  )
}


.predict_brma_estimate_draws <- function(context, location_state, scale_state) {

  object               <- context[["object"]]
  conditioning_depth   <- context[["conditioning_depth"]]
  new_data             <- context[["new_data"]]
  outcome_data         <- context[["outcome_data"]]
  outcome_type         <- context[["outcome_type"]]
  posterior_samples    <- context[["posterior_samples"]]
  K                    <- context[["K"]]
  random_mv            <- context[["random_mv"]]
  mu_samples           <- location_state[["mu"]]
  fixed_mu             <- location_state[["fixed_mu"]]
  blup_vi              <- location_state[["blup_vi"]]
  blup_bias_offset     <- location_state[["blup_bias_offset"]]
  tau_within_samples   <- scale_state[["within"]]
  tau_between_samples  <- scale_state[["between"]]

  if (conditioning_depth == "marginal") {
    if (random_mv) {
      true_effects_samples <- fixed_mu +
        .predict_brma_mv_new_effect_random_draws(
          object            = object,
          data              = new_data,
          posterior_samples = posterior_samples
        )
    } else {
      effect_mean <- fixed_mu
      if (context[["is_multilevel"]]) {
        effect_mean <- effect_mean + .evaluate.brma.cluster_effects(
          fit               = object[["fit"]],
          tau_between       = tau_between_samples,
          cluster           = context[["fit_data"]][["cluster"]],
          same_data         = FALSE,
          effect_direction  = context[["effect_direction"]],
          posterior_samples = posterior_samples
        )
      }
      true_effects_samples <- if (outcome_type == "norm") {
        .evaluate.brma.true_effects.norm(
          mu_samples  = effect_mean,
          tau_within  = tau_within_samples,
          yi          = outcome_data[["yi"]],
          sei         = sqrt(blup_vi),
          same_data   = FALSE
        )
      } else {
        .evaluate.brma.true_effects.glmm(
          fit               = object[["fit"]],
          mu_samples        = effect_mean,
          tau_within        = tau_within_samples,
          same_data         = FALSE,
          K                 = K,
          posterior_samples = posterior_samples
        )
      }
    }
  } else if (conditioning_depth == "cluster") {
    true_effects_samples <- if (outcome_type == "norm") {
      .evaluate.brma.true_effects.norm(
        mu_samples  = mu_samples,
        tau_within  = tau_within_samples,
        yi          = outcome_data[["yi"]],
        sei         = sqrt(blup_vi),
        same_data   = FALSE
      )
    } else {
      .evaluate.brma.true_effects.glmm(
        fit               = object[["fit"]],
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        same_data         = FALSE,
        K                 = K,
        posterior_samples = posterior_samples
      )
    }
  } else if (outcome_type != "norm") {
    true_effects_samples <- mu_samples
  } else if (random_mv) {
    true_effects_samples <- fixed_mu +
      .predict_brma_mv_random_posterior_draws(
        object            = object,
        mu_samples        = fixed_mu,
        posterior_samples = posterior_samples,
        bias_offset       = blup_bias_offset
      )
  } else if (context[["is_multilevel"]] &&
             !context[["is_weightfunction"]]) {
    true_effects_samples <- fixed_mu +
      .evaluate.brma.multilevel_posterior.norm(
        mu_samples  = fixed_mu,
        tau_within  = tau_within_samples,
        tau_between = tau_between_samples,
        yi          = outcome_data[["yi"]],
        vi          = blup_vi,
        cluster     = context[["fit_data"]][["cluster"]],
        bias_offset = blup_bias_offset
      )
  } else if (context[["is_known_v"]]) {
    true_effects_samples <- .evaluate.brma.known_v_posterior.norm(
      mu_samples  = fixed_mu,
      tau_within  = tau_within_samples,
      yi          = outcome_data[["yi"]],
      known_V     = .data_known_v_data(new_data),
      bias_offset = blup_bias_offset
    )
  } else {
    true_effects_samples <- .evaluate.brma.true_effects_posterior.norm(
      mu_samples  = if (context[["is_multilevel"]]) mu_samples else fixed_mu,
      tau_within  = tau_within_samples,
      yi          = outcome_data[["yi"]],
      sei         = sqrt(blup_vi),
      bias_offset = blup_bias_offset
    )
  }

  return(true_effects_samples)
}


.predict_brma_response <- function(context, location_state, scale_state) {

  object            <- context[["object"]]
  as_measure        <- context[["as_measure"]]
  bias_adjusted     <- context[["bias_adjusted"]]
  conditioning_depth <- context[["conditioning_depth"]]
  new_data          <- context[["new_data"]]
  known_V_new       <- context[["known_V_new"]]
  priors            <- context[["priors"]]
  is_weightfunction <- context[["is_weightfunction"]]
  is_known_v        <- context[["is_known_v"]]
  outcome_type      <- context[["outcome_type"]]
  posterior_samples <- context[["posterior_samples"]]
  outcome_data      <- context[["outcome_data"]]
  K                 <- context[["K"]]
  mu_samples        <- location_state[["mu"]]
  tau_within_samples <- scale_state[["within"]]

  # different model types have different output structures
  if (is.element(outcome_type, c("bin", "pois"))) {

    true_effects_samples <- .predict_brma_estimate_draws(
      context        = context,
      location_state = location_state,
      scale_state    = scale_state
    )

    if (as_measure) {
      .check_glmm_response_as_measure(outcome_type, outcome_data)
    }

    if (outcome_type == "bin") {

      if (conditioning_depth == "estimate") {
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
        mu_samples     = true_effects_samples,
        logit_baserate = logit_baserate,
        n1i            = outcome_data[["n1i"]],
        n2i            = outcome_data[["n2i"]]
      )

    } else if (outcome_type == "pois") {

      if (conditioning_depth == "estimate") {
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
        mu_samples = true_effects_samples,
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

    if (bias_adjusted || !is_weightfunction) {

      true_effects_samples <- .predict_brma_estimate_draws(
        context        = context,
        location_state = location_state,
        scale_state    = scale_state
      )

      zero_tau <- matrix(
        0,
        nrow = nrow(true_effects_samples),
        ncol = ncol(true_effects_samples)
      )
      if (is_known_v) {
        outcome_samples <- .outcome_rng.norm_known_v(
          mu_samples = true_effects_samples,
          tau_within = zero_tau,
          known_V    = if (!is.null(known_V_new)) known_V_new else .data_known_v_data(new_data)
        )
      } else {
        outcome_samples <- .outcome_rng.norm(
          mu_samples = true_effects_samples,
          tau_within = zero_tau,
          sei        = outcome_data[["sei"]]
        )
      }

    } else {

      response_mu  <- mu_samples
      response_tau <- tau_within_samples
      if (conditioning_depth == "marginal" && context[["is_multilevel"]]) {
        response_mu <- response_mu + .evaluate.brma.cluster_effects(
          fit               = object[["fit"]],
          tau_between       = scale_state[["between"]],
          cluster           = context[["fit_data"]][["cluster"]],
          same_data         = FALSE,
          effect_direction  = context[["effect_direction"]],
          posterior_samples = posterior_samples
        )
      }
      if (conditioning_depth == "estimate") {
        true_effects_samples <- .predict_brma_estimate_draws(
          context        = context,
          location_state = location_state,
          scale_state    = scale_state
        )
        response_mu  <- true_effects_samples
        response_tau <- matrix(
          0,
          nrow = nrow(true_effects_samples),
          ncol = ncol(true_effects_samples)
        )
      }

      selection_context <- .selection_context(
        object            = object,
        posterior_samples = posterior_samples,
        newdata           = new_data
      )

      outcome_samples <- .outcome_rng.selnorm(
        mu_samples        = response_mu,
        tau_within        = response_tau,
        sei               = outcome_data[["sei"]],
        selection_context = selection_context
      )

    }

    # rename samples
    colnames(outcome_samples) <- paste0("yi[", seq_len(K), "]")
  }

  return(.predict_brma_finalize(
    context    = context,
    samples    = outcome_samples,
    title      = .response_prediction_title(outcome_type, as_measure),
    parameters = .conditional_effect_parameters(object),
    effect     = !is.element(outcome_type, c("bin", "pois")) || as_measure
  ))

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
                                          parameters, posterior_samples = NULL,
                                          quiet = FALSE) {

  if (!conditional) {
    return(samples)
  }

  if (is.list(samples) && !is.matrix(samples)) {
    metadata <- attr(samples, "brma_mv_prediction_target", exact = TRUE)
    out <- lapply(samples, function(component_samples) {
      .condition_prediction_samples(
        object            = object,
        samples           = component_samples,
        conditional       = conditional,
        parameters        = parameters,
        posterior_samples = posterior_samples,
        quiet             = quiet
      )
    })
    if (!is.null(metadata)) {
      attr(out, "brma_mv_prediction_target") <- metadata
    }
    return(.new_brma_samples_list(out))
  }

  keep <- .conditional_parameter_rows(
    object            = object,
    parameters        = parameters,
    posterior_samples = posterior_samples,
    rule              = "OR"
  )
  sample_matrix <- as.matrix(samples)[keep, , drop = FALSE]

  if (nrow(sample_matrix) == 0L) {
    stop("No posterior samples remain after conditioning.", call. = FALSE)
  }

  if (!quiet) {
    message("Conditional posterior samples flattened to one chain.")
  }

  metadata <- attr(samples, "brma_mv_prediction_target", exact = TRUE)
  prediction_samples <- attr(samples, "prediction_samples", exact = TRUE)
  if (!is.null(prediction_samples)) {
    prediction_samples <- prediction_samples[keep, , drop = FALSE]
  }
  out <- .new_brma_samples(
    samples            = sample_matrix,
    n_chains           = 1L,
    n_iter             = nrow(sample_matrix),
    title              = paste("Conditional", attr(samples, "title")),
    probs              = attr(samples, "probs"),
    data               = attr(samples, "data"),
    effect_transform   = attr(samples, "effect_transform"),
    prediction_samples = prediction_samples
  )
  if (!is.null(metadata)) {
    attr(out, "brma_mv_prediction_target") <- metadata
  }

  return(out)
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
    indicator <- .extract_posterior_indicator(
      posterior_samples = posterior_samples,
      parameter         = parameter,
      prior             = prior
    )
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
    indicator <- .extract_posterior_indicator(
      posterior_samples = posterior_samples,
      parameter         = parameter,
      prior             = prior
    )
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

  return(.extract_posterior_indicator(
    posterior_samples = posterior_samples,
    parameter         = parameter
  ))
}
