#' @title Extract method for Robust Bayesian Meta-Analysis Fits
#'
#' @description \code{residuals.RoBMA} extract residuals based on the RoBMA model.
#' Only available for normal-normal models estimated using the spike-and-slab
#' algorithm (i.e., \code{algorithm = "ss"}).
#'
#' @inheritParams summary.RoBMA
#' @inheritParams pooled_effect
#' @inheritParams predict.RoBMA
#'
#' @examples \dontrun{
#' require(metafor)
#' dat <- escalc(measure = "OR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
#'               data = dat.bcg)
#'
#' # fit meta-regression
#' robma_dat <- data.frame(
#'   logOR  = dat$yi,
#'   se     = sqrt(dat$vi),
#'   ablat  = dat$ablat,
#'   alloc  = dat$alloc
#' )
#'
#' fit <- NoBMA.reg(~ ablat + alloc, data = robma_dat,
#'                  seed = 1, algorithm = "ss", parallel = TRUE)
#'
#' residuals(fit)
#'
#' }
#'
#' @return \code{pooled_effect} returns a list of tables of class 'BayesTools_table'.
#' @seealso [predict.RoBMA()]
#' @export
residuals.RoBMA <- function(object, conditional = FALSE, output_scale = NULL, probs = c(.025, .975), as_samples = FALSE, ...){

  # obtain predictions on the model scale
  # data are already on the model scale
  # scale the residuals back to the outcome scale
  # (residuals are differences, must be scaled instead of transformed)

  if(object[["add_info"]][["algorithm"]] != "ss")
    stop("Predictions can only be computed for spike and slab models.")
  if(inherits(object, "BiBMA") || inherits(object, "BiBMA.reg"))
    stop("The true effects can only be computed for normal-normal (NoBMA / RoBMA) models.")

  # get the model fitting scale
  if (is.BiBMA(object)) {
    model_scale <- "logOR"
  } else {
    model_scale <- object$add_info[["effect_measure"]]
  }
  if(is.null(output_scale)){
    output_scale <- object$add_info[["output_scale"]]
  }else if(object$add_info[["output_scale"]] == "y" & .transformation_var(output_scale) != "y"){
    stop("Models estimated using the general effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale)
  }

  # get the prediction (automatically checks all the input)
  preds <- predict.RoBMA(object, newdata = NULL, type = "terms",
                         conditional = conditional, output_scale = .transformation_invar(model_scale),
                         incorporate_publication_bias = TRUE, as_samples = TRUE)

  # get the data: dispatch between meta-regression / meta-analysis input
  data <- .get_outcome_data(object)

  # compute the residuals
  resids <- lapply(1:nrow(data), function(i){
    data[["y"]][i] - preds[[i]]
  })

  # transform estimates
  resids <- lapply(resids, function(x) {
    .scale(x, from = model_scale, to = output_scale)
  })
  names(resids) <- sapply(seq_along(resids), function(x) paste0("residual[", x, "]"))

  # return only samples if requested
  if(as_samples){
    return(resids)
  }

  # obtain estimates tables
  estimates <- BayesTools::ensemble_estimates_table(
    samples    = resids,
    parameters = names(resids),
    probs      = probs,
    title      = if(conditional) "Conditional residuals:" else "Residuals:",
    footnotes  = c(.scale_note_simple(object$add_info[["prior_scale"]], output_scale))
  )

  # create the output object
  output <- list(
    call       = object[["call"]],
    title      = .object_title(object),
    estimates  = estimates,
    footnotes  = c(.scale_note_simple(object$add_info[["prior_scale"]], output_scale))
  )


  class(output) <- "summary.RoBMA"
  attr(output, "type") <- "ensemble"

  return(output)
}


#' @title Extract Residuals from brma Object
#'
#' @description \code{residuals.brma} extracts residuals based on the brma model.
#' Residuals are computed as the difference between observed effect sizes and
#' predicted values from the model.
#'
#' @param object a fitted brma object
#' @param bias_adjusted whether to adjust predictions for publication bias.
#' Defaults to \code{TRUE}, returning bias-corrected residuals.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param as_samples whether posterior samples should be returned instead of
#' a summary table. Defaults to \code{FALSE}.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' Residuals are computed as \code{yi - predicted_mu}, where \code{predicted_mu}
#' is obtained from \code{predict.brma(type = "terms")}. The \code{type = "terms"}
#' prediction gives the fixed effects component (the mean effect size at the
#' predictor levels), which corresponds to the fitted values in standard
#' meta-analysis.
#'
#' For models with publication bias adjustment (selection models), the
#' \code{bias_adjusted} argument controls whether predictions account for
#' publication bias. When \code{bias_adjusted = TRUE} (default), residuals
#' reflect deviations from the bias-corrected effect size.
#'
#' @return If \code{as_samples = FALSE}, returns a \code{brma.predict} object
#' containing a summary table with residual estimates and credible intervals.
#' If \code{as_samples = TRUE}, returns a matrix of posterior samples for
#' residuals.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat, seed = 1)
#'
#' # extract residuals
#' residuals(fit)
#'
#' # get posterior samples
#' resid_samples <- residuals(fit, as_samples = TRUE)
#' }
#'
#' @seealso [predict.brma()]
#' @export
residuals.brma <- function(object,
                           bias_adjusted = TRUE,
                           probs = c(.025, .975),
                           as_samples = FALSE,
                           ...) {

  # input validation
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(as_samples, "as_samples")

  # check model type - only normal outcome models supported
  outcome_type <- .outcome_type(object)
  if (outcome_type != "norm") {
    stop("Residuals are only available for normal outcome models (not binomial or Poisson).", call. = FALSE)
  }

  # get predictions using predict.brma
  # type = "terms" gives fixed effects component (fitted values)
  preds <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "terms",
    bias_adjusted = bias_adjusted,
    as_samples    = TRUE,
    quiet         = TRUE
  )

  # extract observed data
  outcome_data <- object[["data"]][["outcome"]]
  K            <- nrow(outcome_data)

  # compute residuals: yi - predicted_mu
  # preds is S x K matrix of posterior samples
  residuals_matrix <- matrix(NA_real_, nrow = nrow(preds), ncol = K)
  for (k in seq_len(K)) {
    residuals_matrix[, k] <- outcome_data$yi[k] - preds[, k]
  }

  # set column names
  colnames(residuals_matrix) <- paste0("residual[", seq_len(K), "]")

  # return samples if requested
  if (as_samples) {
    return(residuals_matrix)
  }

  # create summary table
  residuals_table <- BayesTools::ensemble_estimates_table(
    samples    = asplit(residuals_matrix, 2),
    parameters = colnames(residuals_matrix),
    probs      = probs,
    title      = "Residuals:"
  )

  # create output object
  output <- list(
    summary = residuals_table,
    data    = object[["data"]]
  )
  class(output) <- "brma.predict"

  return(output)
}

