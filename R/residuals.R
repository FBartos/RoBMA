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

