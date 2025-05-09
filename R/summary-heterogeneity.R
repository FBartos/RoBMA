

#' @title Summarizes heterogeneity of a RoBMA model
#'
#' @description Computes the prediction interval, the absolute
#' heterogeneity  (tau, tau^2), and relative measures of heterogeneity
#' (I^2, H^2) for a fitted RoBMA object.
#'
#' @param type whether to show the overall RoBMA results (\code{"ensemble"})
#' or a detailed summary of the individual models (\code{"individual"}).
#' Can be abbreviated to first letters.
#' @inheritParams summary.RoBMA
#'
#' @details
#' The `conditional` argument allows for computing the conditional prediction interval based
#' on models assuming the presence of the effect and the conditional heterogeneity estimates
#' tau, tau^2, I^2, and H^2 assuming the presence of the heterogeneity.
#'
#' @return \code{summary_heterogeneity} returns a list of tables of class 'BayesTools_table'.
#'
#' @export
summary_heterogeneity <- function(object, type = "ensemble", conditional = FALSE,
                                  output_scale = NULL, probs = c(.025, .975),
                                  short_name = FALSE, remove_spike_0 = FALSE){

  # apply version changes to RoBMA object
  object <- .update_object(object)

  if(is.BiBMA(object))
    stop("The 'summary_heterogeneity' function is not available for Binomial meta-analytic models.")

  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_char(type, "type")
  BayesTools::check_bool(short_name, "short_name")
  BayesTools::check_bool(remove_spike_0, "remove_spike_0")

  # transform the data and posterior to scale used for fitting the model (to compute the metrics)
  model_scale <- object$add_info[["effect_measure"]]
  if(is.null(output_scale)){
    output_scale <- object$add_info[["output_scale"]]
  }else if(object$add_info[["output_scale"]] == "y" & .transformation_var(output_scale, estimation = FALSE) != "y"){
    stop("Models estimated using the generall effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale, estimation = FALSE)
  }

  if(model_scale != object$add_info[["output_scale"]]){
    object <- .transform_posterior(object, object$add_info[["output_scale"]], model_scale)
  }

  # v_tilde for I^2 and H^2 statistic
  if(.is_regression(object)){
    w <- 1/object[["data"]][["outcome"]][["se"]]^2
  }else{
    w <- 1/object[["data"]][["se"]]^2
  }
  v_tilde <- ((length(w) - 1) * sum(w)) / (sum(w)^2 - sum(w^2))


  if(substr(type,1,1) == "e" || object[["add_info"]][["algorithm"]] == "ss"){

    # prediction intervals (computed on the model scale -- returned on the prior scale)
    if(object[["add_info"]][["algorithm"]] == "ss"){
      PI <- list(
        "averaged"    = .compute_model_predictions(object[["model"]], model_scale),
        "conditional" = .compute_model_predictions(object[["model"]], model_scale, conditional = TRUE)
      )
    }else{
      PI <- .compute_ensemble_predictions(object, model_scale)
    }

    # obtain estimates tables
    estimates <- BayesTools::ensemble_estimates_table(
      samples    = list(
        PI    = .transform_mu(PI[["averaged"]], from = object$add_info[["prior_scale"]], to = output_scale),
        tau   = .scale(object$RoBMA$posteriors$tau,   model_scale, output_scale),
        tau2  = .scale(object$RoBMA$posteriors$tau^2, model_scale, output_scale),
        I2    = .compute_I2(object$RoBMA$posteriors$tau, v_tilde),
        H2    = .compute_H2(object$RoBMA$posteriors$tau, v_tilde)
      ),
      parameters = c("PI","tau", "tau2", "I2", "H2"),
      probs      = probs,
      title      = "Model-averaged heterogeneity estimates:",
      footnotes  = c(.heterogeneity_scale_note(model_scale, object$add_info[["prior_scale"]], output_scale, .is_regression(object))),
      warnings   = .collect_errors_and_warnings(object)
    )

    # deal with possibly empty table in case of no alternative models
    if(is.null(object$RoBMA[["posteriors_conditional"]])){
      estimates_conditional                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
      colnames(estimates_conditional)          <- c("Mean", "Median", probs)
      class(estimates_conditional)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_conditional))
      attr(estimates_conditional, "type")      <- rep("estimate", ncol(estimates_conditional))
      attr(estimates_conditional, "rownames")  <- TRUE
      attr(estimates_conditional, "title")     <- "Conditional estimates:"
      attr(estimates_conditional, "warnings")  <- .collect_errors_and_warnings(object)
    }else{
      estimates_conditional <- BayesTools::ensemble_estimates_table(
        samples    = list(
          PI    = .transform_mu(PI[["conditional"]], from = object$add_info[["prior_scale"]], to = output_scale),
          tau   = .scale(object$RoBMA$posteriors_conditional$tau,   model_scale, output_scale),
          tau2  = .scale(object$RoBMA$posteriors_conditional$tau^2, model_scale, output_scale),
          I2    = .compute_I2(object$RoBMA$posteriors_conditional$tau, v_tilde),
          H2    = .compute_H2(object$RoBMA$posteriors_conditional$tau, v_tilde)
        ),
        parameters = c("PI", "tau", "tau2", "I2", "H2"),
        probs      = probs,
        title      = "Conditional heterogeneity estimates:",
        footnotes  = c(.heterogeneity_scale_note(model_scale, object$add_info[["prior_scale"]], output_scale, .is_regression(object))),
        warnings   = .collect_errors_and_warnings(object)
      )
    }

    # create the output object
    output <- list(
      call       = object[["call"]],
      title      = .object_title(object),
      estimates  = estimates
    )

    if(conditional){
      output$estimates_conditional <- estimates_conditional
    }

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "ensemble"

    return(output)

  }else if(substr(type,1,1) == "i"){

    output <- list(
      call       = object[["call"]],
      title      = .object_title(object),
      models     = list()
    )

    for(i in seq_along(object[["models"]])){

      summary  <- BayesTools::model_summary_table(
        model             = object[["models"]][[i]],
        remove_parameters = "gamma",
        short_name        = short_name,
        remove_spike_0    = remove_spike_0
      )

      if(BayesTools::is.prior.point(object[["models"]][[i]]$priors[["tau"]])){
        tau <- object[["models"]][[i]]$priors[["tau"]]$parameters[["location"]]
      }else{
        tau <- suppressWarnings(coda::as.mcmc(object[["models"]][[i]][["fit"]]))[,"tau"]
      }

      PI <- .compute_model_predictions(object[["models"]][[i]], model_scale, object$add_info[["seed"]])

      estimates <- BayesTools::ensemble_estimates_table(
        samples    = list(
          PI    = .transform_mu(PI, from = object[["models"]][[i]][["prior_scale"]], to = output_scale),
          tau   = .scale(tau, object[["models"]][[i]][["prior_scale"]], output_scale),
          tau2  = .scale(.scale(tau, object[["models"]][[i]][["prior_scale"]], model_scale)^2, model_scale, output_scale),
          I2    = .compute_I2(.scale(tau, object[["models"]][[i]][["prior_scale"]], model_scale), v_tilde),
          H2    = .compute_H2(.scale(tau, object[["models"]][[i]][["prior_scale"]], model_scale), v_tilde)
        ),
        parameters = c("PI", "tau", "tau2", "I2", "H2"),
        probs      = probs,
        title      = "Heterogeneity estimates:",
        footnotes  = c(.heterogeneity_scale_note(model_scale, object$add_info[["prior_scale"]], output_scale, .is_regression(object))),
        warnings   = .collect_errors_and_warnings(object)
      )

      output[["models"]][[i]] <- list(
        summary   = summary,
        estimates = estimates
      )

    }

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "individual"

    return(output)

  }else{
    stop(paste0("Unknown summary type: '", type, "'."))
  }
}


# heterogeneity calculation functions according to metafor description in metafor::print.rma function
.compute_v_tilde <- function(se){
  w       <- 1/se^2
  v_tilde <- ((length(w) - 1) * sum(w)) / (sum(w)^2 - sum(w^2))
  return(v_tilde)
}
.compute_I2      <- function(tau, v_tilde){
  100 * tau^2 / (tau^2 + v_tilde)
}
.compute_H2      <- function(tau, v_tilde){
  (tau^2 + v_tilde) / v_tilde
}

.compute_model_predictions    <- function(model, model_scale, seed = NULL, conditional = FALSE){

  if(model[["prior_scale"]] != model[["output_scale"]])
    stop("The prior_scale does not match the output_scale. (Individual models' MCMC samples must have been transformed earlier.)")

  posterior_samples <- try(suppressWarnings(coda::as.mcmc(model[["fit"]])))

  if(.is_model_regression(model)){
    if(BayesTools::is.prior.point(model$priors$terms[["intercept"]])){
      mu <- model$priors$terms[["intercept"]]$parameters[["location"]]
    }else{
      mu <- posterior_samples[,"mu_intercept"]
    }
  }else{
    if(BayesTools::is.prior.point(model$priors[["mu"]])){
      mu <- model$priors[["mu"]]$parameters[["location"]]
    }else{
      mu <- posterior_samples[,"mu"]
    }
  }

  if(BayesTools::is.prior.point(model$priors[["tau"]])){
    tau <- model$priors[["tau"]]$parameters[["location"]]
  }else{
    tau <- posterior_samples[,"tau"]
  }

  mu  <- .scale(mu,  model[["prior_scale"]], model_scale)
  tau <- .scale(tau, model[["prior_scale"]], model_scale)

  if(is.null(seed)){
    set.seed(1)
  }else{
    set.seed(seed)
  }

  # these needs to be simulated as the posteriors can be non-normal
  predictions <- stats::rnorm(n = max(length(mu), length(tau))*10, mean = mu, sd = tau)
  predictions <- .scale(predictions, model_scale, model[["prior_scale"]])

  if(conditional){

    if(.is_model_regression(model)) {
      mu_is_null   <- attr(model$priors$terms[["intercept"]], "components") == "null"
      mu_indicator <- posterior_samples[,"mu_intercept_indicator"]
    } else {
      mu_is_null   <- attr(model$priors[["mu"]], "components") == "null"
      mu_indicator <- posterior_samples[,"mu_indicator"]
    }

    predictions <- predictions[mu_indicator %in% which(!mu_is_null)]
  }

  return(predictions)
}
.compute_ensemble_predictions <- function(object, model_scale, n_samples = 100000){

  model_predictions      <- lapply(object[["models"]], function(model) .compute_model_predictions(model, model_scale, object$add_info[["seed"]]))
  post_probs             <- sapply(object[["models"]], function(model) model$inference[["post_prob"]])
  post_probs_conditional <- object$RoBMA$inference_conditional$Effect[["post_probs"]]

  if(is.null(object$add_info[["seed"]])){
    set.seed(1)
  }else{
    set.seed(object$add_info[["seed"]])
  }

  prediction_interval             <- .mix_predictions(model_predictions, post_probs, n_samples)
  prediction_interval_conditional <- .mix_predictions(model_predictions, post_probs_conditional, n_samples)

  return(list(
    averaged    = prediction_interval,
    conditional = prediction_interval_conditional
  ))
}
.mix_predictions              <- function(model_predictions, post_probs, n_samples){

  predictions <- c()
  sample_ind  <- NULL
  models_ind  <- NULL
  for (i in seq_along(model_predictions)[round(post_probs * n_samples) > 1]) {

    temp_ind <- sample(length(model_predictions[[i]]), round(n_samples * post_probs[i]), replace = TRUE)
    predictions <- c(predictions, model_predictions[[i]][temp_ind])

    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  attr(predictions, "sample_ind") <- sample_ind
  attr(predictions, "models_ind") <- models_ind

  return(predictions)
}
.heterogeneity_scale_note <- function(model_scale, prior_scale, output_scale, regression){
  return(sprintf(
    "The prediction interval %1$s(PI) is summarized on the %2$s scale.\nThe absolute heterogeneity (tau, tau^2) is summarized on the %3$s scale.\nThe relative heterogeneity indicies (I^2 and H^2) were computed on the %4$s scale.",
    if(regression) "for the average effect " else "",
    .transformation_names(output_scale, estimation = FALSE),
    switch(
      output_scale,
      "OR" = .transformation_names("logOR", estimation = FALSE),
      "r"  = .transformation_names("z", estimation = FALSE),
      .transformation_names(output_scale, estimation = FALSE)
    ),
    .transformation_names(model_scale, estimation = FALSE)))
}
