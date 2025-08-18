#' @title Transform RoBMA object into z-curve
#'
#' @description
#' This function transforms the estimated RoBMA model into a
#' z-curve object that can be further summarized and plotted.
#' Only available for normal-normal models estimated using the spike-and-slab
#' algorithm (i.e., \code{algorithm = "ss"}).
#'
#' @param x A RoBMA object
#' @param significance_level Significance level used for computation of z-curve
#' estimates.
#' @param max_samples Maximum number of samples from the posterior distribution
#' that will be used for estimating z-curve estimates.
#'
#' @return \code{as_zcurve} returns a list of tables of class 'zcurve_RoBMA'.
#'
#' @export
as_zcurve <- function(x, significance_level = stats::qnorm(0.975), max_samples = 1000){

  if(x[["add_info"]][["algorithm"]] != "ss")
    stop("Predictions can only be computed for spike and slab models.")
  if(inherits(x, "BiBMA") || inherits(x, "BiBMA.reg"))
    stop("The true effects can only be computed for normal-normal (NoBMA / RoBMA) models.")
  BayesTools::check_real(significance_level, "significance_level", lower = 0)
  BayesTools::check_int(max_samples, "max_samples", lower = 10)

  # compute the main estimates and CIs
  # plotting densities are generated only if figures are requested
  EDR             <- .zcurve_EDR(x, z_threshold = significance_level, max_samples = max_samples, conditional = FALSE)
  EDR_conditional <- try(.zcurve_EDR(x, z_threshold = significance_level, max_samples = max_samples, conditional = TRUE))


  # compute data summaries
  if(inherits(x, "RoBMA.reg") || inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg")){
    outcome_data <- x$data[["outcome"]]
  }else{
    outcome_data <- x[["data"]]
  }

  # store new z-curve objects
  x[["coefficients"]] <- c("EDR" = mean(EDR))
  x[["zcurve"]]       <- list(
    estimates = list(
      EDR = EDR
    ),
    estimates_conditional = list(
      EDR = EDR_conditional
    ),
    data = list(
      significance_level = significance_level,
      z                  = outcome_data[["y"]] / outcome_data[["se"]],
      N_significant      = sum(abs(outcome_data[["y"]] / outcome_data[["se"]]) > significance_level),
      N_observed         = nrow(outcome_data)
    )
  )

  # add class
  class(x) <- c("zcurve_RoBMA", class(x))

  return(x)
}


#' @title Summarize fitted zcurve_RoBMA object
#'
#' @description \code{summary.zcurve_RoBMA} creates summary tables for a
#' zcurve_RoBMA object.
#'
#' @inheritParams summary.RoBMA
#'
#' @return \code{summary.RoBMA} returns a list of tables of class 'BayesTools_table'.
#'
#' @seealso [as_zcurve()]
#' @export
summary.zcurve_RoBMA <- function(object, conditional = FALSE, probs = c(.025, .975), ...){

  # compute summary for the z-curve object
  # most functions are based on the zcurve package

  # get proportion of significant results
  obs_proportion <- stats::prop.test(object$zcurve$data[["N_significant"]], object$zcurve$data[["N_observed"]], conf.level = 0.95)
  info_text      <- sprintf("Fitted using %1$i observation, %2$i significant (ODR = %3$.2f, 95 CI [%4$.2f, %5$.2f]).",
                            object$zcurve$data[["N_observed"]],
                            object$zcurve$data[["N_significant"]],
                            obs_proportion$estimate,
                            obs_proportion$conf.int[1],
                            obs_proportion$conf.int[2]
  )

  estimates <- cbind.data.frame(
    "EDR"           = object$zcurve$estimates[["EDR"]],
    "Soric FDR"     = .get_Soric_FDR(object$zcurve$estimates[["EDR"]], stats::pnorm(object$zcurve$data[["significance_level"]], lower.tail = FALSE) * 2),
    "File Drawer R" = .get_file_drawer_R(object$zcurve$estimates[["EDR"]]),
    "Expected N"    = .get_expected_N(object$zcurve$estimates[["EDR"]], object$zcurve$data[["N_significant"]]),
    "Missing N"     = .get_missing_N(object$zcurve$estimates[["EDR"]], object$zcurve$data[["N_significant"]], object$zcurve$data[["N_observed"]])
  )

  estimates <- BayesTools::ensemble_estimates_table(
    samples    = estimates,
    parameters = names(estimates),
    probs      = probs,
    title      = "Model-averaged estimates:",
    footnotes  = info_text,
    warnings   = .collect_errors_and_warnings(object)
  )

  if(conditional){

    estimates_conditional <- cbind.data.frame(
      "EDR"           = object$zcurve$estimates_conditional[["EDR"]],
      "Soric FDR"     = .get_Soric_FDR(object$zcurve$estimates_conditional[["EDR"]], stats::pnorm(object$zcurve$data[["significance_level"]], lower.tail = FALSE) * 2),
      "File Drawer R" = .get_file_drawer_R(object$zcurve$estimates_conditional[["EDR"]]),
      "Expected N"    = .get_expected_N(object$zcurve$estimates_conditional[["EDR"]], object$zcurve$data[["N_significant"]]),
      "Missing N"     = .get_missing_N(object$zcurve$estimates_conditional[["EDR"]], object$zcurve$data[["N_significant"]], object$zcurve$data[["N_observed"]])
    )

    estimates_conditional <- BayesTools::ensemble_estimates_table(
      samples    = estimates_conditional,
      parameters = names(estimates_conditional),
      probs      = probs,
      title      = "Conditional estimates:",
      footnotes  = info_text,
      warnings   = .collect_errors_and_warnings(object)
    )

  }

  # create the output object
  output <- list(
    call       = object[["call"]],
    title      = "Z-curve",
    estimates  = estimates
  )

  if(conditional){
    output$estimates_conditional <- estimates_conditional
  }

  class(output) <- "summary.zcurve_RoBMA"
  return(output)
}


#' @title Plots fitted zcurve_RoBMA object
#'
#' @description \code{plot.zcurve_RoBMA} creates a z-curve figure for a
#' zcurve_RoBMA object.
#'
#' @inheritParams as_zcurve
#' @inheritParams plot.RoBMA
#' @inheritParams summary.RoBMA
#' @param z_sequence Sequence determining the plotting range and binning of z-values. Defaults to \code{seq(-6, 6, 0.2)}.
#' @param plot_fit Whether the model fit should be included in the figure. Defaults to \code{TRUE}.
#' @param plot_extrapolation Whether the model extrapolation should be included in the figure. Defaults to \code{TRUE}.
#' @param plot_CI Whether the credible intervals should be included in the figure. Defaults to \code{TRUE}.
#'
#' @return \code{plot.zcurve_RoBMA} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [as_zcurve()]
#' @export
plot.zcurve_RoBMA <- function(x, conditional = FALSE, plot_type = "base",
                              probs = c(.025, .975), max_samples = 500,
                              z_sequence = seq(-6, 6, 0.2),
                              plot_fit = TRUE, plot_extrapolation = TRUE, plot_CI = TRUE, ...){

  # plots the z-curve object
  # most functions are based on the zcurve package
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_int(max_samples, "max_samples", lower = 10)
  BayesTools::check_bool(plot_fit, "plot_fit")
  BayesTools::check_bool(plot_extrapolation, "plot_extrapolation")
  BayesTools::check_bool(plot_CI, "plot_CI")

  # extract the data
  dots <- list()
  z    <- x$zcurve$data[["z"]]
  z_in_range <- z > min(z_sequence) & z < max(z_sequence)

  # create z-curve histogram data
  z_hist <- graphics::hist(z[z_in_range], breaks = z_sequence, plot = FALSE)

  # adjust histogram for values outside of plotting range
  z_hist$density <- z_hist$density * mean(z_in_range)

  # compute the expected density under the RoBMA model
  z_densities_fit           <- .zcurve_densities(x, z_sequence = z_sequence, max_samples = max_samples, conditional = conditional, extrapolate = FALSE)
  z_densities_extrapolation <- .zcurve_densities(x, z_sequence = z_sequence, max_samples = max_samples, conditional = conditional, extrapolate = TRUE)

  # create plotting objects
  df_hist <- data.frame(
    x       = z_hist$mids,
    density = z_hist$density,
    breaks  = diff(z_hist$breaks)
  )
  df_fit <- data.frame(
    x     = z_sequence,
    y     = colMeans(z_densities_fit),
    y_lCI = apply(z_densities_fit, 2, stats::quantile, probs = 0.025),
    y_uCI = apply(z_densities_fit, 2, stats::quantile, probs = 0.975)
  )
  df_extrapolation <- data.frame(
    x     = z_sequence,
    y     = colMeans(z_densities_extrapolation),
    y_lCI = apply(z_densities_extrapolation, 2, stats::quantile, probs = 0.025),
    y_uCI = apply(z_densities_extrapolation, 2, stats::quantile, probs = 0.975)
  )

  # allow data return for JASP
  if(isTRUE(dots[["as_data"]])){
    return(list(
      hist          = df_hist,
      fit           = df_fit,
      extrapolation = df_extrapolation
    ))
  }


  if(plot_type == "ggplot"){

    # Initialize ggplot with histogram as density
    out <- ggplot2::ggplot() +
      ggplot2::geom_col(
        ggplot2::aes(
          x = df_hist$x,
          y = df_hist$density),
        fill  = NA,
        color = "black",
        width = df_hist$breaks
      ) +
      ggplot2::labs(x = "z", y = "Density")

    # Add extrapolation fit
    if(plot_extrapolation){
      if(plot_CI){
        out <- out + ggplot2::geom_ribbon(
          ggplot2::aes(
            x    = df_extrapolation$x,
            ymin = df_extrapolation$y_lCI,
            ymax = df_extrapolation$y_uCI),
          fill = scales::alpha("blue", 0.2)
        )
      }
      out <- out + ggplot2::geom_line(
        ggplot2::aes(
          x = df_extrapolation$x,
          y = df_extrapolation$y),
        color     = "blue",
        linewidth = 1
      )
    }

    # Add extrapolated region
    if(plot_fit){
      if(plot_CI){
        out <- out + ggplot2::geom_ribbon(
          ggplot2::aes(
            x    = df_fit$x,
            ymin = df_fit$y_lCI,
            ymax = df_fit$y_uCI),
          fill = scales::alpha("gray40", 0.4)
        )
      }
      out <- out + ggplot2::geom_line(
        ggplot2::aes(
          x = df_fit$x,
          y = df_fit$y),
        color     = "gray40",
        linewidth = 1,
        linetype  = 2
      )
    }

  }else if(plot_type == "base"){

    if(!is.null(dots[["ylim"]])){
      ylim <- dots[["ylim"]]
    }else{
      ylim <- c(0, max(z_hist$density))
      if(plot_fit){
        if(plot_CI){
          ylim <- range(c(ylim, df_extrapolation$y_uCI))
        }else{
          ylim <- range(c(ylim, df_extrapolation$y))
        }
      }
      if(plot_extrapolation){
        if(plot_CI){
          ylim <- range(c(ylim, df_fit$y_uCI))
        }else{
          ylim <- range(c(ylim, df_fit$y))
        }
      }
    }


    graphics::plot(z_hist, freq = FALSE, las = 1, density = 0, angle = 0,
                   border = "black", xlab = "z", main = "", ylim = ylim)

    if(plot_extrapolation){
      if(plot_CI){
        graphics::polygon(
          c(df_extrapolation$x, rev(df_extrapolation$x)),
          c(df_extrapolation$y_lCI, rev(df_extrapolation$y_uCI)),
          border = NA, col = scales::alpha("blue", .20))
      }
      graphics::lines(df_extrapolation$x, df_extrapolation$y, lwd = 2, col = "blue")
    }

    if(plot_fit){
      if(plot_CI){
        graphics::polygon(
          c(df_fit$x, rev(df_fit$x)),
          c(df_fit$y_lCI, rev(df_fit$y_uCI)),
          border = NA, col = scales::alpha("gray40", .40))
      }
      graphics::lines(df_fit$x, df_fit$y, lwd = 2, col = "gray40", lty = 2)
    }

  }


  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(out)
  }
}

.plot_zcurve_range <- function(z, z_min, z_max, z_step){

  z_seq <- seq(from = z_min, to = z_max, by = z_step)

  # include all thresholds larger then min(z)
  if(min(z) > z_min){
    is_larger <- min(which(z_seq > min(z)))
    if(is_larger > 1){
      z_seq <- z_seq[(is_larger-1):length(z_seq)]
    }
  }

  # include all thresholds smaller then max(z)
  if(max(z) < z_max){
    is_lower <- max(which(z_seq < max(z)))
    if(is_lower < length(z_seq)){
      z_seq <- z_seq[1:(is_lower+1)]
    }
  }

  return(z_seq)
}


#' @title Prints summary object for zcurve_RoBMA method
#'
#' @param x a summary of a zcurve_RoBMA object
#' @param ... additional arguments
#'
#'
#' @return \code{summary.zcurve_RoBMA} invisibly returns the print statement.
#'
#' @seealso [as_zcurve()]
#' @export
print.summary.zcurve_RoBMA <- function(x, ...){

  cat("Call:\nas_zcurve: ")
  print(x[["call"]])

  cat("\n")
  cat(x[["title"]])

  for(type in c("estimates", "estimates_conditional")){
    if(!is.null(x[[type]])){
      cat("\n")
      print(x[[type]])
    }
  }
}


# imported from zcurve
.get_Soric_FDR     <- function(EDR, sig_level){
  ((1/EDR) - 1)*(sig_level/(1-sig_level))
}
.get_file_drawer_R <- function(EDR){
  (1-EDR)/EDR
}
.get_expected_N    <- function(EDR, N_sig){
  .get_file_drawer_R(EDR)*N_sig + N_sig
}
.get_missing_N     <- function(EDR, N_sig, N_obs){
  .get_expected_N(EDR, N_sig) - N_obs
}


#' @title Prints a fitted zcurve_RoBMA object
#'
#' @param x a fitted zcurve_RoBMA object.
#' @param ... additional arguments.
#'
#'
#' @return \code{print.zcurve_RoBMA} invisibly returns the print statement.
#'
#' @seealso [as_zcurve()]
#' @export
print.zcurve_RoBMA <- function(x, ...){
  cat("Call:\nas_zcurve: ")
  print(x$call)
  cat("\nEstimates:\n")
  print(stats::coef(x))
}


.zcurve_EDR       <- function(x, z_threshold, max_samples = 1000, conditional = FALSE){

  # get the model fitting scale
  if (is.BiBMA(x)) {
    model_scale <- "logOR"
  } else {
    model_scale <- x$add_info[["effect_measure"]]
  }

  # extract posterior samples (and obtain conditional indicator)
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["model"]][["fit"]]))
  priors            <- x[["priors"]]

  # subset the samples to speed up the computation
  if(conditional){
    # if conditional output is to be provided, condition first and then subset
    # otherwise there might be no conditional samples left

    # select the indicator
    if(inherits(x, "RoBMA.reg") || inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg")){
      mu_indicator <- posterior_samples[,"mu_intercept_indicator"]
      mu_is_null   <- attr(x[["model"]]$priors$terms[["intercept"]], "components") == "null"
    }else{
      mu_indicator <- posterior_samples[,"mu_indicator"]
      mu_is_null   <- attr(x[["model"]]$priors$mu, "components") == "null"
    }
    mu_indicator <- mu_indicator %in% which(!mu_is_null)

    selected_samples_ind <- unique(round(seq(from = 1, to = sum(mu_indicator), length.out = max_samples)))
    selected_samples_ind <- seq_len(nrow(posterior_samples))[mu_indicator][selected_samples_ind]
    n_samples            <- length(selected_samples_ind)
    posterior_samples    <- posterior_samples[selected_samples_ind,,drop=FALSE]

  }else{

    selected_samples_ind <- unique(round(seq(from = 1, to = nrow(posterior_samples), length.out = max_samples)))
    n_samples            <- length(selected_samples_ind)
    posterior_samples    <- posterior_samples[selected_samples_ind,,drop=FALSE]
  }

  # dispatch between meta-regression / meta-analysis input
  if(inherits(x, "RoBMA.reg") || inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg")){
    newdata.predictors <- do.call(cbind.data.frame, x$data[["predictors"]])
    newdata.outcome    <- x$data[["outcome"]]
  }else{
    newdata.outcome <- x[["data"]]
  }

  # obtain the (study-specific) mu estimate
  # meta-regression and meta-analysis separately
  if(inherits(x, "RoBMA.reg") || inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg")){

    mu_samples  <- t(BayesTools::JAGS_evaluate_formula(
      fit         = x$model$fit,
      formula     = x$formula,
      parameter   = "mu",
      data        = newdata.predictors,
      prior_list  = attr(x$model$fit, "prior_list")
    ))[selected_samples_ind,,drop=FALSE]

  }else{

    mu_samples   <- matrix(posterior_samples[,"mu"], ncol = nrow(newdata.outcome), nrow = n_samples)

  }

  # transform the samples to the model fitting scale (the data are at the model fitting scale)
  mu_samples  <- .scale(mu_samples,  x$add_info[["output_scale"]], model_scale)

  # adjust effect samples and data to match original observed data direction
  # (to undo fitting adjustment for one-sided weightfunction)
  if(!is.null(x$add_info[["effect_direction"]]) && x$add_info[["effect_direction"]] == "negative"){
    mu_samples            <- -mu_samples
    newdata.outcome[,"y"] <- -newdata.outcome[,"y"] # (this one is never used but we flipp it as well for consistency)
  }

  # add conditional warnings
  if(conditional){
    if(sum(mu_indicator) < 100)
      warning(gettextf("There is only a very small number of posterior samples (%1s) assuming presence of the effect. The resulting estimates are not reliable.",
                       sum(mu_indicator)), call. = FALSE, immediate. = TRUE)
    if(sum(mu_indicator) <= 2)
      stop("Less or equal to 2 posterior samples assuming presence of the effects. The estimates could not be computed.")
  }

  # predicting responses without selection models does not require incorporating the between-study random-effects
  # (the marginalized and non-marginalized parameterization are equivalent)
  tau_samples <- posterior_samples[,"tau"]
  tau_samples <- .scale(tau_samples, x$add_info[["output_scale"]], model_scale)

  # compute the proportion of estimates larger than threshold under the population parameters
  outcome_thresholds <- rep(NA, nrow(mu_samples))

  for(j in 1:nrow(mu_samples)){
    # create containers for temporal samples from the posterior distribution
    temp_thresholds <- rep(NA, ncol(mu_samples))

    # compute the densities and threshold for each observation
    for(i in seq_len(ncol(mu_samples))){
      temp_thresholds[i] <-
        stats::pnorm(z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE) +
        stats::pnorm(-z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
    }

    # store the results
    outcome_thresholds[j] <- mean(temp_thresholds)
  }

  return(outcome_thresholds)
}
.zcurve_densities <- function(x, z_sequence, max_samples = 1000, conditional = FALSE, extrapolate = FALSE){

  # two modes of the function
  # extrapolate = FALSE: produces z-density of the fitted model as is == fit assessment
  #                      - requires incoporating publication bias indicies (i.e., PET, PEESE, weighted likelihoods)
  # extrapolate = TRUE: produces z-density of the expected unbiased distribution
  #                     - assummes that selection and SE dependent bias would not come into play
  #                     - whereever weighted likelihoods are involved, they need to be inverse-weighted to extrapolate to missing estimates

  # get the model fitting scale
  if (is.BiBMA(x)) {
    model_scale <- "logOR"
  } else {
    model_scale <- x$add_info[["effect_measure"]]
  }

  # extract posterior samples (and obtain conditional indicator)
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["model"]][["fit"]]))
  priors            <- x[["priors"]]

  # subset the samples to speed up the computation
  if(conditional){
    # if conditional output is to be provided, condition first and then subset
    # otherwise there might be no conditional samples left

    # select the indicator
    if(inherits(x, "RoBMA.reg") || inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg")){
      mu_indicator <- posterior_samples[,"mu_intercept_indicator"]
      mu_is_null   <- attr(x[["model"]]$priors$terms[["intercept"]], "components") == "null"
    }else{
      mu_indicator <- posterior_samples[,"mu_indicator"]
      mu_is_null   <- attr(x[["model"]]$priors$mu, "components") == "null"
    }
    mu_indicator <- mu_indicator %in% which(!mu_is_null)

    selected_samples_ind <- unique(round(seq(from = 1, to = sum(mu_indicator), length.out = max_samples)))
    selected_samples_ind <- seq_len(nrow(posterior_samples))[mu_indicator][selected_samples_ind]
    n_samples            <- length(selected_samples_ind)
    posterior_samples    <- posterior_samples[selected_samples_ind,,drop=FALSE]

  }else{

    selected_samples_ind <- unique(round(seq(from = 1, to = nrow(posterior_samples), length.out = max_samples)))
    n_samples            <- length(selected_samples_ind)
    posterior_samples    <- posterior_samples[selected_samples_ind,,drop=FALSE]
  }

  # dispatch between meta-regression / meta-analysis input
  if(inherits(x, "RoBMA.reg") || inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg")){
    newdata.predictors <- do.call(cbind.data.frame, x$data[["predictors"]])
    newdata.outcome    <- x$data[["outcome"]]
  }else{
    newdata.outcome <- x[["data"]]
  }

  # obtain the (study-specific) mu estimate
  # meta-regression and meta-analysis separately
  if(inherits(x, "RoBMA.reg") || inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg")){

    mu_samples  <- t(BayesTools::JAGS_evaluate_formula(
      fit         = x$model$fit,
      formula     = x$formula,
      parameter   = "mu",
      data        = newdata.predictors,
      prior_list  = attr(x$model$fit, "prior_list")
    ))[selected_samples_ind,,drop=FALSE]

  }else{

    mu_samples   <- matrix(posterior_samples[,"mu"], ncol = nrow(newdata.outcome), nrow = n_samples)

  }

  # transform the samples to the model fitting scale (the data are at the model fitting scale)
  mu_samples  <- .scale(mu_samples,  x$add_info[["output_scale"]], model_scale)

  # add conditional warnings
  if(conditional){
    if(sum(mu_indicator) < 100)
      warning(gettextf("There is only a very small number of posterior samples (%1s) assuming presence of the effect. The resulting estimates are not reliable.",
                       sum(mu_indicator)), call. = FALSE, immediate. = TRUE)
    if(sum(mu_indicator) <= 2)
      stop("Less or equal to 2 posterior samples assuming presence of the effects. The estimates could not be computed.")
  }

  # adjust effect samples and data to match original observed data direction
  # (to undo fitting adjustment for one-sided weightfunction)
  if(!is.null(x$add_info[["effect_direction"]]) && x$add_info[["effect_direction"]] == "negative"){
    newdata.outcome[,"y"] <- -newdata.outcome[,"y"]
    z_sequence            <- -z_sequence
    mu_samples            <- -mu_samples
  }

  # extract the list of priors
  priors_bias  <- priors[["bias"]]
  if(!extrapolate){
    # add PET/PEESE adjustment
    if(any(sapply(priors_bias, is.prior.PET))){
      PET_samples <- posterior_samples[,"PET"]
      # PET is scale invariant (no-scaling needed)
    }else{
      PET_samples <- rep(0, n_samples)
    }
    if(any(sapply(priors_bias, is.prior.PEESE))){
      PEESE_samples <- posterior_samples[,"PEESE"]
      # PEESE scales with the inverse
      PEESE_samples <- .scale(PEESE_samples, model_scale, x$add_info[["output_scale"]])
    }else{
      PEESE_samples <- rep(0, n_samples)
    }

    for(i in seq_len(ncol(mu_samples))){
      mu_samples[,i] <- mu_samples[,i] + PET_samples * newdata.outcome[i,"se"] + PEESE_samples * newdata.outcome[i,"se"]^2
    }
  }

  if(inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg") || (length(priors[["bias"]]) == 1 && is.prior.none(priors[["bias"]][[1]]))){

    # predicting responses without selection models does not require incorporating the between-study random-effects
    # (the marginalized and non-marginalized parameterization are equivalent)
    tau_samples <- posterior_samples[,"tau"]
    tau_samples <- .scale(tau_samples, x$add_info[["output_scale"]], model_scale)

    # compute the density ad specified support
    # outcome_lower     <- rep(NA, nrow(mu_samples))
    # outcome_higher    <- rep(NA, nrow(mu_samples))
    outcome_densities <- matrix(NA, nrow = nrow(mu_samples), ncol = length(z_sequence))

    for(j in 1:nrow(mu_samples)){
      # create containers for temporal samples from the posterior distribution
      # temp_lower     <- rep(NA, ncol(mu_samples))
      # temp_higher    <- rep(NA, ncol(mu_samples))
      temp_densities <- matrix(NA, nrow = ncol(mu_samples), ncol = length(z_sequence))

      # compute the densities and threshold for each observation
      for(i in seq_len(ncol(mu_samples))){
        # temp_lower[i]      <- stats::pnorm(z_sequence[1]                  * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
        # temp_higher[i]     <- stats::pnorm(z_sequence[length(z_sequence)] * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE)
        temp_densities[i,] <- stats::dnorm(z_sequence * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2)) * newdata.outcome[i,"se"]
        # the density needs to be transformed due to the support change
      }

      # store the results
      # outcome_lower[j]      <- mean(temp_lower)
      # outcome_higher[j]     <- mean(temp_higher)
      outcome_densities[j,] <- colMeans(temp_densities)
    }

  }else{

    # required for study ids / crit_x values in selection models
    fit_data <- .fit_data_ss(
      data             = newdata.outcome,
      priors           = priors,
      effect_direction = x$add_info[["effect_direction"]],
      prior_scale      = x$add_info[["prior_scale"]],
      weighted         = FALSE,
      weighted_type    = FALSE,
      multivariate     = .is_multivariate(x)
    )

    # predicting response requires incorporating the between-study random effects if selection models are present
    # (we use approximate selection likelihood which samples the true study effects instead of marginalizing them)
    if(.is_multivariate(x)){

      tau_samples <- posterior_samples[,"tau"]
      rho_samples <- posterior_samples[,"rho"]
      # deal with computer precision errors from JAGS
      rho_samples[rho_samples>1] <- 1
      rho_samples[rho_samples<0] <- 0
      # tau_within  = tau * sqrt(rho)
      # tau_between = tau * sqrt(1-rho)
      tau_within_samples  <- tau_samples * sqrt(rho_samples)
      tau_between_samples <- tau_samples * sqrt(1-rho_samples)
      gamma_samples       <- posterior_samples[,grep("gamma", colnames(posterior_samples)),drop = FALSE]

      tau_between_samples <- .scale(tau_between_samples, x$add_info[["output_scale"]], model_scale)
      tau_within_samples  <- .scale(tau_within_samples,  x$add_info[["output_scale"]], model_scale)

      # incorporate within study heterogeneity into the predictor
      # either estimated for prediction on the same data or integrated over for new data
      for(i in seq_len(nrow(newdata.outcome))){
        mu_samples[,i] <- mu_samples[,i] + gamma_samples[,fit_data$study_ids[i]] * tau_within_samples
      }

      # tau_between samples work as tau for the final sampling step
      tau_samples <- tau_between_samples

    }else{

      tau_samples  <- posterior_samples[,"tau"]
      tau_samples  <- .scale(tau_samples,  x$add_info[["output_scale"]], model_scale)

    }

    # selection models are sampled separately for increased efficiency
    bias_indicator           <- posterior_samples[,"bias_indicator"]
    weightfunction_indicator <- bias_indicator %in% which(sapply(priors[["bias"]], is.prior.weightfunction))

    # compute the density add specified support
    # outcome_lower     <- rep(NA, nrow(mu_samples))
    # outcome_higher    <- rep(NA, nrow(mu_samples))
    outcome_densities <- matrix(NA, nrow = nrow(mu_samples), ncol = length(z_sequence))

    for(j in 1:nrow(mu_samples)){
      # create containers for temporal samples from the posterior distribution
      # temp_lower     <- rep(NA, ncol(mu_samples))
      # temp_higher    <- rep(NA, ncol(mu_samples))
      temp_densities <- matrix(NA, nrow = ncol(mu_samples), ncol = length(z_sequence))

      for(i in seq_len(ncol(mu_samples))){

        # sample normal models/PET/PEESE
        if(!weightfunction_indicator[j]){
          # temp_lower[i]      <- stats::pnorm(z_sequence[1]                  * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
          # temp_higher[i]     <- stats::pnorm(z_sequence[length(z_sequence)] * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE)
          temp_densities[i,] <- stats::dnorm(z_sequence * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2)) * newdata.outcome[i,"se"]
        }

        # sample selection models
        if(weightfunction_indicator[j]){
          # # if including probability > & < than plotting range, those would also need to be re-standardized
          # temp_lower[i]  <- .pwnorm_fast.ss(
          #   q          = z_sequence[1] * newdata.outcome[i,"se"],
          #   mean       = mu_samples[j,i],
          #   sd         = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
          #   omega      = posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
          #   crit_x     = fit_data$crit_y[, i, drop=FALSE],
          #   lower.tail = TRUE)
          # temp_higher[i] <- .pwnorm_fast.ss(
          #   q          = z_sequence[length(z_sequence)] * newdata.outcome[i,"se"],
          #   mean       = mu_samples[j,i],
          #   sd         = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
          #   omega      = posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
          #   crit_x     = fit_data$crit_y[, i, drop=FALSE],
          #   lower.tail = FALSE)


          # the density needs to be transformed due to the support change
          # (+ extrapolation by removing standardizing constant in case of selection models)
          if(extrapolate){
            temp_consts <- .dwnorm_fast.ss(
              x      = 0,
              mean   = mu_samples[j,i],
              sd     = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
              omega  = posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
              crit_x = fit_data$crit_y[, i],
              attach_constant = TRUE
            )
            temp_densities[i,] <- stats::dnorm(z_sequence * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2)) * newdata.outcome[i,"se"] / attr(temp_consts, "constant")
          }else{
            temp_out <- .dwnorm_fast.ss(
              x      = z_sequence * newdata.outcome[i,"se"],
              mean   = mu_samples[j,i],
              sd     = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
              omega  = matrix(posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
                              nrow = length(z_sequence), ncol = length(grep("omega", colnames(posterior_samples))), byrow = TRUE),
              crit_x = fit_data$crit_y[, i],
              attach_constant = TRUE
            )
            temp_densities[i,] <- temp_out * newdata.outcome[i,"se"]
          }

        }
      }

      # store the results
      # outcome_lower[j]      <- mean(temp_lower)
      # outcome_higher[j]     <- mean(temp_higher)
      outcome_densities[j,] <- colMeans(temp_densities)
    }
  }

  return(outcome_densities)
}
