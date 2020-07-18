#' @title Prints a fitted RoBMA object
#'
#' @param x a fitted RoBMA object.
#' @param ... additional arguments.
#' @export  print.RoBMA
#' @rawNamespace S3method(print, RoBMA)
#' @seealso [RoBMA()]
print.RoBMA <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimates:\n")
  print(stats::coef(x))
}


#' @title Summarize fitted RoBMA object
#'
#' @description \code{summary.RoBMA} creates a numerical
#' summary of the RoBMA object.
#'
#' @param object a fitted RoBMA object.
#' @param type whether to show the overall RoBMA results (\code{"ensemble"}),
#' an overview of the individual models (\code{"models"}), or detailed summary
#' for the individual models (\code{"individual"}).
#' @param conditional show the conditional estimates (assuming that the
#' alternative is true). Defaults to \code{FALSE}. Only available for
#' \code{type == "conditional"}.
#' @param diagnostics show the maximum R-hat and minimum ESS for the main
#' parameters in each of the models. Only available for \code{type = "ensemble"}.
#' @param include_theta whether the estimated random effects should be included
#' either in the summaries.
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .50, .975)}
#' @param logBF show log of the BFs. Defaults to \code{FALSE}.
#' @param BF01 show BF in support of the null hypotheses. Defaults to
#' \code{FALSE}.
#' @param digits_estimates a number of decimals for rounding the estimates.
#' Defaults to \code{3}.
#' @param digits_BF a number of decimals for rounding the BFs. Defaults to \code{3}.
#' @param ... additional arguments
#'
#' @return summary of a RoBMA object
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' # summary can provide many details about the model
#' summary(fit)
#'
#' # note that the summary function contains additional arguments
#' # that allow to obtain a specific output, i.e, the conditional estimates
#' # (assuming that the non-null models are true) can be obtained
#' summary(fit, conditional = TRUE)
#'
#' # overview of the models and their prior and posterior probability, marginal likelihood,
#' # and inclusion Bayes factor:
#' summary(fit, type = "models")
#'
#' # and the model diagnostics overview, containing maximum R-hat and minimum ESS across parameters
#' # but see '?diagnostics' for diagnostics plots for individual model parameters
#' summary(fit, type = "models", diagnostics = TRUE)
#'
#' # summary of individual models and their parameters can be further obtained by
#' summary(fit, type = "individual")
#'
#' }
#' @note See [diagnostics()] for visual convergence checks of the individual models.
#' @method summary RoBMA
#' @export summary.RoBMA
#' @rawNamespace S3method(summary, RoBMA)
#' @seealso [RoBMA()] [diagnostics()]
summary.RoBMA       <- function(object, type = if(diagnostics) "models" else "ensemble", conditional = FALSE, diagnostics = FALSE, include_theta = FALSE,
                                probs = c(.025, .975), logBF = FALSE, BF01 = FALSE,
                                digits_estimates = 3, digits_BF = 3,...){

  # print diagnostics if all models fail to converge
  if(!any(object$add_info$converged)){
    if(substr(type,1,1) != "m" & !diagnostics)warning("All models failed to converge. Model diagnostics were printed instead.")
    type        <- "models"
    diagnostics <- TRUE
  }


  if(substr(type,1,1) == "e"){

    ### model estimates
    # compute quantiles
    if(!is.null(probs)){
      if(length(probs) != 0){

        if(!is.numeric(probs) | !is.vector(probs))stop("The passed probabilities 'probs' must be a numeric vector.")
        if(!(all(probs > 0) & all(probs < 1)))stop("The passed probabilities 'probs' must be higher than 0 and lower than 1.")

        # average median estimates
        averaged_q <- rbind(
          mu     = if(length(object$RoBMA$samples$averaged$mu)  == 0) c(NA, NA) else unname(stats::quantile(object$RoBMA$samples$averaged$mu,  probs)),
          tau    = if(length(object$RoBMA$samples$averaged$tau) == 0) c(NA, NA) else unname(stats::quantile(object$RoBMA$samples$averaged$tau, probs)))
        if(ncol(object$RoBMA$samples$averaged$omega) != 0)averaged_q <- rbind(averaged_q, matrix(apply(object$RoBMA$samples$averaged$omega, 2, stats::quantile, probs = probs), ncol = length(probs), byrow = TRUE))
        if(include_theta & nrow(object$RoBMA$samples$averaged$theta) != 0)averaged_q <- rbind(averaged_q, matrix(apply(object$RoBMA$samples$averaged$theta, 2, stats::quantile, probs = probs), ncol = length(probs), byrow = TRUE))


        # conditional mean estimates
        conditional_q <- rbind(
          mu     = if(length(object$RoBMA$samples$conditional$mu)  == 0) c(NA, NA) else unname(stats::quantile(object$RoBMA$samples$conditional$mu,  probs)),
          tau    = if(length(object$RoBMA$samples$conditional$tau) == 0) c(NA, NA) else unname(stats::quantile(object$RoBMA$samples$conditional$tau, probs)))
        if(ncol(object$RoBMA$samples$averaged$omega) != 0)conditional_q <- rbind(conditional_q, matrix(apply(object$RoBMA$samples$conditional$omega, 2, stats::quantile, probs = probs), ncol = length(probs), byrow = TRUE))
        if(include_theta & nrow(object$RoBMA$samples$conditional$theta) != 0)conditional_q <- rbind(conditional_q, matrix(apply(object$RoBMA$samples$conditional$theta, 2, stats::quantile, probs = probs), ncol = length(probs), byrow = TRUE))


        colnames(averaged_q)    <- probs
        colnames(conditional_q) <- probs

      }else{
        averaged_q    <- NULL
        conditional_q <- NULL
      }
    }else{
      averaged_q    <- NULL
      conditional_q <- NULL
    }

    # averaged mean estimates
    averaged_e <- c(
      mu     = if(length(object$RoBMA$samples$averaged$mu)  == 0) NA else base::mean(object$RoBMA$samples$averaged$mu),
      tau    = if(length(object$RoBMA$samples$averaged$tau) == 0) NA else base::mean(object$RoBMA$samples$averaged$tau),
      if(ncol(object$RoBMA$samples$averaged$omega) == 0) NULL else apply(object$RoBMA$samples$averaged$omega, 2, base::mean),
      if(include_theta){if(nrow(object$RoBMA$samples$averaged$theta) == 0) NULL else apply(object$RoBMA$samples$averaged$theta, 2, base::mean)}
    )

    # average median estimates
    averaged_m <- c(
      mu     = if(length(object$RoBMA$samples$averaged$mu)  == 0) NA else stats::median(object$RoBMA$samples$averaged$mu),
      tau    = if(length(object$RoBMA$samples$averaged$tau) == 0) NA else stats::median(object$RoBMA$samples$averaged$tau),
      if(ncol(object$RoBMA$samples$averaged$omega) == 0) NULL else apply(object$RoBMA$samples$averaged$omega, 2, stats::median),
      if(include_theta){if(nrow(object$RoBMA$samples$averaged$theta) == 0) NULL else apply(object$RoBMA$samples$averaged$theta, 2, stats::median)}
    )

    # conditional mean estimates
    conditional_e <- c(
      mu     = if(length(object$RoBMA$samples$conditional$mu)  == 0) NA else base::mean(object$RoBMA$samples$conditional$mu),
      tau    = if(length(object$RoBMA$samples$conditional$tau) == 0) NA else base::mean(object$RoBMA$samples$conditional$tau),
      if(ncol(object$RoBMA$samples$conditional$omega) == 0) NULL else apply(object$RoBMA$samples$conditional$omega, 2, base::mean),
      if(include_theta){if(nrow(object$RoBMA$samples$conditional$theta) == 0) NULL else apply(object$RoBMA$samples$conditional$theta, 2, base::mean)}
    )

    # conditional median estimates
    conditional_m <- c(
      mu     = if(length(object$RoBMA$samples$conditional$mu)  == 0) NA else stats::median(object$RoBMA$samples$conditional$mu),
      tau    = if(length(object$RoBMA$samples$conditional$tau) == 0) NA else stats::median(object$RoBMA$samples$conditional$tau),
      if(ncol(object$RoBMA$samples$conditional$omega) == 0) NULL else apply(object$RoBMA$samples$conditional$omega, 2, stats::median),
      if(include_theta){if(nrow(object$RoBMA$samples$conditional$theta) == 0) NULL else apply(object$RoBMA$samples$conditional$theta, 2, stats::median)}
    )

    # create estimates tables
    averaged_tab    <- cbind.data.frame(Mean = averaged_m,    Median = averaged_m,    averaged_q)
    conditional_tab <- cbind.data.frame(Mean = conditional_e, Median = conditional_m, conditional_q)

    ### model types overview
    mm_mu       <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "mu"))[object$add_info$converged]
    mm_tau      <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "tau"))[object$add_info$converged]
    mm_omega    <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "omega"))[object$add_info$converged]

    # number of model types
    models_n    <- c(
      mu    = sum(mm_mu),
      tau   = sum(mm_tau),
      omega = sum(mm_omega)
    )

    # extract model weights
    prior_weights_all   <- sapply(object$models, function(m)m$prior_odds)[object$add_info$converged]
    # standardize model weights
    prior_weights_all   <- prior_weights_all / sum(prior_weights_all)
    # conditional model weights
    models_prior <- c(
      mu    <- sum(prior_weights_all[mm_mu]),
      tau   <- sum(prior_weights_all[mm_tau]),
      omega <- sum(prior_weights_all[mm_omega])
    )

    # conditional model posteriors
    posterior_weights_all <- object$RoBMA$posterior_prob$all[object$add_info$converged]
    models_posteriors     <- c(
      mu    <- sum(posterior_weights_all[mm_mu]),
      tau   <- sum(posterior_weights_all[mm_tau]),
      omega <- sum(posterior_weights_all[mm_omega])
    )

    # BF
    BF <- unlist(object$RoBMA$BF)
    BF[is.nan(BF)] <- NA
    if(BF01){
      BF <- 1/BF
    }else{
      BF <- BF
    }
    if(logBF){
      BF <- log(BF)
    }

    # create overview table
    overview_tab <- cbind.data.frame(models_n, models_prior, models_posteriors, BF)
    rownames(overview_tab) <- c("Effect", "Heterogeneity", "Pub. bias")
    colnames(overview_tab) <- c("Models", "Prior prob.", "Post. prob.",
                                paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")"))


    ### return results
    res <- list(
      call     = object$call,
      overview = overview_tab,
      averaged = averaged_tab,
      add_info = list(
        weight_type      = if(is.null(.get_omega_mapping(object$models, cuts_only = TRUE))) "none" else if(any(sapply(object$models, function(m)m$priors$omega$distribution) == "one.sided"))"one-sided" else "two-sided",
        n_models         = length(object$models),
        effect_size      = object$add_info$effect_size,
        mu_transform      = if(object$add_info$effect_size %in% c("r","OR"))object$add_info$mu_transform,
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        study_names      = object$add_info$study_names,
        type             = "ensemble",
        failed           = sum(!object$add_info$converged)
      )

    )

    if(conditional){
      res$conditional <- conditional_tab
    }

  }else if(substr(type,1,1) == "m"){

    priors_mu      <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$mu, silent = TRUE))
    priors_tau     <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$tau, silent = TRUE))
    priors_omega   <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$omega, silent = TRUE))
    prior_odds     <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_odds)
    prior_prob     <- prior_odds / sum(prior_odds)
    marg_lik       <- sapply(1:length(object$models), function(i)object$models[[i]]$marg_lik$logml)
    posterior_prob <- bridgesampling::post_prob(marg_lik, prior_prob = prior_prob)
    BF             <- sapply(1:length(object$models), function(i){
      temp_mm <- rep(FALSE, length(object$models))
      temp_mm[i] <- TRUE
      .inclusion_BF(prior_prob, posterior_prob, temp_mm)
    })

    BF[is.nan(BF)] <- NA
    if(BF01){
      BF <- 1/BF
    }else{
      BF <- BF
    }
    if(logBF){
      BF <- log(BF)
    }

    overview_tab <- cbind.data.frame(priors_mu, priors_tau, priors_omega, prior_prob, posterior_prob, marg_lik, BF,
                                     stringsAsFactors = FALSE)
    rownames(overview_tab) <- NULL
    colnames(overview_tab) <- c("Prior mu","Prior tau","Prior omega", "Prior prob.", "Post. prob.", "log(MargLik)",
                                paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")"))


    # add the summary model diagnostics
    if(diagnostics){

      diagnostics_tab <- overview_tab[,1:3]

      # extract max(R-hat) & min(ESS)
      diag_sum <- sapply(1:length(object$models), function(i){

        temp_x <- object$models[[i]]$fit

        if(length(temp_x) == 0 | any(class(object$models[[i]]$fit) %in% c("simpleError","error"))){

          return(c(NA, NA, NA))

        }else{

          s.x <- object$models[[i]]$fit_summary
          names.omegas <- rownames(s.x)[grepl("omega", rownames(s.x))]
          names.omegas <- names.omegas[-length(names.omegas)] # remove the last one since it's fixed
          if(include_theta){
            names.thetas <- rownames(s.x)[grepl("theta", rownames(s.x))]
          }
          s.x <- s.x[rownames(s.x) %in% c("mu", "tau", names.omegas, if(include_theta)names.thetas), ]

          if(length(dim(s.x)) == 2){
            return(c(
              MCerr = max(s.x[, 7]),
              Rhat  = max(s.x[,11]),
              ESS   = min(s.x[, 9])
            ))
          }else{
            return(c(
              MCerr = max(s.x[ 7]),
              Rhat  = max(s.x[11]),
              ESS   = min(s.x[ 9])
            ))
          }

        }
       })

      diagnostics_tab$"max(MCMC error)"  <- diag_sum[1,]
      diagnostics_tab$"min(ESS)"         <- diag_sum[3,]
      diagnostics_tab$"max(Rhat)"        <- diag_sum[2,]
    }


    res <- list(
      call     = object$call,
      overview = overview_tab,
      add_info = list(
        n_models         = length(object$models),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        effect_size      = object$add_info$effect_size,
        mu_transform      = if(object$add_info$effect_size %in% c("r","OR"))object$add_info$mu_transform,
        study_names      = object$add_info$study_names,
        type             = "models"
      )
    )

    if(diagnostics){
      res$diagnostics <- diagnostics_tab
    }


  }else if(substr(type, 1, 1) == "i"){

    overview_tabs <- list()

    prior_odds     <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_odds)
    prior_prob     <- prior_odds / sum(prior_odds)
    posterior_prob <- object$RoBMA$posterior_prob$all
    marg_lik       <- sapply(1:length(object$models), function(i)object$models[[i]]$marg_lik$logml)
    BF             <- sapply(1:length(object$models), function(i){
      temp_mm <- rep(FALSE, length(object$models))
      temp_mm[i] <- TRUE
      .inclusion_BF(prior_prob, posterior_prob, temp_mm)
    })

    BF[is.nan(BF)] <- NA
    if(BF01){
      BF <- 1/BF
    }else{
      BF <- BF
    }
    if(logBF){
      BF <- log(BF)
    }


    for(i in 1:length(object$models)){

      if(length(object$models[[i]]$fit) == 0 |  any(class(object$models[[i]]$fit) %in% c("simpleError","error"))){

        s.x <- NULL

      }else{

        s.x <- object$models[[i]]$fit_summary

        names.omegas <- rownames(s.x)[grepl("omega", rownames(s.x))]
        if(include_theta){
          names.thetas <- rownames(s.x)[grepl("theta", rownames(s.x))]
        }

        s.x <- s.x[rownames(s.x) %in% c("mu", "tau", names.omegas, if(include_theta)names.thetas), ]

        if(length(dim(s.x)) == 0){
          s.x <- data.frame(matrix(s.x, ncol = 11))
          rownames(s.x) <- rownames(object$models[[i]]$fit_summary)[rownames(object$models[[i]]$fit_summary) %in% c("mu", "tau", names.omegas, if(include_theta)names.thetas)]
        }else{
          s.x <- data.frame(s.x)
        }
        s.x <- s.x[,c(4, 5, 1, 2, 3, 7, 8, 9, 11)]


        colnames(s.x) <- c("Mean", "SD", ".025", "Median", ".975", "MCMC error", "Error % of SD", "ESS", "Rhat")

        if(object$models[[i]]$priors$omega$distribution != "point"){
          # rename omegas
          rownames(s.x)[grepl("omega", rownames(s.x))] <- paste0(
            "omega[",c(0, rev(object$models[[i]]$priors$omega$parameters$steps)), ",",
            c(rev(object$models[[i]]$priors$omega$parameters$steps), 1), "]"
          )
          # and order them
          s.x[c(1:nrow(s.x))[grepl("omega", rownames(s.x))],] <- s.x[rev(c(1:nrow(s.x))[grepl("omega", rownames(s.x))]),]
        }

        if(object$add_info$effect_size %in%  c("r", "OR")){
          s.x[grepl("mu",    rownames(s.x)),1:6] <- .transform(s.x[grepl("mu",    rownames(s.x)),1:6], object$add_info$effect_size, object$add_info$mu_transform)
          s.x[grepl("theta", rownames(s.x)),1:6] <- .transform(s.x[grepl("theta", rownames(s.x)),1:6], object$add_info$effect_size, object$add_info$mu_transform)
        }
      }


      overview_tabs[[i]] <- list(
        priors     = list(
          mu          = object$models[[i]]$priors$mu,
          tau         = object$models[[i]]$priors$tau,
          omega       = object$models[[i]]$priors$omega
        ),
        tab        = s.x,
        prior_prob = prior_prob[i],
        marg_lik   = marg_lik[i],
        posterior_prob = posterior_prob[i],
        BF         = BF[i],
        add_info   = list(
          weight_type   = if(object$models[[i]]$priors$omega$distribution == "point") "none" else if(object$models[[i]]$priors$omega$distribution == "one.sided")"one-sided" else "two-sided",
          effect_size   = object$add_info$effect_size,
          mu_transform  = if(object$add_info$effect_size %in% c("r","OR"))object$add_info$mu_transform,
          study_names   = object$add_info$study_names,
          BF_type       = paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")")
        )
      )
    }


    res <- list(
      call     = object$call,
      overview = overview_tabs,
      add_info = list(
        n_models         = length(object$models),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        effect_size      = object$add_info$effect_size,
        mu_transform      = if(object$add_info$effect_size %in% c("r","OR"))object$add_info$mu_transform,
        BF_type          = paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")"),
        study_names      = object$add_info$study_names,
        type             = "individual"
      )

    )

  }

  class(res) <- "summary.RoBMA"
  return(res)
}

#' @title Prints summary object for RoBMA method
#'
#' @param x a summary of a RoBMA object
#' @param ... additional arguments
#' @method print.summary RoBMA
#' @export print.summary.RoBMA
#' @rawNamespace S3method(print, summary.RoBMA)
#' @seealso [RoBMA()]
print.summary.RoBMA <- function(x, ...){

  # format the output before printing
  if(x$add_info$type == "ensemble"){

    overview <- x$overview
    overview$Models <- paste0(overview$Models,"/",  x$add_info$n_models - x$add_info$failed)
    overview[,2:3]  <- format(round(overview[,2:3], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    overview[,4]    <- format(round(overview[,4],   x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)

    # round the results (a loop is required to deal with NAs)
    averaged <- x$averaged
    for(i in 1:ncol(averaged)){
      averaged[,i] <- format(round(averaged[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    }
    if(!is.null(x$conditional)){
      conditional <- x$conditional
      for(i in 1:ncol(conditional)){
        conditional[,i] <- format(round(conditional[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }
    }

  }else if(x$add_info$type == "models"){

    overview <- x$overview
    overview[,4:6]  <- format(round(overview[,4:6], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    overview[,7]    <- format(round(overview[,7],   x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)

    if(!is.null(x$diagnostics)){
      diagnostics     <- x$diagnostics
      diagnostics[,4] <- format(round(diagnostics[,4], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      diagnostics[,5] <- round(diagnostics[,5])
      diagnostics[,6] <- format(round(diagnostics[,6], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    }
  }else if(x$add_info$type == "individual"){

    overview <- x$overview

    for(i in 1:length(overview)){

      if(!is.null(overview[[i]]$tab)){
        overview[[i]]$tab[,c(1:6,9)] <- format(round(overview[[i]]$tab[,c(1:6,9)], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
        overview[[i]]$tab[,7]        <- format(round(overview[[i]]$tab[,7], 1), nsmall = 1)
        overview[[i]]$tab[,8]        <- round(overview[[i]]$tab[,8])
      }

      overview[[i]]$prior_prob     <- format(round(overview[[i]]$prior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$marg_lik       <- format(round(overview[[i]]$marg_lik, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$posterior_prob <- format(round(overview[[i]]$posterior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$BF             <- format(round(overview[[i]]$BF, x$add_info$digits_BF), nsmall = x$add_info$digits_BF)
    }
  }


  cat("Call:\n")
  print(x$call)
  cat("\n")


  if(x$add_info$type == "ensemble"){

    cat("Robust Bayesian Meta-Analysis\n")
    print(overview, quote = FALSE, right = TRUE)
    cat("\n")

    cat("Model-averaged estimates\n")
    print(averaged, quote = FALSE, right = TRUE)
    if(x$add_info$effect_size == "r"){
      cat(paste0("(Tau is on ", if(x$add_info$mu_transform == "cohens_d") "Cohen's d" else if(x$add_info$mu_transform == "fishers_z") "Fisher's z", " scale.)\n"))
    }else if(x$add_info$effect_size == "OR"){
      cat(paste0("(Tau is on ", if(x$add_info$mu_transform == "log_OR") "log(OR)" else if(x$add_info$mu_transform == "cohens_d") "Cohen's d", " scale.)\n"))
    }
    if(x$add_info$weight_type != "none")cat(paste0("(Estimated omegas correspond to ", x$add_info$weight_type, " p-values)\n"))
    if(x$add_info$failed != 0)cat(paste0("\033[0;31m",x$add_info$failed, ifelse(x$add_info$failed == 1, " model", " models"), " failed to converge and ",ifelse(x$add_info$failed == 1, "was", "were")," omited from the summary.\033[0m\n"))


    if(!is.null(x$conditional)){
      cat("\n")
      cat("Conditional estimates\n")
      print(conditional, quote = FALSE, right = TRUE)
      if(x$add_info$effect_size == "r"){
        cat(paste0("(Tau is on ", if(x$add_info$mu_transform == "cohens_d") "Cohen's d" else if(x$add_info$mu_transform == "fishers_z") "Fisher's z", " scale.)\n"))
      }else if(x$add_info$effect_size == "OR"){
        cat(paste0("(Tau is on ", if(x$add_info$mu_transform == "log_OR") "log(OR)" else if(x$add_info$mu_transform == "cohens_d") "Cohen's d", " scale.)\n"))
      }
      if(x$add_info$weight_type != "none")cat(paste0("(Estimated omegas correspond to ", x$add_info$weight_type, " p-values)\n"))
    }

  }else if(x$add_info$type == "models"){

    cat("Robust Bayesian Meta-Analysis\n")
    print(overview, quote = FALSE, right = TRUE)

    if(!is.null(x$diagnostics)){
      cat("\n")
      cat("Models diagnostics overview\n")
      print(diagnostics, quote = FALSE, right = TRUE)
    }

  }else if(x$add_info$type == "individual"){

    cat("Individual Models Summary\n\n")

    for(i in 1:length(overview)){

      cat(paste0("Model: ", i, "\n"))
      cat(paste0("Prior prob.:\t",         overview[[i]]$prior_prob,     "\t\t",
                                              "Prior mu:\t",print(overview[[i]]$priors$mu, silent = T),"\n"))
      cat(paste0("log(MargLik):\t",        overview[[i]]$marg_lik,       ifelse(nchar(as.character(overview[[i]]$marg_lik)) < 8, "\t\t", "\t"),
                                              "Prior tau:\t",print(overview[[i]]$priors$tau, silent = T),"\n"))
      cat(paste0("Post. prob.:\t",        overview[[i]]$posterior_prob, "\t\t",
                                              "Prior omega:\t",print(overview[[i]]$priors$omega, silent = T),"\n"))
      cat(paste0(x$add_info$BF_type,":\t", overview[[i]]$BF,        "\t\t","\n"))

      cat("\nModel Coefficients:\n")
      print(overview[[i]]$tab, quote = FALSE, right = TRUE)
      if(any(rownames(overview[[i]]$tab) == "tau")){
        if(x$add_info$effect_size == "r"){
          cat(paste0("(Tau is on ", if(x$add_info$mu_transform == "cohens_d") "Cohen's d" else if(x$add_info$mu_transform == "fishers_z") "Fisher's z", " scale.)\n"))
        }else if(x$add_info$effect_size == "OR"){
          cat(paste0("(Tau is on ", if(x$add_info$mu_transform == "log_OR") "log(OR)" else if(x$add_info$mu_transform == "cohens_d") "Cohen's d", " scale.)\n"))
        }
      }
      if(any(grepl("omega", rownames(overview[[i]]$tab)))){
        cat(paste0("(Estimated omegas correspond to ", overview[[i]]$add_info$weight_type, " p-values)\n"))
      }
      if(i != length(overview))cat("\n\n")
    }
  }

}

#' @title Reports whether x is a RoBMA object
#'
#' @param x an object to test
#' @export is.RoBMA
is.RoBMA            <- function(x){
  inherits(x, "RoBMA")
}

