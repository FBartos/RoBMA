
### basic plotting functions ----
#' @title Plots brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) distribution a brma object.
#'
#' @param x a fitted RoBMA object
#' @param parameter a parameter to be plotted. Defaults to
#' \code{"mu"} (for the effect size). The additional options
#' are \code{"tau"} (for the heterogeneity),
#' \code{"weightfunction"} (for the estimated weightfunction),
#' or \code{"PET-PEESE"} (for the PET-PEESE regression).
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
plot.brma  <- function(
    x, parameter, parameter_mods, parameter_scale,
    prior = FALSE, plot_type = "base", dots_prior = NULL, ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter       = parameter,
    parameter_mods  = parameter_mods,
    parameter_scale = parameter_scale,
    object          = x
  )

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = parameter
  )

  ### set up plotting arguments
  n_levels   <- .get_samples_n_levels(samples, parameter)
  dots       <- .set_dots_plot(..., n_levels = n_levels)
  dots_prior <- .set_dots_prior(dots_prior)

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- parameter
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- TRUE
  args$show_figures             <- NULL

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


#' @title Plots Weight Function of brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) weight function of a brma object.
#'
#' @param x a fitted RoBMA object
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
plot_weightfunction <- function(x, ...)  UseMethod("plot_weightfunction")

#' @export
#' @rdname plot_weightfunction
plot_weightfunction.brma  <- function(
    x, rescale_p_values = TRUE, show_data = TRUE,
    prior = FALSE, plot_type = "base", dots_prior = NULL, ...) {

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(rescale_p_values, "rescale_p_values")
  BayesTools::check_bool(show_data, "show_data")

  if (!.is_weightfunction(x)) {
    stop("'plot_weightfunction' is available only for models with a weightfunction component.", call. = FALSE)
  }

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = "omega"
  )

  ### set up plotting arguments
  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- "weightfunction"
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- FALSE
  args$rescale_x                <- rescale_p_values

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


#' @title Plots Weight Function of brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) weight function of a brma object.
#'
#' @param x a fitted RoBMA object
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
plot_PETPEESE <- function(x, ...)  UseMethod("plot_PETPEESE")

#' @export
#' @rdname plot_PETPEESE
plot_PETPEESE.brma  <- function(
    x, show_data = TRUE,
    prior = FALSE, plot_type = "base", dots_prior = NULL, ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(show_data, "show_data")

  if (!(.is_PET(x) || .is_PEESE(x))) {
    stop("'plot_PETPEESE' is available only for models with a PET or PEESE component.", call. = FALSE)
  }

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = c("mu", if (.is_PET(x)) "PET", if (.is_PEESE(x)) "PEESE")
  )

  ### set up plotting arguments
  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)

  # set plotting range (make sure all effect sizes are within the plotting range)
  data_outcome <- x[["data"]][["outcome"]]
  if (is.null(dots[["ylim"]])) {
    dots[["ylim"]] <- range(pretty(c(0, range(data_outcome$yi))))
  }

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- "PETPEESE"
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- FALSE

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # include data
  if (show_data) {
    if (plot_type == "ggplot") {
      plot <- plot + ggplot2::geom_point(
        ggplot2::aes(
          x  = data_outcome$sei,
          y  = data_outcome$yi),
        size  = if(!is.null(dots[["cex"]])) dots[["cex"]] else 2,
        shape = if(!is.null(dots[["pch"]])) dots[["pch"]] else 16
      )
    } else if(plot_type == "base") {
      graphics::points(data_outcome$sei, data_outcome$yi, cex = if(!is.null(dots[["cex"]])) dots[["cex"]] else 1, pch = if(!is.null(dots[["pch"]])) dots[["pch"]] else 16)
    }
  }

  # return the plots
  if (plot_type == "base") {
    return(invisible(plot))
  } else if (plot_type == "ggplot") {
    return(plot)
  }
}



.clean_PET_PEESE_samples <- function(samples, parameter) {
  for (i in seq_along(attr(samples[["PET"]], "prior_list"))) {
    if (is.prior.PET(attr(samples[["PET"]], "prior_list")[[i]])) {
      class(attr(samples[["PET"]], "prior_list")[[i]])   <- class(attr(samples[["PET"]], "prior_list")[[i]])[!class(attr(samples[["PET"]], "prior_list")[[i]]) %in% "prior.PET"]
    } else if (!is.prior.point(attr(samples[["PET"]], "prior_list")[[i]])){
      attr(samples[["PET"]], "prior_list")[[i]] <- prior("point", parameter = list("location" = 0), prior_weights = attr(samples[["PET"]], "prior_list")[[i]][["prior_weights"]])
    }
  }
  for (i in seq_along(attr(samples[["PEESE"]], "prior_list"))) {
    if (is.prior.PEESE(attr(samples[["PEESE"]], "prior_list")[[i]])) {
      class(attr(samples[["PEESE"]], "prior_list")[[i]])   <- class(attr(samples[["PEESE"]], "prior_list")[[i]])[!class(attr(samples[["PEESE"]], "prior_list")[[i]]) %in% "prior.PEESE"]
    } else if (!is.prior.point(attr(samples[["PEESE"]], "prior_list")[[i]])) {
      attr(samples[["PEESE"]], "prior_list")[[i]] <- prior("point", parameter = list("location" = 0), prior_weights = attr(samples[["PEESE"]], "prior_list")[[i]][["prior_weights"]])
    }
  }
  return(samples)
}




### basic MCMC diagnostics functions ----

#' @export
diagnostic_plots <- function(x, ...) UseMethod("diagnostic_plots")

#' @export
#' @rdname diagnostic_plots
diagnostic_plots.brma <- function(x, parameter, parameter_mods, parameter_scale, type, plot_type = "base", lags = 30, ...){

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_char(type, "type", allow_values = c("autocorrelation", "trace", "density"))

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter       = parameter,
    parameter_mods  = parameter_mods,
    parameter_scale = parameter_scale,
    object          = x
  )

  # a message with info about multiple plots
  if(plot_type == "base" && parameter == "omega") {
    message("Multiple plots will be produced. See '?layout' for help with setting multiple plots.")
  }

  dots  <- .set_dots_diagnostics(..., type = type, chains = x[["fit_control"]][["chains"]])

  # get the parameter name
  args                   <- dots
  args$fit               <- x[["fit"]]
  args$parameter         <- parameter
  args$type              <- type
  args$plot_type         <- plot_type
  args$lags              <- lags
  args$transform_factors <- TRUE
  args$short_name        <- FALSE
  args$parameter_names   <- FALSE
  args$formula_prefix    <- FALSE

  plots <- do.call(BayesTools::JAGS_diagnostics, args)

  # return the plots
  if(plot_type == "base"){
    return(invisible(plots))
  }else if(plot_type == "ggplot"){
    return(plots)
  }
}

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_autocorrelation <- function(x, ...) UseMethod("diagnostic_plots_autocorrelation")

#' @rdname diagnostic_plots
diagnostic_plots_autocorrelation.brma <- function(x, parameter = NULL, plot_type = "base", lags = 30, ...) {
  diagnostic_plots(x = x, parameter = parameter, type = "autocorrelation", plot_type = plot_type, lags = 30, ...)
}

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_trace <- function(x, ...) UseMethod("diagnostic_plots_trace")

#' @rdname diagnostic_plots
diagnostic_plots_trace.brma           <- function(x, parameter = NULL, plot_type = "base", ...) {
  diagnostic_plots(x = x, parameter = parameter, type = "trace", plot_type = plot_type, ...)
}

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_density <- function(x, ...) UseMethod("diagnostic_plots_density")

#' @rdname diagnostic_plots
diagnostic_plots_density.brma         <- function(x, parameter = NULL, plot_type = "base", ...) {
  diagnostic_plots(x = x, parameter = parameter, type = "density", plot_type = plot_type, ...)
}

### prediction functions ----



### funnel plot functions ----
#' @export
funnel <- function(x, ...) UseMethod("funnel")

#' @title Funnel Plot for brma Object
#'
#' @description \code{funnel.brma} creates a funnel plot for
#' a \code{"brma"} object. This function visualizes the sampling
#' distribution of observed effect sizes, incorporating heterogeneity
#' and publication bias if present in the model.
#'
#' @param x a fitted brma object
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param incorporate_heterogeneity Whether heterogeneity should be incorporated
#' into the sampling distribution. Defaults to \code{TRUE}.
#' @param incorporate_publication_bias Whether publication bias should be incorporated
#' into the sampling distribution. Defaults to \code{TRUE}.
#' @param max_samples Maximum number of samples from the posterior distribution
#' that will be used for estimating the funnel plot under publication bias.
#' Defaults to \code{500}.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function.
#'
#' @details
#' The \code{funnel.brma} function differs from the corresponding
#' \link[metafor]{funnel} function in two regards: 1) The heterogeneity is
#' by default incorporated into the model and 2) the sampling distribution
#' under publication bias is approximated by sampling from the estimated
#' weighted normal distribution under the pooled effect.
#'
#' The sampling distribution is drawn under the mean effect size and heterogeneity
#' estimates (the uncertainty about those values is not incorporated into the figure).
#'
#' @return \code{funnel.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat, seed = 1)
#'
#' # create funnel plot
#' funnel(fit)
#'
#' # funnel plot without publication bias adjustment
#' funnel(fit, incorporate_publication_bias = FALSE)
#'
#' # funnel plot without heterogeneity
#' funnel(fit, incorporate_heterogeneity = FALSE)
#'
#' # ggplot2 version
#' funnel(fit, plot_type = "ggplot")
#' }
#'
#' @export
#' @rdname funnel
funnel.brma <- function(x,
                        plot_type = "base",
                        incorporate_heterogeneity = TRUE,
                        incorporate_publication_bias = TRUE,
                        max_samples = 500,
                        ...) {

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(incorporate_heterogeneity, "incorporate_heterogeneity")
  BayesTools::check_bool(incorporate_publication_bias, "incorporate_publication_bias")
  BayesTools::check_int(max_samples, "max_samples", lower = 10)

  # check model type - only normal outcome models supported
  outcome_type <- .outcome_type(x)
  if (outcome_type != "norm") {
    stop("Funnel plots are only available for normal outcome models (not binomial or Poisson).", call. = FALSE)
  }

  ### extract structural information
  is_mods           <- .is_mods(x)
  is_multilevel     <- .is_multilevel(x)
  is_weightfunction <- .is_weightfunction(x)
  effect_direction  <- .effect_direction(x)

  ### compute residuals: observed yi - predicted mu (type="terms")
  # predictions are on the original scale
  res_samples <- predict.brma(
    object        = x,
    newdata       = NULL,
    type          = "terms",
    bias_adjusted = TRUE,
    as_samples    = TRUE,
    quiet         = TRUE
  )

  # compute residuals for each observation: yi - mu_predicted
  outcome_data <- x[["data"]][["outcome"]]
  K            <- nrow(outcome_data)

  # res_samples is S x K matrix, compute residuals
  residuals_matrix <- matrix(NA_real_, nrow = nrow(res_samples), ncol = K)
  for (k in seq_len(K)) {
    residuals_matrix[, k] <- outcome_data$yi[k] - res_samples[, k]
  }

  # use mean residuals for visualization
  res <- colMeans(residuals_matrix)

  # extract standard errors
  sei <- outcome_data$sei

  ### set up plotting range
  se_range    <- pretty(c(0, max(sei)))
  se_sequence <- seq(0, max(se_range), length.out = 21)

  ### extract posterior samples and priors
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["fit"]]))
  priors            <- x[["priors"]]

  ### compute funnel boundaries
  if (!incorporate_publication_bias || !is_weightfunction) {

    # compute contours without publication bias (normal distribution)

    if (incorporate_heterogeneity) {
      # extract tau samples
      tau_result <- .evaluate.brma.tau(
        fit           = x[["fit"]],
        scale_data    = x[["data"]][["scale"]],
        scale_formula = if (.is_scale(x)) .create_fit_formula_list(data = x[["data"]], "scale") else NULL,
        scale_priors  = priors[["scale"]],
        is_scale      = .is_scale(x),
        is_multilevel = is_multilevel,
        K             = 1L
      )
      # for funnel plot, use total tau (sqrt of sum of squared components)
      tau_within_samples  <- tau_result[["tau_within"]][, 1]
      tau_between_samples <- tau_result[["tau_between"]][, 1]
      tau_samples         <- sqrt(tau_within_samples^2 + tau_between_samples^2)
    } else {
      tau_samples <- rep(0, nrow(posterior_samples))
    }

    # compute normal quantiles
    ci_left  <- sapply(se_sequence, function(se) stats::qnorm(0.025, mean = 0, sd = mean(sqrt(se^2 + tau_samples^2))))
    ci_right <- sapply(se_sequence, function(se) stats::qnorm(0.975, mean = 0, sd = mean(sqrt(se^2 + tau_samples^2))))

  } else {

    # compute contours with publication bias (weighted normal distribution)

    # subsample for efficiency
    selected_samples_ind <- unique(round(seq(from = 1, to = nrow(posterior_samples), length.out = max_samples)))
    n_samples            <- length(selected_samples_ind)
    posterior_samples    <- posterior_samples[selected_samples_ind, , drop = FALSE]

    # get pooled mu samples
    mu_samples <- pooled_effect.brma(
      object        = x,
      bias_adjusted = FALSE,
      as_samples    = TRUE,
      quiet         = TRUE
    )
    mu_samples <- mu_samples[selected_samples_ind, 1]

    # get tau samples
    if (incorporate_heterogeneity) {
      tau_result <- .evaluate.brma.tau(
        fit           = x[["fit"]],
        scale_data    = x[["data"]][["scale"]],
        scale_formula = if (.is_scale(x)) .create_fit_formula_list(data = x[["data"]], "scale") else NULL,
        scale_priors  = priors[["scale"]],
        is_scale      = .is_scale(x),
        is_multilevel = is_multilevel,
        K             = 1L
      )
      tau_within_samples  <- tau_result[["tau_within"]][selected_samples_ind, 1]
      tau_between_samples <- tau_result[["tau_between"]][selected_samples_ind, 1]

      # for multilevel models, incorporate within-study heterogeneity into mu
      if (is_multilevel) {
        mu_samples <- mu_samples + stats::rnorm(n_samples) * tau_within_samples
        tau_samples <- tau_between_samples
      } else {
        tau_samples <- sqrt(tau_within_samples^2 + tau_between_samples^2)
      }
    } else {
      tau_samples <- rep(0, n_samples)
    }

    # extract omega samples for weight function
    omega_samples <- posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE]

    # get critical values for weight function
    steps <- BayesTools::weightfunctions_mapping(
      priors[["bias"]][sapply(priors[["bias"]], is.prior.weightfunction)],
      cuts_only = TRUE,
      one_sided = TRUE
    )
    steps <- rev(steps)[c(-1, -length(steps))]

    # compute quantiles at specific standard errors
    ci_left  <- rep(NA_real_, length(se_sequence))
    ci_right <- rep(NA_real_, length(se_sequence))

    for (j in seq_along(se_sequence)) {

      # compute critical y values for this se
      crit_y <- .get_cutoffs(
        y                = mu_samples,
        se               = rep(se_sequence[j], n_samples),
        prior            = list(distribution = "one.sided", parameters = list(steps = steps)),
        original_measure = rep("yi", n_samples),
        effect_measure   = "yi"
      )

      # sample from weighted normal for each posterior sample
      temp_lower  <- rep(NA_real_, n_samples)
      temp_higher <- rep(NA_real_, n_samples)

      for (i in seq_len(n_samples)) {
        temp_lower[i] <- .qwnorm_fast.ss(
          p      = 0.025,
          mean   = mu_samples[i],
          sd     = sqrt(tau_samples[i]^2 + se_sequence[j]^2),
          omega  = omega_samples[i, , drop = FALSE],
          crit_x = crit_y[i, ]
        ) - mu_samples[i]
        temp_higher[i] <- .qwnorm_fast.ss(
          p      = 0.975,
          mean   = mu_samples[i],
          sd     = sqrt(tau_samples[i]^2 + se_sequence[j]^2),
          omega  = omega_samples[i, , drop = FALSE],
          crit_x = crit_y[i, ]
        ) - mu_samples[i]
      }

      ci_left[j]  <- mean(temp_lower)
      ci_right[j] <- mean(temp_higher)
    }
  }

  ### create plotting objects
  x_range <- pretty(c(ci_left, ci_right))
  df_points <- data.frame(
    x = res,
    y = sei
  )
  df_funnel <- data.frame(
    x = c(rev(ci_left),     ci_right),
    y = c(rev(se_sequence), se_sequence)
  )
  df_funnel_edge1 <- data.frame(
    x = ci_left,
    y = se_sequence
  )
  df_funnel_edge2 <- data.frame(
    x = ci_right,
    y = se_sequence
  )
  df_background <- data.frame(
    x = c(min(x_range),  max(x_range),  max(x_range),  min(x_range)),
    y = c(min(se_range), min(se_range), max(se_range), max(se_range))
  )

  ### get additional plotting options
  dots <- list(...)

  ### create the plot
  if (plot_type == "ggplot") {

    out <- ggplot2::ggplot() +
      ggplot2::geom_polygon(
        mapping = ggplot2::aes(
          x = df_background$x,
          y = df_background$y),
        fill    = "grey"
      ) +
      ggplot2::geom_polygon(
        mapping = ggplot2::aes(
          x = df_funnel$x,
          y = df_funnel$y),
        fill    = "white"
      ) +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = c(0, 0),
          y = range(se_range)),
        linetype = "dotted"
      ) +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = df_funnel_edge1$x,
          y = df_funnel_edge1$y),
        linetype = "dotted"
      ) +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = df_funnel_edge2$x,
          y = df_funnel_edge2$y),
        linetype = "dotted"
      ) +
      ggplot2::geom_point(
        mapping = ggplot2::aes(
          x = df_points$x,
          y = df_points$y),
        fill  = "black",
        shape = if (is.null(dots[["shape"]])) 21 else dots[["shape"]],
        size  = if (is.null(dots[["size"]]))  2  else dots[["size"]]
      )

    out <- out +
      ggplot2::scale_x_continuous(breaks = x_range, limits = range(x_range), name = gettext("Residual")) +
      ggplot2::scale_y_reverse(breaks = rev(se_range), limits = rev(range(se_range)), name = gettext("Standard Error"))

  } else if (plot_type == "base") {

    # Set up the plot area with reversed y-axis
    graphics::plot(
      NA, NA,
      xlim = range(x_range),
      ylim = rev(range(se_range)),
      xlab = gettext("Residual"),
      ylab = gettext("Standard Error"),
      type = "n",
      axes = FALSE
    )
    graphics::axis(1, at = x_range)
    graphics::axis(2, at = rev(se_range), las = 1)

    # Draw background polygon (grey)
    graphics::polygon(df_background$x, df_background$y, col = "grey", border = NA)

    # Draw funnel polygon (white)
    graphics::polygon(df_funnel$x, df_funnel$y, col = "white", border = NA)

    # Vertical dotted line at x = 0
    graphics::lines(c(0, 0), range(se_range), lty = "dotted")

    # Funnel edges (dotted lines)
    graphics::lines(df_funnel_edge1$x, df_funnel_edge1$y, lty = "dotted")
    graphics::lines(df_funnel_edge2$x, df_funnel_edge2$y, lty = "dotted")

    # Determine point shape and size
    point_shape <- if (is.null(dots[["shape"]])) 21 else dots[["shape"]]
    point_size  <- if (is.null(dots[["size"]]))  1  else dots[["size"]]

    # Plot points with black fill
    graphics::points(df_points$x, df_points$y, pch = point_shape, bg = "black", cex = point_size)
  }

  # return the plots
  if (plot_type == "base") {
    return(invisible())
  } else if (plot_type == "ggplot") {
    return(out)
  }
}

### Common plotting helper functions ----
# Internal helper function for selecting and validating a plot parameter
#
# Processes parameter selection for plotting based on three mutually exclusive arguments:
# parameter, parameter_mods, parameter_scale. Only one should be specified.
#
# @param parameter Base parameter name ("mu", "tau", "rho")
# @param parameter_mods Moderator parameter for location regression
# @param parameter_scale Moderator parameter for scale regression
# @param object A brma object
#
# @return A character string with the validated parameter name (with appropriate prefix).
#   Attributes:
#   - parameter_samples: The JAGS samples name
#   - type: One of "base", "mods", "scale" indicating the parameter source
.check_and_select_plot_parameter <- function(parameter, parameter_mods, parameter_scale, object) {

  # Check for model characteristics
  is_mods           <- .is_mods(object)
  is_scale          <- .is_scale(object)
  is_multilevel     <- .is_multilevel(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)

  # Determine which arguments were provided
  has_parameter       <- !missing(parameter)       && !is.null(parameter)
  has_parameter_mods  <- !missing(parameter_mods)  && !is.null(parameter_mods)
  has_parameter_scale <- !missing(parameter_scale) && !is.null(parameter_scale)

  # Count how many were specified

  n_specified <- sum(c(has_parameter, has_parameter_mods, has_parameter_scale))

  # Error if more than one specified
  if (n_specified > 1) {
    stop("Only one of 'parameter', 'parameter_mods', or 'parameter_scale' can be specified.", call. = FALSE)
  }

  # Default behavior if none specified
  if (n_specified == 0) {
    if (is_mods) {
      # Default to intercept for meta-regression models
      return("mu_intercept")
    } else {
      # Default to mu for simple models
      return("mu")
    }
  }

  # Process based on which argument was specified
  if (has_parameter) {

    # Validate base parameters
    BayesTools::check_char(parameter, "parameter")

    # Handle mu parameter
    if (parameter == "mu") {
      if (is_mods) {
        # mu becomes mu_intercept when mods are present
        return("mu_intercept")
      } else {
        return("mu")
      }
    }

    # Handle tau parameter
    if (parameter == "tau") {
      if (is_scale) {
        # tau becomes log_tau_intercept when scale is present
        return("log_tau_intercept")
      } else {
        return("tau")
      }
    }

    # Handle rho parameter (only for multilevel models)
    if (parameter == "rho") {
      if (!is_multilevel) {
        stop("The 'rho' parameter is only available for multilevel models.", call. = FALSE)
      }
      return("rho")
    }

    ### Handle publication bias parameters
    if (parameter == "PET") {
      if (!is_PET) {
        stop("The 'PET' parameter is only available for PET models.", call. = FALSE)
      }
      return("PET")
    }

    if (parameter == "PEESE") {
      if (!is_PEESE) {
        stop("The 'PEESE' parameter is only available for PEESE models.", call. = FALSE)
      }
      return("PEESE")
    }

    if (parameter == "omega" || parameter == "weightfunction") {
      if (!is_weightfunction) {
        stop("The 'omega'/'weightfunction' parameter is only available for selection models.", call. = FALSE)
      }
      return("omega")
    }

    # Unknown base parameter
    # Build list of valid parameters dynamically
    valid_params <- c("'mu'", "'tau'")
    if (is_multilevel)     valid_params <- c(valid_params, "'rho'")
    if (is_PET)            valid_params <- c(valid_params, "'PET'")
    if (is_PEESE)          valid_params <- c(valid_params, "'PEESE'")
    if (is_weightfunction) valid_params <- c(valid_params, "'omega'/'weightfunction'")

    stop(sprintf(
      "Unknown parameter '%s'. Valid base parameters are: %s.",
      parameter,
      paste(valid_params, collapse = ", ")
    ), call. = FALSE)

  } else if (has_parameter_mods) {

    # Validate that mods are present
    if (!is_mods) {
      stop("The 'parameter_mods' argument can only be used for models with moderators.", call. = FALSE)
    }

    BayesTools::check_char(parameter_mods, "parameter_mods")

    # Get available terms from the mods formula
    mods_formula    <- attr(object[["data"]][["mods"]], "formula")
    available_terms <- c("intercept", attr(stats::terms(mods_formula), "term.labels"))

    # Validate the specified term exists
    if (!parameter_mods %in% available_terms) {
      stop(sprintf(
        "The specified 'parameter_mods' term '%s' is not available. Available terms are: %s.",
        parameter_mods,
        paste0("'", available_terms, "'", collapse = ", ")
      ), call. = FALSE)
    }

    # Prefix with "mu_"
    return(paste0("mu_", parameter_mods))

  } else if (has_parameter_scale) {

    # Validate that scale is present
    if (!is_scale) {
      stop("The 'parameter_scale' argument can only be used for location-scale models.", call. = FALSE)
    }

    BayesTools::check_char(parameter_scale, "parameter_scale")

    # Get available terms from the scale formula
    scale_formula <- attr(object[["data"]][["scale"]], "formula")
    available_terms <- c("intercept", attr(stats::terms(scale_formula), "term.labels"))

    # Validate the specified term exists
    if (!parameter_scale %in% available_terms) {
      stop(sprintf(
        "The specified 'parameter_scale' term '%s' is not available. Available terms are: %s.",
        parameter_scale,
        paste0("'", available_terms, "'", collapse = ", ")
      ), call. = FALSE)
    }

    # Prefix with "log_tau_"
    return(paste0("log_tau_", parameter_scale))
  }
}

.set_dots_plot        <- function(..., n_levels = 1) {

  dots <- list(...)
  if (is.null(dots[["col"]]) & n_levels == 1) {
    dots[["col"]]      <- "black"
  }else if (is.null(dots[["col"]]) & n_levels > 1) {
    dots[["col"]]      <- grDevices::palette.colors(n = n_levels + 1, palette = "Okabe-Ito")[-1]
  }
  if (is.null(dots[["col.fill"]])) {
    dots[["col.fill"]] <- "#4D4D4D4C" # scales::alpha("grey30", .30)
  }

  return(dots)
}
.set_dots_prior       <- function(dots_prior) {

  if (is.null(dots_prior)) {
    dots_prior <- list()
  }

  if (is.null(dots_prior[["col"]])) {
    dots_prior[["col"]]      <- "grey60"
  }
  if (is.null(dots_prior[["lty"]])) {
    dots_prior[["lty"]]      <- 1
  }
  if (is.null(dots_prior[["col.fill"]])) {
    dots_prior[["col.fill"]] <- "#B3B3B34C" # scales::alpha("grey70", .30)
  }

  return(dots_prior)
}
.set_dots_diagnostics <- function(..., type, chains) {

  dots <- list(...)
  if (is.null(dots[["col"]])) {
    dots[["col"]] <- if(type == "autocorrelation") "black" else rev(scales::viridis_pal()(chains))
  }

  return(dots)
}
.get_samples_n_levels <- function(samples, parameter) {
  if (inherits(samples[[parameter]], "mixed_posteriors.factor")) {
    if (attr(samples[[parameter]],"orthonormal") || attr(samples[[parameter]],"meandif")) {
      n_levels <- length(attr(samples[[parameter]],"level_names"))
    } else if (attr(samples[[parameter]],"treatment")) {
      n_levels <- length(attr(x$add_info[["predictors"]],"level_names")) - 1
    }
  } else {
    n_levels <- 1
  }
}
