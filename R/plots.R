#' @title Plots a fitted RoBMA object
#'
#' @description \code{plot.RoBMA} allows to visualize
#' different \code{"RoBMA"} object parameters in various
#' ways. See \code{type} for the different model types.
#'
#' @param x a fitted RoBMA object
#' @param parameter a parameter to be plotted. Either
#' \code{"mu"}, \code{"tau"}, \code{"theta"}, or
#' \code{"omega"}. A bivariate plot for model-averaged
#' estimates of "mu" and "tau" can be obtained by
#' \code{c("mu","tau")} if \code{type = "averaged"}. In
#' addition a forest plot with the original estimates can
#' be obtained by \code{"forest"} or added to the theta
#' estimates by \code{c("theta", "forest")}.
#' @param type what type of estimates should be plotted.
#' Options are \code{"averaged"} for the model-averaged
#' estimates, \code{"conditional"} for the conditional
#' estimates, or \code{"individual"} for the
#' individual models estimates. The options
#' \code{c("individual", "conditional")} can be supplied
#' together to show only coditional individual models.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot2"} for plotting. The later
#' requires \pkg{ggplot2} package to be installed.
#' @param mean whether the mean should be plotted.
#' @param median whether the median should be plotted.
#' @param CI width of the confidence intervals.
#' @param prior add prior density to the plot. Only available
#' for \code{type = "averaged"} or \code{type = "conditional"}.
#' Defaults to \code{FALSE}.
#' @param order either (1) ordering of the studies for
#' \code{parameter = "theta"} or \code{parameter = "forest"}.
#' Defaults to \code{NULL} - ordering as supplied to the fitting
#' function. However, studies can be ordered either
#' \code{"ascending"} or \code{"descending"} by effect size, or
#' by \code{"alphabetical"} by labels.
#' Or (2) ordering models for \code{type = "individual"}. The
#' default orders models according to their number. However,
#' models can be ordered either \code{"ascending"} or
#' \code{"descending"} by posterior model probability
#' \code{c("ascending", "prob")}, or marginal likelihood
#' \code{c("descending", "marglik")}
#' by marginal likelihood.
#' @param digits_estimates number of decimals to be displays for
#' \code{parameter = "theta"}, \code{parameter = "forest"}, and
#' \code{type = "individual"} plot.
#' @param show_figures which figures should be returned in case when
#' multiple plots are generated. Useful when
#' \code{parameter = "omega", type = "individual"} which generates
#' a figure for each weights cut-off. Defaults to \code{-1} which
#' omits the first weight. Set to \code{NULL} to show all figures
#' or to \code{c(1,3)} to show only the first and third one.
#' @param weights whether the weights or weight function should
#' be returned. Only applicable when \code{parameter = "omega"}.
#' Defaults to \code{FALSE} - the weight function is plotted.
#' @param ... additional arguments to be passed to
#' \link[graphics]{par} if \code{plot_type = "base"}. Especially
#' useful for \code{parameter == "theta"},
#' \code{parameter == "forest"} or \code{type = "individual"}
#' where automatic margins might cut out parts of the labels.
#'
#' @examples \donttest{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # plot function allows to visualize the results of a fitted RoBMA object, for example,
#' # the model-averaged mean parameter estimate
#' plot(fit, parameter = "mu")
#'
#' # or show both the prior and posterior distribution
#' plot(fit, parameter = "mu", prior = TRUE)
#'
#' # condtional plots might by obtained by specifying
#' plot(fit, parameter = "mu", type = "conditional")
#'
#' # plotting function also allows to visualize the weight function (or individual weights by adding 'weights = TRUE')
#' plot(fit, parameter = "omega")
#'
#' # or the forest plot (the estimated study effects can be shown by setting 'parameter = "theta"')
#' plot(fit, parameter = "forest")
#'
#' # it is also possible to compare the individual model estimates, and order them by the posterior probability
#' plot(fit, parameter = "mu", type = "individual", order = "prob")
#'
#' }
#' @export  plot.RoBMA
#' @rawNamespace S3method(plot, RoBMA)
#' @seealso [RoBMA()]
plot.RoBMA <- function(x, parameter,
                       type = "averaged", plot_type = "base",
                       mean = TRUE, median = FALSE, CI = .95, prior = FALSE,
                       order = NULL, digits_estimates = 2,
                       show_figures = if(parameter == "omega" & (weights | any(type %in% "individual")) ) -1, weights = FALSE, ...){

  ### settings
  # deal with misspecified arguments
  if(length(parameter) == 1){
    if(parameter %in% c("weight", "weights"))parameter <- "omega"
  }

  for(i in 1:length(type)){
    if(substr(type[i], 1, 1) == "c")type[i] <- "conditional"
    if(substr(type[i], 1, 1) == "i")type[i] <- "individual"
    if(substr(type[i], 1, 1) == "m")type[i] <- "individual" # for 'models'
    if(substr(type[i], 1, 1) == "a")type[i] <- "averaged"
  }

  if(plot_type == "ggplot2")plot_type <- "ggplot"


  if(sum(x$add_info$converged) == 0)stop("There is no converged model in the ensemble.")
  if(length(type) > 2)stop("The passed type is not supported for plotting. See '?plot.RoBMA' for more details.")
  if(length(type) == 2)if(!all(type %in% c("conditional","individual")))stop("The passed type is not supported for plotting. See '?plot.RoBMA' for more details.")
  if(!all(type %in% c("averaged","conditional","individual")))stop("The passed type is not supported for plotting. See '?plot.RoBMA' for more details.")
  if(!plot_type %in% c("base", "ggplot"))stop("The passed plot_type is not supported for plotting. See '?plot.RoBMA' for more details.")

  # check availability of ggplot
  if(plot_type == "ggplot"){
    if(!try(requireNamespace("ggplot2")))stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  }


  ### plotting
  if(any(type %in% "individual")){

    if(any(parameter %in% c("theta", "forest")))stop("Individual plot is not supported for theta parameter (forest plot).")
    output <- .plot.RoBMA_ind(x, parameter, type, plot_type, mean, median, CI, order, digits_estimates, show_figures, ...)

  }else if(length(parameter) == 2){

    # option for creating model-averaged bivariate plots of mu and tau
    if(all(parameter %in% c("mu", "tau"))){

      if(type != "averaged")stop("Bivariate plot is available only for averaged parameters.")
      output <- .plot.RoBMA_biv(x, parameter, plot_type, ...)

    }else if(all(parameter %in% c("theta", "forest"))){

      # creating forrest plot with both observed and estimated effect sizes
      output <- .plot.RoBMA_theta(x, parameter, type, plot_type, mean, median, CI, order, digits_estimates, ...)

    }else{
      stop("The parameters combination is not supported for plotting. See '?plot.RoBMA' for more details.")
    }

  }else if(parameter %in% c("mu", "tau") | (parameter ==  "omega" & weights)){
    output <- .plot.RoBMA_par(x, parameter, type, plot_type, mean, median, CI, prior, show_figures, ...)
  }else if(parameter == "omega" & !weights){
    output <- .plot.RoBMA_weightf(x, type, plot_type, mean, median, CI, prior, ...)
  }else if(parameter %in% c("theta", "forest")){
    output <- .plot.RoBMA_theta(x, parameter, type, plot_type, mean, median, CI, order, digits_estimates, ...)
  }else{
    stop("The parameter is not supported for plotting. See '?plot.RoBMA' for more details.")
  }


  # return the plots
  if(plot_type == "base"){
    return(invisible(NULL))
  }else if(plot_type == "ggplot"){
    return(output)
  }
}

.plot.RoBMA_par     <- function(fit, par, type, plot_type, mean, median, CI, prior, show_figures, ...){

  # get the parameter values & density
  y  <- fit$RoBMA$samples[[type]][[par]]
  if(length(y) == 0)stop("The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")

  # make sure that the parameter is in a matrix
  if(length(dim(y)) == 0)y <- matrix(y, ncol = 1)

  # obtain model information
  fit_models   <- fit$models[fit$add_info$converged]
  model_priors <- sapply(fit_models, function(m)m$priors[[par]], simplify = FALSE)
  prior_odds   <- sapply(fit_models, function(m)m$prior_odds)
  if(type == "conditional"){
    mm_par       <- sapply(fit_models, function(m)!.is_parameter_null(m$priors, par))
    model_priors <- model_priors[mm_par]
    prior_odds   <- prior_odds[mm_par]
  }
  prior_prob   <- prior_odds / sum(prior_odds)
  if(par == "omega"){
    spikes_loc <- 1
  }else{
    spikes_loc   <- unique(unlist(sapply(model_priors, function(p)if(p$distribution == "point") p$parameters$location)))
  }
  if(par %in% c("mu","tau")){
    trunc_lower  <- unlist(sapply(model_priors, function(p)p$truncation$lower))
    trunc_upper  <- unlist(sapply(model_priors, function(p)p$truncation$upper))
    trunc_lower  <- if(!all(is.null(trunc_lower))) min(trunc_lower) else -Inf
    trunc_upper  <- if(!all(is.null(trunc_upper))) max(trunc_upper) else  Inf
  }

  # get the indicies
  if(par == "omega"){
    omega_ind <- .get_omega_mapping(fit$models)
    if(type == "conditional")omega_ind <- omega_ind[mm_par]
    if(is.null(omega_ind))stop("The ensemble does not cointain any converged model adjusting for publication bias.")
  }


  # transform values to correlations if needed
  if(fit$add_info$effect_size == "r" & par == "mu"){
    if(fit$add_info$mu_transform == "cohens_d"){
      spikes_loc <- psych::d2r(spikes_loc)
      trunc_lower = if(!is.infinite(trunc_lower)) psych::d2r(trunc_lower) else -Inf
      trunc_upper = if(!is.infinite(trunc_upper)) psych::d2r(trunc_upper) else  Inf
    }else if(fit$add_info$mu_transform == "fishers_z"){
      spikes_loc <- psych::fisherz2r(spikes_loc)
      trunc_lower = if(!is.infinite(trunc_lower)) psych::fisherz2r(trunc_lower) else -Inf
      trunc_upper = if(!is.infinite(trunc_upper)) psych::fisherz2r(trunc_upper) else  Inf
    }
  }


  # get the estimated densitites
  res_CI     <- NULL
  res_den    <- list()
  res_spikes <- list()

  for(i in 1:ncol(y)){
    # create a temporary y
    temp_y   <- y[,i]

    # get out the posterior on spikes
    res_spikes[[i]] <- list()
    if(length(spikes_loc) > 0){
      for(j in 1:length(spikes_loc)){
        res_spikes[[i]][[j]] <- data.frame(
          x = rep(spikes_loc[j], 2),
          y = c(0, sum(temp_y == spikes_loc[j])/nrow(y))
        )
        temp_y  <- temp_y[temp_y != spikes_loc[j]]
      }
    }

    # compute the density of the rest
    if(length(temp_y) > 0){
      if(par %in% c("mu", "tau")){

        with_trunc <- list()
        if(!is.infinite(trunc_lower))with_trunc$from <- trunc_lower
        if(!is.infinite(trunc_upper))with_trunc$to   <- trunc_upper

        temp_den <- do.call(stats::density, c(list(x = temp_y), with_trunc))

        # truncate the posterior density appropriatelly
        no_support <- .get_no_support(fit_models, par)
        if(!is.null(no_support)){
          for(n in 1:length(no_support)){
            temp_den$y[temp_den$x > no_support[[n]]$lower & temp_den$x < no_support[[n]]$upper] <- 0
          }
        }

      }else if(par == "omega"){
        temp_den <- stats::density(temp_y, from = 0, to = 1)
      }

      # rescale the density by the proportion of included samples
      temp_den$y   <- temp_den$y * (length(temp_y) / length(y[,i]))

      # add the start and end points
      temp_den$x <- c(temp_den$x[1], temp_den$x, temp_den$x[length(temp_den$x)])
      temp_den$y <- c(0,             temp_den$y, 0)

      res_den[[i]] <- data.frame(x = temp_den$x, y = temp_den$y)
    }else{
      res_den[[i]] <- list()
    }

    # the CIs are based on all samples
    res_CI       <- rbind(res_CI, c(lCI = stats::quantile(y[,i], (1-CI)/2), uCI = stats::quantile(y[,i], 1-(1-CI)/2)))
  }

  # get the prior densities
  if(prior){

    # get the plotting data
    prior_den    <- list()
    prior_spikes <- list()

    # get range
    if(par != "omega"){
      x_range <- NULL
      for(j in 1:length(model_priors)){
        if(model_priors[[j]]$distribution == "point"){
          x_range <- rbind(x_range, rep(model_priors[[j]]$parameters$location,2))
        }else{
          x_range <- rbind(x_range, .plot.prior_data(
            prior        = model_priors[[j]],
            mu_transform = if(fit$add_info$effect_size == "r" & par == "mu")fit$add_info$mu_transform,
            samples      = 100)$x_range)
        }
      }
      x_range <- c(min(c(x_range[,1], res_den[[1]]$x)), max(c(x_range[,2], res_den[[1]]$x)))
    }else{
      x_range <- c(0, 1)
    }

    # get the data
    temp_prior_dat <- list()
    for(j in 1:length(model_priors)){
      temp_prior_dat[[j]] <- .plot.prior_data(
        prior        = model_priors[[j]],
        x_range      = x_range,
        mu_transform = if(fit$add_info$effect_size == "r" & par == "mu")fit$add_info$mu_transform,
        weights      = TRUE)
    }

    for(i in 1:ncol(y)){

      # get densities across this range
      temp_den       <- NULL
      temp_spikes    <- NULL

      for(j in 1:length(model_priors)){

        # select the correct one in case of omega
        if(par == "omega"){
          # the indicies are reversed since .plot.prior_data generates the coefficients in reverse order
          temp_par_ind <- if(model_priors[[j]]$distribution == "point") 1 else -1* omega_ind[[j]][i] + max(omega_ind[[j]]) + 1
        }else{
          temp_par_ind <- 1
        }

        # add either spike or density
        if(temp_prior_dat[[j]]$prob[[temp_par_ind]]){
          temp_spikes$x <- c(temp_spikes$x, temp_prior_dat[[j]]$df[[temp_par_ind]]$x[2])
          temp_spikes$y <- c(temp_spikes$y, temp_prior_dat[[j]]$df[[temp_par_ind]]$y[2] * prior_prob[j])
        }else{
          temp_den$x <- cbind(temp_den$x, temp_prior_dat[[j]]$df[[temp_par_ind]]$x)
          temp_den$y <- cbind(temp_den$y, temp_prior_dat[[j]]$df[[temp_par_ind]]$y * prior_prob[j])
        }


      }

      # agregate
      prior_spikes[[i]] <- list()
      if(!is.null(temp_spikes)){
        for(j in 1:length(unique(temp_spikes$x))){
          prior_spikes[[i]][[j]] <- data.frame(
            x = rep(unique(temp_spikes$x[j]), 2),
            y = c(0, sum(temp_spikes$y[temp_spikes$x == unique(temp_spikes$x)[j]]))
          )
        }
      }
      if(!is.null(temp_den)){
        if(length(dim(temp_den$y)) == 0)temp_den$y <- matrix(temp_den$y, ncol = 1)
        prior_den[[i]] <- data.frame(
          x = temp_den$x[,1],
          y = apply(temp_den$y, 1, sum)
        )
      }else{
        prior_den[[i]] <- list()
      }
    }
  }


  # get parameter names
  if(par == "mu"){
    if(fit$add_info$effect_size == "r")par_names <- list(bquote(~rho~.(paste0("(",type,")"))))
    if(fit$add_info$effect_size == "d")par_names <- list(bquote("Cohen's"~italic(d)~.(paste0("(",type,")"))))
    if(fit$add_info$effect_size == "y")par_names <- list(bquote(~mu~.(paste0("(",type,")"))))
  }else if(par == "tau"){
    if(fit$add_info$effect_size == "r"){
      if(fit$add_info$mu_transform == "cohens_d"){
        par_names <- list(bquote(~tau~"(Cohen's"~italic(d)~ "scale;"~.(paste0(type,")"))))
      }else if(fit$add_info$mu_transform == "fishers_z"){
        par_names <- list(bquote(~tau~"(Fisher's"~italic(z)~ "scale;"~.(paste0(type,")"))))
      }
    }else{
      par_names <- list(bquote(~tau~.(paste0("(",type,")"))))
    }
  }else if(par == "omega"){

    summary_info <- summary(fit, conditional = type == "conditional")
    sum_all <- summary_info[[type]]

    par_names <- sapply(rownames(sum_all)[grepl(par, rownames(sum_all))], function(x){
      bquote(~omega[~.(substr(x,6,nchar(x)))]~.(paste0("(",summary_info$add_info$weight_type, ")")))
    })
  }


  # create appropriate plot
  plots <- list()
  if(is.null(show_figures) | length(res_den) == 1){
    plots_ind <- c(1:ncol(y))
  }else{
    plots_ind <- c(1:ncol(y))[show_figures]
  }
  # a message with info about muliple plots
  if(plot_type == "base" & length(plots_ind) > 1)cat(paste0(length(plots_ind), " plots will be produced. See '?layout' for help with setting multiple plots."))


  for(i in plots_ind){


    # get parameter range
    if(par %in% c("mu", "tau")){

      # get range of all x
      temp_xlim <- range(c(
        if(!is.null(res_den[[i]]))res_den[[i]]$x,
        if(length(res_spikes[[i]]) != 0)unlist(sapply(1:length(res_spikes[[i]]), function(j)res_spikes[[i]][[j]]$x[res_spikes[[i]][[j]]$y > 0], simplify = FALSE)),
        if(prior)if(!is.null(prior_den[[i]]))prior_den[[i]]$x,
        if(prior)if(length(prior_spikes[[i]]) != 0)unlist(sapply(1:length(prior_spikes[[i]]), function(j)prior_spikes[[i]][[j]]$x[prior_spikes[[i]][[j]]$y > 0], simplify = FALSE)),
        0
      )
      )

      # deal with collapsed range due to only point posterior density
      if(temp_xlim[1] == temp_xlim[2]){
        temp_xlim[1] <- temp_xlim[1] - .5
        temp_xlim[2] <- temp_xlim[2] + .5
      }

      if(par == "tau" & temp_xlim[1] < 0)temp_xlim[1] <- 0

    }else{
      temp_xlim <- c(0, 1)
    }

    # range of y can be defined only on the spikes, but we wanna wanna use spikes for scalling
    temp_all_y <- c(
      if(!is.null(res_den[[i]]))res_den[[i]]$y,
      if(prior)if(!is.null(prior_den[[i]]))prior_den[[i]]$y
    )
    if(length(temp_all_y) == 0 | all(temp_all_y == 0)){
      temp_ylim <- c(0, 1)
    }else{
      temp_ylim <- c(0, max(temp_all_y))
    }


    # rescale the spike probabilities
    any_spikes <- FALSE
    if(length(res_spikes[[i]]) != 0){
      for(j in 1:length(length(res_spikes[[i]]))){
        res_spikes[[i]][[j]]$y <- res_spikes[[i]][[j]]$y * temp_ylim[2]
        if(res_spikes[[i]][[j]]$y[2] > 0)any_spikes <- TRUE
      }
    }
    if(prior){
      if(length(prior_spikes[[i]]) != 0){
        prior_spikes[[i]][[j]]$y <- prior_spikes[[i]][[j]]$y * temp_ylim[2]
        if(prior_spikes[[i]][[j]]$y[2] > 0)any_spikes <- TRUE
      }
    }
    any_density <- FALSE
    if(!is.null(res_den[[i]]$x))any_density <- TRUE
    if(prior)if(!is.null(prior_den[[i]]$x))any_density <- TRUE


    if(plot_type == "base"){

      # set up margins
      if(length(list(...)) == 0){
        graphics::par(mar = c(4,4,1,if(any_spikes) 4 else 1))
      }else{
        graphics::par(list(...))
      }

      graphics::plot(NA, bty = "n", las = 1, xlim = temp_xlim, ylim = temp_ylim, xlab = "", ylab = "", main = "", type = "l",
                     yaxt = if(!any_density) "n")
      graphics::mtext(par_names[[i]], side = 1, line = 2.5, cex = 1.25)
      # add axis
      if(any_density){
        graphics::mtext("Density",      side = 2, line = 2.5, cex = 1.25)
        if(any_spikes)graphics::axis(4, at = seq(0, temp_ylim[2], length.out = 5), labels = c("0", ".25", ".50", ".75", "1"), las = 1)
        if(any_spikes)graphics::mtext("Probability",  side = 4, line = 2.5, cex = 1.25)
      }else{
        if(any_spikes)graphics::axis(2, at = seq(0, temp_ylim[2], length.out = 5), labels = c("0", ".25", ".50", ".75", "1"), las = 1)
        if(any_spikes)graphics::mtext("Probability",  side = 2, line = 2.5, cex = 1.25)
      }

      # maybe add CI shading and mean/median later
      if(FALSE){
        if(CI)graphics::polygon(
          x = c(res_CI[i,1], res_CI[i,1], res_den[[i]]$x[res_den[[i]]$x >= res_CI[i,1] & res_den[[i]]$x <= res_CI[i,2]], res_CI[i,2], res_CI[i,2]),
          y = c(0, res_den[[i]]$y[which.max(res_den[[i]]$x >= res_CI[i,1])], res_den[[i]]$y[res_den[[i]]$x >= res_CI[i,1] & res_den[[i]]$x <= res_CI[i,2]], res_den[[i]]$y[which.min(res_den[[i]]$x <= res_CI[i,2])], 0),
          col = "grey80"
        )
        if(mean)graphics::lines(rep(mean(y[,i]),2), c(0,res_den[[i]]$y[which.max(res_den[[i]]$x > mean(y[,i]))]) ,lwd = 2)
        if(median)graphics::lines(rep(median(y[,i]),2), c(0,res_den[[i]]$y[which.max(res_den[[i]]$x > median(y[,i]))]) ,lwd = 2)
      }

      # add densities
      if(prior)if(!is.null(prior_den[[i]]))graphics::lines(prior_den[[i]]$x, prior_den[[i]]$y, lwd = 2, lty = 2, col = "grey50")
      if(!is.null(res_den[[i]]))graphics::lines(res_den[[i]]$x, res_den[[i]]$y, lwd = 2, lty = 1)


      # add spikes
      if(prior){
        if(length(prior_spikes[[i]]) != 0){
          if(prior_spikes[[i]][[j]]$y[2] > 0)graphics::arrows(x0 = prior_spikes[[i]][[j]]$x[1], y0 = prior_spikes[[i]][[j]]$y[1], y1 = prior_spikes[[i]][[j]]$y[2], lwd = 3, lty = 1, col = "grey50")
        }
      }
      if(length(res_spikes[[i]]) != 0){
        for(j in 1:length(length(res_spikes[[i]]))){
          if(res_spikes[[i]][[j]]$y[2] > 0)graphics::arrows(x0 = res_spikes[[i]][[j]]$x[1], y0 = res_spikes[[i]][[j]]$y[1], y1 = res_spikes[[i]][[j]]$y[2], lwd = 3, lty = 1)
        }
      }


      plots <- NULL

    }else if(plot_type == "ggplot"){


      temp_plot <- ggplot2::ggplot()

      # maybe later
      if(FALSE){
        if(CI)temp_plot <- temp_plot + ggplot2::geom_polygon(data = data.frame(
          xx = c(res_CI[i,1], res_CI[i,1], res_den[[i]]$x[res_den[[i]]$x >= res_CI[i,1] & res_den[[i]]$x <= res_CI[i,2]], res_CI[i,2], res_CI[i,2]),
          yy = c(0, res_den[[i]]$y[which.max(res_den[[i]]$x >= res_CI[i,1])], res_den[[i]]$y[res_den[[i]]$x >= res_CI[i,1] & res_den[[i]]$x <= res_CI[i,2]], res_den[[i]]$y[which.min(res_den[[i]]$x <= res_CI[i,2])], 0)
        ),
        ggplot2::aes_string(
          x = "xx",
          y = "yy"
        ),
        fill = "grey80"
        )
        if(mean)temp_plot   <- temp_plot + ggplot2::geom_line(data = data.frame(
          x_e = rep(mean(y[,i]),2),
          y_e = c(0,res_den[[i]]$y[which.max(res_den[[i]]$x > mean(y[,i]))])
        ),
        ggplot2::aes_string(x = "x_e", y = "y_e"), size = 1.25)
        if(median)temp_plot <- temp_plot + ggplot2::geom_line(data = data.frame(
          x_m = rep(median(y[,i]),2),
          y_m = c(0,res_den[[i]]$y[which.max(res_den[[i]]$x > median(y[,i]))])
        ),
        ggplot2::aes_string(x = "x_m", y = "y_m"), size = 1.25)
      }

      # add densities
      if(prior)if(length(prior_den[[i]]) != 0)temp_plot <- temp_plot + ggplot2::geom_line(
        data  = prior_den[[i]],
        ggplot2::aes_string(x = "x", y = "y"),
        color = "grey50", linetype = 2, size = 1.25)
      if(length(res_den[[i]]) != 0)temp_plot <- temp_plot + ggplot2::geom_line(
        data  = res_den[[i]],
        ggplot2::aes_string(x = "x", y = "y"),
        color = "black", linetype = 1, size = 1.25)

      # add spikes
      if(prior){
        if(length(prior_spikes[[i]]) != 0){
          if(prior_spikes[[i]][[j]]$y[2] > 0)temp_plot <- temp_plot + ggplot2::geom_segment(
            data = data.frame(
              x       = prior_spikes[[i]][[j]]$x[1],
              y_start = prior_spikes[[i]][[j]]$y[1],
              y_end   = prior_spikes[[i]][[j]]$y[2]),
            ggplot2::aes_string(x = "x", xend = "x", y = "y_start", yend = "y_end"),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm")),
            color = "grey50", size = 1.25)
        }
      }
      if(length(res_spikes[[i]]) != 0){
        for(j in 1:length(length(res_spikes[[i]]))){
          if(res_spikes[[i]][[j]]$y[2] > 0)temp_plot <- temp_plot + ggplot2::geom_segment(
            data = data.frame(
              x       = res_spikes[[i]][[j]]$x[1],
              y_start = res_spikes[[i]][[j]]$y[1],
              y_end   = res_spikes[[i]][[j]]$y[2]),
            ggplot2::aes_string(x = "x", xend = "x", y = "y_start", yend = "y_end"),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm")),
            color = "black", size = 1.25)
        }
      }

      # add axis
      if(any_density & any_spikes){
        temp_plot <- temp_plot + ggplot2::scale_y_continuous(
          name     = "Density",
          limits = range(pretty(temp_ylim)),
          breaks = pretty(temp_ylim),
          labels = pretty(temp_ylim),
          sec.axis = ggplot2::sec_axis(
            ~ .,
            name   = "Probability",
            breaks = seq(0, temp_ylim[2], length.out = 5),
            labels = c("0", ".25", ".50", ".75", "1"))
        )
      }else if(any_density){
        temp_plot <- temp_plot + ggplot2::scale_y_continuous(
          name   = "Density",
          limits = range(pretty(temp_ylim)),
          breaks = pretty(temp_ylim),
          labels = pretty(temp_ylim))
      }else if(any_spikes){
        temp_plot <- temp_plot + ggplot2::scale_y_continuous(
          name   = "Probability",
          limits = temp_ylim,
          breaks = seq(0, temp_ylim[2], length.out = 5),
          labels = c("0", ".25", ".50", ".75", "1"))
      }

      temp_plot <- temp_plot + ggplot2::scale_x_continuous(
        name   = par_names[[i]],
        limits = range(pretty(temp_xlim)),
        breaks = pretty(temp_xlim),
        labels = pretty(temp_xlim))

      plots <- c(plots, list(temp_plot))
    }
  }


  if(length(plots) == 1)plots <- plots[[1]]

  return(plots)
}
.plot.RoBMA_biv     <- function(fit, par, plot_type, ...){

  # get the parameter values & density
  y1     <- fit$RoBMA$samples[["averaged"]][[par[1]]]
  y2     <- fit$RoBMA$samples[["averaged"]][[par[2]]]


  par_names <- list()
  for(p in 1:length(par)){
    if(par[p] == "mu"){
      if(fit$add_info$effect_size == "r")par_names[[p]] <- bquote(~rho~"(averaged)")
      if(fit$add_info$effect_size == "d")par_names[[p]] <- bquote("Cohen's"~italic(d)~"(averaged)")
      if(fit$add_info$effect_size == "y")par_names[[p]] <- bquote(~mu~"(averaged)")
    }else if(par[p] == "tau"){
      if(fit$add_info$effect_size == "r"){
        if(fit$add_info$mu_transform == "cohens_d"){
          par_names[[p]] <- bquote(~tau~"(Cohen's"~italic(d)~ "scale; averaged)")
        }else if(fit$add_info$mu_transform == "fishers_z"){
          par_names[[p]] <- bquote(~tau~"(Fisher's"~italic(z)~ "scale; averaged)")
        }
      }else{
        par_names[[p]] <- bquote(~tau~"(averaged)")
      }
    }
  }



  # create appropriate plot
  if(plot_type == "base"){

    # set up margins
    if(length(list(...)) == 0){
      graphics::par(mar = c(4,4,2,2))
    }else{
      graphics::par(list(...))
    }

    graphics::plot(y1, y2, bty = "n", las = 1, pch = 16, col = scales::alpha("black", 0.4),
                   xlab = "", ylab = "", main = "")
    graphics::mtext(par_names[[1]], side = 1, line = 2.5, cex = 1.25)
    graphics::mtext(par_names[[2]], side = 2, line = 2.5, cex = 1.25)

    out <- NULL

  }else if(plot_type == "ggplot"){

    out <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = y1, y = y2), alpha = .4)
    out <- out + ggplot2::scale_x_continuous(
      name   = par_names[[1]],
      limits = range(pretty(range(y1))),
      breaks = pretty(range(y1)),
      labels = pretty(range(y1)))
    out <- out + ggplot2::scale_y_continuous(
      name   = par_names[[2]],
      limits = range(pretty(range(y2))),
      breaks = pretty(range(y2)),
      labels = pretty(range(y2)))

  }

  return(out)
}
.plot.RoBMA_weightf <- function(fit, type, plot_type, mean, median, CI, prior, ...){

  # deal with only null models
  if(all(sapply(fit$models,function(m).is_parameter_null(m$priors, "omega"))) & type == "conditional")stop("The ensemble cointains no non-null model adjusting for publication bias.")
  if(all(sapply(fit$models,function(m)m$priors$omega$distribution == "point")))stop("The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")

  all_cuts  <- .get_omega_mapping(fit$models, cuts_only = TRUE)

  # get the x-axis coordinate order
  coord_order <- sort(rep(1:(length(all_cuts)-1),2), decreasing = TRUE)
  x           <- all_cuts[rev(c(1, sort(rep(2:(length(all_cuts)-1), 2)), length(all_cuts)))]

  # get the y-axis coordinates
  if(mean){
    y_mean <- apply(fit$RoBMA$samples[[type]]$omega, 2, mean)[coord_order]
  }
  if(median){
    y_median <- apply(fit$RoBMA$samples[[type]]$omega, 2, median)[coord_order]
  }
  if(CI){
    y_lCI <- apply(fit$RoBMA$samples[[type]]$omega, 2, stats::quantile, probs = (1-CI)/2)[coord_order]
    y_uCI <- apply(fit$RoBMA$samples[[type]]$omega, 2, stats::quantile, probs = 1-(1-CI)/2)[coord_order]
  }


  if(prior){
    prior_samples <- .plot.prior_data_joint_omega(fit$models[fit$add_info$converged], type)
    if(is.null(prior_samples))prior_samples <- matrix(1, ncol = 1, nrow = 1)

    if(mean){
      prior_mean <- apply(prior_samples, 2, mean)[coord_order]
    }
    if(median){
      prior_median <- apply(prior_samples, 2, median)[coord_order]
    }
    if(CI){
      prior_lCI <- apply(prior_samples, 2, stats::quantile, probs = (1-CI)/2)[coord_order]
      prior_uCI <- apply(prior_samples, 2, stats::quantile, probs = 1-(1-CI)/2)[coord_order]
    }
  }


  # axis labels
  x_labes <- trimws(rev(all_cuts), which = "both", whitespace = "0")
  y_labes <- trimws(seq(0,1,.1),   which = "both", whitespace = "0")
  x_labes[length(x_labes)] <- 0 # fix the omitted 0
  y_labes[1] <- 0

  #main_text <- paste0("Estimated ", ifelse(type == "conditional", "conditional ","averaged "), "Weights")
  x_text <- bquote(italic(p)*.(paste0("-value (",ifelse(any(sapply(fit$models, function(m)m$priors$omega$distribution) == "one.sided"),"one-sided", "two-sided"),")")))
  y_text <- bquote("Publication prob. ("*omega*";"~.(ifelse(type == "conditional", "conditional)","averaged)")))


  # create appropriate plot
  if(plot_type == "base"){

    # set up margins
    if(length(list(...)) == 0){
      graphics::par(mar = c(4,4,2,2))
    }else{
      graphics::par(list(...))
    }

    graphics::plot(NA, type = "n", xlim = c(0, 1), ylim = c(0,1), xaxt = "n", yaxt = "n", bty = "n",
                   xlab = "", ylab = "", main = "")
    graphics::axis(1, at = rev(all_cuts), labels = x_labes)
    graphics::axis(2, at = seq(0,1,.1),   labels = y_labes, las = 1)
    graphics::mtext(x_text, side = 1, line = 2.5, cex = 1.25)
    graphics::mtext(y_text, side = 2, line = 2.5, cex = 1.25)

    if(CI)graphics::polygon(
      x = c(x, rev(x)),
      y = c(y_lCI, rev(y_uCI)),
      col = "grey80", border = NA
    )
    if(prior){
      if(CI){
        graphics::lines(x, prior_lCI, lwd = 1, lty = 2, col = "grey50")
        graphics::lines(x, prior_uCI, lwd = 1, lty = 2, col = "grey50")
      }
    }

    if(mean)graphics::lines(x, y_mean,     lwd = 2)
    if(median)graphics::lines(x, y_median, lwd = 2)
    if(prior){
      if(mean)graphics::lines(x, prior_mean,     lwd = 2, lty = 2, col = "grey50")
      if(median)graphics::lines(x, prior_median, lwd = 2, lty = 2, col = "grey50")
    }

    out <- NULL

  }else if(plot_type == "ggplot"){

    out <- ggplot2::ggplot()
    if(CI)out     <- out + ggplot2::geom_polygon(
      ggplot2::aes(
        x  = c(x, rev(x)),
        y  = c(round(y_lCI,5), rev(round(y_uCI,5)))), # rounding is reuqired for nummerical imprecission issues
      fill = "grey80"
    )
    if(prior){
      if(CI)out <- out + ggplot2::geom_path(ggplot2::aes(x = x, y = prior_lCI), size = 1, linetype = 2, color = "grey50")
      if(CI)out <- out + ggplot2::geom_path(ggplot2::aes(x = x, y = prior_uCI), size = 1, linetype = 2, color = "grey50")
    }

    if(mean)out   <- out + ggplot2::geom_path(ggplot2::aes(x = x, y = y_mean),   size = 1.25)
    if(median)out <- out + ggplot2::geom_path(ggplot2::aes(x = x, y = y_median), size = 1.25)
    if(prior){
      if(mean)out   <- out + ggplot2::geom_path(ggplot2::aes(x = x, y = prior_mean),   size = 1.25, linetype = 2, color = "grey50")
      if(median)out <- out + ggplot2::geom_path(ggplot2::aes(x = x, y = prior_median), size = 1.25, linetype = 2, color = "grey50")
    }
    out <- out + ggplot2::scale_x_continuous(x_text, breaks = rev(all_cuts), labels = x_labes, limits = c(0,1))
    out <- out + ggplot2::scale_y_continuous(y_text, breaks = seq(0,1,.1),   labels = y_labes, limits = c(0,1))

  }

  return(out)

}
.plot.RoBMA_theta   <- function(fit, par, type, plot_type, mean, median, CI, order, digits_estimates, ...){

  # CIs must be set
  if(is.null(CI))CI <- .95
  if(!median)mean   <- TRUE
  if(length(par) == 2){
    mean   <- TRUE
    median <- FALSE
  }else if(par == "forest"){
    mean   <- TRUE
    median <- FALSE
  }

  # get the parameter values
  mu     <- fit$RoBMA$samples[[type]][["mu"]]
  if(length(mu) == 0)stop("The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameters or all of the models did not converge.")
  mu_est <- if(median) stats::median(mu) else if(mean) base::mean(mu)
  mu_CI  <- unname(stats::quantile(mu, probs = c((1-CI)/2, 1-(1-CI)/2)))

  ## obtain the plotting information depending on the plot settings
  # get the original studies mean and CI
  if(any(par == "forest")){
    results_orig <- .get_effect_and_ci(fit$add_info, .95)
  }
  # get the estimated theta
  if(any(par == "theta")){
    theta_est   <- fit$RoBMA$samples[[type]][["theta"]]
    results_est <- cbind.data.frame(
      lCI  = apply(theta_est, 2, stats::quantile, prob = (1-CI)/2),
      est  = apply(theta_est, 2, if(median) stats::median else if(mean) base::mean),
      uCI  = apply(theta_est, 2, stats::quantile, prob = 1-(1-CI)/2),
      name = if(!is.null(fit$add_info$study_names)) fit$add_info$study_names else paste0("Study ", 1:ncol(theta_est))
    )
  }

  if(length(par) == 1){
    if(par == "forest"){
      results <- results_orig
    }else if(par == "theta"){
      results <- results_est
    }
    if(!is.null(order)){
      if(order == "ascending"){
        results <- results[order(results$est, decreasing = TRUE),]
      }else if(order == "descending"){
        results <- results[order(results$est, decreasing = FALSE),]
      }else if(order == "alphabetical"){
        results <- results[order(results$name, decreasing = TRUE),]
      }
    }
  }else{
    if(!is.null(order)){
      if(order == "ascending"){
        results_est  <- results_est[order(results_est$est,  decreasing = TRUE),]
        results_orig <- results_orig[order(results_est$est, decreasing = TRUE),]
      }else if(order == "descending"){
        results_est  <- results_est[order(results_est$est,  decreasing = FALSE),]
        results_orig <- results_orig[order(results_est$est, decreasing = FALSE),]
      }else if(order == "alphabetical"){
        results_est  <- results_est[order(results_est$name,  decreasing = TRUE),]
        results_orig <- results_orig[order(results_est$name, decreasing = TRUE),]
      }
    }
  }


  # general plotting info
  if(length(par) == 2){
    y_at <- seq(.025, 1, length.out = nrow(results_orig)*3 + 6)
    y_at <- y_at[-c(1:3 ,(length(y_at) - 2):length(y_at))]
    y_at <- y_at - (y_at[1] - y_at[2])/2
    y_at <- y_at[1:length(y_at) %% 3 != 1]
  }else{
    y_at <- seq(.025, 1, length.out = nrow(results) + 4)
    y_at <- y_at[-c(1, 2 ,length(y_at) - 1, length(y_at))]
  }


  # estimate information
  est_info <- paste0(ifelse(median,"Median", "Mean"),
                     if(CI)paste0(" [", round(CI*100), "% CI]"))

  if(fit$add_info$effect_size == "r")x_name <- bquote(~rho~.(paste0("(",type,")")))
  if(fit$add_info$effect_size == "d")x_name <- bquote("Cohen's"~italic(d)~.(paste0("(",type,")")))
  if(fit$add_info$effect_size == "y")x_name <- bquote(~mu~.(paste0("(",type,")")))

  if(length(par) == 2){
    est_info <- c(est_info, sapply(1:nrow(results_orig), function(i){
      c(paste0(
        if(round(results_orig$est[i], digits_estimates) >= 0)" ",
        format(round(results_orig$est[i], digits_estimates), nsmall = digits_estimates),
        if(CI)paste0(" [",format(round(results_orig$lCI[i], digits_estimates), nsmall = digits_estimates),", ",
                     format(round(results_orig$uCI[i], digits_estimates), nsmall = digits_estimates),"]")
      ),
      paste0(
        if(round(results_est$est[i], digits_estimates) >= 0)" ",
        format(round(results_est$est[i], digits_estimates), nsmall = digits_estimates),
        if(CI)paste0(" [",format(round(results_est$lCI[i], digits_estimates), nsmall = digits_estimates),", ",
                     format(round(results_est$uCI[i], digits_estimates), nsmall = digits_estimates),"]")
      )
      )

    }))
    est_info <- c(est_info, paste0(
      if(round(mu_est, digits_estimates) >= 0)" ",
      format(round(mu_est, digits_estimates), nsmall = digits_estimates),
      if(CI)paste0(" [",format(round(mu_CI[1], digits_estimates), nsmall = digits_estimates),", ",
                   format(round(mu_CI[2], digits_estimates), nsmall = digits_estimates),"]")))
  }else{
    est_info <- c(est_info, sapply(1:nrow(results), function(i){
      paste0(
        if(round(results$est[i], digits_estimates) >= 0)" ",
        format(round(results$est[i], digits_estimates), nsmall = digits_estimates),
        if(CI)paste0(" [",format(round(results$lCI[i], digits_estimates), nsmall = digits_estimates),", ",
                     format(round(results$uCI[i], digits_estimates), nsmall = digits_estimates),"]")
      )
    }))
    est_info <- c(est_info, paste0(
      if(round(mu_est, digits_estimates) >= 0)" ",
      format(round(mu_est, digits_estimates), nsmall = digits_estimates),
      if(CI)paste0(" [",format(round(mu_CI[1], digits_estimates), nsmall = digits_estimates),", ",
                   format(round(mu_CI[2], digits_estimates), nsmall = digits_estimates),"]")))
  }


  if(length(par) == 2){
    x_range <- range(c(results_orig$lCI, results_est$lCI, results_orig$uCI, results_est$uCI,0))
  }else{
    x_range <- range(c(results$lCI, results$uCI, 0))
  }
  x_range[1] <- x_range[1] - (x_range[2] - x_range[1])*.05
  x_range[2] <- x_range[2] + (x_range[2] - x_range[1])*.05

  # create appropriate plot
  if(plot_type == "base"){

    # set up margins
    if(length(list(...)) == 0){
      graphics::par(mar = c(4,max(sapply(as.character(if(length(par) == 2) results_est$name else results$name), nchar))*12/25,
                            2,max(sapply(as.character(est_info), nchar))*10/25))
    }else{
      graphics::par(list(...))
    }


    graphics::plot(NA, bty = "n", las = 1, xlab = "", ylab = "", main = "", yaxt = "n", ylim = c(0,1), xlim = x_range)
    if(length(par) == 2){
      graphics::axis(2,
                     at     = c(0.025, sapply(1:(length(y_at)/2), function(i)mean(y_at[(2*i-1):(2*i)])), 1),
                     labels = rev(c("Study:", as.character(results_est$name), "Overall")),
                     las = 1, col = NA)
    }else{
      graphics::axis(2, at = c(0.025, y_at, 1), labels = rev(c("Study:",as.character(results$name), "Overall")), las = 1, col = NA)
    }
    graphics::mtext(x_name, side = 1, line = 2.5, cex = 1.25)
    graphics::axis(4, at = c(0.025, y_at, 1), labels = rev(est_info), las = 1, col = NA, hadj = 0)
    graphics::lines(c(0, 0), c(0, 1), lty = 3)
    if(length(par) == 2){
      for(i in 1:nrow(results_orig)){
        graphics::lines(c(results_orig$lCI[i],results_orig$uCI[i]), rep(rev(y_at)[2*i-1], 2))
        graphics::points(results_orig$est[i], rev(y_at)[2*i-1], pch = 15)

        graphics::lines(c(results_est$lCI[i],results_est$uCI[i]), rep(rev(y_at)[2*i], 2), col = "grey50")
        graphics::points(results_est$est[i], rev(y_at)[2*i], pch = 15, col = "grey50")
      }
    }else{
      for(i in 1:nrow(results)){
        graphics::lines(c(results$lCI[i],results$uCI[i]), rep(rev(y_at)[i], 2))
        graphics::points(results$est[i], rev(y_at)[i], pch = 15)
      }
    }

    graphics::polygon(
      x = c(mu_CI[1], mu_est, mu_CI[2], mu_est),
      y = c(.025, .01, .025, .04),
      col = "black"
    )

    out <- NULL

  }else if(plot_type == "ggplot"){

    out <- ggplot2::ggplot()

    if(length(par) == 2){
      out <- out + ggplot2::geom_errorbarh(ggplot2::aes(
        xmin = as.vector(sapply(1:(length(y_at)/2), function(i)c(results_orig$lCI[i], results_est$lCI[i]))),
        xmax = as.vector(sapply(1:(length(y_at)/2), function(i)c(results_orig$uCI[i], results_est$uCI[i]))),
        y    = rev(y_at)),
        color = rep(c("black","grey50"), length(y_at)/2))
      out   <- out + ggplot2::geom_point(ggplot2::aes(
        x = as.vector(sapply(1:(length(y_at)/2), function(i)c(results_orig$est[i], results_est$est[i]))),
        y = rev(y_at)),
        shape = 15,
        color = rep(c("black","grey50"), length(y_at)/2))
    }else{
      out <- out + ggplot2::geom_errorbarh(ggplot2::aes(
        xmin = results$lCI,
        xmax = results$uCI,
        y    = rev(y_at)))
      out   <- out + ggplot2::geom_point(ggplot2::aes(x = results$est, y = rev(y_at)), shape = 15)
    }

    out <- out + ggplot2::geom_polygon(
      ggplot2::aes(
        x = c(mu_CI[1], mu_est, mu_CI[2], mu_est),
        y = c(.025, .01, .025, .04)),
      fill = "black",
    )
    if(length(par) == 2){
      out <- out + ggplot2::scale_y_continuous("",
                                               breaks   = c(0.025, sapply(1:(length(y_at)/2), function(i)mean(y_at[(2*i-1):(2*i)])), 1),
                                               labels   = rev(c("Study:", as.character(results_est$name), "Overall")),
                                               limits   = c(0, 1),
                                               sec.axis = ggplot2::sec_axis(
                                                 ~ .,
                                                 breaks = c(0.025, y_at, 1),
                                                 labels = rev(est_info))
      )
    }else{
      out <- out + ggplot2::scale_y_continuous("",
                                               breaks   = c(0.025, y_at, 1),
                                               labels   = rev(c("Study:",as.character(results$name), "Overall")),
                                               limits   = c(0, 1),
                                               sec.axis = ggplot2::sec_axis(
                                                 ~ .,
                                                 breaks = c(0.025, y_at, 1),
                                                 labels = rev(est_info))
      )
    }

    out <- out + ggplot2::geom_line(ggplot2::aes(x = c(0,0), y = c(0,1)), linetype = "dotted")
    out <- out + ggplot2::scale_x_continuous(
      name   = x_name,
      limits = range(pretty(x_range)),
      breaks = pretty(x_range),
      labels = pretty(x_range))
    out <- out + ggplot2::theme(axis.title.y      = ggplot2::element_blank(),
                                axis.line.y       = ggplot2::element_blank(),
                                axis.ticks.y      = ggplot2::element_blank(),
                                axis.text.y       = ggplot2::element_text(
                                  hjust = 0,
                                  color = "black"),
                                axis.text.y.right = ggplot2::element_text(
                                  hjust = 1,
                                  color = if(length(par) == 2)c("black",rep(c("grey50","black"), length(y_at)/2),"black")))
  }


  return(out)
}
.plot.RoBMA_ind     <- function(fit, par, type, plot_type, mean, median, CI, order, digits_estimates, show_figures, ...){

  # only 95% CI is supported and only one estimates is set
  CI        <- .95
  omega_ind <- .get_omega_mapping(fit$models)
  if(!median) mean  <- TRUE else mean <- FALSE

  # select the type of individual plots
  type <- type[type != "individual"]
  if(length(type) == 0)type <- "averaged"

  # get the parameter values
  summary_info <- summary(fit, probs = c((1-CI)/2, 1-(1-CI)/2), conditional = type == "conditional")
  sum_all <- summary_info[[type]]
  sum_ind <- summary(fit, type = "individual", probs = c((1-CI)/2, 1-(1-CI)/2))[["overview"]]

  # the overall estimates
  par_est <- sum_all[grepl(par, rownames(sum_all)) ,c("Mean", "Median")[c(mean, median)]]
  par_CI  <- sum_all[grepl(par, rownames(sum_all)) ,3:4]

  # individual model estimates
  if(par == "omega"){
    if(is.null(omega_ind))stop("The ensemble cointains no non-null model adjusting for publication bias.")
    mod_est   <- vector(mode = "list", length = ncol(do.call(rbind,omega_ind)))
    mod_CI    <- vector(mode = "list", length = ncol(do.call(rbind,omega_ind)))
  }else{
    mod_est   <- NULL
    mod_CI    <- NULL
  }
  mod_names   <- list()
  mod_marglik <- NULL
  mod_prior   <- NULL

  if(type == "conditional"){
    model_ind <- c(1:length(sum_ind))[sapply(fit$models, function(m)!.is_parameter_null(m$priors, par))]
  }else if(type == "averaged"){
    model_ind <- c(1:length(sum_ind))
  }
  if(length(model_ind) == 0)stop("The ensemble contains no non-null model with the specified parameter.")

  for(i in model_ind){

    if(fit$models[[i]]$priors[[par]]$distribution != "point"){

      # add the estimates values
      if(par == "omega"){

        for(j in 1:length(omega_ind[[i]])){
          mod_est[[j]] <- c(mod_est[[j]], sum_ind[[i]]$tab[rev(rownames(sum_ind[[i]]$tab)[grepl("omega",rownames(sum_ind[[i]]$tab))])[omega_ind[[i]][j]], c("Mean", "Median")[c(mean, median)]])
          mod_CI[[j]]  <- rbind(mod_CI[[j]], unlist(sum_ind[[i]]$tab[rev(rownames(sum_ind[[i]]$tab)[grepl("omega",rownames(sum_ind[[i]]$tab))])[omega_ind[[i]][j]], c(3, 5)]))
        }

      }else{

        mod_est <- c(mod_est,    sum_ind[[i]]$tab[rownames(sum_ind[[i]]$tab) == par, c("Mean", "Median")[c(mean, median)]])
        mod_CI  <- rbind(mod_CI, unlist(sum_ind[[i]]$tab[rownames(sum_ind[[i]]$tab) == par, c(3, 5)]))

      }

    }else{

      # add the estimates null values
      if(par == "omega"){
        for(j in 1:ncol(do.call(rbind,omega_ind))){
          mod_est[[j]] <- c(mod_est[[j]],    1)
          mod_CI[[j]]  <- rbind(mod_CI[[j]], c(1,1))
        }
      }else{
        mod_est <- c(mod_est,        fit$models[[i]]$priors[[par]]$parameters$location)
        mod_CI  <- rbind(mod_CI, rep(fit$models[[i]]$priors[[par]]$parameters$location, 2))
      }

    }

    # marglik information
    mod_marglik    <- c(mod_marglik, fit$models[[i]]$marg_lik$logml)
    mod_prior      <- c(mod_prior,   fit$models[[i]]$prior_odds)

    # prior names
    mod_names <- c(mod_names, list(
      list(
        mu    = print(fit$models[[i]]$priors$mu,    silent = TRUE, plot = TRUE),
        tau   = print(fit$models[[i]]$priors$tau,   silent = TRUE, plot = TRUE),
        omega = print(fit$models[[i]]$priors$omega, silent = TRUE, plot = TRUE))
    ))
    #}

  }

  # compute posterior prob
  if(all(mod_prior == 0) | length(mod_prior) == 0)stop("The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameters or all of the models did not converge.")
  mod_prior <- mod_prior/sum(mod_prior)
  mod_post  <- bridgesampling::post_prob(mod_marglik, prior_prob = mod_prior)

  # merge the model results together
  if(par == "omega"){
    mod_res   <- vector(mode = "list", length = ncol(do.call(rbind,omega_ind)))
    for(i in 1:length(mod_res)){
      mod_res[[i]] <- data.frame(
        est     = mod_est[[i]],
        lCI     = mod_CI[[i]][,1],
        uCI     = mod_CI[[i]][,2]
      )
    }
  }else{
    mod_res <- list(data.frame(
      est     = mod_est,
      lCI     = mod_CI[,1],
      uCI     = mod_CI[,2]
    ))
  }

  # order the values
  for(i in 1:length(mod_res)){
    if(!is.null(order)){
      if(any(order == "marglik")){
        if(any(order == "ascending")){
          mod_res[[i]] <- mod_res[[i]][order(mod_marglik, decreasing = TRUE),]
        }else{
          mod_res[[i]] <- mod_res[[i]][order(mod_marglik, decreasing = FALSE),]
        }
      }else if(any(order == "prob")){
        if(any(order == "ascending")){
          mod_res[[i]] <- mod_res[[i]][order(mod_post, decreasing = TRUE),]
        }else{
          mod_res[[i]] <- mod_res[[i]][order(mod_post, decreasing = FALSE),]
        }
      }else{
        if(any(order == "descending")){
          mod_res[[i]] <- mod_res[[i]][nrow(mod_res[[i]]):1,]
        }
      }
    }
  }
  mod_marglik <- mod_marglik[as.numeric(rownames(mod_res[[1]]))]
  mod_post    <- mod_post[as.numeric(rownames(mod_res[[1]]))]


  # general plotting info
  y_at <- seq(.05, 1, length.out = nrow(mod_res[[1]]) + 2)
  y_at <- y_at[-c(1, length(y_at))]
  est_name  <- paste0("Overall (", if(type == "averaged") "Model-Averaged" else if(type == "conditional") "Conditional", ")")
  mod_names <- unlist(sapply(as.numeric(rownames(mod_res[[1]])), function(i){
    c(bquote(~mu~"~"~.(mod_names[[i]]$mu)),
      bquote(~tau~"~"~.(mod_names[[i]]$tau)),
      bquote(~omega~"~"~.(mod_names[[i]]$omega)))
  }))

  # x axis names
  if(par == "mu"){
    if(fit$add_info$effect_size == "r")par_names <- list(bquote(~rho))
    if(fit$add_info$effect_size == "d")par_names <- list(bquote("Cohen's"~italic(d)))
    if(fit$add_info$effect_size == "y")par_names <- list(bquote(~mu))
  }else if(par == "tau"){
    if(fit$add_info$effect_size == "r"){
      if(fit$add_info$mu_transform == "cohens_d"){
        par_names <- list(bquote(~tau~"(Cohen's"~italic(d)~ "scale)"))
      }else if(fit$add_info$mu_transform == "fishers_z"){
        par_names <- list(bquote(~tau~"(Fisher's"~italic(z)~ "scale"))
      }
    }else{
      par_names <- list(bquote(~tau))
    }
  }else if(par == "omega"){
    par_names <- sapply(rownames(sum_all)[grepl(par, rownames(sum_all))], function(x){
      bquote(~omega[~.(substr(x,6,nchar(x)))]~.(paste0("(",summary_info$add_info$weight_type, ")")))
    })
  }


  # estimate information
  est_info <- vector(mode = "list", length = length(mod_res))
  for(i in 1:length(mod_res)){
    est_info[[i]] <- paste0(ifelse(median,"Median", "Mean"),
                            paste0(" [", round(CI*100), "% CI]\nPost. prob. (Prior prob.)"))
    est_info[[i]] <- c(est_info[[i]], sapply(1:nrow(mod_res[[1]]), function(j){
      paste0(
        format(round(mod_res[[i]]$est[j], digits_estimates), nsmall = digits_estimates),
        paste0(" [",format(round(mod_res[[i]]$lCI[j], digits_estimates), nsmall = digits_estimates),", ",
               format(round(mod_res[[i]]$uCI[j], digits_estimates), nsmall = digits_estimates),"]\n",
               format(round(mod_post[j], digits_estimates), nsmall = digits_estimates),
               " (",format(round(mod_prior[j], digits_estimates), nsmall = digits_estimates),")")
      )
    }))
    est_info[[i]] <- c(est_info[[i]], paste0(
      format(round(par_est[i], digits_estimates), nsmall = digits_estimates),
      if(CI)paste0(" [",format(round(par_CI[i,1], digits_estimates), nsmall = digits_estimates),", ",
                   format(round(par_CI[i,2], digits_estimates), nsmall = digits_estimates),"]")))
  }


  # create appropriate plot
  plots <- list()
  if(is.null(show_figures) | length(mod_res) == 1){
    plots_ind <- c(1:length(mod_res))
  }else{
    plots_ind <- c(1:length(mod_res))[show_figures]
  }
  # a message with info about muliple plots
  if(plot_type == "base" & length(plots_ind) > 1)cat(paste0(length(plots_ind), " plots will be produced. See '?layout' for help with setting multiple plots."))

  if(par == "omega"){
    x_range <- c(0, 1)
  }else{
    x_range <- c(min(mod_res[[i]]$lCI), max(mod_res[[i]]$uCI))
    x_range[1] <- x_range[1] - (x_range[2] - x_range[1])*.05
    x_range[2] <- x_range[2] + (x_range[2] - x_range[1])*.05

    if(all(x_range == 0)){
      if(par == "mu")x_range <- c(-.5, .5)
      if(par == "tau")x_range <- c(0, .5)
    }
  }

  for(i in plots_ind){
    if(plot_type == "base"){

      # set up margins
      if(length(list(...)) == 0){
        graphics::par(mar = c(4,max(sapply(as.character(mod_names), nchar))*1/5,
                              0,max(sapply(as.character(est_info[[i]]), nchar))*2/7))
      }else{
        graphics::par(list(...))
      }

      graphics::plot(NA, bty = "n", las = 1, xlab = "", ylab = "", main = "", yaxt = "n", ylim = c(0,1.05), xlim = x_range)

      y_at_lab <- sort(rep(y_at, 3))
      y_at_lab <- y_at_lab + c(-1,0,1)*(y_at[2] - y_at[1])/4

      for(j in 1:length(y_at_lab)){
        graphics::axis(2, at = rev(y_at_lab)[j], labels = mod_names[[j]], las = 1, col = NA)
      }
      graphics::axis(2, at = c(.05, 1), labels = c(est_name, "Model:"), las = 1, col = NA)
      graphics::mtext(par_names[[i]], side = 1, line = 2.5, cex = 1.25)
      graphics::axis(4, at = c(0.05, y_at, 1), labels = rev(est_info[[i]]), las = 1, col = NA, hadj = 0)


      if(par != "omega")graphics::lines(c(0, 0), c(0, 1), lty = 3)
      if(par == "omega")graphics::lines(c(1, 1), c(0, 1), lty = 3)

      for(j in 1:nrow(mod_res[[i]])){
        graphics::lines(c(mod_res[[i]]$lCI[j],mod_res[[i]]$uCI[j]), rep(rev(y_at)[j], 2))
        graphics::points(mod_res[[i]]$est[j], rev(y_at)[j], pch = 15, cex = .5 + mod_post[j])
      }

      graphics::polygon(
        x = c(par_CI[i, 1], par_est[i], par_CI[i, 2], par_est[i]),
        y = c(.05, .025, .05, .075),
        col = "black"
      )

      plots <- NULL

    }else if(plot_type == "ggplot"){

      temp_plot <- ggplot2::ggplot(mod_res[[i]])
      temp_plot <- temp_plot + ggplot2::geom_errorbarh(ggplot2::aes_string(
        xmin = "lCI",
        xmax = "uCI",
        y    = rev(y_at)))
      temp_plot <- temp_plot + ggplot2::geom_point(ggplot2::aes_string(
        x = "est",
        y = rev(y_at)),
        shape = 15, size = 1 + 5*mod_post)
      temp_plot <- temp_plot + ggplot2::geom_polygon(
        data = data.frame(
          xx = c(par_est[i], par_CI[i, 1],  par_est[i], par_CI[i, 2]),
          yy = c(.025, .05, .075, .05)
        ),
        ggplot2::aes_string(
          x = "xx",
          y = "yy"),
        fill = "black"
      )

      y_at_lab <- sort(rep(y_at, 3))
      y_at_lab <- y_at_lab + c(-1,0,1)*(y_at[2] - y_at[1])/4

      temp_plot <- temp_plot + ggplot2::scale_y_continuous("",
                                                           breaks   = c(0.05, y_at_lab, 1.00),
                                                           labels   = rev(c("Model:",mod_names, est_name)),
                                                           limits   = c(0, 1.00),
                                                           sec.axis = ggplot2::sec_axis(
                                                             ~ .,
                                                             breaks = c(0.05, y_at, 1.00),
                                                             labels = rev(est_info[[i]]))
      )

      if(par != "omega")temp_plot <- temp_plot + ggplot2::geom_line(
        data = data.frame(
          xl = c(0,0),
          yl = c(0,1)
        ),
        ggplot2::aes_string(x = "xl", y = "yl"),
        linetype = "dotted")
      if(par == "omega")temp_plot <- temp_plot + ggplot2::geom_line(
        data = data.frame(
          xl = c(1,1),
          yl = c(0,1)
        ),
        ggplot2::aes_string(x = "xl", y = "yl"),
        linetype = "dotted")

      temp_plot <- temp_plot + ggplot2::scale_x_continuous(
        name   = par_names[[i]],
        limits = if(par == "omega") c(-.01, 1.01) else range(pretty(x_range)),
        breaks = pretty(x_range),
        labels = pretty(x_range))
      temp_plot <- temp_plot + ggplot2::theme(axis.title.y      = ggplot2::element_blank(),
                                              axis.line.y       = ggplot2::element_blank(),
                                              axis.ticks.y      = ggplot2::element_blank(),
                                              axis.text.y       = ggplot2::element_text(
                                                hjust = 0,
                                                color = "black"),
                                              axis.text.y.right = ggplot2::element_text(hjust = 1))

      plots <- c(plots, list(temp_plot))
    }
  }

  if(length(plots) == 1)plots <- plots[[1]]

  return(plots)
}
