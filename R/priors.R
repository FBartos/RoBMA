#' @title Creates a RoBMA prior
#'
#' @description \code{prior} creates a prior distribution for fitting
#' a RoBMA model. The prior can be visualized by a \code{plot} function.
#'
#' @param distribution name of the prior distribution. The
#' possible options are
#' \describe{
#'   \item{\code{"point"}}{for a point density characterized by a
#'   \code{location} parameter.}
#'   \item{\code{"normal"}}{for a normal distribution characterized
#'   by a \code{mean} and \code{sd} parameters.}
#'   \item{\code{"cauchy"}}{for a Cauchy distribution characterized
#'   by a \code{location} and \code{scale} parameters. Internally
#'   converted into a generalized t-distribution with \code{df = 1}.}
#'   \item{\code{"t"}}{for a generalized t-distribution characterized
#'   by a \code{location}, \code{scale}, and \code{df} parameters.}
#'   \item{\code{"gamma"}}{for a gamma distribution characterized
#'   by either \code{shape} and \code{rate}, or \code{shape} and
#'   \code{scale} parameters. The later is internally converted to
#'   the \code{shape} and \code{rate} parametrization}
#'   \item{\code{"invgamma"}}{for an inverse-gamma distribution
#'   characterized by a \code{shape} and \code{scale} parameters. The
#'   JAGS part uses a 1/gamma distribution with a shape and rate
#'   parameter.}
#'   \item{\code{"two.sided"}}{for a two-sided weight function
#'   characterized by a vector \code{steps} and vector \code{alpha}
#'   parameters. The \code{alpha} parameter determines an alpha
#'   parameter of Dirichlet distribution which which cumulative sum
#'   is used for the weights omega.}
#'   \item{\code{"one.sided"}}{for a one-sided weight function
#'   characterized by a either a vector \code{steps} and vector
#'   \code{alpha} parameter, leading to a monotonic one-sided
#'   function, or by a a vector \code{steps}, vector \code{alpha1},
#'   and vector \code{alpha2} parameters leading non-monotonic
#'   one-sided weight function. The \code{alpha} / \code{alpha1} and
#'   \code{alpha2} parameters determine an alpha parameter of
#'   Dirichlet distribution which which cumulative sum is used for
#'   the weights omega.}
#'   \item{\code{"uniform"}}{for a uniform distribution defined on a
#'   range from \code{a} to \code{b}}S
#' }
#' @param parameters list of appropriate parameters for a given
#' \code{distribution}.
#' @param truncation list with two elements, \code{lower} and
#' \code{upper}, that define the lower and upper truncation of the
#' distribution. Defaults to \code{list(lower = -Inf, upper = Inf)}.
#' The lower truncation point is automatically set to 0 if it is
#' specified outside of the support of distributions defined only for
#' positive numbers.
#' @param prior_odds prior odds associated with a given distribution.
#' [RoBMA()] creates models corresponding to all combinations of prior
#' distributions for each of the model parameters (mu, tau, omega), and
#' sets the model priors odds to the product of its prior distributions.
#'
#' @examples # create a standart normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # create a half-normal standart normal prior distribution
#' p2 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1),
#' truncation = list(lower = 0, upper = Inf))
#'
#'
#' @export  prior
#' @rawNamespace S3method(print, RoBMA.prior)
#' @seealso [plot.RoBMA.prior()], \link[stats]{Normal}, \link[stats]{Cauchy},
#' \link[extraDistr]{LocationScaleT}, \link[stats]{GammaDist}, \link[extraDistr]{InvGamma}.
prior <- function(distribution, parameters, truncation = list(lower = -Inf, upper = Inf), prior_odds = 1){

  # general input check
  if(!(is.vector(distribution) & length(distribution) == 1 & is.character(distribution)))stop("Argument 'distribution' is incorectly specified. It must be a character vector of length one.")
  if(!is.list(parameters))stop("Argument 'parameters' must be a list.")
  if(!all(sapply(parameters, function(par)(is.vector(par) & is.numeric(par)))))stop("Elements of the 'parameter' argument must be numeric vectors.")
  if(!(is.list(truncation) & length(truncation) == 2))stop("Argument 'truncation' must be a list of length two.")
  if(!(is.vector(prior_odds) & is.numeric(prior_odds) & length(prior_odds) == 1))stop("Argument 'prior_odds' must be a numeric vector of length two.")
  if(is.null(names(truncation)))names(truncation) <- c("lower", "upper")
  if(truncation$lower >= truncation$upper)stop("The lower truncation point needs to be lower than the upper truncation point.")
  if(prior_odds <= 0)stop("Argument 'prior_odds' must be positive.")

  distribution <- tolower(distribution)


  # create an output object
  output <- list()

  # check whether the values are appropriate for the individual distribution, whether all parameters are included etc...
  if(distribution %in% c("norm", "normal")){

    if(length(parameters) != 2)stop("Normal prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("mean", "sd")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("mean", "sd")], sep = ", ", collapse = ""), " are not supported for a normal distribution."))
    }else{
      names(parameters) <- c("mean", "sd")
    }
    if(parameters$sd <= 0)stop("Parameter 'sd' must be positive.")

    # add the values to the output
    output$distribution <- "normal"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("t", "student", "student-t", "student t")){

    if(length(parameters) != 3)stop("Student-t prior distribution requires 3 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("df", "location", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("df", "location", "scale")], sep = ", ", collapse = ""), " are not supported for a student-t distribution."))
    }else{
      if(is.null(names(parameters)))names(parameters) <- c("location", "scale", "df")
    }
    if(parameters$df < 1)stop("Parameter 'df' must be positive.")
    if(!parameters$scale > 0)stop("Parameter 'scale' must be positive.")

    # add the values to the output
    output$distribution <- "t"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("cauchy")){

    if(length(parameters) != 2)stop("Cauchy prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("location", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("location", "scale")], sep = ", ", collapse = ""), " are not supported for a Cauchy distribution."))
    }else{
      if(is.null(names(parameters)))names(parameters) <- c("location", "scale")
    }
    if(!parameters$scale > 0)stop("Parameter 'scale' must be positive.")

    # add the values to the output
    output$distribution <- "t"
    output$parameters   <- list(df = 1, location = parameters$location, scale = parameters$scale) # pass as t-distribution
    output$truncation   <- truncation

  }else if(distribution %in% c("invgamma", "inversegamma", "inverse-gamma", "inverse gamma")){

    if(length(parameters) != 2)stop("Inverse gamma prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("shape", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("shape", "scale")], sep = ", ", collapse = ""), " are not supported for an inverse gamma distribution."))
    }else{
      names(parameters) <- c("shape", "scale")
    }
    if(!parameters$shape > 0)stop("Parameter 'shape' must be positive.")
    if(!parameters$scale > 0)stop("Parameter 'scale' must be positive.")

    if(truncation$lower == -Inf)truncation$lower <- 0 # change the defaul lower truncation
    if(truncation$lower < 0)stop("Lower truncation point must be non-negative.")
    # add the values to the output
    output$distribution <- "invgamma"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("gamma")){

    if(length(parameters) != 2)stop("Gamma prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("shape", "rate", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("shape", "rate", "scale")], sep = ", ", collapse = ""), " are not supported for a gamma distribution."))
    }else{
      names(parameters) <- c("shape", "rate")
    }
    if(!is.null(parameters$scale)){
      parameters$rate  <- 1/parameters$scale
      parameters$scale <- NULL
    }
    if(!parameters$shape > 0)stop("Parameter 'shape' must be positive.")
    if(!parameters$rate > 0)stop("Parameter 'rate' must be positive.")

    if(truncation$lower == -Inf)truncation$lower <- 0 # change the defaul lower truncation
    if(truncation$lower < 0)stop("Lower truncation point must be non-negative.")
    # add the values to the output
    output$distribution <- "gamma"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("point", "spike")){

    if(length(parameters) != 1)stop("Point prior distribution requires 1 parameter.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("location")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("location")], sep = ", ", collapse = ""), " are not supported for an inverse gamma distribution."))
    }else{
      names(parameters) <- c("location")
    }

    # add the values to the output
    output$distribution <- "point"
    output$parameters   <- parameters

  }else if(distribution %in% c("two.sided", "twosided", "two-sided", "two sided")){

    if(length(parameters) != 2)stop("Two-sided step function requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("alpha", "steps")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("alpha", "steps")], sep = ", ", collapse = ""), " are not supported for a two-sided weight function."))
    }else{
      names(parameters) <- c("steps", "alpha")
    }
    if(!all(parameters$alpha > 0))stop("Parameters 'alpha' must be positive.")
    if(!all(parameters$steps < 1 & parameters$steps > 0))stop("Parameters 'steps' must be higer than 0 and lower than 1.")
    if(!(all(parameters$steps == cummax(parameters$steps))))stop("Parameters 'steps' must be monotonically increasing.")
    if(length(parameters$steps) != length(parameters$alpha) - 1)stop("The parameter alpha needs to have one more argument then there are steps.")

    # reverse the ordering of alpha and weights - to corespond to ordering on t-statistics
    parameters$steps <- rev(parameters$steps)
    parameters$alpha <- rev(parameters$alpha)

    # add the values to the output
    output$distribution <- "two.sided"
    output$parameters   <- parameters

  }else if(distribution %in% c("one.sided", "onesided", "one-sided", "one sided")){

    if(!length(parameters) %in% c(2,3))stop("One-sided step function requires 2 or 3 parameters.")
    if(length(parameters) == 2){

      if(!is.null(names(parameters))){
        if(!all(names(parameters) %in% c("alpha", "steps")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("alpha", "steps")], sep = ", ", collapse = ""), " are not supported for a one-sided monotonic weight function."))
      }else{
        names(parameters) <- c("steps", "alpha")
      }
      if(!all(parameters$alpha > 0))stop("Parameters 'alpha' must be positive.")
      if(!all(parameters$steps < 1 & parameters$steps > 0))stop("Parameters 'steps' must be higer than 0 and lower than 1.")
      if(!(all(parameters$steps == cummax(parameters$steps))))stop("Parameters 'steps' must be monotonically increasing.")
      if(length(parameters$steps) != length(parameters$alpha) - 1)stop("The parameter alpha needs to have one more argument then there are steps.")

      # reverse the ordering of alpha and weights - to corespond to ordering on t-statistics
      parameters$steps <- rev(parameters$steps)
      parameters$alpha <- rev(parameters$alpha)

    }else if(length(parameters) == 3){

      if(!is.null(names(parameters))){
        if(!all(names(parameters) %in% c("alpha1", "alpha2", "steps")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("alpha1", "alpha2", "steps")], sep = ", ", collapse = ""), " are not supported for a one-sided non-monotonic weight function."))
      }else{
        names(parameters) <- c("steps","alpha1", "alpha2")
      }
      if(!all(parameters$alpha1 > 0))stop("Parameters 'alpha1' must be positive.")
      if(!all(parameters$alpha2 > 0))stop("Parameters 'alpha2' must be positive.")
      if(!all(parameters$steps < 1 & parameters$steps > 0))stop("Parameters 'steps' must be higer than 0 and lower than 1.")
      if(!(all(parameters$steps == cummax(parameters$steps))))stop("Parameters 'steps' must be monotonically increasing.")
      if(sum(parameters$steps <= .5) != length(parameters$alpha1) - 1)stop("The parameter alpha1 needs to have one more argument then there are steps <= .5.")
      if(sum(parameters$steps > .5) != length(parameters$alpha2) - 1)stop("The parameter alpha2 needs to have one more argument then there are steps > .5.")

      # reverse the ordering of alpha and weights - to corespond to ordering on t-statistics
      parameters$steps  <- rev(parameters$steps)
      parameters$alpha1 <- rev(parameters$alpha1)
    }

    # add the values to the output
    output$distribution <- "one.sided"
    output$parameters   <- parameters

  }else if(distribution %in% c("uniform", "unif")){

    if(length(parameters) != 2)stop("Uniform prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("a", "b")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("a", "b")], sep = ", ", collapse = ""), " are not supported for a normal distribution."))
    }else{
      names(parameters) <- c("a", "b")
    }
    if(parameters$a >= parameters$b)stop("Parameter 'a' must be lower than parameter 'b'.")
    if(!is.infinite(truncation$lower))if(truncation$lower != parameters$a)stop("Lower truncation must correspond to the parameter 'a'.")
    if(!is.infinite(truncation$upper))if(truncation$upper != parameters$b)stop("Upper truncation must correspond to the parameter 'b'.")

    # add the values to the output
    output$distribution <- "uniform"
    output$parameters   <- parameters
    output$truncation   <- list(lower = parameters$a, upper = parameters$b)

  }

  output$prior_odds <- prior_odds
  class(output) <- "RoBMA.prior"

  return(output)
}


#' @title Prints a RoBMA.prior object
#'
#' @param x a RoBMA prior
#' @param ... additional arguments
#' \describe{
#'   \item{silent}{to silently return the print message.}
#'   \item{plot}{to return \link[base]{bquote} formatted
#'   prior name for plotting.}
#'  }
#' @export  print.RoBMA.prior
#' @rawNamespace S3method(print, RoBMA.prior)
#' @seealso [prior()]
print.RoBMA.prior <- function(x, ...){

  dots <- list(...)
  silent <- if(is.null(dots$silent)) FALSE else as.logical(dots$silent)
  plot   <- if(is.null(dots$plot))   FALSE else as.logical(dots$plot)
  if(plot)silent <- TRUE

  name <- switch(x$distribution,
                 "normal"    = "Normal",
                 "t"         = "gen. Student-t",
                 "gamma"     = "Gamma",
                 "invgamma"  = "InvGamma",
                 "point"     = "Spike",
                 "two.sided" = "Two-sided",
                 "one.sided" = "One-sided",
                 "uniform"   = "Uniform")


  if(x$distribution == "normal"){
    parameters <- c(x$parameters$mean, x$parameters$sd)
  }else if(x$distribution == "t"){
    if(x$parameters$df == 1){
      name <- "Cauchy"
      parameters <- c(x$parameters$location, x$parameters$scale)
    }else{
      parameters <- c(x$parameters$location, x$parameters$scale, x$parameters$df)
    }
  }else if(x$distribution == "gamma"){
    parameters <- c(x$parameters$shape, x$parameters$rate)
  }else if(x$distribution == "invgamma"){
    parameters <- c(x$parameters$shape, x$parameters$scale)
  }else if(x$distribution == "point"){
    parameters <- c(x$parameters$location)
  }else if(x$distribution == "two.sided"){
    parameters <- c(
      paste0("(",paste(x$parameters$steps, collapse = ", "),")"),
      paste0("(",paste(x$parameters$alpha, collapse = ", "),")")
    )
  }else if(x$distribution == "one.sided"){
    if(all(names(x$parameters) %in% c("steps", "alpha1", "alpha2"))){
      parameters <- c(
        paste0("(",paste(x$parameters$steps, collapse = ", "),")"),
        paste0("(",paste(x$parameters$alpha1, collapse = ", "),")"),
        paste0("(",paste(x$parameters$alpha2, collapse = ", "),")")
      )
    }else{
      parameters <- c(
        paste0("(",paste(x$parameters$steps, collapse = ", "),")"),
        paste0("(",paste(x$parameters$alpha, collapse = ", "),")")
      )
    }

  }else if(x$distribution == "uniform"){
    parameters <- c(x$parameters$a, x$parameters$b)
  }

  if(plot){
    if(!x$distribution %in% c("point","one.sided","two.sided","uniform")){
      output <- bquote(italic(.(name))*.(paste0("(",paste(parameters, collapse = ", "),")"))[~"["~.(x$truncation$lower)*", "*.(x$truncation$upper)~"]"])
    }else{
      output <- bquote(italic(.(name))*.(paste0("(",paste(parameters, collapse = ", "),")")))
    }

  }else{
    output <- paste0(name,"(",paste(parameters, collapse = ", "),")")
    if(!x$distribution %in% c("point","one.sided","two.sided","uniform")){
      output <- paste0(output, "[", x$truncation$lower, ", ", x$truncation$upper, "]")
    }
  }


  if(!silent)cat(output)
  return(invisible(output))
}


#' @title Plots a RoBMA.prior object
#'
#' @param x a RoBMA prior
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot2"} for plotting. The later
#' requires \pkg{ggplot2} package to be installed.
#' @param mu_transform whether and how should the prior distribution
#' be transformed. If the prior distribution is constructed for
#' effect sizes supplied as correlations, the prior for mu paraparameter
#' is not defined on the correlation scale directly, but transformed into
#' it. Only possible if the \code{par_name = "mu"}. Defaults to
#' \code{NULL}. Other options are \code{"cohens_d"} and \code{"fishers_z"}.
#' @param weights whether the weights or weight function should
#' be returned. Only applicable for priors on the omega parameter.
#' Defaults to \code{FALSE} - the weight function is plotted.
#' @param show_figures which figures should be returned in case of
#' multiple plots are generated. Useful when priors for the omega
#' parameter are plotted and \code{weights = TRUE}.
#' @param samples how many samples should be drawn for the
#' density plot (applies only for the weight functions, other
#' prior distributions are plotted using the pdf). Defaults
#' to \code{10000}.
#' @param points how many points should be used for drawing the
#' density plot. Defaults to \code{1000}.
#' @param par_name a name of parameter to be included in the x-axis
#' label
#' @param ... additional arguments
#' @export  plot.RoBMA.prior
#' @rawNamespace S3method(plot, RoBMA.prior)
#' @seealso [prior()]
plot.RoBMA.prior <- function(x, plot_type = "base", mu_transform = NULL,
                             show_figures = -1, weights = FALSE, par_name = NULL,
                             samples = 1e6, points = 1000, ...){

  # check input
  if(plot_type == "ggplot2")plot_type <- "ggplot"
  if(!plot_type %in% c("base", "ggplot"))stop("The passed plot_type is not supported for plotting. See '?plot.RoBMA' for more details.")

  # check availability of ggplot
  if(plot_type == "ggplot"){
    if(!try(requireNamespace("ggplot2")))stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  }


  # get plotting data
  if(is.null(par_name))par_name <- ""
  if(par_name == "mu" & !is.null(mu_transform))par_name <- "rho"
  plot_data <- .plot.prior_data(x, samples, points, weights, par_name, mu_transform)


  # do the plotting
  plots <- list()
  if(is.null(show_figures) | length(plot_data$df) == 1){
    plots_ind <- c(1:length(plot_data$df))
  }else{
    plots_ind <- c(1:length(plot_data$df))[show_figures]
  }

  for(i in plots_ind){

    if(plot_type == "base"){

      # a message with info about muliple plots
      if(i == 1 & length(plot_data$df) > 1)cat(paste0(length(plot_data$df), " plots will be produced. See '?layout' for help with setting multiple plots."))

      # set up margins
      if(length(list(...)) == 0){
        graphics::par(mar = c(4, 4, 1, 1))
      }else{
        graphics::par(list(...))
      }

      graphics::plot(NA, bty = "n", las = 1, xlab = "", ylab = "", main = "", ylim = plot_data$y_range[[i]], xlim = plot_data$x_range,
           xaxt = if((x$distribution %in% c("two.sided", "one.sided") | par_name == "omega") & !weights) "n")

      if((x$distribution %in% c("two.sided", "one.sided") | par_name == "omega") & !weights){
        graphics::axis(1,
                       at     = c(0, sort(x$parameters$steps), 1),
                       labels = c(0, trimws(c(sort(x$parameters$steps), 1), which = "both", whitespace = "0")))
      }
      graphics::mtext(plot_data$x_lab[[i]], side = 1, line = 2.5, cex = 1.25)
      graphics::mtext(plot_data$y_lab[[i]], side = 2, line = 2.5, cex = 1.25)

      if(!is.null(plot_data$df[[i]]$uCI)){
        graphics::polygon(
          x = c(plot_data$df[[i]]$x,   rev(plot_data$df[[i]]$x)),
          y = c(plot_data$df[[i]]$lCI, rev(plot_data$df[[i]]$uCI)),
          col = "grey80", border = NA
        )
      }
      graphics::lines(plot_data$df[[i]]$x, plot_data$df[[i]]$y, lwd = 2)


      plots <- NULL

    }else if(plot_type == "ggplot"){

      temp_plot <- ggplot2::ggplot(plot_data$df[[i]])

      if(!is.null(plot_data$df[[i]]$uCI)){
        temp_plot <- temp_plot + ggplot2::geom_polygon(
          data = data.frame(
            xx = c(plot_data$df[[i]]$x,   rev(plot_data$df[[i]]$x)),
            yy = c(plot_data$df[[i]]$lCI, rev(plot_data$df[[i]]$uCI))
          ),
          ggplot2::aes_string(
            x = "xx",
            y = "yy"),
          fill = "grey80"
        )
      }

      temp_plot <- temp_plot + ggplot2::geom_line(ggplot2::aes_string(x = "x", y = "y"), size = 1.25)

      if((x$distribution %in% c("two.sided", "one.sided") | par_name == "omega") & !weights){
        temp_plot <- temp_plot + ggplot2::scale_x_continuous(
          name   = plot_data$x_lab[[i]],
          breaks = c(0, sort(x$parameters$steps), 1),
          labels = c(0, trimws(c(sort(x$parameters$steps), 1), which = "both", whitespace = "0")),
          limits = plot_data$x_range)
      }else{
        temp_plot <- temp_plot + ggplot2::scale_x_continuous(
          name   = plot_data$x_lab[[i]],
          limits = range(pretty(range(c(0,plot_data$x_range)))),
          breaks = pretty(range(c(0,plot_data$x_range))),
          labels = pretty(range(c(0,plot_data$x_range))))
      }

      temp_plot <- temp_plot + ggplot2::scale_y_continuous(
        name   = plot_data$y_lab[[i]],
        limits = range(pretty(plot_data$y_range[[i]])),
        breaks = pretty(plot_data$y_range[[i]]),
        labels = pretty(plot_data$y_range[[i]]))

      plots <- c(plots, list(temp_plot))
    }
  }

  if(length(plots) == 1)plots <- plots[[1]]


  # return the plots
  if(plot_type == "base"){
    return(invisible(NULL))
  }else if(plot_type == "ggplot"){
    return(plots)
  }
}

### helper functions for generating data required for plotting
.plot.prior_data <- function(prior, samples = 1e6, points = 1000, weights = FALSE, par_name = "", mu_transform = NULL, x_range = NULL){

  df      <- list()
  names   <- list()
  prob    <- list()
  x_lab   <- list()
  y_lab   <- list()
  y_range <- list()

  x_range_passed <- !is.null(x_range)

  # prepare density truncation if needed
  if(!is.null(mu_transform)){
    if(prior$distribution == "point"){
      with_trunc <- list(
        from = -1,
        to   = 1
      )
    }else{
      if(mu_transform == "cohens_d"){
        with_trunc <- list(
          from = if(!is.infinite(prior$truncation$lower)) psych::d2r(prior$truncation$lower) else -1,
          to   = if(!is.infinite(prior$truncation$upper)) psych::d2r(prior$truncation$upper) else  1
        )
      }else if(mu_transform == "fishers_z"){
        with_trunc <- list(
          from = if(!is.infinite(prior$truncation$lower)) psych::fisherz2r(prior$truncation$lower) else -1,
          to   = if(!is.infinite(prior$truncation$upper)) psych::fisherz2r(prior$truncation$upper) else  1
        )
      }
    }
  }else{
    with_trunc <- list(
      from = prior$truncation$lower,
      to   = prior$truncation$upper
    )
  }


  if(prior$distribution == "one.sided" & all(names(prior$parameters) %in% c("alpha1", "alpha2", "steps"))){

    x_range <- c(0,1)
    if(weights){
      names   <- sapply(1:(length(prior$parameters$steps)+1), function(i){
        bquote(~omega[~.(paste0("[", c(1, prior$parameters$steps)[i], ",",
                                c(prior$parameters$steps, 0)[i],"]"
        ))])
      })
    }else{
      names[[1]] <- bquote(~omega~"~"~.(print(prior, plot = TRUE)))
    }

    # using sampling aproximation for obtaining the values
    temp_out <- .plot_prior_data_omega(prior, samples, points, weights, return_samples = FALSE)
    df       <- temp_out$df
    prob     <- temp_out$prob

  }else if(prior$distribution == "two.sided" | (prior$distribution == "one.sided" & all(names(prior$parameters) %in% c("alpha", "steps")))){

    x_range <- c(0,1)
    if(weights){
      names   <- sapply(1:length(prior$parameters$alpha), function(i){
        bquote(~omega[~.(paste0("[", c(1, prior$parameters$steps)[i], ",",
                                c(prior$parameters$steps, 0)[i],"]"
        ))])
      })
    }else{
      names[[1]] <- bquote(~omega~"~"~.(print(prior, plot = TRUE)))
    }

    # analytical solution for the default cutoffs
    if(all(prior$parameters$alpha == 1) & length(prior$parameters$alpha) <= 3 & weights){

      if(length(prior$parameters$alpha) == 2){
        prob <- list(TRUE, FALSE)
        df   <- list(
          data.frame(
            x = c(1, 1),
            y = c(0, 1)
          ),
          data.frame(
            x = seq(0, 1, length.out = points),
            y = rep(1, points)
          )
        )
      }else if(length(prior$parameters$alpha) == 3){
        prob <- list(TRUE, FALSE, FALSE)
        df   <- list(
          data.frame(
            x = c(1, 1),
            y = c(0, 1)
          ),
          data.frame(
            x = seq(0, 1, length.out = points),
            y = seq(0, 1, length.out = points) * 2
          ),
          data.frame(
            x = seq(0, 1, length.out = points),
            y = -1*seq(0, 1, length.out = points)*2 + 2
          )
        )
      }

    }else{

      # sampling for more difficult weight-functions
      temp_out <- .plot_prior_data_omega(prior, samples, points, weights, return_samples = FALSE)
      df   <- temp_out$df
      prob <- temp_out$prob

    }

  }else if(prior$distribution == "point"){

    if(par_name == "omega"){
      df[[1]]  <- data.frame(
        x   = c(0, 1),
        y   = c(1, 1),
        lCI = c(1, 1),
        uCI = c(1, 1)
      )
      prob[[1]]  <- FALSE
      names[[1]] <- bquote(~omega~"~"~.(print(prior, plot = TRUE)))
    }else{
      prob[[1]] <- TRUE
      x_range   <- c(max(c(prior$truncation$lower, prior$parameters$location - .5)),
                     min(c(prior$truncation$upper, prior$parameters$location + .5)))
      df[[1]]   <- data.frame(
        x = rep(prior$parameters$location, 2),
        y = c(0, 1)
      )

      if(!is.null(mu_transform)){
        if(mu_transform == "cohens_d"){
          df[[1]]$x   <- psych::d2r(df[[1]]$x)
        }else if(mu_transform == "fishers_z"){
          df[[1]]$x   <- psych::fisherz2r(df[[1]]$x)
        }
      }

    }


  }else if(prior$distribution == "normal"){

    prob[[1]] <- FALSE
    if(is.null(mu_transform)){
      # analytical for non-transformed parameters
      if(is.null(x_range)){
        x_range <- c(
          if(!is.infinite(prior$truncation$lower)) prior$truncation$lower else prior$parameters$mean - 3*prior$parameters$sd,
          if(!is.infinite(prior$truncation$upper)) prior$truncation$upper else prior$parameters$mean + 3*prior$parameters$sd
        )
      }
      temp_x <- seq(x_range[1], x_range[2], length.out = points)
      df[[1]] <- data.frame(
        x = temp_x,
        y = exp(stats::dnorm(temp_x, mean = prior$parameters$mean, sd = prior$parameters$sd, log = TRUE) -
                  log(
                    stats::pnorm(prior$truncation$upper, prior$parameters$mean, prior$parameters$sd, lower.tail = TRUE, log.p = FALSE) -
                      stats::pnorm(prior$truncation$lower, prior$parameters$mean, prior$parameters$sd, lower.tail = TRUE, log.p = FALSE)
                  ))
      )
    }else{
      # sampling for transformed parameters
      x_range <- c(with_trunc$from, with_trunc$to)
      temp_x  <- NULL
      while(length(temp_x) < samples){
        new_temp_x <- stats::rnorm(samples, mean = prior$parameters$mean, sd = prior$parameters$sd)
        new_temp_x <- new_temp_x[new_temp_x > prior$truncation$lower & new_temp_x < prior$truncation$upper]
        temp_x     <- c(temp_x, new_temp_x)
      }
      # transform
      if(mu_transform == "cohens_d"){
        temp_x   <- psych::d2r(temp_x)
      }else if(mu_transform == "fishers_z"){
        temp_x   <- psych::fisherz2r(temp_x)
        temp_x[is.nan(temp_x)] <- ifelse(temp_x[is.nan(temp_x)] > 0, 1, -1)
      }
      temp_den <- stats::density(temp_x, from = with_trunc$from, to = with_trunc$to, n = points)
      temp_den$y[c(1, length(temp_den$y))] <- 0
      df[[1]] <- data.frame(
        x = temp_den$x,
        y = temp_den$y
      )
    }


  }else if(prior$distribution == "t"){

    prob[[1]] <- FALSE
    if(is.null(mu_transform)){
      # analytical for non-transformed parameters
      if(is.null(x_range)){
        x_range <- c(
          if(!is.infinite(prior$truncation$lower)) prior$truncation$lower else prior$parameters$location - 3*prior$parameters$scale,
          if(!is.infinite(prior$truncation$upper)) prior$truncation$upper else prior$parameters$location + 3*prior$parameters$scale
        )
      }
      temp_x <- seq(x_range[1], x_range[2], length.out = points)
      df[[1]] <- data.frame(
        x = temp_x,
        y = exp(extraDistr::dlst(temp_x, df = prior$parameters$df, mu = prior$parameters$location, sigma = prior$parameters$scale, log = TRUE) -
                  log(
                    extraDistr::plst(prior$truncation$upper, df = prior$parameters$df, mu = prior$parameters$location, sigma = prior$parameters$scale, lower.tail = TRUE, log.p = FALSE) -
                      extraDistr::plst(prior$truncation$lower, df = prior$parameters$df, mu = prior$parameters$location, sigma = prior$parameters$scale, lower.tail = TRUE, log.p = FALSE)
                  ))
      )
    }else{
      # sampling for transformed parameters
      x_range <- c(with_trunc$from, with_trunc$to)
      temp_x  <- NULL
      while(length(temp_x) < samples){
        new_temp_x <- extraDistr::rlst(samples, df = prior$parameters$df, mu = prior$parameters$location, sigma = prior$parameters$scale)
        new_temp_x <- new_temp_x[new_temp_x > prior$truncation$lower & new_temp_x < prior$truncation$upper]
        temp_x     <- c(temp_x, new_temp_x)
      }
      # transform
      if(mu_transform == "cohens_d"){
        temp_x   <- psych::d2r(temp_x)
      }else if(mu_transform == "fishers_z"){
        temp_x   <- psych::fisherz2r(temp_x)
        temp_x[is.nan(temp_x)] <- ifelse(temp_x[is.nan(temp_x)] > 0, 1, -1)
      }
      temp_den <- stats::density(temp_x, from = with_trunc$from, to = with_trunc$to, n = points)
      temp_den$y[c(1, length(temp_den$y))] <- 0
      df[[1]] <- data.frame(
        x = temp_den$x,
        y = temp_den$y
      )
    }

  }else if(prior$distribution == "gamma"){

    prob[[1]] <- FALSE
    if(is.null(mu_transform)){
      # analytical for non-transformed parameters
      if(is.null(x_range)){
        x_range <- c(
          if(!is.infinite(prior$truncation$lower)) prior$truncation$lower else stats::qgamma(.01, shape = prior$parameters$shape, rate = prior$parameters$rate),
          if(!is.infinite(prior$truncation$upper)) prior$truncation$upper else stats::qgamma(.95, shape = prior$parameters$shape, rate = prior$parameters$rate)
        )
      }
      if(prior$truncation$lower == 0 & x_range[1] < .5)x_range[1] <- 0
      temp_x <- seq(x_range[1], x_range[2], length.out = points)
      df[[1]] <- data.frame(
        x = temp_x,
        y = exp(stats::dgamma(temp_x, shape = prior$parameters$shape, rate = prior$parameters$rate, log = TRUE)  -
                  log(
                    stats::pgamma(prior$truncation$upper, shape = prior$parameters$shape, rate = prior$parameters$rate, lower.tail = TRUE, log.p = FALSE) -
                      stats::pgamma(prior$truncation$lower, shape = prior$parameters$shape, rate = prior$parameters$rate, lower.tail = TRUE, log.p = FALSE)
                  ))
      )
      if(df[[1]]$x[1] == 0)df[[1]]$y[1] <- 0
    }else{
      # sampling for transformed parameters
      x_range <- c(with_trunc$from, with_trunc$to)
      temp_x  <- NULL
      while(length(temp_x) < samples){
        new_temp_x <- stats::rgamma(samples, shape = prior$parameters$shape, rate = prior$parameters$rate)
        new_temp_x <- new_temp_x[new_temp_x > prior$truncation$lower & new_temp_x < prior$truncation$upper]
        temp_x     <- c(temp_x, new_temp_x)
      }
      # transform
      if(mu_transform == "cohens_d"){
        temp_x   <- psych::d2r(temp_x)
      }else if(mu_transform == "fishers_z"){
        temp_x   <- psych::fisherz2r(temp_x)
        temp_x[is.nan(temp_x)] <- ifelse(temp_x[is.nan(temp_x)] > 0, 1, -1)
      }
      temp_den <- stats::density(temp_x, from = with_trunc$from, to = with_trunc$to, n = points)
      temp_den$y[c(1, length(temp_den$y))] <- 0
      df[[1]] <- data.frame(
        x = temp_den$x,
        y = temp_den$y
      )
    }

  }else if(prior$distribution == "invgamma"){

    prob[[1]] <- FALSE
    if(is.null(mu_transform)){
      # analytical for non-transformed parameters
      if(is.null(x_range)){
        x_range <- c(
          if(!is.infinite(prior$truncation$lower)) prior$truncation$lower else extraDistr::qinvgamma(.01, alpha = prior$parameters$shape, beta = prior$parameters$scale),
          if(!is.infinite(prior$truncation$upper)) prior$truncation$upper else extraDistr::qinvgamma(.95, alpha = prior$parameters$shape, beta = prior$parameters$scale)
        )
      }
      if(prior$truncation$lower == 0 & x_range[1] < .5)x_range[1] <- 0
      temp_x <- seq(x_range[1], x_range[2], length.out = points)
      df[[1]] <- data.frame(
        x = temp_x,
        y = exp(extraDistr::dinvgamma(temp_x, alpha = prior$parameters$shape, beta = prior$parameters$scale, log = TRUE) -
                  log(
                    extraDistr::pinvgamma(prior$truncation$upper, alpha = prior$parameters$shape, beta = prior$parameters$scale, lower.tail = TRUE, log.p = FALSE) -
                      extraDistr::pinvgamma(prior$truncation$lower, alpha = prior$parameters$shape, beta = prior$parameters$scale, lower.tail = TRUE, log.p = FALSE)
                  ))
      )
      if(df[[1]]$x[1] == 0)df[[1]]$y[1] <- 0

    }else{
      # sampling for transformed parameters
      x_range <- c(with_trunc$from, with_trunc$to)
      temp_x  <- NULL
      while(length(temp_x) < samples){
        new_temp_x <- extraDistr::rinvgamma(samples, alpha = prior$parameters$shape, beta = prior$parameters$scale)
        new_temp_x <- new_temp_x[new_temp_x > prior$truncation$lower & new_temp_x < prior$truncation$upper]
        temp_x     <- c(temp_x, new_temp_x)
      }
      # transform
      if(mu_transform == "cohens_d"){
        temp_x   <- psych::d2r(temp_x)
      }else if(mu_transform == "fishers_z"){
        temp_x   <- psych::fisherz2r(temp_x)
        temp_x[is.nan(temp_x)] <- ifelse(temp_x[is.nan(temp_x)] > 0, 1, -1)
      }
      temp_den <- stats::density(temp_x, from = with_trunc$from, to = with_trunc$to, n = points)
      temp_den$y[c(1, length(temp_den$y))] <- 0
      df[[1]] <- data.frame(
        x = temp_den$x,
        y = temp_den$y
      )
    }
  }else if(prior$distribution == "uniform"){

    prob[[1]] <- FALSE
    if(is.null(mu_transform)){
      # analytical for non-transformed parameters
      if(is.null(x_range)){
        x_range <- c(
          prior$parameters$a,
          prior$parameters$b
        )
      }
      temp_x <- seq(x_range[1], x_range[2], length.out = points)
      df[[1]] <- data.frame(
        x = temp_x,
        y = exp(stats::dunif(temp_x, min = prior$parameters$a, max = prior$parameters$b, log = TRUE))
      )
      # add the 'truncation'
      if(x_range[1] == prior$parameters$a)df[[1]] <- rbind(c(x = prior$parameters$a, y = 0), df[[1]])
      if(x_range[2] == prior$parameters$b)df[[1]] <- rbind(df[[1]], c(x = prior$parameters$b, y = 0))
    }else{
      # sampling for transformed parameters
      x_range <- c(with_trunc$from, with_trunc$to)
      temp_x  <- stats::runif(samples, min = prior$parameters$a, max = prior$parameters$b)

      # transform
      if(mu_transform == "cohens_d"){
        temp_x   <- psych::d2r(temp_x)
      }else if(mu_transform == "fishers_z"){
        temp_x   <- psych::fisherz2r(temp_x)
        temp_x[is.nan(temp_x)] <- ifelse(temp_x[is.nan(temp_x)] > 0, 1, -1)
      }
      temp_den <- stats::density(temp_x, from = with_trunc$from, to = with_trunc$to, n = points)
      temp_den$y[c(1, length(temp_den$y))] <- 0
      df[[1]] <- data.frame(
        x = temp_den$x,
        y = temp_den$y
      )
    }
  }


  # add vertical edges to truncations, and correct all values above truncation for forced range variables
  if(prior$distribution %in% c("normal", "t", "gamma", "invgamma", "uniform")){
    if(!is.infinite(prior$truncation$lower)){
      if(any(df[[1]]$x < with_trunc$from)){
        df[[1]]$y[df[[1]]$x < with_trunc$from] <- 0
      }else if(!x_range_passed){
        df[[1]] <- rbind(c(with_trunc$from, 0), df[[1]])
      }
    }
    if(!is.infinite(with_trunc$to)){
      if(any(df[[1]]$x < with_trunc$to)){
        df[[1]]$y[df[[1]]$x > with_trunc$to] <- 0
      }else if(!x_range_passed){
        df[[1]] <- rbind(df[[1]], c(with_trunc$to, 0))
      }
    }
  }


  # add name if not weights
  if(length(names) == 0){
    if(par_name == ""){
      names[[1]] <- bquote(~"~"~.(print(prior, plot = TRUE)))
    }else if(par_name == "mu"){
      names[[1]] <- bquote(~mu~"~"~.(print(prior, plot = TRUE)))
    }else if(par_name == "tau"){
      names[[1]] <- bquote(~tau~"~"~.(print(prior, plot = TRUE)))
    }else if(par_name == "rho"){
      names[[1]] <- bquote(~rho~"~"~.(print(prior, plot = TRUE)))
    }
  }

  # compute some suplementary plotting info
  for(i in 1:length(df)){

    if((prior$distribution %in% c("two.sided", "one.sided") | par_name == "omega") & !weights){
      if(par_name == "omega"){
        x_lab[[i]] <- bquote(italic(p)*"-value")
      }else{
        x_lab[[i]] <- bquote(italic(p)*.(paste0("-value (",if(prior$distribution == "two.sided") "two-sided" else "one-sided",")")))
      }
      y_lab[[i]] <- names[[i]]

    }else{

      x_lab[[i]] <- names[[i]]
      y_lab[[i]] <- if(!prob[[i]]) "Density" else "Probability"

    }

    if(prob[[i]]){
      y_range[[i]] <- c(0, 1)
    }else if((prior$distribution %in% c("two.sided", "one.sided") | par_name == "omega") & !weights){
      y_range[[i]] <- c(0, 1)
    }else{
      y_range[[i]] <- c(0, max(df[[i]]$y[!is.nan(df[[i]]$y)]))
    }

  }


  return(list(
    df      = df,
    x_range = x_range,
    y_range = y_range,
    x_lab   = x_lab,
    y_lab   = y_lab,
    names   = names,
    prob    = prob
  ))
}
.plot_prior_data_omega <- function(prior, samples, points, weights, return_samples = FALSE){

  prob <- list()
  if(!weights){
    temp_mean <- NULL
    temp_lCI  <- NULL
    temp_uCI  <- NULL
    df        <- list()
  }else{
    if(all(names(prior$parameters) %in% c("alpha1", "alpha2", "steps"))){
      df  <- vector(mode = "list", length = length(prior$parameters$alpha1) + length(prior$parameters$alpha2) - 1)
    }else if(all(names(prior$parameters) %in% c("alpha", "steps"))){
      df  <- vector(mode = "list", length = length(prior$parameters$alpha))
    }
  }

  if(prior$distribution == "two.sided" | (prior$distribution == "one.sided" & all(names(prior$parameters) %in% c("alpha", "steps")))){

    # simulate data from dirichlet distribution
    eta <- NULL
    for(a in prior$parameters$alpha){
      eta <- cbind(eta, stats::rgamma(samples, shape = 1, rate = a))
    }
    eta <- eta / apply(eta, 1, sum)

    # transform to the cumulative sum
    omega <- t(apply(eta, 1, cumsum))

    # return those samples if requested
    if(return_samples){
      return(omega)
    }

    # create either weights or weight-function summaries
    for(i in 1:(length(prior$parameters$alpha)-1)){

      if(weights){
        prob[[length(prior$parameters$alpha) - (i-1)]] <- FALSE
        temp_d    <- stats::density(omega[,i], from = 0, to = 1)
        df[[length(prior$parameters$alpha) - (i-1)]]   <- data.frame(
          x = temp_d$x,
          y = temp_d$y
        )

      }else{
        temp_mean <- c(temp_mean, mean(omega[,i]))
        temp_lCI  <- c(temp_lCI,  stats::quantile(omega[,i], .025))
        temp_uCI  <- c(temp_uCI,  stats::quantile(omega[,i], .975))
      }


    }

    # add the last constant weight
    if(weights){
      prob[[1]] <- TRUE
      df[[1]]   <- data.frame(
        x = c(1, 1),
        y = c(1, 1)
      )
    }else{
      df[[1]]  <- data.frame(
        x   = c(0, prior$parameters$steps[sort(rep(1:length(prior$parameters$steps),2), decreasing = TRUE)], 1),
        y   = c(1, 1, temp_mean[sort(rep(1:length(temp_mean),2), decreasing = TRUE)]),
        lCI = c(1, 1, temp_lCI[sort(rep(1:length(temp_lCI),2),   decreasing = TRUE)]),
        uCI = c(1, 1, temp_uCI[sort(rep(1:length(temp_uCI),2),   decreasing = TRUE)])
      )
      prob[[1]] <- FALSE
    }

  }else if(prior$distribution == "one.sided" & all(names(prior$parameters) %in% c("alpha1", "alpha2", "steps"))){

    # simulate data from dirichlet distribution
    eta1  <- NULL
    eta2  <- NULL
    for(a in prior$parameters$alpha1){
      eta1 <- cbind(eta1, stats::rgamma(samples, shape = 1, rate = a))
    }
    for(a in prior$parameters$alpha2){
      eta2 <- cbind(eta2, stats::rgamma(samples, shape = 1, rate = a))
    }
    if(length(dim(eta1)) == 0)eta1 <- matrix(eta1, ncol = 1)
    if(length(dim(eta2)) == 0)eta2 <- matrix(eta2, ncol = 1)

    eta1 <- eta1/apply(eta1, 1, sum)
    eta2 <- eta2/apply(eta2, 1, sum) * (1-eta1[,1])

    # assign proper coefficients to the correct place
    omega <- matrix(nrow = samples, ncol = length(prior$parameters$alpha1) + length(prior$parameters$alpha2) - 1)
    for(j1 in 1:length(prior$parameters$alpha1)){
      omega[,length(prior$parameters$alpha2) - 1 + j1] = apply(matrix(eta1[,1:j1], ncol = j1), 1, sum)
    }
    for(j2 in 2:length(prior$parameters$alpha2)){
      omega[,j2-1] = apply(matrix(eta2[,j2:length(prior$parameters$alpha2)], ncol = length(prior$parameters$alpha2) - j2 + 1), 1, sum) + eta1[,1]
    }

    # return those samples if requested
    if(return_samples){
      return(omega)
    }

    # create either weights or weight-function summaries
    for(i in 1:(ncol(omega)-1)){

      if(weights){
        prob[[ncol(omega) - (i-1)]] <- FALSE
        temp_d    <- stats::density(omega[,i], from = 0, to = 1, n = points)
        df[[ncol(omega) - (i-1)]] <- data.frame(
          x = temp_d$x,
          y = temp_d$y
        )

      }else{
        temp_mean <- c(temp_mean, mean(omega[,i]))
        temp_lCI  <- c(temp_lCI,  stats::quantile(omega[,i], .025))
        temp_uCI  <- c(temp_uCI,  stats::quantile(omega[,i], .975))
      }


    }

    if(weights){
      prob[[1]] <- TRUE
      df[[1]]   <- data.frame(
        x = c(1, 1),
        y = c(1, 1)
      )
    }else{
      df[[1]]  <- data.frame(
        x   = c(0, prior$parameters$steps[sort(rep(1:length(prior$parameters$steps),2), decreasing = TRUE)], 1),
        y   = c(1, 1, temp_mean[sort(rep(1:length(temp_mean),2), decreasing = TRUE)]),
        lCI = c(1, 1, temp_lCI[sort(rep(1:length(temp_lCI),2),   decreasing = TRUE)]),
        uCI = c(1, 1, temp_uCI[sort(rep(1:length(temp_uCI),2),   decreasing = TRUE)])
      )
      prob[[1]] <- FALSE
    }
  }

  return(list(
    prob = prob,
    df   = df
  ))
}
.plot.prior_data_joint_omega <- function(models, type, samples = 1e6){

  omega_ind  <- .get_omega_mapping(models)
  priors     <- sapply(models, function(m)m$priors$omega, simplify = FALSE)
  prior_odds <- sapply(models, function(m)m$prior_odds)
  if(type == "conditional"){
    mm_par     <- sapply(models, function(m)!.is_parameter_null(m$priors, "omega"))
    priors     <- priors[mm_par]
    prior_odds <- prior_odds[mm_par]
    omega_ind  <- omega_ind[mm_par]
  }
  prior_prob   <- prior_odds / sum(prior_odds)

  # return NULL if only null models are specified
  if(all(sapply(omega_ind, is.null)))return(NULL)

  omega_samples <- matrix(nrow = 0, ncol = ncol(do.call(rbind,omega_ind)))

  for(i in 1:length(priors)){

    if(round(samples * prior_prob[i]) < 1)next
    if(priors[[i]]$distribution == "point"){

      omega_samples <- rbind(omega_samples, matrix(1, nrow = round(samples * prior_prob[i]), ncol = ncol(omega_samples)))

    }else{

      temp_out <- .plot_prior_data_omega(priors[[i]], round(samples * prior_prob[i]), NA, TRUE, TRUE)
      temp_sam <- NULL
      for(j in omega_ind[[i]]){
        temp_sam <- cbind(temp_sam, temp_out[,j])
      }

      omega_samples <- rbind(omega_samples, temp_sam)
    }
  }

  return(omega_samples)
}

#' @title Prints summary of \code{"RoBMA"} model implied by the specified priors
#'
#' @description \code{check_setup} prints summary of \code{"RoBMA"} model
#' implied by the specified priors. It is useful for checking the ensemble
#' configuration prior to fitting all of the models.
#'
#' @inheritParams RoBMA
#' @param models should the models details be printed
#' @export check_setup
#' @seealso [RoBMA()], [prior()]
check_setup <- function(priors_mu    = prior(distribution = "normal",   parameters = list(mean = 0, sd = 1)),
                        priors_tau   = prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15)),
                        priors_omega = list(
                          prior(distribution = "two.sided", parameters = list(alpha = c(1, 1),     steps = c(.05)),      prior_odds = 1/2),
                          prior(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),  steps = c(.05, .10)), prior_odds = 1/2)
                        ),
                        priors_mu_null    = prior(distribution = "point", parameters = list(location = 0)),
                        priors_tau_null   = prior(distribution = "point", parameters = list(location = 0)),
                        priors_omega_null = prior(distribution = "point", parameters = list(location = 1)),
                        models = FALSE){

  object <- list()
  object$priors  <- list(
    mu    = .set_parameter_priors(priors_mu_null,    priors_mu,    "mu"),
    tau   = .set_parameter_priors(priors_tau_null,   priors_tau,   "tau"),
    omega = .set_parameter_priors(priors_omega_null, priors_omega, "omega")
  )
  object$models  <- .get_models(object$priors)


  ### model types overview
  mm_mu       <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "mu"))
  mm_tau      <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "tau"))
  mm_omega    <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "omega"))

  # number of model types
  models_n    <- c(
    mu    = sum(mm_mu),
    tau   = sum(mm_tau),
    omega = sum(mm_omega)
  )

  # extract model weights
  prior_weights_all   <- sapply(object$models, function(m)m$prior_odds)
  # standardize model weights
  prior_weights_all   <- prior_weights_all / sum(prior_weights_all)
  # conditional model weights
  models_prior <- c(
    mu    <- sum(prior_weights_all[mm_mu]),
    tau   <- sum(prior_weights_all[mm_tau]),
    omega <- sum(prior_weights_all[mm_omega])
  )

  # create overview table
  overview_tab <- cbind.data.frame(models_n, models_prior)
  rownames(overview_tab) <- c("Effect", "Heterogeneity", "Pub. bias")
  colnames(overview_tab) <- c("Models", "Prior prob.")

  overview      <- overview_tab
  overview[,1]  <- paste0(overview[,1],"/",length(object$models))
  overview[,2]  <- format(round(overview[,2], 3), nsmall = 3)

  output <- list(
    overview = overview_tab,
    add_info = list(
      n_models = length(object$models)
    )
  )


  ### model details
  if(models){
    priors_mu      <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$mu, silent = TRUE))
    priors_tau     <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$tau, silent = TRUE))
    priors_omega   <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$omega, silent = TRUE))
    prior_odds     <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_odds)
    prior_prob     <- prior_odds / sum(prior_odds)

    models_tab <- cbind.data.frame(priors_mu, priors_tau, priors_omega, prior_prob, stringsAsFactors = FALSE)
    rownames(models_tab) <- NULL
    colnames(models_tab) <- c("Prior mu","Prior tau","Prior omega", "Prior prob.")

    models_print      <- models_tab
    models_print[,4]  <- format(round(models_print[,4], 3), nsmall = 3)

    output$models <- models_tab
  }



  cat("Robust Bayesian Meta-Analysis (Set-Up)\n")
  print(overview, quote = FALSE, right = TRUE)

  if(models){
    cat("\nModels Overview\n")
    print(models_print, quote = FALSE, right = TRUE)
  }


  return(invisible(output))

}

