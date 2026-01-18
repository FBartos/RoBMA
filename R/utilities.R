#' @title Options for the RoBMA package
#'
#' @description A placeholder object and functions for the RoBMA package.
#' (adapted from the runjags R package).
#'
#' @param name the name of the option to get the current value of - for a list of
#' available options, see details below.
#' @param ... named option(s) to change - for a list of available options, see
#' details below.
#'
#' @details
#' The available options are:
#' \describe{
#'   \item{\code{max_cores}}{number of cores to use for parallel computing (default to all available cores - 1)}
#'   \item{\code{check_scaling}}{whether to check scaling of predictors (default \code{TRUE})}
#'   \item{\code{silent}}{whether to suppress output (default \code{FALSE})}
#'   \item{\code{autocompute.loo}}{whether to automatically compute LOO (default \code{FALSE})}
#'   \item{\code{autocompute.waic}}{whether to automatically compute WAIC (default \code{FALSE})}
#'   \item{\code{autocompute.marglik}}{whether to automatically compute marginal likelihood (default \code{FALSE})}
#'   \item{\code{default_UISD.effect}}{default scaling of the unit information standard deviation for the effect size parameter (default \code{0.5})}
#'   \item{\code{default_UISD.heterogeneity}}{default scaling of the unit information standard deviation for the heterogeneity parameter (default \code{0.25})}
#'   \item{\code{default_UISD.mods}}{default scaling of the unit information standard deviation for the moderators (default \code{0.25})}
#'   \item{\code{default_UISD.scale}}{default scaling of the unit information standard deviation for the scale parameter (default \code{0.5})}
#'   \item{\code{default_informed_priors.mods}}{default scaling of informed priors for moderators (default \code{0.5})}
#'   \item{\code{default_informed_priors.scale}}{default scaling of informed priors for the scale parameter (default \code{0.5})}
#'   \item{\code{default_bias_weightfunction.alpha}}{default alpha for the weightfunction (default \code{1})}
#'   \item{\code{default_bias_PET.scale}}{default scale for the PET (default \code{1})}
#'   \item{\code{default_bias_PEESE.scale}}{default scale for the PEESE (default \code{5})}
#' }
#'
#' @return The current value of all available RoBMA options (after applying any
#' changes specified) is returned invisibly as a named list.
#'
#' @export RoBMA.options
#' @export RoBMA.get_option
#' @name RoBMA_options
#' @aliases RoBMA_options RoBMA.options RoBMA.get_option
NULL


#' @rdname RoBMA_options
RoBMA.options    <- function(...){

	opts <- list(...)

	for(i in seq_along(opts)){

	  if(!names(opts)[i] %in% names(RoBMA.private))
	    stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep=""))

	  assign(names(opts)[i], opts[[i]] , envir = RoBMA.private)
	}

	return(invisible(RoBMA.private$options))
}

#' @rdname RoBMA_options
RoBMA.get_option <- function(name){

	if(length(name)!=1)
	  stop("Only 1 option can be retrieved at a time")

	if(!name %in% names(RoBMA.private))
	  stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))

	# Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
	return(eval(RoBMA.private[[name]]))
}

# export the function directly to suppress import warnings
.runjags__findjags <- function() runjags::findjags()

# adapted from the runjags package version 2.2.0
RoBMA.private <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("defaultoptions",  list(
  jagspath = expression(.runjags__findjags()),
  envir    = RoBMA.private))

assign("options",         RoBMA.private$defaultoptions,   envir = RoBMA.private)
assign("RoBMA_version",   "notset",                       envir = RoBMA.private)
assign("min_jags_major",  4,                              envir = RoBMA.private)
assign("max_jags_major",  4,                              envir = RoBMA.private)
assign("max_cores",       parallel::detectCores(logical = TRUE) - 1,  envir = RoBMA.private)
assign("check_scaling",   TRUE,                                       envir = RoBMA.private) # TODO: remove?

assign("silent", FALSE,  envir = RoBMA.private)

assign("autocompute.loo",     FALSE, envir = RoBMA.private)
assign("autocompute.waic",    FALSE, envir = RoBMA.private)
assign("autocompute.marglik", FALSE, envir = RoBMA.private)


### default scaling of unit information prior for different parameters
assign("default_UISD.effect",          1/2,  envir = RoBMA.private)
assign("default_UISD.heterogeneity",   1/4,  envir = RoBMA.private)
assign("default_UISD.mods",            1/4,  envir = RoBMA.private)
assign("default_UISD.scale",           1/2,  envir = RoBMA.private)

### default scaling of informed priors for moderators
assign("default_informed_priors.mods",  1/2,  envir = RoBMA.private)
assign("default_informed_priors.scale", 1/2,  envir = RoBMA.private)

### default setting for the publication bias priors
assign("default_bias_weightfunction.alpha",  1,  envir = RoBMA.private)
assign("default_bias_PET.scale",             1,  envir = RoBMA.private)
assign("default_bias_PEESE.scale",           5,  envir = RoBMA.private)


# check and fix number of threads (sometimes bugs out during installation)
.check_max_cores <- function(){
  if(RoBMA.private$max_cores > parallel::detectCores(logical = TRUE) - 1){
    RoBMA.options(max_cores = parallel::detectCores(logical = TRUE) - 1)
  }
}
