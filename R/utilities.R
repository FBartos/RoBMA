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



# adapted from the runjags package version 2.2.0
RoBMA.private <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("defaultoptions",  list(
  jagspath = expression(runjags::findjags()),
  envir    = RoBMA.private))

assign("options",         RoBMA.private$defaultoptions,   envir = RoBMA.private)
assign("RoBMA_version",   "notset",                       envir = RoBMA.private)
assign("min_jags_major",  4,                              envir = RoBMA.private)
assign("max_jags_major",  4,                              envir = RoBMA.private)
assign("max_cores",       parallel::detectCores(logical = TRUE) - 1,  envir = RoBMA.private)

# check proper BayesTools package version
.check_BayesTools <- function(){

  RoBMA.version      <- try(utils::packageVersion("RoBMA"))
  BayesTools.version <- try(utils::packageVersion("BayesTools"))

  if(inherits(RoBMA.version, "try-error") | inherits(BayesTools.version, "try-error")){
    return(invisible(FALSE))
  }

  if(is.null(RoBMA.version) | is.null(BayesTools.version)){
    return(invisible(FALSE))
  }

  BayesTools_required <- switch(
    paste0(RoBMA.version, collapse = "."),
    "2.1.1" = c("0.1.3", "0.1.3"),
    "2.1.2" = c("0.1.3", "0.1.3"),
    "2.2.0" = c("0.1.3", "0.1.3"),
    "2.2.1" = c("0.2.3", "999.999.999"),
    "2.2.2" = c("0.2.3", "999.999.999"),
    "2.2.3" = c("0.2.3", "999.999.999"),
    "2.3.0" = c("0.2.3", "999.999.999"),
    "2.3.1" = c("0.2.3", "999.999.999"),
    stop("New RoBMA version needs to be defined in '.check_BayesTools' function!")
  )

  min_OK <- sum(as.numeric(strsplit(BayesTools_required[1], ".", fixed = TRUE)[[1]]) * c(1e9, 1e6, 1e3)) <=
    sum(unlist(BayesTools.version) * c(1e9, 1e6, 1e3))
  max_OK <- sum(as.numeric(strsplit(BayesTools_required[2], ".", fixed = TRUE)[[1]]) * c(1e9, 1e6, 1e3)) >=
    sum(unlist(BayesTools.version) * c(1e9, 1e6, 1e3))

  if(min_OK && max_OK){
    return(invisible(TRUE))
  }else{
    warning(sprintf(
      "RoBMA version %1$s requires BayesTools version higher or equal %2$s and lower or equal %3$s.",
      paste0(RoBMA.version, collapse = "."),
      BayesTools_required[1],
      BayesTools_required[2]
    ), call.= FALSE)
    return(invisible(FALSE))
  }
}

