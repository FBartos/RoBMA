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

	if(length(opts) > 0){

		options    <- RoBMA.private$options
		recognised <- pmatch(names(opts), names(options))

		if(any(is.na(recognised))){
			warning(paste("Igoring unmatched or ambiguous option(s): ", paste(names(opts)[is.na(recognised)], collapse=", ")))
      opts <- opts[!is.na(recognised)]
		}

		optnames <- names(options)[recognised[!is.na(recognised)]]

		if(length(optnames)>0) for(i in 1:length(optnames)){
			options[optnames[i]] <- opts[[i]]
		}

		assign("options", options , envir = RoBMA.private)
	}

	return(invisible(RoBMA.private$options))
}

#' @rdname RoBMA_options
RoBMA.get_option <- function(name){

	if(length(name)!=1)
	  stop("Only 1 option can be retrieved at a time")

	opt <- pmatch(name,names(RoBMA.private$options))

	if(is.na(opt))
	  stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))

	# Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
	return(eval(RoBMA.private$options[[opt]]))
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
