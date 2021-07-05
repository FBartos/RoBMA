# addapted from the runjags package version 2.2.0
.onLoad <- function(libname, pkgname){

  RoBMA.private$RoBMA_version <- utils::packageDescription(pkgname, fields='Version')

  # Get and save the library location, getting rid of any trailing / caused by r_arch being empty:
  modloc <- gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
  if(!file.exists(file.path(modloc, paste('RoBMA', .Platform$dynlib.ext, sep=''))))
    modloc <- ''

  RoBMA.private$modulelocation <- modloc

  setopts <- mget('.RoBMA.options', envir=.GlobalEnv, ifnotfound = list(.RoBMA.options = NULL))[[1]]
  if(!is.null(setopts)){
    if(!is.list(setopts)){
      warning('Ignoring invalid (non-list) specification for .RoBMA.options on loading the RoBMA package', call.=FALSE)
    }else{
      newopts <- do.call('RoBMA.options', args = setopts)
    }
  }

  # To ensure that cleanup.jags is run when R is quit:
  reg.finalizer(RoBMA.private, .onDetach, onexit = TRUE)

}

.onAttach <- function(libname, pkgname){

  # This will be run after load if the package is attached:
  setopts <- mget('.RoBMA.options', envir=.GlobalEnv, ifnotfound=list(.RoBMA.options=NULL))[[1]]
  if(!is.null(setopts) && !RoBMA.getOption('silent.RoBMA')){
    packageStartupMessage(paste('Attaching RoBMA (version ', RoBMA.private$RoBMA_version, ') and setting user-specified options', sep=''))
  }

  packageStartupMessage("This is a preview of the 2.0 version of the RoBMA package.")
  # 1.2.0 message: "Please, note the following changes in version 1.2.0 (see NEWS for more details):\n- all models are now estimated using the likelihood of effects sizes (instead of t-statistics)\n- parametrization of random effect models changed (the studies' true effects are marginalized out of the likelihood)"
}

.onDetach <- function(libpath){

  # Just in case it is not always safe to try and access an element of an env that is in the process of being deleted (when R quits):
  if(!is.null(RoBMA.private$dynlibname)){
    .dynunloadmodule()
  }

  all.folders <- try(runjags::runjags.getOption('full.cleanup'), silent = TRUE)

  if(class(all.folders)=='try-error'){
    all.folders <- FALSE
  }

  try(runjags::cleanup.jags(all.folders = all.folders, silent = TRUE), silent = TRUE)
}

.dynunloadmodule <- function(){

  if(is.null(runjagsprivate$dynlibname)){
    warning('Unable to load the dynlib as it has not been loaded')
    invisible(FALSE)
  }

  slibpath <- system.file("libs", paste(.Platform$r_arch, if(.Platform$r_arch!="") "/" else "", if(.Platform$OS.type=="unix") "runjags.so" else "runjags.dll", sep=""), package="RoBMA")
  swcat("Unloading shared library from:  ", slibpath, "\n", sep="")
  success <- try(dyn.unload(slibpath))

  if(inherits(success, 'try-error'))
    warning("The internal dynlib could not be unloaded")

  RoBMA.private$dynlibname <- NULL
  return(invisible(TRUE))
}
