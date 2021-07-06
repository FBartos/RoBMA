# adapted from the runjags package version 2.2.0
.onLoad <- function(libname, pkgname){

  RoBMA.private$RoBMA_version <- utils::packageDescription(pkgname, fields='Version')

  # Get and save the library location, getting rid of any trailing / caused by r_arch being empty:
  modloc <- gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
  if(!file.exists(file.path(modloc, paste('RoBMA', .Platform$dynlib.ext, sep='')))){
    modloc <- ''
    warning('The RoBMA module could not be loaded.', call. =FALSE)
  }else{
    rjags::load.module("RoBMA", path = modloc)
    library.dynam("RoBMA", "RoBMA", libname)
  }

  RoBMA.private$modulelocation <- modloc

  setopts <- mget('.RoBMA.options', envir=.GlobalEnv, ifnotfound = list(.RoBMA.options = NULL))[[1]]
  if(!is.null(setopts)){
    if(!is.list(setopts)){
      warning('Ignoring invalid (non-list) specification for .RoBMA.options on loading the RoBMA package', call.=FALSE)
    }else{
      newopts <- do.call('RoBMA.options', args = setopts)
    }
  }

}

.onAttach <- function(libname, pkgname){

  packageStartupMessage("This is a preview of the 2.0 version of the RoBMA package.")
  # 1.2.0 message: "Please, note the following changes in version 1.2.0 (see NEWS for more details):\n- all models are now estimated using the likelihood of effects sizes (instead of t-statistics)\n- parametrization of random effect models changed (the studies' true effects are marginalized out of the likelihood)"
}

.onUnload <- function(libpath){

  # Just in case it is not always safe to try and access an element of an env that is in the process of being deleted (when R quits):
  if(!is.null(RoBMA.private$modulelocation)){
    rjags::unload.module("RoBMA")
    library.dynam.unload("RoBMA", libpath)
  }
}
