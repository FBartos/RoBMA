# adapted from the runjags package version 2.2.0
.onLoad <- function(libname, pkgname){

  requireNamespace("runjags")
  requireNamespace("mvtnorm")

  RoBMA.private$RoBMA_version <- utils::packageDescription(pkgname, fields = 'Version')

  # Get and save the library location, getting rid of any trailing / caused by r_arch being empty:
  module_location <- gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
  if(!file.exists(file.path(module_location, paste('RoBMA', .Platform$dynlib.ext, sep='')))){
    module_location <- NULL
    warning('The RoBMA module could not be loaded.', call. = FALSE)
  }else{
    rjags::load.module("RoBMA", path = module_location, quiet = TRUE)
    if(!"RoBMA" %in% rjags::list.modules()){
      warning('The RoBMA module could not be loaded.', call. = FALSE)
    }
  }

  RoBMA.private$module_location <- module_location
  RoBMA.private$lib_name        <- libname

  setopts <- mget('.RoBMA.options', envir=.GlobalEnv, ifnotfound = list(.RoBMA.options = NULL))[[1]]
  if(!is.null(setopts)){
    if(!is.list(setopts)){
      warning('Ignoring invalid (non-list) specification for .RoBMA.options on loading the RoBMA package', call.=FALSE)
    }else{
      newopts <- do.call('RoBMA.options', args = setopts)
    }
  }

  # Check if correct version of BayesTools is installed
  .check_BayesTools()

  # Check and fix number of threads (sometimes bugs out during installation)
  .check_max_cores()

}

.onAttach <- function(libname, pkgname){

  packageStartupMessage("RoBMA version 3.3 now features spike-and-slab style model-averaging via the \'algorithm = \"ss\"\' argument.\nSee \'vignette(\"FastRoBMA\", package = \"RoBMA\")' for more details (\'algorithm = \"ss\"\' argument will become the default setting in the future major release of the package).")

}

.onUnload <- function(libpath){

  # tricking the dyn.library unload
  if(!is.null(RoBMA.private$lib_name)){
    library.dynam("RoBMA", "RoBMA", RoBMA.private$lib_name)
  }

  # Just in case it is not always safe to try and access an element of an env that is in the process of being deleted (when R quits):
  if(!is.null(RoBMA.private$module_location)){
    rjags::unload.module("RoBMA")
  }
}

.load_RoBMA_module <- function(pkgname = "RoBMA"){

  if(is.null(RoBMA.private$module_location) || (!is.null(RoBMA.private$module_location) && RoBMA.private$module_location == "")){
    libnames         <- .libPaths()
    module_locations <- sapply(libnames, function(libname) gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else "")))
    sapply(module_locations, function(module_location) rjags::load.module("RoBMA", path = module_location))
  }else{
    rjags::load.module("RoBMA", path = RoBMA.private$module_location)
  }

}
