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

  # Check and fix number of threads (sometimes bugs out during installation)
  .check_max_cores()

  # Register S3 methods for posterior package (if available)
  .register_posterior_methods()

  # Register S3 methods for loo package (if available)
  .register_loo_methods()

}

.onAttach <- function(libname, pkgname){

  packageStartupMessage("Welcome to RoBMA 4.0")

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


# ---------------------------------------------------------------------------- #
# S3 method registration for posterior package
# ---------------------------------------------------------------------------- #

# Register as_draws methods with the posterior package when it is available.
# This allows posterior::as_draws(brma_object) to work correctly.
.register_posterior_methods <- function() {

  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(invisible(NULL))
  }

  # Register methods for brma class
  .s3_register("posterior::as_draws",        "brma")
  .s3_register("posterior::as_draws_array",  "brma")
  .s3_register("posterior::as_draws_df",     "brma")
  .s3_register("posterior::as_draws_list",   "brma")
  .s3_register("posterior::as_draws_matrix", "brma")
  .s3_register("posterior::as_draws_rvars",  "brma")

  # Register methods for brma_samples class
  .s3_register("posterior::as_draws",        "brma_samples")
  .s3_register("posterior::as_draws_array",  "brma_samples")
  .s3_register("posterior::as_draws_df",     "brma_samples")
  .s3_register("posterior::as_draws_list",   "brma_samples")
  .s3_register("posterior::as_draws_matrix", "brma_samples")
  .s3_register("posterior::as_draws_rvars",  "brma_samples")

  invisible(NULL)
}

# Register loo and loo_compare methods with the loo package when it is available.
# This allows loo::loo(brma_object) and loo::loo_compare(brma_object) to work correctly.
.register_loo_methods <- function() {

  if (!requireNamespace("loo", quietly = TRUE)) {
    return(invisible(NULL))
  }

  # Register methods for brma class
  .s3_register("loo::loo",         "brma")
  .s3_register("loo::loo_compare", "brma")
  .s3_register("loo::loo_compare", "loo")

  invisible(NULL)
}

# S3 method registration helper (adapted from vctrs package)
# This is a standard pattern used by tidyverse and other packages to
# conditionally register S3 methods with another package's generics.
.s3_register <- function(generic, class, method = NULL) {

  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]

  caller <- parent.frame()

  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) {
      top
    } else {
      caller
    }
  }
  get_method <- function(method, env) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = env)
    } else {
      method
    }
  }

  register <- function(...) {
    envir <- asNamespace(package)
    method_fn <- get_method(method, get_method_env())
    registerS3method(generic, class, method_fn, envir = envir)
  }

  setHook(packageEvent(package, "onLoad"), register)

  if (isNamespaceLoaded(package)) {
    register()
  }

  invisible()
}
