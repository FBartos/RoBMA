#' @importFrom graphics hist lines
#' @importFrom Matrix Cholesky
#' @importFrom stats AIC BIC coef cooks.distance dfbetas fitted hatvalues influence logLik model.matrix nobs plogis predict qlogis qqnorm residuals rstandard rstudent terms update
#' @importFrom utils capture.output getFromNamespace getS3method tail
NULL

# The Matrix import loads the registered CHOLMOD C-callables used by r-known-v.

# Saved at installation to check the runtime ABI before loading our module.
.RoBMA_jags_build <- rjags::pkg.compile.flags("ALL")

.onLoad <- function(libname, pkgname) {

  RoBMA.private[["selection_runtime_loading"]] <- TRUE
  RoBMA.private[["selection_runtime_initialized"]] <- FALSE
  on.exit(RoBMA.private[["selection_runtime_loading"]] <- FALSE, add = TRUE)
  requireNamespace("BayesTools")
  requireNamespace("mvtnorm")

  RoBMA.private$RoBMA_version   <- utils::packageDescription(pkgname, fields = "Version")
  RoBMA.private$module_location <- .RoBMA_module_location(libname, pkgname)
  RoBMA.private$lib_name        <- libname

  source_only <- is.null(RoBMA.private$module_location) &&
    dir.exists(file.path(getNamespaceInfo(pkgname, "path"), "src"))
  if (!source_only) {
    .load_RoBMA_module(pkgname = pkgname)
    if (!.check_RoBMA_native_routines(pkgname = pkgname, warn = FALSE)) {
      stop("The RoBMA native routines are missing from its JAGS module. Reinstall RoBMA.", call. = FALSE)
    }
  }

  setopts <- mget(".RoBMA.options", envir = .GlobalEnv, ifnotfound = list(.RoBMA.options = NULL))[[1]]
  if (!is.null(setopts)) {
    if (!is.list(setopts)) {
      warning("Ignoring invalid non-list .RoBMA.options on loading RoBMA.",
              call. = FALSE)
    } else {
      do.call("RoBMA.options", args = setopts)
    }
  }
  .native_threads_configure(RoBMA.private[["native_threads"]])

  if (.selection_runtime_available()) {
    .selection_runtime_configure(.selection_runtime_settings(capacity_bytes = 0))
    RoBMA.private[["selection_runtime_initialized"]] <- FALSE
  }
  .check_max_cores()
  .register_posterior_methods()
  .register_loo_methods()
  .hypothesis_density_backends_register()
}

.onAttach <- function(libname, pkgname) {

  packageStartupMessage(paste0(
    "Welcome to RoBMA ", utils::packageVersion(pkgname), ".\n",
    "See `vignette('v00-introduction', package = 'RoBMA')` for introduction to the package."
  ))
}

.onUnload <- function(libpath) {

  # Run finalizers for unused models while their samplers are still registered.
  gc()
  tryCatch(
    {
      if ("RoBMA" %in% rjags::list.modules()) {
        rjags::unload.module("RoBMA", quiet = TRUE)
      }
    },
    error = function(e) NULL
  )
  tryCatch(
    .hypothesis_density_backends_unregister(),
    error = function(e) NULL
  )
}

.RoBMA_module_location <- function(libname, pkgname) {

  arch             <- if (.Platform$r_arch != "") .Platform$r_arch else ""
  module_location  <- file.path(libname, pkgname, "libs", arch)
  module_file      <- file.path(module_location, paste0("RoBMA", .Platform$dynlib.ext))

  if (file.exists(module_file)) {
    return(normalizePath(module_location, winslash = "/", mustWork = TRUE))
  }

  source_location <- .RoBMA_source_module_location(pkgname)
  if (!is.null(source_location)) {
    return(source_location)
  }

  dll_path <- .RoBMA_loaded_dll_path(pkgname)
  if (is.null(dll_path)) {
    return(NULL)
  }

  return(dirname(dll_path))
}

.RoBMA_source_module_location <- function(pkgname) {

  source_path <- tryCatch(getNamespaceInfo(pkgname, "path"), error = function(e) NULL)
  if (is.null(source_path)) {
    return(NULL)
  }

  arch <- if (.Platform$r_arch != "") .Platform$r_arch else ""
  paths <- c(
    file.path(source_path, "src", arch),
    file.path(source_path, "src")
  )
  paths <- unique(paths[nzchar(paths)])

  for (path in paths) {
    module_file <- file.path(path, paste0("RoBMA", .Platform$dynlib.ext))
    if (file.exists(module_file)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
  }

  return(NULL)
}

.RoBMA_loaded_dll_path <- function(pkgname) {

  dlls <- getLoadedDLLs()
  if (!pkgname %in% names(dlls)) {
    return(NULL)
  }

  dll_path <- dlls[[pkgname]][["path"]]
  if (is.null(dll_path) || !file.exists(dll_path)) {
    return(NULL)
  }

  return(normalizePath(dll_path, winslash = "/", mustWork = TRUE))
}

.load_RoBMA_module <- function(pkgname = "RoBMA", path = RoBMA.private$module_location,
                               quiet = TRUE) {

  if (pkgname %in% rjags::list.modules()) {
    return(invisible(TRUE))
  }
  if (is.null(path)) {
    stop("The RoBMA JAGS module is missing. Reinstall RoBMA using JAGS 5 or later.", call. = FALSE)
  }

  # RoBMA's factories must be registered after the BayesTools factories.
  if (!BayesTools::BayesTools_load_JAGS_module(quiet = quiet)) {
    stop("The BayesTools JAGS module is required to load RoBMA.", call. = FALSE)
  }
  libname <- RoBMA.private$lib_name
  installed_dll <- file.path(libname, pkgname, "libs", .Platform$r_arch,
                             paste0(pkgname, .Platform$dynlib.ext))
  if (file.exists(installed_dll)) {
    rjags::load.pkg.module(pkgname, libname, .RoBMA_jags_build, quiet = quiet)
  } else {
    # pkgload places the compiled module under src/.
    rjags::load.pkg.module(pkgname, libname, .RoBMA_jags_build, check.only = TRUE)
    rjags::load.module(pkgname, path = path, quiet = quiet)
  }
  invisible(TRUE)
}

.check_RoBMA_native_routines <- function(pkgname = "RoBMA", warn = TRUE) {

  required_symbols <- c(
    "RoBMA_selnorm_cache_control",
    "RoBMA_selnorm_sampler_control",
    "RoBMA_selnorm_sampler_report",
    "RoBMA_selnorm_kernel_loglik_matrix",
    "RoBMA_glmm_binom_marginal_loglik",
    "RoBMA_glmm_pois_marginal_loglik"
  )

  loaded <- vapply(
    required_symbols,
    is.loaded,
    FUN.VALUE = logical(1),
    PACKAGE   = pkgname
  )

  if (!all(loaded) && warn) {
    warning(
      "RoBMA native routines are not loaded from the package DLL. ",
      "Compiled likelihood helpers will be unavailable.",
      call. = FALSE
    )
  }

  return(invisible(all(loaded)))
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
