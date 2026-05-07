#' @importFrom graphics hist lines
#' @importFrom stats coef cooks.distance dfbetas fitted hatvalues influence logLik model.matrix nobs plogis predict qlogis qqnorm residuals rstandard rstudent terms update vcov
#' @importFrom utils capture.output getS3method
NULL

.onLoad <- function(libname, pkgname) {

  requireNamespace("BayesTools")
  requireNamespace("runjags")
  requireNamespace("mvtnorm")

  .check_bayestools_forward_api()

  RoBMA.private$RoBMA_version   <- utils::packageDescription(pkgname, fields = "Version")
  RoBMA.private$module_location <- .RoBMA_module_location(libname, pkgname)
  RoBMA.private$lib_name        <- libname

  .load_RoBMA_module(pkgname = pkgname, warn = TRUE)
  .load_RoBMA_native_routines(pkgname = pkgname, libname = libname, warn = TRUE)
  .check_RoBMA_native_routines(pkgname = pkgname)

  setopts <- mget(".RoBMA.options", envir = .GlobalEnv, ifnotfound = list(.RoBMA.options = NULL))[[1]]
  if (!is.null(setopts)) {
    if (!is.list(setopts)) {
      warning("Ignoring invalid non-list .RoBMA.options on loading RoBMA.",
              call. = FALSE)
    } else {
      do.call("RoBMA.options", args = setopts)
    }
  }

  .check_max_cores()
  .register_posterior_methods()
  .register_loo_methods()
}

.onAttach <- function(libname, pkgname) {

  packageStartupMessage(paste0(
    "Welcome to RoBMA ", utils::packageVersion(pkgname), ".\n",
    "See `vignette('v00-introduction', package = 'RoBMA')` for introduction to the package."
  ))
}

.check_bayestools_forward_api <- function() {

  required <- c(
    "formula_add_intercept",
    "plot_transformed_prior",
    "JAGS_formula_design"
  )
  missing <- required[!vapply(
    required,
    function(x) exists(x, envir = asNamespace("BayesTools"), inherits = FALSE),
    logical(1)
  )]

  posterior_plot_args         <- names(formals(BayesTools::plot_posterior))
  marginal_plot_args          <- names(formals(BayesTools::plot_marginal))
  missing_posterior_plot_args <- setdiff(
    c("data", "show_data", "dots_data"),
    posterior_plot_args
  )
  missing_marginal_plot_args  <- setdiff(
    c("legend", "legend_title", "legend_labels", "legend_position"),
    marginal_plot_args
  )

  if (length(missing) > 0 ||
      length(missing_posterior_plot_args) > 0 ||
      length(missing_marginal_plot_args) > 0) {
    details <- character(0)
    if (length(missing) > 0) {
      details <- c(details, paste0("missing functions: ", paste(missing, collapse = ", ")))
    }
    if (length(missing_posterior_plot_args) > 0) {
      details <- c(
        details,
        paste0(
          "BayesTools::plot_posterior() missing arguments: ",
          paste(missing_posterior_plot_args, collapse = ", ")
        )
      )
    }
    if (length(missing_marginal_plot_args) > 0) {
      details <- c(
        details,
        paste0(
          "BayesTools::plot_marginal() missing arguments: ",
          paste(missing_marginal_plot_args, collapse = ", ")
        )
      )
    }

    stop(
      "RoBMA requires a BayesTools build with the forward APIs (",
      paste(details, collapse = "; "),
      ").",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.onUnload <- function(libpath) {

  tryCatch(
    {
      if ("RoBMA" %in% rjags::list.modules()) {
        rjags::unload.module("RoBMA")
      }
    },
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
                               quiet = TRUE, warn = FALSE) {

  loaded <- tryCatch("RoBMA" %in% rjags::list.modules(), error = function(e) FALSE)
  if (loaded) {
    return(TRUE)
  }

  if (is.null(path) || !dir.exists(path)) {
    if (warn) {
      warning(
        "RoBMA JAGS module was not found in the installed package library. ",
        "Model fitting requires this module; reinstall RoBMA after installing JAGS >= 4.3.1.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  load_error <- NULL
  tryCatch(
    rjags::load.module("RoBMA", path = path, quiet = quiet),
    error = function(e) load_error <<- conditionMessage(e)
  )

  loaded <- tryCatch("RoBMA" %in% rjags::list.modules(), error = function(e) FALSE)
  if (!loaded && warn) {
    message <- paste0(
      "RoBMA JAGS module failed to load from '", path, "'. ",
      "Model fitting requires this module; reinstall RoBMA after installing JAGS >= 4.3.1."
    )
    if (!is.null(load_error)) {
      message <- paste0(message, " rjags error: ", load_error)
    }
    warning(message, call. = FALSE)
  }

  return(loaded)
}

.load_RoBMA_native_routines <- function(pkgname = "RoBMA", libname = RoBMA.private$lib_name,
                                        warn = FALSE) {

  if (isTRUE(.check_RoBMA_native_routines(pkgname = pkgname, warn = FALSE))) {
    return(TRUE)
  }

  load_error <- NULL
  tryCatch(
    library.dynam("RoBMA", pkgname, libname),
    error = function(e) load_error <<- conditionMessage(e)
  )

  if (!isTRUE(.check_RoBMA_native_routines(pkgname = pkgname, warn = FALSE)) &&
      !is.null(RoBMA.private$module_location)) {
    module_file <- file.path(RoBMA.private$module_location, paste0("RoBMA", .Platform$dynlib.ext))
    if (file.exists(module_file)) {
      tryCatch(
        dyn.load(module_file),
        error = function(e) load_error <<- conditionMessage(e)
      )
    }
  }

  loaded <- .check_RoBMA_native_routines(pkgname = pkgname, warn = FALSE)
  if (!isTRUE(loaded) && warn) {
    message <- paste0(
      "RoBMA native routines failed to load from the package DLL. ",
      "Compiled likelihood helpers will be unavailable; reinstall RoBMA after installing JAGS >= 4.3.1."
    )
    if (!is.null(load_error)) {
      message <- paste0(message, " R loader error: ", load_error)
    }
    warning(message, call. = FALSE)
  }

  return(isTRUE(loaded))
}

.check_RoBMA_native_routines <- function(pkgname = "RoBMA", warn = TRUE) {

  required_symbols <- c(
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
