
.onLoad <- function(libname, pkgname){

  # load runjags
  requireNamespace("BayesTools")
  requireNamespace("runjags")

  hereIsTheModule <- file.path(libname, pkgname)
  path <- file.path(hereIsTheModule, paste0("libs", Sys.getenv("R_ARCH")))
  tryCatch(rjags::load.module("RoBMA", path = path), error = function(e) warning(sprintf("The RoBMA module couldn't be loaded from %s. libname: %s, pkgname: %s.\n", path, libname, pkgname)))

}

.onAttach <- function(libname, pkgname){

  packageStartupMessage("This is a preview of the 2.0 version of the RoBMA package.")
  # 1.2.0 message: "Please, note the following changes in version 1.2.0 (see NEWS for more details):\n- all models are now estimated using the likelihood of effects sizes (instead of t-statistics)\n- parametrization of random effect models changed (the studies' true effects are marginalized out of the likelihood)"

}

# .onUnload <- function(libpath) {
#
#   tryCatch(rjags::unload.module("RoBMA"), error = function(e) warning(sprintf("The RoBMA module couldn't be unloaded.")))
# }
