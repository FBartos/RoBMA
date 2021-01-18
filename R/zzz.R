
.onLoad <- function(libname, pkgname){

  # load runjags
  requireNamespace("runjags")

  hereIsTheModule <- file.path(libname, pkgname)
  path <- file.path(hereIsTheModule, paste0("libs", Sys.getenv("R_ARCH")))
  tryCatch(rjags::load.module("RoBMA", path = path), error = function(e) warning(sprintf("The RoBMA module couldn't be loaded from %s. libname: %s, pkgname: %s.\n", path, libname, pkgname)))

  packageStartupMessage("Parametrization of random effect models changed and all models are estimated using likelihood of effects sizes since version 1.2.0, see NEWS for more details.")
}
