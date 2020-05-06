.onLoad <- function(libname, pkgname) {
  # load runjags
  requireNamespace("runjags")

  hereIsTheModule <- sub("/DESCRIPTION", '', sub("/Meta.*", '', attr(packageDescription("RoBMA"), "file")))
  # load the JAGS module - it always fails during the instalation process for some reason
  tryCatch(rjags::load.module("RoBMA", path = paste0(hereIsTheModule, "/libs", Sys.getenv("R_ARCH")) ), error = function(e)cat("The RoBMA module couldn't be loaded.\n"))

}
