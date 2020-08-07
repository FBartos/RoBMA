
.onLoad <- function(libname, pkgname){

  # load runjags
  requireNamespace("runjags")

  hereIsTheModule <- file.path(libname, pkgname)
  path <- file.path(hereIsTheModule, paste0("libs", Sys.getenv("R_ARCH")))
  tryCatch(rjags::load.module("RoBMA", path = path), error = function(e) warning(sprintf("The RoBMA module couldn't be loaded from %s. libname: %s, pkgname: %s.\n", path, libname, pkgname)))

}
