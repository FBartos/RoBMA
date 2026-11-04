# Build the development DLL with the package's release flags before a runner
# loads it through devtools/pkgload.
#
# pkgbuild::compile_dll(), which devtools::load_all() and devtools::test() call,
# compiles with its debug flags by default (-UNDEBUG -g -O0 appended after the
# package's own -O2, so -O0 wins). Every timing a runner records under
# load_all() is therefore a timing of unoptimized native code, and work that
# moves computation from optimized library calls into package C++ (the
# selection tail kernel, the quadrature rules) looks slower there while being
# several times faster in the installed build. The scenario and test runners
# source this file so that the DLL they load is the one users get.
#
# The DLL is rebuilt with the release flags when the sources changed or when
# the last build recorded here was not a release build (for example after an
# interactive load_all() in an IDE). Set ROBMA_DEBUG_DLL=TRUE to keep the
# debug build for a debugging session.

ensure_optimized_dll <- function(path = ".", quiet = TRUE) {

  path <- normalizePath(path, mustWork = TRUE)
  if (identical(toupper(Sys.getenv("ROBMA_DEBUG_DLL", "")), "TRUE")) {
    return(invisible("debug"))
  }
  if (!requireNamespace("pkgbuild", quietly = TRUE)) {
    return(invisible("unavailable"))
  }

  dll    <- file.path(path, "src", paste0("RoBMA", .Platform$dynlib.ext))
  marker <- file.path(path, "src", ".optimized-dll")
  # A release DLL is one this helper built and that nothing rebuilt since:
  # any later compile (an IDE load_all(), a debug run) leaves a newer DLL.
  release <- file.exists(dll) && file.exists(marker) &&
    file.info(dll)[["mtime"]] <= file.info(marker)[["mtime"]] + 1

  if (!release) {
    # A rebuild must start from clean objects: make only recompiles sources
    # newer than their objects, so debug objects would otherwise survive into
    # the "release" DLL and only the link step would run.
    pkgbuild::clean_dll(path)
  }
  pkgbuild::compile_dll(path, force = !release, debug = FALSE, quiet = quiet)
  # compile_dll() rebuilds only changed objects; the marker follows the DLL.
  if (file.exists(dll)) {
    writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), marker)
    Sys.setFileTime(marker, file.info(dll)[["mtime"]] + 1)
  }
  invisible(if (release) "kept" else "rebuilt")
}
