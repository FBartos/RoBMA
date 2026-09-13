if (interactive()) {

  library(devtools)
  library(testthat)
  library(vdiffr)

  # R restores .RData after .Rprofile. Reload development helpers afterward
  # so saved functions and absolute paths cannot override this checkout.
  .First <- local({
    project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    previous_first <- get0(".First", envir = .GlobalEnv, inherits = FALSE)

    function() {
      if (is.null(previous_first)) {
        rm(".First", envir = .GlobalEnv)
      } else {
        assign(".First", previous_first, envir = .GlobalEnv)
        previous_first()
      }
      source(file.path(project_root, ".dev", "test-tests.R"), local = .GlobalEnv)
    }
  })
}
