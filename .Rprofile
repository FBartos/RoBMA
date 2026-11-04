if (interactive()) {

  library(devtools)
  library(testthat)
  library(vdiffr)

  # R restores .RData after .Rprofile. An active binding captures any saved
  # .First while keeping this wrapper in place until startup invokes it.
  local({
    project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    state <- new.env(parent = emptyenv())
    state$had_previous <- exists(".First", envir = .GlobalEnv, inherits = FALSE)
    if (state$had_previous) {
      state$previous <- get(".First", envir = .GlobalEnv, inherits = FALSE)
      rm(".First", envir = .GlobalEnv)
    }

    first <- function() {
      had_previous <- state$had_previous
      previous     <- if (had_previous) state$previous else NULL
      rm(".First", envir = .GlobalEnv)
      if (had_previous) {
        assign(".First", previous, envir = .GlobalEnv)
        if (is.function(previous)) {
          previous()
        }
      }
      source(file.path(project_root, ".dev", "test-tests.R"), local = .GlobalEnv)
    }

    makeActiveBinding(".First", function(value) {
      if (missing(value)) {
        return(first)
      }
      state$had_previous <- TRUE
      state$previous     <- value
      invisible(NULL)
    }, .GlobalEnv)
  })
}
