.robma_user_tests_root <- local({

  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    ofile <- frame[["ofile"]]
    if (is.null(ofile)) "" else ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  editor_path <- ""
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()[["path"]],
      error = function(error) ""
    )
  }

  paths <- c(editor_path, rev(ofiles))
  paths <- paths[nzchar(paths)]
  candidates <- c(
    dirname(normalizePath(paths, winslash = "/", mustWork = FALSE)),
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  )
  candidates <- unique(candidates)
  for (candidate in candidates) {
    repeat {
      description <- file.path(candidate, "DESCRIPTION")
      if (file.exists(description)) {
        package <- tryCatch(
          read.dcf(description, fields = "Package")[[1L]],
          error = function(error) ""
        )
        if (identical(package, "RoBMA")) {
          return(candidate)
        }
      }

      parent <- dirname(candidate)
      if (identical(parent, candidate)) {
        break
      }
      candidate <- parent
    }
  }

  stop(
    "Cannot locate the RoBMA package root from the active editor or working directory.",
    call. = FALSE
  )
})

source(
  file.path(.robma_user_tests_root, ".dev", "test-tests.R"),
  local = TRUE
)

.robma_user_test_results <- test_tests(
  refit           = TRUE,
  stop_on_failure = TRUE
)
