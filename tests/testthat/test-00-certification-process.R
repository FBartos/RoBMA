.review_certification_process <- function() {

  path <- testthat::test_path("..", "..", "tools", "test-profile.R")
  skip_if_not(file.exists(path), "Developer profile tool requires a source checkout.")
  definitions <- Filter(function(expr) is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
    identical(expr[[2L]], as.name(".run_certification_process")), parse(path))
  stopifnot(length(definitions) == 1L)
  environment <- new.env(parent = parent.frame())
  eval(definitions[[1L]], environment)
  environment$.run_certification_process
}


# A process that writes one piece of its output to the file of its stream on every wait() and
# ends after the last piece, the way a child redirected to files grows them while it runs.
.review_certification_fake_process <- function(status = 0L, failure = NULL,
    cleanup_failure = NULL) {

  state <- new.env(parent = emptyenv())
  state$alive <- TRUE
  state$pieces <- list(
    list(stream = "stdout", text = strrep("0123456789", 100L)),
    list(stream = "stdout", text = strrep("0123456789", 100L)),
    list(stream = "stderr", text = "retained stderr\n"),
    list(stream = "stdout", text = paste0(strrep("0123456789", 100L), "\n")))
  state$expected_output <- paste0(
    vapply(Filter(function(piece) piece$stream == "stdout", state$pieces),
           function(piece) piece$text, character(1L)), collapse = "")
  state$expected_error <- "retained stderr\n"
  state$events <- character()
  state$files_existed <- logical()
  state$constructor <- NULL
  append_bytes <- function(path, text) {
    connection <- file(path, open = "ab")
    writeBin(charToRaw(text), connection)
    close(connection)
  }
  process <- list(
    is_alive = function() state$alive,
    wait = function(timeout) {
      state$events <- c(state$events, paste0("wait:", timeout))
      if (!is.null(failure)) stop(failure)
      if (length(state$pieces) > 0L) {
        piece <- state$pieces[[1L]]
        state$pieces <- state$pieces[-1L]
        append_bytes(state$constructor[[piece$stream]], piece$text)
        if (length(state$pieces) == 0L) state$alive <- FALSE
      }
      invisible(NULL)
    },
    kill_tree = function(close_connections = TRUE) {
      state$events <- c(state$events, paste0("kill:", close_connections))
      if (!close_connections && !is.null(cleanup_failure)) stop(cleanup_failure)
      state$alive <- FALSE
      invisible(integer())
    },
    get_exit_status = function() status
  )
  constructor <- list(new = function(...) {
    state$constructor <- list(...)
    state$files_existed <- file.exists(c(state$constructor$stdout, state$constructor$stderr))
    process
  })
  list(state = state, constructor = constructor)
}


test_that("certification subprocesses retain complete success and nonzero output", {

  helper <- .review_certification_process()
  for (status in c(0L, 2L)) {
    local({
      fake <- .review_certification_fake_process(status)
      testthat::local_mocked_bindings(process = fake$constructor, .package = "processx")
      message_output <- capture.output(
        output <- capture.output(result <- helper("owned-test-command", "owned-test-arg", 3)),
        type = "message")
      expect_identical(result$status, status)
      expect_false(result$timeout)
      expect_identical(result$stdout, fake$state$expected_output)
      expect_identical(result$stderr, fake$state$expected_error)
      expect_identical(paste0(output, collapse = "\n"), sub("\n$", "", result$stdout))
      expect_identical(paste0(message_output, collapse = "\n"), sub("\n$", "", result$stderr))
      # One wait per piece, then the tree is ended before the final wait and the closing cleanup.
      expect_identical(fake$state$events,
        c(rep("wait:1000", 4L), "kill:FALSE", "wait:1000", "kill:TRUE"))
      expect_true(fake$state$constructor$cleanup_tree)
      expect_true(fake$state$constructor$windows_hide_window)
      # The child writes to two distinct files, which exist while it runs and are removed after.
      files <- c(fake$state$constructor$stdout, fake$state$constructor$stderr)
      expect_false(anyDuplicated(files) > 0L)
      expect_true(all(fake$state$files_existed))
      expect_false(any(file.exists(files)))
    })
  }
})


test_that("certification subprocess interruptions and cleanup failures unwind", {

  helper <- .review_certification_process()
  interrupt <- structure(list(message = "owned fixture interrupt", call = NULL),
    class = c("interrupt", "condition"))
  cleanup_error <- structure(list(message = "owned cleanup failure", call = NULL),
    class = c("owned_cleanup_error", "error", "condition"))
  for (kind in c("interrupt", "cleanup")) {
    local({
      fake <- .review_certification_fake_process(
        failure = if (kind == "interrupt") interrupt else NULL,
        cleanup_failure = if (kind == "cleanup") cleanup_error else NULL)
      testthat::local_mocked_bindings(process = fake$constructor, .package = "processx")
      caught <- tryCatch(
        capture.output(capture.output(helper("owned-test-command", character(), 3)),
          type = "message"),
        interrupt = identity, error = identity)
      expect_identical(caught, if (kind == "interrupt") interrupt else cleanup_error)
      expect_false(fake$state$alive)
      expect_identical(tail(fake$state$events, 1L), "kill:TRUE")
      expect_false(any(file.exists(
        c(fake$state$constructor$stdout, fake$state$constructor$stderr))))
    })
  }
})


test_that("a real certification child returns its output and exit status through files", {

  skip_on_cran()
  skip_if_not_installed("processx")
  helper  <- .review_certification_process()
  rscript <- file.path(R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  child <- tempfile("owned-certification-child-", fileext = ".R")
  on.exit(unlink(child), add = TRUE)
  writeLines(c(
    "writeLines(sprintf('line %d', 1:20000))",
    "message('a message on stderr')",
    "quit(save = 'no', status = 3L)"
  ), child)

  message_output <- capture.output(
    output <- capture.output(result <- helper(rscript, c("--vanilla", child), 60)),
    type = "message")

  # Windows text-mode files carry "\r\n".
  lines <- sprintf("line %d", 1:20000)
  expect_identical(result$status, 3L)
  expect_false(result$timeout)
  expect_identical(gsub("\r", "", result$stdout, fixed = TRUE),
    paste0(paste(lines, collapse = "\n"), "\n"))
  expect_match(result$stderr, "a message on stderr", fixed = TRUE)
  # The command echoed by processx precedes the child's lines.
  expect_identical(gsub("\r", "", tail(output, length(lines)), fixed = TRUE), lines)
  expect_true(any(grepl("a message on stderr", message_output, fixed = TRUE)))
})


test_that("a real certification child past its limit is killed with its grandchild", {

  skip_on_cran()
  skip_if_not_installed("processx")
  helper  <- .review_certification_process()
  rscript <- file.path(R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  directory <- normalizePath(tempfile("owned-certification-child-"),
    winslash = "/", mustWork = FALSE)
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  pid_file   <- file.path(directory, "grandchild.pid")
  grandchild <- file.path(directory, "grandchild.R")
  child      <- file.path(directory, "child.R")
  writeLines(c(
    sprintf('writeLines(as.character(Sys.getpid()), "%s")', pid_file),
    "Sys.sleep(120)"
  ), grandchild)
  writeLines(c(
    "writeLines('child started')",
    "flush(stdout())",
    sprintf('system2("%s", c("--vanilla", "%s"), wait = FALSE)', rscript, grandchild),
    "Sys.sleep(120)"
  ), child)

  started <- proc.time()[["elapsed"]]
  output <- capture.output(result <- helper(rscript, c("--vanilla", child), 3))
  elapsed <- proc.time()[["elapsed"]] - started

  expect_true(result$timeout)
  expect_false(identical(result$status, 0L))
  expect_lt(elapsed, 10)
  expect_match(result$stdout, "child started", fixed = TRUE)
  expect_true(any(grepl("child started", output, fixed = TRUE)))
  expect_true(file.exists(pid_file), info = "The grandchild started before the limit.")
  if (file.exists(pid_file)) {
    pid <- as.integer(readLines(pid_file, warn = FALSE))
    expect_false(pid %in% ps::ps_pids())
  }
})
