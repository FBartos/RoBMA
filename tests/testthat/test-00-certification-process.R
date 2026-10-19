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


.review_certification_fake_process <- function(status = 0L, failure = NULL,
    cleanup_failure = NULL) {

  state <- new.env(parent = emptyenv())
  state$alive <- TRUE
  state$output <- paste0(paste(rep("0123456789", 300L), collapse = ""), "\n")
  state$error <- "retained stderr\n"
  state$expected_output <- state$output
  state$expected_error <- state$error
  state$events <- character()
  state$closed <- FALSE
  state$constructor <- NULL
  read <- function(field, n) {
    value <- substr(state[[field]], 1L, n)
    state[[field]] <- substring(state[[field]], nchar(value) + 1L)
    if (!nzchar(state$output) && !nzchar(state$error)) state$alive <- FALSE
    value
  }
  process <- list(
    is_alive = function() state$alive,
    poll_io = function(timeout) {
      state$events <- c(state$events, paste0("poll:", timeout))
      if (!is.null(failure)) stop(failure)
      c(output = if (nzchar(state$output)) "ready" else "closed",
        error = if (nzchar(state$error)) "ready" else "closed", process = "ready")
    },
    read_output = function(n) read("output", n),
    read_error = function(n) read("error", n),
    kill_tree = function(close_connections = TRUE) {
      state$events <- c(state$events, paste0("kill:", close_connections))
      if (!close_connections && !is.null(cleanup_failure)) stop(cleanup_failure)
      state$alive <- FALSE
      if (close_connections) state$closed <- TRUE
      invisible(integer())
    },
    wait = function(timeout) {
      state$events <- c(state$events, paste0("wait:", timeout))
      invisible(NULL)
    },
    get_exit_status = function() status
  )
  constructor <- list(new = function(...) {
    state$constructor <- list(...)
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
      expect_true(fake$state$closed)
      expect_true(match("kill:FALSE", fake$state$events) < match("wait:1000", fake$state$events))
      expect_identical(tail(fake$state$events, 1L), "kill:TRUE")
      expect_true(fake$state$constructor$cleanup_tree)
      expect_true(fake$state$constructor$windows_hide_window)
      expect_identical(fake$state$constructor$env, Sys.getenv())
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
      caught <- tryCatch(capture.output(helper("owned-test-command", character(), 3)),
        interrupt = identity, error = identity)
      expect_identical(caught, if (kind == "interrupt") interrupt else cleanup_error)
      expect_true(fake$state$closed)
      expect_false(fake$state$alive)
      expect_identical(tail(fake$state$events, 1L), "kill:TRUE")
    })
  }
})
