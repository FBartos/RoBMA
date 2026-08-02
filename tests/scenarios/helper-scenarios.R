.scenario_helpers_dir <- local({

  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    ofile <- frame[["ofile"]]
    if (is.null(ofile)) "" else ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  candidates <- c(
    dirname(ofiles),
    getwd(),
    file.path(getwd(), "tests", "scenarios"),
    file.path(getwd(), "scenarios")
  )
  candidates <- unique(normalizePath(
    candidates,
    winslash = "/",
    mustWork = FALSE
  ))
  available <- vapply(candidates, function(candidate) {
    file.exists(file.path(candidate, "helper-scenarios.R"))
  }, logical(1))
  if (!any(available)) {
    stop("Cannot locate tests/scenarios/helper-scenarios.R.", call. = FALSE)
  }

  candidates[which(available)[1L]]
})

.scenario_visual_writer <- file.path(
  dirname(.scenario_helpers_dir),
  "testthat",
  "helper-visual-writer.R"
)
if (!exists(".write_canonical_svg", mode = "function")) {
  source(.scenario_visual_writer, local = TRUE)
}

.scenario_state <- new.env(parent = emptyenv())


.scenario_is_true <- function(value) {

  if (length(value) != 1L || is.na(value)) {
    return(FALSE)
  }

  return(tolower(as.character(value)) %in% c("1", "true", "yes", "y"))
}


.scenario_on_ci <- function() {

  return(
    .scenario_is_true(Sys.getenv("CI")) ||
      .scenario_is_true(Sys.getenv("GITHUB_ACTIONS"))
  )
}


.scenario_validate_flag <- function(value, name) {

  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }

  return(value)
}


.scenario_validate_name <- function(name, type = "artifact") {

  valid <- is.character(name) && length(name) == 1L && !is.na(name) &&
    grepl("^[a-z0-9][a-z0-9_-]*$", name)
  if (!valid) {
    stop(
      "Scenario ", type, " names must use lowercase letters, numbers, ",
      "underscores, and hyphens.",
      call. = FALSE
    )
  }

  return(name)
}


.scenario_config <- function() {

  if (!exists("config", envir = .scenario_state, inherits = FALSE)) {
    stop("Call scenario_start() before using scenario helpers.", call. = FALSE)
  }

  return(get("config", envir = .scenario_state, inherits = FALSE))
}


# Start one scenario. An environment variable set by tools/test-scenario.R can
# force regeneration even when the scenario file keeps its local flag FALSE.
scenario_start <- function(name, regenerate = NULL,
                           root = .scenario_helpers_dir,
                           create_missing = NULL, width = 150L) {

  name <- .scenario_validate_name(name, "scenario")

  caller_regenerate <- get0(
    "REGENERATE_SCENARIO_FILES",
    envir      = parent.frame(),
    inherits   = TRUE,
    ifnotfound = FALSE
  )
  caller_regenerate <- .scenario_validate_flag(
    caller_regenerate,
    "REGENERATE_SCENARIO_FILES"
  )
  env_regenerate <- .scenario_is_true(
    Sys.getenv("ROBMA_SCENARIO_REGENERATE")
  )

  if (is.null(regenerate)) {
    regenerate <- caller_regenerate
  }
  regenerate <- .scenario_validate_flag(regenerate, "regenerate")
  regenerate <- regenerate || env_regenerate

  if (is.null(create_missing)) {
    create_missing <- !.scenario_on_ci() || regenerate
  }
  create_missing <- .scenario_validate_flag(create_missing, "create_missing")

  if (!is.numeric(width) || length(width) != 1L || is.na(width) ||
      width < 20 || width != as.integer(width)) {
    stop("'width' must be one integer of at least 20.", call. = FALSE)
  }

  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  config <- list(
    name           = name,
    root           = root,
    cache_dir      = file.path(root, "cache", name),
    results_dir    = file.path(root, "results", name),
    regenerate     = regenerate,
    create_missing = create_missing,
    width          = as.integer(width)
  )
  assign("config", config, envir = .scenario_state)

  if (regenerate) {
    message("Regenerating fits and locked outputs for scenario '", name, "'.")
  }

  return(invisible(config))
}


.scenario_fit_path <- function(name) {

  config <- .scenario_config()
  name   <- .scenario_validate_name(name, "fit")

  return(file.path(config[["cache_dir"]], paste0(name, ".rds")))
}


# Load an existing cached fit without evaluating code, or evaluate and cache
# the fitting expression when the cache is missing or regeneration is enabled.
scenario_fit <- function(name, code) {

  config <- .scenario_config()
  path   <- .scenario_fit_path(name)
  expr   <- substitute(code)

  if (file.exists(path) && !config[["regenerate"]]) {
    return(tryCatch(
      readRDS(path),
      error = function(error) {
        stop(
          "Cached scenario fit '", name,
          "' cannot be read. Regenerate the scenario cache.",
          call. = FALSE
        )
      }
    ))
  }

  fit <- eval(expr, envir = parent.frame())
  dir.create(config[["cache_dir"]], recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, path)

  return(fit)
}


.scenario_write_lines <- function(lines, path) {

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(lines, connection, useBytes = TRUE)

  return(invisible(path))
}


# Capture a deliberately printed result and lock it to a tracked text file.
scenario_text <- function(name, code) {

  config <- .scenario_config()
  name   <- .scenario_validate_name(name, "text")
  path   <- file.path(config[["results_dir"]], paste0(name, ".txt"))
  expr   <- substitute(code)
  value  <- NULL

  old_options <- options(width = config[["width"]])
  on.exit(options(old_options), add = TRUE)
  output <- capture.output(
    value <- eval(expr, envir = parent.frame()),
    type = "output"
  )

  if (!file.exists(path) && !config[["create_missing"]] &&
      !config[["regenerate"]]) {
    testthat::fail(paste0(
      "Missing locked scenario text 'results/", config[["name"]], "/",
      basename(path), "'. Create and review it locally."
    ))
    return(invisible(value))
  }

  if (config[["regenerate"]] || !file.exists(path)) {
    .scenario_write_lines(output, path)
    message(
      if (config[["regenerate"]]) "Regenerated" else "Added",
      " scenario text: ", path
    )
  } else {
    expected <- readLines(path, warn = FALSE, encoding = "UTF-8")
    testthat::expect_equal(
      output,
      expected,
      info = paste0("Scenario text '", config[["name"]], "/", name, "' changed.")
    )
  }

  return(invisible(value))
}


.scenario_snapshot_context <- function() {

  get_snapshotter <- get0(
    "get_snapshotter",
    envir    = asNamespace("testthat"),
    inherits = FALSE
  )
  snapshotter <- get_snapshotter()
  if (is.null(snapshotter) || is.null(snapshotter$file)) {
    stop(
      "scenario_plot() must run through tools/test-scenario.R or testthat::test_file().",
      call. = FALSE
    )
  }

  return(snapshotter)
}


# Render a plot expression on vdiffr's device and compare it with a tracked SVG.
scenario_plot <- function(name, code) {

  testthat::skip_if_not_installed("vdiffr")

  config      <- .scenario_config()
  name        <- .scenario_validate_name(name, "plot")
  snapshotter <- .scenario_snapshot_context()
  if (!identical(snapshotter$file, config[["name"]])) {
    stop(
      "Scenario '", config[["name"]], "' must be stored in 'test-",
      config[["name"]], ".R' so vdiffr uses the expected snapshot directory.",
      call. = FALSE
    )
  }

  expected <- file.path(
    snapshotter$snap_dir,
    config[["name"]],
    paste0(name, ".svg")
  )
  if (!file.exists(expected) && !config[["create_missing"]] &&
      !config[["regenerate"]]) {
    testthat::fail(paste0(
      "Missing locked scenario plot '_snaps/", config[["name"]], "/",
      basename(expected), "'. Create and review it locally."
    ))
    return(invisible(NULL))
  }

  expr     <- substitute(code)
  eval_env <- parent.frame()
  figure <- function() {

    value <- eval(expr, envir = eval_env)
    if (is.function(value)) {
      value <- value()
    }
    if (inherits(value, "ggplot")) {
      print(value)
    }
    return(invisible(value))
  }

  writer <- function(plot, file, title = "") {

    .write_canonical_svg(plot, file, title)
    if (config[["regenerate"]]) {
      dir.create(dirname(expected), recursive = TRUE, showWarnings = FALSE)
      copied <- file.copy(file, expected, overwrite = TRUE)
      if (!copied) {
        stop("Failed to regenerate scenario plot: ", expected, call. = FALSE)
      }
    }
  }

  vdiffr::expect_doppelganger(
    title  = name,
    fig    = figure,
    writer = writer
  )

  return(invisible(NULL))
}
