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


.scenario_is_interactive <- function() {

  return(interactive())
}


.scenario_validate_flag <- function(value, name) {

  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }

  return(value)
}


.scenario_validate_name <- function(name, type = "artifact") {

  valid <- is.character(name) && length(name) == 1L && !is.na(name) &&
    grepl("^[a-z0-9][a-z0-9._-]*$", name) && !endsWith(name, ".")
  if (!valid) {
    stop(
      "Scenario ", type, " names must use lowercase letters, numbers, ",
      "underscores, hyphens, and internal periods.",
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


.scenario_fit_paths <- function(name) {

  config <- .scenario_config()
  name   <- .scenario_validate_name(name, "fit")

  return(list(
    fit  = file.path(config[["cache_dir"]], paste0(name, ".rds")),
    call = file.path(config[["cache_dir"]], paste0(name, ".call.txt"))
  ))
}


# Load an existing cached fit only when its fitting call still matches.
# Otherwise evaluate the expression and replace the fit and call cache.
scenario_fit <- function(name, code) {

  config   <- .scenario_config()
  paths    <- .scenario_fit_paths(name)
  expr     <- substitute(code)
  fit_call <- deparse(expr, width.cutoff = 500L)

  cached_call <- if (file.exists(paths[["call"]])) {
    tryCatch(
      readLines(paths[["call"]], warn = FALSE, encoding = "UTF-8"),
      error = function(error) NULL
    )
  } else {
    NULL
  }
  call_matches <- identical(cached_call, fit_call)

  if (file.exists(paths[["fit"]]) && !config[["regenerate"]] && call_matches) {
    return(tryCatch(
      readRDS(paths[["fit"]]),
      error = function(error) {
        stop(
          "Cached scenario fit '", name,
          "' cannot be read. Regenerate the scenario cache.",
          call. = FALSE
        )
      }
    ))
  }
  if (file.exists(paths[["fit"]]) && !config[["regenerate"]]) {
    reason <- if (is.null(cached_call)) {
      "Missing cached call for"
    } else {
      "Fit call changed for"
    }
    message(
      reason, " scenario fit '", name, "'; refitting."
    )
  }

  fit <- eval(expr, envir = parent.frame())
  dir.create(config[["cache_dir"]], recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, paths[["fit"]])
  .scenario_write_lines(fit_call, paths[["call"]])

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
  if (.scenario_is_interactive() && length(output) > 0L) {
    writeLines(output)
  }

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


.scenario_evaluate_plot <- function(expr, envir, restore_par = TRUE) {

  if (restore_par) {
    old_par <- graphics::par(no.readonly = TRUE)
    old_par[["new"]] <- NULL
    on.exit(graphics::par(old_par), add = TRUE)
  }

  value <- eval(expr, envir = envir)
  if (is.function(value)) {
    value <- value()
  }
  if (inherits(value, "ggplot")) {
    print(value)
  }

  return(invisible(value))
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

    # vdiffr closes its temporary device, so its par settings cannot leak.
    .scenario_evaluate_plot(expr, eval_env, restore_par = FALSE)
  }
  if (.scenario_is_interactive()) {
    .scenario_evaluate_plot(expr, eval_env)
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
