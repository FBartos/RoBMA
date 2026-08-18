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
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", name) && !endsWith(name, ".")
  if (!valid) {
    stop(
      "Scenario ", type, " names must use letters, numbers, ",
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


# Start one scenario. Direct interactive execution shows and updates output;
# test_scenario(), tools/test-scenario.R, and non-interactive execution compare
# quietly unless their explicit controls request otherwise.
scenario_start <- function(name, regenerate = NULL, refit = NULL, update = NULL,
                           update_timings = NULL, show_output = NULL,
                           root = .scenario_helpers_dir,
                           create_missing = NULL, width = 150L) {

  name <- .scenario_validate_name(name, "scenario")
  direct_interactive <- .scenario_is_interactive() &&
    !.scenario_is_true(Sys.getenv("ROBMA_SCENARIO_RUNNER"))
  env_regenerate <- .scenario_is_true(
    Sys.getenv("ROBMA_SCENARIO_REGENERATE")
  )

  if (is.null(regenerate)) {
    regenerate <- FALSE
  }
  regenerate <- .scenario_validate_flag(regenerate, "regenerate")
  regenerate <- regenerate || env_regenerate

  env_refit <- .scenario_is_true(Sys.getenv("ROBMA_SCENARIO_REFIT"))
  if (is.null(refit)) {
    refit <- FALSE
  }
  refit <- .scenario_validate_flag(refit, "refit")
  refit <- refit || env_refit || regenerate

  env_update <- .scenario_is_true(Sys.getenv("ROBMA_SCENARIO_UPDATE"))
  if (is.null(update)) {
    update <- direct_interactive
  }
  update <- .scenario_validate_flag(update, "update")
  update <- update || env_update || regenerate

  env_update_timings <- .scenario_is_true(
    Sys.getenv("ROBMA_SCENARIO_UPDATE_TIMINGS")
  )
  if (is.null(update_timings)) {
    update_timings <- direct_interactive
  }
  update_timings <- .scenario_validate_flag(update_timings, "update_timings")
  update_timings <- update_timings || env_update_timings || update

  if (is.null(show_output)) {
    show_output <- direct_interactive
  }
  show_output <- .scenario_validate_flag(show_output, "show_output")

  if (is.null(create_missing)) {
    create_missing <- direct_interactive || update
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
    snapshots_dir  = file.path(root, "_snaps", name),
    timings_dir    = file.path(root, "timings"),
    regenerate     = regenerate,
    refit          = refit,
    update         = update,
    update_timings = update_timings,
    show_output    = show_output,
    create_missing = create_missing,
    direct         = direct_interactive,
    width          = as.integer(width)
  )
  assign("config", config, envir = .scenario_state)
  assign(
    "timing_baseline",
    .scenario_read_timings(.scenario_timing_paths()[["expected"]]),
    envir = .scenario_state
  )
  assign("timings", .scenario_empty_timings(), envir = .scenario_state)
  assign(
    "timing_unavailable",
    .scenario_empty_timing_unavailable(),
    envir = .scenario_state
  )
  assign("timing_warned", character(), envir = .scenario_state)

  if (refit && update) {
    message("Refitting models and updating locked outputs for scenario '",
            name, "'.")
  } else if (refit) {
    message("Refitting models for scenario '", name, "'.")
  } else if (update) {
    message("Updating locked outputs for scenario '", name, "'.")
  }

  return(invisible(config))
}


.scenario_fit_paths <- function(name) {

  config <- .scenario_config()
  name   <- .scenario_validate_name(name, "fit")

  return(list(
    fit    = file.path(config[["cache_dir"]], paste0(name, ".rds")),
    call   = file.path(config[["cache_dir"]], paste0(name, ".call.txt")),
    timing = file.path(config[["cache_dir"]], paste0(name, ".timing.rds"))
  ))
}


# Load an existing cached fit only when its fitting call still matches.
# Otherwise evaluate the expression and replace the fit and call cache.
scenario_fit <- function(name, code, cache_version = NULL) {

  config   <- .scenario_config()
  paths    <- .scenario_fit_paths(name)
  expr     <- substitute(code)
  fit_call <- deparse(expr, width.cutoff = 500L)
  if (!is.null(cache_version)) {
    if (!is.numeric(cache_version) || length(cache_version) != 1L ||
        is.na(cache_version) || cache_version < 1L ||
        cache_version != as.integer(cache_version)) {
      stop("'cache_version' must be NULL or one positive integer.",
           call. = FALSE)
    }
    fit_call <- c(
      paste0("cache_version: ", as.integer(cache_version)),
      fit_call
    )
  }

  cached_call <- if (file.exists(paths[["call"]])) {
    tryCatch(
      readLines(paths[["call"]], warn = FALSE, encoding = "UTF-8"),
      error = function(error) NULL
    )
  } else {
    NULL
  }
  call_matches <- identical(cached_call, fit_call)

  if (file.exists(paths[["fit"]]) && !config[["refit"]] && call_matches) {
    fit <- tryCatch(
      readRDS(paths[["fit"]]),
      error = function(error) {
        stop(
          "Cached scenario fit '", name,
          "' cannot be read. Regenerate the scenario cache.",
          call. = FALSE
        )
      }
    )
    fit_timing <- if (file.exists(paths[["timing"]])) {
      tryCatch(
        suppressWarnings(readRDS(paths[["timing"]])),
        error = function(error) NULL
      )
    } else {
      NULL
    }
    valid_timing <- is.list(fit_timing) &&
      identical(fit_timing[["version"]], 1L) &&
      identical(fit_timing[["fit_call"]], fit_call) &&
      is.numeric(fit_timing[["elapsed"]]) &&
      length(fit_timing[["elapsed"]]) == 1L &&
      is.finite(fit_timing[["elapsed"]]) &&
      fit_timing[["elapsed"]] >= 0
    if (valid_timing) {
      .scenario_register_timing(
        "fit", name, fit_timing[["elapsed"]], fit_timing[["provenance"]]
      )
    } else {
      .scenario_mark_timing_unavailable(
        "fit", name,
        "the cached fit predates valid timing metadata; rerun with 'refit = TRUE'"
      )
    }
    return(fit)
  }
  if (file.exists(paths[["fit"]]) && !config[["refit"]]) {
    reason <- if (is.null(cached_call)) {
      "Missing cached call for"
    } else {
      "Fit call changed for"
    }
    message(
      reason, " scenario fit '", name, "'; refitting."
    )
  }

  started <- .scenario_clock()
  fit     <- eval(expr, envir = parent.frame())
  elapsed <- max(0, .scenario_clock() - started)
  provenance <- .scenario_timing_provenance()
  dir.create(config[["cache_dir"]], recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, paths[["fit"]])
  saveRDS(
    list(
      version    = 1L,
      fit_call   = fit_call,
      elapsed    = elapsed,
      provenance = provenance
    ),
    paths[["timing"]]
  )
  .scenario_write_lines(fit_call, paths[["call"]])
  .scenario_register_timing("fit", name, elapsed, provenance)

  return(fit)
}


.scenario_write_lines <- function(lines, path) {

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(lines, connection, useBytes = TRUE)

  return(invisible(path))
}


.scenario_candidate_path <- function(path) {

  return(sub("(\\.[^.]+)$", ".new\\1", path))
}


.scenario_clock <- function() {

  return(unname(proc.time()[["elapsed"]]))
}


.scenario_timing_provenance <- function() {

  system <- Sys.info()[["sysname"]]
  if (is.null(system) || is.na(system)) {
    system <- "unknown"
  }

  return(list(
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform  = paste(system, R.version$platform, sep = "/")
  ))
}


.scenario_empty_timings <- function() {

  return(data.frame(
    type      = character(),
    name      = character(),
    elapsed   = numeric(),
    r_version = character(),
    platform  = character(),
    stringsAsFactors = FALSE
  ))
}


.scenario_empty_timing_unavailable <- function() {

  return(data.frame(
    type   = character(),
    name   = character(),
    reason = character(),
    stringsAsFactors = FALSE
  ))
}


.scenario_timing_key <- function(timings) {

  return(paste(timings[["type"]], timings[["name"]], sep = "/"))
}


.scenario_order_timings <- function(timings) {

  type_order <- match(timings[["type"]], c("fit", "text", "plot"))
  return(timings[order(type_order, timings[["name"]]), , drop = FALSE])
}


.scenario_timing_paths <- function() {

  config   <- .scenario_config()
  expected <- file.path(config[["timings_dir"]], paste0(config[["name"]], ".tsv"))

  return(list(
    expected  = expected,
    candidate = .scenario_candidate_path(expected)
  ))
}


.scenario_read_timings <- function(path) {

  if (!file.exists(path)) {
    return(.scenario_empty_timings())
  }
  timings <- tryCatch(
    utils::read.delim(
      path,
      header           = TRUE,
      stringsAsFactors = FALSE,
      check.names      = FALSE
    ),
    error = function(error) {
      stop("Cannot read scenario timing baseline: ", path, ".", call. = FALSE)
    }
  )
  required <- names(.scenario_empty_timings())
  if (!identical(names(timings), required)) {
    stop("Scenario timing baseline has invalid columns: ", path, ".",
         call. = FALSE)
  }

  timings[["type"]]      <- as.character(timings[["type"]])
  timings[["name"]]      <- as.character(timings[["name"]])
  timings[["elapsed"]]   <- as.numeric(timings[["elapsed"]])
  timings[["r_version"]] <- as.character(timings[["r_version"]])
  timings[["platform"]]  <- as.character(timings[["platform"]])
  valid <- !is.na(timings[["type"]]) & !is.na(timings[["name"]]) &
    timings[["type"]] %in% c("fit", "text", "plot") &
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", timings[["name"]]) &
    !endsWith(timings[["name"]], ".") &
    is.finite(timings[["elapsed"]]) & timings[["elapsed"]] >= 0 &
    !is.na(timings[["r_version"]]) & nzchar(timings[["r_version"]]) &
    !is.na(timings[["platform"]]) & nzchar(timings[["platform"]])
  if (any(!valid) || anyDuplicated(.scenario_timing_key(timings))) {
    stop("Scenario timing baseline has invalid rows: ", path, ".",
         call. = FALSE)
  }

  rownames(timings) <- NULL
  return(.scenario_order_timings(timings))
}


.scenario_write_timings <- function(timings, path) {

  timings <- .scenario_order_timings(timings)
  lines <- c(
    paste(names(timings), collapse = "\t"),
    if (nrow(timings) > 0L) {
      apply(timings, 1L, function(row) {
        row[["elapsed"]] <- sprintf("%.9f", as.numeric(row[["elapsed"]]))
        paste(row, collapse = "\t")
      })
    }
  )
  .scenario_write_lines(lines, path)
  return(invisible(path))
}


.scenario_current_timings <- function() {

  if (!exists("timings", envir = .scenario_state, inherits = FALSE)) {
    return(.scenario_empty_timings())
  }

  return(get("timings", envir = .scenario_state, inherits = FALSE))
}


.scenario_timing_issues <- function(final = TRUE) {

  current     <- .scenario_current_timings()
  baseline    <- get("timing_baseline", envir = .scenario_state, inherits = FALSE)
  unavailable <- get(
    "timing_unavailable", envir = .scenario_state, inherits = FALSE
  )
  issues <- data.frame(id = character(), text = character(), stringsAsFactors = FALSE)
  add_issue <- function(id, text) {

    issues <<- rbind(
      issues,
      data.frame(id = id, text = text, stringsAsFactors = FALSE)
    )
  }

  baseline_key <- .scenario_timing_key(baseline)
  current_key  <- .scenario_timing_key(current)
  matched      <- match(current_key, baseline_key)
  comparable   <- !is.na(matched)
  regressed <- comparable &
    current[["elapsed"]] > 1.20 * baseline[["elapsed"]][matched]
  for (i in which(regressed)) {
    old <- baseline[["elapsed"]][matched[[i]]]
    new <- current[["elapsed"]][[i]]
    increase <- if (old == 0) Inf else 100 * (new / old - 1)
    add_issue(
      paste0("call:", current_key[[i]]),
      paste0(
        current_key[[i]], ": ", sprintf("%.3f", new), " s vs ",
        sprintf("%.3f", old), " s (+", sprintf("%.1f", increase),
        "%; threshold 20%)"
      )
    )
  }

  if (final && nrow(baseline) == 0L) {
    add_issue(
      "baseline:missing",
      paste0(
        "timing baseline is unavailable; review '",
        chartr("\\", "/", .scenario_timing_paths()[["candidate"]]),
        "' and update timings"
      )
    )
  }
  if (final && nrow(unavailable) > 0L) {
    for (i in seq_len(nrow(unavailable))) {
      key <- paste(unavailable[["type"]][[i]], unavailable[["name"]][[i]],
                   sep = "/")
      add_issue(
        paste0("unavailable:", key),
        paste0(key, " timing is unavailable because ",
               unavailable[["reason"]][[i]], ".")
      )
    }
  }

  observed_key <- c(current_key, .scenario_timing_key(unavailable))
  if (final && nrow(baseline) > 0L) {
    added   <- setdiff(observed_key, baseline_key)
    missing <- setdiff(baseline_key, observed_key)
    if (length(added) > 0L || length(missing) > 0L) {
      parts <- character()
      if (length(added) > 0L) {
        parts <- c(parts, paste0("added: ", paste(added, collapse = ", ")))
      }
      if (length(missing) > 0L) {
        parts <- c(parts, paste0("missing: ", paste(missing, collapse = ", ")))
      }
      add_issue(
        "structure",
        paste0("timing entries changed (", paste(parts, collapse = "; "), ")")
      )
    }
  }

  same_entries <- nrow(baseline) > 0L && nrow(unavailable) == 0L &&
    setequal(current_key, baseline_key)
  if (final && same_entries) {
    baseline_elapsed <- baseline[["elapsed"]][match(current_key, baseline_key)]
    percentage_change <- ifelse(
      baseline_elapsed == 0,
      ifelse(current[["elapsed"]] == 0, 0, Inf),
      100 * (current[["elapsed"]] / baseline_elapsed - 1)
    )
    average_change <- mean(percentage_change)
    if (average_change > 5) {
      add_issue(
        "average",
        paste0(
          "average timing regression: ", sprintf("%.1f", average_change),
          "% (unweighted mean across ", length(percentage_change),
          " calls; threshold 5%)"
        )
      )
    }
  }

  if (nrow(issues) > 0L) {
    issues <- issues[!duplicated(issues[["id"]]), , drop = FALSE]
  }
  rownames(issues) <- NULL
  return(issues)
}


.scenario_warn_timing_issues <- function(issues) {

  if (nrow(issues) == 0L) {
    return(invisible(issues))
  }
  warned <- get("timing_warned", envir = .scenario_state, inherits = FALSE)
  issues <- issues[!issues[["id"]] %in% warned, , drop = FALSE]
  if (nrow(issues) == 0L) {
    return(invisible(issues))
  }
  assign(
    "timing_warned", c(warned, issues[["id"]]), envir = .scenario_state
  )
  warning(
    "Scenario timing warnings for '", .scenario_config()[["name"]],
    "':\n- ", paste(issues[["text"]], collapse = "\n- "),
    call. = FALSE
  )

  return(invisible(issues))
}


.scenario_merge_timings <- function(baseline, current) {

  if (nrow(current) == 0L) {
    return(baseline)
  }
  keep <- !.scenario_timing_key(baseline) %in% .scenario_timing_key(current)
  return(.scenario_order_timings(rbind(
    baseline[keep, , drop = FALSE], current
  )))
}


.scenario_register_timing <- function(type, name, elapsed,
                                      provenance = .scenario_timing_provenance()) {

  type <- match.arg(type, c("fit", "text", "plot"))
  name <- .scenario_validate_name(name, type)
  if (!is.numeric(elapsed) || length(elapsed) != 1L ||
      !is.finite(elapsed) || elapsed < 0) {
    stop("Scenario elapsed time must be one finite non-negative number.",
         call. = FALSE)
  }
  valid_provenance <- is.list(provenance) &&
    is.character(provenance[["r_version"]]) &&
    length(provenance[["r_version"]]) == 1L &&
    !is.na(provenance[["r_version"]]) &&
    nzchar(provenance[["r_version"]]) &&
    is.character(provenance[["platform"]]) &&
    length(provenance[["platform"]]) == 1L &&
    !is.na(provenance[["platform"]]) &&
    nzchar(provenance[["platform"]])
  if (!valid_provenance) {
    provenance <- .scenario_timing_provenance()
  }

  current <- .scenario_current_timings()
  row <- data.frame(
    type      = type,
    name      = name,
    elapsed   = as.numeric(elapsed),
    r_version = provenance[["r_version"]],
    platform  = provenance[["platform"]],
    stringsAsFactors = FALSE
  )
  keep    <- .scenario_timing_key(current) != .scenario_timing_key(row)
  current <- .scenario_order_timings(rbind(current[keep, , drop = FALSE], row))
  assign("timings", current, envir = .scenario_state)
  unavailable <- get(
    "timing_unavailable", envir = .scenario_state, inherits = FALSE
  )
  unavailable <- unavailable[
    .scenario_timing_key(unavailable) != .scenario_timing_key(row),
    , drop = FALSE
  ]
  assign("timing_unavailable", unavailable, envir = .scenario_state)

  config <- .scenario_config()
  if (config[["direct"]]) {
    baseline <- get(
      "timing_baseline", envir = .scenario_state, inherits = FALSE
    )
    complete <- nrow(baseline) > 0L &&
      setequal(.scenario_timing_key(current), .scenario_timing_key(baseline))
    .scenario_warn_timing_issues(.scenario_timing_issues(final = complete))
    paths <- .scenario_timing_paths()
    .scenario_write_timings(current, paths[["candidate"]])
    if (config[["update_timings"]] && nrow(unavailable) == 0L) {
      .scenario_write_timings(
        .scenario_merge_timings(baseline, current), paths[["expected"]]
      )
      unlink(paths[["candidate"]])
    }
  }

  return(invisible(elapsed))
}


.scenario_mark_timing_unavailable <- function(type, name, reason) {

  type <- match.arg(type, c("fit", "text", "plot"))
  name <- .scenario_validate_name(name, type)
  unavailable <- get(
    "timing_unavailable", envir = .scenario_state, inherits = FALSE
  )
  row <- data.frame(
    type = type, name = name, reason = reason, stringsAsFactors = FALSE
  )
  current <- .scenario_current_timings()
  current <- current[
    .scenario_timing_key(current) != .scenario_timing_key(row),
    , drop = FALSE
  ]
  assign("timings", current, envir = .scenario_state)
  keep <- .scenario_timing_key(unavailable) != .scenario_timing_key(row)
  assign(
    "timing_unavailable",
    rbind(unavailable[keep, , drop = FALSE], row),
    envir = .scenario_state
  )
  if (.scenario_config()[["direct"]]) {
    .scenario_warn_timing_issues(.scenario_timing_issues(final = TRUE))
  }

  return(invisible(NULL))
}


.scenario_finalize_timing <- function(scenario = NULL, allow_update = TRUE) {

  if (!exists("config", envir = .scenario_state, inherits = FALSE)) {
    return(invisible(NULL))
  }
  config <- .scenario_config()
  if (!is.null(scenario) && !identical(config[["name"]], scenario)) {
    return(invisible(NULL))
  }
  current <- .scenario_current_timings()
  unavailable <- get(
    "timing_unavailable", envir = .scenario_state, inherits = FALSE
  )
  if (nrow(current) == 0L && nrow(unavailable) == 0L) {
    return(invisible(NULL))
  }

  paths <- .scenario_timing_paths()
  if (nrow(current) > 0L) {
    .scenario_write_timings(current, paths[["candidate"]])
  }
  issues <- .scenario_timing_issues(final = TRUE)
  .scenario_warn_timing_issues(issues)

  can_update <- config[["update_timings"]] && allow_update &&
    nrow(current) > 0L && nrow(unavailable) == 0L
  if (can_update) {
    existed <- file.exists(paths[["expected"]])
    .scenario_write_timings(current, paths[["expected"]])
    unlink(paths[["candidate"]])
    message(
      if (existed) "Updated" else "Added",
      " scenario timings: ", paths[["expected"]]
    )
  } else if (config[["update_timings"]] && !allow_update) {
    warning(
      "Scenario timing baseline for '", config[["name"]],
      "' was not updated because the scenario did not finish.",
      call. = FALSE
    )
  }

  return(invisible(current))
}


.scenario_plot_filename <- function(name) {

  name <- gsub("[^a-z0-9]", "-", tolower(name))
  name <- gsub("--+", "-", name)
  name <- gsub("^-|-$", "", name)

  return(paste0(name, ".svg"))
}


.scenario_copy_file <- function(from, to, type) {

  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(from, to, overwrite = TRUE)
  if (!copied) {
    stop("Failed to write scenario ", type, ": ", to, call. = FALSE)
  }

  return(invisible(to))
}


.scenario_text_candidates <- function(root, scenarios) {

  changes <- lapply(scenarios, function(scenario) {
    candidates <- list.files(
      file.path(root, "results", scenario),
      pattern    = "[.]new[.]txt$",
      full.names = TRUE
    )
    expected   <- sub("[.]new[.]txt$", ".txt", candidates)
    keep       <- file.exists(expected)
    candidates <- candidates[keep]
    expected   <- expected[keep]

    data.frame(
      name = paste0(
        scenario, "/",
        sub("[.]new[.]txt$", "", basename(candidates))
      ),
      expected  = expected,
      candidate = candidates,
      stringsAsFactors = FALSE
    )
  })
  changes <- do.call(rbind, changes)
  if (is.null(changes)) {
    changes <- data.frame(
      name      = character(),
      expected  = character(),
      candidate = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(changes) <- NULL

  return(changes[order(changes[["name"]]), , drop = FALSE])
}


.scenario_snapshot_review <- function(path) {

  required <- c("shiny", "diffviewer")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) {
    message(
      "Interactive scenario snapshot review is unavailable; install: ",
      paste(required[!available], collapse = ", "), "."
    )
    return(FALSE)
  }

  testthat::snapshot_review(path = path)
  return(TRUE)
}


.scenario_review_text_snapshots <- function(root, scenarios) {

  changes <- .scenario_text_candidates(root, scenarios)
  if (nrow(changes) == 0L) {
    return(invisible(changes))
  }

  message(
    "Scenario text snapshot mismatches:\n- ",
    paste(changes[["name"]], collapse = "\n- ")
  )
  if (!.scenario_is_interactive()) {
    message("Cached candidates were kept for later interactive review.")
    return(invisible(changes))
  }

  review_root <- tempfile("robma-scenario-review-")
  on.exit(unlink(review_root, recursive = TRUE), add = TRUE)
  staged_expected  <- character(nrow(changes))
  staged_candidate <- character(nrow(changes))
  for (i in seq_len(nrow(changes))) {
    parts <- strsplit(changes[["name"]][[i]], "/", fixed = TRUE)[[1L]]
    directory <- file.path(review_root, "_snaps", parts[[1L]])
    staged_expected[[i]] <- file.path(
      directory,
      basename(changes[["expected"]][[i]])
    )
    staged_candidate[[i]] <- .scenario_candidate_path(staged_expected[[i]])
    .scenario_copy_file(
      changes[["expected"]][[i]], staged_expected[[i]], "review snapshot"
    )
    .scenario_copy_file(
      changes[["candidate"]][[i]], staged_candidate[[i]], "review candidate"
    )
  }

  reviewed <- .scenario_snapshot_review(review_root)
  if (!reviewed) {
    return(invisible(changes))
  }

  changes[["status"]] <- "skipped"
  for (i in seq_len(nrow(changes))) {
    if (file.exists(staged_candidate[[i]])) {
      next
    }

    staged <- readLines(
      staged_expected[[i]], warn = FALSE, encoding = "UTF-8"
    )
    old <- readLines(
      changes[["expected"]][[i]], warn = FALSE, encoding = "UTF-8"
    )
    new <- readLines(
      changes[["candidate"]][[i]], warn = FALSE, encoding = "UTF-8"
    )
    if (identical(staged, new)) {
      .scenario_copy_file(
        changes[["candidate"]][[i]],
        changes[["expected"]][[i]],
        "accepted text snapshot"
      )
      changes[["status"]][[i]] <- "accepted"
      unlink(changes[["candidate"]][[i]])
    } else if (identical(staged, old)) {
      changes[["status"]][[i]] <- "rejected"
      unlink(changes[["candidate"]][[i]])
    } else {
      stop(
        "Reviewed scenario text has an unexpected value: ",
        changes[["name"]][[i]], ".",
        call. = FALSE
      )
    }
  }

  message(
    "Scenario text review complete: ",
    paste0(changes[["name"]], " [", changes[["status"]], "]",
           collapse = ", "), "."
  )
  return(invisible(changes))
}


# Capture a result as R would print it at the top level and lock the output to a
# tracked text file.
scenario_text <- function(name, code) {

  set.seed(1)
  config <- .scenario_config()
  name   <- .scenario_validate_name(name, "text")
  path   <- file.path(config[["results_dir"]], paste0(name, ".txt"))
  candidate <- .scenario_candidate_path(path)
  expr   <- substitute(code)
  result <- NULL

  old_options <- options(width = config[["width"]])
  on.exit(options(old_options), add = TRUE)
  started <- .scenario_clock()
  output <- capture.output(
    {
      result <- suppressMessages(withVisible(
        eval(expr, envir = parent.frame())
      ))
      if (result[["visible"]]) {
        print(result[["value"]])
      }
    },
    type = "output"
  )
  elapsed <- max(0, .scenario_clock() - started)
  .scenario_register_timing("text", name, elapsed)
  value <- result[["value"]]
  if (isTRUE(config[["show_output"]]) && length(output) > 0L) {
    writeLines(output)
  }

  if (!file.exists(path) && !config[["create_missing"]] &&
      !config[["update"]]) {
    testthat::fail(paste0(
      "Missing locked scenario text 'results/", config[["name"]], "/",
      basename(path), "'. Create and review it locally."
    ))
    return(invisible(value))
  }

  if (config[["update"]] || !file.exists(path)) {
    .scenario_write_lines(output, path)
    unlink(candidate)
    message(
      if (config[["update"]]) "Updated" else "Added",
      " scenario text: ", path
    )
  } else {
    expected <- readLines(path, warn = FALSE, encoding = "UTF-8")
    if (identical(output, expected)) {
      unlink(candidate)
      testthat::succeed()
    } else {
      .scenario_write_lines(output, candidate)
      testthat::expect_equal(
        output,
        expected,
        info = paste0(
          "Scenario text '", config[["name"]], "/", name,
          "' changed. The candidate is cached at ", candidate,
          " and will be offered for review after test_scenario() completes."
        )
      )
    }
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
    return(NULL)
  }

  return(snapshotter)
}


.scenario_evaluate_plot <- function(expr, envir, restore_par = TRUE) {

  if (restore_par) {
    old_par <- graphics::par(no.readonly = TRUE)
    old_par[["new"]] <- NULL
    graphics::par(new = FALSE)
    on.exit({
      graphics::par(old_par)
      graphics::par(new = FALSE)
    }, add = TRUE)
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


# Draw when interactive output is enabled, and compare with a tracked SVG.
scenario_plot <- function(name, code) {

  set.seed(1)
  config      <- .scenario_config()
  name        <- .scenario_validate_name(name, "plot")
  expr        <- substitute(code)
  eval_env    <- parent.frame()
  show_output <- isTRUE(config[["show_output"]])
  snapshotter <- .scenario_snapshot_context()
  figure_elapsed <- numeric()

  figure <- function() {

    set.seed(1)
    started   <- .scenario_clock()
    completed <- FALSE
    on.exit({
      if (completed) {
        figure_elapsed <<- c(
          figure_elapsed,
          max(0, .scenario_clock() - started)
        )
      }
    }, add = TRUE)
    # vdiffr closes its temporary device, so its par settings cannot leak.
    .scenario_evaluate_plot(expr, eval_env, restore_par = FALSE)
    completed <- TRUE
    return(invisible(NULL))
  }
  on.exit({
    if (length(figure_elapsed) > 0L) {
      .scenario_register_timing("plot", name, sum(figure_elapsed))
    }
  }, add = TRUE)

  if (is.null(snapshotter)) {
    if (!show_output && !.scenario_is_interactive()) {
      stop(
        "scenario_plot() must run through tools/test-scenario.R or testthat::test_file().",
        call. = FALSE
      )
    }

    set.seed(1)
    .scenario_evaluate_plot(expr, eval_env)

    testthat::skip_if_not_installed("vdiffr")
    expected  <- file.path(
      config[["snapshots_dir"]],
      .scenario_plot_filename(name)
    )
    candidate <- .scenario_candidate_path(expected)
    generated <- tempfile("robma-scenario-", fileext = ".svg")
    on.exit(unlink(generated), add = TRUE)
    .write_canonical_svg(figure, generated, name)

    if (!file.exists(expected) && !config[["create_missing"]] &&
        !config[["update"]]) {
      .scenario_copy_file(generated, candidate, "plot candidate")
      testthat::fail(paste0(
        "Missing locked scenario plot '_snaps/", config[["name"]], "/",
        basename(expected), "'. The new snapshot is at ", candidate, "."
      ))
      return(invisible(NULL))
    }

    if (config[["update"]] || !file.exists(expected)) {
      .scenario_copy_file(generated, expected, "plot")
      unlink(candidate)
      message(
        if (config[["update"]]) "Updated" else "Added",
        " scenario plot: ", expected
      )
      testthat::succeed()
    } else if (testthat::compare_file_text(expected, generated)) {
      unlink(candidate)
      testthat::succeed()
    } else {
      .scenario_copy_file(generated, candidate, "plot candidate")
      testthat::fail(paste0(
        "Scenario plot '", config[["name"]], "/", name,
        "' changed. The old snapshot was kept; the new snapshot is at ",
        candidate, "."
      ))
    }

    return(invisible(NULL))
  }

  testthat::skip_if_not_installed("vdiffr")

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
    .scenario_plot_filename(name)
  )
  if (!file.exists(expected) && !config[["create_missing"]] &&
      !config[["update"]]) {
    testthat::fail(paste0(
      "Missing locked scenario plot '_snaps/", config[["name"]], "/",
      basename(expected), "'. Create and review it locally."
    ))
    return(invisible(NULL))
  }

  if (show_output) {
    set.seed(1)
    .scenario_evaluate_plot(expr, eval_env)
  }

  writer <- function(plot, file, title = "") {

    .write_canonical_svg(plot, file, title)
    if (config[["update"]]) {
      .scenario_copy_file(file, expected, "plot")
      unlink(.scenario_candidate_path(expected))
    }
  }

  vdiffr::expect_doppelganger(
    title  = name,
    fig    = figure,
    writer = writer
  )

  return(invisible(NULL))
}


# Compare estimates against a reference on a common difference scale.
scenario_agreement_plot <- function(reference, estimate, main = "",
                                    reference_label = "metafor",
                                    estimate_label = "RoBMA") {

  if (!is.numeric(reference) || !is.numeric(estimate)) {
    stop("'reference' and 'estimate' must be numeric vectors.", call. = FALSE)
  }
  if (length(reference) != length(estimate)) {
    stop("'reference' and 'estimate' must have the same length.",
         call. = FALSE)
  }
  if (!is.character(reference_label) || length(reference_label) != 1L ||
      is.na(reference_label) || !nzchar(reference_label)) {
    stop("'reference_label' must be one non-empty string.", call. = FALSE)
  }
  if (!is.character(estimate_label) || length(estimate_label) != 1L ||
      is.na(estimate_label) || !nzchar(estimate_label)) {
    stop("'estimate_label' must be one non-empty string.", call. = FALSE)
  }

  difference <- estimate - reference
  keep       <- is.finite(reference) & is.finite(estimate)
  if (!any(keep)) {
    stop("Agreement plot is unavailable: no finite reference-estimate pairs.",
         call. = FALSE)
  }
  if (any(!keep)) {
    warning(
      sum(!keep), " non-finite reference-estimate pairs were removed.",
      immediate. = TRUE
    )
  }
  band <- if (sum(keep) == 1L) {
    0
  } else {
    0.1 * stats::sd(reference[keep])
  }
  y_limit <- max(2 * band, 1.05 * abs(difference[keep]))

  graphics::plot(
    reference[keep], difference[keep],
    main = main, xlab = reference_label,
    ylab = paste0(estimate_label, " - ", reference_label),
    ylim = c(-y_limit, y_limit), type = "n"
  )
  graphics::rect(
    graphics::par("usr")[[1L]], -band,
    graphics::par("usr")[[2L]],  band,
    col = "grey90", border = NA
  )
  graphics::abline(h = 0, lty = 2, col = "grey50")
  graphics::points(reference[keep], difference[keep], pch = 19, cex = 0.7)

  return(invisible(NULL))
}


.scenario_artifacts_in_file <- function(path) {

  artifacts <- list(
    text    = character(),
    plot    = character(),
    dynamic = c(text = FALSE, plot = FALSE)
  )
  expressions <- parse(path, keep.source = FALSE)

  walk <- function(node) {

    if (!is.call(node)) {
      return(invisible(NULL))
    }

    call_name <- if (is.symbol(node[[1L]])) {
      as.character(node[[1L]])
    } else {
      ""
    }
    type <- switch(
      call_name,
      scenario_text = "text",
      scenario_plot = "plot",
      NULL
    )
    if (!is.null(type)) {
      artifact_name <- node[[2L]]
      if (is.character(artifact_name) && length(artifact_name) == 1L &&
          !is.na(artifact_name)) {
        artifacts[[type]] <<- c(artifacts[[type]], artifact_name)
      } else {
        artifacts[["dynamic"]][[type]] <<- TRUE
      }
    }

    children <- as.list(node)[-1L]
    for (i in seq_along(children)) {
      if (identical(children[[i]], quote(expr = ))) {
        next
      }
      walk(children[[i]])
    }
    return(invisible(NULL))
  }

  for (expression in expressions) {
    walk(expression)
  }
  artifacts[["text"]] <- unique(artifacts[["text"]])
  artifacts[["plot"]] <- unique(artifacts[["plot"]])

  return(artifacts)
}


.scenario_orphan_artifacts <- function(path) {

  filename <- basename(path)
  if (!grepl("^test-.*\\.[Rr]$", filename)) {
    stop("Scenario files must be named 'test-<scenario>.R'.", call. = FALSE)
  }

  scenario_name <- sub("^test-(.*)\\.[Rr]$", "\\1", filename)
  scenario_root <- dirname(path)
  artifacts     <- .scenario_artifacts_in_file(path)
  orphans       <- character()

  specifications <- list(
    text = list(
      directory = file.path(scenario_root, "results", scenario_name),
      extension = ".txt",
      candidate = "\\.new\\.txt$"
    ),
    plot = list(
      directory = file.path(scenario_root, "_snaps", scenario_name),
      extension = ".svg",
      candidate = "\\.new\\.svg$"
    )
  )

  for (type in names(specifications)) {
    if (artifacts[["dynamic"]][[type]]) {
      next
    }
    specification <- specifications[[type]]
    existing <- list.files(
      specification[["directory"]],
      pattern    = paste0("\\", specification[["extension"]], "$"),
      full.names = FALSE
    )
    existing <- existing[!grepl(specification[["candidate"]], existing)]
    expected <- if (type == "plot") {
      vapply(artifacts[[type]], .scenario_plot_filename, character(1))
    } else {
      paste0(artifacts[[type]], specification[["extension"]])
    }
    unused   <- setdiff(existing, expected)
    if (length(unused) > 0L) {
      directory <- if (type == "text") "results" else "_snaps"
      orphans <- c(
        orphans,
        file.path(directory, scenario_name, unused)
      )
    }
  }

  return(chartr("\\", "/", sort(orphans)))
}


.scenario_report_orphans <- function(path) {

  orphans <- .scenario_orphan_artifacts(path)
  if (length(orphans) > 0L) {
    message(
      "Orphaned scenario artifacts are not referenced by ", basename(path),
      ":\n- ", paste(orphans, collapse = "\n- ")
    )
  }

  return(invisible(orphans))
}


.scenario_list_files <- function(root = .scenario_helpers_dir) {

  files <- list.files(
    root,
    pattern    = "^test-[a-z0-9]([a-z0-9._-]*[a-z0-9_-])?\\.[Rr]$",
    full.names = TRUE
  )
  scenario_names <- sub(
    "^test-([a-z0-9]([a-z0-9._-]*[a-z0-9_-])?)\\.[Rr]$",
    "\\1",
    basename(files)
  )
  names(files) <- scenario_names

  return(files)
}


.scenario_filter_files <- function(files, filter = NULL) {

  if (is.null(filter)) {
    return(files)
  }
  if (!is.character(filter) || length(filter) != 1L || is.na(filter)) {
    stop("'filter' must be NULL or one regular expression.", call. = FALSE)
  }

  selected <- tryCatch(
    grepl(filter, names(files)),
    error = function(error) {
      stop("Invalid scenario 'filter': ", conditionMessage(error),
           call. = FALSE)
    }
  )
  files <- files[selected]
  if (length(files) == 0L) {
    stop("No scenarios matched 'filter = ", encodeString(filter, quote = "'"),
         ".", call. = FALSE)
  }

  return(files)
}


.scenario_test_file <- function(path, reporter, stop_on_failure) {

  return(testthat::test_file(
    path,
    reporter        = reporter,
    package         = "RoBMA",
    load_package    = "none",
    stop_on_failure = stop_on_failure
  ))
}


# Run maintainer scenarios from an interactive development session. `filter`
# is a regular expression matched against scenario names.
test_scenario <- function(filter = NULL, reporter = "progress", refit = FALSE,
                          update = FALSE, update_timings = FALSE,
                          regenerate = FALSE,
                          load_package = TRUE, stop_on_failure = FALSE,
                          root = .scenario_helpers_dir) {

  refit           <- .scenario_validate_flag(refit, "refit")
  update          <- .scenario_validate_flag(update, "update")
  update_timings  <- .scenario_validate_flag(
    update_timings,
    "update_timings"
  )
  regenerate      <- .scenario_validate_flag(regenerate, "regenerate")
  load_package    <- .scenario_validate_flag(load_package, "load_package")
  stop_on_failure <- .scenario_validate_flag(
    stop_on_failure,
    "stop_on_failure"
  )
  refit  <- refit || regenerate
  update <- update || regenerate
  update_timings <- update_timings || update

  root         <- normalizePath(root, winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(
    file.path(root, "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  files <- .scenario_list_files(root)
  if (length(files) == 0L) {
    stop("No scenarios are defined.", call. = FALSE)
  }
  files <- .scenario_filter_files(files, filter)

  if (load_package) {
    devtools::load_all(project_root, quiet = TRUE)
  }

  environment_names <- c(
    "AGENT",
    "NOT_CRAN",
    "ROBMA_SCENARIO_REGENERATE",
    "ROBMA_SCENARIO_REFIT",
    "ROBMA_SCENARIO_UPDATE",
    "ROBMA_SCENARIO_UPDATE_TIMINGS",
    "ROBMA_SCENARIO_RUNNER"
  )
  old_environment <- Sys.getenv(environment_names, unset = NA_character_)
  old_directory   <- getwd()
  on.exit({
    setwd(old_directory)
    for (i in seq_along(environment_names)) {
      name <- environment_names[[i]]
      if (is.na(old_environment[[i]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old_environment[[i]]), name))
      }
    }
  }, add = TRUE)
  on.exit(
    .scenario_review_text_snapshots(root, names(files)),
    add = TRUE
  )

  Sys.setenv(
    AGENT                             = "1",
    NOT_CRAN                          = "true",
    ROBMA_SCENARIO_REGENERATE         = if (regenerate) "TRUE" else "FALSE",
    ROBMA_SCENARIO_REFIT              = if (refit) "TRUE" else "FALSE",
    ROBMA_SCENARIO_UPDATE             = if (update) "TRUE" else "FALSE",
    ROBMA_SCENARIO_UPDATE_TIMINGS     = if (update_timings) "TRUE" else "FALSE",
    ROBMA_SCENARIO_RUNNER             = "TRUE"
  )
  setwd(project_root)

  results <- vector("list", length(files))
  names(results) <- names(files)
  for (i in seq_along(files)) {
    path          <- files[[i]]
    scenario_name <- names(files)[[i]]
    completed     <- FALSE
    message("Running scenario '", scenario_name, "'.")
    results[[i]] <- tryCatch(
      {
        result    <- .scenario_test_file(path, reporter, stop_on_failure)
        completed <- TRUE
        result
      },
      finally = {
        .scenario_finalize_timing(scenario_name, allow_update = completed)
        .scenario_report_orphans(path)
      }
    )
  }
  return(invisible(results))
}
