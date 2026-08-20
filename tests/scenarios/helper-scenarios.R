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

if (!exists(".scenario_state", inherits = TRUE)) {
  .scenario_state <- new.env(parent = emptyenv())
}


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
    update_timings <- FALSE
  }
  update_timings <- .scenario_validate_flag(update_timings, "update_timings")
  update_timings <- update_timings || env_update_timings

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
    cached_phases <- if (is.list(fit_timing)) {
      fit_timing[["phases"]]
    } else {
      NULL
    }
    valid_phases <- is.numeric(cached_phases) &&
      !is.null(names(cached_phases)) &&
      !anyNA(names(cached_phases)) &&
      !any(!nzchar(names(cached_phases))) &&
      !anyDuplicated(names(cached_phases)) &&
      all(names(cached_phases) %in% c("model", "loo", "marglik")) &&
      all(is.finite(cached_phases)) &&
      all(cached_phases >= 0) &&
      (length(cached_phases) == 0L ||
       identical(names(cached_phases)[[1L]], "model"))
    valid_timing <- is.list(fit_timing) &&
      identical(fit_timing[["version"]], 2L) &&
      identical(fit_timing[["fit_call"]], fit_call) &&
      is.numeric(fit_timing[["elapsed"]]) &&
      length(fit_timing[["elapsed"]]) == 1L &&
      is.finite(fit_timing[["elapsed"]]) &&
      fit_timing[["elapsed"]] >= 0 &&
      valid_phases
    if (valid_timing) {
      .scenario_register_fit_timings(
        name       = name,
        elapsed    = fit_timing[["elapsed"]],
        phases     = fit_timing[["phases"]],
        provenance = fit_timing[["provenance"]]
      )
    } else {
      .scenario_mark_timing_unavailable(
        "fit", name,
        paste(
          "the cached fit predates valid timing metadata with post-fit phases;",
          "rerun with",
          "'refit = TRUE'"
        )
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

  started    <- .scenario_clock()
  evaluation <- .scenario_evaluate_fit(expr, envir = parent.frame())
  elapsed    <- max(0, .scenario_clock() - started)
  fit        <- evaluation[["fit"]]
  phases     <- evaluation[["phases"]]
  if (length(phases) > 0L) {
    phases <- c(model = max(0, elapsed - sum(phases)), phases)
  }
  provenance <- .scenario_timing_provenance()
  dir.create(config[["cache_dir"]], recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, paths[["fit"]])
  saveRDS(
    list(
      version    = 2L,
      fit_call   = fit_call,
      elapsed    = elapsed,
      phases     = phases,
      provenance = provenance
    ),
    paths[["timing"]]
  )
  .scenario_write_lines(fit_call, paths[["call"]])
  .scenario_register_fit_timings(name, elapsed, phases, provenance)

  return(fit)
}


.scenario_write_lines <- function(lines, path) {

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(lines, connection, useBytes = TRUE)

  return(invisible(path))
}


.scenario_evaluate_fit <- function(expr, envir) {

  phase_functions <- c(loo = "add_loo", marglik = "add_marglik")
  collector <- new.env(parent = emptyenv())
  collector[["elapsed"]] <- stats::setNames(numeric(), character())
  bindings <- vector("list", length(phase_functions))
  names(bindings) <- names(phase_functions)

  for (phase in names(phase_functions)) {
    function_name <- phase_functions[[phase]]
    if (!exists(function_name, envir = envir, mode = "function",
                inherits = TRUE)) {
      next
    }
    bindings[[phase]] <- list(
      name      = function_name,
      local     = exists(function_name, envir = envir, inherits = FALSE),
      value     = if (exists(function_name, envir = envir, inherits = FALSE)) {
        get(function_name, envir = envir, inherits = FALSE)
      } else {
        NULL
      },
      original  = get(function_name, envir = envir, mode = "function",
                      inherits = TRUE)
    )
  }
  on.exit({
    for (binding in rev(bindings)) {
      if (is.null(binding)) {
        next
      }
      if (binding[["local"]]) {
        assign(binding[["name"]], binding[["value"]], envir = envir)
      } else if (exists(binding[["name"]], envir = envir, inherits = FALSE)) {
        rm(list = binding[["name"]], envir = envir)
      }
    }
  }, add = TRUE)

  for (phase in names(bindings)) {
    binding <- bindings[[phase]]
    if (is.null(binding)) {
      next
    }
    wrapper <- local({
      this_phase <- phase
      original   <- binding[["original"]]
      function(...) {

        started <- .scenario_clock()
        on.exit({
          phase_elapsed <- max(0, .scenario_clock() - started)
          elapsed <- collector[["elapsed"]]
          if (!this_phase %in% names(elapsed)) {
            elapsed[[this_phase]] <- 0
          }
          elapsed[[this_phase]] <- elapsed[[this_phase]] + phase_elapsed
          collector[["elapsed"]] <- elapsed
        }, add = TRUE)
        original(...)
      }
    })
    assign(binding[["name"]], wrapper, envir = envir)
  }

  fit <- eval(expr, envir = envir)
  return(list(fit = fit, phases = collector[["elapsed"]]))
}


.scenario_register_fit_timings <- function(name, elapsed, phases,
                                            provenance) {

  .scenario_register_timing("fit", name, elapsed, provenance)
  phase_types <- c(
    model   = "fit_model",
    loo     = "fit_loo",
    marglik = "fit_marglik"
  )
  for (phase in names(phases)) {
    .scenario_register_timing(
      phase_types[[phase]], name, phases[[phase]], provenance
    )
  }

  return(invisible(elapsed))
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


.scenario_timing_types <- function() {

  return(c(
    "fit", "fit_model", "fit_loo", "fit_marglik", "time", "text", "plot"
  ))
}


.scenario_order_timings <- function(timings) {

  type_order <- match(timings[["type"]], .scenario_timing_types())
  timings <- timings[order(type_order, timings[["name"]]), , drop = FALSE]
  rownames(timings) <- NULL
  return(timings)
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
    timings[["type"]] %in% .scenario_timing_types() &
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
  assessable   <- current[["elapsed"]] >= 0.5
  comparable   <- !is.na(matched) & assessable
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
    missing <- setdiff(baseline_key, observed_key)
    if (length(missing) > 0L) {
      add_issue(
        "structure",
        paste0("timing entries are missing: ", paste(missing, collapse = ", "))
      )
    }
  }

  same_entries <- nrow(baseline) > 0L && nrow(unavailable) == 0L &&
    setequal(current_key, baseline_key)
  if (final && same_entries) {
    split_fit_names <- current[["name"]][current[["type"]] == "fit_model"]
    average_rows <- assessable & !(current[["type"]] == "fit" &
      current[["name"]] %in% split_fit_names)
    average_current <- current[average_rows, , drop = FALSE]
    if (nrow(average_current) > 0L) {
      average_key <- .scenario_timing_key(average_current)
      baseline_elapsed <- baseline[["elapsed"]][match(
        average_key, baseline_key
      )]
      percentage_change <- ifelse(
        baseline_elapsed == 0,
        ifelse(average_current[["elapsed"]] == 0, 0, Inf),
        100 * (average_current[["elapsed"]] / baseline_elapsed - 1)
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


.scenario_merge_timings <- function(baseline, current,
                                    accept_slower = FALSE) {

  if (nrow(current) == 0L) {
    return(baseline)
  }
  baseline_key <- .scenario_timing_key(baseline)
  current_key  <- .scenario_timing_key(current)
  matched      <- match(current_key, baseline_key)
  replace <- is.na(matched) | accept_slower
  comparable <- !is.na(matched)
  replace[comparable] <- replace[comparable] |
    current[["elapsed"]][comparable] <
      baseline[["elapsed"]][matched[comparable]]
  replacement_key <- current_key[replace]
  keep <- !baseline_key %in% replacement_key

  return(.scenario_order_timings(rbind(
    baseline[keep, , drop = FALSE], current[replace, , drop = FALSE]
  )))
}


.scenario_timing_candidate_needed <- function(baseline, current, unavailable,
                                              final = FALSE) {

  baseline_key <- .scenario_timing_key(baseline)
  current_key  <- .scenario_timing_key(current)
  matched      <- match(current_key, baseline_key)
  slower <- !is.na(matched) &
    current[["elapsed"]] > baseline[["elapsed"]][matched]
  missing <- final && length(setdiff(
    baseline_key,
    c(current_key, .scenario_timing_key(unavailable))
  )) > 0L

  return(any(slower) || nrow(unavailable) > 0L || missing)
}


.scenario_register_timing <- function(type, name, elapsed,
                                      provenance = .scenario_timing_provenance()) {

  type <- match.arg(type, .scenario_timing_types())
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
    updated <- .scenario_merge_timings(
      baseline,
      current,
      accept_slower = config[["update_timings"]]
    )
    if (!identical(updated, baseline)) {
      .scenario_write_timings(updated, paths[["expected"]])
      assign("timing_baseline", updated, envir = .scenario_state)
    }
    if (.scenario_timing_candidate_needed(
      updated, current, unavailable, final = FALSE
    )) {
      .scenario_write_timings(current, paths[["candidate"]])
    } else {
      unlink(paths[["candidate"]])
    }
  }

  return(invisible(elapsed))
}


.scenario_mark_timing_unavailable <- function(type, name, reason) {

  type <- match.arg(type, .scenario_timing_types())
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
  issues <- .scenario_timing_issues(final = TRUE)
  .scenario_warn_timing_issues(issues)

  baseline <- get(
    "timing_baseline", envir = .scenario_state, inherits = FALSE
  )
  if (allow_update && nrow(current) > 0L) {
    existed <- file.exists(paths[["expected"]])
    replace_all <- config[["update_timings"]] && nrow(unavailable) == 0L
    updated <- if (replace_all) {
      current
    } else {
      .scenario_merge_timings(
        baseline,
        current,
        accept_slower = config[["update_timings"]]
      )
    }
    if (!identical(updated, baseline)) {
      .scenario_write_timings(updated, paths[["expected"]])
      assign("timing_baseline", updated, envir = .scenario_state)
      message(
        if (existed) "Updated" else "Added",
        " scenario timings: ", paths[["expected"]]
      )
    }
    baseline <- updated
  } else if (config[["update_timings"]] && !allow_update) {
    warning(
      "Scenario timing baseline for '", config[["name"]],
      "' was not updated because the scenario did not finish.",
      call. = FALSE
    )
  }

  keep_candidate <- !allow_update || .scenario_timing_candidate_needed(
    baseline, current, unavailable, final = TRUE
  )
  if (nrow(current) > 0L && keep_candidate) {
    .scenario_write_timings(current, paths[["candidate"]])
  } else {
    unlink(paths[["candidate"]])
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
      type     = rep("table", length(candidates)),
      scenario = rep(scenario, length(candidates)),
      name = if (length(candidates) == 0L) character() else paste0(
        scenario, "/", sub("[.]new[.]txt$", "", basename(candidates))
      ),
      expected  = expected,
      candidate = candidates,
      stringsAsFactors = FALSE
    )
  })
  changes <- do.call(rbind, changes)
  if (is.null(changes)) {
    changes <- data.frame(
      type      = character(),
      scenario  = character(),
      name      = character(),
      expected  = character(),
      candidate = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(changes) <- NULL

  return(changes[order(changes[["scenario"]], changes[["name"]]), , drop = FALSE])
}


.scenario_plot_candidates <- function(root, scenarios) {

  changes <- lapply(scenarios, function(scenario) {
    candidates <- list.files(
      file.path(root, "_snaps", scenario),
      pattern    = "[.]new[.]svg$",
      full.names = TRUE
    )
    expected   <- sub("[.]new[.]svg$", ".svg", candidates)
    keep       <- file.exists(expected)
    candidates <- candidates[keep]
    expected   <- expected[keep]

    data.frame(
      type      = rep("figure", length(candidates)),
      scenario  = rep(scenario, length(candidates)),
      name      = if (length(candidates) == 0L) character() else paste0(
        scenario, "/", sub("[.]new[.]svg$", "", basename(candidates))
      ),
      expected  = expected,
      candidate = candidates,
      stringsAsFactors = FALSE
    )
  })
  changes <- do.call(rbind, changes)
  if (is.null(changes)) {
    changes <- data.frame(
      type      = character(),
      scenario  = character(),
      name      = character(),
      expected  = character(),
      candidate = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(changes) <- NULL

  return(changes[order(changes[["scenario"]], changes[["name"]]), , drop = FALSE])
}


.scenario_snapshot_candidates <- function(root, scenarios,
                                          types = c("table", "figure")) {

  types <- match.arg(types, c("table", "figure"), several.ok = TRUE)
  changes <- list()
  if ("table" %in% types) {
    changes[["table"]] <- .scenario_text_candidates(root, scenarios)
  }
  if ("figure" %in% types) {
    changes[["figure"]] <- .scenario_plot_candidates(root, scenarios)
  }
  changes <- do.call(rbind, changes)
  rownames(changes) <- NULL

  return(changes[order(
    changes[["scenario"]], match(changes[["type"]], c("table", "figure")),
    changes[["name"]]
  ), , drop = FALSE])
}


.scenario_snapshot_review <- function(path, files = NULL) {

  required <- c("shiny", "diffviewer")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) {
    message(
      "Interactive snapshot review is unavailable; install: ",
      paste(required[!available], collapse = ", "), "."
    )
    return(FALSE)
  }

  testthat::snapshot_review(files = files, path = path)
  return(TRUE)
}


.scenario_files_identical <- function(x, y) {

  return(identical(
    unname(tools::md5sum(x)),
    unname(tools::md5sum(y))
  ))
}


.scenario_review_snapshots <- function(root, scenarios,
                                       types = c("table", "figure")) {

  changes <- .scenario_snapshot_candidates(root, scenarios, types)
  if (nrow(changes) == 0L) {
    message("No cached scenario snapshots to review.")
    return(invisible(changes))
  }

  message(
    "Scenario snapshot mismatches:\n- ",
    paste0(
      "[", changes[["type"]], "] ", changes[["name"]],
      collapse = "\n- "
    )
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
    directory <- file.path(
      review_root, "_snaps", changes[["scenario"]][[i]]
    )
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

    if (.scenario_files_identical(
      staged_expected[[i]], changes[["candidate"]][[i]]
    )) {
      .scenario_copy_file(
        changes[["candidate"]][[i]],
        changes[["expected"]][[i]],
        paste("accepted", changes[["type"]][[i]], "snapshot")
      )
      changes[["status"]][[i]] <- "accepted"
      unlink(changes[["candidate"]][[i]])
    } else if (.scenario_files_identical(
      staged_expected[[i]], changes[["expected"]][[i]]
    )) {
      changes[["status"]][[i]] <- "rejected"
      unlink(changes[["candidate"]][[i]])
    } else {
      stop(
        "Reviewed scenario snapshot has an unexpected value: ",
        changes[["name"]][[i]], ".",
        call. = FALSE
      )
    }
  }

  message(
    "Scenario snapshot review complete: ",
    paste0(
      changes[["name"]], " [", changes[["type"]], ":",
      changes[["status"]], "]", collapse = ", "
    ), "."
  )
  return(invisible(changes))
}


.scenario_review_text_snapshots <- function(root, scenarios) {

  return(.scenario_review_snapshots(root, scenarios, types = "table"))
}


review_scenario_snapshots <- function(filter = NULL, root = NULL) {

  has_last_run <- exists(
    "last_scenario_run", envir = .scenario_state, inherits = FALSE
  )
  if (is.null(root) && has_last_run) {
    last_run  <- get(
      "last_scenario_run", envir = .scenario_state, inherits = FALSE
    )
    root      <- last_run[["root"]]
    scenarios <- last_run[["scenarios"]]
  } else {
    if (is.null(root)) {
      root <- .scenario_helpers_dir
    }
    root      <- normalizePath(root, winslash = "/", mustWork = TRUE)
    scenarios <- names(.scenario_list_files(root))
  }
  files <- stats::setNames(rep("", length(scenarios)), scenarios)
  files <- .scenario_filter_files(files, filter)

  return(.scenario_review_snapshots(root, names(files)))
}


review_test_snapshots <- function(files = NULL,
                                  root = file.path(
                                    dirname(.scenario_helpers_dir), "testthat"
                                  )) {

  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!.scenario_is_interactive()) {
    message("Cached test snapshots require an interactive session for review.")
    return(invisible(FALSE))
  }

  return(invisible(.scenario_snapshot_review(root, files = files)))
}


# Evaluate an ordinary scenario computation and record its wall time.
scenario_time <- function(name, code) {

  .scenario_config()
  name    <- .scenario_validate_name(name, "time")
  expr    <- substitute(code)
  started <- .scenario_clock()
  result  <- withVisible(eval(expr, envir = parent.frame()))
  elapsed <- max(0, .scenario_clock() - started)
  .scenario_register_timing("time", name, elapsed)

  if (result[["visible"]]) {
    return(result[["value"]])
  }
  return(invisible(result[["value"]]))
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

    testthat::skip_if_not_installed("vdiffr")
    expected  <- file.path(
      config[["snapshots_dir"]],
      .scenario_plot_filename(name)
    )
    candidate <- .scenario_candidate_path(expected)
    generated <- tempfile("robma-scenario-", fileext = ".svg")
    on.exit(unlink(generated), add = TRUE)
    .write_canonical_svg(figure, generated, name)

    if (show_output) {
      set.seed(1)
      .scenario_evaluate_plot(expr, eval_env)
    }

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

  if (show_output) {
    set.seed(1)
    .scenario_evaluate_plot(expr, eval_env)
  }

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


.scenario_ex_string <- function(value, name, allow_null = FALSE) {

  if (allow_null && is.null(value)) {
    return(invisible(NULL))
  }
  if (!is.character(value) || length(value) != 1L ||
      is.na(value) || !nzchar(value)) {
    requirement <- if (allow_null) "NULL or one non-empty string" else "one non-empty string"
    stop("'", name, "' must be ", requirement, ".", call. = FALSE)
  }
  return(invisible(value))
}


.scenario_ex_unavailable <- function(source, parameter) {

  error <- simpleError(paste0(
    source, " estimate '", parameter, "' is unavailable."
  ))
  class(error) <- c("scenario_estimate_unavailable", class(error))
  stop(error)
}


# Extract one coefficient, variance component, or correlation from metafor.
.scenario_ex_m_one <- function(fit, parameter, statistic, coefficient_table) {

  .scenario_ex_string(parameter, "parameter")
  .scenario_ex_string(statistic, "statistic")
  unavailable <- function() {
    .scenario_ex_unavailable("Metafor", parameter)
  }
  statistic_unavailable <- function() {
    stop(
      "Statistic '", statistic, "' is unavailable for metafor estimate '",
      parameter, "'.",
      call. = FALSE
    )
  }

  if (!is.null(coefficient_table) && nrow(coefficient_table) > 0L) {
    normalize_name <- function(x) tolower(gsub("[^[:alnum:]]", "", x))
    coefficient_name <- if (parameter %in% c("mu", "intercept")) "intrcpt" else parameter
    rows <- which(normalize_name(rownames(coefficient_table)) == normalize_name(coefficient_name))
    if (length(rows) > 1L) {
      stop("Metafor estimate '", parameter, "' is ambiguous.", call. = FALSE)
    }
    if (length(rows) == 1L) {
      statistic_name <- switch(
        tolower(statistic),
        mean = "estimate", sd = "se", ci_0.025 = "ci.lb", ci_0.975 = "ci.ub",
        statistic
      )
      column <- match(tolower(statistic_name), tolower(names(coefficient_table)))
      if (is.na(column)) {
        statistic_unavailable()
      }
      return(unname(coefficient_table[[column]][[rows]]))
    }
  }

  parts <- regmatches(
    parameter,
    regexec("^(sigma|tau|gamma|rho|phi)(?:\\[([^]]+)\\])?(\\^2)?$", parameter, perl = TRUE)
  )[[1L]]
  if (length(parts) == 0L) {
    unavailable()
  }
  if (!tolower(statistic) %in% c("estimate", "mean")) {
    statistic_unavailable()
  }

  family_name <- parts[[2L]]
  index_name  <- parts[[3L]]
  squared     <- nzchar(parts[[4L]])
  is_variance <- family_name %in% c("sigma", "tau", "gamma")
  field_name  <- if (is_variance) paste0(family_name, "2") else family_name
  values      <- fit[[field_name]]
  if (is.null(values) || length(values) == 0L || (!is_variance && squared)) {
    unavailable()
  }

  if (identical(index_name, "total")) {
    if (!identical(family_name, "sigma")) {
      unavailable()
    }
    value <- sum(values)
  } else {
    labels <- switch(
      family_name,
      sigma = fit[["s.names"]],
      tau   = fit[["g.levels.f"]],
      gamma = fit[["h.levels.f"]],
      NULL
    )
    if (is.list(labels)) {
      labels <- labels[[1L]]
    }
    if (!nzchar(index_name)) {
      if (length(values) != 1L) {
        stop(
          "Metafor estimate '", parameter,
          "' is ambiguous; specify its index or label.",
          call. = FALSE
        )
      }
      index <- 1L
    } else if (grepl("^[0-9]+$", index_name)) {
      index <- as.integer(index_name)
    } else {
      index <- match(index_name, as.character(labels))
    }
    if (is.na(index) || index < 1L || index > length(values)) {
      unavailable()
    }
    value <- values[[index]]
  }

  if (is_variance && !squared) {
    value <- sqrt(value)
  }
  return(unname(value))
}


# Extract one named estimate from long-form RoBMA summary output.
.scenario_ex_r_one <- function(parameter, component, statistic, tables) {

  .scenario_ex_string(parameter, "parameter")
  .scenario_ex_string(component, "component", allow_null = TRUE)
  .scenario_ex_string(statistic, "statistic")
  parameter_aliases <- switch(
    parameter,
    mu        = c("mu", "intercept"),
    intercept = c("intercept", "mu"),
    parameter
  )

  found_parameter <- FALSE
  for (table in tables) {
    if (is.null(table) || !"parameter" %in% names(table) ||
        (!is.null(component) && !"component" %in% names(table))) {
      next
    }
    matches <- table[["parameter"]] %in% parameter_aliases
    if (!is.null(component)) {
      matches <- matches & table[["component"]] == component
    }
    rows            <- which(matches)
    found_parameter <- found_parameter || length(rows) > 0L
    if (length(rows) > 1L) {
      stop(
        "RoBMA estimate '", parameter,
        "' is ambiguous; specify 'component'.",
        call. = FALSE
      )
    }
    if (length(rows) == 1L && statistic %in% names(table)) {
      return(unname(table[[statistic]][[rows]]))
    }
  }

  if (found_parameter) {
    stop(
      "Statistic '", statistic, "' is unavailable for RoBMA estimate '",
      parameter, "'.",
      call. = FALSE
    )
  }
  .scenario_ex_unavailable("RoBMA", parameter)
}


.scenario_ex_values <- function(fit, parameter, component, statistic, extractor) {

  if (!is.character(parameter) || length(parameter) == 0L ||
      anyNA(parameter) || any(!nzchar(parameter))) {
    stop("'parameter' must contain one or more non-empty strings.", call. = FALSE)
  }
  n_parameters <- length(parameter)
  if (!is.character(statistic) || !length(statistic) %in% c(1L, n_parameters) ||
      anyNA(statistic) || any(!nzchar(statistic))) {
    stop(
      "'statistic' must contain one string or one string per parameter.",
      call. = FALSE
    )
  }
  statistic <- rep(statistic, length.out = n_parameters)

  if (is.null(component)) {
    component <- rep(list(NULL), n_parameters)
  } else {
    if (!is.character(component) || !length(component) %in% c(1L, n_parameters)) {
      stop(
        "'component' must be NULL, one string, or one string per parameter.",
        call. = FALSE
      )
    }
    component <- rep(component, length.out = n_parameters)
    component <- lapply(component, function(value) if (is.na(value)) NULL else value)
  }

  values <- vapply(seq_len(n_parameters), function(i) {
    tryCatch(
      extractor(fit, parameter[[i]], component[[i]], statistic[[i]]),
      scenario_estimate_unavailable = function(error) NA_real_
    )
  }, numeric(1L))
  if (n_parameters == 1L) {
    return(unname(values))
  }
  parameter_names <- names(parameter)
  if (is.null(parameter_names)) {
    parameter_names <- parameter
  } else {
    unnamed <- !nzchar(parameter_names)
    parameter_names[unnamed] <- parameter[unnamed]
  }
  names(values) <- parameter_names
  return(values)
}


# Extract one or more metafor estimates. Unavailable selectors return NA.
ex_m <- function(fit, parameter, statistic = "estimate") {

  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required to extract metafor estimates.", call. = FALSE)
  }
  coefficient_table <- as.data.frame(stats::coef(summary(fit)))
  return(.scenario_ex_values(
    fit, parameter, NULL, statistic,
    function(fit, parameter, component, statistic) {
      .scenario_ex_m_one(fit, parameter, statistic, coefficient_table)
    }
  ))
}


# Extract one or more RoBMA estimates from long-form summary output. Use
# 'component' when a parameter occurs in more than one summary component.
ex_r <- function(fit, parameter, component = NULL, statistic = "Mean") {

  tables <- list(
    as.data.frame(summary(fit, include_mcmc_diagnostics = FALSE)),
    tryCatch(
      as.data.frame(summary_heterogeneity(fit)),
      error = function(e) NULL
    )
  )
  return(.scenario_ex_values(
    fit, parameter, component, statistic,
    function(fit, parameter, component, statistic) {
      .scenario_ex_r_one(parameter, component, statistic, tables)
    }
  ))
}


# Extract pooled-effect statistics aligned with metafor::predict().
ex_p <- function(fit) {

  output     <- data.frame(pooled_effect(fit))
  statistics <- c("Mean", "CI_0.025", "CI_0.975", "PI_0.025", "PI_0.975")
  values     <- as.data.frame(t(output[, statistics, drop = FALSE]), check.names = FALSE)
  names(values)    <- output[["parameter"]]
  rownames(values) <- c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")
  return(values)
}


# Dispatch by model implementation. Parameter vectors return named vectors;
# model lists return model-named vectors or parameter-column data frames.
ex <- function(fit, parameter, component = NULL, statistic = NULL, ...) {

  UseMethod("ex")
}


ex.rma <- function(fit, parameter, component = NULL, statistic = NULL, ...) {

  if (is.null(statistic)) {
    statistic <- "estimate"
  }
  return(ex_m(fit, parameter, statistic))
}


ex.brma <- function(fit, parameter, component = NULL, statistic = NULL, ...) {

  if (is.null(statistic)) {
    statistic <- "Mean"
  }
  return(ex_r(fit, parameter, component, statistic))
}


ex.list <- function(fit, parameter, component = NULL, statistic = NULL, ...) {

  if (length(fit) == 0L) {
    stop("'fit' must contain one or more models.", call. = FALSE)
  }
  model_names <- names(fit)
  if (is.null(model_names)) {
    model_names <- paste("model", seq_along(fit))
  } else {
    unnamed <- !nzchar(model_names)
    model_names[unnamed] <- paste("model", which(unnamed))
  }
  values <- lapply(fit, function(model) {
    ex(
      model,
      parameter = parameter,
      component = component,
      statistic = statistic,
      ...
    )
  })
  if (length(parameter) == 1L) {
    values <- unlist(values, use.names = FALSE)
    names(values) <- model_names
    return(values)
  }
  values <- as.data.frame(do.call(rbind, values), check.names = FALSE)
  rownames(values) <- model_names
  return(values)
}


ex.default <- function(fit, parameter, component = NULL, statistic = NULL, ...) {

  stop(
    "'fit' must be a metafor model, a RoBMA model, or a list of models.",
    call. = FALSE
  )
}


# Compare marginal diagnostics from a reference fit and a brma.mv fit.
plot_marginal_diagnostics <- function(fit_reference, fit_brma) {

  reference_values <- list(
    "Residuals"      = as.numeric(stats::residuals(fit_reference)),
    "Rstandard"      = stats::rstandard(fit_reference)[["z"]],
    "Hat values"     = as.numeric(stats::hatvalues(fit_reference)),
    "Cooks distance" = stats::cooks.distance(fit_reference),
    "DFBETAS"        = unlist(stats::dfbetas(fit_reference))
  )
  brma_values <- list(
    "Residuals"      = stats::residuals(
      fit_brma,
      type               = "outcome",
      conditioning_depth = "marginal"
    ),
    "Rstandard"      = suppressWarnings(stats::rstandard(
      fit_brma,
      conditioning_depth = "marginal"
    ))[["z"]],
    "Hat values"     = suppressWarnings(stats::hatvalues(fit_brma)),
    "Cooks distance" = stats::cooks.distance(fit_brma),
    "DFBETAS"        = unlist(stats::dfbetas(fit_brma))
  )

  graphics::par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  for (diagnostic in names(reference_values)) {
    scenario_agreement_plot(
      reference_values[[diagnostic]],
      brma_values[[diagnostic]],
      main = diagnostic
    )
  }

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
  assign(
    "last_scenario_run",
    list(root = root, scenarios = names(files)),
    envir = .scenario_state
  )

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
    review_scenario_snapshots(),
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
