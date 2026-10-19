context("Scenario cleanup regressions")

source(testthat::test_path("..", "scenarios", "helper-scenarios.R"))


test_that("scenario extraction preserves computation and refit failures", {

  local_mocked_s3_method(
    "summary", "scenario_cleanup_fit",
    function(object, ...) data.frame(
      component = "location", parameter = "intercept", Mean = 0.3
    )
  )
  helper_env <- environment(ex_r)
  had_local <- exists("summary_heterogeneity", helper_env, inherits = FALSE)
  if (had_local) old <- get("summary_heterogeneity", helper_env, inherits = FALSE)
  on.exit({
    if (had_local) {
      assign("summary_heterogeneity", old, helper_env)
    } else {
      rm("summary_heterogeneity", envir = helper_env)
    }
  }, add = TRUE)
  assign("summary_heterogeneity", function(object, ...) data.frame(
    component = "study", parameter = "tau", Mean = 0.2
  ), helper_env)
  fit <- structure(list(), class = "scenario_cleanup_fit")
  expect_equal(ex_r(fit, "intercept"), 0.3)
  expect_true(is.na(ex_r(fit, "missing")))
  expect_true(is.na(ex_r(fit, "tau", component = "observation")))

  for (failure_class in c("scenario_computation_error", "BayesTools_refit_required")) {
    failure <- structure(
      list(message = "deliberate heterogeneity failure", call = NULL),
      class = c(failure_class, "error", "condition")
    )
    assign("summary_heterogeneity", function(object, ...) stop(failure), helper_env)
    caught <- tryCatch(ex_r(fit, "intercept"), error = identity)
    expect_identical(caught, failure)
    caught <- tryCatch(ex_r(fit, "missing"), error = identity)
    expect_identical(caught, failure)
  }
})


test_that("empty shared diagnostics fail before graphics setup", {

  helper_env <- environment(plot_marginal_diagnostics)
  old <- get(".scenario_marginal_diagnostic_values", helper_env, inherits = FALSE)
  on.exit(assign(".scenario_marginal_diagnostic_values", old, helper_env), add = TRUE)
  assign(".scenario_marginal_diagnostic_values", function(fit) fit, helper_env)
  failure <- tryCatch(
    plot_marginal_diagnostics(list(a = 1), list(b = 2)),
    error = identity
  )
  expect_s3_class(failure, "error")
  expect_identical(conditionMessage(failure),
                   "'fit_reference' and 'fit_brma' have no shared marginal diagnostics.")
})


test_that("unmatched scenario filters name the argument clearly", {

  files <- c(alpha = "test-alpha.R")
  failure <- tryCatch(.scenario_filter_files(files, "missing"), error = identity)
  expect_s3_class(failure, "error")
  expect_identical(conditionMessage(failure),
                   "No scenarios matched 'filter' = 'missing'.")
  expect_identical(.scenario_filter_files(files, "alpha"), files)
})


test_that("managed missing rendering fails before evaluating its plot", {

  root <- tempfile("robma-scenario-cleanup-render-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  old_runner <- Sys.getenv("ROBMA_SCENARIO_RUNNER", unset = NA_character_)
  on.exit({
    if (is.na(old_runner)) Sys.unsetenv("ROBMA_SCENARIO_RUNNER") else
      Sys.setenv(ROBMA_SCENARIO_RUNNER = old_runner)
  }, add = TRUE)
  Sys.setenv(ROBMA_SCENARIO_RUNNER = "TRUE")
  scenario_start("unit", root = root, show_output = FALSE)
  helper_env <- environment(scenario_plot)
  had_local <- exists("requireNamespace", helper_env, inherits = FALSE)
  if (had_local) old <- get("requireNamespace", helper_env, inherits = FALSE)
  on.exit({
    if (had_local) assign("requireNamespace", old, helper_env) else
      rm("requireNamespace", envir = helper_env)
  }, add = TRUE)
  real_require_namespace <- base::requireNamespace
  missing_dependency <- function(package, ...) {
    if (identical(package, "vdiffr")) return(FALSE)
    real_require_namespace(package, ...)
  }
  assign("requireNamespace", missing_dependency, helper_env)
  # Exercise the installed testthat checker's actual missing-package decision
  # on the old path as well as the helper's new managed-rendering guard.
  missing_check <- testthat::skip_if_not_installed
  environment(missing_check) <- list2env(list(requireNamespace = missing_dependency),
    parent = environment(missing_check))
  testthat::local_mocked_bindings(skip_if_not_installed = missing_check,
    .package = "testthat")
  artifact_state <- function() {
    files <- sort(list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE))
    tools::md5sum(files[!dir.exists(files)])
  }
  before <- artifact_state()
  evaluated <- FALSE
  failure <- tryCatch(scenario_plot("figure", {
    evaluated <- TRUE
  }), error = identity, skip = identity)
  expect_s3_class(failure, "error")
  expect_identical(conditionMessage(failure),
                   "Package 'vdiffr' is required to render managed scenario plots; this run is incomplete.")
  expect_false(evaluated)
  expect_identical(artifact_state(), before)
})


test_that("actual Assink coarse-fit bodies restore the previous sampler on every exit", {

  fit_names <- c("fit_bselmodel.mv_V_cond_coarse", "fit_bselmodel.mv_V_CCI_coarse")
  bodies <- list()
  scenario_file <- testthat::test_path("..", "scenarios", "test-assink2016.R")
  if (!file.exists(scenario_file) && file.exists(testthat::test_path("..", "..", ".git"))) {
    stop("The actual Assink scenario source is missing from the developer checkout.",
         call. = FALSE)
  }
  skip_if_not(file.exists(scenario_file),
    "Actual scenario definitions are excluded from package archives.")
  source_calls <- parse(scenario_file, keep.source = FALSE)
  for (test_call in source_calls) {
    if (!is.call(test_call) ||
        !identical(test_call[[1L]], quote(testthat::test_that))) next
    for (statement in as.list(test_call[[3L]])[-1L]) {
      if (!is.call(statement) ||
          !identical(statement[[1L]], as.name("<-"))) next
      fit_call <- statement[[3L]]
      if (!is.call(fit_call) ||
          !identical(fit_call[[1L]], as.name("scenario_fit"))) next
      if (is.character(fit_call[[2L]]) && fit_call[[2L]] %in% fit_names) {
        bodies[[fit_call[[2L]]]] <- fit_call[[3L]]
      }
    }
  }
  expect_setequal(names(bodies), fit_names)

  # Evaluate the real scenario expression, with only expensive work replaced.
  # No constructor arguments are forced and no model is fitted.
  for (body in bodies) {
    for (previous in c("default", "coarse_corrected")) {
      for (stage in c("success", "fit_error", "marglik_error", "loo_interrupt")) {
        state <- new.env(parent = emptyenv())
        state$sampler <- previous
        state$observed <- NULL
        env <- new.env(parent = baseenv())
        env$RoBMA.options <- function(...) {

          requested_options <- list(...)
          if (length(requested_options)) state$sampler <- requested_options[["selection.sampler"]]
          return(list(selection.sampler = state$sampler))
        }
        env$bselmodel.mv <- function(...) {

          state$observed <- state$sampler
          if (stage == "fit_error") stop("deliberate fit failure")
          return(list(marker = "fit"))
        }
        env$add_marglik <- function(fit) {

          if (stage == "marglik_error") stop("deliberate marginal-likelihood failure")
          return(fit)
        }
        env$add_loo <- function(fit) {

          if (stage == "loo_interrupt") {
            stop(structure(list(message = "deliberate interrupt", call = NULL),
                           class = c("interrupt", "condition")))
          }
          return(fit)
        }
        evaluate <- eval(call("function", pairlist(), body), env)
        result <- tryCatch(evaluate(), error = identity, interrupt = identity)
        expect_identical(state$observed, "coarse_corrected")
        expect_identical(state$sampler, previous)
        if (stage == "success") {
          expect_identical(result, list(marker = "fit"))
        } else if (stage == "loo_interrupt") {
          expect_s3_class(result, "interrupt")
        } else {
          expect_s3_class(result, "error")
        }
      }
    }
  }
})
