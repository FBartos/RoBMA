context("Session cache of hypothesis plans on fitted models")

source(testthat::test_path("common-functions.R"))


# Fitted models with a location, a scale or a moderator component.
.plan_cache_fit_names <- function() {

  utils::head(intersect(
    c("dat.lehmann2018_RoBMA_mods", "konstantopoulos2011_3lvl",
      "dat.lehmann2018-3PSM"),
    list_fits()
  ), 1L)
}


# The number of plans `expr` builds.
.plan_cache_builds <- function(expr) {

  builds   <- 0L
  original <- .hypothesis_plan
  value    <- testthat::with_mocked_bindings(
    expr,
    .hypothesis_plan = function(...) {
      builds <<- builds + 1L
      original(...)
    },
    .package = "RoBMA"
  )

  list(value = value, builds = builds)
}


test_that("hypothesis() and hypothesis_quantities() build the plan of a statement once per fit", {

  fit_names <- .plan_cache_fit_names()
  skip_if(length(fit_names) == 0L, "No cached fitted model is active.")

  withr::defer(.hypothesis_plan_cache_clear())
  for (fit_name in fit_names) {
    fit <- load_fit(fit_name, validate = FALSE)
    .hypothesis_plan_cache_clear()
    metadata <- .brma_parameter_catalog_metadata(fit)
    entries  <- metadata[["entries"]]
    entries  <- entries[entries[["component"]] != "bias", , drop = FALSE]
    entry    <- as.list(entries[1L, setdiff(names(entries), "aliases"),
                                drop = FALSE])
    root      <- .hypothesis_quantities_reference_root(metadata, entry)
    statement <- paste0("`", root, "` = 0")
    plan_of   <- function(...) {
      .hypothesis_plans(fit, statement, component = entry[["component"]],
                        metadata = metadata, ...)[[1L]]
    }

    first <- .plan_cache_builds(plan_of())
    expect_identical(first$builds, 1L, info = fit_name)
    again <- .plan_cache_builds(plan_of())
    expect_identical(again$builds, 0L, info = fit_name)
    expect_identical(again$value, first$value, info = fit_name)
    # Another argument of the plan is another plan of the statement.
    expect_identical(.plan_cache_builds(plan_of(n_samples = 500L))$builds, 1L,
                     info = fit_name)
    expect_identical(.plan_cache_builds(plan_of(n_samples = 500L))$builds, 0L,
                     info = fit_name)
    expect_identical(.plan_cache_builds(plan_of(standardized = TRUE))$builds, 1L,
                     info = fit_name)

    # hypothesis_quantities() plans every statement of the fit once, the
    # statement above among them, and once again for the next call.
    .hypothesis_plan_cache_clear()
    quantities <- .plan_cache_builds(suppressWarnings(hypothesis_quantities(fit)))
    expect_gt(quantities$builds, 1L)
    expect_identical(.plan_cache_builds(plan_of())$builds, 0L, info = fit_name)
    repeated <- .plan_cache_builds(suppressWarnings(hypothesis_quantities(fit)))
    expect_identical(repeated$builds, 0L, info = fit_name)
    expect_identical(repeated$value, quantities$value, info = fit_name)

    # And so does hypothesis(), which reads the plan the method above kept
    # (whether it evaluates the statement or stops with the plan's refusal).
    call <- .plan_cache_builds(suppressWarnings(tryCatch(
      hypothesis(fit, statement, component = entry[["component"]],
                 density_method = "KDE"),
      error = function(e) e
    )))
    expect_identical(call$builds, 0L, info = fit_name)
  }
})


test_that("the results of a cache are those of a call without one", {

  fit_names <- .plan_cache_fit_names()
  skip_if(length(fit_names) == 0L, "No cached fitted model is active.")

  withr::defer(.hypothesis_plan_cache_clear())
  for (fit_name in fit_names) {
    fit <- load_fit(fit_name, validate = FALSE)
    metadata <- .brma_parameter_catalog_metadata(fit)
    entries  <- metadata[["entries"]]
    entries  <- utils::head(entries[entries[["component"]] != "bias", ,
                                    drop = FALSE], 3L)
    run <- function(cached) {
      .hypothesis_plan_cache_clear()
      testthat::with_mocked_bindings(
        {
          # Each statement twice: the second time from the cache.
          statements <- lapply(seq_len(nrow(entries)), function(i) {
            entry     <- as.list(entries[i, setdiff(names(entries), "aliases"),
                                         drop = FALSE])
            statement <- paste0(
              "`", .hypothesis_quantities_reference_root(metadata, entry),
              "` > 0"
            )
            call <- function() {
              suppressWarnings(tryCatch(
                hypothesis(fit, statement, component = entry[["component"]],
                           density_method = "KDE", seed = 1),
                error = function(e) e
              ))
            }
            list(call(), call())
          })
          list(statements = statements,
               quantities = suppressWarnings(hypothesis_quantities(fit)))
        },
        .hypothesis_plan_cache = if (cached) {
          .hypothesis_plan_cache
        } else {
          function(object = NULL) new.env(parent = emptyenv())
        },
        .package = "RoBMA"
      )
    }
    cached   <- run(TRUE)
    uncached <- run(FALSE)
    expect_identical(cached, uncached, info = fit_name)
  }
})
