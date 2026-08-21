context("Test profiles")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-profile-cases.R"))


.with_test_profile_env <- function(code) {

  variables <- c("ROBMA_TEST_PROFILE", "ROBMA_TEST_ACTIVE_FITS")
  old <- Sys.getenv(variables, unset = NA_character_)
  on.exit({
    for (variable in variables) {
      if (is.na(old[[variable]])) {
        Sys.unsetenv(variable)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old[[variable]]), variable))
      }
    }
  }, add = TRUE)

  force(code)
}


test_that("quiet LLM reporter counts skips without printing each reason", {

  output_path <- tempfile("robma-quiet-reporter-", fileext = ".txt")
  on.exit(unlink(output_path), add = TRUE)
  reporter <- quiet_llm_reporter(file = output_path)
  skip <- structure(
    list(message = "routine certification skip"),
    class = c("expectation_skip", "expectation", "condition")
  )
  success <- structure(
    list(),
    class = c("expectation_success", "expectation", "condition")
  )
  failure <- tryCatch(
    testthat::expectation("failure", "deliberate failure"),
    error = identity
  )

  reporter$add_result("profiles", "quiet output", skip)
  reporter$add_result("profiles", "quiet output", success)
  reporter$add_result("profiles", "quiet output", failure)
  reporter$end_reporter()
  output <- readLines(output_path, warn = FALSE)

  expect_equal(reporter$n_skip, 1L)
  expect_equal(reporter$n_ok, 1L)
  expect_equal(reporter$n_fail, 1L)
  expect_true(any(grepl("deliberate failure", output, fixed = TRUE)))
  expect_true(any(grepl("SKIP 1", output, fixed = TRUE)))
  expect_false(any(grepl("routine certification skip", output, fixed = TRUE)))
})


test_that("interactive test runner dispatches filtered and comprehensive profiles", {

  runner_env <- new.env(parent = globalenv())
  source(
    testthat::test_path("..", "..", ".dev", "test-tests.R"),
    local = runner_env
  )
  calls <- list()
  assign(
    ".run_robma_test_profile",
    function(profile, clean = FALSE, filter = NULL) {

      calls[[length(calls) + 1L]] <<- list(
        profile = profile,
        clean   = clean,
        filter  = filter
      )
      return(invisible(TRUE))
    },
    envir = runner_env
  )
  assign(
    "review_test_snapshots",
    function(...) invisible(data.frame()),
    envir = runner_env
  )

  runner_env$test_tests(
    filter          = "interpret",
    refit           = TRUE,
    load_package    = FALSE,
    stop_on_failure = TRUE,
    root            = testthat::test_path()
  )
  expect_identical(
    vapply(calls, `[[`, character(1), "profile"),
    c("refresh-standard", "filter")
  )
  expect_identical(vapply(calls, `[[`, logical(1), "clean"), c(TRUE, FALSE))
  expect_identical(calls[[2L]][["filter"]], "interpret")

  calls <- list()
  runner_env$test_tests(
    filter       = "input-data",
    load_package = FALSE,
    root         = testthat::test_path()
  )
  expect_identical(
    vapply(calls, `[[`, character(1), "profile"),
    "filter"
  )
  expect_false(calls[[1L]][["clean"]])

  calls <- list()
  runner_env$test_tests(
    load_package = FALSE,
    root         = testthat::test_path()
  )
  expect_identical(
    vapply(calls, `[[`, character(1), "profile"),
    c("refresh-standard", "standard", "certification")
  )
  expect_false(any(vapply(calls, `[[`, logical(1), "clean")))
  expect_error(
    runner_env$test_tests(
      update_timings = TRUE,
      load_package   = FALSE,
      root           = testthat::test_path()
    ),
    "Ordinary tests have no timing baselines"
  )
})


test_that("certification cases partition expensive evidence", {

  cases   <- certification_cases()
  catalog <- fit_catalog()

  expect_identical(CERTIFICATION_CASE_TIMEOUT_SECONDS, 3600)
  expect_identical(names(cases), c(
    "numerical-kernels",
    "normal-models",
    "glmm-models",
    "multivariate-core",
    "multivariate-extended",
    "multivariate-singular",
    "multivariate-parity-cs",
    "multivariate-parity-nested",
    "multivariate-parity-har",
    "multivariate-parity-treatment",
    "iwmde-qcmde"
  ))

  test_files <- list.files(
    testthat::test_path(),
    pattern    = "^test-.*\\.[Rr]$",
    full.names = FALSE
  )
  test_stems <- sub("^test-", "", sub("\\.[Rr]$", "", test_files))
  for (name in names(cases)) {
    case <- cases[[name]]
    expect_true(nzchar(case[["description"]]), info = name)
    expect_true(nzchar(case[["test_filter"]]), info = name)
    expect_true(any(grepl(case[["test_filter"]], test_stems)), info = name)
    expect_true(all(case[["fit_sources"]] %in% test_files), info = name)
    required_tests <- case[["required_tests"]]
    expect_s3_class(required_tests, "data.frame")
    expect_gt(nrow(required_tests), 0L)
    expect_identical(names(required_tests), c("file", "test"), info = name)
    expect_false(anyDuplicated(required_tests) > 0L, info = name)
    expect_true(all(required_tests[["file"]] %in% test_files), info = name)
    for (i in seq_len(nrow(required_tests))) {
      required_file <- required_tests[["file"]][[i]]
      required_test <- required_tests[["test"]][[i]]
      required_stem <- sub(
        "^test-", "", sub("\\.[Rr]$", "", required_file)
      )
      expect_true(grepl(case[["test_filter"]], required_stem), info = name)
      source_lines <- readLines(
        testthat::test_path(required_file),
        warn = FALSE
      )
      declaration <- paste0("test_that(\"", required_test, "\"")
      expect_true(
        any(grepl(declaration, source_lines, fixed = TRUE)),
        info = paste(name, required_file, required_test)
      )
    }
    fit_names <- certification_case_fit_names(name, catalog = catalog)
    fit_sources <- catalog[["source_file"]][match(fit_names, catalog[["name"]])]
    expect_true(all(fit_sources %in% case[["fit_sources"]]), info = name)
  }

  assigned <- unique(unlist(
    lapply(
      names(cases),
      certification_case_fit_names,
      catalog = catalog
    ),
    use.names = FALSE
  ))
  certification_fits <- catalog[["name"]][
    catalog[["profile"]] == "certification"
  ]
  expect_setequal(certification_fits, intersect(certification_fits, assigned))

  smoke_group <- c(
    "brma.mv_latent",
    "brma.mv_whitened",
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_random"
  )
  extended_names <- certification_case_fit_names("multivariate-extended")
  glmm_names     <- certification_case_fit_names("glmm-models")
  iwmde_names    <- certification_case_fit_names("iwmde-qcmde")
  expect_true(all(smoke_group %in% extended_names))
  expect_true(all(smoke_group %in% iwmde_names))
  expect_true(all(c(
    "nielweise2008_glmm_effect_null",
    "nielweise2008_glmm",
    "dat.lehmann2018-3PSM_effect_null",
    "dat.lehmann2018-3PSM"
  ) %in% iwmde_names))
  expect_true("brma.mv_block_mvn_random_scale" %in% extended_names)
  expect_true("bcg_BMA.glmm_3lvl_location_scale" %in% glmm_names)
})


test_that("certification evidence requires one executed passing result", {

  required <- .required_tests("test-required.R", "required evidence")
  successful <- data.frame(
    file    = c("test-required.R", "test-unrelated.R"),
    test    = c("required evidence", "expected unrelated skip"),
    skipped = c(FALSE, TRUE),
    passed  = c(2L, 0L),
    stringsAsFactors = FALSE
  )

  expect_true(isTRUE(validate_certification_evidence(
    successful,
    required,
    "synthetic"
  )))
  expect_error(
    validate_certification_evidence(
      successful[0L, , drop = FALSE],
      required,
      "synthetic"
    ),
    "observed 0"
  )

  renamed <- successful
  renamed[["test"]][[1L]] <- "renamed evidence"
  expect_error(
    validate_certification_evidence(renamed, required, "synthetic"),
    "observed 0"
  )

  skipped <- successful
  skipped[["skipped"]][[1L]] <- TRUE
  expect_error(
    validate_certification_evidence(skipped, required, "synthetic"),
    "skipped required test"
  )

  zero <- successful
  zero[["passed"]][[1L]] <- 0L
  expect_error(
    validate_certification_evidence(zero, required, "synthetic"),
    "no passing expectations"
  )

  duplicated <- rbind(successful[1L, ], successful[1L, ])
  expect_error(
    validate_certification_evidence(duplicated, required, "synthetic"),
    "observed 2"
  )
})


test_that("active fit selection is explicit and fails closed", {

  .with_test_profile_env({
    Sys.setenv(ROBMA_TEST_PROFILE = "certification")
    Sys.unsetenv("ROBMA_TEST_ACTIVE_FITS")
    all_fits <- active_fit_catalog()[["name"]]
    selected <- all_fits[1:2]

    Sys.setenv(ROBMA_TEST_ACTIVE_FITS = paste(selected, collapse = ","))
    expect_identical(active_fit_catalog()[["name"]], selected)

    Sys.setenv(ROBMA_TEST_ACTIVE_FITS = "__none__")
    expect_equal(nrow(active_fit_catalog()), 0L)

    Sys.setenv(ROBMA_TEST_ACTIVE_FITS = "not-a-catalog-fit")
    expect_error(active_fit_catalog(), "Unknown active cached fit")
  })
})


test_that("certification case fit filters select their source files", {

  for (name in certification_case_names()) {
    case       <- certification_case(name)
    fit_filter <- certification_case_fit_filter(name)
    if (length(case[["fit_sources"]]) == 0L) {
      expect_null(fit_filter, info = name)
    } else {
      stems <- sub("^test-", "", sub("\\.[Rr]$", "", case[["fit_sources"]]))
      expect_true(all(grepl(fit_filter, stems)), info = name)
    }
  }
})


test_that("fit-source prerequisites may use inactive cached fits", {

  helper_env    <- environment(skip_if_missing_fits)
  old_list_fits <- get("list_fits", envir = helper_env, inherits = TRUE)
  observed      <- new.env(parent = emptyenv())
  on.exit(assign("list_fits", old_list_fits, envir = helper_env), add = TRUE)
  assign(
    "list_fits",
    function(validate = TRUE, active_only = TRUE, ...) {
      observed[["validate"]]    <- validate
      observed[["active_only"]] <- active_only
      "dependency"
    },
    envir = helper_env
  )

  expect_no_error(skip_if_missing_fits("dependency", active_only = FALSE))
  expect_true(observed[["validate"]])
  expect_false(observed[["active_only"]])
})


test_that("raw lazy caches expose only active-profile fits", {

  helper_env    <- environment(lazy_fits)
  old_list_fits <- get("list_fits", envir = helper_env, inherits = TRUE)
  observed      <- list()
  on.exit(assign("list_fits", old_list_fits, envir = helper_env), add = TRUE)
  assign(
    "list_fits",
    function(validate = TRUE, active_only = TRUE, ...) {
      observed[[length(observed) + 1L]] <<- c(
        validate    = validate,
        active_only = active_only
      )
      "active"
    },
    envir = helper_env
  )

  fits  <- lazy_fits(c("active", "inactive"), validate = FALSE)
  infos <- lazy_infos(c("inactive", "active"), validate = FALSE)

  expect_identical(names(fits), "active")
  expect_identical(names(infos), "active")
  expect_length(observed, 2L)
  for (call in observed) {
    expect_false(call[["validate"]])
    expect_true(call[["active_only"]])
  }
})
