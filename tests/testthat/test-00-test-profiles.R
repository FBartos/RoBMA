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
  iwmde_names    <- certification_case_fit_names("iwmde-qcmde")
  expect_true(all(smoke_group %in% extended_names))
  expect_true(all(smoke_group %in% iwmde_names))
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
