context("CRAN smoke test")

# The only tests CRAN runs (tests/testthat.R): the package and its native code
# load, every registered method resolves, and the JAGS module compiles and
# samples. Full behavior is covered by the test profiles.

on_cran <- get("on_cran", envir = asNamespace("testthat"), inherits = FALSE)

skip_if_not(on_cran(), "CRAN smoke test only")

test_that("native routines and S3 methods are registered", {

  expect_true(isTRUE(validate_native_test_symbols()))

  # columns: generic, class, method, and the package of a foreign generic;
  # methods for generics of suggested packages that are not installed (such as
  # metafor in CRAN's check without suggested packages) are skipped
  methods  <- getNamespaceInfo("RoBMA", "S3methods")
  owners   <- unique(stats::na.omit(methods[, 4L]))
  missing_owners <- owners[!vapply(owners, requireNamespace, logical(1), quietly = TRUE)]
  methods  <- methods[!methods[, 4L] %in% missing_owners, , drop = FALSE]
  resolved <- vapply(seq_len(nrow(methods)), function(i) {
    owner <- if (is.na(methods[i, 4L])) "RoBMA" else methods[i, 4L]
    !is.null(getS3method(methods[i, 1L], methods[i, 2L], optional = TRUE,
                         envir = asNamespace(owner)))
  }, logical(1))
  expect_identical(methods[!resolved, 3L], character())

  # exports include re-exported BayesTools functions
  exports <- getNamespaceExports("RoBMA")
  missing <- exports[!vapply(exports, function(name) {
    !inherits(try(getExportedValue("RoBMA", name), silent = TRUE), "try-error")
  }, logical(1))]
  expect_identical(missing, character())
})

test_that("package loads and fits minimal brma and selection models", {

  yi  <- c(-0.12, 0.05, 0.18, 0.31)
  sei <- c(0.20, 0.18, 0.22, 0.19)
  fit_settings <- list(
    yi = yi, sei = sei, measure = "GEN",
    prior_unit_information_sd = 1,
    sample = 100, burnin = 50, adapt = 500,
    chains = 1, seed = 1, silent = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )

  fit <- do.call(brma, fit_settings)
  expect_s3_class(fit, "brma.norm")
  expect_s3_class(summary(fit), "summary.brma")

  # the selection likelihood is a distribution of the RoBMA JAGS module
  fit_selection <- do.call(bselmodel, c(fit_settings, list(
    prior_bias = prior_weightfunction("one-sided", .025, wf_fixed(c(1, .5)))
  )))
  expect_s3_class(fit_selection, "bselmodel")
  expect_s3_class(summary(fit_selection), "summary.brma")
  expect_true("RoBMA" %in% rjags::list.modules())
  expect_match(paste(fit_selection$fit$model, collapse = "\n"), "dselnorm_", fixed = TRUE)
})
