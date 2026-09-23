context("Diagnostic notes")

test_that("DFBETAS note objects print and preserve data-frame contract", {

  x <- data.frame(mu = c(NaN, NaN))
  x <- RoBMA:::.diagnostic_with_note(
    x,
    class = "dfbetas.brma",
    note  = RoBMA:::.diagnostic_zero_variance_note(
      diagnostic = "DFBETAS",
      parameters = "mu"
    )
  )

  expect_s3_class(x, "data.frame")
  expect_s3_class(x, "dfbetas.brma")
  expect_true(all(is.nan(x[["mu"]])))
  expect_true(any(grepl("Note: DFBETAS could not be computed", capture.output(print(x)))))
})


test_that("COVRATIO note objects print and preserve numeric values", {

  x <- RoBMA:::.diagnostic_with_note(
    rep(NaN, 2),
    class = "covratio.brma",
    note  = RoBMA:::.diagnostic_zero_variance_note(
      diagnostic = "COVRATIO",
      parameters = "mu",
      variance   = "posterior"
    )
  )

  expect_s3_class(x, "covratio.brma")
  expect_true(all(is.nan(as.numeric(x))))
  expect_true(any(grepl("Note: COVRATIO could not be computed", capture.output(print(x)))))
})


test_that("COVRATIO exclusion notes print and preserve numeric values", {

  x <- RoBMA:::.diagnostic_with_note(
    c(.9, 1.1),
    class = "covratio.brma",
    note  = RoBMA:::.diagnostic_excluded_zero_variance_note(
      diagnostic = "COVRATIO",
      parameters = "rho",
      variance   = "posterior"
    )
  )

  expect_s3_class(x, "covratio.brma")
  expect_equal(as.numeric(x), c(.9, 1.1))
  expect_true(any(grepl("Note: COVRATIO excluded parameter", capture.output(print(x)))))
})


test_that("Diagnostic label helpers preserve study labels and repair invalid row names", {

  object <- brma.norm(
    yi        = c(.10, .20, .30),
    sei       = c(.10, .20, .30),
    slab      = c("Study A", "Study A", "Study C"),
    only_data = TRUE
  )
  object[["data"]][["outcome"]][["slab"]][3] <- NA

  labels <- RoBMA:::.diagnostic_study_labels(object)

  expect_equal(labels, c("Study A", "Study A.1", "3"))

  out <- RoBMA:::.diagnostic_set_rownames(data.frame(x = 1:3), object)
  expect_equal(rownames(out), labels)

  vec <- RoBMA:::.diagnostic_set_names(1:3, object)
  expect_equal(names(vec), labels)
})


test_that("diagnostic plot wrappers expose the full selector API", {

  wrapper_formals <- list(
    autocorrelation = names(formals(plot_diagnostic_autocorrelation.brma)),
    trace           = names(formals(plot_diagnostic_trace.brma)),
    density         = names(formals(plot_diagnostic_density.brma))
  )

  for (formals_i in wrapper_formals) {
    expect_true(all(c(
      "parameter_mods", "parameter_scale", "component", "type", "lags"
    ) %in% formals_i))
  }
})
