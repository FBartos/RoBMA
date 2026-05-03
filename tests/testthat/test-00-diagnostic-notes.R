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
