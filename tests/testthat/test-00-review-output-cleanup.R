test_that("native thread option failures name NA and preserve atomic updates", {

  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  for (value in list(0, -1, .5, Inf, c(1, 2), "two")) {
    expect_error(RoBMA.options(max_cores = 2L, native_threads = value),
      "Option 'native_threads' must be NA or a positive integer.", fixed = TRUE)
    expect_identical(RoBMA.options(), previous)
  }
  expect_silent(RoBMA.options(native_threads = NA_integer_))
  expect_identical(RoBMA.get_option("native_threads"), NA_integer_)
})
