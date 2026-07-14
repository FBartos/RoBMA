context("fit cache source hash portability")

source(testthat::test_path("common-functions.R"))


test_that("native source hashes ignore line-ending representation", {

  file_lf      <- tempfile("robma-native-lf-", fileext = ".cc")
  file_crlf    <- tempfile("robma-native-crlf-", fileext = ".cc")
  file_changed <- tempfile("robma-native-changed-", fileext = ".cc")
  on.exit(unlink(c(file_lf, file_crlf, file_changed)), add = TRUE)

  writeBin(charToRaw("int value() { return 1; }\n"), file_lf)
  writeBin(charToRaw("int value() { return 1; }\r\n"), file_crlf)
  writeBin(charToRaw("int value() { return 2; }\n"), file_changed)

  expect_identical(
    .fit_cache_source_file_md5(file_lf),
    .fit_cache_source_file_md5(file_crlf)
  )
  expect_false(identical(
    .fit_cache_source_file_md5(file_lf),
    .fit_cache_source_file_md5(file_changed)
  ))
})
