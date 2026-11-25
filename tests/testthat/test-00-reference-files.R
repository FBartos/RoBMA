context("Human-reviewed text and table references")

source(testthat::test_path("common-functions.R"))


test_that("changed reference text is retained as a review candidate", {

  root <- tempfile("robma-reference-test-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(root, recursive = TRUE)
  reference <- file.path(root, "output.txt")
  candidate <- file.path(root, "output.new.txt")
  writeLines("old output", reference)

  expect_failure(
    test_reference_text(
      text      = "new output",
      filename  = "output.txt",
      print_dir = root
    ),
    "Reference candidate cached"
  )
  expect_equal(readLines(reference, warn = FALSE), "old output")
  expect_equal(readLines(candidate, warn = FALSE), "new output")

  expect_no_error(test_reference_text(
    text      = "old output",
    filename  = "output.txt",
    print_dir = root
  ))
  expect_false(file.exists(candidate))
})


test_that("changed reference tables use the same candidate workflow", {

  root <- tempfile("robma-reference-test-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(root, recursive = TRUE)
  table <- data.frame(value = 1:2)
  reference <- file.path(root, "table.txt")
  candidate <- file.path(root, "table.new.txt")
  writeLines("old table", reference)

  expect_failure(
    test_reference_table(
      table     = table,
      filename  = "table.txt",
      print_dir = root
    ),
    "Reference candidate cached"
  )
  expect_equal(readLines(reference, warn = FALSE), "old table")
  expect_equal(
    readLines(candidate, warn = FALSE),
    capture_output_lines(table, print = TRUE, width = 150)
  )
})
