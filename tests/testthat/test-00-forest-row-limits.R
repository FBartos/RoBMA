test_that("forest automatic limits retain custom study and summary rows", {

  for (rows in list(c(-10, -12, -14), -10)) {
    ylim <- .forest_default_ylim(
      k = 3, row = -1, dots = list(rows = rows), predstyle = "line"
    )
    study_rows <- if (length(rows) == 1L) rows - 0:2 else rows
    expect_true(all(study_rows > ylim[1L] & study_rows < ylim[2L]))
  }

  for (style in c("line", "bar", "shade", "dist")) {
    ylim <- .forest_default_ylim(k = 3, row = 20, dots = list(), predstyle = style)
    expect_true(20 > ylim[1L] && 20 < ylim[2L])
    expect_true(all(1:3 > ylim[1L] & 1:3 < ylim[2L]))
  }

  expect_equal(.forest_default_ylim(3, -1, list(), "line"), c(-2, 6))
  expect_equal(.forest_default_ylim(3, -1, list(), "bar"), c(-3, 6))
})
