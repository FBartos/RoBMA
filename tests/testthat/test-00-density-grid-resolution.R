test_that("adaptive density grids stop at the floating-point resolution limit", {

  bounds <- c(1, 1 + 4 * .Machine$double.eps)
  expect_error(
    .iwmde_complete_display_grid(bounds, bounds, n_points = 20L),
    "too few representable values"
  )
  grid <- .iwmde_complete_display_grid(bounds, bounds, n_points = 5L)
  expect_identical(grid, 1 + (0:4) * .Machine$double.eps)
})
