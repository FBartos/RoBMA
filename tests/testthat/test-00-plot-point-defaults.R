context("Plot point defaults")

.expect_metafor_point_style <- function(dots) {

  expect_identical(dots[["pch"]], 21)
  expect_identical(dots[["col"]], "black")
  expect_identical(dots[["bg"]], "#A6A6A6")
  expect_identical(dots[["cex"]], 1)
  expect_identical(dots[["size"]], 2)
}

test_that("point plotting defaults use metafor-style markers", {

  .expect_metafor_point_style(.plot_point_style_defaults(list()))
  .expect_metafor_point_style(.set_dots_funnel(list()))
  .expect_metafor_point_style(.set_dots_qqnorm(list()))
  .expect_metafor_point_style(.set_dots_radial())
  .expect_metafor_point_style(.set_dots_regplot())
})

test_that("point plotting defaults respect overrides", {

  custom <- list(
    pch  = 19,
    col  = "blue",
    bg   = "lightblue",
    cex  = 1.5,
    size = 3
  )

  expect_equal(.plot_point_style_defaults(custom)[names(custom)], custom)
  expect_equal(.set_dots_funnel(custom)[names(custom)], custom)
  expect_equal(.set_dots_qqnorm(custom)[names(custom)], custom)
  expect_equal(do.call(.set_dots_radial, custom)[names(custom)], custom)
  expect_equal(do.call(.set_dots_regplot, custom)[names(custom)], custom)
})
