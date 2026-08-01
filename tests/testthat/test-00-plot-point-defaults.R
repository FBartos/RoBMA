context("Plot point defaults")

test_that("regplot point sizes use relative precision without an offset", {

  for (scale in c(1e-200, 1e200)) {
    expect_equal(
      .regplot_normalized_precision(scale * c(1, sqrt(2), 2)),
      c(1, 1 / 3, 0)
    )
  }

  expect_identical(.regplot_normalized_precision(rep(3, 4)), numeric(4))
})

test_that("plot affine transforms retain non-zero represented slopes", {

  raw <- 2^-400 * c(-2, -1, 1, 2)

  large <- .plot_brma_affine_sample_transform(raw, raw * 2^200)
  small <- .plot_brma_affine_sample_transform(raw / 2^-400, raw)

  expect_equal(large[["slope"]], 2^200)
  expect_equal(small[["slope"]], 2^-400)
})

test_that("regplot variance conversion rejects invalid values", {

  expect_equal(
    .regplot_variance_sd(matrix(c(0, 1, 4)), "Test variance"),
    matrix(c(0, 1, 2))
  )
  expect_error(
    .regplot_variance_sd(matrix(c(1, -1)), "Test variance"),
    "finite and non-negative"
  )
})

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
