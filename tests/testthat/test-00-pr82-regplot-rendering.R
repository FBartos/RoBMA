test_that("regplot clips bands without discarding vertices", {

  testthat::skip_if_not_installed("ggplot2")
  band <- data.frame(
    x = c(-2, 0, 2, 2, 0, -2), y = c(-1, 0, 1, 3, 2, 1),
    lower = c(-1, 0, 1, 1, 0, -1), upper = c(1, 2, 3, 3, 2, 1),
    group = "all", group_id = 1
  )
  data <- list(
    points = data.frame(x = 0, y = 1, size = 1, group = "all", group_id = 1),
    pred = NULL, ci = band, pi = NULL, si = NULL, refline = NULL,
    mod_type = "continuous", groups = "all", xlim = c(-1, 1),
    ylim = c(0, 2), xlab = "x", ylab = "y"
  )
  plot <- .regplot_plot_ggplot(data, .set_dots_regplot())
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[1L]]$x, band$x)
  expect_equal(built$data[[1L]]$y, band$y)
  expect_identical(plot$coordinates$limits$x, data$xlim)
  expect_identical(plot$coordinates$limits$y, data$ylim)
})

test_that("regplot supports zero standard errors in cosmetic point sizes", {

  expect_equal(.regplot_normalized_precision(c(0, 0.1, 0.2)), c(1, 0, 0))
  expect_equal(.regplot_normalized_precision(c(0, 0)), c(0, 0))
  expect_error(.regplot_normalized_precision(c(-0.1, 0.2)), "non-negative", fixed = TRUE)
})

test_that("explicit grouped line colors override the automatic palette", {

  dots <- .set_dots_regplot(lcol = "red")
  expect_identical(
    .regplot_palette(c("a", "b"), dots$lcol, attr(dots, "lcol_supplied")),
    c(a = "red", b = "red")
  )
  defaults <- .set_dots_regplot()
  expect_length(unique(.regplot_palette(c("a", "b"), defaults$lcol,
                                         attr(defaults, "lcol_supplied"))), 2L)
})

test_that("regplot jitter preserves the caller's random stream", {

  withr::local_seed(200)
  before <- .Random.seed
  kind <- RNGkind()
  first <- .regplot_jitter_values(10, 0.2)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  expect_identical(.regplot_jitter_values(10, 0.2), first)
  rm(".Random.seed", envir = .GlobalEnv)
  .regplot_jitter_values(10, 0.2)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
