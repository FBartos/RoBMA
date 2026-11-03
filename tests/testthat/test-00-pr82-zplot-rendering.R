test_that("zplot respects explicit vertical limits in both renderers", {

  object <- list(zplot = list(data = list(z = rep(0.1, 10))), priors = NULL)
  for (backend in c("base", "ggplot")) {
    dots <- .get_dots_hist_zplot(list(ylim = c(0, 0.2)), backend, max_density = 2)
    expect_equal(dots$ylim, c(0, 0.2))
  }
  expect_equal(.get_dots_hist_zplot(list(.zplot_auto_ymax = 3), "base", 2)$ylim,
               c(0, 3))
  testthat::skip_if_not_installed("ggplot2")
  plot <- plot.zplot_brma(object, plot_type = "ggplot", plot_fit = FALSE,
                          plot_thresholds = FALSE, ylim = c(0, 0.2))
  expect_equal(plot$coordinates$limits$y, c(0, 0.2))
  expect_gt(max(ggplot2::ggplot_build(plot)$data[[1L]]$y), 0.2)
})

test_that("specific zplot line arguments override shared graphical arguments", {

  testthat::skip_if_not_installed("ggplot2")
  seen_alpha <- NULL
  testthat::local_mocked_bindings(
    lines.zplot_brma = function(..., alpha) {
      seen_alpha <<- alpha
      data.frame(x = c(-1, 1), y = c(0.1, 0.2), y_lCI = 0, y_uCI = 0.3)
    },
    .package = "RoBMA"
  )
  object <- list(zplot = list(data = list(z = c(-1, 0, 1))), priors = NULL)
  plot <- plot.zplot_brma(
    object, plot_type = "ggplot", plot_ci = FALSE, plot_thresholds = FALSE,
    col = "black", alpha = 0.4, dots_fit = list(col = "red", alpha = 0.25)
  )
  expect_equal(seen_alpha, 0.25)
  expect_identical(plot$layers[[2L]]$aes_params$colour, "red")
})

test_that("coarse histogram bins retain every selection threshold", {

  prior <- BayesTools::prior_weightfunction(
    "one-sided", c(0.05, 0.1), weights = BayesTools::wf_fixed(c(1, 0.7, 0.5))
  )
  priors <- list(outcome = list(bias = prior))
  thresholds <- stats::qnorm(c(0.05, 0.1), lower.tail = FALSE)
  bins <- .zplot_bins(from = 0, to = 3, by = 1.5, length.out = NULL,
    thresholds = .zplot_threshold(priors))
  expect_true(all(thresholds %in% bins))
  expect_identical(range(bins), c(0, 3))
  expect_true(all(diff(bins) > 0))
})
test_that("zplot display grids share the signed selection cutoffs", {

  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    weights = BayesTools::wf_fixed(c(1, .2)))
  priors <- list(outcome = list(bias = prior))
  positive <- .zplot_threshold(priors, "positive")
  negative <- .zplot_threshold(priors, "negative")
  expect_equal(negative, -positive, tolerance = 1e-14)
  for (thresholds in list(positive, negative)) {
    histogram <- .zplot_bins(-6, 6, .5, NULL,
      type = "hist", thresholds = thresholds)
    density <- .zplot_bins(-6, 6, .05, NULL,
      type = "dens", thresholds = thresholds)
    expect_true(all(thresholds %in% histogram))
    expect_true(all(thresholds %in% density))
  }
  prior <- BayesTools::prior_weightfunction("two-sided", .05,
    weights = BayesTools::wf_fixed(c(1, .2)))
  priors <- list(outcome = list(bias = prior))
  expect_equal(.zplot_threshold(priors, "positive"),
    .zplot_threshold(priors, "negative"), tolerance = 1e-14)
  expect_identical(.zplot_threshold(NULL, "negative"), numeric())
})

test_that("zplot histograms of an empty plotting range have zero density", {

  object <- list(zplot = list(data = list(z = c(-1, 0, 1))), priors = NULL)
  expect_message(
    data <- hist.zplot_brma(object, from = 20, to = 25, as_data = TRUE),
    "3 z-statistics are out of the plotting range"
  )
  expect_identical(data$density, rep(0, nrow(data)))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_message(hist.zplot_brma(object, from = 20, to = 25), "out of the plotting range")
  expect_equal(graphics::par("usr")[3:4], c(0, 1) + c(-0.04, 0.04))
})

test_that("zplot histograms honour explicit horizontal limits", {

  object <- list(zplot = list(data = list(z = c(-1, 0, 1))), priors = NULL)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  hist.zplot_brma(object, xlim = c(-2, 2))
  expect_equal(graphics::par("usr")[1:2], c(-2, 2) + c(-0.16, 0.16))
  plot.zplot_brma(object, plot_fit = FALSE, plot_thresholds = FALSE, xlim = c(-3, 3))
  expect_equal(graphics::par("usr")[1:2], c(-3, 3) + c(-0.24, 0.24))
  expect_error(hist.zplot_brma(object, xlim = c(2, -2)), "'xlim' must be an increasing")

  testthat::skip_if_not_installed("ggplot2")
  plot <- hist.zplot_brma(object, plot_type = "ggplot", xlim = c(-2, 2))
  expect_equal(plot$coordinates$limits$x, c(-2, 2))
  plot <- plot.zplot_brma(object, plot_type = "ggplot", plot_fit = FALSE,
                          plot_thresholds = FALSE, xlim = c(-2, 2))
  expect_equal(plot$coordinates$limits$x, c(-2, 2))
})

test_that("zplot grids keep the requested range and bin count", {

  expect_length(.zplot_bins(-6, 6, .5, 10L, type = "hist", thresholds = numeric()), 11L)
  expect_length(.zplot_bins(-6, 6, .5, 10L, type = "dens", thresholds = numeric()), 10L)
  expect_identical(.zplot_bins(-6, 6, .5, NULL, thresholds = numeric()), seq(-6, 6, .5))
  for (type in c("hist", "dens")) {
    bins <- .zplot_bins(-6, 6, .7, NULL, type = type, thresholds = numeric())
    expect_identical(range(bins), c(-6, 6))
    expect_true(all(diff(bins) > 0))
  }

  object <- list(zplot = list(data = list(z = c(-1, 0, 1))), priors = NULL)
  expect_identical(nrow(hist.zplot_brma(object, length.out = 10, as_data = TRUE)), 10L)
  expect_error(hist.zplot_brma(object, from = 3, to = -3), "'from' smaller than 'to'")
  expect_error(hist.zplot_brma(object, by = 0), "'by' must be a positive")
  expect_error(hist.zplot_brma(object, length.out = 0), "length.out")
  expect_error(plot.zplot_brma(object, plot_fit = FALSE, by.hist = -1),
               "'by.hist' must be a positive")
})
