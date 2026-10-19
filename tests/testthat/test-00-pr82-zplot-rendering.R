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
  bins <- .zplot_bins(priors, from = 0, to = 3, by = 1.5, length.out = NULL)
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
    histogram <- .zplot_bins(priors, -6, 6, .5, NULL,
      type = "hist", thresholds = thresholds)
    density <- .zplot_bins(priors, -6, 6, .05, NULL,
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
