context("(02) qCMDE bound rows at a finite support boundary")

source(testthat::test_path("common-functions.R"))


test_that("PET qCMDE density lines on smaller grids are not rejected", {

  skip_on_cran()
  # PET rows are a Gaussian kernel times a truncated Cauchy prior on (0, Inf).
  # Their symmetric Gaussian-bound interval crosses 0; clamping it there put
  # the lower end at about 1e-8, so a log-chart range too wide for 50
  # normalization points and a validation change above 5% that rejected these
  # lines. Each end now sits where its own side's bound meets the target.
  fits <- c("dat.lehmann2018-PET", "dat.lehmann2018_RoBMA")
  skip_if_missing_fits(fits)
  for (name in fits) {
    fit  <- load_fit(name)
    plot <- plot(fit, parameter = "PET", density_method = "qCMDE",
      density_control = list(n_points = 40L), plot_type = "ggplot")
    expect_s3_class(plot, "ggplot")
    line <- which(vapply(plot[["layers"]], function(layer) {
      inherits(layer[["geom"]], "GeomLine")
    }, logical(1)))
    expect_length(line, 1L)
    curve <- plot[["layers"]][[line]][["data"]]
    expect_true(all(is.finite(curve[["y"]]) & curve[["y"]] >= 0))
    expect_true(all(curve[["x"]] > 0))
  }
})
