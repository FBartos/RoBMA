context("(02) qCMDE bound rows at a finite support boundary")

source(testthat::test_path("common-functions.R"))


test_that("40-point PET density curves agree with finer normalization", {

  skip_on_cran()
  # PET rows have a Gaussian kernel times a truncated Cauchy prior on (0, Inf).
  # Compare all displayed ordinates with the finer 400-point normalization.
  # The PET-only fit needs a finer normalization grid to resolve the
  # conditional ordinates of the current posterior draws.
  budgets <- c("dat.lehmann2018-PET" = 100L, "dat.lehmann2018_RoBMA" = 50L)
  skip_if_missing_fits(names(budgets))
  for (name in names(budgets)) {
    fit  <- load_fit(name)
    plot <- plot(fit, parameter = "PET", density_method = "qCMDE",
      density_control = list(n_points = 40L, normalization_points = budgets[[name]]),
      plot_type = "ggplot")
    expect_s3_class(plot, "ggplot")
    line <- which(vapply(plot[["layers"]], function(layer) {
      inherits(layer[["geom"]], "GeomLine")
    }, logical(1)))
    expect_length(line, 1L)
    curve <- plot[["layers"]][[line]][["data"]]
    expect_true(all(is.finite(curve[["y"]]) & curve[["y"]] >= 0))
    expect_true(all(curve[["x"]] > 0))
    reference <- plot(fit, parameter = "PET", density_method = "qCMDE",
      density_control = list(n_points = 40L, normalization_points = 400L), plot_type = "ggplot")
    reference <- reference[["layers"]][[line]][["data"]]
    expect_equal(curve[["x"]], reference[["x"]], tolerance = 1e-12)
    expect_lt(max(abs(curve[["y"]] / reference[["y"]] - 1)), 0.001)
  }
})
