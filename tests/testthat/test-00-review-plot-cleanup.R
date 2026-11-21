test_that("posterior aggregate bias and composite density requests refuse at admission", {

  object <- bselmodel(yi = c(-.1, 0, .1), sei = rep(.1, 3), measure = "GEN",
    prior_bias = BayesTools::prior_weightfunction("one-sided", .025,
      BayesTools::wf_fixed(c(1, .2))), prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)
  expect_error(plot.brma(object, parameter = "bias"), paste0(
    "Posterior density plots are unavailable for aggregate 'bias'. ",
    "Use 'parameter = \"PET\"', 'parameter = \"PEESE\"', or ",
    "'parameter = \"omega\"' to plot a supported component."), fixed = TRUE)
  for (method in c("qCMDE", "IWMDE")) {
    error <- tryCatch(plot.brma(object, parameter = "omega", density_method = method), error = identity)
    expect_s3_class(error, "RoBMA_density_method_bias_target")
    expect_s3_class(error, "RoBMA_density_method_unavailable")
    expect_identical(conditionMessage(error), paste0(
      "qCMDE/IWMDE density plots are unavailable for composite ",
      "weightfunction targets. Use 'density_method = \"KDE\"'."))
  }
  for (argument in c("residual", "type", "unit", "conditioning_depth")) {
    args <- stats::setNames(list(FALSE), argument)
    expect_error(do.call(bfunnel.brma, c(list(x = object), args)), argument, fixed = TRUE)
  }
  # Aggregate prior plotting remains supported by the existing prior renderer.
  expect_s3_class(plot_prior(object, parameter = "bias", plot_type = "ggplot"), "ggplot")
})

test_that("factor-cell plots honor color aliases and transformed default labels", {

  object <- review_factor_replay(c("a", "b"), "treatment", distribution = "normal", measure = "OR")
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  args <- list(x = object, parameter = "g[b]", plot_type = "ggplot", prior = TRUE,
    transform = "EXP", col = "purple", color = "red", colour = "orange",
    dots_prior = list(col = "green", color = "blue", colour = "yellow"))
  plotted <- do.call(plot.brma, args)
  expect_match(ggplot2::get_labs(plotted)$x, "(odds ratio)", fixed = TRUE)
  line_layers <- which(vapply(plotted$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1L)))
  built <- ggplot2::ggplot_build(plotted)
  colors <- unique(unlist(lapply(built$data[line_layers], `[[`, "colour")))
  expect_true(all(c("purple", "green") %in% colors))
  args$par_name <- "My cell"
  expect_identical(ggplot2::get_labs(do.call(plot.brma, args))$x, "My cell")
  args$col <- NULL
  args$color <- NULL
  args$colour <- "red"
  args$dots_prior <- list(colour = "blue")
  plotted <- do.call(plot.brma, args)
  line_layers <- which(vapply(plotted$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1L)))
  built <- ggplot2::ggplot_build(plotted)
  colors <- unique(unlist(lapply(built$data[line_layers], `[[`, "colour")))
  expect_true(all(c("red", "blue") %in% colors))
})

test_that("categorical regplot admission precedes numeric evaluation validation", {

  object <- review_factor_replay(c("a", "b"), "treatment", distribution = "normal")
  expect_error(regplot(object, mod = "g", at = "b", as_data = TRUE),
    "'at' is only available for continuous moderators.", fixed = TRUE)
})

test_that("factor posterior levels select the same curve in plots and overlays", {

  levels <- c("May", "Jun", "Jul")
  first  <- review_factor_replay(levels, "orthonormal")
  second <- review_factor_replay(levels, "orthonormal", independent = TRUE)
  class(first)  <- setdiff(class(first), "only_priors.brma")
  class(second) <- setdiff(class(second), "only_priors.brma")
  colors <- .plot_level_palette(length(levels))
  full <- plot(first, "g", plot_type = "ggplot")
  one  <- plot(first, "g", level = "Jul", plot_type = "ggplot")
  separate <- plot(first, "g", level = "Jul", individual = TRUE,
    plot_type = "ggplot")
  expect_s3_class(full, "ggplot")
  expect_s3_class(one, "ggplot")
  expect_named(separate, "Jul")

  plotted_curve <- function(plot, color) {
    rows <- ggplot2::ggplot_build(plot)$data
    curves <- lapply(rows, function(layer) {
      layer[layer$colour == color, c("x", "y"), drop = FALSE]
    })
    result <- do.call(rbind, curves)
    rownames(result) <- NULL
    result
  }
  july <- plotted_curve(full, colors[3L])
  expect_gt(nrow(july), 1L)
  expect_identical(plotted_curve(one, colors[3L]), july)
  expect_identical(plotted_curve(separate[["Jul"]], colors[3L]), july)
  expect_equal(nrow(plotted_curve(one, colors[1L])), 0L)

  expect_error(lines(second, "g", level = c("Jun", "Jul"), individual = TRUE, plot_type = "ggplot"),
    "'lines()' adds to one plot. Use 'level' to select the overlay or 'plot(..., individual = TRUE)' for separate figures.",
    fixed = TRUE)
  overlay <- one + lines(second, "g", level = "Jul", individual = FALSE,
    plot_type = "ggplot", lty = 2)
  line_types <- unique(unlist(lapply(ggplot2::ggplot_build(overlay)$data, `[[`, "linetype")))
  expect_length(line_types, 2L)
  second_plot <- plot(second, "g", level = "Jul", plot_type = "ggplot")
  expected <- rbind(july, plotted_curve(second_plot, colors[3L]))
  expect_equal(plotted_curve(overlay, colors[3L]), expected)
  expect_error(plot(first, "g[Jul]", level = "Jul"),
    "Select the whole factor term", fixed = TRUE)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(first, "g", level = "Jul", prior = TRUE)
  axes <- graphics::par("usr")
  lines(second, "g", level = "Jul", lty = 2)
  expect_identical(graphics::par("usr"), axes)
})

test_that("precomputed posterior plots estimate only selected factor levels", {

  fit <- review_factor_replay(c("May", "Jun", "Jul"), "orthonormal")
  class(fit) <- setdiff(class(fit), "only_priors.brma")
  estimate <- .mock_iwmde_estimate_success()
  visited <- character()
  testthat::local_mocked_bindings(
    .iwmde_estimate = function(..., metadata) {
      cell <- metadata[["level"]]
      visited <<- c(visited, cell)
      if (identical(cell, "May")) {
        return(list(diagnostics = list(density = list(status = "unavailable",
          reason = "May has no eligible rows."))))
      }
      estimate(..., metadata = metadata)
    }, .package = "RoBMA")

  for (method in c("qCMDE", "IWMDE")) {
    visited <- character()
    selected <- plot(fit, "g", level = "Jul", density_method = method, plot_type = "ggplot")
    expect_s3_class(selected, "ggplot")
    expect_identical(visited, "Jul")
    expect_error(plot(fit, "g", level = c("Jul", "May"), density_method = method),
      class = "RoBMA_density_plot_error")
  }
})

test_that("a selected level with unit atomic mass needs no density estimate", {

  dat <- data.frame(yi = (1:6) / 10, sei = .1, g = factor(rep(c("a", "b", "c"), 2)))
  setup <- review_compile_object(brma(
    yi = yi, sei = sei, mods = ~ g, data = dat, measure = "GEN", prior_unit_information_sd = 1,
    prior_mods = list(g = BayesTools::prior_factor("point", list(.2), contrast = "treatment")),
    prior_heterogeneity = BayesTools::prior("point", list(.1)), only_priors = TRUE, silent = TRUE
  ))
  fit <- review_attach_replay(setup, cbind(mu_intercept = seq(-.5, .5, length.out = 30),
    "mu_g[1]" = .2, "mu_g[2]" = .2, tau = .1))
  class(fit) <- setdiff(class(fit), "only_priors.brma")
  testthat::local_mocked_bindings(.iwmde_estimate_density = function(...) {
    stop("An atomic level has no continuous component to estimate.")
  }, .package = "RoBMA")
  reference <- ggplot2::ggplot_build(plot(fit, "g", level = "c", plot_type = "ggplot"))$data
  for (method in c("qCMDE", "IWMDE")) {
    plotted <- plot(fit, "g", level = "c", density_method = method, plot_type = "ggplot")
    expect_equal(ggplot2::ggplot_build(plotted)$data, reference)
  }
})
