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
