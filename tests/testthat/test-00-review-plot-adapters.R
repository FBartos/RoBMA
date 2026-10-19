.review_plot_adapter_replay <- function(selected = FALSE) {

  args <- list(yi = c(.1, .2, .3, .4), sei = rep(.2, 4L), measure = "GEN",
    prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    only_priors = TRUE, silent = TRUE)
  if (selected) args$prior_bias <- BayesTools::prior_weightfunction("one-sided",
    .025, BayesTools::wf_fixed(c(1, .2)))
  object <- do.call(if (selected) bselmodel else brma, args)
  replay <- .plot_range_point_replay(object)
  stopifnot(identical(class(replay$object), class(object)))
  replay$object
}

test_that("admitted forest adapters preserve presentation and explicit overrides", {

  fit <- .review_plot_adapter_replay()
  slab <- paste0("Custom study ", seq_len(nobs(fit)))
  prepared <- forest.brma(fit, addpred = TRUE, predstyle = "shade", row = -4,
    mlab = "Prepared", predlim = c(-2, 2), border = "red", constarea = TRUE,
    slab = slab, as_data = TRUE)
  expect_identical(prepared$studies$slab, slab)
  expect_identical(prepared$forest_args$slab, slab)
  expect_identical(attr(prepared$forest_args$x, "slab"), slab)
  expect_identical(prepared$addpoly_args$rows, -4)
  expect_identical(prepared$addpoly_args$predstyle, "shade")
  expect_identical(prepared$addpoly_args$mlab, "Prepared")
  expect_identical(prepared$addpoly_args$predlim, c(-2, 2))
  expect_identical(prepared$addpoly_args$border, "red")
  expect_identical(prepared$addpoly_args$constarea, TRUE)
  bare <- forest.brma(fit, addfit = FALSE, as_data = TRUE)
  expect_identical(bare$plot_defaults$addfit, FALSE)
  expect_null(bare$forest_args$ylim)
  custom <- forest.brma(fit, addfit = FALSE, ylim = c(-5, 30), as_data = TRUE)
  expect_identical(custom$forest_args$ylim, c(-5, 30))
  .with_temp_plot_device({
    direct <- forest.brma(fit, addpred = TRUE, predstyle = "shade", row = -4,
      mlab = "Prepared", predlim = c(-2, 2), border = "red", constarea = TRUE, slab = slab)
    replay <- forest.metafor_forest.brma(prepared)
    expect_equal(replay, direct)
    expect_silent(forest.metafor_forest.brma(prepared, predstyle = "bar", row = -2,
      mlab = "Override", border = "blue"))
    expect_equal(forest.metafor_forest.brma(bare), forest.brma(fit, addfit = FALSE))
    legacy <- as_metafor_forest(fit)
    expect_null(legacy$plot_defaults)
    expect_equal(forest.metafor_forest.brma(legacy), forest.brma(fit))
  })
})

test_that("admitted selected forest and funnels retain reference provenance", {

  fit <- .review_plot_adapter_replay(selected = TRUE)
  posterior <- .get_posterior_samples(fit$fit)
  withr::local_seed(481)
  before <- .Random.seed
  pooled <- pooled_effect(fit, bias_adjusted = FALSE, .posterior_samples = posterior)
  expect_identical(.Random.seed, before)
  explicit <- predict(fit, newdata = data.frame(sei = .1), type = "estimate",
    bias_adjusted = FALSE, quiet = TRUE, .posterior_samples = posterior)
  expect_identical(attr(explicit, "RoBMA_target", exact = TRUE)$event, "new_publication")
  expect_false(identical(.Random.seed, before))
  prepared <- forest.brma(fit, addpred = TRUE, bias_adjusted = FALSE, as_data = TRUE)
  expect_identical(attr(prepared, "RoBMA_target", exact = TRUE)$reported_target,
    "fixed-location average-expanded-design reference")
  specified <- forest.brma(fit, addpred = TRUE, bias_adjusted = FALSE,
    newdata = data.frame(sei = .1), as_data = TRUE)
  expect_identical(attr(specified, "RoBMA_target", exact = TRUE)$event, "new_publication")
  for (fun in list(funnel, bfunnel)) {
    actual <- fun(fit, refline = .1, ylim = c(0, .27), xlim = c(-1, 1),
      max_samples = 24L, as_data = TRUE)
    expect_identical(actual$refline, data.frame(x = c(.1, .1), y = c(0, .27)))
    target <- attr(actual, "RoBMA_target", exact = TRUE)
    expect_identical(target$reference_law, "independently_published_new_estimate_product_rule")
    expect_identical(target$selection_model, .data_selection_model(fit$data))
    expect_true(actual$estimand %in% c("plugin", "posterior_predictive"))
    plotted <- fun(fit, refline = .1, ylim = c(0, .27), xlim = c(-1, 1),
      max_samples = 24L, plot_type = "ggplot")
    reference_layer <- which(vapply(plotted$layers, function(layer) {
      identical(layer$data, actual$refline)
    }, logical(1L)))
    expect_length(reference_layer, 1L)
    if (length(reference_layer) == 1L) {
      reference <- ggplot2::ggplot_build(plotted)$data[[reference_layer]]
      expect_identical(sum(is.finite(reference$y)), 2L)
      expect_equal(reference$y, -actual$refline$y)
    }
    .with_temp_plot_device(expect_null(fun(fit, refline = .1,
      ylim = c(0, .27), xlim = c(-1, 1), max_samples = 24L, plot_type = "base")))
  }
})

test_that("admitted known-V QQ data retain actual LOO deletion context", {

  fit <- add_loo(.review_cleanup_plain_replay(known_v = TRUE))
  qq <- qqnorm.brma(fit, as_data = TRUE)
  expected <- .selection_deletion_retained_context(fit$data, fit$priors)
  expect_identical(expected, "remaining_data")
  expect_identical(qq$unit, "estimate")
  expect_null(qq$conditioning_depth)
  expect_identical(qq$retained_context, expected)
  expect_false("conditioning_depth" %in% names(qq$points))
  expect_identical(qq$points$retained_context, rep(expected, nobs(fit)))
  standard <- qqnorm.brma(fit, type = "rstandard", as_data = TRUE)
  expect_identical(standard$conditioning_depth, "marginal")
  expect_null(standard$retained_context)
})
