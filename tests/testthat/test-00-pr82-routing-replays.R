# Attach a synthetic posterior to a constructor-only object: compile its
# formulas, then finalize the fit through the BayesTools fit contract. The
# draws are the matrix 'samples', or the result of calling 'samples' on the
# default draws (every column of the compiled priors zero, point priors at their
# value); without 'samples', the default draws.
.replay_object <- function(object, samples = NULL) {

  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- lapply(names(args[["formula_list"]]), function(parameter) {
    BayesTools::JAGS_formula(
      formula = args[["formula_list"]][[parameter]], parameter = parameter,
      data = args[["formula_data_list"]][[parameter]],
      prior_list = args[["formula_prior_list"]][[parameter]],
      formula_scale = args[["formula_scale_list"]][[parameter]],
      prior_random = args[["formula_random_prior_list"]][[parameter]],
      random_effects_compile = args[["formula_random_effects_compile_list"]][[parameter]]
    )
  })
  names(compiled) <- names(args[["formula_list"]])
  design <- lapply(compiled, `[[`, "formula_design")
  priors <- list()
  for (item in compiled) priors <- c(priors, item[["prior_list"]])

  if (is.null(samples) || is.function(samples)) {
    columns <- unique(unlist(lapply(names(priors), function(parameter) {
      BayesTools:::.prior_linear_prior_columns(parameter, priors[[parameter]])
    })))
    defaults <- matrix(0, 2L, length(columns), dimnames = list(NULL, columns))
    for (name in names(priors)) {
      if (BayesTools::is.prior.point(priors[[name]])) {
        defaults[, name] <- mean(priors[[name]])
      }
    }
    samples <- if (is.function(samples)) samples(defaults) else defaults
  }

  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "list")
  )
  attr(fit, "formula_design") <- design
  attr(fit, "formula_scale")  <- lapply(compiled, `[[`, "formula_scale")
  attr(fit, "prior_list")     <- c(priors, .create_fit_priors(object[["data"]],
                                                                object[["priors"]]))
  object[["formula_design"]] <- design
  object[["fit"]]            <- as_bayestools_fit(fit)

  list(object = object, samples = samples)
}

# A selected multivariate model whose estimate-level source is integrated and
# whose study-level source is conditioned, with both standard deviations fixed
# at one: the row variances of the two sources are x^2 and z^2.
.selected_source_replay <- function(dat) {

  object <- bselmodel.mv(
    yi = yi, vi = vi,
    random = list(estimate = ~ id(0 + x | esid), context = ~ id(0 + z | study)),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    standardize_continuous_predictors = FALSE,
    prior_heterogeneity = BayesTools::prior_random(
      estimate = BayesTools::random_block(sd = BayesTools::prior("point", list(1))),
      context  = BayesTools::random_block(sd = BayesTools::prior("point", list(1)))
    ),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .025, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model()
    ),
    only_priors = TRUE, silent = TRUE
  )
  .replay_object(object)
}

test_that("a selected outcome funnel needs one common source split", {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4), vi = rep(.04, 4L), esid = factor(1:4),
    study = factor(c("a", "a", "b", "b")), x = c(1, 1, 0, 0), z = c(0, 0, 1, 1)
  )

  # The total heterogeneity is common (1 on every row) but the integrated and
  # conditioned parts swap between the first and the last two rows.
  for (rows in list(1:4, c(3L, 4L, 1L, 2L))) {
    replay <- .selected_source_replay(dat[rows, ])
    common <- .funnel_common_heterogeneity(replay[["object"]], replay[["samples"]])
    expect_false(common[["sources"]][["common"]])
    expect_false(common[["common"]])
    expect_true(common[["total_common"]])
    expect_null(common[["tau"]])
    # the refusal names the split, not study-specific heterogeneity
    expect_error(
      bfunnel(replay[["object"]], as_data = TRUE),
      paste0("outcome contours require one common marginal heterogeneity ",
             "distribution across studies. The total heterogeneity is common, ",
             "but the split of its selected sources into integrated and ",
             "conditioned parts varies across rows; use a LOO-PIT residual ",
             "funnel"),
      fixed = TRUE
    )
  }
  expect_error(
    .funnel_stop_no_common_heterogeneity("funnel()", list(total_common = FALSE)),
    "This model has study-specific marginal heterogeneity", fixed = TRUE
  )

  # A split that applies to every row keeps the outcome contour.
  dat[["x"]] <- 1
  dat[["z"]] <- 1
  replay <- .selected_source_replay(dat)
  common <- .funnel_common_heterogeneity(replay[["object"]], replay[["samples"]])
  expect_true(common[["sources"]][["common"]])
  expect_true(common[["common"]])
  expect_equal(common[["sources"]][["integrated"]], c(1, 1), tolerance = 1e-12)
  expect_equal(common[["sources"]][["conditioned"]], c(1, 1), tolerance = 1e-12)
  expect_equal(common[["tau"]], rep(sqrt(2), 2L), tolerance = 1e-12)
})


# Nested random intercepts, `study` (sampled) and `esid_study` (marginalized),
# with draws from `draws`, a function of the default draws.
.regplot_nested_replay <- function(scale, draws) {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4), study = factor(c("a", "a", "b", "b")),
    esid = factor(1:4), x = c(0, 1, 0, 1)
  )
  object <- brma.mv(
    yi = yi, V = rep(.04, 4L), random = ~ 1 | study / esid, scale = scale,
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE
  )
  replay <- .replay_object(object, draws)
  replay[["dat"]] <- dat
  replay
}

# The regplot grid, standardized by the fitted predictor as the scale
# regression sees it.
.regplot_grid <- function(dat) {

  grid <- data.frame(x = c(-1, .5, 2))
  list(
    newdata = grid,
    z       = (grid[["x"]] - mean(dat[["x"]])) / stats::sd(dat[["x"]])
  )
}

test_that("regplot counts the allocated variance of a common scale once", {

  replay <- .regplot_nested_replay(~ x, function(draws) {
    draws[, "mu__xRE_ALLOCx_heterogeneity__weight[1]"] <- .2
    draws[, "mu__xRE_ALLOCx_heterogeneity__weight[2]"] <- .8
    draws[, "log_tau_intercept"] <- c(.5, .8)
    draws[, "log_tau_x"]         <- c(.4, -.3)
    draws
  })
  object  <- replay[["object"]]
  samples <- replay[["samples"]]
  grid    <- .regplot_grid(replay[["dat"]])
  tau     <- samples[, "log_tau_intercept"] * exp(outer(samples[, "log_tau_x"], grid[["z"]]))

  # The marginalized block owns its allocation share; the sampled block the rest.
  marginalized <- .regplot_mv_marginalized_variance_samples(
    object, grid[["newdata"]], samples
  )
  expect_equal(unname(marginalized), .2 * tau^2, tolerance = 1e-12)
  total <- .regplot_mv_random_sd_samples(object, grid[["newdata"]], samples)
  expect_equal(unname(total^2), unname(tau^2), tolerance = 1e-12)

  # The pointwise fallback sums the same block-owned variances.
  testthat::local_mocked_bindings(
    .brma_mv_random_effects_marginal_vcov = function(...) {
      stop("JAGS_bridgesampling() row-indexed external SD source 'tau'.", call. = FALSE)
    },
    .package = "RoBMA"
  )
  fallback <- .regplot_mv_random_sd_samples(object, grid[["newdata"]], samples)
  expect_equal(unname(fallback^2), unname(tau^2), tolerance = 1e-12)
})

test_that("regplot evaluates a component-specific scale on the prediction grid", {

  replay <- .regplot_nested_replay(list(esid = ~ x), function(draws) {
    draws[, "mu__xREx__study_intercept"]    <- c(.3, .6)
    draws[, "log_tau_esid_study_intercept"] <- c(.5, .8)
    draws[, "log_tau_esid_study_x"]         <- c(.4, -.3)
    draws
  })
  object  <- replay[["object"]]
  samples <- replay[["samples"]]
  grid    <- .regplot_grid(replay[["dat"]])
  tau_esid <- samples[, "log_tau_esid_study_intercept"] *
    exp(outer(samples[, "log_tau_esid_study_x"], grid[["z"]]))
  study_variance <- matrix(samples[, "mu__xREx__study_intercept"]^2, 2L, 3L)

  marginalized <- .regplot_mv_marginalized_variance_samples(
    object, grid[["newdata"]], samples
  )
  expect_equal(unname(marginalized), unname(tau_esid^2), tolerance = 1e-12)
  total <- .regplot_mv_random_sd_samples(object, grid[["newdata"]], samples)
  expect_equal(unname(total^2), unname(study_variance + tau_esid^2), tolerance = 1e-12)

  testthat::local_mocked_bindings(
    .brma_mv_random_effects_marginal_vcov = function(...) {
      stop("JAGS_bridgesampling() row-indexed external SD source 'tau'.", call. = FALSE)
    },
    .package = "RoBMA"
  )
  fallback <- .regplot_mv_random_sd_samples(object, grid[["newdata"]], samples)
  expect_equal(unname(fallback^2), unname(study_variance + tau_esid^2), tolerance = 1e-12)
})


# Reuse the compiler/finalized-fit replay above for actual public zplot calls.
# These repeated point-prior states test display routing, not MCMC accuracy.
.review_display_replay <- function(selected = FALSE) {

  yi <- seq(-5.5, 5.5, length.out = 24L)
  point <- BayesTools::prior("point", list(location = 0))
  if (selected) {
    object <- bselmodel.mv(yi = yi, V = diag(24L), random = NULL,
      measure = "GEN", prior_effect = point,
      prior_bias = BayesTools::prior_weightfunction("one-sided", .025,
        BayesTools::wf_fixed(c(1, .2))), prior_unit_information_sd = 1,
      only_priors = TRUE, silent = TRUE)
    samples <- matrix(0, 24L, 1L, dimnames = list(NULL, "mu"))
  } else {
    object <- brma(yi = yi, sei = rep(1, 24L), measure = "GEN",
      prior_effect = point,
      prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
      prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
    samples <- cbind(mu = rep(0, 24L), tau = rep(.1, 24L))
  }
  .replay_object(object, samples)[["object"]]
}


test_that("actual ordinary heterogeneity summaries identify their component", {

  result <- summary_heterogeneity(.review_display_replay())
  frame <- as.data.frame(result)
  expect_identical(result[["component"]], "heterogeneity")
  expect_identical(frame[["component"]], rep("heterogeneity", 4L))
  expect_identical(frame[["parameter"]], c("tau", "tau2", "I2", "H2"))
  expect_identical(data.frame(result), frame)
  expect_equal(unname(frame[["Mean"]]), c(.1, .01, 100 * .01 / 1.01, 1.01),
    tolerance = 1e-12)
})


test_that("zplot public controls distinguish explicit, stored and legacy inputs", {

  skip_if_not_installed("ggplot2")
  invalid_message <- "'integration_control' must be created by set_selection_likelihood_control()."
  legacy_message <- paste0(
    "'integration_control' is unavailable for this legacy zplot object. ",
    "Supply 'integration_control = set_selection_likelihood_control()' ",
    "or recreate the object with 'as_zplot()'.")
  control <- set_selection_likelihood_control()
  override <- set_selection_likelihood_control(seed = 203L)
  for (selected in c(FALSE, TRUE)) {
    object <- .review_display_replay(selected)
    for (invalid in list(NULL, list())) {
      expect_error(as_zplot(object, integration_control = invalid),
        invalid_message, fixed = TRUE)
    }
    z <- as_zplot(object, integration_control = control)
    expect_identical(z[["zplot"]][["data"]][["integration_control"]], control)
    stored <- lines(z, plot_type = "ggplot", as_data = TRUE, max_samples = Inf)
    expect_equal(lines(z, plot_type = "ggplot", as_data = TRUE,
      max_samples = Inf, integration_control = override), stored, tolerance = 1e-12)
    for (invalid in list(NULL, list())) {
      expect_error(lines(z, as_data = TRUE, integration_control = invalid),
        invalid_message, fixed = TRUE)
      expect_error(plot(z, plot_type = "ggplot", integration_control = invalid),
        invalid_message, fixed = TRUE)
      expect_error(plot(z, plot_type = "ggplot",
        dots_fit = list(integration_control = invalid)), invalid_message, fixed = TRUE)
      corrupt <- z
      corrupt[["zplot"]][["data"]]["integration_control"] <- list(invalid)
      expect_error(lines(corrupt, as_data = TRUE), invalid_message, fixed = TRUE)
      expect_error(plot(corrupt, plot_type = "ggplot"), invalid_message, fixed = TRUE)
    }
    legacy <- z
    legacy[["zplot"]][["data"]][["integration_control"]] <- NULL
    expect_s3_class(plot(legacy, plot_type = "ggplot", plot_fit = FALSE,
      plot_extrapolation = FALSE), "ggplot")
    expect_s3_class(plot(legacy, plot_type = "ggplot", plot_extrapolation = TRUE,
      dots_fit = list(integration_control = control),
      dots_extrapolation = list(integration_control = control)), "ggplot")
    for (both in c(FALSE, TRUE)) {
      expect_s3_class(plot(legacy, plot_type = "ggplot",
        plot_extrapolation = both, integration_control = control), "ggplot")
      if (selected) {
        expect_error(plot(legacy, plot_type = "ggplot", plot_extrapolation = both),
          legacy_message, fixed = TRUE)
      } else {
        expect_s3_class(plot(legacy, plot_type = "ggplot",
          plot_extrapolation = both), "ggplot")
      }
    }
    if (selected) {
      error <- tryCatch(lines(legacy, as_data = TRUE), error = identity)
      expect_identical(conditionMessage(error), legacy_message)
      expect_null(conditionCall(error))
    } else {
      expect_equal(lines(legacy, as_data = TRUE, max_samples = Inf), stored,
        tolerance = 0)
    }
    expect_equal(lines(legacy, as_data = TRUE, max_samples = Inf,
      integration_control = control), stored, tolerance = 0)
  }
})


test_that("standalone ggplot histograms scale to actual added density layers", {

  skip_if_not_installed("ggplot2")
  z <- as_zplot(.review_display_replay())
  bars <- hist(z, plot_type = "ggplot", as_data = TRUE, plot_thresholds = FALSE)
  density <- lines(z, plot_type = "ggplot", as_data = TRUE, plot_ci = FALSE)
  expect_gt(max(density[["y"]]), max(bars[["density"]]))
  histogram <- hist(z, plot_type = "ggplot", plot_thresholds = FALSE)
  composed <- histogram + lines(z, plot_type = "ggplot", plot_ci = FALSE)
  expect_null(histogram[["coordinates"]][["limits"]][["y"]])
  built <- ggplot2::ggplot_build(composed)
  expect_gte(built[["layout"]][["panel_params"]][[1L]][["y.range"]][[2L]],
    max(density[["y"]]))
  explicit <- hist(z, plot_type = "ggplot", ylim = c(0, .1),
    plot_thresholds = FALSE) + lines(z, plot_type = "ggplot", plot_ci = FALSE)
  expect_identical(explicit[["coordinates"]][["limits"]][["y"]], c(0, .1))
  combined <- plot(z, plot_type = "ggplot", plot_ci = FALSE,
    plot_thresholds = FALSE)
  expect_equal(combined[["coordinates"]][["limits"]][["y"]],
    c(0, max(density[["y"]]) * 1.05), tolerance = 0)
  expect_identical(.get_dots_hist_zplot(list(), "base", max(bars[["density"]]))[["ylim"]],
    c(0, max(bars[["density"]])))
  expect_error(hist(z, plot_type = "ggplot", add = TRUE), paste0(
    "'add = TRUE' is unavailable for 'plot_type = \"ggplot\"'. ",
    "Combine 'hist()' and 'lines()' layers with '+'."), fixed = TRUE)
})


test_that("zplot component styles override global aliases in actual layers", {

  skip_if_not_installed("ggplot2")
  z <- as_zplot(.review_display_replay(TRUE))
  layer <- function(plot, geom, occurrence = 1L) {
    indices <- which(vapply(plot[["layers"]], function(layer) {
      inherits(layer[["geom"]], geom)
    }, logical(1L)))
    ggplot2::layer_data(plot, indices[[occurrence]])
  }
  styles <- list(
    list(col = "red", lty = 2, lwd = 2),
    list(color = "red", linetype = 2, linewidth = 2),
    list(color = "red", col = "orange", linetype = 2, lty = 5,
      linewidth = 2, lwd = 5)
  )
  for (style in styles) {
    plot <- plot(z, plot_type = "ggplot", plot_ci = FALSE,
      plot_extrapolation = TRUE, color = "black", linetype = 1, linewidth = .5,
      dots_fit = style,
      dots_extrapolation = list(col = "blue", lty = 3, lwd = 3),
      dots_thresholds = list(col = "green", lty = 4, lwd = 4))
    fit <- layer(plot, "GeomLine", 1L)
    extrapolation <- layer(plot, "GeomLine", 2L)
    threshold <- layer(plot, "GeomVline")
    expect_identical(unique(fit[["colour"]]), "red")
    expect_equal(unique(fit[["linetype"]]), 2)
    expect_equal(unique(fit[["linewidth"]]), 2)
    expect_identical(unique(extrapolation[["colour"]]), "blue")
    expect_equal(unique(extrapolation[["linetype"]]), 3)
    expect_equal(unique(extrapolation[["linewidth"]]), 3)
    expect_identical(unique(threshold[["colour"]]), "green")
    expect_equal(unique(threshold[["linetype"]]), 4)
    expect_equal(unique(threshold[["linewidth"]]), 4)
  }
  global <- hist(z, plot_type = "ggplot", color = "purple", linetype = 5,
    linewidth = 2)
  threshold <- layer(global, "GeomVline")
  expect_identical(unique(threshold[["colour"]]), "purple")
  expect_equal(unique(threshold[["linetype"]]), 5)
  expect_equal(unique(threshold[["linewidth"]]), 2)
})


test_that("zplot reconversion and summary wording retain the public result", {

  z <- as_zplot(.review_display_replay())
  again <- as_zplot(z)
  expect_identical(sum(class(again) == "zplot_brma"), 1L)
  printed <- capture.output(value <- withVisible(print(z)))
  expect_s3_class(value[["value"]], "summary.zplot_brma")
  expect_false(value[["visible"]])
  expect_true(any(grepl("conditioning depth: marginal.", printed, fixed = TRUE)))
  expect_false(any(grepl("conditioning depth: marginal).", printed, fixed = TRUE)))
})
