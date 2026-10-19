source(testthat::test_path("common-functions.R"))

test_that("spike-and-slab plotting preserves the parent probability measure", {

  object <- brma(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  # Isolate frontend dispatch here; compiler-backed random selectors have
  # their own fitted checks below and in the plotting verification gallery.
  for (probability in c(0, .1, 1)) {
    parent <- BayesTools::prior_spike_and_slab(
      prior("normal", list(mean = 0, sd = 1), list(lower = 0, upper = Inf)),
      prior("point", list(location = probability)))
    testthat::local_mocked_bindings(
      .select_plot_prior_parameter = function(...) list(
        prior = parent, label = "tau", source = "random"),
      .package = "RoBMA")
    plot <- plot_prior(object, "tau", plot_type = "ggplot")
    lines <- Filter(function(x) inherits(x$geom, "GeomLine"), plot$layers)
    atoms <- Filter(function(x) inherits(x$geom, "GeomSegment"), plot$layers)
    if (probability > 0) {
      expect_length(lines, 1L)
      curve <- lines[[1L]]$data
      keep <- curve$x > 0
      expect_true(any(keep))
      expect_true(all(is.finite(curve$y[keep])))
      expect_equal(curve$y[keep] / (2 * stats::dnorm(curve$x[keep])),
        rep(probability, sum(keep)), tolerance = 1e-12)
    }
    if (probability < 1) {
      expect_length(atoms, 1L)
      if (probability == 0) {
        # A pure point uses the primary probability axis. Read the producer's
        # probability data, rather than treating an arrow height as mass.
        producer <- get(".plot_data_prior_list.simple", asNamespace("BayesTools"))
        data <- producer(list(parent), x_seq = c(.1, .5, 1), x_range = NULL,
          x_range_quant = NULL, n_points = 500, n_samples = 10000,
          force_samples = FALSE, individual = FALSE, transformation = NULL,
          transformation_arguments = NULL, transformation_settings = FALSE)
        points <- Filter(function(x) inherits(x, "density.prior.point"), data)
        expect_length(points, 1L)
        expect_identical(points[[1L]]$y, 1)
        expect_identical(plot$scales$get_scales("y")$name, "Probability")
      } else {
        expect_true(is.finite(plot$bt_scale_y2_state$scale_y2))
        expect_equal(atoms[[1L]]$data$yend / plot$bt_scale_y2_state$scale_y2,
          1 - probability, tolerance = 1e-12)
      }
    }
    if (probability == 0) expect_length(lines, 0L)
    if (probability == 1) expect_length(atoms, 0L)
    .with_temp_plot_device(expect_silent(plot_prior(object, "tau")))
  }
})

# plot_prior() and print_prior() of the random-effect selectors of fitted
# objects. A quantity that is exactly one declared prior under the identity
# transform shows that prior; any other quantity shows the prior law
# BayesTools derives from the declared priors (parameter_prior_density(), the
# law plot(prior = TRUE) overlays) without a posterior. The expectations are
# analytic identities of the fits' own declared priors.

test_that("plot_prior and print_prior show the declared prior of a random SD", {

  skip_on_cran()
  fit   <- load_fit("brma.mv_block_mvn_random")
  prior <- attr(fit[["fit"]], "prior_list")[["mu__xREx__estimate_intercept"]]
  expect_true(BayesTools::is.prior(prior))

  plot  <- plot_prior(fit, "tau", plot_type = "ggplot")
  curve <- .random_prior_curve(plot)
  expect_gt(nrow(curve), 500L)
  .expect_random_prior_curve(
    curve, exp(BayesTools::lpdf(prior, curve[["x"]])), 1e-8,
    support = unlist(prior[["truncation"]])
  )
  expect_equal(
    .random_prior_curve(plot_prior(
      fit, "tau", plot_type = "ggplot", standardized_coefficients = FALSE
    )),
    curve
  )
  expect_equal(
    .random_prior_curve(plot_prior(
      fit, "tau", plot_type = "ggplot", component = "random"
    )),
    curve
  )

  # base plots run and return nothing
  .with_temp_plot_device(expect_null(plot_prior(fit, "tau")))

  printed <- capture.output(returned <- print_prior(fit, "tau"))
  expect_identical(returned, prior)
  expect_identical(printed[[1L]], "tau:")
  expect_true(all(
    paste0("  ", print(prior, silent = TRUE)) %in% printed
  ))
  expect_silent(print_prior(fit, "tau", silent = TRUE))

  # several selectors, and the whole specification of print_prior(fit, "random")
  plots <- plot_prior(fit, c("tau", "tau2"), plot_type = "ggplot")
  expect_named(plots, c("tau", "tau2"))
  expect_s3_class(
    print_prior(fit, "random", silent = TRUE), "prior_random"
  )
})

test_that("plot_prior draws the prior of a variance as the square of the SD prior", {

  skip_on_cran()
  fit   <- load_fit("brma.mv_block_mvn_random")
  prior <- attr(fit[["fit"]], "prior_list")[["mu__xREx__estimate_intercept"]]

  # the density of v = sigma^2 is f(sqrt(v)) / (2 sqrt(v))
  plot  <- plot_prior(fit, "tau2", plot_type = "ggplot")
  curve <- .random_prior_curve(plot)
  expect_gt(nrow(curve), 500L)
  .expect_random_prior_curve(
    curve,
    exp(BayesTools::lpdf(prior, sqrt(curve[["x"]]))) / (2 * sqrt(curve[["x"]])),
    1e-6
  )
  .with_temp_plot_device(expect_null(plot_prior(fit, "tau2")))

  # the derivation line and the declared source prior
  printed <- capture.output(returned <- print_prior(fit, "tau2"))
  expect_identical(printed[[1L]], "tau2 = tau^2")
  expect_identical(printed[[2L]], "  tau:")
  expect_true(all(
    paste0("    ", print(prior, silent = TRUE)) %in% printed
  ))
  expect_s3_class(returned, "prior_linear_density")
  expect_silent(print_prior(fit, "tau2", silent = TRUE))

  # xlim, n_points and par_name select the plotted range, grid and label
  limited <- plot_prior(
    fit, "tau2", plot_type = "ggplot", xlim = c(0.05, 0.4), n_points = 40,
    par_name = "variance"
  )
  curve <- .random_prior_curve(limited)
  expect_gte(min(curve[["x"]]), 0.05)
  expect_lte(max(curve[["x"]]), 0.4)
  expect_lte(nrow(curve), 45L)
  .expect_random_prior_curve(
    curve,
    exp(BayesTools::lpdf(prior, sqrt(curve[["x"]]))) / (2 * sqrt(curve[["x"]])),
    1e-6, minimum = 30L
  )
  expect_match(
    paste(format(ggplot2::get_labs(limited)[["x"]]), collapse = ""),
    "variance", fixed = TRUE
  )
  .with_temp_plot_device(expect_null(
    plot_prior(fit, "tau2", xlim = c(0.05, 0.4), n_points = 40, col = "red")
  ))

  # the evaluation grid of a prior law is not selected by x_seq
  expect_error(plot_prior(fit, "tau2", x_seq = c(0.1, 0.2)), "x_seq")
})

test_that("random selectors of gated allocations show prior laws with their point masses", {

  skip_on_cran()
  fit <- load_fit("BMA.mv_random_components")
  prior_list <- attr(fit[["fit"]], "prior_list")
  total   <- prior_list[["mu__xRE_ALLOCx_heterogeneity__allocation_sd"]]
  weights <- prior_list[["mu__xRE_ALLOCx_heterogeneity__weight"]]
  selection <- function(selector) {
    .brma_parameter_select_entry(fit, selector)[["selection"]]
  }
  # two components with independent inclusion gates of probability 0.5

  # the total SD is 0 when no component is included (probability 0.25) and
  # the scale prior otherwise: the mixed law has the atom and is drawn with
  # an arrow at 0
  law <- BayesTools::parameter_prior_density(fit[["fit"]], selection("tau_total"))
  expect_equal(law[["points"]][["x"]], 0)
  expect_equal(law[["points"]][["p"]], 0.25, tolerance = 1e-12)
  layers <- .random_prior_layers(plot_prior(fit, "tau_total", plot_type = "ggplot"))
  expect_length(layers[["arrows"]], 1L)
  expect_equal(layers[["arrows"]][[1L]][["x"]], 0)
  expect_gt(layers[["arrows"]][[1L]][["yend"]], 0)
  printed <- capture.output(print_prior(fit, "tau_total"))
  expect_match(
    printed[[1L]], "tau_total = sigma_total * sqrt(sum of I_c * w[c] over the components c)",
    fixed = TRUE
  )
  expect_true(any(grepl("^  w:$", printed)))
  expect_true(any(grepl("^  p_study:$", printed)))
  expect_true(any(grepl("^  p_observation:$", printed)))
  expect_true(any(grepl("Spike(0.5)", printed, fixed = TRUE)))
  expect_true(all(paste0("    ", print(weights, silent = TRUE)) %in% printed))

  # the share of a component is a mixture of the atoms at 0 and 1 and the
  # Beta(1, 1) marginal, given an active component: each of the three
  # configurations (only this component, only the other one, both) has
  # probability 1/3
  law <- BayesTools::parameter_prior_density(fit[["fit"]], selection("tau2_prop(study)"))
  expect_equal(law[["points"]][["x"]], c(0, 1))
  expect_equal(law[["points"]][["p"]], rep(1 / 3, 2), tolerance = 1e-12)
  expect_equal(law[["density"]][["mass"]], 1 / 3, tolerance = 1e-12)
  printed <- capture.output(print_prior(fit, "tau2_prop(study)"))
  expect_identical(
    printed[[1L]],
    paste0(
      "tau2_prop(study) = I_study * w[1] / sum of I_c * w[c] over the ",
      "components c, the share of the total variance of allocation ",
      "'heterogeneity' of the active components"
    )
  )
  expect_identical(
    grep(":$", printed, value = TRUE), c("  w:", "  p_study:", "  p_observation:")
  )
  layers <- .random_prior_layers(plot_prior(fit, "tau2_prop(study)", plot_type = "ggplot"))
  expect_length(layers[["arrows"]], 2L)
  expect_equal(
    sort(vapply(layers[["arrows"]], function(arrow) arrow[["x"]], numeric(1))),
    c(0, 1)
  )
  curve <- .random_prior_curve(plot_prior(fit, "tau2_prop(study)", plot_type = "ggplot"))
  .expect_random_prior_curve(
    curve, rep(1 / 3, nrow(curve)), 1e-6, support = c(0, 1)
  )

  # a component SD is drawn the same way, with the gate priors in its derivation
  layers <- .random_prior_layers(plot_prior(fit, "study: tau", plot_type = "ggplot"))
  expect_length(layers[["arrows"]], 1L)
  expect_equal(layers[["arrows"]][[1L]][["x"]], 0)
  printed <- capture.output(print_prior(fit, "study: tau"))
  expect_match(
    printed[[1L]], "study: tau = sigma_total * I_study * sqrt(w[1])", fixed = TRUE
  )
  expect_true(any(grepl("^  p_study:$", printed)))
  .with_temp_plot_device(expect_null(plot_prior(fit, "tau2_prop(study)")))
  .with_temp_plot_device(expect_null(plot_prior(fit, "tau_total", xlim = c(0, 2))))
})

test_that("random selectors of a nested gated allocation derive from its allocation chain", {

  skip_on_cran()
  fit <- load_fit("lane_nested_total")
  prior_list <- attr(fit[["fit"]], "prior_list")
  total <- prior_list[["mu__xRE_ALLOCx_heterogeneity__allocation_sd"]]
  included <- 0.5

  # the total is the scale prior times the gate: the atom at 0 and the scaled
  # scale prior, whose derivation does not involve the split weights
  law <- BayesTools::parameter_prior_density(
    fit[["fit"]], .brma_parameter_select_entry(fit, "tau_total")[["selection"]]
  )
  expect_equal(law[["points"]][["x"]], 0)
  expect_equal(law[["points"]][["p"]], 1 - included, tolerance = 1e-12)
  plot <- plot_prior(fit, "tau_total", plot_type = "ggplot")
  curve <- .random_prior_curve(plot)
  .expect_random_prior_curve(
    curve, included * exp(BayesTools::lpdf(total, curve[["x"]])),
    1e-6, support = unlist(total[["truncation"]])
  )
  expect_length(.random_prior_layers(plot)[["arrows"]], 1L)
  printed <- capture.output(print_prior(fit, "tau_total"))
  expect_identical(printed[[1L]], "tau_total = sigma_total * I_component_1")
  expect_identical(grep(":$", printed, value = TRUE), c("  sigma_total:", "  p_component_1:"))

  # a component SD follows the gate of the root and the split weights
  printed <- capture.output(print_prior(fit, "esid_study: tau"))
  expect_identical(
    printed[[1L]], "esid_study: tau = sigma_total * I_component_1 * sqrt(w[1])"
  )
  expect_identical(
    grep(":$", printed, value = TRUE), c("  sigma_total:", "  p_component_1:", "  w:")
  )
  .with_temp_plot_device(expect_null(plot_prior(fit, "esid_study: tau")))

  # a share of the split depends on the split weights alone
  printed <- capture.output(print_prior(fit, "tau2_prop(study)"))
  expect_identical(
    printed[[1L]],
    "tau2_prop(study) = w[2], the share of the total variance of allocation 'component_1_split'"
  )
  expect_identical(grep(":$", printed, value = TRUE), "  w:")
})

test_that("random selectors without a prior density are refused with the target classes", {

  skip_on_cran()
  fit <- load_fit("lane_random_correlation_direct")

  expect_error(
    plot_prior(fit, "tau(x)"),
    "prior of the random-effect quantity 'tau(x)' is unavailable",
    fixed = TRUE, class = "RoBMA_hypothesis_target"
  )
  for (call in list(
    quote(plot_prior(fit, "tau(x)")),
    quote(print_prior(fit, "tau(x)")),
    quote(plot_prior(fit, "rho(intercept,x)", standardized_coefficients = FALSE))
  )) {
    error <- tryCatch(eval(call), error = identity)
    expect_s3_class(error, "RoBMA_hypothesis_target")
    expect_s3_class(error, "RoBMA_hypothesis_unavailable")
  }

  # the correlation has its exact prior on the fitted scale
  expect_error(
    plot_prior(fit, "rho(intercept,x)", standardized_coefficients = FALSE),
    "standardized_coefficients = TRUE", fixed = TRUE
  )
  .with_temp_plot_device(expect_null(plot_prior(fit, "rho(intercept,x)")))
  expect_gt(
    nrow(.random_prior_curve(plot_prior(fit, "rho(intercept,x)", plot_type = "ggplot"))),
    100L
  )
})

test_that("only_priors objects point to the random-effect specification", {

  skip_on_cran()
  data <- data.frame(
    yi    = c(0.10, 0.25, 0.15, 0.30, 0.05),
    sei   = c(0.2, 0.25, 0.22, 0.28, 0.17),
    study = c("s1", "s1", "s2", "s2", "s3")
  )
  random <- brma.mv(
    yi = yi, V = diag(data[["sei"]]^2), random = ~ 1 | study, data = data,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE
  )

  for (call in list(
    quote(plot_prior(random, "tau")),
    quote(print_prior(random, "tau", silent = TRUE))
  )) {
    error <- tryCatch(eval(call), error = identity)
    expect_identical(class(error), c("simpleError", "error", "condition"))
    expect_match(
      conditionMessage(error), "The specified parameter 'tau' is not available."
    )
    expect_match(
      conditionMessage(error),
      "print_prior(object, parameter = \"random\")", fixed = TRUE
    )
  }
  expect_s3_class(
    print_prior(random, parameter = "random", silent = TRUE), "prior_random"
  )

  # objects without random-effect priors give the original message only
  plain <- BMA(
    yi = yi, sei = sei, data = data, measure = "SMD", only_priors = TRUE
  )
  error <- tryCatch(print_prior(plain, "tau2", silent = TRUE), error = identity)
  expect_match(conditionMessage(error), "is not available.")
  expect_false(grepl("random", conditionMessage(error), fixed = TRUE))
})

test_that("allocation derivations retain cycle, unknown-parent, and supported-depth causes", {

  allocation <- function(name, parent = NULL) list(label = name,
    parent = if (is.null(parent)) NULL else list(allocation = parent))
  spec <- list(formula_parameter = "mu")
  definitions <- lapply(seq_len(17L), function(i) allocation(as.character(i),
    if (i == 1L) NULL else as.character(i - 1L)))
  design <- list(mu = list(parameter = "mu", random_allocations = definitions))
  # Resolving the root on the last supported iteration is valid.
  expect_length(.brma_random_allocation_chain(design, spec, definitions[[16L]]), 16L)
  error <- tryCatch(.brma_random_allocation_chain(design, spec, definitions[[17L]]), error = identity)
  expect_s3_class(error, "RoBMA_random_allocation_depth_unavailable")
  expect_s3_class(error, "RoBMA_random_prior_derivation_error")
  expect_identical(error$reason, "unsupported_depth")
  expect_identical(error$supported_depth, 16L)
  expect_identical(conditionMessage(error), paste0(
    "Random-effect allocation prior derivation is unavailable: ",
    "the parent chain exceeds the supported depth of 16 allocations."))
  cyclical <- allocation("cycle", "cycle")
  design$mu$random_allocations <- list(cyclical)
  error <- tryCatch(.brma_random_allocation_chain(design, spec, cyclical), error = identity)
  expect_s3_class(error, "RoBMA_random_allocation_cycle")
  expect_s3_class(error, "RoBMA_random_prior_derivation_error")
  expect_identical(error$reason, "parent_cycle")
  expect_identical(conditionMessage(error), "Random-effect allocation metadata contain a parent cycle.")
  error <- tryCatch(.brma_random_allocation_chain(design, spec,
    allocation("leaf", "absent")), error = identity)
  expect_s3_class(error, "RoBMA_random_allocation_parent_unavailable")
  expect_identical(error$reason, "unknown_parent")
  expect_identical(conditionMessage(error), "Random-effect allocation metadata reference unknown parent 'absent'.")
})

test_that("random variance derivations preserve tau2 inside grouping names", {

  dat <- data.frame(yi = c(.1, .2, .3, .4), study = c("a", "a", "b", "b"), esid = 1:4)
  setup <- review_compile_object(brma.mv(yi = yi, V = diag(.04, 4L),
    random = list(tau2_study = ~ id(1 | study), other = ~ id(1 | esid)),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    prior_heterogeneity = BayesTools::prior_random(sd = BayesTools::prior("normal",
      list(mean = 0, sd = 1), list(lower = 0, upper = Inf))),
    only_priors = TRUE, silent = TRUE))
  priors <- c(setup$compiled$mu$prior_list,
    .create_fit_priors(setup$object$data, setup$object$priors))
  # Include the compiler-declared random monitors, rather than attaching
  # prior columns alone to a random-effect model.
  terms <- setup$compiled$mu$formula_design$random_effects
  latent_bases <- vapply(terms, function(term) paste0(term$parameter_stem, "_xRE_Zx"), character(1L))
  latent <- unlist(lapply(terms, function(term) {
    index <- expand.grid(group = seq_len(term$n_groups), column = seq_len(term$n_columns))
    paste0(term$parameter_stem, "_xRE_Zx[", index$group, ",", index$column, "]")
  }), use.names = FALSE)
  additional <- unique(unlist(lapply(setup$compiled, `[[`, "add_parameters")))
  columns <- unique(c(review_prior_columns(priors), setdiff(additional, latent_bases), latent))
  samples <- matrix(.5, 24L, length(columns), dimnames = list(NULL, columns))
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  quantities <- BayesTools::parameter_catalog(object$fit)$quantities
  quantities <- quantities[quantities$quantity == "var", , drop = FALSE]
  available <- .brma_random_parameter_specs(quantities)
  selected <- available$parameter[available$quantity == "var" &
    grepl("tau2_study", available$label, fixed = TRUE)]
  expect_length(selected, 1L)
  label <- available$label[match(selected, available$parameter)]
  text <- capture.output(print_prior(object, parameter = selected))
  expect_identical(text[1L], paste0(label, " = tau^2"))
  expect_match(text[1L], "tau2_study", fixed = TRUE)
  expect_identical(text[2L], "  tau:")
})

test_that("allocation derivation preserves unexpected source errors", {

  error <- structure(list(message = "Synthetic source programming failure.", call = NULL),
    class = c("review_source_programming_error", "error", "condition"))
  testthat::local_mocked_bindings(
    .marginalized_random_effect_allocation_factors = function(...) stop(error),
    .marginalized_random_effect_has_allocation = function(...) TRUE,
    .brma_random_parameter_component_column = function(...) 1L,
    .package = "RoBMA")
  actual <- tryCatch(.brma_random_prior_derivation(
    object = list(fit = list()), spec = list(allocation_derived = TRUE),
    provenance = list(term = list(sd_binding = list(allocations = list(list()))))),
    error = identity)
  expect_identical(actual, error)
})
