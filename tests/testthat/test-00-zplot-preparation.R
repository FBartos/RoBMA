.zplot_preparation_object <- function(kind = "normal") {

  dat <- data.frame(yi = rep(c(-.3, .1, .7), 4L), sei = rep(c(.2, .3, .4), 4L),
    study = rep(seq_len(4L), each = 3L))
  args <- list(yi = dat$yi, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  if (kind == "mv") {
    args$V <- diag(dat$sei^2)
    constructor <- bselmodel.mv
  } else {
    args$sei <- dat$sei
    constructor <- if (kind == "normal") brma else bselmodel
  }
  if (kind == "cluster") args$cluster <- dat$study
  do.call(constructor, args)
}


test_that("shared zplot preparation preserves complete summaries and target metadata", {

  control <- set_selection_likelihood_control(seed = 13L)
  estimates <- list(EDR = seq(.2, .4, length.out = 10L), weights = seq(1, 2, length.out = 10L))
  calls <- list()
  testthat::local_mocked_bindings(
    .zplot_fun.brma = function(object, ...) {
      calls[[length(calls) + 1L]] <<- list(...)
      estimates
    },
    .package = "RoBMA"
  )
  for (kind in c("normal", "selected", "cluster", "mv")) {
    object <- .zplot_preparation_object(kind)
    original <- object
    prepared <- .zplot_prepare(object, 1.5, 10, "marginal", control)
    result <- as_zplot(object, significance_level = 1.5, max_samples = 10,
      integration_control = control)
    expected <- original
    vector_target <- .zplot_vector_selection_target(object)
    expected[["zplot"]] <- list(estimates = list(EDR = estimates$EDR,
      weights = if (vector_target) NULL else estimates$weights),
      data = prepared$data, target = prepared$target)
    class(expected) <- c("zplot_brma", class(original))
    expect_identical(result, expected, info = kind)
    expect_identical(object, original, info = kind)
    expect_false(inherits(prepared, "zplot_brma"))
    expect_false("estimates" %in% names(prepared))
    expect_identical(prepared$data$z, .outcome_data_yi(object) / .outcome_data_sei(object))
    expect_identical(prepared$data$N_significant, sum(abs(prepared$data$z) > 1.5))
    expect_identical(prepared$data$N_observed, length(prepared$data$z))
    expect_identical(prepared$target$missing_count,
      if (vector_target) "unavailable" else "univariate_convention")
    expect_identical(tail(calls, 1L)[[1L]], list(z_threshold = 1.5,
      max_samples = 10L, extrapolate = TRUE, conditioning_depth = "marginal",
      integration_control = control))
    expect_s3_class(summary(result), "summary.zplot_brma")
  }
  expect_length(calls, 4L)
})


test_that("direct zplot retains both renderers and plot arguments without summary work", {

  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("vdiffr")
  control <- set_selection_likelihood_control(seed = 19L)
  summaries <- 0L
  density_calls <- list()
  # Ten deterministic density draws, including a normal and a displaced
  # component, exercise the common plot-data and credible-band calculations.
  density <- function(z) vapply(z, function(value) {
    stats::dnorm(value, mean = rep(c(0, .4), 5L), sd = seq(1, 1.3, length.out = 10L))
  }, numeric(10L))
  record_density <- function(object, z_sequence, max_samples, conditioning_depth,
      integration_control, parallel, cores) {
    density_calls[[length(density_calls) + 1L]] <<- list(
      data = object$zplot$data, target = object$zplot$target,
      max_samples = max_samples, conditioning_depth = conditioning_depth,
      integration_control = integration_control, parallel = parallel, cores = cores,
      grid = z_sequence)
    density(z_sequence)
  }
  testthat::local_mocked_bindings(
    .zplot_fun.brma = function(object, z_threshold = NULL, z_sequence = NULL,
        max_samples, conditioning_depth, integration_control, extrapolate,
        parallel = FALSE, cores = 1L) {
      if (!is.null(z_threshold)) {
        summaries <<- summaries + 1L
        return(list(EDR = rep(.3, 10L), weights = rep(2, 10L)))
      }
      value <- record_density(object, z_sequence, max_samples, conditioning_depth,
        integration_control, parallel, cores)
      if (extrapolate) value / 2 else value
    },
    .zplot_density_pair = function(object, z_sequence, max_samples,
        conditioning_depth, integration_control, parallel, cores) {
      value <- record_density(object, z_sequence, max_samples, conditioning_depth,
        integration_control, parallel, cores)
      list(fitted = value, extrapolated = value / 2)
    },
    .package = "RoBMA"
  )
  render <- function(fun, args, backend) {
    value <- NULL
    svg <- if (backend == "base") {
      file <- tempfile(fileext = ".svg")
      on.exit(unlink(file), add = TRUE)
      vdiffr::write_svg(function() {
        value <<- withVisible(do.call(fun, args))
      }, file, title = "Zplot preparation")
      readLines(file, warn = FALSE)
    } else {
      value <- withVisible(do.call(fun, args))
      NULL
    }
    list(result = value, svg = svg)
  }
  for (kind in c("normal", "selected", "cluster", "mv")) {
    object <- .zplot_preparation_object(kind)
    original <- object
    zobject <- as_zplot(object, significance_level = 1.5, max_samples = 10,
      integration_control = control)
    for (backend in c("base", "ggplot")) {
      for (extrapolate in c(FALSE, TRUE)) {
        args <- list(plot_type = backend, max_samples = 10, plot_extrapolation = extrapolate,
          from = -3, to = 4, length.out.hist = 12L, length.out.lines = 11L,
          probs = c(.1, .9), parallel = TRUE, cores = 2L,
          dots_fit = list(col = "red"), dots_extrapolation = list(col = "blue"))
        # Use the existing public object method as the comparison route.
        reference <- render(plot.zplot_brma, c(list(x = zobject), args), backend)
        reference_density <- tail(density_calls, 1L)[[1L]]
        before <- summaries
        direct <- render(zplot, c(list(object = object,
          significance_level = 1.5, summary_max_samples = 10,
          integration_control = control), args), backend)
        expect_identical(summaries, before)
        expect_identical(tail(density_calls, 1L)[[1L]], reference_density)
        expect_identical(object, original)
        if (backend == "base") {
          expect_identical(direct$result, list(value = NULL, visible = FALSE))
          expect_identical(direct$svg, reference$svg)
        } else {
          expect_identical(ggplot2::ggplot_build(direct$result$value)$data,
            ggplot2::ggplot_build(reference$result$value)$data)
        }
      }
    }
    before <- summaries
    explicit <- zplot(zobject, plot_type = "ggplot", max_samples = 10,
      length.out.lines = 11L)
    expect_s3_class(explicit, "ggplot")
    expect_identical(summaries, before)
  }
})


test_that("direct multivariate zplots match full summaries on ten replay draws", {

  testthat::skip_if_not_installed("ggplot2")
  for (selected in c(FALSE, TRUE)) {
    args <- list(yi = c(.1, -.2), V = matrix(c(.04, .01, .01, .09), 2L), mods = ~ 1,
      data = data.frame(study = c("a", "a")),
      prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
      measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
    if (selected) args$prior_bias <- BayesTools::prior_weightfunction(
      "one-sided", .025, BayesTools::wf_fixed(c(1, .2)))
    setup <- review_compile_object(do.call(if (selected) bselmodel.mv else brma.mv, args))
    samples <- matrix(seq(-.1, .3, length.out = 10L), ncol = 1L,
      dimnames = list(NULL, review_coefficient_name(setup, "intercept")))
    object <- review_attach_replay(setup, samples)
    zobject <- as_zplot(object, significance_level = 1.5, max_samples = 10L)
    reference <- plot(zobject, plot_type = "ggplot", max_samples = 10L,
      length.out.lines = 11L, plot_extrapolation = TRUE)
    direct <- zplot(object, significance_level = 1.5, summary_max_samples = 10L,
      plot_type = "ggplot", max_samples = 10L, length.out.lines = 11L,
      plot_extrapolation = TRUE)
    expect_identical(ggplot2::ggplot_build(direct)$data,
      ggplot2::ggplot_build(reference)$data)
  }
})


test_that("direct zplot retains validation and density failures but skips summary failures", {

  object <- .zplot_preparation_object()
  computations <- 0L
  fail_density <- FALSE
  testthat::local_mocked_bindings(
    .zplot_fun.brma = function(object, z_threshold = NULL, z_sequence = NULL, ...) {
      computations <<- computations + 1L
      if (!is.null(z_threshold)) stop("Summary-only integration failure.", call. = FALSE)
      if (fail_density) stop("Density integration failure.", call. = FALSE)
      matrix(stats::dnorm(z_sequence), 10L, length(z_sequence), byrow = TRUE)
    },
    .package = "RoBMA"
  )
  invalid <- list(list(significance_level = 0), list(significance_level = 39),
    list(summary_max_samples = 5L), list(summary_max_samples = NA_real_),
    list(conditioning_depth = "invalid"), list(conditioning_depth = "cluster"),
    list(integration_control = list(relative_tolerance = -1)),
    list(max_samples = 5L), list(parallel = NA), list(cores = 0L),
    list(probs = c(-1, 1)), list(plot_type = "invalid"))
  for (args in invalid) {
    expect_error(do.call(zplot, c(list(object = object), args)))
  }
  glmm <- object
  attr(glmm$data, "outcome_type") <- "bin"
  class(glmm) <- c("brma.glmm", "brma")
  expect_error(zplot(glmm), "GLMM objects are not supported", fixed = TRUE)
  expect_identical(computations, 0L)
  expect_error(as_zplot(object, max_samples = 10), "Summary-only integration failure.", fixed = TRUE)
  testthat::skip_if_not_installed("ggplot2")
  expect_s3_class(zplot(object, plot_type = "ggplot", max_samples = 10,
    summary_max_samples = Inf, length.out.lines = 11L), "ggplot")
  fail_density <- TRUE
  expect_error(zplot(object, plot_type = "ggplot", max_samples = 10,
    length.out.lines = 11L), "Density integration failure.", fixed = TRUE)
})


test_that("direct zplot preserves registered and global subclass conversions", {

  converted <- plotted <- NULL
  transform <- function(object, ...) {
    converted <<- list(object = object, arguments = list(...))
    structure(list(marker = "custom conversion"), class = "zplot_preparation_result")
  }
  render <- function(x, ...) {
    plotted <<- list(object = x, arguments = list(...))
    "custom plot result"
  }
  namespace <- asNamespace("RoBMA")
  conversion_table <- get(".__S3MethodsTable__.", envir = namespace)
  plot_table <- get(".__S3MethodsTable__.", envir = asNamespace("base"))
  registerS3method("as_zplot", "zplot_preparation_custom", transform, envir = namespace)
  registerS3method("plot", "zplot_preparation_result", render, envir = asNamespace("base"))
  withr::defer(rm(list = "as_zplot.zplot_preparation_custom", envir = conversion_table))
  withr::defer(rm(list = "plot.zplot_preparation_result", envir = plot_table))
  rlang::local_bindings(as_zplot.zplot_preparation_global = transform, .env = globalenv())
  testthat::local_mocked_bindings(
    .zplot_prepare = function(...) stop("Custom conversion must be used."),
    .package = "RoBMA"
  )
  control <- set_selection_likelihood_control(seed = 31L)
  for (subclass in c("zplot_preparation_custom", "zplot_preparation_global")) {
    object <- structure(list(marker = "custom model"),
      class = c("zplot_preparation_no_override", subclass, "brma"))
    result <- zplot(object, significance_level = 1.5, summary_max_samples = 13L,
      conditioning_depth = "estimate", integration_control = control,
      plot_type = "custom backend", custom_argument = 17L)
    expect_identical(result, "custom plot result")
    expect_identical(converted, list(object = object, arguments = list(
      significance_level = 1.5, max_samples = 13L, conditioning_depth = "estimate",
      integration_control = control)))
    expect_identical(plotted, list(object = structure(list(marker = "custom conversion"),
      class = "zplot_preparation_result"), arguments = list(
      plot_type = "custom backend", custom_argument = 17L)))
  }
})


test_that("subclasses using default conversion retain summary-free direct plotting", {

  object <- .zplot_preparation_object()
  class(object) <- c("zplot_preparation_no_override", class(object))
  # An ordinary plot subclass method was also preceded by zplot_brma in the
  # previous as_zplot() route; it must not divert the direct plot either.
  rlang::local_bindings(plot.zplot_preparation_no_override = function(...) {
    stop("Ordinary subclass plotting must not replace zplot plotting.")
  }, .env = globalenv())
  computations <- 0L
  testthat::local_mocked_bindings(
    .zplot_fun.brma = function(object, z_threshold = NULL, z_sequence = NULL, ...) {
      if (!is.null(z_threshold)) stop("Summary integration must not be used.")
      computations <<- computations + 1L
      matrix(stats::dnorm(z_sequence), 10L, length(z_sequence), byrow = TRUE)
    },
    .package = "RoBMA"
  )
  expect_s3_class(zplot(object, plot_type = "ggplot", max_samples = 10L,
    length.out.lines = 11L), "ggplot")
  expect_identical(computations, 1L)
})
