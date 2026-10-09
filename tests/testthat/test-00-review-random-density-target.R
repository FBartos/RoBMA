context("Random-effect density targets follow the requested coefficient scale")

# A one-slope random SD of a standardized predictor with a 'tau' alias: the
# fitted SD coordinate is 'source' on the standardized scale and the original
# scale multiplies it by 1 / sd(x).
.review_random_slope_replay <- function() {

  withr::local_seed(2)
  k   <- 30L
  dat <- data.frame(
    study = rep(sprintf("s%02d", seq_len(10L)), each = 3L),
    x     = stats::rnorm(k)
  )
  dat[["yi"]] <- 0.2 + 0.1 * dat[["x"]] +
    stats::rnorm(10L, 0, 0.1)[as.integer(factor(dat[["study"]]))] +
    stats::rnorm(k, 0, 0.15)
  object <- brma.mv(
    yi = yi, V = diag(rep(.0225, k)), random = ~ id(0 + x | study), data = dat,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_heterogeneity = BayesTools::prior_random(
      sd = BayesTools::prior("normal", list(0, 1), truncation = list(0, Inf))
    ),
    standardize_continuous_predictors = TRUE, only_priors = TRUE, silent = TRUE
  )
  setup    <- review_compile_object(object)
  compiled <- setup[["compiled"]][["mu"]]
  columns  <- review_prior_columns(compiled[["prior_list"]])
  samples  <- matrix(.5, 100L, length(columns), dimnames = list(NULL, columns))
  sd_names <- attr(compiled[["formula_scale"]], "random_effect_sd_leaves")[[1L]][["leaf_names"]]
  samples[, sd_names] <- seq(.2, 1.6, length.out = 100L)
  samples[, review_coefficient_name(setup, "intercept")] <- seq(-.1, .1, length.out = 100L)

  list(object = review_attach_replay(setup, samples), samples = samples)
}

test_that("the random SD density target transforms the draws to the requested scale", {

  replay <- .review_random_slope_replay()
  object <- replay[["object"]]
  mixed <- lapply(c(original = FALSE, standardized = TRUE), function(standardized) {
    as.numeric(.brma_random_parameter_mixed_posterior(
      object, "tau", standardized_coefficients = standardized
    )[[1L]])
  })
  # the scales differ, so the targets must differ as well
  expect_gt(max(abs(mixed[["original"]] - mixed[["standardized"]])), .1)

  for (standardized in c(FALSE, TRUE)) {
    target <- .brma_random_parameter_density_target(
      object, "tau", standardized_coefficients = standardized
    )
    transformed <- BayesTools::parameter_transform_forward(
      replay[["samples"]][, target[["parameter"]]], target[["display_transform"]]
    )
    expect_equal(
      as.numeric(transformed),
      mixed[[if (standardized) "standardized" else "original"]],
      tolerance = 1e-12, info = paste("standardized:", standardized)
    )
  }
  # the default is the original scale
  default <- .brma_random_parameter_density_target(object, "tau")
  expect_identical(default[["display_transform"]],
                   .brma_random_parameter_density_target(
                     object, "tau", standardized_coefficients = FALSE
                   )[["display_transform"]])
})

test_that("random SD point hypotheses plan their density target on the statement scale", {

  replay <- .review_random_slope_replay()
  object <- replay[["object"]]

  for (standardized in c(FALSE, TRUE)) {
    plan <- .hypothesis_plans(object, "tau = .5", standardized = standardized)[[1L]]
    target <- plan[["density_target"]]
    expect_false(is.null(target), info = paste("standardized:", standardized))
    transformed <- BayesTools::parameter_transform_forward(
      replay[["samples"]][, target[["parameter"]]], target[["display_transform"]]
    )
    expect_equal(
      as.numeric(transformed), as.numeric(plan[["draws"]][["posterior"]]),
      tolerance = 1e-12, info = paste("standardized:", standardized)
    )
  }
})

test_that("random SD plots build their density target on the plotted coefficient scale", {

  replay <- .review_random_slope_replay()
  object <- replay[["object"]]
  # the replayed posterior is a fitted posterior, plotted by the posterior method
  class(object) <- setdiff(class(object), "only_priors.brma")

  for (standardized in c(FALSE, TRUE)) {
    attached <- NULL
    stop_after_attach <- function(...) {
      stop(structure(
        class = c("review_plot_target_built", "error", "condition"),
        list(message = "plot target built", call = NULL)
      ))
    }
    outcome <- testthat::with_mocked_bindings(
      tryCatch(
        plot(object, "tau", standardized_coefficients = standardized,
             density_method = "qCMDE"),
        review_plot_target_built = identity
      ),
      .iwmde_check_density_method_supported = function(...) invisible(NULL),
      .plot_brma_attach_iwmde = function(object, samples, parameter,
                                         display_transform = NULL, ...) {
        attached <<- list(parameter = parameter, display_transform = display_transform)
        stop_after_attach()
      },
      .package = "RoBMA"
    )
    info <- paste("standardized:", standardized)
    expect_s3_class(outcome, "review_plot_target_built")
    expect_false(is.null(attached), info = info)

    # the target of the plot is the mixed posterior drawn on the plotted scale
    plotted <- as.numeric(.brma_random_parameter_mixed_posterior(
      object, "tau", standardized_coefficients = standardized
    )[[1L]])
    transformed <- BayesTools::parameter_transform_forward(
      replay[["samples"]][, attached[["parameter"]]], attached[["display_transform"]]
    )
    expect_equal(as.numeric(transformed), plotted, tolerance = 1e-12, info = info)
  }
})
