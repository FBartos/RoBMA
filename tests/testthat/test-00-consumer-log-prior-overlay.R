.consumer_log_scale_replay <- function() {

  dat <- data.frame(yi = c(-.1, .2, .3, .1), sei = rep(.2, 4L))
  object <- brma(yi = yi, sei = sei, scale = ~ 1, data = dat,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
    only_priors = TRUE, silent = TRUE)
  setup <- review_compile_object(object)
  columns <- unique(unlist(lapply(setup$compiled, function(compiled) {
    unlist(lapply(names(compiled$prior_list), function(parameter) {
      BayesTools:::.prior_linear_prior_columns(parameter, compiled$prior_list[[parameter]])
    }))
  })))
  base_priors <- .create_fit_priors(object$data, object$priors)
  stopifnot(identical(names(base_priors), "mu"))
  columns <- c(names(base_priors), columns)
  samples <- matrix(0, 80L, length(columns), dimnames = list(NULL, columns))
  samples[, "mu"] <-
    stats::qnorm((seq_len(80L) - .5) / 80) / 5
  samples[, review_coefficient_name(setup, "intercept", "log_tau")] <-
    exp(stats::qnorm((seq_len(80L) - .5) / 80) / 10)
  object <- review_attach_replay(setup, samples)
  # Only remove the prior-only dispatcher after attaching current posterior draws.
  class(object) <- setdiff(class(object), "only_priors.brma")
  object
}

test_that("posterior scale LOG overlays retain automatic and explicit display ranges", {

  object <- .consumer_log_scale_replay()
  expect_true(inherits(object, "brma"))
  expect_false(inherits(object, "only_priors.brma"))
  original <- serialize(object, NULL)
  for (limits in list(NULL, c(-3, 1))) {
    arguments <- list(x = object, parameter = "intercept", component = "scale",
      standardized_coefficients = TRUE, transform = "LOG", prior = TRUE,
      plot_type = "ggplot")
    if (!is.null(limits)) arguments$xlim <- limits
    graph <- do.call(plot, arguments)
    expect_s3_class(graph, "ggplot")
    # BayesTools draws the prior first, then overlays the posterior.
    prior <- ggplot2::ggplot_build(graph)$data[[1L]]
    expect_true(all(is.finite(prior$x)))
    expect_true(any(prior$x < 0))
    expect_equal(prior$y, 2 * stats::dnorm(exp(prior$x)) * exp(prior$x), tolerance = 1e-12)
    if (is.null(limits)) {
      expect_lt(min(prior$x), -5)
    } else {
      expect_equal(range(prior$x), limits, tolerance = 1e-12)
    }
  }
  expect_identical(serialize(object, NULL), original)
})
