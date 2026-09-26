exp_affine_test_sample <- function(mixture = FALSE) {

  condition_event <- structure(list(
    conditional      = character(),
    conditional_rule = "AND",
    families         = list(),
    condition_key    = "<averaged>"
  ), class = "BayesTools_condition_event")
  atoms <- BayesTools::posterior_atom_attribute(
    source = "product_space_structure"
  )
  sample_class <- c("mixed_posteriors.simple", "mixed_posteriors")
  if (isTRUE(mixture)) {
    sample_class <- c("mixed_posteriors.mixture", sample_class)
  }

  with_draw_metadata(
    c(.15, .25, .35),
    class     = sample_class,
    condition = list(
      conditional              = character(),
      condition_key            = "<averaged>",
      resolved_condition_event = condition_event,
      averaged                 = TRUE
    ),
    atoms     = atoms
  )
}


exp_affine_test_prior_density <- function() {

  structure(list(
    points = data.frame(x = numeric(), p = numeric())
  ), class = c("prior_linear_density", "prior_density"))
}


test_that("atom-free averaged exp-affine posteriors are certified", {

  sample <- exp_affine_test_sample(mixture = TRUE)

  expect_null(.hypothesis_plan_exp_affine_certify(
    sample        = sample,
    prior_density = exp_affine_test_prior_density()
  ))

  atomic <- sample
  BayesTools::posterior_metadata(atomic, "atoms") <- BayesTools::posterior_atom_attribute(
    point_masses = data.frame(x = 0, mass = .1)
  )
  refusal <- .hypothesis_plan_exp_affine_certify(
    sample        = atomic,
    prior_density = exp_affine_test_prior_density()
  )
  expect_match(refusal[["reason"]], "posterior is atom-free", fixed = TRUE)
})


test_that("exp-affine targets are tested on their own scale with the certified prior density", {

  density <- exp_affine_test_prior_density()
  # The target's marginal posterior (BayesTools::marginal_posterior() in the
  # plan): its draws with the prior density and the exact support of the map.
  target  <- with_draw_metadata(
    c(.15, .25, .35),
    class         = c("marginal_posterior.simple", "marginal_posterior"),
    prior_density = density,
    support       = BayesTools::posterior_support_attribute(c(0, Inf))
  )
  captured <- NULL
  testthat::local_mocked_bindings(
    hypothesis_BF = function(posterior, hypothesis, ...) {
      captured <<- list(posterior = posterior, hypothesis = hypothesis, ...)
      structure(
        data.frame(BF = 2, check.names = FALSE),
        hypothesis_ast = hypothesis, raw_BF = 2,
        class = c("BayesTools_hypothesis_BF", "data.frame")
      )
    },
    .package = "BayesTools"
  )
  plan <- list(
    parameter   = "log_tau_intercept",
    point       = TRUE,
    conditional = FALSE,
    draws       = list(posterior = target),
    targets     = list()
  )
  hypothesis <- BayesTools::hypothesis_parse("log_tau_intercept = 0.2")
  .hypothesis_plan_execute_scalar(
    plans           = list(plan),
    hypothesis      = hypothesis,
    object          = list(),
    logBF           = FALSE,
    BF01            = FALSE,
    seed            = 1,
    density_method  = "KDE",
    density_control = NULL,
    columns         = "all"
  )

  # No prior draws: the posterior carries the prior density, and the
  # statement keeps its original-scale values.
  expect_null(captured[["prior"]])
  expect_identical(captured[["density_method"]], "KDE")
  expect_identical(
    BayesTools::posterior_metadata(captured[["posterior"]], "prior_density"),
    density
  )
  expect_equal(as.numeric(captured[["posterior"]]), c(.15, .25, .35))
  expect_identical(captured[["hypothesis"]], hypothesis)
})


test_that("exp-affine coefficients refuse qCMDE/IWMDE with the original-scale classes", {

  skip_on_cran()
  # The original-scale scale intercept of a scale formula with a standardized
  # continuous predictor is an exp(affine) map of the fitted coefficients.
  # plot() and hypothesis() refuse qCMDE/IWMDE for it with the same
  # density-method classes; hypothesis() keeps its own classes first.
  set.seed(1)
  k   <- 30L
  dat <- data.frame(x = stats::rnorm(k), sei = stats::runif(k, 0.1, 0.3))
  dat[["yi"]] <- stats::rnorm(k, 0.2 + 0.1 * dat[["x"]], dat[["sei"]])
  fit <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ x, scale = ~ x, data = dat,
    measure = "SMD", chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  plan <- .hypothesis_plans(fit, "intercept = 0.5", component = "scale")[[1L]]
  expect_identical(plan[["kind"]], "exp_affine")

  classes <- c("RoBMA_density_method_original_scale",
               "RoBMA_density_method_unavailable")
  for (method in c("qCMDE", "IWMDE")) {
    error <- tryCatch(
      plot(fit, "intercept", component = "scale", density_method = method,
           density_control = list(n_points = 20, samples = 50)),
      error = identity
    )
    expect_identical(class(error), c(classes, "error", "condition"), info = method)
    expect_identical(
      conditionMessage(error),
      paste0(
        "qCMDE/IWMDE does not support the fitted nonlinear joint transform ",
        "for 'log_tau_intercept'. Use density_method = 'KDE' or ",
        "standardized_coefficients = TRUE."
      ),
      info = method
    )
    error <- tryCatch(
      hypothesis(fit, "intercept = 0.5", component = "scale",
                 density_method = method),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_method", "RoBMA_hypothesis_unavailable", classes,
        "error", "condition"),
      info = method
    )
  }
  # The normal approximation is no qCMDE/IWMDE request: its refusal keeps the
  # hypothesis classes only.
  error <- tryCatch(
    hypothesis(fit, "intercept = 0.5", component = "scale",
               density_method = "normal"),
    error = identity
  )
  expect_identical(
    class(error),
    c("RoBMA_hypothesis_method", "RoBMA_hypothesis_unavailable", "error",
      "condition")
  )
  kde <- suppressWarnings(hypothesis(
    fit, "intercept = 0.5", component = "scale", density_method = "KDE"
  ))
  expect_true(is.finite(attr(kde, "raw_BF")))
})
