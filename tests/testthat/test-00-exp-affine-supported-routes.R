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

  sample  <- exp_affine_test_sample()
  density <- exp_affine_test_prior_density()
  target  <- .hypothesis_plan_exp_affine_posterior(
    sample        = sample,
    prior_density = density,
    support       = c(0, Inf),
    parameter     = "log_tau_intercept"
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
