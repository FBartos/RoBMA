exp_affine_test_samples <- function(mixture = FALSE) {

  condition_event <- structure(list(
    conditional      = character(),
    conditional_rule = "AND",
    families         = list(),
    condition_key    = "<averaged>"
  ), class = "BayesTools_condition_event")
  atoms <- structure(list(
    declared                = TRUE,
    locations               = matrix(
      numeric(), nrow = 0L, ncol = 1L,
      dimnames = list(NULL, "log_tau_intercept")
    ),
    mass                     = numeric(),
    source                   = "product_space_structure",
    component_probabilities = c(.4, .6)
  ), class = c("BayesTools_posterior_atoms", "list"))
  sample_class <- c("mixed_posteriors.simple", "mixed_posteriors")
  if (isTRUE(mixture)) {
    sample_class <- c("mixed_posteriors.mixture", sample_class)
  }
  posterior <- structure(
    c(.15, .25, .35),
    class                    = sample_class,
    conditional              = character(),
    condition_key            = "<averaged>",
    resolved_condition_event = condition_event,
    posterior_atoms          = atoms
  )
  prior_density <- structure(list(
    points = data.frame(x = numeric(), p = numeric())
  ), class = c("prior_linear_density", "prior_density"))

  structure(
    list(log_tau_intercept = posterior),
    prior_densities = list(log_tau_intercept = prior_density)
  )
}


test_that("atom-free averaged exp-affine posteriors are certified", {

  samples <- exp_affine_test_samples(mixture = TRUE)

  expect_invisible(
    .hypothesis_brma_exp_affine_certify(
      samples     = samples,
      target      = "log_tau_intercept",
      conditional = FALSE
    )
  )

  atomic_samples <- samples
  atoms <- attr(
    atomic_samples[["log_tau_intercept"]],
    "posterior_atoms",
    exact = TRUE
  )
  atoms[["locations"]] <- matrix(
    0, nrow = 1L, ncol = 1L,
    dimnames = list("location", "log_tau_intercept")
  )
  atoms[["mass"]] <- .1
  attr(atomic_samples[["log_tau_intercept"]], "posterior_atoms") <- atoms
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      samples     = atomic_samples,
      target      = "log_tau_intercept",
      conditional = FALSE
    ),
    "posterior is atom-free"
  )
})


test_that("exp-affine point diagnostics are restored to the fitted scale", {

  samples <- exp_affine_test_samples()
  target <- list(
    target = "log_tau_intercept",
    route  = list(type = "exp_affine")
  )
  original <- BayesTools::hypothesis_parse("log_tau_intercept = 0.2")

  testthat::local_mocked_bindings(
    transform_prior_samples = function(...) {
      cbind(log_tau_intercept = c(.1, .2, .4))
    },
    hypothesis_BF = function(hypothesis, ...) {
      parsed <- list(list(
        input = "log_tau_intercept = -1.6094379124341003",
        left  = list(
          type  = "point",
          label = "log_tau_intercept = -1.6094379124341003",
          value = log(.2)
        ),
        right = list(
          type  = "not_point",
          label = "log_tau_intercept != -1.6094379124341003",
          value = log(.2)
        ),
        explicit = FALSE
      ))
      structure(data.frame(
        Alternative = "log_tau_intercept = -1.6094379124341003",
        Null        = "log_tau_intercept != -1.6094379124341003",
        BF          = 2,
        BF_error    = 3,
        prior       = 4,
        posterior   = 8,
        method      = "kernel Savage-Dickey",
        check.names = FALSE
      ), hypothesis_ast = hypothesis, parsed = parsed, raw_BF = 2)
    },
    .package = "BayesTools"
  )

  result <- .hypothesis_brma_exp_affine_kde(
    object      = list(fit = structure(list(), class = "BayesTools_fit")),
    samples     = samples,
    hypothesis  = original,
    parameter   = "log_tau_intercept",
    target_info = target,
    conditional = FALSE,
    logBF       = FALSE,
    BF01        = FALSE,
    seed        = 1,
    n_samples   = 3L,
    columns     = "all"
  )

  expect_equal(result[["prior"]], 4 / .2)
  expect_equal(result[["posterior"]], 8 / .2)
  expect_equal(result[["BF"]], 2)
  expect_equal(result[["BF_error"]], 3)
  expect_equal(attr(result, "raw_BF"), 2)
})
