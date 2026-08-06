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
