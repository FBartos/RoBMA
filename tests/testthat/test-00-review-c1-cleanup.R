test_that("bridge extraction distinguishes a missing raw upstream result", {

  # A supplied scalar is a real standardized result with no raw bridge object.
  object <- brma(yi = c(.1, -.2), sei = c(.2, .3), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE)
  object[["marglik"]] <- BayesTools::bridgesampling_object(-2)
  error <- tryCatch(bridge_sampler(object), error = identity)
  expect_s3_class(error, "RoBMA_bridge_unavailable")
  expect_identical(error[["reason"]], "missing_upstream")
  expect_identical(error[["logml"]], -2)
  expect_null(error[["call"]])
  expect_identical(conditionMessage(error), paste0(
    "The stored marginal-likelihood result does not contain a raw ",
    "bridge-sampling object. Recompute it with 'add_marglik()'."
  ))
})


test_that("bridge diagnostics name the declared sampling target", {

  for (sampling in c("integrate", "condition")) {
    object <- bselmodel(yi = c(.4, -.2), sei = c(.8, 1.1),
      cluster = c("a", "a"), measure = "GEN",
      prior_unit_information_sd = 1, only_priors = TRUE,
      selection = BayesTools::selection_model(
        known_sampling_variance = sampling, other_random_effects = "integrate"))
    priors <- .create_fit_priors(object[["data"]], object[["priors"]])
    setup <- .marglik_sampling_latent_setup(object[["data"]], object[["priors"]], priors)
    expected <- if (sampling == "integrate") {
      "full sampling error integrated before selection normalization"
    } else "full sampling error retained through Gaussian augmentation"
    expect_identical(setup[["diagnostics"]][["target"]], expected)
    if (sampling == "condition") {
      cluster <- .marglik_cluster_effects_setup(object[["data"]], object[["priors"]], priors)
      expect_false(cluster[["marginalized"]])
      expect_identical(cluster[["diagnostics"]][["reason"]],
        "cluster effects remain in the retained-sampling selection context")
    }
  }
})
