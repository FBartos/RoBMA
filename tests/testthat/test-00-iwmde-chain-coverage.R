context("qCMDE/IWMDE refuse a fitted chain without eligible rows with a typed condition")

# A model-averaged effect (a normal slab or a spike at zero) replayed with two
# chains: the first visits the slab, the second never does.
.chain_coverage_replay <- function() {

  dat <- data.frame(yi = c(.1, .25, .3, .4, .2, .15, .05, .35), sei = .2)
  object <- BMA(
    yi = dat[["yi"]], sei = dat[["sei"]], measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_effect_null = BayesTools::prior("spike", list(0)),
    prior_heterogeneity = BayesTools::prior(
      "normal", list(0, 1), truncation = list(0, Inf)
    ),
    only_priors = TRUE, silent = TRUE
  )
  withr::local_seed(2905)
  n <- 300L
  chain <- function(slab) {
    cbind(
      mu           = if (slab) stats::rnorm(n, .3, .1) else rep(0, n),
      mu_indicator = if (slab) 2 else 1,
      tau          = abs(stats::rnorm(n, .15, .05)) + .02
    )
  }

  review_attach_replay(review_compile_object(object), list(chain(TRUE), chain(FALSE)))
}

test_that("a chain that never visits the conditioned target has a typed refusal with a remedy", {

  object <- .chain_coverage_replay()
  expect_identical(length(object[["fit"]][["mcmc"]]), 2L)

  for (method in c("qCMDE", "IWMDE")) {
    failure <- tryCatch(
      suppressWarnings(hypothesis(
        object, "mu = 0.3", density_method = method,
        density_control = list(n_points = 20, samples = 100), seed = 1
      )),
      error = identity
    )
    expect_s3_class(failure, "RoBMA_density_chain_coverage")
    # the public ordinate class and the construction classes are kept
    expect_s3_class(failure, "RoBMA_density_ordinate_error")
    expect_s3_class(failure, "iwmde_construction_error")
    expect_equal(failure[["stage"]], "chain-coverage validation")
    expect_equal(unname(failure[["chain_coverage"]][["eligible_chain_counts"]]), c(300L, 0L))
    expect_match(
      conditionMessage(failure), "no eligible continuous row remains in fitted chain(s): 2.",
      fixed = TRUE
    )
    # the message ends with the remedies
    expect_match(
      conditionMessage(failure),
      "Use density_method = 'KDE' to estimate the density from the posterior draws, or run longer chains and check convergence.$"
    )
  }

  # KDE has no chain gate
  expect_no_error(suppressWarnings(hypothesis(object, "mu = 0.3", density_method = "KDE")))
})

test_that("the density-only routes keep the remedy of a chain-coverage refusal", {

  construction_failure <- function(...) {
    tryCatch(
      .iwmde_stop_construction_failure("qCMDE", "mu", integer(), "stage", ...),
      error = identity
    )
  }
  coverage <- construction_failure(
    detail = "no eligible continuous row remains in fitted chain(s): 2",
    condition_class = "RoBMA_density_chain_coverage", remedy = "Remedy."
  )
  other <- construction_failure(detail = "invalid matrix")

  testthat::local_mocked_bindings(.iwmde_estimate = function(failure) stop(failure), .package = "RoBMA")
  wrapped <- tryCatch(.iwmde_estimate_density(coverage), error = identity)
  expect_identical(
    class(wrapped)[1:3],
    c("RoBMA_density_construction_error", "RoBMA_density_chain_coverage", "iwmde_construction_error")
  )
  # the remedy is not repeated by the generic KDE advice
  expect_identical(conditionMessage(wrapped), conditionMessage(coverage))

  # other construction failures still receive the generic advice
  wrapped <- tryCatch(.iwmde_estimate_density(other), error = identity)
  expect_false(inherits(wrapped, "RoBMA_density_chain_coverage"))
  expect_match(conditionMessage(wrapped), "Use density_method = 'KDE'", fixed = TRUE)
})
