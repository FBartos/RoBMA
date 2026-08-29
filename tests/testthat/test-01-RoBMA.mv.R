source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("RoBMA.mv")


.robma_mv_fit_input <- function() {

  data <- data.frame(
    yi    = c(0.05, 0.16, 0.22, 0.34, 0.08, 0.19),
    study = factor(rep(c("s1", "s2", "s3"), each = 2L)),
    obs   = factor(seq_len(6L)),
    x     = rep(c(0, 1), 3L)
  )
  V <- kronecker(
    diag(3L),
    matrix(c(0.010, 0.004, 0.004, 0.014), nrow = 2L)
  )

  list(data = data, V = V)
}


.robma_mv_fit_bias_priors <- function() {

  list(
    BayesTools::prior_weightfunction(
      "one-sided",
      steps   = 0.025,
      weights = BayesTools::wf_cumulative(c(1, 1))
    ),
    BayesTools::prior_PET(
      "cauchy",
      parameters = list(0, 1),
      truncation = list(0, Inf)
    ),
    BayesTools::prior_PEESE(
      "cauchy",
      parameters = list(0, 1),
      truncation = list(0, Inf)
    )
  )
}


test_that("RoBMA.mv fits the complete multivariate product space", {

  input <- .robma_mv_fit_input()
  fit <- RoBMA.mv(
    yi = yi, V = input[["V"]], mods = ~ x,
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    data = input[["data"]], measure = "GEN",
    prior_bias = .robma_mv_fit_bias_priors(),
    prior_unit_information_sd = 1,
    selection_control = set_selection_likelihood_control(
      points_per_scramble = 256L,
      scrambles            = 8L,
      relative_tolerance   = 0.02,
      seed                 = 47L
    ),
    chains = 2, sample = 300, burnin = 150, adapt = 500,
    seed = 196, silent = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
  fit <- suppressWarnings(add_loo(fit, unit = "estimate"))
  save_fit(
    "RoBMA.mv_exact_product_space",
    fit,
    info = input
  )

  expect_identical(
    class(fit),
    c("RoBMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma")
  )
  expect_identical(fit[["selection_likelihood"]][["type"]], "exact")
  expect_true(.is_PET(fit))
  expect_true(.is_PEESE(fit))
  expect_true(.is_weightfunction(fit))
  expect_length(summary(fit)[["inclusion_random"]][["post_prob"]], 2L)
  expect_s3_class(fit[["loo"]][["estimate"]], "loo")
})
