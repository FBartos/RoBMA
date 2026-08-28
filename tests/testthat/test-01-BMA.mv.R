source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("BMA.mv")


.bma_mv_fit_data <- function() {

  data <- data.frame(
    yi    = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05, 0.11, 0.16),
    study = factor(rep(c("s1", "s2", "s3", "s4"), each = 2L)),
    obs   = factor(seq_len(8L)),
    x     = rep(c(0, 1), 4L)
  )
  V <- matrix(0, nrow = 8L, ncol = 8L)
  block <- matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2L)
  for (study in seq_len(4L)) {
    index <- ((study - 1L) * 2L + 1L):(study * 2L)
    V[index, index] <- block
  }

  list(data = data, V = V)
}


test_that("BMA.mv fits a product space over fixed and random components", {

  input <- .bma_mv_fit_data()
  fit <- BMA.mv(
    yi = yi, V = input[["V"]], mods = ~ x,
    random = list(study = ~ 1 | study, observation = ~ 1 | obs),
    data = input[["data"]], measure = "GEN",
    prior_heterogeneity = list(
      BayesTools::prior(
        "normal", parameters = list(mean = 0, sd = 0.25),
        truncation = list(0, Inf), prior_weights = 3
      ),
      BayesTools::prior(
        "normal", parameters = list(mean = 0, sd = 0.75),
        truncation = list(0, Inf), prior_weights = 2
      )
    ),
    prior_heterogeneity_null = BayesTools::prior(
      "spike", parameters = list(location = 0), prior_weights = 5
    ),
    prior_unit_information_sd = 1,
    chains = 2, sample = 300, burnin = 150, adapt = 500,
    seed = 193, silent = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit(
    "BMA.mv_random_components",
    fit,
    info = input
  )

  expect_identical(
    class(fit),
    c("BMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma")
  )
  expect_length(
    fit[["formula_design"]][["mu"]][["random_allocations"]][[1L]][["inclusion"]],
    2L
  )
  expect_length(summary(fit)[["inclusion_random"]][["post_prob"]], 2L)
  expect_s3_class(fit[["loo"]][["estimate"]], "loo")
})
