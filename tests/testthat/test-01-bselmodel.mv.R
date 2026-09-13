source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("bselmodel.mv")


test_that("bselmodel.mv fits the marginal product selection target", {

  data <- data.frame(
    yi    = c(0.05, 0.16, 0.22, 0.34, 0.08, 0.19),
    study = factor(rep(c("s1", "s2", "s3"), each = 2L)),
    x     = rep(c(0, 1), 3L)
  )
  V <- matrix(0, nrow = 6L, ncol = 6L)
  block <- matrix(c(0.010, 0.004, 0.004, 0.014), nrow = 2L)
  for (study in seq_len(3L)) {
    rows <- ((study - 1L) * 2L + 1L):(study * 2L)
    V[rows, rows] <- block
  }

  fit <- bselmodel.mv(
    yi = yi, V = V, mods = ~ x, random = ~ 1 | study,
    data = data, measure = "GEN",
    prior_bias = BayesTools::prior_weightfunction(
      side    = "one-sided",
      steps   = 0.025,
      weights = BayesTools::wf_cumulative(c(1, 1)),
      model   = BayesTools::selection_model(
        other_random_effects    = "integrate",
        known_sampling_variance = "integrate",
        group            = "study"
      )
    ),
    prior_unit_information_sd = 1,
    selection_control = set_selection_likelihood_control(
      points_per_scramble = 256L,
      scrambles            = 8L,
      relative_tolerance   = 0.02,
      seed                 = 41L
    ),
    chains = 2, sample = 300, burnin = 150, adapt = 500,
    seed = 194, silent = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
  fit <- suppressWarnings(add_loo(fit, unit = "estimate"))
  fit <- add_marglik(
    fit,
    repetitions = 1,
    maxiter     = 2000,
    silent      = TRUE
  )
  save_fit(
    "bselmodel.mv_marg_random",
    fit,
    info = list(data = data, V = V)
  )

  expect_identical(
    class(fit),
    c("bselmodel.mv", "bselmodel", "brma.mv", "brma.norm", "brma")
  )
  model <- .data_selection_model(fit[["data"]])
  expect_identical(
    model[c("other_random_effects", "known_sampling_variance")],
    list(other_random_effects = "integrate", known_sampling_variance = "integrate")
  )
  expect_identical(model[["groups"]][["group_index"]], rep(1:3, each = 2L))
  expect_identical(model[["branches"]][[1L]][["weight_rule"]], "product")
  expect_s3_class(fit[["loo"]][["estimate"]], "loo")
  expect_true(is.finite(as.numeric(logml(fit))))
})
