context("Multilevel selection model for the IWMDE location fast path")

source(testthat::test_path("common-functions.R"))

skip_on_cran()
skip_if_not_installed("metadat")
skip_refit_if_cached("iwmde-multilevel-selection")


test_that("multilevel selection models integrating their random effects fit", {

  skip_if_not_certification(
    "This fit certifies the IWMDE multilevel weightfunction location fast path."
  )
  # The default selection model conditions on the cluster random effects, so
  # IWMDE evaluates their local conditional likelihood. Integrating them
  # makes every row's likelihood marginal: formula coefficients then take
  # the multilevel selected-normal location fast path, which
  # test-02-iwmde-fast-paths.R compares with the scalar evaluation.
  data(dat.lehmann2018, package = "metadat")
  fit <- bselmodel(
    yi = yi, vi = vi, mods = ~ Preregistered, cluster = Full_Citation,
    data = dat.lehmann2018, measure = "SMD",
    selection = selection_model(other_random_effects = "integrate"),
    chains = 2, sample = 1000, burnin = 500, adapt = 500,
    seed = 1, silent = TRUE
  )
  fit <- add_marglik(fit)
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018-3PSM_3lvl_mods_marginal", fit)

  expect_s3_class(fit, "bselmodel")
  expect_identical(
    .data_selection_model(fit[["data"]])[c(
      "estimate_random_effects", "other_random_effects", "known_sampling_variance"
    )],
    list(
      estimate_random_effects = "integrate",
      other_random_effects    = "integrate",
      known_sampling_variance = "integrate"
    )
  )
})
